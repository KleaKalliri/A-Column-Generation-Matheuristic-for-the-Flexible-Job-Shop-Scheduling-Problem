// =============================================================================
// branch_and_price.cpp  –  Branch-and-Price for FJSSP
//
// Implements branch_and_price.h.
// Mirrors branch_and_price_machine_based.py.
// =============================================================================

#include "branch_and_price.h"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/sat_parameters.pb.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sat = operations_research::sat;

static constexpr double DINF = std::numeric_limits<double>::infinity();
static constexpr int    IINF = std::numeric_limits<int>::max() / 2;

// ===========================================================================
//  CP-SAT solver
// ===========================================================================

CPSatResult cp_sat_solve(const Instance& inst,
                          const BranchAndPriceNode& node,
                          double time_limit_sec,
                          double ub,
                          const Schedule* hint_schedule)
{
    CPSatResult result;
    int J = inst.num_jobs;

    sat::CpModelBuilder model;
    int horizon = 0;
    for (int j = 0; j < J; ++j)
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            for (auto& [m, pt] : inst.proc_time[j][o])
                horizon += pt;
        }

    // Determine eligible machines per op under node constraints
    auto get_eligible = [&](int j, int o) -> std::vector<int> {
        Op jop{j, o};
        std::vector<int> result;

        for (int m : inst.elig_mach[j][o]) {
            // required_machine
            auto it_req = node.required_machine.find(jop);
            if (it_req != node.required_machine.end() && it_req->second != m) continue;
            // forbidden
            auto it_for = node.forbidden_machines.find(jop);
            if (it_for != node.forbidden_machines.end() && it_for->second.count(m)) continue;
            // allowed
            auto it_all = node.allowed_machines.find(jop);
            if (it_all != node.allowed_machines.end() && !it_all->second.count(m)) continue;

            result.push_back(m);
        }
        return result;
    };

    // Minimum permitted processing time at the current B&P node.
    std::vector<std::vector<int>> min_pt_allowed(J);
    std::vector<std::vector<int>> job_earliest(J);
    std::vector<std::vector<int>> job_remaining(J);

    for (int j = 0; j < J; ++j) {
        const int n = inst.num_ops[j];

        min_pt_allowed[j].resize(n, IINF);
        job_earliest[j].resize(n, 0);
        job_remaining[j].resize(n + 1, 0);

        for (int o = 0; o < n; ++o) {
            const auto eligs = get_eligible(j, o);

            if (eligs.empty()) {
                return result;
            }

            for (int m : eligs) {
                min_pt_allowed[j][o] = std::min(
                    min_pt_allowed[j][o],
                    inst.proc_time[j][o].at(m)
                );
            }

            if (min_pt_allowed[j][o] == IINF) {
                return result;
            }
        }

        // Earliest start from preceding operations.
        for (int o = 1; o < n; ++o) {
            job_earliest[j][o] =
                job_earliest[j][o - 1]
                + min_pt_allowed[j][o - 1];
        }

        // Minimum work from operation o through the final operation.
        for (int o = n - 1; o >= 0; --o) {
            job_remaining[j][o] =
                min_pt_allowed[j][o]
                + job_remaining[j][o + 1];
        }
    }

    // Variables
    // interval_vars[(j,o,m)] = optional interval on machine m
    std::map<std::tuple<int,int,int>, sat::IntervalVar> intervals;
    std::map<std::tuple<int,int,int>, sat::BoolVar>     presence;
    std::map<std::pair<int,int>, sat::IntVar>           start_vars;
    std::map<std::pair<int,int>, sat::IntVar>           end_vars;

    for (int j = 0; j < J; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            Op jop{j, o};
            auto eligs = get_eligible(j, o);
            if (eligs.empty()) return result;  // infeasible

            // int slb = 0, sub = horizon;
            // {
            //     auto it = node.start_lb.find(jop);
            //     if (it != node.start_lb.end()) slb = it->second;
            // }
            // {
            //     auto it = node.start_ub.find(jop);
            //     if (it != node.start_ub.end()) sub = it->second;
            // }

            int slb = job_earliest[j][o];

            auto it_lb = node.start_lb.find(jop);
            if (it_lb != node.start_lb.end()) {
                slb = std::max(slb, it_lb->second);
            }

            // The current operation and all successors must fit before H.
            int sub = H - job_remaining[j][o];

            auto it_ub = node.start_ub.find(jop);
            if (it_ub != node.start_ub.end()) {
                sub = std::min(sub, it_ub->second);
            }

            if (slb > sub) {
                return result;
            }

            auto sv = model.NewIntVar({slb, sub});
            // end[j,o] = sv + sum_m presence[j,o,m] * pt[j,o,m]
            // (exactly one machine is selected, so this linear expr equals sv + pt_assigned)
            auto ev = model.NewIntVar({slb, horizon});
            start_vars[{j,o}] = sv;
            end_vars[{j,o}]   = ev;

            // Optional intervals for each eligible machine
            std::vector<sat::BoolVar> mach_pres;
            sat::LinearExpr pt_sum;   // sum_m pv_m * pt_m
            for (int m : eligs) {
                int pt = inst.proc_time[j][o].at(m);
                auto pv = model.NewBoolVar();
                auto iv = model.NewOptionalFixedSizeIntervalVar(sv, (int64_t)pt, pv);
                intervals[{j,o,m}] = iv;
                presence[{j,o,m}]  = pv;
                mach_pres.push_back(pv);
                // Accumulate: presence_var * proc_time
                pt_sum += sat::LinearExpr::Term(pv, (int64_t)pt);
            }

            // Exactly one machine assigned
            model.AddExactlyOne(mach_pres);

            // ev = sv + sum_m pv_m * pt_m  =>  ev - sv = pt_sum
            sat::LinearExpr lhs = sat::LinearExpr(ev) - sat::LinearExpr(sv);
            model.AddEquality(lhs, pt_sum);
        }
    }

    // NoOverlap per machine
    for (int m : inst.all_machines) {
        std::vector<sat::IntervalVar> mach_intervals;
        for (int j = 0; j < J; ++j)
            for (int o = 0; o < inst.num_ops[j]; ++o)
                if (intervals.count({j,o,m}))
                    mach_intervals.push_back(intervals.at({j,o,m}));
        if (!mach_intervals.empty())
            model.AddNoOverlap(mach_intervals);
    }

    // Job precedence
    for (int j = 0; j < J; ++j) {
        for (int o = 0; o < inst.num_ops[j] - 1; ++o) {
            model.AddLessOrEqual(end_vars[{j,o}], start_vars[{j,o+1}]);
        }
    }

    // Order constraints
    // for (auto& [m, oc_set] : node.order_constraints) {
    //     for (auto& [A, B] : oc_set) {
    //         auto [j1,o1] = A;
    //         auto [j2,o2] = B;
    //         // Check if both are on machine m
    //         auto ita = presence.find({j1,o1,m});
    //         auto itb = presence.find({j2,o2,m});
    //         if (ita != presence.end() && itb != presence.end()) {
    //             // end(A) <= start(B) when both are on m
    //             // Since order constraint means A before B on same machine:
    //             model.AddLessOrEqual(end_vars[{j1,o1}], start_vars[{j2,o2}])
    //                 .OnlyEnforceIf({ita->second, itb->second});
    //         }
    //     }
    // }

    // Makespan
    sat::IntVar cmax = model.NewIntVar({0, horizon});
    for (int j = 0; j < J; ++j) {
        int last_o = inst.num_ops[j] - 1;
        model.AddLessOrEqual(end_vars[{j, last_o}], cmax);
    }

    // If an incumbent UB is available, only search for strict improvements.
    // For integer makespan values, cmax < ceil(ub) is equivalent to cmax <= ub - 1.
    if (std::isfinite(ub)) {
        int64_t ub_cutoff = static_cast<int64_t>(std::ceil(ub));
        model.AddLessOrEqual(cmax, ub_cutoff - 1);
    }

    // Optional warm-start hint. Provide operation start/end hints and a
    // one-hot machine-assignment hint for optional intervals.
    if (hint_schedule != nullptr) {
        for (int j = 0; j < J; ++j) {
            for (int o = 0; o < inst.num_ops[j]; ++o) {
                Op jop{j, o};
                auto it_h = hint_schedule->find(jop);
                if (it_h == hint_schedule->end()) continue;

                auto [hint_m, hint_s, hint_e] = it_h->second;

                auto it_sv = start_vars.find({j, o});
                if (it_sv != start_vars.end())
                    model.AddHint(it_sv->second, hint_s);

                auto it_ev = end_vars.find({j, o});
                if (it_ev != end_vars.end())
                    model.AddHint(it_ev->second, hint_e);

                auto eligs = get_eligible(j, o);
                for (int m : eligs) {
                    auto it_p = presence.find({j, o, m});
                    if (it_p == presence.end()) continue;
                    model.AddHint(it_p->second, (m == hint_m) ? 1 : 0);
                }
            }
        }
    }

    model.Minimize(cmax);

    sat::CpSolverResponse response;
    {
        sat::SatParameters params;
        params.set_max_time_in_seconds(time_limit_sec);
        params.set_log_search_progress(false);
        response = sat::SolveWithParameters(model.Build(), params);
    }

    if (response.status() == sat::CpSolverStatus::OPTIMAL ||
        response.status() == sat::CpSolverStatus::FEASIBLE)
    {
        int64_t cms = sat::SolutionIntegerValue(response, cmax);
        result.makespan = (double)cms;
        result.found    = true;

        for (int j = 0; j < J; ++j) {
            for (int o = 0; o < inst.num_ops[j]; ++o) {
                int64_t sv = sat::SolutionIntegerValue(response, start_vars[{j,o}]);
                // Find assigned machine
                for (int m : get_eligible(j, o)) {
                    if (presence.count({j,o,m})) {
                        bool pv = sat::SolutionBooleanValue(response, presence[{j,o,m}]);
                        if (pv) {
                            int pt = inst.proc_time[j][o].at(m);
                            result.schedule[{j,o}] = {m, (int)sv, (int)(sv+pt)};
                            break;
                        }
                    }
                }
            }
        }
    }

    return result;
}

// ===========================================================================
//  Propagation rules
// ===========================================================================

static int compute_earliest_start(int j, int o, const Instance& inst,
                                   const std::map<Op, int>& start_lb)
{
    // Critical path lower bound: sum of min processing times of preceding ops
    int lb = 0;
    for (int k = 0; k < o; ++k) {
        int mn = IINF;
        for (auto& [m, pt] : inst.proc_time[j][k]) mn = std::min(mn, pt);
        lb += mn;
    }
    auto it = start_lb.find({j, o});
    if (it != start_lb.end()) lb = std::max(lb, it->second);
    return lb;
}

/// Transitive closure of order constraints (Floyd-Warshall style on (j,o) pairs)
// static bool transitive_closure(
//     std::map<int, std::set<std::pair<Op,Op>>>& order_constraints)
// {
//     bool changed = true;
//     while (changed) {
//         changed = false;
//         for (auto& [m, oc] : order_constraints) {
//             std::vector<std::pair<Op,Op>> new_pairs;
//             for (auto& [A, B] : oc) {
//                 for (auto& [C, D] : oc) {
//                     if (B == C && !oc.count({A, D})) {
//                         new_pairs.push_back({A, D});
//                     }
//                 }
//             }
//             // Check cycles: if A->B and B->A, infeasible
//             for (auto& [A, B] : oc) {
//                 if (oc.count({B, A})) return false;  // cycle
//             }
//             for (auto& p : new_pairs) {
//                 if (oc.insert(p).second) changed = true;
//             }
//         }
//     }
//     return true;
// }

bool propagate_node(BranchAndPriceNode& node, const Instance& inst, double best_ub)
{
    bool changed = true;
    int  iters   = 0;
    const int MAX_PROP_ITER = 30;

    while (changed && iters++ < MAX_PROP_ITER) {
        changed = false;
        int J = inst.num_jobs;

        // R1: OMA domain filtering
        //   - Intersect allowed_machines with elig_mach
        //   - Remove forbidden_machines from domain
        //   - If exactly one machine left, fix it
        for (int j = 0; j < J; ++j) {
            for (int o = 0; o < inst.num_ops[j]; ++o) {
                Op jop{j, o};

                // Build domain
                std::set<int> domain(inst.elig_mach[j][o].begin(), inst.elig_mach[j][o].end());

                // Intersect with allowed
                auto it_all = node.allowed_machines.find(jop);
                if (it_all != node.allowed_machines.end()) {
                    std::set<int> inter;
                    for (int m : it_all->second)
                        if (domain.count(m)) inter.insert(m);
                    if (inter != domain) { domain = inter; changed = true; }
                    if (domain.empty()) return false;
                    node.allowed_machines[jop] = domain;
                }

                // Remove forbidden
                auto it_for = node.forbidden_machines.find(jop);
                if (it_for != node.forbidden_machines.end()) {
                    for (int m : it_for->second) {
                        if (domain.erase(m)) changed = true;
                    }
                    if (domain.empty()) return false;
                }

                // Check required
                auto it_req = node.required_machine.find(jop);
                if (it_req != node.required_machine.end()) {
                    if (!domain.count(it_req->second)) return false;
                    domain = {it_req->second};
                }

                // If single machine in domain, fix it
                if (domain.size() == 1) {
                    int m = *domain.begin();
                    auto it = node.required_machine.find(jop);
                    if (it == node.required_machine.end()) {
                        node.required_machine[jop] = m;
                        changed = true;
                    }
                }
            }
        }

        // R2: Time-window domain pruning (forward + backward pass)
        if (best_ub < 1e29) {
            for (int j = 0; j < J; ++j) {
                int n_ops = inst.num_ops[j];

                // Helper: effective allowed domain for op (j,o)
                auto get_domain = [&](int jj, int oo) -> std::set<int> {
                    std::set<int> dom(inst.elig_mach[jj][oo].begin(), inst.elig_mach[jj][oo].end());
                    auto ita = node.allowed_machines.find({jj, oo});
                    if (ita != node.allowed_machines.end()) {
                        std::set<int> inter;
                        for (int m : ita->second) if (dom.count(m)) inter.insert(m);
                        dom = inter;
                    }
                    auto itf = node.forbidden_machines.find({jj, oo});
                    if (itf != node.forbidden_machines.end())
                        for (int m : itf->second) dom.erase(m);
                    auto itr = node.required_machine.find({jj, oo});
                    if (itr != node.required_machine.end())
                        dom = {itr->second};
                    return dom;
                };

                // Forward pass: earliest start per op
                std::vector<int> es(n_ops, 0);
                for (int o = 1; o < n_ops; ++o) {
                    auto dom = get_domain(j, o - 1);
                    if (dom.empty()) return false;
                    int min_pt = IINF;
                    for (int m : dom) min_pt = std::min(min_pt, inst.proc_time[j][o-1].at(m));
                    es[o] = es[o - 1] + min_pt;
                }

                // Backward pass: latest finish per op
                std::vector<int> lf(n_ops, (int)std::floor(best_ub));
                for (int o = n_ops - 2; o >= 0; --o) {
                    auto dom = get_domain(j, o + 1);
                    if (dom.empty()) return false;
                    int min_pt = IINF;
                    for (int m : dom) min_pt = std::min(min_pt, inst.proc_time[j][o+1].at(m));
                    lf[o] = lf[o + 1] - min_pt;
                }

                // Update start_lb from forward pass; prune machines that cannot fit
                for (int o = 0; o < n_ops; ++o) {
                    Op jop{j, o};

                    // Update start_lb
                    auto it_lb = node.start_lb.find(jop);
                    int cur_lb = (it_lb != node.start_lb.end()) ? it_lb->second : 0;
                    int new_lb = std::max(cur_lb, es[o]);
                    if (new_lb > cur_lb) {
                        node.start_lb[jop] = new_lb;
                        changed = true;
                    }

                    // Prune machines where es[o] + pt > lf[o]
                    auto dom = get_domain(j, o);
                    if (dom.empty()) return false;
                    std::set<int> pruned;
                    for (int m : dom) {
                        int pt = inst.proc_time[j][o].at(m);
                        if (es[o] + pt > lf[o] + 1e-9)
                            pruned.insert(m);
                    }
                    if (!pruned.empty()) {
                        std::set<int> new_dom;
                        for (int m : dom) if (!pruned.count(m)) new_dom.insert(m);
                        if (new_dom.empty()) return false;
                        node.allowed_machines[jop] = new_dom;
                        changed = true;
                        // Singleton -> fix required_machine
                        if (new_dom.size() == 1) {
                            int sole = *new_dom.begin();
                            if (!node.required_machine.count(jop)) {
                                node.required_machine[jop] = sole;
                                changed = true;
                            }
                        }
                    }
                }
            }
        } else {
            // No UB yet: only do forward start_lb update
            for (int j = 0; j < J; ++j) {
                int lb = 0;
                for (int o = 0; o < inst.num_ops[j]; ++o) {
                    Op jop{j, o};
                    auto it = node.start_lb.find(jop);
                    int cur_lb = (it != node.start_lb.end()) ? it->second : 0;
                    int new_lb = std::max(cur_lb, lb);
                    if (new_lb > cur_lb) { node.start_lb[jop] = new_lb; changed = true; }
                    std::set<int> domain(inst.elig_mach[j][o].begin(), inst.elig_mach[j][o].end());
                    auto itf = node.forbidden_machines.find(jop);
                    if (itf != node.forbidden_machines.end())
                        for (int m : itf->second) domain.erase(m);
                    auto ita = node.allowed_machines.find(jop);
                    if (ita != node.allowed_machines.end()) {
                        std::set<int> inter;
                        for (int m : ita->second) if (domain.count(m)) inter.insert(m);
                        domain = inter;
                    }
                    if (domain.empty()) return false;
                    int min_pt = IINF;
                    for (int m : domain) min_pt = std::min(min_pt, inst.proc_time[j][o].at(m));
                    lb = new_lb + min_pt;
                }
            }
        }

        // R7: CP lower bounds (job critical path)
        for (int j = 0; j < J; ++j) {
            int lb = 0;
            for (int o = 0; o < inst.num_ops[j]; ++o) {
                Op jop{j, o};
                int cur = 0;
                auto it = node.start_lb.find(jop);
                if (it != node.start_lb.end()) cur = it->second;
                lb = std::max(lb, cur);
                int min_pt = IINF;
                for (int m : inst.elig_mach[j][o]) {
                    auto it2 = node.forbidden_machines.find(jop);
                    bool forb = (it2 != node.forbidden_machines.end() && it2->second.count(m));
                    if (!forb) min_pt = std::min(min_pt, inst.proc_time[j][o].at(m));
                }
                if (min_pt == IINF) return false;
                lb += min_pt;
            }
            // lb is a lower bound on makespan for this job
            // Use to update node.cp_lower_bound
            if ((double)lb > node.cp_lower_bound) {
                node.cp_lower_bound = (double)lb;
                changed = true;
            }
        }
    }

    return true;
}

// ===========================================================================
//  Column feasibility check
// ===========================================================================

bool column_is_feasible_at_node(const Column& col, int machine,
                                  const BranchAndPriceNode& node,
                                  const Instance& inst)
{
    // Check each op in the column
    for (auto& [jop, s] : col.s_vals) {
        auto [j, o] = jop;
        int e = col.e_vals.at(jop);

        // required_machine
        auto it_req = node.required_machine.find(jop);
        if (it_req != node.required_machine.end() && it_req->second != machine)
            return false;

        // forbidden_machines
        auto it_for = node.forbidden_machines.find(jop);
        if (it_for != node.forbidden_machines.end() && it_for->second.count(machine))
            return false;

        // allowed_machines
        auto it_all = node.allowed_machines.find(jop);
        if (it_all != node.allowed_machines.end() && !it_all->second.count(machine))
            return false;

        // start bounds
        auto it_lb = node.start_lb.find(jop);
        if (it_lb != node.start_lb.end() && s < it_lb->second) return false;

        auto it_ub = node.start_ub.find(jop);
        if (it_ub != node.start_ub.end() && s > it_ub->second) return false;
    }

    // Order constraints: only applied when check_order=true (not for inherited pool)
    // if (check_order) {
    //     auto it_oc = node.order_constraints.find(machine);
    //     if (it_oc != node.order_constraints.end()) {
    //         for (auto& [A, B] : it_oc->second) {
    //             if (col.a.count(A) && col.a.count(B)) {
    //                 int eA = col.e_vals.at(A);
    //                 int sB = col.s_vals.at(B);
    //                 if (eA > sB) return false;
    //             }
    //         }
    //     }
    // }

    return true;
}

// ===========================================================================
//  B&P solver implementation
// ===========================================================================

BranchAndPriceSolver::BranchAndPriceSolver(const Instance& inst)
    : inst_(inst), horizon_(compute_horizon(inst))
{}

bool BranchAndPriceSolver::update_ub(double makespan, const Schedule& sched)
{
    if (makespan < best_ub_ - BOUND_EPS) {
        best_ub_       = makespan;
        best_schedule_ = sched;
        std::cout << "[B&P] New incumbent: " << makespan << "\n";
        return true;
    }
    return false;
}

bool BranchAndPriceSolver::is_integral(
    const BranchAndPriceNode&,
    const std::map<std::pair<int,int>, double>& y_values,
    const std::map<int, std::vector<Column>>& cols)
{
    // 1. Every non-zero y must be near 1
    for (auto& [mp, y] : y_values)
        if (y > LP_POS_TOL && y < 1.0 - INT_TOL)
            return false;

    // 2. Each op (j,o) must be covered by exactly one y=1 column
    for (int j = 0; j < inst_.num_jobs; ++j)
        for (int o = 0; o < inst_.num_ops[j]; ++o) {
            double coverage = 0.0;
            for (auto& [mp, y] : y_values) {
                auto [m, p] = mp;
                if (y > 1.0 - INT_TOL && cols.at(m)[p].a.count({j, o}))
                    coverage += y;
            }
            if (std::abs(coverage - 1.0) > INT_TOL) return false;
        }
    return true;
}

Schedule BranchAndPriceSolver::extract_schedule(
    const BranchAndPriceNode&,
    const std::map<std::pair<int,int>, double>& y_values,
    const std::map<int, std::vector<Column>>& cols)
{
    Schedule sched;
    for (auto& [mp, y] : y_values) {
        auto [m, p] = mp;
        if (y > 1.0 - INT_TOL) {
            const Column& col = cols.at(m)[p];
            for (auto& [jop, s] : col.s_vals) {
                int e = col.e_vals.at(jop);
                sched[jop] = {m, s, e};
            }
        }
    }
    return sched;
}

std::optional<BranchDecision>
BranchAndPriceSolver::select_branch(
    const BranchAndPriceNode& node,
    const std::map<std::pair<int,int>, double>& y_values,
    const std::map<int, std::vector<Column>>& cols)
{

    // Priority 2: OMA branching — fallback when no fractional order pair exists
    // For each op (j,o), compute mass[m] = sum of y[m,p] where (j,o) is in col p on m
    std::map<Op, std::map<int, double>> op_mass;
    for (auto& [mp, y] : y_values) {
        if (y < LP_POS_TOL) continue;
        auto [m, p] = mp;
        auto& col = cols.at(m)[p];
        for (auto& [jop, _] : col.s_vals) {
            op_mass[jop][m] += y;
        }
    }

    Op best_op = {-1, -1};
    double best_score = -1.0;
    int    best_domain = INT_MAX;  // MRV tie-break: smallest remaining domain

    // Helper: effective domain size for op (j,o)
    auto domain_size = [&](const Op& jop) -> int {
        auto [jj, oo] = jop;
        std::set<int> dom(inst_.elig_mach[jj][oo].begin(), inst_.elig_mach[jj][oo].end());
        auto ita = node.allowed_machines.find(jop);
        if (ita != node.allowed_machines.end()) {
            std::set<int> inter;
            for (int mm : ita->second) if (dom.count(mm)) inter.insert(mm);
            dom = inter;
        }
        auto itf = node.forbidden_machines.find(jop);
        if (itf != node.forbidden_machines.end())
            for (int mm : itf->second) dom.erase(mm);
        return (int)dom.size();
    };

    for (auto& [jop, mass_map] : op_mass) {
        // Already fixed?
        if (node.required_machine.count(jop)) continue;

        // Score = total fractional weight across all machines
        double score = 0.0;
        for (auto& [m, wt] : mass_map)
            score += std::min(wt, 1.0 - wt);

        int dom = domain_size(jop);
        if (score > best_score + 1e-9 ||
            (score > best_score - 1e-9 && dom < best_domain)) {
            best_score  = score;
            best_domain = dom;
            best_op     = jop;
        }
    }

    if (best_op.first >= 0 && best_score > INT_TOL) {
        // auto& mass_map = op_mass[best_op];
        // int   fix_m    = -1;
        // double fix_w   = -1;
        // for (auto& [m, w] : mass_map) {
        //     if (w > fix_w) { fix_w = w; fix_m = m; }
        // }
        // if (fix_m >= 0)
           return BranchDecision::oma_fix(best_op, fix_m);
    }

    return std::nullopt;
}

std::vector<BranchAndPriceNode>
BranchAndPriceSolver::process_node(BranchAndPriceNode& node)
{
    // Filter inherited pool: OMA/time-bound checks only (no order constraints).
    // Order-violating columns will be excluded naturally by CG pricing; keeping
    // them here ensures the RMP has enough warmstart columns to stay feasible.
    std::map<int, std::vector<Column>> feas_cols;
    for (int m : inst_.all_machines) {
        feas_cols[m] = {};
        for (auto& col : node.cols[m]) {
            if (column_is_feasible_at_node(col, m, node, inst_))
                feas_cols[m].push_back(col);
        }
    }

    // If no columns exist for any machine, add heuristic columns
    bool any_empty = false;
    for (int m : inst_.all_machines)
        if (feas_cols[m].empty()) { any_empty = true; break; }

    if (any_empty) {
        auto init_cols = build_initial_columns(inst_, horizon_);
        for (auto& [m, cols] : init_cols)
            for (auto& col : cols)
                if (column_is_feasible_at_node(col, m, node, inst_))
                    feas_cols[m].push_back(col);
    }

    // ── Repair pricing: cover any op not present in any surviving column ───
    // (more targeted than the any_empty heuristic above; matches Python logic)
    {
        // Build set of covered (j,o) pairs
        std::set<Op> covered;
        for (auto& [m, cols] : feas_cols)
            for (auto& col : cols)
                covered.insert(col.a.begin(), col.a.end());

        // For each op not covered, try to generate a repair column
        for (int j = 0; j < inst_.num_jobs; ++j) {
            for (int o = 0; o < inst_.num_ops[j]; ++o) {
                Op jop{j, o};
                if (covered.count(jop)) continue;

                bool found = false;
                for (int m : inst_.elig_mach[j][o]) {
                    // Respect OMA constraints at this node
                    auto req_it = node.required_machine.find(jop);
                    if (req_it != node.required_machine.end() && req_it->second != m) continue;
                    auto forb_it = node.forbidden_machines.find(jop);
                    if (forb_it != node.forbidden_machines.end() && forb_it->second.count(m)) continue;
                    auto allow_it = node.allowed_machines.find(jop);
                    if (allow_it != node.allowed_machines.end() && !allow_it->second.count(m)) continue;

                    auto repair_cols = repair_pricing_for_op(
                        m, jop, inst_, horizon_,
                        node.required_machine, node.forbidden_machines, node.allowed_machines,
                         node.start_lb, node.start_ub);
                    for (auto& col : repair_cols) {
                        if (column_is_feasible_at_node(col, m, node, inst_)) {
                            feas_cols[m].push_back(col);
                            covered.insert(jop);
                            found = true;
                            break;
                        }
                    }
                    if (found) break;
                }
                if (!found) {
                    // Last resort: build a trivial single-op column on any eligible
                    // machine (OMA constraints only, ignoring order constraints) so
                    // that CG has a warmstart column and can find an order-compatible
                    // solution during pricing.
                    for (int m : inst_.elig_mach[j][o]) {
                        auto req_it = node.required_machine.find(jop);
                        if (req_it != node.required_machine.end() && req_it->second != m) continue;
                        auto forb_it = node.forbidden_machines.find(jop);
                        if (forb_it != node.forbidden_machines.end() && forb_it->second.count(m)) continue;
                        auto allow_it = node.allowed_machines.find(jop);
                        if (allow_it != node.allowed_machines.end() && !allow_it->second.count(m)) continue;

                        int lb_start = 0;
                        auto it_lb = node.start_lb.find(jop);
                        if (it_lb != node.start_lb.end()) lb_start = it_lb->second;
                        int pt = inst_.proc_time[j][o].at(m);
                        Column col = make_column(m, {{j, o, lb_start, lb_start + pt}});
                        feas_cols[m].push_back(col);
                        covered.insert(jop);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    // No machine can provide any column covering this op → truly infeasible
                    node.status = NodeStatus::INFEASIBLE;
                    ++pruned_infeasible_;
                    return {};
                }
            }
        }
    }

    // Build RMP and run column generation
    RMP rmp(inst_, horizon_);
    for (auto& [m, cols] : feas_cols)
        rmp.add_columns(m, cols);

    CGResult cg = column_generation(
        inst_, horizon_, rmp, CG_NODE_MAX_ITER,
        &node.required_machine, &node.forbidden_machines, &node.allowed_machines, &node.start_lb, &node.start_ub,
        true);

    last_cg_iters_  = cg.iterations;
    last_cg_time_   = cg.cg_time;
    last_beam_frac_ = cg.beam_frac;
    last_lp_cmax_   = cg.cmax_value;
    node.lp_cmax    = cg.cmax_value;  // set early so PRUNED lines can report it

    // CG failure (RMP infeasible or numerical failure) => treat as infeasible
    if (!cg.solved_optimally || !std::isfinite(cg.lp_bound)) {
        node.lower_bound = node.cp_lower_bound;
        node.status = NodeStatus::INFEASIBLE;
        ++pruned_infeasible_;
        return {};
    }

    // Check artificial variables before using lp_bound: when arts > 0 the
    // objective is inflated by BIG_M so it is not a valid LP bound.
    bool has_art = false;
    for (int m : inst_.all_machines) {
        if (rmp.get_artificial(m) > LP_POS_TOL) { has_art = true; break; }
    }
    if (!has_art) {
        for (int j = 0; j < inst_.num_jobs && !has_art; ++j)
            for (int o = 0; o < inst_.num_ops[j] && !has_art; ++o)
                if (rmp.get_assign_artificial({j, o}) > LP_POS_TOL)
                    has_art = true;
    }
    if (has_art) {
        node.lower_bound = node.cp_lower_bound;
        node.status = NodeStatus::INFEASIBLE;
        ++pruned_infeasible_;
        return {};
    }

    double lb = std::max(cg.lp_bound, node.cp_lower_bound);
    node.lower_bound = lb;

    // Pruning: Cmax is integer, so the best reachable integer solution from this
    // node is ceil(lp_bound).  If ceil(lp_bound) >= best_ub_ we cannot improve.
    // Equivalent continuous condition: lp_bound > best_ub_ - 1 + INT_TOL
    if (lb > best_ub_ - 1.0 + INT_TOL) {
        node.status = NodeStatus::PRUNED;
        ++pruned_bound_;
        return {};
    }

    // Save column state for children
    node.cols = cg.columns;
    node.y_values = cg.y_values;

    // Check integrality
    if (is_integral(node, cg.y_values, cg.columns)) {
        node.status = NodeStatus::INTEGER;
        ++integer_nodes_;
        auto sched = extract_schedule(node, cg.y_values, cg.columns);
        update_ub(cg.cmax_value, sched);
        return {};
    }

    // Branch
    auto branch_dec = select_branch(node, cg.y_values, cg.columns);
    if (!branch_dec) {
        // Should not happen if not integral; treat as pruned
        node.status = NodeStatus::PRUNED;
        ++pruned_bound_;
        return {};
    }

    // Create child nodes
    std::vector<BranchAndPriceNode> children(2);
    for (int ci = 0; ci < 2; ++ci) {
        children[ci]          = node;  // copy all constraints
        children[ci].node_id  = next_node_id_++;
        children[ci].parent_id= node.node_id;
        children[ci].depth    = node.depth + 1;
        children[ci].status   = NodeStatus::OPEN;
        children[ci].lower_bound = lb;
        children[ci].branch_decision = *branch_dec;
    }

    if (branch_dec->type == BranchDecision::Type::OMA) {
        Op jop   = branch_dec->oma_op;
        int fix_m = branch_dec->oma_machine;

        // Child 0: fix op to machine fix_m
        children[0].required_machine[jop] = fix_m;
        children[0].allowed_machines[jop] = {fix_m};
        for (int m : inst_.all_machines) {
            if (m != fix_m) children[0].forbidden_machines[jop].insert(m);
        }

        // Child 1: forbid machine fix_m
        children[1].forbidden_machines[jop].insert(fix_m);
        auto& dom = children[1].allowed_machines[jop];
        if (dom.empty()) {
            for (int m : inst_.elig_mach[jop.first][jop.second]) dom.insert(m);
        }
        dom.erase(fix_m);
        if (dom.empty()) {
            children[1].status = NodeStatus::INFEASIBLE;
        }
    } 

    // Propagate
    std::vector<BranchAndPriceNode> result;
    for (auto& child : children) {
        if (child.status == NodeStatus::INFEASIBLE) continue;
        // Pre-filter inherited column pool: OMA/time-bound only (no order constraints)
        for (int m : inst_.all_machines) {
            std::vector<Column> kept;
            for (auto& col : child.cols[m])
                if (column_is_feasible_at_node(col, m, child, inst_))
                    kept.push_back(col);
            child.cols[m] = std::move(kept);
        }
        if (!propagate_node(child, inst_, best_ub_)) continue;
        if (child.cp_lower_bound >= best_ub_ - BOUND_EPS) continue;
        result.push_back(std::move(child));
    }

    return result;
}

BranchAndPriceResult BranchAndPriceSolver::solve()
{
    auto t_start = std::chrono::steady_clock::now();

    // Root node
    BranchAndPriceNode root;
    root.node_id    = next_node_id_++;
    root.depth      = 0;
    root.lower_bound= 0.0;

    // Build initial columns for root
    root.cols = build_initial_columns(inst_, horizon_);

    // Run CP-SAT for initial upper bound
    std::cout << "[B&P] Computing initial upper bound via CP-SAT...\n";
    auto cp_res = cp_sat_solve(inst_, root, CP_SAT_INITIAL_TL);
    if (cp_res.found) {
        update_ub(cp_res.makespan, cp_res.schedule);
        std::cout << "[B&P] Initial UB = " << best_ub_ << "\n";
    }

    // Best-first priority queue (min-heap by lower_bound)
    using QElem = std::pair<double, BranchAndPriceNode>;
    auto cmp = [](const QElem& a, const QElem& b){ return a.first > b.first; };
    std::priority_queue<QElem, std::vector<QElem>, decltype(cmp)> pq(cmp);
    pq.push({root.lower_bound, root});

    int nodes_explored = 0;
    double best_lb = 0.0;  // track best lower bound for reporting

    while (!pq.empty() && nodes_explored < MAX_NODES) {
        auto [lb, node] = pq.top();
        pq.pop();

        if (lb >= best_ub_ - BOUND_EPS) { ++pruned_bound_; continue; }

        ++nodes_explored;
        if (nodes_explored % 10 == 0) {
            std::ostringstream line;
            line << "[B&P] Node " << nodes_explored
                 << " | queue_LB=" << std::fixed << std::setprecision(4) << lb
                 << " | UB=" << std::fixed << std::setprecision(4) << best_ub_
                 << " | Queue=" << pq.size()
                 << "\n";
            std::cout << line.str();
        }

        // Periodic CP-SAT heuristic
        if (nodes_explored % CP_SAT_NODE_FREQ == 0) {
            auto cpres = cp_sat_solve(inst_, node, CP_SAT_NODE_TL, best_ub_);
            if (cpres.found) update_ub(cpres.makespan, cpres.schedule);
        }

        auto children = process_node(node);

        // Print a line for every terminal (pruned/infeasible/integer) node
        if (node.status != NodeStatus::OPEN) {
            const char* reason =
                node.status == NodeStatus::INFEASIBLE ? "infeasible"  :
                node.status == NodeStatus::INTEGER    ? "integrality" : "bound";
            std::ostringstream pline;
            pline << "[B&P] Node " << node.node_id
                  << " | LB=" << std::fixed << std::setprecision(4) << node.lower_bound
                  << " | LP_Cmax=" << std::fixed << std::setprecision(4) << node.lp_cmax
                  << " | UB=" << std::fixed << std::setprecision(4) << best_ub_
                  << " | PRUNED (" << reason << ")\n";
            std::cout << pline.str();
        }

        for (auto& child : children) {
            best_lb = std::max(best_lb, child.lower_bound);
            pq.push({child.lower_bound, child});
        }

        if (on_node_cb)
            on_node_cb(node.node_id, node.lower_bound, best_ub_);
    }

    auto t_end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    bool proven_optimal = (pq.empty() && nodes_explored < MAX_NODES);

    BranchAndPriceResult res;
    res.best_bound      = proven_optimal ? best_ub_ : best_lb;
    res.best_obj        = best_ub_;
    res.best_schedule   = best_schedule_;
    res.proven_optimal  = proven_optimal;
    res.nodes_explored    = nodes_explored;
    res.integer_nodes     = integer_nodes_;
    res.pruned_bound      = pruned_bound_;
    res.pruned_infeasible = pruned_infeasible_;
    res.elapsed_seconds   = elapsed;

    return res;
}
