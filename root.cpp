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
#include <sstream>
#include <stdexcept>

// #include <fstream>
// #include "mlp_policy.h"

namespace sat = operations_research::sat;

static constexpr double DINF = std::numeric_limits<double>::infinity();
static constexpr int    IINF = std::numeric_limits<int>::max() / 2;

// Column pool size caps — prune after each CG to control warm-start pool growth.
static constexpr int MAX_COLS_PER_MACHINE_ROOT = 250;
static constexpr int MAX_COLS_PER_MACHINE_NODE = 150;

static inline double to_integer_lb(double lb_cont) {
    if (!std::isfinite(lb_cont)) return lb_cont;
    return std::ceil(lb_cont - 1e-6);  
}

struct SequenceCandidate {
    bool valid = false;
    int machine = -1;
    Op A{-1, -1};
    Op B{-1, -1};
    double ab_mass = 0.0;
    double ba_mass = 0.0;
    double score = -1.0;
};

// ===========================================================================
//  LP-rounding hint for CP-SAT
// ===========================================================================

/// Assign each operation to the eligible machine having the largest
/// CG/LP assignment mass x[j,o,m]. Then build an earliest-start schedule.
static Schedule compute_lp_hint(
    const Instance& inst,
    const BranchAndPriceNode& node,
    const std::map<std::tuple<int, int, int>, double>& x_vals)
{
    // Step 1: choose one machine for every operation.
    std::map<Op, int> assigned;

    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            const Op jop{j, o};

            // Respect a machine fixed by branching.
            auto it_required = node.required_machine.find(jop);

            if (it_required != node.required_machine.end()) {
                assigned[jop] = it_required->second;
                continue;
            }

            int best_machine = -1;
            double best_weight =
                -std::numeric_limits<double>::infinity();

            for (int m : inst.elig_mach[j][o]) {
                auto it_forbidden =
                    node.forbidden_machines.find(jop);

                if (it_forbidden !=
                        node.forbidden_machines.end() &&
                    it_forbidden->second.count(m)) {
                    continue;
                }

                auto it_allowed =
                    node.allowed_machines.find(jop);

                if (it_allowed !=
                        node.allowed_machines.end() &&
                    !it_allowed->second.count(m)) {
                    continue;
                }

                double lp_weight = 0.0;

                auto it_x = x_vals.find({j, o, m});
                if (it_x != x_vals.end()) {
                    lp_weight = it_x->second;
                }

                // Prefer the highest LP mass.
                // Break ties using the shorter processing time.
                if (best_machine < 0 ||
                    lp_weight > best_weight + 1e-12 ||
                    (std::abs(lp_weight - best_weight) <= 1e-12 &&
                     inst.proc_time[j][o].at(m) <
                     inst.proc_time[j][o].at(best_machine)))
                {
                    best_weight = lp_weight;
                    best_machine = m;
                }
            }

            if (best_machine >= 0) {
                assigned[jop] = best_machine;
            }
        }
    }

    // Step 2: construct an earliest-start schedule for the chosen routing.
    std::map<int, int> machine_free;
    std::map<Op, int> operation_end;

    Schedule hint;

    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            const Op jop{j, o};

            auto it_machine = assigned.find(jop);
            if (it_machine == assigned.end()) {
                continue;
            }

            const int machine = it_machine->second;
            const int processing_time =
                inst.proc_time[j][o].at(machine);

            int earliest_start = 0;

            // Job precedence.
            if (o > 0) {
                auto it_previous =
                    operation_end.find({j, o - 1});

                if (it_previous != operation_end.end()) {
                    earliest_start =
                        std::max(
                            earliest_start,
                            it_previous->second
                        );
                }
            }

            // Current B&P node lower bound.
            auto it_lb = node.start_lb.find(jop);

            if (it_lb != node.start_lb.end()) {
                earliest_start =
                    std::max(
                        earliest_start,
                        it_lb->second
                    );
            }

            // Machine availability.
            auto it_machine_free =
                machine_free.find(machine);

            if (it_machine_free != machine_free.end()) {
                earliest_start =
                    std::max(
                        earliest_start,
                        it_machine_free->second
                    );
            }

            const int end =
                earliest_start + processing_time;

            operation_end[jop] = end;
            machine_free[machine] = end;

            hint[jop] = {
                machine,
                earliest_start,
                end
            };
        }
    }

    return hint;
}



// ===========================================================================
//  CP-SAT solver
// ===========================================================================

// CPSatResult cp_sat_solve(const Instance& inst,
//                           const BranchAndPriceNode& node,
//                           double time_limit_sec,
//                           double ub,
//                           const Schedule* hint_schedule)
static CPSatResult cp_sat_solve_impl(
    const Instance& inst,
    const BranchAndPriceNode& node,
    double time_limit_sec,
    double ub,
    const Schedule* hint_schedule,
    bool feasibility_only,
    int target_cmax)
{
    CPSatResult result;
    int J = inst.num_jobs;

    sat::CpModelBuilder model;
    int horizon = 0;
    for (int j = 0; j < J; ++j){
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            for (auto& [m, pt] : inst.proc_time[j][o])
                horizon += pt;
        }
    }

    // In feasibility mode, all operations must finish by target_cmax.
    // In normal optimization mode, retain the original large horizon.
    const int H = feasibility_only
        ? target_cmax
        : horizon;

    if (H < 0) {
        return result;
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
            if (eligs.empty()) {
                return result;
            }

            // Minimum processing time among the machines still allowed
            // at this B&P node.
            int min_pt = IINF;

            for (int m : eligs) {
                min_pt = std::min(
                    min_pt,
                    inst.proc_time[j][o].at(m)
                );
            }

            if (min_pt == IINF) {
                return result;
            }

            // Earliest allowed start.
            int slb = 0;

            auto it_lb = node.start_lb.find(jop);
            if (it_lb != node.start_lb.end()) {
                slb = std::max(slb, it_lb->second);
            }

            // The operation must have enough room to finish by H.
            int sub = H - min_pt;

            // Apply the B&P node start upper bound.
            auto it_ub = node.start_ub.find(jop);
            if (it_ub != node.start_ub.end()) {
                sub = std::min(sub, it_ub->second);
            }

            // No valid start time exists.
            if (slb > sub) {
                return result;
            }

            // Start must lie inside [slb, sub].
            // auto sv = model.NewIntVar(
            //     sat::Domain(
            //         static_cast<int64_t>(slb),
            //         static_cast<int64_t>(sub)
            //     )
            // );
            auto sv = model.NewIntVar({
                static_cast<int64_t>(slb),
                static_cast<int64_t>(sub)
            });

            // End must lie inside [slb + min_pt, H].
            // auto ev = model.NewIntVar(
            //     sat::Domain(
            //         static_cast<int64_t>(slb + min_pt),
            //         static_cast<int64_t>(H)
            //     )
            // );
            auto ev = model.NewIntVar({
                static_cast<int64_t>(slb + min_pt),
                static_cast<int64_t>(H)
            });

            start_vars[{j, o}] = sv;
            end_vars[{j, o}]   = ev;

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

    // Sequence branching constraints.
    // For each machine m and pair (A,B):
    // if A and B are both assigned to m, enforce end(A) <= start(B).
    for (const auto& [machine, arcs] : node.order_constraints) {
        for (const auto& [A, B] : arcs) {
            auto it_pa = presence.find({A.first, A.second, machine});
            auto it_pb = presence.find({B.first, B.second, machine});

            if (it_pa == presence.end() || it_pb == presence.end()) {
                continue;
            }

            auto ct = model.AddLessOrEqual(
                end_vars[A],
                start_vars[B]
            );

            ct.OnlyEnforceIf({it_pa->second, it_pb->second});
        }
    }

    sat::IntVar cmax = model.NewIntVar({
        int64_t{0},
        static_cast<int64_t>(H)
    });

    for (int j = 0; j < J; ++j) {
        int last_o = inst.num_ops[j] - 1;
        model.AddLessOrEqual(end_vars[{j, last_o}], cmax);
    }

    // Prune makespan above current UB
    // if (std::isfinite(ub))
    //     model.AddLessThan(cmax, model.NewConstant((int64_t)ub));

    // Normal optimization mode:
    // ub = 37 means search for Cmax < 37, therefore Cmax <= 36.
    if (!feasibility_only && std::isfinite(ub)) {
        model.AddLessThan(
            cmax,
            model.NewConstant(static_cast<int64_t>(ub))
        );
    }

    // Pure feasibility mode:
    // search for any schedule satisfying Cmax <= target_cmax.
    if (feasibility_only) {
        model.AddLessOrEqual(
            cmax,
            model.NewConstant(static_cast<int64_t>(target_cmax))
        );
    }

        if (hint_schedule != nullptr) {
        for (int j = 0; j < J; ++j) {
            for (int o = 0; o < inst.num_ops[j]; ++o) {
                const Op jop{j, o};

                auto it_h = hint_schedule->find(jop);
                if (it_h == hint_schedule->end()) {
                    continue;
                }

                const auto [hm, hs, he] = it_h->second;

                // Skip if the hinted machine is not eligible at this node
                // (would otherwise zero-hint all eligible presence vars).
                const auto eligible_for_hint = get_eligible(j, o);
                if (std::find(eligible_for_hint.begin(),
                              eligible_for_hint.end(),
                              hm) == eligible_for_hint.end()) {
                    continue;
                }

                // Hint machine assignment.
                for (int m : eligible_for_hint) {
                    auto it_p = presence.find({j, o, m});

                    if (it_p == presence.end()) {
                        continue;
                    }

                    model.AddHint(
                        it_p->second,
                        static_cast<int64_t>(m == hm ? 1 : 0)
                    );
                }

                // Add timing hints when they fit inside the current target horizon.
                bool valid_timing_hint = false;
                {
                    auto pt_it = inst.proc_time[j][o].find(hm);
                    if (pt_it != inst.proc_time[j][o].end()) {
                        const int hinted_pt = pt_it->second;
                        int hinted_lb = 0;
                        int hinted_ub = H - hinted_pt;

                        auto it_lb = node.start_lb.find(jop);
                        if (it_lb != node.start_lb.end())
                            hinted_lb = std::max(hinted_lb, it_lb->second);

                        auto it_ub = node.start_ub.find(jop);
                        if (it_ub != node.start_ub.end())
                            hinted_ub = std::min(hinted_ub, it_ub->second);

                        valid_timing_hint =
                            hs >= hinted_lb &&
                            hs <= hinted_ub &&
                            he == hs + hinted_pt &&
                            he <= H;
                    }
                }

                if (!feasibility_only || valid_timing_hint) {
                    auto it_sv = start_vars.find({j, o});
                    auto it_ev = end_vars.find({j, o});

                    if (it_sv != start_vars.end()) {
                        model.AddHint(
                            it_sv->second,
                            static_cast<int64_t>(hs)
                        );
                    }

                    if (it_ev != end_vars.end()) {
                        model.AddHint(
                            it_ev->second,
                            static_cast<int64_t>(he)
                        );
                    }
                }
            } // closes for o
        }     // closes for j
    }         // closes if hint_schedule

    // Objective only in optimization mode.
    if (!feasibility_only) {
        model.Minimize(cmax);
    }

    sat::CpSolverResponse response;
    {
        sat::SatParameters params;
        if (time_limit_sec > 0.0)
            params.set_max_time_in_seconds(time_limit_sec);
        params.set_log_search_progress(false);
        response = sat::SolveWithParameters(model.Build(), params);
    }

    if (response.status() == sat::CpSolverStatus::OPTIMAL ||
    response.status() == sat::CpSolverStatus::FEASIBLE)
    {
        int actual_cmax = 0;

        for (int j = 0; j < J; ++j) {
            for (int o = 0; o < inst.num_ops[j]; ++o) {
                const int start = static_cast<int>(
                    sat::SolutionIntegerValue(
                        response,
                        start_vars.at({j, o})
                    )
                );

                bool machine_found = false;

                for (int m : get_eligible(j, o)) {
                    auto it = presence.find({j, o, m});

                    if (it == presence.end()) {
                        continue;
                    }

                    if (!sat::SolutionBooleanValue(response, it->second)) {
                        continue;
                    }

                    const int pt = inst.proc_time[j][o].at(m);
                    const int end = start + pt;

                    result.schedule[{j, o}] = {
                        m,
                        start,
                        end
                    };

                    actual_cmax = std::max(actual_cmax, end);
                    machine_found = true;
                    break;
                }

                if (!machine_found) {
                    result.schedule.clear();
                    result.makespan = -1.0;
                    result.found = false;
                    return result;
                }
            }
        }

        result.makespan = static_cast<double>(actual_cmax);
        result.found = true;
    }

    return result;
}

CPSatResult cp_sat_solve(
    const Instance& inst,
    const BranchAndPriceNode& node,
    double time_limit_sec,
    double ub,
    const Schedule* hint_schedule)
{
    return cp_sat_solve_impl(
        inst,
        node,
        time_limit_sec,
        ub,
        hint_schedule,
        /*feasibility_only=*/false,
        /*target_cmax=*/IINF
    );
}

CPSatResult cp_sat_feasibility_solve(
    const Instance& inst,
    const BranchAndPriceNode& node,
    double time_limit_sec,
    int target_cmax,
    const Schedule* hint_schedule)
{
    return cp_sat_solve_impl(
        inst,
        node,
        time_limit_sec,
        /*ub=*/std::numeric_limits<double>::infinity(),
        hint_schedule,
        /*feasibility_only=*/true,
        target_cmax
    );
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

                // Sync fully-reduced domain back to allowed_machines so it is always
                // canonical (forbidden already applied). Rebuild forbidden_machines from
                // elig_mach \ domain to keep the two maps consistent and avoid disagreement
                // between code that checks only one of them.
                node.allowed_machines[jop] = domain;
                node.forbidden_machines[jop].clear();
                for (int mm : inst.elig_mach[j][o])
                    if (!domain.count(mm))
                        node.forbidden_machines[jop].insert(mm);

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
                //
                // The earliest start of every operation must respect both:
                //   1. completion of the previous operation;
                //   2. the explicit start_lb already stored at the node.
                std::vector<int> es(n_ops, 0);

                if (n_ops > 0) {
                    const Op first_operation{j, 0};

                    auto first_lb_it =
                        node.start_lb.find(first_operation);

                    if (first_lb_it != node.start_lb.end()) {
                        es[0] = std::max(
                            0,
                            first_lb_it->second
                        );
                    }

                    // Detect an inconsistent time window immediately.
                    auto first_ub_it =
                        node.start_ub.find(first_operation);

                    if (first_ub_it != node.start_ub.end() &&
                        es[0] > first_ub_it->second)
                    {
                        return false;
                    }
                }

                for (int o = 1; o < n_ops; ++o) {
                    // Processing time of the previous operation under its
                    // current machine domain.
                    const auto previous_domain =
                        get_domain(j, o - 1);

                    if (previous_domain.empty()) {
                        return false;
                    }

                    int previous_min_processing_time = IINF;

                    for (int machine : previous_domain) {
                        previous_min_processing_time =
                            std::min(
                                previous_min_processing_time,
                                inst.proc_time[j][o - 1].at(machine)
                            );
                    }

                    if (previous_min_processing_time == IINF) {
                        return false;
                    }

                    // Lower bound imposed by job precedence.
                    const int precedence_lower_bound =
                        es[o - 1] +
                        previous_min_processing_time;

                    // Explicit lower bound already present at this B&P node.
                    const Op current_operation{j, o};

                    int explicit_lower_bound = 0;

                    auto current_lb_it =
                        node.start_lb.find(current_operation);

                    if (current_lb_it != node.start_lb.end()) {
                        explicit_lower_bound =
                            current_lb_it->second;
                    }

                    // The explicit bound is now propagated immediately to all
                    // following operations through es[o].
                    es[o] =
                        std::max(
                            precedence_lower_bound,
                            explicit_lower_bound
                        );

                    // Check consistency with an explicit upper start bound.
                    auto current_ub_it =
                        node.start_ub.find(current_operation);

                    if (current_ub_it != node.start_ub.end() &&
                        es[o] > current_ub_it->second)
                    {
                        return false;
                    }
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

                    // Update the node with the fully propagated earliest start.
                    auto current_lb_it =
                        node.start_lb.find(jop);

                    const int current_lower_bound =
                        current_lb_it != node.start_lb.end()
                            ? current_lb_it->second
                            : 0;

                    if (es[o] > current_lower_bound) {
                        node.start_lb[jop] = es[o];
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

    // required_machine: if an op is required on THIS machine, it must appear in the column.
    // A column on machine m that omits a required-on-m op is infeasible for that branch.
    for (auto& [req_op, req_m] : node.required_machine) {
        if (req_m == machine && !col.a.count(req_op)) return false;
    }

    auto oc_it = node.order_constraints.find(machine);

    if (oc_it != node.order_constraints.end()) {
        for (const auto& [A, B] : oc_it->second) {
            const bool has_A = col.a.count(A);
            const bool has_B = col.a.count(B);

            // Constraint është conditional:
            // aplikohet vetëm nëse A dhe B janë të dyja në kolonën e kësaj makine.
            if (!has_A || !has_B) {
                continue;
            }

            const int end_A = col.e_vals.at(A);
            const int start_B = col.s_vals.at(B);

            if (end_A > start_B) {
                return false;
            }
        }
    }

    return true;
}

// ===========================================================================
//  Static branching helpers
// ===========================================================================

/// Column signature: ordered list of (job, op_index) pairs in machine-scheduling order.
static std::vector<std::pair<int,int>> get_col_signature(const Column& col) {
    std::vector<std::pair<int,int>> sig;
    sig.reserve(col.ops.size());
    for (auto& [j, o, s, e] : col.ops)
        sig.push_back({j, o});
    return sig;
}

/// x[(j,o,m)] = sum_{p} a[j,o,m,p] * y[m,p]
static std::map<std::tuple<int,int,int>, double>
compute_assignment_masses(const std::map<std::pair<int,int>, double>& y_values,
                           const std::map<int, std::vector<Column>>& cols)
{
    std::map<std::tuple<int,int,int>, double> x;
    for (auto& [mp, y] : y_values) {
        if (y < LP_POS_TOL) continue;
        auto [m, p] = mp;
        auto it = cols.find(m);
        if (it == cols.end() || p >= (int)it->second.size()) continue;
        for (auto& [jop, s] : it->second[p].s_vals) {
            auto [j, o] = jop;
            x[{j, o, m}] += y;
        }
    }
    return x;
}

/// Returns true when every operation has LP mass ≈ 1 on exactly one machine.
/// Populates selected_machine_out.
static bool is_routing_fixed_lp(
    const BranchAndPriceNode& node,
    const std::map<std::tuple<int,int,int>, double>& x_vals,
    const Instance& inst,
    std::map<Op, int>& selected_machine_out)
{
    selected_machine_out.clear();
    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            Op jop{j, o};
            auto req_it = node.required_machine.find(jop);
            if (req_it != node.required_machine.end()) {
                selected_machine_out[jop] = req_it->second;
                continue;
            }
            int    best_m   = -1;
            double best_val = 0.0;
            for (int m : inst.elig_mach[j][o]) {
                double v = 0.0;
                auto it = x_vals.find({j, o, m});
                if (it != x_vals.end()) v = it->second;
                if (v > best_val) { best_val = v; best_m = m; }
            }
            if (best_val >= 1.0 - 10 * LP_POS_TOL && best_m >= 0)
                selected_machine_out[jop] = best_m;
            else
                return false;
        }
    }
    return true;
}

/// Extract an integer schedule from fixed routing + machine orders via longest-path.
/// Returns {schedule, makespan} or {{}, DINF} on failure (cycle or missing data).
static std::pair<Schedule, double>
extract_integer_schedule_from_routing(
    const std::map<Op, int>& selected_machine,
    // const std::map<int, std::vector<Op>>& machine_orders,
    const Instance& inst)
{
    std::vector<Op> all_ops;
    for (int j = 0; j < inst.num_jobs; ++j)
        for (int o = 0; o < inst.num_ops[j]; ++o)
            all_ops.push_back({j, o});

    std::map<Op, std::vector<std::pair<Op, int>>> adj;
    std::map<Op, int> in_deg;
    for (auto& jop : all_ops) { adj[jop] = {}; in_deg[jop] = 0; }

    // Job-precedence arcs: (j,o) -> (j,o+1) with weight = pt[j,o,selected_m]
    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0; o + 1 < inst.num_ops[j]; ++o) {
            Op from{j, o}, to{j, o + 1};
            auto it = selected_machine.find(from);
            if (it == selected_machine.end()) return {{}, DINF};
            int m = it->second;
            auto pt_it = inst.proc_time[j][o].find(m);
            if (pt_it == inst.proc_time[j][o].end()) return {{}, DINF};
            adj[from].push_back({to, pt_it->second});
            in_deg[to]++;
        }
    }

    // Longest-path via Kahn's topological sort
    std::map<Op, int> dist;
    for (auto& jop : all_ops) dist[jop] = 0;

    std::queue<Op> tq;
    for (auto& jop : all_ops)
        if (in_deg[jop] == 0) tq.push(jop);

    std::vector<Op> topo;
    while (!tq.empty()) {
        Op u = tq.front(); tq.pop();
        topo.push_back(u);
        for (auto& [v, wt] : adj[u]) {
            dist[v] = std::max(dist[v], dist[u] + wt);
            if (--in_deg[v] == 0) tq.push(v);
        }
    }

    if ((int)topo.size() != (int)all_ops.size())
        return {{}, DINF};  // cycle detected

    Schedule sched;
    double cmax = 0.0;
    for (auto& jop : all_ops) {
        auto it = selected_machine.find(jop);
        if (it == selected_machine.end()) return {{}, DINF};
        int m  = it->second;
        int s  = dist[jop];
        int pt = inst.proc_time[jop.first][jop.second].at(m);
        sched[jop] = {m, s, s + pt};
        cmax = std::max(cmax, (double)(s + pt));
    }
    return {sched, cmax};
}

static bool validate_schedule_result(
    const Schedule& sched,
    const Instance& inst)
{
    // ---------------------------------------------------------
    // 1. The schedule must contain every operation exactly once.
    // ---------------------------------------------------------
    int total_operations = 0;

    for (int j = 0; j < inst.num_jobs; ++j) {
        total_operations += inst.num_ops[j];
    }

    if (static_cast<int>(sched.size()) != total_operations) {
        std::cerr
            << "[Schedule validation] Wrong number of operations:"
            << " expected=" << total_operations
            << " found=" << sched.size()
            << "\n";

        return false;
    }

    // Machine intervals:
    // machine -> (start, end, operation)
    std::map<
        int,
        std::vector<
            std::tuple<int, int, Op>
        >
    > machine_intervals;

    // ---------------------------------------------------------
    // 2. Check every operation individually.
    // ---------------------------------------------------------
    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            const Op operation{j, o};

            auto schedule_it =
                sched.find(operation);

            if (schedule_it == sched.end()) {
                std::cerr
                    << "[Schedule validation] Missing operation:"
                    << " job=" << j
                    << " op=" << o
                    << "\n";

                return false;
            }

            const auto [machine, start, end] =
                schedule_it->second;

            // Start and end times must be valid.
            if (start < 0 || end <= start) {
                std::cerr
                    << "[Schedule validation] Invalid time interval:"
                    << " job=" << j
                    << " op=" << o
                    << " machine=" << machine
                    << " start=" << start
                    << " end=" << end
                    << "\n";

                return false;
            }

            // Machine must be eligible.
            const auto& eligible_machines =
                inst.elig_mach[j][o];

            if (std::find(
                    eligible_machines.begin(),
                    eligible_machines.end(),
                    machine)
                == eligible_machines.end())
            {
                std::cerr
                    << "[Schedule validation] Ineligible machine:"
                    << " job=" << j
                    << " op=" << o
                    << " machine=" << machine
                    << "\n";

                return false;
            }

            // Processing time must exist for this machine.
            auto processing_time_it =
                inst.proc_time[j][o].find(machine);

            if (processing_time_it ==
                inst.proc_time[j][o].end())
            {
                std::cerr
                    << "[Schedule validation] Missing processing time:"
                    << " job=" << j
                    << " op=" << o
                    << " machine=" << machine
                    << "\n";

                return false;
            }

            const int expected_processing_time =
                processing_time_it->second;

            if (end != start + expected_processing_time) {
                std::cerr
                    << "[Schedule validation] Wrong processing time:"
                    << " job=" << j
                    << " op=" << o
                    << " machine=" << machine
                    << " start=" << start
                    << " end=" << end
                    << " expected_end="
                    << start + expected_processing_time
                    << "\n";

                return false;
            }

            machine_intervals[machine].push_back({
                start,
                end,
                operation
            });
        }
    }

    // ---------------------------------------------------------
    // 3. Check job precedence:
    //
    // end(j,o) <= start(j,o+1)
    // ---------------------------------------------------------
    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0;
             o + 1 < inst.num_ops[j];
             ++o)
        {
            const auto& current =
                sched.at({j, o});

            const auto& next =
                sched.at({j, o + 1});

            const int current_end =
                std::get<2>(current);

            const int next_start =
                std::get<1>(next);

            if (current_end > next_start) {
                std::cerr
                    << "[Schedule validation] Job precedence violation:"
                    << " job=" << j
                    << " op=" << o
                    << " ends=" << current_end
                    << " but op=" << o + 1
                    << " starts=" << next_start
                    << "\n";

                return false;
            }
        }
    }

    // ---------------------------------------------------------
    // 4. Check machine no-overlap.
    // ---------------------------------------------------------
    for (auto& [machine, intervals] :
         machine_intervals)
    {
        std::sort(
            intervals.begin(),
            intervals.end(),
            [](const auto& first,
               const auto& second)
            {
                if (std::get<0>(first) !=
                    std::get<0>(second))
                {
                    return std::get<0>(first) <
                           std::get<0>(second);
                }

                return std::get<1>(first) <
                       std::get<1>(second);
            }
        );

        for (std::size_t index = 1;
             index < intervals.size();
             ++index)
        {
            const auto& previous =
                intervals[index - 1];

            const auto& current =
                intervals[index];

            const int previous_end =
                std::get<1>(previous);

            const int current_start =
                std::get<0>(current);

            if (current_start < previous_end) {
                const Op previous_operation =
                    std::get<2>(previous);

                const Op current_operation =
                    std::get<2>(current);

                std::cerr
                    << "[Schedule validation] Machine overlap:"
                    << " machine=" << machine
                    << " previous=(j="
                    << previous_operation.first
                    << ",o="
                    << previous_operation.second
                    << ",end="
                    << previous_end
                    << ") current=(j="
                    << current_operation.first
                    << ",o="
                    << current_operation.second
                    << ",start="
                    << current_start
                    << ")\n";

                return false;
            }
        }
    }

    return true;
}

static double compute_schedule_makespan(
    const Schedule& schedule)
{
    double makespan = 0.0;

    for (const auto& [operation, assignment] :
         schedule)
    {
        const int end =
            std::get<2>(assignment);

        makespan =
            std::max(
                makespan,
                static_cast<double>(end)
            );
    }

    return makespan;
}

/// Validate a column: no duplicate ops, correct processing times, machine compatible.
static bool validate_column_impl(const Column& col, int machine, const Instance& inst)
{
    std::set<Op> seen;
    int prev_e = -1;
    for (auto& [j, o, s, e] : col.ops) {
        Op jop{j, o};
        if (seen.count(jop)) return false;
        seen.insert(jop);
        auto& elig = inst.elig_mach[j][o];
        if (!std::count(elig.begin(), elig.end(), machine)) return false;
        auto pt_it = inst.proc_time[j][o].find(machine);
        if (pt_it == inst.proc_time[j][o].end()) return false;
        if (e != s + pt_it->second || s < 0) return false;
        if (s < prev_e) return false;  // overlap on machine
        prev_e = e;
    }
    return true;
}

// ===========================================================================
//  B&P solver implementation
// ===========================================================================

BranchAndPriceSolver::BranchAndPriceSolver(const Instance& inst)
    : inst_(inst), horizon_(compute_horizon(inst))
{}

void BranchAndPriceSolver::set_initial_ub(double makespan, const Schedule& sched)
{
    update_ub(makespan, sched);
}

bool BranchAndPriceSolver::update_ub(
    double reported_makespan,
    const Schedule& schedule)
{
    // Never accept an invalid schedule as an incumbent.
    if (!validate_schedule_result(
            schedule,
            inst_))
    {
        std::cerr
            << "[B&P] Rejected invalid incumbent schedule."
            << " Reported makespan="
            << reported_makespan
            << "\n";

        return false;

        // During debugging, instead of return false, you can use:
        //
        // throw std::runtime_error(
        //     "Invalid schedule passed to update_ub()."
        // );
    }

    // Compute the makespan directly from the validated schedule.
    const double actual_makespan =
        compute_schedule_makespan(
            schedule
        );

    if (std::abs(
            actual_makespan -
            reported_makespan)
        > BOUND_EPS)
    {
        std::cerr
            << "[B&P] Makespan mismatch:"
            << " reported="
            << reported_makespan
            << " actual="
            << actual_makespan
            << ". Using actual value.\n";
    }

    if (actual_makespan <
        best_ub_ - BOUND_EPS)
    {

        best_ub_ =
            actual_makespan;

        best_schedule_ =
            schedule;

        // Store this accepted incumbent schedule as a source of future columns.
        incumbent_schedules_.push_back(schedule);

        // Keep only the last few incumbents to avoid growing the pool too much.
        static constexpr std::size_t MAX_INCUMBENT_SCHEDULES = 5;

        if (incumbent_schedules_.size() > MAX_INCUMBENT_SCHEDULES) {
            incumbent_schedules_.erase(
                incumbent_schedules_.begin()
            );
        }

        horizon_ =
            static_cast<int>(
                std::ceil(best_ub_)
            );

        std::cout
            << "[B&P] New incumbent: UB="
            << actual_makespan
            << "  horizon="
            << horizon_
            << " incumbent_pool="
            << incumbent_schedules_.size()
            << "\n";

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


struct OMACandidate {
    bool valid = false;
    Op op{-1, -1};
    int machine = -1;
    double x = 0.0;
    double score = -1.0;
};

static std::map<std::tuple<int, int, int>, double>
compute_oma_masses(
    const std::map<int, std::vector<Column>>& cols,
    const std::map<std::pair<int, int>, double>& y_values)
{
    std::map<std::tuple<int, int, int>, double> x;

    for (const auto& [machine, machine_cols] : cols) {
        for (int p = 0; p < static_cast<int>(machine_cols.size()); ++p) {
            auto it_y = y_values.find({machine, p});
            if (it_y == y_values.end()) {
                continue;
            }

            const double y = it_y->second;
            if (y <= LP_POS_TOL) {
                continue;
            }

            const Column& col = machine_cols[p];

            for (const Op& op : col.a) {
                x[{op.first, op.second, machine}] += y;
            }
        }
    }

    return x;
}


OMACandidate find_best_oma_candidate(
    const Instance& inst,
    const BranchAndPriceNode& node,
    const std::map<int, std::vector<Column>>& cols,
    const std::map<std::pair<int, int>, double>& y_values)
{
    OMACandidate best;

    const auto x_vals =
        compute_oma_masses(cols, y_values);

    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            Op op{j, o};

            // If already fixed to one machine, no OMA branching needed.
            if (node.required_machine.count(op)) {
                continue;
            }

            for (int machine : inst.elig_mach[j][o]) {
                // Skip machines already forbidden.
                auto it_forbidden =
                    node.forbidden_machines.find(op);

                if (it_forbidden != node.forbidden_machines.end() &&
                    it_forbidden->second.count(machine))
                {
                    continue;
                }

                // If allowed_machines exists, respect it.
                auto it_allowed =
                    node.allowed_machines.find(op);

                if (it_allowed != node.allowed_machines.end() &&
                    !it_allowed->second.count(machine))
                {
                    continue;
                }

                double x = 0.0;

                auto it_x =
                    x_vals.find({j, o, machine});

                if (it_x != x_vals.end()) {
                    x = it_x->second;
                }

                // We branch only on fractional OMA mass.
                if (x <= LP_POS_TOL || x >= 1.0 - LP_POS_TOL) {
                    continue;
                }

                // Balanced branching is best when x is close to 0.5.
                // Score in [0,1], maximum at x=0.5.
                const double balance_score =
                    4.0 * x * (1.0 - x);

                // Prefer operations with more routing flexibility.
                // const double flexibility_bonus =
                //     1.0 + 0.05 *
                //     static_cast<double>(
                //         std::max(
                //             0,
                //             static_cast<int>(inst.elig_mach[j][o].size()) - 1
                //         )
                //     );

                // const double score =
                //     balance_score * flexibility_bonus;

                const double score = balance_score;

                if (score > best.score) {
                    best.valid = true;
                    best.op = op;
                    best.machine = machine;
                    best.x = x;
                    best.score = score;
                }
            }
        }
    }

    return best;
}
 

static SequenceCandidate find_best_sequence_candidate(
    const BranchAndPriceNode& node,
    const std::map<int, std::vector<Column>>& cols,
    const std::map<std::pair<int, int>, double>& y_values)
{
    SequenceCandidate best;

    for (const auto& [machine, machine_cols] : cols) {
        std::map<std::pair<Op, Op>, double> before_mass;

        for (int p = 0; p < static_cast<int>(machine_cols.size()); ++p) {
            auto it_y = y_values.find({machine, p});
            if (it_y == y_values.end()) {
                continue;
            }

            const double y = it_y->second;
            if (y <= LP_POS_TOL) {
                continue;
            }

            const Column& col = machine_cols[p];

            std::vector<Op> ordered_ops;
            ordered_ops.reserve(col.ops.size());

            for (const auto& entry : col.ops) {
                const int j = std::get<0>(entry);
                const int o = std::get<1>(entry);
                ordered_ops.push_back({j, o});
            }

            for (int i = 0; i < static_cast<int>(ordered_ops.size()); ++i) {
                for (int k = i + 1; k < static_cast<int>(ordered_ops.size()); ++k) {
                    Op A = ordered_ops[i];
                    Op B = ordered_ops[k];

                    // Mos branch mbi operacione të të njëjtit job,
                    // sepse job precedence e përcakton tashmë rendin.
                    if (A.first == B.first) {
                        continue;
                    }

                    before_mass[{A, B}] += y;
                }
            }
        }

        for (const auto& [pair_ab, ab_mass] : before_mass) {
            Op A = pair_ab.first;
            Op B = pair_ab.second;

            auto it_ba = before_mass.find({B, A});
            if (it_ba == before_mass.end()) {
                continue;
            }

            const double ba_mass = it_ba->second;

            if (ab_mass <= LP_POS_TOL || ba_mass <= LP_POS_TOL) {
                continue;
            }

            // Mos branch përsëri mbi të njëjtin çift në këtë makinë.
            auto it_oc = node.order_constraints.find(machine);
            if (it_oc != node.order_constraints.end()) {
                if (it_oc->second.count({A, B}) ||
                    it_oc->second.count({B, A}))
                {
                    continue;
                }
            }

            const double total = ab_mass + ba_mass;
            if (total <= LP_POS_TOL) {
                continue;
            }

            const double balance =
                1.0 -
                std::abs(ab_mass - ba_mass) /
                std::max(1e-9, total);

            const double score =
                std::min(1.0, total) * balance;

            if (score > best.score) {
                best.valid = true;
                best.machine = machine;
                best.A = A;
                best.B = B;
                best.ab_mass = ab_mass;
                best.ba_mass = ba_mass;
                best.score = score;
            }
        }
    }

    return best;
}

std::optional<BranchDecision>
BranchAndPriceSolver::select_branch(
    const BranchAndPriceNode& node,
    const std::map<std::pair<int,int>, double>& y_values,
    const std::map<int, std::vector<Column>>& cols)
{

    // const OMACandidate oma =
    //     find_best_oma_candidate(
    //         inst_,
    //         node,
    //         node.cols,
    //         node.y_values
    //     );

    // const SequenceCandidate seq =
    //     find_best_sequence_candidate(
    //         node,
    //         node.cols,
    //         node.y_values
    //     );

    const OMACandidate oma =
        find_best_oma_candidate(
            inst_,
            node,
            cols,
            y_values
        );

    const SequenceCandidate seq =
        find_best_sequence_candidate(
            node,
            cols,
            y_values
        );

    // static constexpr double OMA_WEIGHT = 1.15;
    // static constexpr double SEQ_WEIGHT = 0.9;

    static constexpr double OMA_WEIGHT = 1.15;
    static constexpr double SEQ_WEIGHT = 0.9;

    double oma_score =
        oma.valid ? OMA_WEIGHT * oma.score : -1.0;

    double seq_score =
        seq.valid ? SEQ_WEIGHT * seq.score : -1.0;


    const double same_machine_mass =
        seq.valid
            ? seq.ab_mass + seq.ba_mass
            : 0.0;

    const bool strong_sequence_candidate =
        seq.valid &&
        same_machine_mass >= 0.75 &&
        seq.score >= 0.60;

    const bool choose_sequence =
        strong_sequence_candidate &&
        (
            !oma.valid ||
            (
                node.depth >= 2 &&
                seq_score >= 0.80 * oma_score
            )
        );

    // Optional dynamic boost:
    // If no operation is fixed yet, OMA is usually more useful.

    // if (node.required_machine.empty() && oma.valid) {
    //     oma_score *= 1.20;
    // }

    // If routing is already significantly restricted, sequence becomes stronger.
    if (!node.required_machine.empty() && seq.valid) {
        seq_score *= 1.05;
    }

    if (oma_score >= seq_score && oma.valid) {
        std::cout
            << "  Branch OMA: J" << oma.op.first
            << "o" << oma.op.second
            << " on M" << oma.machine
            << " x=" << oma.x
            << " raw_score=" << oma.score
            << " weighted_score=" << oma_score
            << " seq_best=" << seq_score
            << "\n";

        return BranchDecision::oma(
            oma.op,
            oma.machine
        );
    }

    if (seq.valid) {
        std::cout
            << "  Branch SEQUENCE: machine=" << seq.machine
            << " left: J" << seq.A.first << "o" << seq.A.second
            << " before J" << seq.B.first << "o" << seq.B.second
            << " | right: J" << seq.B.first << "o" << seq.B.second
            << " before J" << seq.A.first << "o" << seq.A.second
            << " ab_mass=" << seq.ab_mass
            << " ba_mass=" << seq.ba_mass
            << " raw_score=" << seq.score
            << " weighted_score=" << seq_score
            << " oma_best=" << oma_score
            << "\n";

        return BranchDecision::sequence(
            seq.machine,
            seq.A,
            seq.B
        );
    }

    return std::nullopt;
}

/// Prune a per-machine column pool to at most max_per_machine columns.
/// Always retains LP-active columns (y > 1e-6).
/// Fills the remaining budget from the tail of each machine's list
/// (most recently generated columns, which had the best reduced costs).
static void prune_column_pool(
    std::map<int, std::vector<Column>>& cols,
    const std::map<std::pair<int,int>, double>& y_values,
    int max_per_machine)
{
    for (auto& [m, mcols] : cols) {
        int n = (int)mcols.size();
        if (n <= max_per_machine) continue;

        // Mark LP-active columns (y > 1e-6).
        std::vector<bool> active(n, false);
        int n_active = 0;
        for (int ci = 0; ci < n; ++ci) {
            auto it = y_values.find({m, ci});
            if (it != y_values.end() && it->second > 1e-6) {
                active[ci] = true;
                ++n_active;
            }
        }

        int budget = std::max(0, max_per_machine - n_active);

        // Keep all active, then fill budget from the end (most recent).
        std::vector<Column> kept;
        kept.reserve(n_active + budget);
        for (int ci = 0; ci < n; ++ci)
            if (active[ci]) kept.push_back(mcols[ci]);
        int taken = 0;
        for (int ci = n - 1; ci >= 0 && taken < budget; --ci)
            if (!active[ci]) { kept.push_back(mcols[ci]); ++taken; }
        mcols = std::move(kept);
    }
}

std::vector<BranchAndPriceNode>
BranchAndPriceSolver::process_node(BranchAndPriceNode& node)
{
    const double remaining_at_node_start =
        BAP_GLOBAL_TIME_LIMIT -
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solve_start_
        ).count();
    if (remaining_at_node_start <= 0.0) {
        std::cout << "  [B&P] global time limit reached before node processing\n";
        return {};
    }

    // ── Step 1: Filter inherited column pool ──────────────────────────────────
    // Enforce OMA, time-bounds, order constraints, and banned/fixed-column rules.
    int cols_before = 0, removed_oma = 0, removed_banned = 0;
    std::map<int, std::vector<Column>> feas_cols;
    for (int m : inst_.all_machines) {
        feas_cols[m] = {};
        for (auto& col : node.cols[m]) {
            ++cols_before;
            if (!column_is_feasible_at_node(col, m, node, inst_)) {
                ++removed_oma; continue;
            }
            // if (!column_is_feasible_at_node(col, m, node, inst_, /*check_order=*/true)) {
            //     ++removed_order; continue;
            // }
            auto sig = get_col_signature(col);
            if (node.banned_column_signatures.count(sig)) {
                ++removed_banned; continue;
            }
            auto it_fix = node.fixed_machine_columns.find(m);
            if (it_fix != node.fixed_machine_columns.end() && sig != it_fix->second) {
                ++removed_banned; continue;
            }
            feas_cols[m].push_back(col);
        }
    }
    int cols_after = 0;
    for (int m : inst_.all_machines) cols_after += (int)feas_cols[m].size();

    // If filtering leaves a small pool, inject additional branching-feasible
    // initial columns before running CG.
    static constexpr int MIN_COLS_AFTER_FILTER = 0;
    if (cols_after < MIN_COLS_AFTER_FILTER) {
        auto init_cols = build_initial_columns(inst_, horizon_);
        int added = 0;
        for (auto& [m, cols] : init_cols) {
            for (auto& col : cols) {
                if (!column_is_feasible_at_node(col, m, node, inst_))
                    continue;
                // if (!column_is_feasible_at_node(col, m, node, inst_))
                //     continue;
                auto sig = get_col_signature(col);
                if (node.banned_column_signatures.count(sig))
                    continue;
                auto it_fix = node.fixed_machine_columns.find(m);
                if (it_fix != node.fixed_machine_columns.end() && sig != it_fix->second)
                    continue;
                feas_cols[m].push_back(col);
                ++added;
            }
        }
        cols_after = 0;
        for (int m : inst_.all_machines) cols_after += (int)feas_cols[m].size();
        std::cout << "  low-column fallback: added=" << added
                  << " -> after_filter=" << cols_after << "\n";
    }

    // ── Per-node debug header ─────────────────────────────────────────────────
    {
        int n_req  = (int)node.required_machine.size();
        int n_forb = 0;
        for (auto& [op, ms] : node.forbidden_machines) n_forb += (int)ms.size();

        // int n_ord  = 0;
        // // for (auto& [m, oc] : node.order_constraints) n_ord += (int)oc.size();

        int n_ord = 0;
        for (auto& [m, oc] : node.order_constraints) {
            n_ord += static_cast<int>(oc.size());
        }

        std::cout << "[Node " << node.node_id << " depth=" << node.depth << "]\n"
                  << "  UB=" << std::fixed << std::setprecision(4) << best_ub_
                  << "  T=" << horizon_ << "\n"
                //   << "  required_OMA=" << n_req
                //   << " forbidden_OMA=" << n_forb

                // << "  required_machine=" << n_req
                // << " forbidden_machine=" << n_forb
                // << " order_constraints=" << n_ord

                //   << " order_constraints=" << n_ord << "\n"
                  << "  columns before_filter=" << cols_before
                  << " after_filter=" << cols_after;
                //   << " (removed OMA=" << removed_oma
                //   << " order=" << removed_order
                //   << " banned/fixed=" << removed_banned << ")\n";
    }

    // ── Step 2: Ensure every machine has at least one column and every op is covered ──

    // If any machine has an empty pool after filtering, refill from heuristics
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

    // Local helper: check whether a column contains every mandatory operation.
    auto column_contains_all =
        [](const Column& column,
        const std::set<Op>& mandatory_operations) -> bool
    {
        for (const Op& operation : mandatory_operations) {
            if (!column.a.count(operation)) {
                return false;
            }
        }

        return true;
    };

    // Local helper: avoid inserting an exact duplicate column.
    auto column_already_exists =
        [&](int machine,
            const Column& candidate) -> bool
    {
        auto machine_it =
            feas_cols.find(machine);

        if (machine_it == feas_cols.end()) {
            return false;
        }

        for (const Column& existing :
            machine_it->second)
        {
            if (existing.ops == candidate.ops) {
                return true;
            }
        }

        return false;
    };

    // Build the mandatory-operation set for every machine.
    std::map<int, std::set<Op>>
        mandatory_by_machine;

    for (const auto& [operation, machine] :
        node.required_machine)
    {
        mandatory_by_machine[machine].insert(
            operation
        );
    }

    // Repair every machine that has one or more mandatory operations.
    for (const auto& [machine, mandatory_operations] :
        mandatory_by_machine)
    {
        if (mandatory_operations.empty()) {
            continue;
        }

        // First check whether the inherited/initial pool already contains
        // a valid column with the complete mandatory set.
        bool joint_column_exists = false;

        auto machine_columns_it =
            feas_cols.find(machine);

        if (machine_columns_it !=
            feas_cols.end())
        {
            for (const Column& column :
                machine_columns_it->second)
            {
                if (!column_contains_all(
                        column,
                        mandatory_operations))
                {
                    continue;
                }

                if (!column_is_feasible_at_node(
                        column,
                        machine,
                        node,
                        inst_))
                {
                    continue;
                }

                joint_column_exists = true;
                break;
            }
        }

        if (joint_column_exists) {
            continue;
        }

        std::cout
            << "  [Mandatory repair] machine="
            << machine
            << " required_ops="
            << mandatory_operations.size()
            << "\n";

        // Generate columns containing ALL mandatory operations simultaneously.
        auto repair_columns =
            repair_pricing_for_machine(
                machine,
                mandatory_operations,
                inst_,
                horizon_,
                node.required_machine,
                node.forbidden_machines,
                node.allowed_machines,
                node.start_lb,
                node.start_ub,
                /*top_k=*/10
            );

        int added_columns = 0;

        for (Column& column : repair_columns) {
            // Defensive verification: the returned column must contain
            // the full mandatory set.
            if (!column_contains_all(
                    column,
                    mandatory_operations))
            {
                continue;
            }

            // Verify all branching restrictions and time bounds.
            if (!column_is_feasible_at_node(
                    column,
                    machine,
                    node,
                    inst_))
            {
                continue;
            }

            // Verify that the complete column lies inside the current horizon.
            bool inside_horizon = true;

            for (const auto& [operation, end] :
                column.e_vals)
            {
                if (end > horizon_) {
                    inside_horizon = false;
                    break;
                }
            }

            if (!inside_horizon) {
                continue;
            }

            if (column_already_exists(
                    machine,
                    column))
            {
                continue;
            }

            feas_cols[machine].push_back(
                std::move(column)
            );

            ++added_columns;
        }

        if (added_columns > 0) {
            std::cout
                << "  [Mandatory repair] added="
                << added_columns
                << " joint columns for machine="
                << machine
                << "\n";
        } else {
            // Do not prune here. The exact Phase-I pricing must decide
            // whether the node is truly infeasible.
            std::cout
                << "  [Mandatory repair] no joint column generated"
                << " for machine=" << machine
                << "; leaving recovery to exact Phase-I pricing.\n";
        }
    }


    {
        std::set<Op> covered;

        for (const auto& [machine, columns] :
            feas_cols)
        {
            for (const Column& column : columns) {
                covered.insert(
                    column.a.begin(),
                    column.a.end()
                );
            }
        }

        for (int j = 0;
            j < inst_.num_jobs;
            ++j)
        {
            for (int o = 0;
                o < inst_.num_ops[j];
                ++o)
            {
                const Op operation{j, o};

                if (covered.count(operation)) {
                    continue;
                }

                bool found = false;

                // If the operation is fixed by branching, try only its
                // required machine. Otherwise try all eligible machines.
                std::vector<int> candidate_machines;

                auto required_it =
                    node.required_machine.find(
                        operation
                    );

                if (required_it !=
                    node.required_machine.end())
                {
                    candidate_machines.push_back(
                        required_it->second
                    );
                } else {
                    candidate_machines =
                        inst_.elig_mach[j][o];
                }

                for (int machine :
                    candidate_machines)
                {
                    auto forbidden_it =
                        node.forbidden_machines.find(
                            operation
                        );

                    if (forbidden_it !=
                            node.forbidden_machines.end() &&
                        forbidden_it->second.count(machine))
                    {
                        continue;
                    }

                    auto allowed_it =
                        node.allowed_machines.find(
                            operation
                        );

                    if (allowed_it !=
                            node.allowed_machines.end() &&
                        !allowed_it->second.count(machine))
                    {
                        continue;
                    }

                    // auto repair_columns =
                    //     repair_pricing_for_op(
                    //         machine,
                    //         operation,
                    //         inst_,
                    //         horizon_,
                    //         node.required_machine,
                    //         node.forbidden_machines,
                    //         node.allowed_machines,
                    //         node.start_lb,
                    //         node.start_ub
                    //     );

                    auto repair_columns =
                        repair_pricing_for_op(
                            machine,
                            operation,
                            inst_,
                            horizon_,
                            node.required_machine,
                            node.forbidden_machines,
                            node.allowed_machines,
                            node.order_constraints,
                            node.start_lb,
                            node.start_ub
                        );

                    for (Column& column :
                        repair_columns)
                    {
                        // The target operation must appear.
                        if (!column.a.count(operation)) {
                            continue;
                        }

                        // This also checks that if the machine has mandatory
                        // operations, the column contains all of them.
                        if (!column_is_feasible_at_node(
                                column,
                                machine,
                                node,
                                inst_))
                        {
                            continue;
                        }

                        bool inside_horizon = true;

                        for (const auto& [column_op, end] :
                            column.e_vals)
                        {
                            if (end > horizon_) {
                                inside_horizon = false;
                                break;
                            }
                        }

                        if (!inside_horizon) {
                            continue;
                        }

                        if (!column_already_exists(
                                machine,
                                column))
                        {
                            // Record every operation covered by the new column,
                            // not only the target operation.
                            for (const Op& covered_operation :
                                column.a)
                            {
                                covered.insert(
                                    covered_operation
                                );
                            }

                            feas_cols[machine].push_back(
                                std::move(column)
                            );
                        } else {
                            // A duplicate valid column still means the target
                            // operation is already covered in the pool.
                            covered.insert(operation);
                        }

                        found = true;
                        break;
                    }

                    if (found) {
                        break;
                    }
                }

                if (!found) {
                    // Do not build an invalid single-op fallback and do not
                    // prune the node here. Assignment artificial variables
                    // will keep the RMP solvable, and exact Phase-I pricing
                    // must determine true infeasibility.
                    std::cout
                        << "  [Coverage repair] no initial column found"
                        << " for op(j=" << j
                        << ",o=" << o
                        << "); leaving it to exact Phase-I pricing.\n";
                }
            }
        }
    }

    {
        std::vector<Schedule> schedules_to_inject =
            incumbent_schedules_;

        if (schedules_to_inject.empty() &&
            !best_schedule_.empty())
        {
            schedules_to_inject.push_back(
                best_schedule_
            );
        }

        if (!schedules_to_inject.empty()) {
            auto inc_cols =
                build_initial_columns(
                    inst_,
                    horizon_,
                    schedules_to_inject
                );

            int injected_cols = 0;

            for (auto& [m, mcols] : inc_cols) {
                for (auto& col : mcols) {
                    if (column_is_feasible_at_node(
                            col,
                            m,
                            node,
                            inst_))
                    {
                        // Avoid exact duplicates already present in feas_cols.
                        bool dup = false;

                        for (auto& existing : feas_cols[m]) {
                            if (existing.ops == col.ops) {
                                dup = true;
                                break;
                            }
                        }

                        if (!dup) {
                            feas_cols[m].push_back(col);
                            ++injected_cols;
                        }
                    }
                }
            }

        }
    }

    // ── Step 3: Run Column Generation ─────────────────────────────────────────
    RMP rmp(inst_, horizon_);
    for (auto& [m, cols] : feas_cols)
        rmp.add_columns(m, cols);

    const double cg_time_limit_base =
        (node.depth == 0) ? CG_TIME_LIMIT_ROOT : CG_TIME_LIMIT_NODE;
    const double cg_time_limit = std::min(
        cg_time_limit_base,
        std::max(0.0, BAP_GLOBAL_TIME_LIMIT -
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - solve_start_
            ).count())
    );

    const bool near_ub = std::isfinite(best_ub_) &&
        (node.lp_cmax > best_ub_ - 1.0 - 0.5);  // brenda 0.5 nga prune threshold
    const int effective_stall =
        (node.depth >= 5 && near_ub) ? 12 : CG_STALL_LIMIT;  // 12 vs 6

    CGResult cg = column_generation(
        inst_, horizon_, rmp, CG_NODE_MAX_ITER,
        &node.required_machine,
        &node.forbidden_machines,
        &node.allowed_machines,
        &node.order_constraints,
        &node.start_lb,
        &node.start_ub,
        /*use_beam_first=*/true,
        cg_time_limit,
        best_ub_,
        effective_stall);

    last_cg_iters_  = cg.iterations;
    last_cg_time_   = cg.cg_time;
    last_beam_frac_ = cg.beam_frac;
    last_lp_cmax_   = cg.cmax_value;
    node.lp_cmax    = cg.cmax_value;

    // Accumulate global column and pricing-split counters.
    int node_cols = 0;
    for (auto& [m, cols] : cg.columns)
        node_cols += (int)cols.size();
    total_cols_generated_ += node_cols;
    total_pricing_time_ += cg.pricing_time;
    total_master_time_  += cg.master_time;
    sum_beam_frac_        += cg.beam_frac;
    ++cg_node_count_;

    // ── Step 4: CG and Phase-I handling ───────────────────────────────────────

    // Exact infeasibility certificate from Phase-I pricing.
    if (cg.phase1_infeasible_certified) {
        // node.lower_bound = node.cp_lower_bound;

        node.lower_bound = to_integer_lb(node.cp_lower_bound);

        node.status = NodeStatus::INFEASIBLE;
        ++pruned_infeasible_;

        std::cout
            << "  PRUNED: exact Phase-I pricing certified infeasibility"
            << " artificial_sum=" << cg.artificial_sum
            << "\n";

        return {};
    }

    if (!cg.solved_optimally) {
        std::cout
            << "  CG stopped before convergence after "
            << cg.iterations
            << " iterations"
            << " time=" << cg.cg_time
            << "s\n";
    }

    // ── Step 5: Check artificial variables ────────────────────────────────────
    auto check_art = [&]() -> bool {
        for (int m : inst_.all_machines) {
            if (rmp.get_artificial(m) > LP_POS_TOL) {
                return true;
            }
        }

        for (int j = 0; j < inst_.num_jobs; ++j) {
            for (int o = 0; o < inst_.num_ops[j]; ++o) {
                if (rmp.get_assign_artificial({j, o}) > LP_POS_TOL) {
                    return true;
                }
            }
        }

        return false;
    };

    const bool has_art = check_art();

    // A positive artificial without a timeout should not be used
    // for branching or bound pruning either.
    if (has_art) {
        node.lower_bound = node.cp_lower_bound;
        node.status = NodeStatus::INFEASIBLE;
        ++pruned_infeasible_;

        heuristic_cg_prune_occurred_ = true;

        std::cout
            << "  PRUNED HEURISTICALLY: positive artificial variables remain"
            << " artificial_sum=" << cg.artificial_sum
            << "\n";

        return {};
    }

    double lb = node.cp_lower_bound;

    if (cg.solved_optimally && cg.phase1_feasible && !has_art && std::isfinite(cg.lp_bound)) {
        lb = std::max(node.cp_lower_bound, cg.lp_bound);
    }
    node.lower_bound = to_integer_lb(lb);

    {
        // const double gap_pct = (std::isfinite(best_ub_) && best_ub_ > 1e-9)
        //     ? (best_ub_ - lb) / best_ub_ * 100.0
        //     : 100.0;

        const double gap_pct = (std::isfinite(best_ub_) && best_ub_ > 1e-9)
            ? (best_ub_ - node.lower_bound) / best_ub_ * 100.0
            : 100.0;
            
        std::cout
            << "  LB=" << std::fixed << std::setprecision(4) << lb
            << "  LB int=" << std::fixed << std::setprecision(4) << node.lower_bound
            << "  UB=" << std::setprecision(4) << best_ub_
            << "  Gap=" << std::setprecision(2) << gap_pct << "%"
            << "  Cols=" << node_cols
            << "  P/M=" << std::setprecision(1) << cg.beam_frac * 100.0 << "%"
            << "  CG_iters=" << cg.iterations
            << "\n";
    }

    // ── Step 5d: one CP-SAT call after CG ──────────────────────────────────
    int next_stall_count = node.consecutive_ub_hits;

    // Convert the CG solution into operation-machine LP masses.
    auto x_vals_step5 = compute_assignment_masses(
        cg.y_values,
        cg.columns
    );

    // Build a pure CG/LP-rounding hint.
    // This function should now have only three arguments after RL removal.
    Schedule lp_hint = compute_lp_hint(
        inst_,
        node,
        x_vals_step5
    );

    // Compute the makespan of the CG-derived hint for diagnostics.
    int hint_cmax = horizon_;

    BranchAndPriceNode cg_node = node;
    {
        // Helper: original node-constrained domain for op (j,o).
        auto get_node_domain = [&](int j, int o) -> std::set<int> {
            const Op jop{j, o};
            auto it_r = node.required_machine.find(jop);
            if (it_r != node.required_machine.end()) return {it_r->second};
            std::set<int> dom(inst_.elig_mach[j][o].begin(),
                              inst_.elig_mach[j][o].end());
            auto it_a = node.allowed_machines.find(jop);
            if (it_a != node.allowed_machines.end()) {
                std::set<int> inter;
                for (int m : it_a->second) if (dom.count(m)) inter.insert(m);
                dom = inter;
            }
            auto it_f = node.forbidden_machines.find(jop);
            if (it_f != node.forbidden_machines.end())
                for (int m : it_f->second) dom.erase(m);
            return dom;
        };

        // ── Phase 1: collect CG-supported machines per op ────────────────────
        std::map<Op, std::set<int>> cg_op_machines;
        for (auto& [mp, y] : cg.y_values) {
            if (y <= LP_POS_TOL) continue;
            auto [m, p] = mp;
            auto it_c = cg.columns.find(m);
            if (it_c == cg.columns.end() || p >= (int)it_c->second.size()) continue;
            for (auto& [jop, _s] : it_c->second[p].s_vals)
                cg_op_machines[jop].insert(m);
        }

        // ── CG confidence per operation-machine assignment ─────────────────────
        // x[j,o,m] = LP mass of assigning operation (j,o) to machine m.
        // If one machine has high LP mass, we can fix that operation for CP-SAT.
        std::map<Op, int>    cg_best_machine;
        std::map<Op, double> cg_best_mass;
        std::map<Op, double> cg_total_mass;

        for (const auto& [jom, x] : x_vals_step5) {
            auto [j, o, m] = jom;
            const Op jop{j, o};

            if (x <= LP_POS_TOL) {
                continue;
            }

            cg_total_mass[jop] += x;

            auto it_best = cg_best_mass.find(jop);
            if (it_best == cg_best_mass.end() || x > it_best->second) {
                cg_best_mass[jop] = x;
                cg_best_machine[jop] = m;
            }
        }

        // Tuning parameter.
        // Higher = safer, fixes fewer ops.
        // Lower  = more aggressive, fixes more ops.
        static constexpr double CG_FIX_CONFIDENCE = 0.85;

        for (int j = 0; j < inst_.num_jobs; ++j) {
            for (int o = 0; o < inst_.num_ops[j]; ++o) {
                const Op jop{j, o};
                if (cg_node.required_machine.count(jop)) continue;

                auto it_cg = cg_op_machines.find(jop);
                if (it_cg == cg_op_machines.end()) continue;

                std::set<int> n_dom = get_node_domain(j, o);
                std::set<int> new_dom;

                // Keep machines supported by positive CG columns.
                for (int m : n_dom)
                    if (it_cg->second.count(m)) new_dom.insert(m);

                // Also keep the machine used by the current incumbent,
                // provided it remains feasible at this B&P node.
                auto inc_it = best_schedule_.find(jop);
                if (inc_it != best_schedule_.end()) {
                    const int inc_mach = std::get<0>(inc_it->second);
                    if (n_dom.count(inc_mach))
                        new_dom.insert(inc_mach);
                }

                // Safe fallback: never create an empty domain.
                if (new_dom.empty()) {
                    new_dom = n_dom;
                }

                // ── Stronger CG fixing ────────────────────────────────────────────────
                // If CG puts a large majority of the assignment mass on one machine,
                // force CP-SAT to use that machine for this operation.
                bool fixed_by_cg_confidence = false;

                auto it_best_m = cg_best_machine.find(jop);
                auto it_best_x = cg_best_mass.find(jop);
                auto it_tot_x  = cg_total_mass.find(jop);

                if (it_best_m != cg_best_machine.end() &&
                    it_best_x != cg_best_mass.end() &&
                    it_tot_x  != cg_total_mass.end())
                {
                    const int best_m = it_best_m->second;

                    const double confidence =
                        it_best_x->second / std::max(1e-9, it_tot_x->second);

                    // Only fix if the best machine is still feasible at this B&P node.
                    if (confidence >= CG_FIX_CONFIDENCE && n_dom.count(best_m)) {
                        new_dom.clear();
                        new_dom.insert(best_m);
                        fixed_by_cg_confidence = true;
                    }
                }

                cg_node.allowed_machines[jop] = new_dom;

                cg_node.forbidden_machines[jop].clear();
                for (int m : inst_.elig_mach[j][o]) {
                    if (!new_dom.count(m)) {
                        cg_node.forbidden_machines[jop].insert(m);
                    }
                }

                if (new_dom.size() == 1) {
                    cg_node.required_machine[jop] = *new_dom.begin();
                }
            }
        }

        // ── Phase 2: score ops for selective relaxation ──────────────────────
        struct OpCandidate {
            Op     jop;
            double score;
            int    extra_machine;  // fastest excluded machine
        };
        std::vector<OpCandidate> candidates;

        for (int j = 0; j < inst_.num_jobs; ++j) {
            for (int o = 0; o < inst_.num_ops[j]; ++o) {
                const Op jop{j, o};
                if (cg_node.required_machine.count(jop)) continue;

                std::set<int> n_dom = get_node_domain(j, o);
                auto it_cg = cg_op_machines.find(jop);
                if (it_cg == cg_op_machines.end()) continue;

                // CG-restricted domain (already intersected in Phase 1).
                std::set<int> cg_dom;
                for (int m : n_dom) if (it_cg->second.count(m)) cg_dom.insert(m);
                if (cg_dom.empty()) continue;

                // Excluded: in node domain but not in CG domain.
                std::vector<int> excluded;
                for (int m : n_dom)
                    if (!cg_dom.count(m)) excluded.push_back(m);
                if (excluded.empty()) continue;

                // ── Criticality ───────────────────────────────────────────────
                int min_pt_cg = IINF;
                for (int m : cg_dom) {
                    auto pt_it = inst_.proc_time[j][o].find(m);
                    if (pt_it != inst_.proc_time[j][o].end())
                        min_pt_cg = std::min(min_pt_cg, pt_it->second);
                }
                if (min_pt_cg == IINF) continue;

                int es = 0;
                {
                    auto it = node.start_lb.find(jop);
                    if (it != node.start_lb.end()) es = it->second;
                }
                double slack       = (double)(horizon_ - (es + min_pt_cg));
                double criticality = std::max(0.0, std::min(1.0,
                    1.0 - slack / std::max(1, horizon_)));

                // ── LP uncertainty ─────────────────────────────────────────────
                double max_x = 0.0;
                for (int m : cg_dom) {
                    auto it_x = x_vals_step5.find({j, o, m});
                    if (it_x != x_vals_step5.end())
                        max_x = std::max(max_x, it_x->second);
                }
                double lp_uncertainty = std::max(0.0, 1.0 - max_x);

                // ── Processing-time gain ───────────────────────────────────────
                // Fastest excluded machine vs. fastest CG machine.
                int min_pt_excl = IINF;
                int best_extra  = -1;
                for (int m : excluded) {
                    auto pt_it = inst_.proc_time[j][o].find(m);
                    if (pt_it != inst_.proc_time[j][o].end()
                            && pt_it->second < min_pt_excl) {
                        min_pt_excl = pt_it->second;
                        best_extra  = m;
                    }
                }
                double pt_gain = (min_pt_excl < min_pt_cg && min_pt_cg > 0)
                    ? (double)(min_pt_cg - min_pt_excl) / min_pt_cg
                    : 0.0;

                // ── Domain loss ────────────────────────────────────────────────
                double domain_loss =
                    1.0 - (double)cg_dom.size() / std::max(1, (int)n_dom.size());

                // ── Composite score ────────────────────────────────────────────
                double score = 0.50 * criticality
                             + 0.25 * lp_uncertainty
                             + 0.20 * pt_gain
                             + 0.05 * domain_loss;

                if (best_extra < 0) best_extra = excluded.front();
                candidates.push_back({jop, score, best_extra});
            }
        }

        // ── Phase 3: relax top-K operations ──────────────────────────────────
        int total_ops = 0;
        for (int j = 0; j < inst_.num_jobs; ++j) total_ops += inst_.num_ops[j];
        const int relax_k = std::max(2, total_ops / 5);

        std::sort(candidates.begin(), candidates.end(),
                  [](const OpCandidate& a, const OpCandidate& b){
                      return a.score > b.score;
                  });

        int relaxed = 0;
        for (auto& cand : candidates) {
            if (relaxed >= relax_k) break;
            const Op& jop = cand.jop;

            auto it_allow = cg_node.allowed_machines.find(jop);
            if (it_allow == cg_node.allowed_machines.end()) continue;

            it_allow->second.insert(cand.extra_machine);

            auto it_forb = cg_node.forbidden_machines.find(jop);
            if (it_forb != cg_node.forbidden_machines.end())
                it_forb->second.erase(cand.extra_machine);

            // Undo singleton fix if domain now has 2+ machines.
            if (it_allow->second.size() > 1) {
                auto it_req = cg_node.required_machine.find(jop);
                if (it_req != cg_node.required_machine.end())
                    cg_node.required_machine.erase(it_req);
            }

            ++relaxed;
        }

        // if (relaxed > 0)
        //     std::cout << "  [CP-SAT LNR] relaxed " << relaxed
        //               << "/" << (int)candidates.size()
        //               << " ops (K=" << relax_k << ")\n";
    }

    // Build the CP-SAT warm-start hint: incumbent schedule (filtered to ops
    // that are still machine-feasible at cg_node) merged with the CG/LP-rounding
    // hint for any operations not covered by the incumbent.
    Schedule cp_hint;
    for (const auto& [jop, mse] : best_schedule_) {
        const int machine = std::get<0>(mse);
        bool machine_allowed = true;

        auto it_req = cg_node.required_machine.find(jop);
        if (it_req != cg_node.required_machine.end() &&
            it_req->second != machine)
            machine_allowed = false;

        auto it_forb = cg_node.forbidden_machines.find(jop);
        if (it_forb != cg_node.forbidden_machines.end() &&
            it_forb->second.count(machine))
            machine_allowed = false;

        auto it_allow = cg_node.allowed_machines.find(jop);
        if (it_allow != cg_node.allowed_machines.end() &&
            !it_allow->second.count(machine))
            machine_allowed = false;

        if (machine_allowed)
            cp_hint[jop] = mse;
    }
    // Fill missing operations from the CG/LP-rounding hint.
    for (const auto& [jop, mse] : lp_hint)
        if (!cp_hint.count(jop))
            cp_hint[jop] = mse;

    // Execute CP-SAT at most once after CG, unless the node-level CP-SAT
    // ablation flag has disabled this improvement call.  The initial UB call
    // is independent and remains enabled.
    if (node_cp_sat_enabled_ && std::isfinite(best_ub_)) {
        const int incumbent_target =
            static_cast<int>(std::ceil(best_ub_)) - 1;

        // Respect both the current B&P horizon and the improving UB target.
        const int target_cmax = std::min(
            horizon_,
            incumbent_target
        );

        if (target_cmax >= 0) {
            const double elapsed_before_cp =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - solve_start_
                ).count();
            const double remaining_global_time =
                BAP_GLOBAL_TIME_LIMIT - elapsed_before_cp;
            const double node_cp_time = std::min(
                CP_SAT_NODE_TL,
                std::max(0.0, remaining_global_time)
            );

            if (node_cp_time <= 0.0) {
                std::cout
                    << "  [CP-SAT after CG] skipped: global time limit reached\n";
            } else {

            std::cout
                << "  [CP-SAT after CG] searching for Cmax <= "
                << target_cmax
                << " TL=" << node_cp_time
                << "s, T=" << horizon_
                << "\n";

            const auto cp_start =
                std::chrono::steady_clock::now();

            CPSatResult sat = cp_sat_feasibility_solve(
                inst_,
                cg_node,
                node_cp_time,
                target_cmax,
                cp_hint.empty() ? nullptr : &cp_hint
            );

            const double cp_elapsed =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - cp_start
                ).count();

            if (sat.found) {
                const bool improved = update_ub(
                    sat.makespan,
                    sat.schedule
                );

                std::cout
                    << "  [CP-SAT after CG] found Cmax="
                    << sat.makespan
                    << (improved ? "  ** new UB **" : "")
                    << " elapsed="
                    << std::fixed
                    << std::setprecision(1)
                    << cp_elapsed
                    << "s\n";

                next_stall_count = 0;
            } else {
                std::cout
                    << "  [CP-SAT after CG] no solution found for Cmax <= "
                    << target_cmax
                    << " elapsed="
                    << std::fixed
                    << std::setprecision(1)
                    << cp_elapsed
                    << "s\n";

                next_stall_count =
                    node.consecutive_ub_hits + 1;
            }
            }
        } else {
            next_stall_count = 0;
        }
    } else if (node_cp_sat_enabled_) {
        next_stall_count = 0;
    }

    // No-solution streak is tracked for diagnostics only under current pruning policy.

    // ── Step 6: Pruning by LP Cmax ─────────────────────────────────────────────
    // Only prune by LP_Cmax when CG is fully converged and no artificial variables remain.
    // If CG timed out or artificials are positive, the LP value is not reliable enough
    // to close the node.
    bool reliable_lp_for_pruning = cg.solved_optimally && !has_art;

    if (reliable_lp_for_pruning &&
        cg.cmax_value > best_ub_ - 1.0 + INT_TOL) {
        node.status = NodeStatus::PRUNED;
        ++pruned_bound_;
        std::cout << "  PRUNED (LP_Cmax > UB-1)\n";
        return {};
    }

    if (!reliable_lp_for_pruning) {
        std::cout << "  Skipping LP_Cmax prune (CG not converged or artificials positive)\n";
    }

    // Save column state for children (pruned to cap to control pool growth)
    node.cols     = cg.columns;
    node.y_values = cg.y_values;
    {
        int cap = (node.depth == 0) ? MAX_COLS_PER_MACHINE_ROOT : MAX_COLS_PER_MACHINE_NODE;
        prune_column_pool(node.cols, cg.y_values, cap);
    }


    if (is_integral(
        node,
        cg.y_values,
        cg.columns))
    {
        node.status =
            NodeStatus::INTEGER;

        ++integer_nodes_;

        Schedule schedule =
            extract_schedule(
                node,
                cg.y_values,
                cg.columns
            );

        const double actual_makespan =
            compute_schedule_makespan(
                schedule
            );

        const bool accepted =
            update_ub(
                actual_makespan,
                schedule
            );

        if (accepted) {
            std::cout
                << "  INTEGER (y binary):"
                << " validated Cmax="
                << actual_makespan
                << "\n";
        } else {
            std::cout
                << "  INTEGER (y binary):"
                << " schedule rejected or did not improve UB"
                << " reported_LP_Cmax="
                << cg.cmax_value
                << "\n";
        }

        return {};
    }

    // ── Step 8: Compute LP assignment masses x[j,o,m] ────────────────────────
    auto x_vals = compute_assignment_masses(cg.y_values, cg.columns);

    // Count fractional OMA candidates for debug
    int n_frac_oma = 0;
    for (int j = 0; j < inst_.num_jobs; ++j)
        for (int o = 0; o < inst_.num_ops[j]; ++o) {
            Op jop{j, o};
            if (node.required_machine.count(jop)) continue;
            for (int m : inst_.elig_mach[j][o]) {
                double v = 0.0;
                auto it = x_vals.find({j, o, m});
                if (it != x_vals.end()) v = it->second;
                if (v > LP_POS_TOL && v < 1.0 - LP_POS_TOL) { ++n_frac_oma; break; }
            }
        }
    // std::cout << "  fractional_OMA_candidates=" << n_frac_oma << "\n";
    std::cout << "  sequence_branching_enabled\n";

    // ── Steps 9-10: Branch selection (OMA first, ORDER second) ────────────────
    auto branch_dec =
        select_branch(
            node,
            cg.y_values,
            cg.columns
        );

    if (branch_dec) {
        std::vector<BranchAndPriceNode> children(2);

        for (int ci = 0; ci < 2; ++ci) {
            children[ci] = node;
            children[ci].node_id = next_node_id_++;
            children[ci].parent_id = node.node_id;
            children[ci].depth = node.depth + 1;
            children[ci].status = NodeStatus::OPEN;
            children[ci].lower_bound = lb;
            children[ci].branch_decision = *branch_dec;
            children[ci].consecutive_ub_hits = next_stall_count;
        }

        BranchAndPriceNode& left = children[0];
        BranchAndPriceNode& right = children[1];

        if (branch_dec->type == BranchDecision::Type::OMA) {
            Op op = branch_dec->oma_op;
            int machine = branch_dec->oma_machine;

            // Left child: fix operation on this machine.
            left.required_machine[op] = machine;
            left.allowed_machines[op].clear();
            left.allowed_machines[op].insert(machine);
            left.forbidden_machines[op].clear();

            // Right child: forbid this machine.
            right.forbidden_machines[op].insert(machine);

            // If right child has only one feasible machine left,
            // automatically fix the operation to that machine.
            std::set<int> remaining;

            for (int m : inst_.elig_mach[op.first][op.second]) {
                auto it_forb =
                    right.forbidden_machines.find(op);

                const bool forbidden =
                    it_forb != right.forbidden_machines.end() &&
                    it_forb->second.count(m);

                if (!forbidden) {
                    remaining.insert(m);
                }
            }

            if (remaining.empty()) {
                right.status = NodeStatus::INFEASIBLE;
            }
            else if (remaining.size() == 1) {
                const int only_machine = *remaining.begin();

                right.required_machine[op] = only_machine;
                right.allowed_machines[op].clear();
                right.allowed_machines[op].insert(only_machine);
                right.forbidden_machines[op].clear();
            }

            std::cout
                << "  Apply OMA branch:"
                << " left fixes J" << op.first << "o" << op.second
                << " on M" << machine
                << " | right forbids M" << machine
                << "\n";
        }
        else if (branch_dec->type == BranchDecision::Type::SEQUENCE) {
            int machine = branch_dec->order_machine;
            Op A = branch_dec->order_pair.first;
            Op B = branch_dec->order_pair.second;

            left.order_constraints[machine].insert({A, B});
            right.order_constraints[machine].insert({B, A});

            std::cout
                << "  Apply SEQUENCE branch:"
                << " machine=" << machine
                << " left J" << A.first << "o" << A.second
                << " before J" << B.first << "o" << B.second
                << " | right J" << B.first << "o" << B.second
                << " before J" << A.first << "o" << A.second
                << "\n";
        }

        // Pre-filter children's column pools and propagate constraints.
        std::vector<BranchAndPriceNode> result;

        for (auto& child : children) {
            if (child.status == NodeStatus::INFEASIBLE) {
                ++pruned_infeasible_;
                continue;
            }

            for (int m : inst_.all_machines) {
                std::vector<Column> kept;

                for (const Column& col : child.cols[m]) {
                    if (column_is_feasible_at_node(
                            col,
                            m,
                            child,
                            inst_))
                    {
                        kept.push_back(col);
                    }
                }

                child.cols[m] = std::move(kept);
            }

            if (!propagate_node(child, inst_, best_ub_)) {
                ++pruned_infeasible_;
                continue;
            }

            if (to_integer_lb(child.cp_lower_bound) >= best_ub_ - BOUND_EPS) {
                ++pruned_bound_;
                continue;
            }

            result.push_back(std::move(child));
        }

        return result;
    }

    // ── No OMA branching candidate ─────────────────────────────────────────
    ++no_branch_nodes_;

    std::cout
        << "  No branching candidate found. "
        << (node_cp_sat_enabled_
            ? "CP-SAT was already executed once after CG."
            : "node-level CP-SAT is disabled.")
        << " [no_branch_nodes="
        << no_branch_nodes_
        << "]\n";

    // Do not execute CP-SAT again in this node.
    node.status = NodeStatus::OPEN;
    return {};

    // ── Step 12: No structural branch found ───────────────────────────────────
    // Do not prune here under current policy; node simply yields no children.
    std::cout << "  No structural branch found. Node closed without prune by policy.\n";
    node.status = NodeStatus::OPEN;
    return {};
}

BranchAndPriceResult BranchAndPriceSolver::solve()
{
    auto t_start = std::chrono::steady_clock::now();
    solve_start_ = t_start;

    // Root node
    BranchAndPriceNode root;
    root.node_id    = next_node_id_++;
    root.depth      = 0;
    root.lower_bound= 0.0;

    // ── Initial UB via CP-SAT (skipped if already set by caller via set_initial_ub)
    if (!std::isfinite(best_ub_)) {
        std::cout << "[B&P] Computing initial upper bound via CP-SAT ("
                  << CP_SAT_INITIAL_TL << "s)...\n";
        auto cp_res = cp_sat_solve(inst_, root, CP_SAT_INITIAL_TL);
        if (cp_res.found) {
            update_ub(cp_res.makespan, cp_res.schedule);
        } else {
            std::cout << "[B&P] CP-SAT found no solution; using pessimistic horizon = "
                      << horizon_ << "\n";
        }
    }


    // ── Set horizon to T = UB (CG uses UB, not UB-1) ───────────────────────
    if (std::isfinite(best_ub_)) {
        int tight = static_cast<int>(best_ub_);
        if (tight < horizon_) {
            horizon_ = tight;
            std::cout << "[B&P] Horizon set to T=" << horizon_
                      << " (= UB = " << best_ub_ << ")\n";
        }
    }


    // Propagate the root node before building its initial column pool.
    //
    // This initializes:
    //   - allowed/forbidden machine domains;
    //   - start_lb/start_ub time windows;
    //   - singleton required-machine decisions;
    //   - root.cp_lower_bound.
    if (!propagate_node(
            root,
            inst_,
            best_ub_))
    {
        throw std::runtime_error(
            "Root node became infeasible during propagation."
        );
    }

    root.lower_bound = to_integer_lb(std::max(root.lower_bound, root.cp_lower_bound));


    // std::vector<Schedule> cpsat_init_schedules = incumbent_schedules_;

    // if (cpsat_init_schedules.empty() &&
    //     !best_schedule_.empty())
    // {
    //     cpsat_init_schedules.push_back(
    //         best_schedule_
    //     );
    // }

    std::vector<Schedule> cpsat_init_schedules;

    if (!best_schedule_.empty()) {
        cpsat_init_schedules.push_back(best_schedule_);
    }

    // Build initial columns for root using all stored incumbent schedules.
    root.cols =
        build_initial_columns(
            inst_,
            horizon_,
            cpsat_init_schedules
        );

    // =========================================================
    //  ROOT-NODE CG DEBUG: run CG, print results, then exit.
    // =========================================================
    static constexpr bool DEBUG_ROOT_CG = false;
    if (DEBUG_ROOT_CG) {
        std::cout << "\n=== Root-node Column Generation (debug) ===\n";

        RMP rmp(inst_, horizon_);
        for (auto& [m, cols] : root.cols)
            rmp.add_columns(m, cols);

        int init_cols = 0;
        for (int m : inst_.all_machines)
            init_cols += rmp.num_columns(m);
        std::cout << "  Initial columns loaded: " << init_cols << "\n";

        auto t0 = std::chrono::steady_clock::now();

        CGResult cg = column_generation(
            inst_,
            horizon_,
            rmp,
            CG_NODE_MAX_ITER,
            &root.required_machine,
            &root.forbidden_machines,
            &root.allowed_machines,
            &root.order_constraints,
            &root.start_lb,
            &root.start_ub,
            /*use_beam_first=*/true,
            CG_TIME_LIMIT_ROOT,
            best_ub_
        );

        double cg_sec = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  LP bound  : " << cg.lp_bound   << "\n";
        std::cout << "  LP Cmax   : " << cg.cmax_value << "\n";
        std::cout << "  Iterations: " << cg.iterations << "\n";
        std::cout << "  CG time   : " << cg_sec        << " s\n";

        int total_cols = 0;
        for (int m : inst_.all_machines)
            total_cols += (int)cg.columns[m].size();
        std::cout << "  Total cols: " << total_cols    << "\n";
        std::cout << "  Beam frac : " << std::fixed << std::setprecision(1)
                  << cg.beam_frac * 100.0 << " %\n";

        // Print chosen columns (y > LP_POS_TOL)
        std::cout << "\n  Chosen columns (y > 0):\n";
        std::cout << std::fixed << std::setprecision(6);
        for (auto& [mp, y] : cg.y_values) {
            if (y < LP_POS_TOL) continue;
            auto [m, p] = mp;
            const Column& col = cg.columns.at(m)[p];
            std::cout << "    y[m=" << m << ", p=" << p << "] = " << y
                      << "  ops: [";
            // Collect ops sorted by start time for readable output
            std::vector<std::pair<int,Op>> by_start;
            for (auto& [jop, s] : col.s_vals)
                by_start.push_back({s, jop});
            std::sort(by_start.begin(), by_start.end());
            bool first = true;
            for (auto& [s, jop] : by_start) {
                if (!first) std::cout << ", ";
                int e = col.e_vals.at(jop);
                std::cout << "(j=" << jop.first << ",o=" << jop.second
                          << ",s=" << s << ",e=" << e << ")";
                first = false;
            }
            std::cout << "]\n";
        }
        std::cout << "===========================================\n\n";

        // Verify constraint (3): assignment sums
        std::cout << "  Assignment constraint check (should all be ~1.0):\n";
        bool assign_ok = true;
        for (int j = 0; j < inst_.num_jobs; ++j) {
            for (int o = 0; o < inst_.num_ops[j]; ++o) {
                double sum = 0.0;
                for (auto& [mp, y] : cg.y_values) {
                    auto [m, p] = mp;
                    const Column& col = cg.columns.at(m)[p];
                    if (col.s_vals.count({j, o})) sum += y;
                }
                if (std::abs(sum - 1.0) > 1e-4) {
                    std::cout << "  *** VIOLATION: op(j=" << j << ",o=" << o
                              << ") sum=" << sum << "\n";
                    assign_ok = false;
                }
            }
        }
        if (assign_ok) std::cout << "  All assignment constraints satisfied.\n";

        // Verify constraint (2): machine coverage sums
        std::cout << "\n  Machine coverage check (should all be ~1.0):\n";
        for (int m : inst_.all_machines) {
            double sum = 0.0;
            for (auto& [mp, y] : cg.y_values) {
                if (mp.first == m) sum += y;
            }
            std::cout << "  machine " << m << ": sum_y = "
                      << std::fixed << std::setprecision(6) << sum << "\n";
        }

        std::cout << "\n===========================================\n\n";
        // std::exit(0);
    }
    // =========================================================

    // DFS + Best-Bound hybrid:
    //   Phase 1 (DFS)  — always dive depth-first for the first DFS_DIVE_NODES nodes,
    //                    regardless of initial UB, to find/improve integer solutions fast.
    //   Phase 2 (BB)   — after the dive budget is spent, drain DFS stack into a
    //                    best-bound min-heap and continue with best-first selection.
    static constexpr int DFS_DIVE_NODES = 50;  // nodes to explore in DFS phase

    using QElem = std::pair<double, BranchAndPriceNode>;
    auto cmp = [](const QElem& a, const QElem& b){ return a.first > b.first; };
    std::priority_queue<QElem, std::vector<QElem>, decltype(cmp)> bfs_pq(cmp);
    std::vector<BranchAndPriceNode> dfs_stack;

    // Always start with DFS dive.
    bool dfs_phase = true;
    dfs_stack.push_back(root);

    auto open_size = [&]() -> std::size_t {
        return dfs_phase ? dfs_stack.size() : bfs_pq.size();
    };

    // int nodes_explored = 0;
    // double best_lb = 0.0;
    int nodes_explored = 0;

    // Start from the propagated root lower bound,
    // not from a hard-coded zero.
    double best_lb =
        root.lower_bound;

    auto has_open = [&]{ return dfs_phase ? !dfs_stack.empty() : !bfs_pq.empty(); };

    auto compute_global_open_lb = [&]() -> double {
        double global_lb = std::numeric_limits<double>::infinity();

        if (dfs_phase) {
            for (const auto& open_node : dfs_stack) {
                global_lb = std::min(
                    global_lb,
                    open_node.lower_bound
                );
            }
        } else {
            if (!bfs_pq.empty()) {
                global_lb = bfs_pq.top().first;
            }
        }

        if (!std::isfinite(global_lb)) {
            // No open nodes left. For plotting, close the curve at UB.
            global_lb = best_ub_;
        }

        // A valid global LB curve should not be above UB in the plot.
        if (std::isfinite(best_ub_)) {
            global_lb = std::min(global_lb, best_ub_);
        }

        return global_lb;
    };

    while (has_open() &&
           nodes_explored < MAX_NODES &&
           std::chrono::duration<double>(
               std::chrono::steady_clock::now() - solve_start_
           ).count() < BAP_GLOBAL_TIME_LIMIT) {
        // Switch to best-bound once DFS dive budget is exhausted.
        if (dfs_phase && nodes_explored >= DFS_DIVE_NODES) {
            for (auto& n : dfs_stack)
                bfs_pq.push({n.lower_bound, n});
            dfs_stack.clear();
            dfs_phase = false;
            std::cout << "[B&P] DFS→Best-Bound switch after " << nodes_explored
                      << " nodes (UB=" << std::fixed << std::setprecision(4)
                      << best_ub_ << ")\n";
        }

        BranchAndPriceNode node;
        double lb;
        if (dfs_phase) {
            node = dfs_stack.back(); dfs_stack.pop_back();
            lb   = node.lower_bound;
        } else {
            auto [l, n] = bfs_pq.top(); bfs_pq.pop();
            lb = l; node = n;
        }

        // if (lb >= best_ub_ - BOUND_EPS) { ++pruned_bound_; continue; }

        // if (lb >= best_ub_ - BOUND_EPS) {
        //     ++pruned_bound_;

        //     if (on_node_cb) {
        //         on_node_cb(
        //             node.node_id,
        //             lb,
        //             best_ub_
        //         );
        //     }

        //     continue;
        // }

        if (lb >= best_ub_ - BOUND_EPS) {
            ++pruned_bound_;

            if (on_node_cb) {
                const double global_lb =
                    compute_global_open_lb();

                on_node_cb(
                    node.node_id,
                    global_lb,
                    best_ub_
                );
            }

            continue;
        }

        ++nodes_explored;
        if (nodes_explored % 10 == 0) {
            std::ostringstream line;
            line << "[B&P] Node " << nodes_explored
                 << " | queue_LB=" << std::fixed << std::setprecision(4) << lb
                 << " | UB=" << std::fixed << std::setprecision(4) << best_ub_
                 << " | Queue=" << open_size()
                 << (dfs_phase ? " [DFS]" : " [BB]")
                 << "\n";
            std::cout << line.str();
        }

        auto children = process_node(node);

        // Print a line for every terminal (pruned/infeasible/integer) node
        if (node.status != NodeStatus::OPEN) {
            const char* reason =
                node.status == NodeStatus::INFEASIBLE ? "infeasible"  :
                node.status == NodeStatus::INTEGER    ? "integrality" : "bound";
            std::ostringstream pline;
            pline << "[B&P] Node " << node.node_id
                  << " | LP_Cmax=" << std::fixed << std::setprecision(4) << node.lp_cmax
                  << " | UB=" << std::fixed << std::setprecision(4) << best_ub_
                  << " | PRUNED (" << reason << ")\n";
            std::cout << pline.str();
        }

        for (auto& child : children) {
            best_lb = std::max(best_lb, child.lower_bound);
            if (dfs_phase) dfs_stack.push_back(child);
            else           bfs_pq.push({child.lower_bound, child});
        }

        // if (on_node_cb)
        //     on_node_cb(node.node_id, node.lower_bound, best_ub_);

        if (on_node_cb) {
            const double global_lb =
                compute_global_open_lb();

            on_node_cb(
                node.node_id,
                global_lb,
                best_ub_
            );
        }
    }

    auto t_end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    const bool global_time_limit_reached =
        elapsed >= BAP_GLOBAL_TIME_LIMIT;

    // bool proven_optimal = (bfs_pq.empty() && dfs_stack.empty() && nodes_explored < MAX_NODES);

    const bool proven_optimal =
    bfs_pq.empty() &&
    dfs_stack.empty() &&
    nodes_explored < MAX_NODES &&
    !global_time_limit_reached &&
    !heuristic_cg_prune_occurred_;

    BranchAndPriceResult res;
    res.best_bound      = proven_optimal ? best_ub_ : best_lb;
    res.best_obj        = best_ub_;
    res.best_schedule   = best_schedule_;
    res.proven_optimal  = proven_optimal;
    res.nodes_explored    = nodes_explored;
    res.integer_nodes     = integer_nodes_;
    res.pruned_bound      = pruned_bound_;
    res.pruned_infeasible  = pruned_infeasible_;
    res.elapsed_seconds    = elapsed;
    res.total_pricing_time = total_pricing_time_;
    res.total_master_time  = total_master_time_;
    res.prec_graph_calls   = prec_graph_calls_;
    res.mean_prec_graph_ms = (prec_graph_calls_ > 0)
        ? (prec_graph_total_time_ / prec_graph_calls_) * 1000.0 : 0.0;
    res.total_cols_generated = total_cols_generated_;
    res.avg_beam_frac        = (cg_node_count_ > 0)
        ? sum_beam_frac_ / cg_node_count_ : 0.0;

    std::cout << "[B&P] no_branch_nodes (select_branch found nothing): "
              << no_branch_nodes_ << "\n";
    std::cout << "[B&P] Total columns generated (tree): " << total_cols_generated_ << "\n";
    std::cout << "[B&P] Avg P/M (pricing fraction):     "
              << std::fixed << std::setprecision(1)
              << res.avg_beam_frac * 100.0 << "%\n";
    if (global_time_limit_reached) {
        std::cout << "[B&P] Global time limit reached: "
                  << BAP_GLOBAL_TIME_LIMIT << " s\n";
    }

    return res;
}
