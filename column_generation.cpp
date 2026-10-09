// =============================================================================
// column_generation.cpp  –  Machine-based column generation for FJSSP
//
// Implements column_generation.h.
// =============================================================================

#include "column_generation.h"
#include "instances.h"

// CPLEX Concert Technology (IL_STD defined in column_generation.h via ilocplex.h)

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <future>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <iomanip>

// ---------------------------------------------------------------------------
// Helper: infinity
// ---------------------------------------------------------------------------
static constexpr double INF = std::numeric_limits<double>::infinity();
static constexpr int    IINF = std::numeric_limits<int>::max() / 2;

static inline double integer_lb_from_cmax(double lp_cmax)
{
    if (!std::isfinite(lp_cmax)) {
        return lp_cmax;
    }

    return std::ceil(lp_cmax - 1e-6);
}

// ===========================================================================
//  compute_horizon
// ===========================================================================

int compute_horizon(const Instance& inst) {
    int h = 0;
    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            int best = 0;
            for (auto& [m, pt] : inst.proc_time[j][o])
                best = std::max(best, pt);
            h += best;
        }
    }
    return h;
}

// ===========================================================================
//  Forward declarations needed by build_initial_columns
// ===========================================================================

static JobsOps build_jobs_ops(
    int machine, const Instance& inst,
    const std::map<Op, int>*           required_machine,
    const std::map<Op, std::set<int>>* forbidden_machines,
    const std::map<Op, std::set<int>>* allowed_machines);

// ===========================================================================
//  Construction heuristics
// ===========================================================================

/// Build a column for each machine using earliest-finish-time dispatch.
/// `job_order` specifies the order in which jobs are dispatched;
/// -1 in job_order means randomize.
static std::map<int, Column>
_heuristic_schedule(const Instance& inst, int horizon,
                    const std::vector<int>& job_order_hint,
                    std::mt19937* rng)
{
    int J = inst.num_jobs;
    // next_op[j] = index of next unscheduled op for job j
    std::vector<int> next_op(J, 0);
    // job_ready[j] = earliest time job j can start its next op
    std::vector<int> job_ready(J, 0);
    // mach_ready[m] = earliest free time for machine m
    std::map<int, int> mach_ready;
    for (int m : inst.all_machines) mach_ready[m] = 0;

    // Columns under construction
    std::map<int, std::vector<SchedEntry>> col_ops;
    for (int m : inst.all_machines) col_ops[m] = {};

    // std::vector<int> order = job_order_hint;
    // std::vector<int> order;
    // order.resize(J);
    // std::iota(order.begin(), order.end(), 0);
    // if (rng) std::shuffle(order.begin(), order.end(), *rng);

    std::vector<int> order = job_order_hint;

    if (order.empty()) {
        order.resize(J);
        std::iota(order.begin(), order.end(), 0);
    }

    if (rng != nullptr) {
        std::shuffle(order.begin(), order.end(), *rng);
    }

    // Schedule all ops in round-robin over job order
    int total_ops = 0;
    for (int j = 0; j < J; ++j) total_ops += inst.num_ops[j];

    for (int iter = 0; iter < total_ops; ++iter) {
        // Find the job/op with earliest finish time
        int best_j = -1, best_o = -1, best_m = -1, best_s = IINF, best_e = IINF;

        for (int j : order) {
            int o = next_op[j];
            if (o >= inst.num_ops[j]) continue;
            for (int m : inst.elig_mach[j][o]) {
                int pt = inst.proc_time[j][o].at(m);
                int s  = std::max(job_ready[j], mach_ready.at(m));
                int e  = s + pt;
                if (e < best_e || (e == best_e && s < best_s)) {
                    best_e = e; best_s = s;
                    best_j = j; best_o = o; best_m = m;
                }
            }
        }
        if (best_j < 0) break;

        col_ops[best_m].push_back({best_j, best_o, best_s, best_e});
        mach_ready[best_m] = best_e;
        job_ready[best_j]  = best_e;
        next_op[best_j]++;
    }

    std::map<int, Column> result;
    for (int m : inst.all_machines)
        result[m] = make_column(m, col_ops[m]);
    return result;
}

/// Dispatching schedule: sort remaining eligible ops by rule and pick best.
/// rule == 0: EFT, rule == 1: SPT (shortest processing time first),
/// rule == 2: LPT, rule == 3: MWKR (most work remaining).
/// When rng is provided, job iteration order is shuffled each step for diversity.
static std::map<int, Column>
_dispatching_schedule(const Instance& inst, int horizon, int rule,
                      std::mt19937* rng = nullptr)
{
    int J = inst.num_jobs;
    std::vector<int> next_op(J, 0);
    std::vector<int> job_ready(J, 0);
    std::map<int, int> mach_ready;
    for (int m : inst.all_machines) mach_ready[m] = 0;

    // Work remaining per job (sum of min proc times of remaining ops)
    auto work_remaining = [&](int j) -> int {
        int w = 0;
        for (int o = next_op[j]; o < inst.num_ops[j]; ++o) {
            int mn = IINF;
            for (auto& [m, pt] : inst.proc_time[j][o]) mn = std::min(mn, pt);
            w += mn;
        }
        return w;
    };

    std::map<int, std::vector<SchedEntry>> col_ops;
    for (int m : inst.all_machines) col_ops[m] = {};

    int total_ops = 0;
    for (int j = 0; j < J; ++j) total_ops += inst.num_ops[j];

    // Base job iteration order; shuffled each step when rng is provided
    std::vector<int> job_order(J);
    std::iota(job_order.begin(), job_order.end(), 0);

    for (int iter = 0; iter < total_ops; ++iter) {
        if (rng) std::shuffle(job_order.begin(), job_order.end(), *rng);

        int best_j = -1, best_o = -1, best_m = -1, best_s = IINF, best_e = IINF;
        double best_key = INF;

        for (int j : job_order) {
            int o = next_op[j];
            if (o >= inst.num_ops[j]) continue;
            for (int m : inst.elig_mach[j][o]) {
                int pt = inst.proc_time[j][o].at(m);
                int s  = std::max(job_ready[j], mach_ready.at(m));
                int e  = s + pt;
                double key;
                switch (rule) {
                    case 0: key = (double)e; break;                // EFT
                    case 1: key = (double)pt; break;               // SPT
                    case 2: key = -(double)pt; break;              // LPT
                    case 3: key = -(double)work_remaining(j); break; // MWKR
                    default: key = (double)e; break;
                }
                if (key < best_key) {
                    best_key = key;
                    best_e = e; best_s = s;
                    best_j = j; best_o = o; best_m = m;
                }
            }
        }
        if (best_j < 0) break;

        col_ops[best_m].push_back({best_j, best_o, best_s, best_e});
        mach_ready[best_m] = best_e;
        job_ready[best_j]  = best_e;
        next_op[best_j]++;
    }

    std::map<int, Column> result;
    for (int m : inst.all_machines)
        result[m] = make_column(m, col_ops[m]);
    return result;
}

/// Canonical key for column deduplication: (machine, ops sorted by start)
static std::vector<SchedEntry>
_col_key(const Column& col) {
    return col.ops; // already sorted by start
}

std::map<int, std::vector<Column>>
build_initial_columns(const Instance& inst, int horizon,
                      const std::vector<Schedule>& cpsat_schedules)
{
    static constexpr int INIT_COL_BUDGET_MULT = 3;  // user-tunable seeding intensity 8

    std::map<int, std::vector<Column>> result;
    for (int m : inst.all_machines) result[m] = {};

    // Track deduplicated keys per machine
    std::map<int, std::set<std::vector<SchedEntry>>> seen;
    for (int m : inst.all_machines) seen[m] = {};

    auto add = [&](const std::map<int, Column>& sched) {
        for (auto& [m, col] : sched) {
            // Discard any column whose operations exceed the horizon
            bool within_horizon = true;
            for (auto& [jop, e] : col.e_vals)
                if (e > horizon) { within_horizon = false; break; }
            if (!within_horizon) continue;
            auto k = _col_key(col);
            if (seen[m].insert(k).second)
                result[m].push_back(col);
        }
    };

    auto make_random_rng = []() {
        std::random_device rd;
        std::seed_seq seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        return std::mt19937(seq);
    };

    // Greedy EFT: fully randomized variants
    {
        std::mt19937 rng = make_random_rng();
        add(_heuristic_schedule(inst, horizon, {}, &rng));
    }
    for (int run = 0; run < 49 * INIT_COL_BUDGET_MULT; ++run) {
        std::mt19937 rng = make_random_rng();
        add(_heuristic_schedule(inst, horizon, {}, &rng));
    }

    // Dispatching rules: fully randomized variants
    for (int rule = 0; rule < 4; ++rule) {
        {
            std::mt19937 rng = make_random_rng();
            add(_dispatching_schedule(inst, horizon, rule, &rng));
        }
        for (int run = 0; run < 12 * INIT_COL_BUDGET_MULT; ++run) {
            std::mt19937 rng = make_random_rng();
            add(_dispatching_schedule(inst, horizon, rule, &rng));
        }
    }

    // Additional randomized EFT variants with diverse seeds (independent RNGs)
    for (int run = 0; run < 300 * INIT_COL_BUDGET_MULT; ++run) {
        std::mt19937 rng = make_random_rng();
        add(_heuristic_schedule(inst, horizon, {}, &rng));
    }

    // 305+. Per-machine DP pricing with zero duals (rc_tol = +inf accepts all).
    // Generates machine-local column diversity independent of global dispatch.
    {
        const int DP_INIT_TOP_K = 50 * INIT_COL_BUDGET_MULT;
        Duals zero_duals;
        for (int machine : inst.all_machines) {
            JobsOps jobs_ops = build_jobs_ops(machine, inst, nullptr, nullptr, nullptr);
            if (jobs_ops.empty()) continue;

            // auto dp_res = dp_pricing(
            //     machine, jobs_ops, inst, zero_duals, horizon,
            //     {}, {}, {},
            //     DP_INIT_TOP_K,
            //     std::numeric_limits<double>::infinity(), 
            //     /*allow_intentional_idle=*/false);

            auto dp_res = dp_pricing(
                machine,
                jobs_ops,
                inst,
                zero_duals,
                horizon,
                {},   // order_constraints
                {},   // required_here
                {},   // start_lb
                {},   // start_ub
                DP_INIT_TOP_K,
                std::numeric_limits<double>::infinity(),
                false
            );

            for (auto& pr : dp_res) {
                bool within_horizon = true;
                for (auto& [jop, e] : pr.col.e_vals)
                    if (e > horizon) { within_horizon = false; break; }
                if (!within_horizon) continue;
                auto k = _col_key(pr.col);
                if (seen[machine].insert(k).second)
                    result[machine].push_back(pr.col);
            }
        }
    }

    // CP-SAT schedules: decompose each full schedule into per-machine columns.
    // Each schedule assigns every operation to exactly one machine; we extract
    // the sub-column for each machine and add it to the pool if it is new.
    for (auto& sched : cpsat_schedules) {
        // Group by machine
        std::map<int, std::vector<SchedEntry>> mach_entries;
        for (auto& [jop, mse] : sched) {
            auto [m, s, e] = mse;
            if (!result.count(m)) continue;  // skip machines not in instance
            mach_entries[m].push_back({jop.first, jop.second, s, e});
        }
        for (auto& [m, entries] : mach_entries) {
            bool within_horizon = true;
            for (auto& [j, o, s, e] : entries)
                if (e > horizon) { within_horizon = false; break; }
            if (!within_horizon) continue;
            if (entries.empty()) continue;
            Column col = make_column(m, entries);
            auto k = _col_key(col);
            if (seen[m].insert(k).second)
                result[m].push_back(col);
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// Repair pricing
// ---------------------------------------------------------------------------

// Forward declaration (defined after RMP)
static JobsOps build_jobs_ops(
    int machine, const Instance& inst,
    const std::map<Op, int>*           required_machine,
    const std::map<Op, std::set<int>>* forbidden_machines,
    const std::map<Op, std::set<int>>* allowed_machines);

std::vector<Column>
repair_pricing_for_op(
    int machine,
    Op required_op,
    const Instance& inst,
    int horizon,
    const std::map<Op, int>& required_machine,
    const std::map<Op, std::set<int>>& forbidden_machines,
    const std::map<Op, std::set<int>>& allowed_machines,
    const std::map<int, std::set<std::pair<Op, Op>>>& order_constraints,
    const std::map<Op, int>& start_lb,
    const std::map<Op, int>& start_ub)
{
    JobsOps jobs_ops = build_jobs_ops(machine, inst,
        &required_machine, &forbidden_machines, &allowed_machines);
    if (jobs_ops.empty()) return {};



    std::set<Op> machine_ops;

    for (const auto& [job, operations] : jobs_ops) {
        for (int operation : operations) {
            machine_ops.insert({
                job,
                operation
            });
        }
    }

    // All operations fixed to this machine plus the target operation.
    std::set<Op> req_here;

    for (const auto& [operation, required_machine_id] :
        required_machine)
    {
        if (required_machine_id == machine) {
            req_here.insert(operation);
        }
    }

    req_here.insert(required_op);

    // Never allow a mandatory operation to be silently omitted.
    for (const Op& operation : req_here) {
        if (!machine_ops.count(operation)) {
            return {};
        }
    }

    // Order constraints scoped to this machine
    // std::set<std::pair<Op, Op>> m_oc;
    // auto oc_it = order_constraints.find(machine);
    // if (oc_it != order_constraints.end()) m_oc = oc_it->second;

    // Run DP with zero duals and rc_tol = +inf: find any valid column containing required_op
    Duals zero_duals;


    std::set<std::pair<Op, Op>> m_oc;

    auto oc_it = order_constraints.find(machine);
    if (oc_it != order_constraints.end()) {
        m_oc = oc_it->second;
    }

    auto results = dp_pricing(
        machine,
        jobs_ops,
        inst,
        zero_duals,
        horizon,
        m_oc,
        req_here,
        start_lb,
        start_ub,
        1,
        std::numeric_limits<double>::infinity(),
        false
    );

    std::vector<Column> cols;
    cols.reserve(results.size());
    for (auto& pr : results)
        cols.push_back(std::move(pr.col));
    return cols;
}

std::vector<Column>
repair_pricing_for_machine(
    int machine,
    const std::set<Op>& mandatory_ops,
    const Instance& inst,
    int horizon,
    const std::map<Op, int>& required_machine,
    const std::map<Op, std::set<int>>& forbidden_machines,
    const std::map<Op, std::set<int>>& allowed_machines,
    const std::map<Op, int>& start_lb,
    const std::map<Op, int>& start_ub,
    int top_k)
{
    if (mandatory_ops.empty()) {
        return {};
    }

    JobsOps jobs_ops = build_jobs_ops(
        machine,
        inst,
        &required_machine,
        &forbidden_machines,
        &allowed_machines
    );

    if (jobs_ops.empty()) {
        return {};
    }

    // Collect every operation that can actually appear on this machine
    // under the current branching restrictions.
    std::set<Op> machine_ops;

    for (const auto& [job, operations] : jobs_ops) {
        for (int op : operations) {
            machine_ops.insert({job, op});
        }
    }

    // Every mandatory operation must still be reachable on this machine.
    // Do not allow dp_pricing() to silently ignore a missing mandatory op.
    for (const Op& mandatory_op : mandatory_ops) {
        if (!machine_ops.count(mandatory_op)) {
            return {};
        }

        auto required_it = required_machine.find(mandatory_op);

        if (required_it == required_machine.end() ||
            required_it->second != machine)
        {
            return {};
        }

        auto forbidden_it =
            forbidden_machines.find(mandatory_op);

        if (forbidden_it != forbidden_machines.end() &&
            forbidden_it->second.count(machine))
        {
            return {};
        }

        auto allowed_it =
            allowed_machines.find(mandatory_op);

        if (allowed_it != allowed_machines.end() &&
            !allowed_it->second.count(machine))
        {
            return {};
        }
    }

    Duals zero_duals;

    auto pricing_results = dp_pricing(
        machine,
        jobs_ops,
        inst,
        zero_duals,
        horizon,
        {}, // order_constraints
        mandatory_ops,
        start_lb,
        start_ub,
        std::max(1, top_k),
        std::numeric_limits<double>::infinity()
    );

    std::vector<Column> result;
    result.reserve(pricing_results.size());

    for (auto& pricing_result : pricing_results) {
        Column& column = pricing_result.col;

        bool valid = true;

        // Verify explicitly that every mandatory operation is present.
        for (const Op& mandatory_op : mandatory_ops) {
            if (!column.a.count(mandatory_op)) {
                valid = false;
                break;
            }
        }

        if (!valid) {
            continue;
        }

        // Verify that the complete column lies inside the current horizon.
        for (const auto& [operation, end] : column.e_vals) {
            if (end > horizon) {
                valid = false;
                break;
            }
        }

        if (!valid) {
            continue;
        }

        result.push_back(std::move(column));
    }

    return result;
}

// ===========================================================================
//  RMP implementation
// ===========================================================================

RMP::RMP(const Instance& inst, int horizon)
    : inst_(inst), horizon_(horizon)
{
    env_   = IloEnv();
    model_ = IloModel(env_);
    cplex_ = IloCplex(model_);
    cplex_.setOut(env_.getNullStream());
    cplex_.setWarning(env_.getNullStream());
    cplex_.setParam(IloCplex::Param::Simplex::Display, 0);
    // Use dual simplex — matches Python/CPLEX behaviour
    cplex_.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);

    build_base_constraints();
}

RMP::~RMP()
{
    env_.end();
}

void RMP::build_base_constraints() {
    // C_max variable: 0 <= C_max
    cmax_var_ = IloNumVar(env_, 0.0, IloInfinity, "Cmax");
    model_.add(cmax_var_);

    // Objective: min C_max
    obj_ = IloMinimize(env_, cmax_var_);
    model_.add(obj_);

    int J = inst_.num_jobs;

    // Artificial variables and machine coverage constraints (2):
    //   sum_{p in P'_m} y_{m,p} + art[m] = 1   for all m
    for (int m : inst_.all_machines) {
        std::string aname = "art_" + std::to_string(m);
        IloNumVar art(env_, 0.0, IloInfinity, aname.c_str());
        model_.add(art);
        art_vars_[m] = art;
        obj_.setLinearCoef(art, BIG_M);

        std::string cname = "mach_" + std::to_string(m);
        // Start: only the artificial variable; y columns added incrementally
        IloRange ctr(env_, 1.0, art, 1.0, cname.c_str());
        model_.add(ctr);
        mach_cov_ctr_[m] = ctr;
    }

    // Op assignment constraints (3): sum_{m,p} a[j,o,m,p]*y[m,p] = 1
    for (int j = 0; j < J; ++j) {
        for (int o = 0; o < inst_.num_ops[j]; ++o) {
            Op jop{j, o};
            std::string aname = "art_assign_" + std::to_string(j) + "_" + std::to_string(o);
            IloNumVar art(env_, 0.0, IloInfinity, aname.c_str());
            model_.add(art);
            obj_.setLinearCoef(art, BIG_M);
            assign_art_vars_[jop] = art;

            std::string cname = "assign_" + std::to_string(j) + "_" + std::to_string(o);
            IloRange ctr(env_, 1.0, art, 1.0, cname.c_str());
            model_.add(ctr);
            op_assign_ctr_[jop] = ctr;
        }
    }

    // Precedence constraints (4): e[j,o] <= s[j,o+1]  =>  e[j,o] - s[j,o+1] <= 0
    for (int j = 0; j < J; ++j) {
        for (int o = 0; o < inst_.num_ops[j] - 1; ++o) {
            Op jop{j, o};
            std::string cname = "prec_" + std::to_string(j) + "_" + std::to_string(o);
            IloRange ctr(env_, -IloInfinity, 0.0, cname.c_str());
            model_.add(ctr);
            prec_ctr_[jop] = ctr;
        }
    }

    // Makespan constraints (5): e[j,last_o] <= C_max  =>  e[j,last_o] - C_max <= 0
    for (int j = 0; j < J; ++j) {
        std::string cname = "mkspan_" + std::to_string(j);
        IloRange ctr(env_, -IloInfinity, 0.0, cname.c_str());
        ctr.setLinearCoef(cmax_var_, -1.0);
        model_.add(ctr);
        mkspan_ctr_[j] = ctr;
    }

    // Initialise per-machine column storage
    for (int m : inst_.all_machines)
        machine_cols_[m] = {};
}

void RMP::add_column_to_lp(int machine, const Column& col, IloNumVar var)
{
    // (2) machine coverage
    mach_cov_ctr_.at(machine).setLinearCoef(var, 1.0);

    int J = inst_.num_jobs;

    // Accumulate coefficients for precedence / makespan constraints first,
    // because a single variable can contribute from two different ops in the
    // same column (e[j,o] to prec[j,o] AND -s[j,o] to prec[j,o-1]).
    std::map<Op, double>  prec_coef;
    std::map<int, double> mksp_coef;

    for (auto& [jop, _s] : col.s_vals) {
        auto [j, o] = jop;
        int e = col.e_vals.at(jop);
        int s = _s;

        // (3) op assignment
        op_assign_ctr_.at(jop).setLinearCoef(var, 1.0);

        // (4) precedence: +e[j,o] to prec_ctr[(j,o)]  if o < nj-1
        if (o < inst_.num_ops[j] - 1)
            prec_coef[{j, o}] += (double)e;

        // (4) precedence: -s[j,o] to prec_ctr[(j,o-1)]  if o > 0
        if (o > 0)
            prec_coef[{j, o - 1}] -= (double)s;

        // (5) makespan: +e[j,last_o]
        if (o == inst_.num_ops[j] - 1)
            mksp_coef[j] += (double)e;
    }

    for (auto& [jop, coef] : prec_coef)
        if (prec_ctr_.count(jop))
            prec_ctr_.at(jop).setLinearCoef(var, coef);

    for (auto& [j, coef] : mksp_coef)
        mkspan_ctr_.at(j).setLinearCoef(var, coef);
}


int RMP::add_columns(
    int machine,
    const std::vector<Column>& new_cols)
{
    auto& machine_columns =
        machine_cols_[machine];

    // Existing signatures on this machine.
    std::set<std::vector<SchedEntry>> existing_keys;

    for (const Column& existing_column :
         machine_columns.cols)
    {
        existing_keys.insert(
            _col_key(existing_column)
        );
    }

    int actually_added = 0;

    for (const Column& column : new_cols) {
        const auto key =
            _col_key(column);

        // The column already exists either in the RMP or earlier
        // in new_cols during this same call.
        if (!existing_keys.insert(key).second) {
            continue;
        }

        const std::string variable_name =
            "y_" +
            std::to_string(machine) +
            "_" +
            std::to_string(
                machine_columns.cols.size()
            );

        IloNumVar variable(
            env_,
            0.0,
            IloInfinity,
            variable_name.c_str()
        );

        model_.add(variable);

        add_column_to_lp(
            machine,
            column,
            variable
        );

        machine_columns.cols.push_back(
            column
        );

        machine_columns.vars.push_back(
            variable
        );

        ++actually_added;
    }

    return actually_added;
}

double RMP::solve(double time_limit_sec) {
    try {
        // The RMP is solved by dual simplex. Apply a per-call limit so that
        // one CPLEX solve cannot overrun the CG/B&P global deadline.
        if (time_limit_sec > 0.0) {
            cplex_.setParam(
                IloCplex::Param::TimeLimit,
                time_limit_sec
            );
        } else {
            cplex_.setParam(
                IloCplex::Param::TimeLimit,
                1e+75
            );
        }

        bool ok = cplex_.solve();
        auto status = cplex_.getStatus();
        optimal_    = ok && (status == IloAlgorithm::Optimal);
        infeasible_ = (status == IloAlgorithm::Infeasible);

        // Do not let a time limit leak into later independent solves.
        cplex_.setParam(IloCplex::Param::TimeLimit, 1e+75);

        if (optimal_)
            return cplex_.getObjValue();
    } catch (IloException&) {
        cplex_.setParam(IloCplex::Param::TimeLimit, 1e+75);
        optimal_    = false;
        infeasible_ = true;
    }
    return INF;
}

bool RMP::is_optimal()    const { return optimal_;    }
bool RMP::is_infeasible() const { return infeasible_; }
double RMP::objective()   const { return cplex_.getObjValue(); }
double RMP::get_cmax()    const { return cplex_.getValue(cmax_var_); }

RMP::MIPResult RMP::solve_mip(double time_limit_sec, double ub_cutoff) {
    MIPResult res;

    // Collect all column (y) variables — skip artificials.
    IloNumVarArray y_vars(env_);
    for (auto& [m, mc] : machine_cols_)
        for (auto& v : mc.vars)
            y_vars.add(v);
    if (y_vars.getSize() == 0) return res;

    // Tighten Cmax upper bound: require strict improvement (Cmax is integer).
    double old_cmax_ub = cmax_var_.getUB();
    if (std::isfinite(ub_cutoff) && ub_cutoff - 1.0 < old_cmax_ub)
        cmax_var_.setUB(ub_cutoff - 1.0);

    // Convert y variables to binary for the MIP solve.
    IloConversion conv(env_, y_vars, ILOBOOL);
    model_.add(conv);

    // Suppress output; set time limit.
    cplex_.setParam(IloCplex::Param::MIP::Display, 0);
    if (time_limit_sec > 0.0)
        cplex_.setParam(IloCplex::Param::TimeLimit, time_limit_sec);

    try {
        bool ok = cplex_.solve();
        auto status = cplex_.getStatus();
        bool feasible = ok && (status == IloAlgorithm::Optimal ||
                               status == IloAlgorithm::Feasible);
        if (feasible) {
            // Only accept solutions where all artificials are driven to zero.
            bool art_ok = true;
            for (auto& [m, art] : art_vars_)
                if (cplex_.getValue(art) > LP_POS_TOL) { art_ok = false; break; }
            if (art_ok)
                for (auto& [jop, art] : assign_art_vars_)
                    if (cplex_.getValue(art) > LP_POS_TOL) { art_ok = false; break; }
            if (art_ok) {
                res.found    = true;
                res.makespan = cplex_.getValue(cmax_var_);
                for (auto& [m, mc] : machine_cols_)
                    for (int p = 0; p < (int)mc.vars.size(); ++p) {
                        double v = cplex_.getValue(mc.vars[p]);
                        if (v > LP_POS_TOL) res.y_values[{m, p}] = v;
                    }
            }
        }
    } catch (IloException&) {}

    // Restore: remove binary conversion, restore Cmax UB, re-solve LP.
    model_.remove(conv);
    conv.end();
    y_vars.end();
    cmax_var_.setUB(old_cmax_ub);
    cplex_.setParam(IloCplex::Param::TimeLimit, 1e+75);
    cplex_.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);
    try {
        bool ok = cplex_.solve();
        auto status = cplex_.getStatus();
        optimal_    = ok && (status == IloAlgorithm::Optimal);
        infeasible_ = (status == IloAlgorithm::Infeasible);
    } catch (IloException&) {
        optimal_ = false; infeasible_ = true;
    }

    return res;
}

double RMP::get_artificial(int machine) const {
    auto it = art_vars_.find(machine);
    if (it == art_vars_.end()) return 0.0;
    return cplex_.getValue(it->second);
}

double RMP::get_assign_artificial(Op jop) const {
    auto it = assign_art_vars_.find(jop);
    if (it == assign_art_vars_.end()) return 0.0;
    return cplex_.getValue(it->second);
}

const std::vector<Column>& RMP::get_columns(int machine) const {
    return machine_cols_.at(machine).cols;
}

int RMP::num_columns(int machine) const {
    auto it = machine_cols_.find(machine);
    if (it == machine_cols_.end()) return 0;
    return static_cast<int>(it->second.cols.size());
}

Duals RMP::get_duals() const {
    Duals d;
    for (auto& [m, ctr] : mach_cov_ctr_)
        d.alpha[m] = cplex_.getDual(ctr);
    for (auto& [jop, ctr] : op_assign_ctr_)
        d.beta[jop] = cplex_.getDual(ctr);
    for (auto& [jop, ctr] : prec_ctr_)
        d.theta[jop] = cplex_.getDual(ctr);
    for (auto& [j, ctr] : mkspan_ctr_)
        d.mu[j] = cplex_.getDual(ctr);
    return d;
}

std::map<std::pair<int,int>, double> RMP::get_y_values() const {
    std::map<std::pair<int,int>, double> y;
    for (auto& [m, mc] : machine_cols_) {
        for (int p = 0; p < (int)mc.vars.size(); ++p) {
            double v = cplex_.getValue(mc.vars[p]);
            if (v > LP_POS_TOL) y[{m, p}] = v;
        }
    }
    return y;
}

void RMP::rebuild_with_feasible_cols(
    const std::map<int, std::vector<Column>>& cols_by_machine)
{
    // Free all existing Concert objects (IloNumVar, IloRange, etc.)
    env_.end();

    // Reinitialise from scratch
    env_   = IloEnv();
    model_ = IloModel(env_);
    cplex_ = IloCplex(model_);
    cplex_.setOut(env_.getNullStream());
    cplex_.setWarning(env_.getNullStream());
    cplex_.setParam(IloCplex::Param::Simplex::Display, 0);
    cplex_.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);

    cmax_var_ = IloNumVar();  // empty handle
    art_vars_.clear();
    assign_art_vars_.clear();
    machine_cols_.clear();
    mach_cov_ctr_.clear();
    op_assign_ctr_.clear();
    prec_ctr_.clear();
    mkspan_ctr_.clear();
    optimal_    = false;
    infeasible_ = false;

    build_base_constraints();

    for (auto& [m, cols] : cols_by_machine)
        add_columns(m, cols);
}

// ===========================================================================
//  Reduced cost computation
// ===========================================================================

double compute_reduced_cost(int machine, const Column& col,
                             const Duals& duals, const Instance& inst)
{
    double rc = 0.0;

    // -alpha[m]
    {
        auto it = duals.alpha.find(machine);
        if (it != duals.alpha.end()) rc -= it->second;
    }

    for (auto& [jop, _s] : col.s_vals) {
        auto [j, o] = jop;
        int s = _s;
        int e = col.e_vals.at(jop);

        // -beta[j,o]
        {
            auto it = duals.beta.find(jop);
            if (it != duals.beta.end()) rc -= it->second;
        }

        // precedence (theta): -theta[j,o] * e[j,o]  if o < nj-1
        if (o < inst.num_ops[j] - 1) {
            auto it = duals.theta.find(jop);
            if (it != duals.theta.end()) rc -= it->second * e;
        }
        // precedence: +theta[j,o-1] * s[j,o]  if o > 0
        if (o > 0) {
            auto it = duals.theta.find({j, o - 1});
            if (it != duals.theta.end()) rc += it->second * s;
        }

        // makespan: -mu[j] * e[j,last_o]
        if (o == inst.num_ops[j] - 1) {
            auto it = duals.mu.find(j);
            if (it != duals.mu.end()) rc -= it->second * e;
        }
    }

    return rc;
}

// ===========================================================================
//  Per-operation weight for pricing
// ===========================================================================

/// w(j, o, s, p) = contribution of scheduling op (j,o) at time [s, s+p] to reduced cost.
static double _weight(int j, int o, int s, int p,
                       const Instance& inst, const Duals& duals)
{
    double w = 0.0;
    int e = s + p;
    Op jop{j, o};
    int nj = inst.num_ops[j];

    // -beta[j,o]
    {
        auto it = duals.beta.find(jop);
        if (it != duals.beta.end()) w -= it->second;
    }

    // theta contribution for this op's END (if o < nj-1):
    //   precedence constraint: -theta[j,o]*e
    if (o < nj - 1) {
        auto it = duals.theta.find(jop);
        if (it != duals.theta.end()) w -= it->second * e;
    }

    // theta contribution for this op's START (if o > 0):
    //   precedence constraint: +theta[j,o-1]*s
    if (o > 0) {
        auto it = duals.theta.find({j, o - 1});
        if (it != duals.theta.end()) w += it->second * s;
    }

    // mu contribution for LAST op's end:
    //   makespan constraint: -mu[j]*e
    if (o == nj - 1) {
        auto it = duals.mu.find(j);
        if (it != duals.mu.end()) w -= it->second * e;
    }

    return w;
}

// ===========================================================================
//  Beam search pricing
// ===========================================================================

struct BeamLabel {
    double score = 0.0;  // accumulated weighted cost
    int    t     = 0;    // time when the machine becomes free

    // Operation identities in machine order.
    // This can be retained for compatibility/debugging.
    std::vector<std::pair<int, int>> ops;

    // Exact operations and times used when this label was evaluated:
    // (job, operation, start, end).
    //
    // These entries must be used directly to create the final column.
    std::vector<SchedEntry> entries;

    // Last operation index placed for every local job.
    std::vector<int> last_o_placed;

    // Whether each local operation has already been placed.
    std::vector<bool> placed;
};

std::vector<PricingResult>
beam_search_pricing(int machine, const JobsOps& jobs_ops,
                    const Instance& inst, const Duals& duals,
                    int horizon,
                    const std::set<std::pair<Op, Op>>& order_constraints,
                    const std::set<Op>&                required_here,
                    const std::map<Op, int>&            start_lb_ops,
                    const std::map<Op, int>&            start_ub_ops,
                    int top_k)
{
    if (jobs_ops.empty()) return {};

    const int n_local_jobs = (int)jobs_ops.size();

    // ── Build local flat tables ───────────────────────────────────────────────
    // Assign consecutive local IDs to every (j,o) visible to this machine.
    std::vector<int> local_start(n_local_jobs + 1, 0);
    for (int ji = 0; ji < n_local_jobs; ++ji)
        local_start[ji + 1] = local_start[ji] + (int)jobs_ops[ji].second.size();
    const int n_local_ops = local_start[n_local_jobs];

    std::vector<int> local_j(n_local_ops), local_o(n_local_ops);
    // job j -> its index in jobs_ops (-1 if absent)
    std::vector<int> job_to_ji(inst.num_jobs, -1);
    for (int ji = 0; ji < n_local_jobs; ++ji) {
        int j = jobs_ops[ji].first;
        if (j < (int)job_to_ji.size()) job_to_ji[j] = ji;
        for (int qi = 0; qi < (int)jobs_ops[ji].second.size(); ++qi) {
            int lid = local_start[ji] + qi;
            local_j[lid] = j;
            local_o[lid] = jobs_ops[ji].second[qi];
        }
    }

    // (j, o) -> local_id (-1 if not present)
    auto op_to_lid = [&](int j, int o) -> int {
        if (j < 0 || j >= (int)job_to_ji.size()) return -1;
        int ji = job_to_ji[j];
        if (ji < 0) return -1;
        for (int qi = 0; qi < (int)jobs_ops[ji].second.size(); ++qi)
            if (jobs_ops[ji].second[qi] == o) return local_start[ji] + qi;
        return -1;
    };

    std::vector<std::vector<int>> successors_local(n_local_ops);

    for (const auto& [A, B] : order_constraints) {
        int lid_A = op_to_lid(A.first, A.second);
        int lid_B = op_to_lid(B.first, B.second);

        if (lid_A >= 0 && lid_B >= 0) {
            successors_local[lid_A].push_back(lid_B);
        }
    }

    // Cache processing times and time-window bounds as flat arrays
    std::vector<int> pt_flat(n_local_ops);
    std::vector<int> s_lb_flat(n_local_ops, 0);
    std::vector<int> s_ub_flat(n_local_ops);
    for (int lid = 0; lid < n_local_ops; ++lid) {
        int j = local_j[lid], o = local_o[lid];
        int pt = inst.proc_time[j][o].at(machine);
        pt_flat[lid]   = pt;
        s_ub_flat[lid] = horizon - pt;
        Op jop{j, o};
        {
            auto it = start_lb_ops.find(jop);
            if (it != start_lb_ops.end()) s_lb_flat[lid] = it->second;
        }
        {
            auto it = start_ub_ops.find(jop);
            if (it != start_ub_ops.end()) s_ub_flat[lid] = it->second;
        }
    }

    // Required ops as flat bool array
    std::vector<bool> req_flat(n_local_ops, false);
    for (auto& req : required_here) {
        int lid = op_to_lid(req.first, req.second);
        if (lid >= 0) req_flat[lid] = true;
    }

    // ── Initial label ─────────────────────────────────────────────────────────
    double alpha_m = 0.0;
    {
        auto it = duals.alpha.find(machine);
        if (it != duals.alpha.end()) alpha_m = it->second;
    }

    BeamLabel init;
    init.score = -alpha_m;
    init.t     = 0;
    init.last_o_placed.assign(n_local_jobs, -1);
    init.placed.assign(n_local_ops, false);

    std::vector<BeamLabel> beam = {init};
    std::vector<BeamLabel> completed;

    for (int depth = 0; depth < n_local_ops && !beam.empty(); ++depth) {
        std::vector<BeamLabel> next_beam;

        for (auto& label : beam) {
            // (new_score, t_new, ji, o, lid)
            // std::vector<std::tuple<double, int, int, int, int>> cands;

            // Candidate fields:
            // (new_score, end_time, job-list index, operation, local id, start_time)
            std::vector<
                std::tuple<
                    double,
                    int,
                    int,
                    int,
                    int,
                    int
                >
            > cands;

            for (int ji = 0; ji < n_local_jobs; ++ji) {
                int j           = jobs_ops[ji].first;
                int last_placed = label.last_o_placed[ji];
                auto& job_ops   = jobs_ops[ji].second;

                for (int qi = 0; qi < (int)job_ops.size(); ++qi) {
                    int o   = job_ops[qi];
                    int lid = local_start[ji] + qi;

                    if (last_placed >= 0 && o <= last_placed) continue;

                    // Për constraint A before B:
                    // ndalo vetëm rastin kur po përpiqemi të vendosim A pasi B është vendosur.
                    // Mos e ndalo kolonën që përmban vetëm B.
                    bool violates_sequence = false;

                    for (int succ_lid : successors_local[lid]) {
                        if (label.placed[succ_lid]) {
                            violates_sequence = true;
                            break;
                        }
                    }

                    if (violates_sequence) {
                        continue;
                    }

                    const int s =
                        std::max(
                            label.t,
                            s_lb_flat[lid]
                        );

                    if (s > s_ub_flat[lid]) {
                        continue;
                    }

                    const int e =
                        s + pt_flat[lid];

                    if (e > horizon) {
                        continue;
                    }

                    const double operation_weight =
                        _weight(
                            j,
                            o,
                            s,
                            pt_flat[lid],
                            inst,
                            duals
                        );

                    const double new_score =
                        label.score + operation_weight;

                    // Store the exact start time used to calculate new_score.
                    cands.push_back({
                        new_score,
                        e,
                        ji,
                        o,
                        lid,
                        s
                    });
                }
            }

            std::sort(
            cands.begin(),
            cands.end()
        );

        int count_for_label = 0;

        for (const auto& [
                new_score,
                end_time,
                ji,
                operation,
                local_id,
                start_time
            ] : cands)
        {
            if (count_for_label++ >= CANDIDATES_PER_LABEL) {
                break;
            }

            BeamLabel new_label = label;

            const int job =
                jobs_ops[ji].first;

            new_label.score =
                new_score;

            new_label.t =
                end_time;

            // Keep operation identity if other code still uses ops.
            new_label.ops.push_back({
                job,
                operation
            });

            // Store exactly the start/end times used by beam expansion.
            new_label.entries.push_back({
                job,
                operation,
                start_time,
                end_time
            });

            new_label.last_o_placed[ji] =
                operation;

            new_label.placed[local_id] =
                true;

            next_beam.push_back(
                std::move(new_label)
            );
        }

            // Keep label as completed (partial schedule is also a candidate column)
            completed.push_back(label);
        }

        // Trim beam to BEAM_WIDTH best by score
        if ((int)next_beam.size() > BEAM_WIDTH) {
            std::partial_sort(next_beam.begin(),
                              next_beam.begin() + BEAM_WIDTH,
                              next_beam.end(),
                              [](const BeamLabel& a, const BeamLabel& b) {
                                  return a.score < b.score;
                              });
            next_beam.resize(BEAM_WIDTH);
        }
        beam = std::move(next_beam);
    }

    for (auto& l : beam) completed.push_back(l);

    // ── Convert labels to columns ─────────────────────────────────────────────
    std::sort(completed.begin(), completed.end(),
              [](const BeamLabel& a, const BeamLabel& b) {
                  return a.score < b.score;
              });

    // std::vector<PricingResult> results;
    // // std::set<std::vector<std::pair<int,int>>> seen_ops;
    // std::set<std::vector<SchedEntry>> seen_columns;
    // int count = 0;

    std::vector<PricingResult> results;

    // Deduplicate using the complete ordered and timed column.
    // SchedEntry = (job, operation, start, end).
    std::set<std::vector<SchedEntry>> seen_columns;

    int count = 0;

    for (auto& lbl : completed) {
        if (lbl.entries.empty()) continue;

        // Required-here check via flat bool array + placed[]
        if (!required_here.empty()) {
            bool all_present = true;
            for (int lid = 0; lid < n_local_ops; ++lid)
                if (req_flat[lid] && !lbl.placed[lid]) { all_present = false; break; }
            if (!all_present) continue;
        }

        Column col = make_column(machine, lbl.entries);

        if (!seen_columns.insert(col.ops).second) {
            continue;
        }

        double rc  = compute_reduced_cost(machine, col, duals, inst);

        if (rc < RC_TOL) {
            results.push_back({col, rc});
            if (++count >= top_k) break;
        }
    }

    // Use the exact times generated during beam expansion.
    // Column col =
    //     make_column(
    //         machine,
    //         lbl.entries
    //     );

    return results;
}

// ===========================================================================
//  Exact DP pricing
// ===========================================================================

// DP state: (ptr_for_each_job, order_mask)
// ptr[i] = minimum op-list index for job i (from jobs_ops)
struct DPKey {
    std::vector<int> ptr;
    int              mask;
    bool operator==(const DPKey& o) const {
        return mask == o.mask && ptr == o.ptr;
    }
};

struct DPKeyHash {
    size_t operator()(const DPKey& k) const {
        size_t h = VecIntHash{}(k.ptr);
        h ^= std::hash<int>{}(k.mask) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

// dp_table[(state, t)] = accumulated_weight (lower = better)
// We store dp_table[state] = vector indexed by t
// But t can be large; use min-value stored per (state,t_quantized).
// For simplicity: dp_table[state][t] with prefix-min.

struct DPBackPtr {
    DPKey prev_state;
    int   prev_t;
    int   ji;  // index in jobs_ops
    int   qi;  // op index within jobs_ops[ji].second
    int   s;   // scheduled start
};

std::vector<PricingResult>
dp_pricing(int machine, const JobsOps& jobs_ops,
           const Instance& inst, const Duals& duals, int horizon,
           const std::set<std::pair<Op, Op>>& order_constraints,
           const std::set<Op>&                required_here,
           const std::map<Op, int>&            start_lb_ops,
           const std::map<Op, int>&            start_ub_ops,
           int top_k,
           double rc_tol,
           bool allow_intentional_idle
        )
{
    if (jobs_ops.empty()) return {};

    int n = (int)jobs_ops.size();

    double alpha_m = 0.0;
    {
        auto it = duals.alpha.find(machine);
        if (it != duals.alpha.end()) alpha_m = it->second;
    }

    // Collect all ops eligible on this machine (from jobs_ops)
    std::set<Op> machine_ops;
    for (auto& [j, ops_list] : jobs_ops)
        for (int o : ops_list) machine_ops.insert({j, o});

    // Build placed_mask tracking.
    // tracked ops = ops in order_constraint pairs (both sides) that are on this machine,
    // plus all required_here ops.
    // op_to_bit[op] = bit index in placed_mask (int32, supports up to 31 tracked ops).
    std::map<Op, int> op_to_bit;
    {
        int bit = 0;
        auto track = [&](const Op& op) {
            if (!op_to_bit.count(op) && machine_ops.count(op)) {
                if (bit < 31) op_to_bit[op] = bit++;
            }
        };

        for (auto& [A, B] : order_constraints)
            if (machine_ops.count(A) && machine_ops.count(B)) { track(A); track(B); }

        std::map<Op, int> forbidden_successor_bits;

        for (auto& [A, B] : order_constraints) {
            if (op_to_bit.count(A) && op_to_bit.count(B)) {
                forbidden_successor_bits[A] |= (1 << op_to_bit.at(B));
            }
        }

        for (auto& req_op : required_here)
            track(req_op);
    }

    // dp_required_bits[op_b] = bitmask of predecessor bits that must all be set
    // in placed_mask before op_b may be placed.
    // std::map<Op, int> dp_required_bits;
    // for (auto& [A, B] : order_constraints)
    //     if (op_to_bit.count(A) && op_to_bit.count(B))
    //         dp_required_bits[B] |= (1 << op_to_bit.at(A));

    std::map<Op, int> forbidden_successor_bits;

    for (const auto& [A, B] : order_constraints) {
        if (op_to_bit.count(A) && op_to_bit.count(B)) {
            forbidden_successor_bits[A] |= (1 << op_to_bit.at(B));
        }
    }

    // required_here_mask: bits for required_here ops; all must be set at terminal states.
    int required_here_mask = 0;
    for (auto& req_op : required_here) {
        auto it = op_to_bit.find(req_op);
        if (it != op_to_bit.end()) required_here_mask |= (1 << it->second);
    }

    // DP: state = (ptr vector, order_bits)
    // ptr[i] = next min op index (in jobs_ops[i].second) to consider
    // When we place qi >= ptr[i] for job i, new_ptr[i] = qi+1

    std::vector<int> init_ptr(n, 0);
    DPKey init_state{init_ptr, 0};

    // dp[(state)] = map<t, accumulated_weight>  (with prefix min)
    using TimeMap = std::unordered_map<int, double>;
    std::unordered_map<DPKey, TimeMap, DPKeyHash> dp;
    dp[init_state][0] = -alpha_m;

    // backpointers: ((state_key, t)) -> DPBackPtr
    // Use a flattened map for efficiency
    using BPKey = std::pair<DPKey, int>;
    struct BPKeyHash {
        size_t operator()(const BPKey& k) const {
            return DPKeyHash{}(k.first) ^ std::hash<int>{}(k.second) * 2654435761u;
        }
    };


    std::unordered_map<BPKey, DPBackPtr, BPKeyHash> bp;

    // Process layers by sum(ptr)
    int max_sum = 0;
    for (auto& [j, ops_list] : jobs_ops) max_sum += (int)ops_list.size();

    for (int total = 0; total < max_sum; ++total) {
        // Collect states with this sum
        std::vector<DPKey> layer_states;
        for (auto& [state, tmap] : dp) {
            int s = 0;
            for (int x : state.ptr) s += x;
            if (s == total) layer_states.push_back(state);
        }

        for (auto& state : layer_states) {
            auto& tmap = dp[state];
            if (tmap.empty()) continue;  
        }

    }

    // Collect terminal states: required_here ops must all have their placed_mask bit set
    std::vector<std::pair<double, std::pair<DPKey, int>>> terminals;
    for (auto& [state, tmap] : dp) {
        if ((state.mask & required_here_mask) != required_here_mask) continue;
        for (auto& [t, w] : tmap)
            terminals.push_back({w, {state, t}});
    }

    std::sort(terminals.begin(), terminals.end(),
        [](const auto& a, const auto& b){ return a.first < b.first; });

    // Backtrack to build columns
    std::vector<PricingResult> results;
    std::set<std::vector<SchedEntry>> seen;

    for (auto& [w, state_t] : terminals) {
        if ((int)results.size() >= top_k) break;

        auto [final_state, final_t] = state_t;

        // Reconstruct column by backtracking
        std::vector<SchedEntry> entries;
        DPKey cur_state = final_state;
        int   cur_t     = final_t;

        int max_bt = max_sum + 5;
        while (max_bt-- > 0) {
            int s_cur = 0;
            for (int x : cur_state.ptr) s_cur += x;
            if (s_cur == 0) break;

            auto it = bp.find({cur_state, cur_t});
            if (it == bp.end()) break;

            auto& bptr = it->second;
            int ji = bptr.ji, qi = bptr.qi, s = bptr.s;
            auto [j, ops_list] = jobs_ops[ji];
            int o  = ops_list[qi];
            int pt = inst.proc_time[j][o].at(machine);
            entries.push_back({j, o, s, s + pt});

            cur_state = bptr.prev_state;
            cur_t     = bptr.prev_t;
        }

        if (entries.empty()) continue;

        Column col = make_column(machine, entries);
        if (seen.count(col.ops)) continue;
        seen.insert(col.ops);

        double rc = compute_reduced_cost(machine, col, duals, inst);
        if (rc < rc_tol)
            results.push_back({col, rc});
    }

    return results;
}

// ===========================================================================
//  Column generation loop
// ===========================================================================

/// Compute reduced cost of a column against given duals.
static double compute_reduced_cost(const Column& col, int machine,
                                    const Duals& d, const Instance& inst)
{
    double rc = 0.0;
    auto it_a = d.alpha.find(machine);
    if (it_a != d.alpha.end()) rc -= it_a->second;

    for (auto& [jop, s] : col.s_vals) {
        auto [j, o] = jop;
        int e  = col.e_vals.at(jop);
        int nj = inst.num_ops[j];

        auto it_b = d.beta.find(jop);
        if (it_b != d.beta.end()) rc -= it_b->second;

        if (o < nj - 1) {
            auto it_t = d.theta.find(jop);
            if (it_t != d.theta.end()) rc -= it_t->second * e;
        }
        if (o > 0) {
            Op prev_op{j, o - 1};
            auto it_t = d.theta.find(prev_op);
            if (it_t != d.theta.end()) rc += it_t->second * s;
        }
        if (o == nj - 1) {
            auto it_m = d.mu.find(j);
            if (it_m != d.mu.end()) rc -= it_m->second * e;
        }
    }
    return rc;
}

/// Build JobsOps for machine m respecting B&P constraints.
static JobsOps build_jobs_ops(
    int machine, const Instance& inst,
    const std::map<Op, int>*           required_machine,
    const std::map<Op, std::set<int>>* forbidden_machines,
    const std::map<Op, std::set<int>>* allowed_machines)
{
    JobsOps result;
    int J = inst.num_jobs;
    for (int j = 0; j < J; ++j) {
        std::vector<int> eligible_ops;
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            Op jop{j, o};

            // Check eligibility
            auto& elig = inst.elig_mach[j][o];
            if (!std::count(elig.begin(), elig.end(), machine)) continue;

            // required_machine: if set and ≠ machine, skip
            if (required_machine) {
                auto it = required_machine->find(jop);
                if (it != required_machine->end() && it->second != machine) continue;
            }
            // forbidden_machines: if machine is forbidden, skip
            if (forbidden_machines) {
                auto it = forbidden_machines->find(jop);
                if (it != forbidden_machines->end() && it->second.count(machine)) continue;
            }
            // allowed_machines: if set, machine must be in it
            if (allowed_machines) {
                auto it = allowed_machines->find(jop);
                if (it != allowed_machines->end() && !it->second.count(machine)) continue;
            }

            eligible_ops.push_back(o);
        }
        if (!eligible_ops.empty())
            result.push_back({j, eligible_ops});
    }
    return result;
}

// CGResult column_generation(
//     const Instance& inst, int horizon, RMP& rmp,
//     int max_iter,
//     const std::map<Op, int>*                         required_machine,
//     const std::map<Op, std::set<int>>*               forbidden_machines,
//     const std::map<Op, std::set<int>>*               allowed_machines,
//     const std::map<int, std::set<std::pair<Op, Op>>>* order_constraints,
//     const std::map<Op, int>*                         start_lb,
//     const std::map<Op, int>*                         start_ub,
//     bool use_beam_first,
//     double cg_time_limit,
//     double incumbent_ub)

CGResult column_generation(
    const Instance& inst,
    int horizon,
    RMP& rmp,
    int max_iter,
    const std::map<Op, int>* required_machine,
    const std::map<Op, std::set<int>>* forbidden_machines,
    const std::map<Op, std::set<int>>* allowed_machines,
    const std::map<int, std::set<std::pair<Op, Op>>>* order_constraints,
    const std::map<Op, int>* start_lb,
    const std::map<Op, int>* start_ub,
    bool use_beam_first,
    double cg_time_limit,
    double incumbent_ub,
    int stall_limit)
{

    CGResult result;

    result.solved_optimally = false;

    result.phase1_feasible = false;
    result.phase1_infeasible_certified = false;

    result.hit_time_limit = false;
    result.hit_iteration_limit = false;

    result.artificial_sum =
        std::numeric_limits<double>::infinity();

    result.lp_bound =
        std::numeric_limits<double>::infinity();

    result.cmax_value =
        std::numeric_limits<double>::infinity();

    result.iterations = 0;
    result.cg_time = 0.0;
    result.beam_frac = 0.0;
    result.pricing_time = 0.0;
    result.master_time = 0.0;


    auto t_cg_start  = std::chrono::steady_clock::now();
    double beam_wall_time = 0.0;    // wall-time for the parallel beam phases

    bool rmp_failed = false;

    // Every RMP solve receives only the remaining time of this CG call.
    // This is required because the outer B&P timer cannot interrupt a
    // CPLEX solve that is already in progress.
    auto remaining_cg_time = [&]() -> double {
        if (cg_time_limit <= 0.0) return 0.0;

        const double elapsed =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_cg_start
            ).count();

        return std::max(0.0, cg_time_limit - elapsed);
    };

    auto solve_rmp_with_remaining_time = [&]() -> double {
        if (cg_time_limit <= 0.0)
            return rmp.solve();

        const double remaining = remaining_cg_time();
        if (remaining <= 0.0)
            return INF;

        return rmp.solve(remaining);
    };

    // Solve initial RMP and account for the complete CPLEX solve time.
    auto t_master_initial = std::chrono::steady_clock::now();
    double obj = solve_rmp_with_remaining_time();
    result.master_time += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_master_initial).count();
    if (!rmp.is_optimal()) {
        if (cg_time_limit > 0.0 &&
            remaining_cg_time() <= 0.0) {
            result.hit_time_limit = true;
        }
        result.lp_bound   = std::numeric_limits<double>::infinity();
        result.cmax_value = std::numeric_limits<double>::infinity();
        result.iterations = 0;
        result.cg_time    = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_cg_start).count();
        result.pricing_time = beam_wall_time;
        return result;
    }

    auto get_artificial_sum = [&]() -> double {
    double sum = 0.0;

    for (int m : inst.all_machines) {
        sum += std::max(
            0.0,
            rmp.get_artificial(m)
        );
    }

    for (int j = 0; j < inst.num_jobs; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            sum += std::max(
                0.0,
                rmp.get_assign_artificial({j, o})
            );
        }
    }

    return sum;
};

auto has_positive_artificials = [&]() -> bool {
    return get_artificial_sum() > LP_POS_TOL;
};

// Initial Phase-I status after the first optimal RMP solve.
result.artificial_sum =
    get_artificial_sum();

result.phase1_feasible =
    result.artificial_sum <= LP_POS_TOL;

Duals prev_raw_duals;  // raw duals from previous iteration (stability centre)
int   total_cols = 0;  // running column count (updated each iter)
for (int m : inst.all_machines) total_cols += rmp.num_columns(m);


static constexpr bool CG_DEBUG = true;

auto count_rmp_cols = [&]() -> int {
    int n = 0;
    for (int m : inst.all_machines)
        n += rmp.num_columns(m);
    return n;
};


// Parallel pricing: always enabled.
const bool use_parallel = true;

// Per-machine pricing result used by parallel helpers below.
struct MachineResult {
    int machine;
    std::vector<PricingResult> found;
};

    // static constexpr int    CG_STALL_LIMIT = 8;
    // static constexpr double CG_DOBJ_TOL    = 1e-3;

    // static constexpr int    CG_STALL_LIMIT = 6;
    static constexpr double CG_DOBJ_TOL    = 1e-3;
    static constexpr double CG_NEAR_PRUNE_MARGIN = 0.10;

    int cg_lb_int_stall_count = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        result.iterations = iter + 1;


        if (cg_time_limit > 0.0 &&
            remaining_cg_time() <= 0.0) {
            result.hit_time_limit = true;
            result.solved_optimally = false;
            result.artificial_sum = get_artificial_sum();
            result.phase1_feasible =
                result.artificial_sum <= LP_POS_TOL;
            break;
        }

        // Get raw (unsmoothed) duals
        Duals raw_duals = rmp.get_duals();


        const double obj_before  = obj;
        const double cmax_before = rmp.get_cmax();
        const double art_before  = get_artificial_sum();
        const int    cols_before = count_rmp_cols();
        

        int beam_returned = 0;
        int beam_negative = 0;
        int beam_kept     = 0;
        int beam_inserted = 0;

        int dp_returned = 0;
        int dp_negative = 0;
        int dp_kept     = 0;
        int dp_inserted = 0;

        double best_rc = INF;
        double worst_kept_rc = INF;

        std::string pricing_mode = "BEAM";

        // Dual smoothing: blend raw duals with previous raw duals
        // Python: alpha = (1-w)*raw + w*prev_raw
        Duals duals = raw_duals;
        if (iter > 0) {
            for (auto& [m, a]  : duals.alpha) {
                auto it = prev_raw_duals.alpha.find(m);
                if (it != prev_raw_duals.alpha.end())
                    a = (1 - STAB_ALPHA) * a + STAB_ALPHA * it->second;
            }
            for (auto& [jo, b] : duals.beta) {
                auto it = prev_raw_duals.beta.find(jo);
                if (it != prev_raw_duals.beta.end())
                    b = (1 - STAB_ALPHA) * b + STAB_ALPHA * it->second;
            }
            for (auto& [jo, t] : duals.theta) {
                auto it = prev_raw_duals.theta.find(jo);
                if (it != prev_raw_duals.theta.end())
                    t = (1 - STAB_ALPHA) * t + STAB_ALPHA * it->second;
            }
            for (auto& [j, mu] : duals.mu) {
                auto it = prev_raw_duals.mu.find(j);
                if (it != prev_raw_duals.mu.end())
                    mu = (1 - STAB_ALPHA) * mu + STAB_ALPHA * it->second;
            }
        }
        // Save raw duals as next iteration's stability centre
        prev_raw_duals = raw_duals;

        // Horizon is fixed for the duration of this CG call (= UB - 1 set by the
        // caller). Do NOT update it from the LP objective here — that would make
        // pricing inconsistent with the global B&P horizon.

        // Price each machine using smoothed duals; verify against raw duals
        int  cols_added_this_iter = 0;
        bool beam_found_any = false;

        // ── Phase 1: beam search (parallel when use_parallel) ────────────────
        auto run_beam = [&](int machine) -> MachineResult {
            JobsOps jobs_ops = build_jobs_ops(machine, inst,
                                               required_machine, forbidden_machines,
                                               allowed_machines);
            if (jobs_ops.empty()) {
                return {
                    machine,
                    {}
                };
            }

            std::set<std::pair<Op,Op>> m_oc;
            if (order_constraints) {
                auto it = order_constraints->find(machine);
                if (it != order_constraints->end()) m_oc = it->second;
            }

            std::set<Op> req_here;
            if (required_machine)
                for (auto& [jop, m2] : *required_machine)
                    if (m2 == machine) req_here.insert(jop);
            std::vector<PricingResult> found;
            if (use_beam_first)
                // found = beam_search_pricing(machine, jobs_ops, inst, duals, horizon,
                //                             req_here,
                //                             start_lb ? *start_lb : std::map<Op,int>{},
                //                             start_ub ? *start_ub : std::map<Op,int>{});

                found = beam_search_pricing(
                    machine,
                    jobs_ops,
                    inst,
                    duals,
                    horizon,
                    m_oc,
                    req_here,
                    start_lb ? *start_lb : std::map<Op,int>{},
                    start_ub ? *start_ub : std::map<Op,int>{}
                );
            return {
                machine,
                std::move(found)
            };
        };

        std::vector<MachineResult> beam_res;
        {
            auto t_b0 = std::chrono::steady_clock::now();
            if (use_parallel) {
                std::vector<std::future<MachineResult>> futs;
                futs.reserve(inst.all_machines.size());
                for (int machine : inst.all_machines)
                    futs.push_back(std::async(std::launch::async, run_beam, machine));
                for (auto& f : futs) {
                    MachineResult mr = f.get();
                    beam_res.push_back(std::move(mr));
                }
            } else {
                for (int machine : inst.all_machines) {
                    MachineResult mr = run_beam(machine);
                    beam_res.push_back(std::move(mr));
                }
            }
            beam_wall_time += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_b0
            ).count();
        }


        // Verify beam results against raw duals; keep the TOP_K_GLOBAL_COLS most negative.
        {
            struct Cand {
                int machine;
                Column col;
                double rc;
            };

            std::vector<Cand> cands;

            for (auto& mr : beam_res) {
                beam_returned += static_cast<int>(mr.found.size());

                for (auto& pr : mr.found) {
                    double rc =
                        compute_reduced_cost(
                            pr.col,
                            mr.machine,
                            raw_duals,
                            inst
                        );

                    best_rc = std::min(best_rc, rc);

                    if (rc < RC_TOL) {
                        cands.push_back({
                            mr.machine,
                            pr.col,
                            rc
                        });
                    }
                }
            }

            beam_negative = static_cast<int>(cands.size());

            if (!cands.empty()) {
                std::sort(
                    cands.begin(),
                    cands.end(),
                    [](const Cand& a, const Cand& b) {
                        return a.rc < b.rc;
                    }
                );

                static constexpr int MAX_PER_MACHINE_PER_ITER = 3; //3

                std::map<int, int> added_per_machine;
                std::vector<Cand> selected;
                selected.reserve(TOP_K_GLOBAL_COLS);

                // cands is already sorted by reduced cost, most negative first.
                for (const Cand& candidate : cands) {
                    const int machine = candidate.machine;

                    if (added_per_machine[machine] >= MAX_PER_MACHINE_PER_ITER) {
                        continue;
                    }

                    selected.push_back(candidate);
                    ++added_per_machine[machine];

                    if (static_cast<int>(selected.size()) >= TOP_K_GLOBAL_COLS) {
                        break;
                    }
                }

                // Fallback: if the per-machine cap selected too few columns,
                // fill the remaining slots with the best remaining candidates.
                if (selected.empty()) {
                    const int keep =
                        std::min(
                            static_cast<int>(cands.size()),
                            TOP_K_GLOBAL_COLS
                        );

                    selected.insert(
                        selected.end(),
                        cands.begin(),
                        cands.begin() + keep
                    );
                }

                beam_kept = static_cast<int>(selected.size());

                if (!selected.empty()) {
                    worst_kept_rc = selected.back().rc;
                }

                std::map<int, std::vector<Column>> by_machine;

                for (const Cand& candidate : selected) {
                    by_machine[candidate.machine].push_back(
                        candidate.col
                    );
                }

                for (auto& [machine, columns] : by_machine) {
                    int added =
                        rmp.add_columns(
                            machine,
                            columns
                        );

                    beam_inserted += added;
                    cols_added_this_iter += added;
                }

                beam_found_any = beam_inserted > 0;
            }
        }


        if (!beam_found_any) {
            pricing_mode = "DP";

            // Beam found nothing — run exact DP.
            // Collect all DP candidates globally; keep TOP_K_GLOBAL_COLS most negative.
            struct DPCand {
                int machine;
                Column col;
                double rc;
            };

            std::vector<DPCand> dp_cands;

            for (int machine : inst.all_machines) {
                JobsOps jobs_ops =
                    build_jobs_ops(
                        machine,
                        inst,
                        required_machine,
                        forbidden_machines,
                        allowed_machines
                    );

                if (jobs_ops.empty()) {
                    continue;
                }

                std::set<Op> req_here;

                if (required_machine) {
                    for (const auto& [jop, m2] : *required_machine) {
                        if (m2 == machine) {
                            req_here.insert(jop);
                        }
                    }
                }

                std::set<std::pair<Op, Op>> m_oc;

                if (order_constraints) {
                    auto it = order_constraints->find(machine);
                    if (it != order_constraints->end()) {
                        m_oc = it->second;
                    }
                }

                auto dp_results = dp_pricing(
                    machine,
                    jobs_ops,
                    inst,
                    raw_duals,
                    horizon,
                    m_oc,
                    req_here,
                    start_lb ? *start_lb : std::map<Op, int>{},
                    start_ub ? *start_ub : std::map<Op, int>{},
                    TOP_K_PER_MACHINE,
                    RC_TOL,
                    true
                );

                dp_returned += static_cast<int>(dp_results.size());

                for (auto& pr : dp_results) {
                    double rc =
                        compute_reduced_cost(
                            pr.col,
                            machine,
                            raw_duals,
                            inst
                        );

                    best_rc = std::min(best_rc, rc);

                    if (rc < RC_TOL) {
                        dp_cands.push_back({
                            machine,
                            pr.col,
                            rc
                        });
                    }
                }
            }

            dp_negative = static_cast<int>(dp_cands.size());

            std::sort(
                dp_cands.begin(),
                dp_cands.end(),
                [](const DPCand& a, const DPCand& b) {
                    return a.rc < b.rc;
                }
            );

            if (static_cast<int>(dp_cands.size()) > TOP_K_GLOBAL_COLS) {
                dp_cands.resize(TOP_K_GLOBAL_COLS);
            }

            dp_kept = static_cast<int>(dp_cands.size());

            if (!dp_cands.empty()) {
                worst_kept_rc = dp_cands.back().rc;
            }

            std::map<int, std::vector<Column>> by_machine;

            for (const auto& candidate : dp_cands) {
                by_machine[candidate.machine].push_back(
                    candidate.col
                );
            }

            for (const auto& [machine, columns] : by_machine) {
                int added =
                    rmp.add_columns(
                        machine,
                        columns
                    );

                dp_inserted += added;
                cols_added_this_iter += added;
            }

            if (dp_inserted == 0) {
                result.artificial_sum =
                    get_artificial_sum();

                result.phase1_feasible =
                    result.artificial_sum <= LP_POS_TOL;

                if (CG_DEBUG) {
                    std::cout
                        << "  [CGDBG-END] it=" << iter
                        << " mode=DP"
                        << " obj=" << std::fixed << std::setprecision(4)
                        << obj_before
                        << " Cmax=" << cmax_before
                        << " art=" << art_before
                        << " cols=" << cols_before
                        << " dp_ret/neg/kept/ins="
                        << dp_returned << "/"
                        << dp_negative << "/"
                        << dp_kept << "/"
                        << dp_inserted
                        << " best_rc=" << best_rc
                        << " reason=no_new_columns"
                        << "\n";
                }

                if (result.artificial_sum > LP_POS_TOL) {
                    result.phase1_feasible = false;
                    result.phase1_infeasible_certified = true;
                    result.solved_optimally = true;

                    std::cout
                        << "  [CG] Exact Phase-I infeasibility certificate:"
                        << " artificial_sum="
                        << result.artificial_sum
                        << "\n";
                } else {
                    result.phase1_feasible = true;
                    result.phase1_infeasible_certified = false;
                    result.solved_optimally = true;

                    std::cout
                        << "  [CG] Exact pricing converged";
                        // << " artificial_sum="
                        // << result.artificial_sum
                        // << "\n";
                }

                break;
            }
        }

        // Re-solve RMP with any newly added columns
        // double new_obj = rmp.solve();
        // if (!rmp.is_optimal()) break;

        // double new_obj = rmp.solve();

        auto t_master0 = std::chrono::steady_clock::now();

        double new_obj = solve_rmp_with_remaining_time();

        double master_sec =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_master0
            ).count();
        result.master_time += master_sec;

        if (!rmp.is_optimal()) {
            if (cg_time_limit > 0.0 &&
                remaining_cg_time() <= 0.0) {
                result.hit_time_limit = true;
            }
            rmp_failed = true;
            result.solved_optimally = false;
            break;
        }


        const double obj_after  = new_obj;
        const double cmax_after = rmp.get_cmax();
        const double art_after  = get_artificial_sum();
        const int    cols_after = count_rmp_cols();

        if (CG_DEBUG) {
            std::cout
                << "  [CGDBG] it=" << iter
                // << " mode=" << pricing_mode
                // << " obj=" << std::fixed << std::setprecision(4)
                // << obj_before << "->" << obj_after
                // << " dObj=" << (obj_before - obj_after)
                << " LP_Cmax=" << cmax_before << "->" << cmax_after
                << " dCmax=" << (cmax_before - cmax_after)
                // << " art=" << art_before << "->" << art_after
                << " cols=" << cols_before << "->" << cols_after
                << " added=" << cols_added_this_iter
                // << " beam_ret/neg/kept/ins="
                // << beam_returned << "/"
                // << beam_negative << "/"
                // << beam_kept << "/"
                // << beam_inserted
                // << " dp_ret/neg/kept/ins="
                // << dp_returned << "/"
                // << dp_negative << "/"
                // << dp_kept << "/"
                // << dp_inserted
                // << " best_rc=" << best_rc
                // << " worst_kept_rc=" << worst_kept_rc
                // << " master_s=" << master_sec
                << "\n";
        }


        // ── Tailing-off stop: LB_int does not improve ─────────────────────────────
        // This is NOT an exact convergence certificate.
        // It only stops CG early when:
        //   1) Phase-I is feasible,
        //   2) integer LB does not change,
        //   3) LP_Cmax improvement is tiny,
        //   4) the node is NOT close to being pruned by bound.
        {
            const bool phase1_ok =
                art_after <= LP_POS_TOL;

            if (phase1_ok) {
                const double lb_int_before =
                    integer_lb_from_cmax(cmax_before);

                const double lb_int_after =
                    integer_lb_from_cmax(cmax_after);

                const double dcmax =
                    std::abs(cmax_before - cmax_after);

                if (lb_int_before == lb_int_after &&
                    dcmax < CG_DOBJ_TOL)
                {
                    ++cg_lb_int_stall_count;
                } else {
                    cg_lb_int_stall_count = 0;
                }

                const bool near_prune =
                    std::isfinite(incumbent_ub) &&
                    cmax_after > incumbent_ub - 1.0 - CG_NEAR_PRUNE_MARGIN;

                if (!near_prune &&
                    cg_lb_int_stall_count >= stall_limit)
                {
                    result.solved_optimally = true;
                    result.hit_time_limit = false;

                    result.artificial_sum =
                        get_artificial_sum();

                    result.phase1_feasible =
                        result.artificial_sum <= LP_POS_TOL;

                    result.lp_bound =
                        rmp.is_optimal()
                            ? rmp.objective()
                            : INF;

                    result.cmax_value =
                        rmp.is_optimal()
                            ? rmp.get_cmax()
                            : INF;

                    if (CG_DEBUG) {
                        std::cout
                            << "  [CGDBG-END] it=" << iter
                            << " reason=tailing_off_not_near_prune"
                            << " LP_Cmax=" << cmax_after
                            << " LB_int=" << lb_int_after
                            << " UB=" << incumbent_ub
                            << " UB_minus_1=" << incumbent_ub - 1.0
                            << " near_prune=" << (near_prune ? 1 : 0)
                            << " stall=" << cg_lb_int_stall_count
                            << " dCmax=" << dcmax
                            << "\n";
                    }

                    break;
                }
            }
        }

        total_cols += cols_added_this_iter;
        obj = new_obj;

        result.artificial_sum =
            get_artificial_sum();

        result.phase1_feasible =
            result.artificial_sum <= LP_POS_TOL;
    }

    // The loop exhausted max_iter without convergence, timeout,
    // or an RMP failure.
    if (!result.solved_optimally &&
        !result.hit_time_limit &&
        !rmp_failed &&
        max_iter > 0 &&
        result.iterations >= max_iter)
    {
        result.hit_iteration_limit = true;
    }

    // Final solve to get accurate objective / solution values
    // rmp.solve();
    // result.lp_bound   = rmp.is_optimal() ? rmp.objective() : obj;
    // result.cmax_value = rmp.get_cmax();
    // result.y_values   = rmp.get_y_values();

    auto t_master_final = std::chrono::steady_clock::now();
    const double final_obj =
    solve_rmp_with_remaining_time();
    result.master_time += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_master_final).count();

    if (rmp.is_optimal()) {
        result.lp_bound =
            rmp.objective();

        result.cmax_value =
            rmp.get_cmax();

        result.y_values =
            rmp.get_y_values();

        // Always return the real final Phase-I state.
        result.artificial_sum =
            get_artificial_sum();

        result.phase1_feasible =
            result.artificial_sum <= LP_POS_TOL;

        // A feasible final RMP invalidates any accidental stale certificate.
        if (result.phase1_feasible) {
            result.phase1_infeasible_certified = false;
        }
    } else {
        result.lp_bound =
            std::isfinite(obj)
                ? obj
                : final_obj;

        result.cmax_value =
            std::numeric_limits<double>::infinity();

        result.y_values.clear();

        result.solved_optimally = false;
    }

    // Collect column lists
    for (int m : inst.all_machines)
        result.columns[m] = rmp.get_columns(m);

    result.cg_time   = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_cg_start).count();
    result.pricing_time = beam_wall_time;
    result.beam_frac = (result.cg_time > 0.0)
        ? beam_wall_time / result.cg_time
        : 0.0;

    return result;
}
