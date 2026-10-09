#pragma once
// =============================================================================
// column_generation.h  –  Machine-based column generation for FJSSP
//
// Mirrors column_generation_machine_based.py.
// Provides:
//   - construction heuristics (greedy EFT, dispatch, random)
//   - RMP class (LP via IBM CPLEX Concert Technology)
//   - beam-search pricing
//   - exact DP pricing (with order-constraint support)
//   - main column_generation() loop
// =============================================================================

#include "types.h"
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

// CPLEX Concert Technology (IL_STD must be defined before including)
#ifndef IL_STD
#  define IL_STD
#endif
#include <ilcplex/ilocplex.h>

// ---------------------------------------------------------------------------
// Heuristic: compute a horizon
// ---------------------------------------------------------------------------

/// Returns sum over all jobs of max processing time across ops and eligible
/// machines — used as a generous upper bound for the time horizon.
int compute_horizon(const Instance& inst);

// ---------------------------------------------------------------------------
// Construction heuristics
// ---------------------------------------------------------------------------

/// Build one feasible schedule column per machine using greedy EFT dispatch.
/// `rand_seed` controls randomization; -1 means deterministic EFT.
std::map<int, Column>
build_initial_heuristic(const Instance& inst, int horizon, int rand_seed = -1);

/// Build multiple sets of initial columns (greedy variants + randomized
/// construction + optional CP-SAT-derived schedules),
/// deduplicate, and return for each machine a list of unique columns.
/// cpsat_schedules: zero or more full schedules produced by CP-SAT; each is
///   decomposed into per-machine columns and added to the pool.
std::map<int, std::vector<Column>>
build_initial_columns(const Instance& inst, int horizon,
                      const std::vector<Schedule>& cpsat_schedules = {});

// ---------------------------------------------------------------------------
// RMP  –  Restricted Master Problem (LP via GLOP)
// ---------------------------------------------------------------------------

/// Per-machine column list with LP variable handles
struct MachineColumns {
    std::vector<Column>      cols;
    std::vector<IloNumVar>   vars;
};

class RMP {
public:
    explicit RMP(const Instance& inst, int horizon);
    ~RMP();

    // Add columns (incremental — does not rebuild existing constraints)
    // void add_columns(int machine, const std::vector<Column>& new_cols);

    // Add columns incrementally and return the number of columns
    // that were actually inserted after duplicate removal.
    int add_columns(
        int machine,
        const std::vector<Column>& new_cols
    );

    // Solve and return objective (or +inf on failure)
    // Solve the current RMP. If time_limit_sec > 0, CPLEX is given that
    // remaining wall-clock limit for this individual solve.
    double solve(double time_limit_sec = 0.0);

    // Solution
    bool is_optimal() const;
    bool is_infeasible() const;
    double objective() const;

    // Dual values from last solve
    Duals get_duals() const;

    // Primal solution for column-based decomposition
    /// y[(m,p)] = weight of column p on machine m
    std::map<std::pair<int, int>, double> get_y_values() const;

    /// C_max value from last solve
    double get_cmax() const;

    // Access column lists
    const std::vector<Column>& get_columns(int machine) const;
    int num_columns(int machine) const;

    // Artificial variable value for machine m (>0 => infeasible assignment)
    double get_artificial(int machine) const;
    // Artificial variable value for op (j,o) assignment (>0 => infeasible)
    double get_assign_artificial(Op jop) const;

    // Remove columns not feasible under B&P node constraints
    // (rebuilds the LP; call only when necessary)
    void rebuild_with_feasible_cols(const std::map<int, std::vector<Column>>& cols_by_machine);

    // Solve the RMP as a binary MIP (y ∈ {0,1}) to obtain a primal integer solution.
    // Returns the best feasible solution found within time_limit_sec (0 = unlimited).
    // The y_values in MIPResult use the same (machine, col_index) keys as get_y_values().
    // Restores continuous LP state after the MIP solve.
    struct MIPResult {
        bool   found    = false;
        double makespan = std::numeric_limits<double>::infinity();
        std::map<std::pair<int,int>, double> y_values;
    };
    MIPResult solve_mip(double time_limit_sec = 0.0,
                        double ub_cutoff = std::numeric_limits<double>::infinity());

private:
    const Instance& inst_;
    int horizon_;

    // CPLEX Concert Technology objects
    IloEnv       env_;
    IloModel     model_;
    IloCplex     cplex_;
    IloObjective obj_;

    // Variables
    IloNumVar             cmax_var_;
    std::map<int, IloNumVar>    art_vars_;
    std::map<Op, IloNumVar>     assign_art_vars_;
    std::map<int, MachineColumns> machine_cols_;

    // Constraints
    std::map<int, IloRange>  mach_cov_ctr_;   // (2)
    std::map<Op, IloRange>   op_assign_ctr_;  // (3)
    std::map<Op, IloRange>   prec_ctr_;       // (4)
    std::map<int, IloRange>  mkspan_ctr_;     // (5)

    bool optimal_   = false;
    bool infeasible_= false;

    void build_base_constraints();
    void add_column_to_lp(int machine, const Column& col, IloNumVar var);
};

// ---------------------------------------------------------------------------
// Pricing oracle  –  generates improving columns (reduced cost < RC_TOL)
// ---------------------------------------------------------------------------

/// Context passed to pricing oracles: which ops are eligible per job on machine m.
/// jobs_ops[i] = (job_idx, list_of_op_indices)
using JobsOps = std::vector<std::pair<int, std::vector<int>>>;

/// A pricing result: new column + its reduced cost
struct PricingResult {
    Column col;
    double rc;
};

/// Beam-search pricing: returns up to TOP_K_PER_MACHINE columns with rc < RC_TOL.
std::vector<PricingResult>
beam_search_pricing(int machine, const JobsOps& jobs_ops,
                    const Instance& inst, const Duals& duals,
                    int horizon,
                    const std::set<std::pair<Op, Op>>& order_constraints = {},
                    const std::set<Op>&                required_here      = {},
                    const std::map<Op, int>&            start_lb_ops      = {},
                    const std::map<Op, int>&            start_ub_ops      = {},
                    int top_k = TOP_K_PER_MACHINE);

/// Exact DP pricing (with optional order constraints).
/// order_constraints: pairs ((j1,o1),(j2,o2)) meaning end(j1,o1) <= start(j2,o2)
/// required_here: set of (j,o) that MUST be included in the column
/// start_lb/ub: time-window bounds on start times
/// rc_tol: reduced-cost threshold (default RC_TOL); pass +inf for repair pricing
std::vector<PricingResult>
dp_pricing(int machine, const JobsOps& jobs_ops,
           const Instance& inst, const Duals& duals, int horizon,
           const std::set<std::pair<Op, Op>>& order_constraints = {},
           const std::set<Op>&                required_here     = {},
           const std::map<Op, int>&            start_lb_ops     = {},
           const std::map<Op, int>&            start_ub_ops     = {},
           int top_k = TOP_K_PER_MACHINE,
           double rc_tol = RC_TOL,
           bool allow_intentional_idle = false
        );

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
    const std::map<Op, int>& start_ub
);

/// Generates columns for one machine that must contain every operation
/// in mandatory_ops. Used after OMA branching when several operations
/// are fixed to the same machine.
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
    int top_k = 5
);

// ---------------------------------------------------------------------------
// Reduced-cost computation (verification)
// ---------------------------------------------------------------------------

double compute_reduced_cost(int machine, const Column& col, const Duals& duals,
                             const Instance& inst);

// ---------------------------------------------------------------------------
// Main column generation loop
// ---------------------------------------------------------------------------

struct CGResult {
    double               lp_bound;
    double               cmax_value;
    std::map<std::pair<int,int>, double> y_values;  // (machine, col_idx) -> weight
    std::map<int, std::vector<Column>>  columns;    // machine -> column list
    bool                 solved_optimally;

    // New Phase-I information
    bool phase1_feasible = false;
    bool phase1_infeasible_certified = false;

    bool hit_time_limit = false;
    bool hit_iteration_limit = false;

    double artificial_sum =
        std::numeric_limits<double>::infinity();

    int                  iterations;
    double               cg_time;    // wall-clock seconds for entire CG call
    double               beam_frac;  // wall-clock fraction of CG spent in beam phases
    double               pricing_time = 0.0; // wall-clock time of beam phases
    double               master_time  = 0.0; // accumulated CPLEX RMP solve time
};



CGResult column_generation(
    const Instance& inst,
    int             horizon,
    RMP&            rmp,
    int             max_iter = CG_MAX_ITER,
    const std::map<Op, int>*                           required_machine  = nullptr,
    const std::map<Op, std::set<int>>*                 forbidden_machines= nullptr,
    const std::map<Op, std::set<int>>*                 allowed_machines  = nullptr,
    const std::map<int, std::set<std::pair<Op, Op>>>* order_constraints = nullptr,
    const std::map<Op, int>*                           start_lb          = nullptr,
    const std::map<Op, int>*                           start_ub          = nullptr,
    bool                                               use_beam_first    = true,
    double                                             cg_time_limit     = 0.0,
    double                                             incumbent_ub      = std::numeric_limits<double>::infinity(),
    int                                                stall_limit       = CG_STALL_LIMIT
);
