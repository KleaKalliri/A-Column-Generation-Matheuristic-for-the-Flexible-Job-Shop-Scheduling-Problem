#pragma once
// =============================================================================
// branch_and_price.h  –  Branch-and-Price for FJSSP
//
// Mirrors branch_and_price_machine_based.py.
// Implements:
//   BranchAndPriceNode  –  a single B&P tree node with all constraint data
//   BranchAndPriceSolver – the main solver with best-first search
// =============================================================================

#include "types.h"
#include "column_generation.h"
// #include "mlp_policy.h"
#include <chrono>
#include <functional>
#include <map>
// #include <memory>
#include <optional>
#include <queue>
#include <set>
#include <string>
#include <tuple>
#include <vector>

// ---------------------------------------------------------------------------
// B&P node
// ---------------------------------------------------------------------------

enum class NodeStatus { OPEN, PRUNED, INTEGER, INFEASIBLE };

struct BranchDecision {
    enum class Type {
        OMA,
        SEQUENCE
    };

    Type type = Type::SEQUENCE;

    // OMA branching:
    // child 0 -> fix operation on machine
    // child 1 -> forbid operation on machine
    Op  oma_op{-1, -1};
    int oma_machine = -1;

    // Sequence branching:
    // child 0 -> A before B on order_machine
    // child 1 -> B before A on order_machine
    int order_machine = -1;
    std::pair<Op, Op> order_pair;

    static BranchDecision oma(Op op, int machine) {
        BranchDecision d;
        d.type = Type::OMA;
        d.oma_op = op;
        d.oma_machine = machine;
        return d;
    }

    static BranchDecision sequence(int machine, Op A, Op B) {
        BranchDecision d;
        d.type = Type::SEQUENCE;
        d.order_machine = machine;
        d.order_pair = {A, B};
        return d;
    }
};

/// A node in the B&P tree.
struct BranchAndPriceNode {
    int        node_id  = 0;
    int        depth    = 0;
    int        parent_id= -1;

    double     lower_bound    = 0.0;
    double     cp_lower_bound = 0.0;
    NodeStatus status         = NodeStatus::OPEN;

    std::optional<BranchDecision> branch_decision;

    // ---- B&P constraint data ----
    /// Fixed machine assignment: (j,o) -> machine
    std::map<Op, int>        required_machine;
    /// Forbidden machines: (j,o) -> set of forbidden machines
    std::map<Op, std::set<int>> forbidden_machines;
    /// Allowed machines: (j,o) -> set of allowed machines
    std::map<Op, std::set<int>> allowed_machines;
    /// Order constraints per machine: machine -> {((j1,o1),(j2,o2)): A before B}
    // std::map<int, std::set<std::pair<Op,Op>>> order_constraints;

    /// Banned column signatures (fallback column branching).
    /// A signature is the ordered list of (job,op) pairs in machine-scheduling order.
    std::set<std::vector<std::pair<int,int>>> banned_column_signatures;

    /// Fixed machine columns (fallback column branching).
    /// machine -> the ONE column signature allowed for that machine; all others banned.
    std::map<int, std::vector<std::pair<int,int>>> fixed_machine_columns;

    std::map<int, std::set<std::pair<Op, Op>>> order_constraints;

    // ---- Time window bounds ----
    std::map<Op, int> start_lb;  // earliest start
    std::map<Op, int> start_ub;  // latest start

    // ---- LP solution state ----
    std::map<int, std::vector<Column>> cols;  // machine -> columns
    std::map<std::pair<int,int>, double> y_values;  // (m,p) -> weight
    double lp_cmax = 0.0;

    // ---- Stall pruning counter ----
    // Counts how many consecutive ancestor nodes (including this one) produced
    // a CP-SAT result equal to the current UB without improving it.
    // Pruned when this reaches 3.
    int consecutive_ub_hits = 0;

    // Comparison for priority queue (min-heap by lower_bound)
    bool operator>(const BranchAndPriceNode& o) const {
        return lower_bound > o.lower_bound;
    }
};

// ---------------------------------------------------------------------------
// CP-SAT interface
// ---------------------------------------------------------------------------

struct CPSatResult {
    double   makespan;
    Schedule schedule;
    bool     found = false;
};

/// Solve FJSSP relaxation via CP-SAT with node constraints (optimization mode).
/// Optional hint_schedule provides an initial warm-start assignment for CP-SAT
/// (e.g. from LP-rounding or an RL policy). Keys are (j,o), values are (m,s,e).
CPSatResult cp_sat_solve(const Instance& inst,
                          const BranchAndPriceNode& node,
                          double time_limit_sec,
                          double ub = std::numeric_limits<double>::infinity(),
                          const Schedule* hint_schedule = nullptr);

/// Solve FJSSP as a pure feasibility problem: find any schedule with
/// Cmax <= target_cmax.  No objective is minimized; CP-SAT stops at the
/// first feasible assignment it finds.  Uses node branching constraints and
/// an optional warm-start hint_schedule derived from CG column masses.
CPSatResult cp_sat_feasibility_solve(const Instance& inst,
                                      const BranchAndPriceNode& node,
                                      double time_limit_sec,
                                      int target_cmax,
                                      const Schedule* hint_schedule = nullptr);

// ---------------------------------------------------------------------------
// Propagation
// ---------------------------------------------------------------------------

/// Applies constraint propagation rules R1-R7 to a node.
/// Returns false if infeasibility is detected.
bool propagate_node(BranchAndPriceNode& node, const Instance& inst, double best_ub = 1e30);

// ---------------------------------------------------------------------------
// Column feasibility check
// ---------------------------------------------------------------------------

bool column_is_feasible_at_node(const Column& col, int machine,
                                  const BranchAndPriceNode& node,
                                  const Instance& inst);

// ---------------------------------------------------------------------------
// Main solver
// ---------------------------------------------------------------------------

struct BranchAndPriceResult {
    double   best_bound;
    double   best_obj;
    Schedule best_schedule;
    bool     proven_optimal;
    int      nodes_explored;
    int      integer_nodes;
    int      pruned_bound;
    int      pruned_infeasible;
    double   elapsed_seconds;
    double   total_pricing_time = 0.0; ///< beam-search wall-time summed over all nodes
    double   total_master_time  = 0.0; ///< CPLEX RMP solve time summed over all nodes
    int      prec_graph_calls    = 0;
    double   mean_prec_graph_ms  = 0.0;  ///< mean wall-time per precedence-graph extraction
    int      total_cols_generated = 0;   ///< total columns generated across all B&P nodes
    double   avg_beam_frac        = 0.0; ///< average pricing/master time split (beam fraction)
};

class BranchAndPriceSolver {
public:
    explicit BranchAndPriceSolver(const Instance& inst);

    /// Call before solve() to supply a known UB (e.g. from CP-SAT run in main).
    /// solve() will skip its internal CP-SAT call if this has been set.
    void set_initial_ub(double makespan, const Schedule& sched);

    /// Enable or disable the CP-SAT improvement call made after column
    /// generation at each B&P node. The initial UB call is unaffected.
    void set_node_cp_sat_enabled(bool enabled) { node_cp_sat_enabled_ = enabled; }

    int    get_horizon() const { return horizon_; }
    double get_best_ub() const { return best_ub_; }

    /// Enable RL training-data logging to a CSV file.
    /// Each node's CP-SAT call logs: state features, LP-rounded action, reward.
    // void set_rl_log(const std::string& path);

    /// Load a trained RL policy from a weights file (written by train_rl_policy.py).
    /// When loaded, cp_sat_solve hints are generated by the MLP instead of
    /// greedy LP-argmax.
    // void set_rl_policy(const std::string& path);

    BranchAndPriceResult solve();

    // Optional: callbacks
    std::function<void(int node_id, double lb, double ub)> on_node_cb;

private:
    const Instance& inst_;
    int             horizon_;

    double best_ub_  = std::numeric_limits<double>::infinity();
    Schedule best_schedule_;
    bool node_cp_sat_enabled_ = true;

    std::vector<Schedule> incumbent_schedules_;

    int next_node_id_       = 0;
    int integer_nodes_      = 0;
    int pruned_bound_       = 0;
    int pruned_infeasible_  = 0;
    int no_branch_nodes_    = 0;  ///< nodes where select_branch found no OMA or ORDER candidate
    int    prec_graph_calls_      = 0;    ///< total precedence-graph extraction invocations
    double prec_graph_total_time_ = 0.0;  ///< accumulated wall-time of all extractions (seconds)
    int    last_cg_iters_  = 0;
    double last_cg_time_   = 0.0;
    double last_beam_frac_ = 0.0;
    double last_lp_cmax_   = 0.0;
    bool heuristic_cg_prune_occurred_ = false;
    int    total_cols_generated_ = 0;   ///< cumulative columns generated across all nodes
    double total_pricing_time_   = 0.0; ///< cumulative parallel beam wall-time (s)
    double total_master_time_    = 0.0; ///< cumulative CPLEX RMP solve time (s)
    double sum_beam_frac_        = 0.0; ///< sum of per-node beam_frac for average computation
    int    cg_node_count_        = 0;   ///< number of nodes where CG ran (denominator for avg)
    std::chrono::steady_clock::time_point solve_start_;

    /// RL data logger — null when logging is disabled.
    // std::shared_ptr<class RLDataLogger> rl_logger_;

    /// Trained RL policy for machine-assignment scoring — null when not loaded.
    // std::shared_ptr<MLPPolicy> mlp_policy_;

    /// Process one B&P node: run CG, check integrality, set up branching.
    /// Returns child nodes (or empty if pruned/integer).
    std::vector<BranchAndPriceNode>
    process_node(BranchAndPriceNode& node);

    /// Select branching variable.
    std::optional<BranchDecision>
    select_branch(const BranchAndPriceNode& node,
                  const std::map<std::pair<int,int>, double>& y_values,
                  const std::map<int, std::vector<Column>>& cols);

    /// Check LP solution integrality.
    bool is_integral(const BranchAndPriceNode& node,
                     const std::map<std::pair<int,int>, double>& y_values,
                     const std::map<int, std::vector<Column>>& cols);

    /// Extract schedule from integral LP solution.
    Schedule extract_schedule(const BranchAndPriceNode& node,
                               const std::map<std::pair<int,int>, double>& y_values,
                               const std::map<int, std::vector<Column>>& cols);

    /// Update upper bound; returns true if improved.
    bool update_ub(double makespan, const Schedule& sched);
};
