#pragma once
// =============================================================================
// types.h  –  Core data structures for machine-based FJSSP Column Generation
//             and Branch-and-Price.
//
// Mirrors the Python types used in:
//   column_generation_machine_based.py  and  branch_and_price_machine_based.py
// =============================================================================

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// Aliases
// ---------------------------------------------------------------------------

/// Operation identity: (job_index, op_index)
using Op = std::pair<int, int>;

/// A single scheduled operation: (job, op_index, start, end)
using SchedEntry = std::tuple<int, int, int, int>;

/// Full schedule: op -> (machine, start, end)
using Schedule = std::map<Op, std::tuple<int, int, int>>;

// ---------------------------------------------------------------------------
// Column  –  a schedule for one machine
// ---------------------------------------------------------------------------

struct Column {
    int machine = -1;

    /// Operations assigned to this machine, sorted by start time.
    /// Each entry: (job, op_index, start, end)
    std::vector<SchedEntry> ops;

    /// Set of (job, op_index) pairs covered by this column.
    std::set<Op> a;

    /// Start times: a[(j,o)] -> start
    std::map<Op, int> s_vals;

    /// End times: a[(j,o)] -> end
    std::map<Op, int> e_vals;

    bool operator==(const Column& o) const {
        return machine == o.machine && ops == o.ops;
    }
};

/// Build a column from a list of (j, o, s, e) tuples on the given machine.
inline Column make_column(int machine, std::vector<SchedEntry> raw_ops) {
    std::sort(raw_ops.begin(), raw_ops.end(),
              [](const SchedEntry& a, const SchedEntry& b) {
                  return std::get<2>(a) < std::get<2>(b);
              });
    Column col;
    col.machine = machine;
    col.ops     = raw_ops;
    for (auto& [j, o, s, e] : raw_ops) {
        col.a.insert({j, o});
        col.s_vals[{j, o}] = s;
        col.e_vals[{j, o}] = e;
    }
    return col;
}

// ---------------------------------------------------------------------------
// Instance  –  processed FJSSP problem data
// ---------------------------------------------------------------------------

struct Instance {
    std::string name;
    int num_jobs     = 0;
    int num_machines = 0;

    /// Number of operations per job: num_ops[j]
    std::vector<int> num_ops;

    /// Eligible machines per operation: elig_mach[j][o] = sorted list of machine IDs
    std::vector<std::vector<std::vector<int>>> elig_mach;

    /// Processing times: proc_time[j][o][m] = processing time
    std::vector<std::vector<std::map<int, int>>> proc_time;

    /// Sorted list of all unique machine IDs
    std::vector<int> all_machines;
};

// ---------------------------------------------------------------------------
// Hash helpers
// ---------------------------------------------------------------------------

struct PairHash {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 * 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

struct VecIntHash {
    size_t operator()(const std::vector<int>& v) const {
        size_t seed = v.size();
        for (int x : v) {
            seed ^= static_cast<size_t>(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// ---------------------------------------------------------------------------
// Dual variables (extracted from LP solution)
// ---------------------------------------------------------------------------

struct Duals {
    std::map<int, double>        alpha;  // machine coverage duals
    std::map<Op, double>         beta;   // operation assignment duals
    std::map<Op, double>         theta;  // precedence duals
    std::map<int, double>        mu;     // makespan duals
};

// ---------------------------------------------------------------------------
// Algorithm constants
// ---------------------------------------------------------------------------

constexpr double RC_TOL               = -1e-4;   // negative-RC threshold
constexpr int    CG_MAX_ITER          = 5000;    // max CG iterations (root)
constexpr double BIG_M                = 1e7;     // artificial variable penalty
constexpr double STAB_ALPHA           = 0.1;     // dual smoothing weight
constexpr double T_BUFFER             = 1.5;     // horizon = T_BUFFER * heuristic cmax
constexpr int    TOP_K_PER_MACHINE    = 8;       // max improving cols per machine per call 8
constexpr int    TOP_K_GLOBAL_COLS    = 20;      // max cols added to RMP per CG iteration 10
constexpr int    BEAM_WIDTH           = 50;      // beam width for heuristic pricing 35
constexpr int    CANDIDATES_PER_LABEL = 10;      // candidates expanded per beam label 15
constexpr double INT_TOL              = 1e-5;    // integrality tolerance
constexpr double BOUND_EPS            = 1e-4;    // pruning tolerance
constexpr double LP_POS_TOL           = 1e-4;    // LP positive-value threshold
constexpr int    MAX_NODES            = 10000;   // B&P node limit
constexpr int    CG_NODE_MAX_ITER     = 2000;    // CG iteration limit per B&P node
constexpr double CG_TIME_LIMIT_ROOT   = 180.0;    // CG wall-clock limit at root node (s) 60
constexpr double CG_TIME_LIMIT_NODE   = 180.0;    // CG wall-clock limit at child nodes (s) 60
constexpr double CP_SAT_INITIAL_TL    = 60.0;     // CP-SAT time limit (initial UB)
constexpr double CP_SAT_NODE_TL       = 30.0;     // CP-SAT time limit at nodes
constexpr int    CP_SAT_NODE_FREQ     = 0;      // run CP-SAT heuristic every N nodes
constexpr double BAP_GLOBAL_TIME_LIMIT = 7200.0;  // global B&P wall-clock limit (s)

constexpr int    CG_STALL_LIMIT       = 10;
