#pragma once
// =============================================================================
// instances.h  –  FJSSP problem instances (C++ equivalent of instances.py)
//
// Two formats supported:
//   machine_pattern: {machine_id: [(op_idx, job_idx, proc_time), ...]}
//   job_op_time:     {machine_id: [(job_idx, op_idx, proc_time), ...]}
//                    (same as machine_pattern but with j/o swapped)
//
// All instances are converted to the canonical Instance struct from types.h.
// =============================================================================

#include "types.h"
#include <stdexcept>
#include <algorithm>
#include <set>
#include <string>
#include <map>
#include <vector>

#include <array>
#include <cstdint>

// ---------------------------------------------------------------------------
// Raw instance description
// ---------------------------------------------------------------------------

/// Entry in a machine_pattern: (op_index, job_index, proc_time) [canonical order]
/// or (job_index, op_index, proc_time) when tuple_order == "job_op_time"
struct RawEntry {
    int first, second, proc_time;
};

enum class TupleOrder { OP_JOB_TIME, JOB_OP_TIME };

struct RawInstance {
    std::string name;
    int num_jobs;
    int num_machines;
    TupleOrder tuple_order = TupleOrder::OP_JOB_TIME;
    /// machine_id -> list of (first, second, proc_time)
    std::map<int, std::vector<RawEntry>> machine_pattern;
};

// ---------------------------------------------------------------------------
// Instance loader
// ---------------------------------------------------------------------------

inline Instance load_raw_instance(const RawInstance& raw) {
    Instance inst;
    inst.name         = raw.name;
    inst.num_jobs     = raw.num_jobs;
    inst.num_machines = raw.num_machines;

    int J = raw.num_jobs;
    std::vector<int> n_ops_per_job(J, 0);
    std::vector<std::vector<std::set<int>>> elig_sets(J);
    std::vector<std::vector<std::map<int, int>>> pt_data(J);

    // First pass: determine num_ops per job
    for (auto& [m, entries] : raw.machine_pattern) {
        for (auto& e : entries) {
            int o = (raw.tuple_order == TupleOrder::OP_JOB_TIME) ? e.first  : e.second;
            int j = (raw.tuple_order == TupleOrder::OP_JOB_TIME) ? e.second : e.first;
            if (o + 1 > n_ops_per_job[j]) n_ops_per_job[j] = o + 1;
        }
    }

    // Size arrays
    for (int j = 0; j < J; ++j) {
        elig_sets[j].resize(n_ops_per_job[j]);
        pt_data[j].resize(n_ops_per_job[j]);
    }

    // Second pass: fill eligibility and processing times
    for (auto& [m, entries] : raw.machine_pattern) {
        for (auto& e : entries) {
            int o = (raw.tuple_order == TupleOrder::OP_JOB_TIME) ? e.first  : e.second;
            int j = (raw.tuple_order == TupleOrder::OP_JOB_TIME) ? e.second : e.first;
            elig_sets[j][o].insert(m);
            pt_data[j][o][m] = e.proc_time;
        }
    }

    // Validate
    for (int j = 0; j < J; ++j) {
        for (int o = 0; o < n_ops_per_job[j]; ++o) {
            if (elig_sets[j][o].empty())
                throw std::runtime_error("Incomplete instance: op (" +
                    std::to_string(j) + "," + std::to_string(o) + ") has no eligible machine");
        }
    }

    // Build elig_mach as sorted vectors
    inst.num_ops  = n_ops_per_job;
    inst.elig_mach.resize(J);
    inst.proc_time.resize(J);
    for (int j = 0; j < J; ++j) {
        inst.elig_mach[j].resize(n_ops_per_job[j]);
        inst.proc_time[j].resize(n_ops_per_job[j]);
        for (int o = 0; o < n_ops_per_job[j]; ++o) {
            inst.elig_mach[j][o] = std::vector<int>(elig_sets[j][o].begin(), elig_sets[j][o].end());
            std::sort(inst.elig_mach[j][o].begin(), inst.elig_mach[j][o].end());
            inst.proc_time[j][o] = pt_data[j][o];
        }
    }

    // Collect all machines
    std::set<int> all_m_set;
    for (auto& [m, _] : raw.machine_pattern) all_m_set.insert(m);
    inst.all_machines = std::vector<int>(all_m_set.begin(), all_m_set.end());
    std::sort(inst.all_machines.begin(), inst.all_machines.end());

    return inst;
}

// Instance load_fjs_instance(const std::string& path);

// ---------------------------------------------------------------------------
// Instance definitions
// ---------------------------------------------------------------------------

inline RawInstance make_small_6j6m() {
    RawInstance r;
    r.name = "small_6jobs_6machines";
    r.num_jobs = 6; r.num_machines = 6;
    r.tuple_order = TupleOrder::OP_JOB_TIME;  // (op, job, time)
    r.machine_pattern = {
        {1, {{1,0,3},{4,1,10},{3,2,9},{1,3,5},{3,3,3},{4,4,3},{3,5,10}}},
        {2, {{2,0,6},{0,1,8},{4,2,1},{0,3,5},{1,4,3},{3,4,4},{0,5,3}}},
        {3, {{0,0,1},{1,1,5},{0,2,5},{2,3,5},{0,4,9},{5,5,1}}},
        {4, {{3,0,7},{4,0,3},{5,1,4},{1,2,4},{3,3,3},{5,4,1},{0,5,3},{1,5,3}}},
        {5, {{5,0,6},{2,1,10},{5,2,7},{4,3,8},{2,4,5},{4,5,4}}},
        {6, {{4,0,3},{3,1,10},{0,2,5},{2,2,8},{4,2,1},{5,3,9},{3,4,4},{2,5,9}}}
    };
    return r;
}

inline RawInstance make_medium_7j7m() {
    RawInstance r;
    r.name = "medium_7jobs_7machines";
    r.num_jobs = 7; r.num_machines = 7;
    r.tuple_order = TupleOrder::OP_JOB_TIME;
    r.machine_pattern = {
        {1, {{0,0,4},{6,1,7},{3,2,6},{5,2,3},{4,3,8},{3,4,5},{2,5,9},{1,6,2}}},
        {2, {{1,0,5},{0,1,8},{6,2,4},{1,3,5},{4,4,3},{5,0,4},{3,5,6},{5,3,7},{2,6,10}}},
        {3, {{2,0,7},{1,1,3},{0,2,6},{0,5,7},{6,3,9},{5,4,4},{4,5,8},{3,6,5},{6,1,2},{6,6,3}}},
        {4, {{3,0,9},{2,1,5},{1,2,7},{0,3,4},{0,0,6},{6,4,6},{5,5,3},{4,6,8}}},
        {5, {{4,0,3},{3,1,8},{2,2,5},{1,3,6},{0,4,9},{6,5,4},{5,6,7},{2,6,5}}},
        {6, {{5,0,6},{4,1,4},{2,3,8},{3,2,5},{1,4,7},{0,5,3},{4,3,5},{5,5,6},{6,6,9}}},
        {7, {{6,0,8},{5,1,6},{2,1,7},{4,2,4},{3,3,7},{2,4,3},{1,5,5},{4,4,9},{0,6,10}}}
    };
    return r;
}

inline RawInstance make_large_8j8m() {
    RawInstance r;
    r.name = "large_8jobs_8machines";
    r.num_jobs = 8; r.num_machines = 8;
    r.tuple_order = TupleOrder::OP_JOB_TIME;
    r.machine_pattern = {
        {1, {{0,0,5},{7,1,8},{6,2,3},{5,3,7},{4,4,4},{3,5,9},{2,6,6},{1,7,2},{1,0,7},{6,3,5}}},
        {2, {{1,0,6},{0,1,9},{7,2,4},{6,3,8},{5,4,3},{4,5,7},{3,6,5},{2,7,10},{2,0,4},{7,3,6}}},
        {3, {{2,0,8},{1,1,4},{0,2,7},{7,3,5},{6,4,9},{5,5,3},{4,6,8},{3,7,6},{3,0,5},{0,3,9}}},
        {4, {{3,0,4},{2,1,7},{1,2,5},{0,3,6},{7,4,8},{6,5,4},{5,6,9},{4,7,3},{4,0,6},{1,3,8}}},
        {5, {{4,0,7},{3,1,5},{2,2,9},{1,3,4},{0,4,6},{7,5,8},{6,6,3},{5,7,7},{5,0,3},{2,3,7}}},
        {6, {{5,0,3},{4,1,8},{3,2,6},{2,3,9},{1,4,5},{0,5,7},{7,6,4},{6,7,8},{6,0,9},{3,3,4}}},
        {7, {{6,0,9},{5,1,3},{4,2,8},{3,3,6},{2,4,7},{1,5,5},{0,6,4},{7,7,6},{7,0,5},{4,3,3}}},
        {8, {{7,0,6},{6,1,7},{5,2,4},{4,3,8},{3,4,5},{2,5,6},{1,6,9},{0,7,3},{0,0,8},{5,3,9}}}
    };
    return r;
}

inline RawInstance make_xlarge_9j9m() {
    RawInstance r;
    r.name = "xlarge_9jobs_9machines";
    r.num_jobs = 9; r.num_machines = 9;
    r.tuple_order = TupleOrder::OP_JOB_TIME;
    r.machine_pattern = {
        {1, {{0,0,6},{8,1,4},{7,2,9},{6,3,3},{5,4,7},{4,5,5},{3,6,8},{2,7,4},{1,8,6},{1,0,8},{7,3,5}}},
        {2, {{1,0,4},{0,1,7},{8,2,5},{7,3,9},{6,4,3},{5,5,8},{4,6,6},{3,7,5},{2,8,9},{2,0,6},{8,3,4}}},
        {3, {{2,0,9},{1,1,5},{0,2,6},{8,3,7},{7,4,4},{6,5,3},{5,6,9},{4,7,8},{3,8,5},{3,0,4},{0,3,8}}},
        {4, {{3,0,5},{2,1,8},{1,2,4},{0,3,6},{8,4,9},{7,5,5},{6,6,7},{5,7,3},{4,8,8},{4,0,7},{1,3,6}}},
        {5, {{4,0,3},{3,1,6},{2,2,8},{1,3,5},{0,4,4},{8,5,9},{7,6,3},{6,7,7},{5,8,6},{5,0,9},{2,3,3}}},
        {6, {{5,0,8},{4,1,3},{3,2,7},{2,3,4},{1,4,6},{0,5,5},{8,6,4},{7,7,9},{6,8,3},{6,0,5},{3,3,7}}},
        {7, {{6,0,4},{5,1,9},{4,2,3},{3,3,8},{2,4,5},{1,5,6},{0,6,7},{8,7,4},{7,8,9},{7,0,3},{4,3,9}}},
        {8, {{7,0,7},{6,1,4},{5,2,9},{4,3,3},{3,4,8},{2,5,6},{1,6,5},{0,7,8},{8,8,4},{8,0,6},{5,3,8}}},
        {9, {{8,0,5},{7,1,8},{6,2,6},{5,3,4},{4,4,9},{3,5,7},{2,6,3},{1,7,6},{0,8,7},{0,0,9},{6,3,6}}}
    };
    return r;
}

inline RawInstance make_xlarge_10j10m() {
    RawInstance r;
    r.name = "xlarge_10jobs_10machines";
    r.num_jobs = 10; r.num_machines = 10;
    r.tuple_order = TupleOrder::OP_JOB_TIME;
    r.machine_pattern = {
        {1,  {{0,0,6},{9,1,4},{8,2,9},{7,3,3},{6,4,7},{5,5,5},{4,6,8},{3,7,4},{2,8,6},{1,9,7},{1,0,8},{8,3,5}}},
        {2,  {{1,0,4},{0,1,7},{9,2,5},{8,3,9},{7,4,3},{6,5,8},{5,6,6},{4,7,5},{3,8,9},{2,9,4},{2,0,6},{9,3,4}}},
        {3,  {{2,0,9},{1,1,5},{0,2,6},{9,3,7},{8,4,4},{7,5,3},{6,6,9},{5,7,8},{4,8,5},{3,9,6},{3,0,4},{0,3,8}}},
        {4,  {{3,0,5},{2,1,8},{1,2,4},{0,3,6},{9,4,9},{8,5,5},{7,6,7},{6,7,3},{5,8,8},{4,9,6},{4,0,7},{1,3,6}}},
        {5,  {{4,0,3},{3,1,6},{2,2,8},{1,3,5},{0,4,4},{9,5,9},{8,6,3},{7,7,7},{6,8,6},{5,9,8},{5,0,9},{2,3,3}}},
        {6,  {{5,0,8},{4,1,3},{3,2,7},{2,3,4},{1,4,6},{0,5,5},{9,6,4},{8,7,9},{7,8,3},{6,9,7},{6,0,5},{3,3,7}}},
        {7,  {{6,0,4},{5,1,9},{4,2,3},{3,3,8},{2,4,5},{1,5,6},{0,6,7},{9,7,4},{8,8,9},{7,9,5},{7,0,3},{4,3,9}}},
        {8,  {{7,0,7},{6,1,4},{5,2,9},{4,3,3},{3,4,8},{2,5,6},{1,6,5},{0,7,8},{9,8,4},{8,9,6},{8,0,6},{5,3,8}}},
        {9,  {{8,0,5},{7,1,8},{6,2,6},{5,3,4},{4,4,9},{3,5,7},{2,6,3},{1,7,6},{0,8,7},{9,9,5},{9,0,9},{6,3,6}}},
        {10, {{9,0,6},{8,1,5},{7,2,8},{6,3,7},{5,4,4},{4,5,9},{3,6,6},{2,7,5},{1,8,8},{0,9,4},{0,0,7},{7,3,5}}}
    };
    return r;
}

inline RawInstance make_xlarge_11j11m() {
    RawInstance r;
    r.name = "xlarge_11jobs_11machines_flexible";
    r.num_jobs = 11;
    r.num_machines = 11;
    r.tuple_order = TupleOrder::OP_JOB_TIME;
    r.machine_pattern = {
        {1,  {{0,0,5},{10,1,9},{9,2,6},{8,3,3},{7,4,7},{6,5,6},{5,6,3},{4,7,7},{3,8,4},{2,9,8},{1,10,7},
              {1,0,4},{3,3,5},{6,6,7},{7,9,7}}},
        {2,  {{1,0,7},{0,1,5},{10,2,3},{9,3,3},{8,4,8},{7,5,8},{6,6,6},{5,7,4},{4,8,4},{3,9,9},{2,10,9},
              {2,1,5},{4,4,6},{6,7,7},{8,10,8}}},
        {3,  {{2,0,9},{1,1,8},{0,2,9},{10,3,8},{9,4,9},{8,5,3},{7,6,9},{6,7,3},{5,8,9},{4,9,3},{3,10,4},
              {3,2,6},{5,5,7},{7,8,8},{9,0,4}}},
        {4,  {{3,0,4},{2,1,4},{1,2,6},{0,3,8},{10,4,3},{9,5,5},{8,6,5},{7,7,7},{6,8,9},{5,9,4},{4,10,6},
              {4,3,7},{6,6,8},{8,9,9},{10,1,5}}},
        {5,  {{4,0,6},{3,1,9},{2,2,5},{1,3,8},{0,4,4},{10,5,7},{9,6,3},{8,7,6},{7,8,9},{6,9,5},{5,10,8},
              {5,4,8},{7,7,9},{9,10,4},{0,2,7}}},
        {6,  {{5,0,8},{4,1,5},{3,2,9},{2,3,6},{1,4,3},{0,5,9},{10,6,6},{9,7,3},{8,8,7},{7,9,4},{6,10,3},
              {6,5,9},{9,8,5},{10,0,6},{1,3,8}}},
        {7,  {{6,0,3},{5,1,8},{4,2,6},{3,3,6},{2,4,4},{1,5,4},{0,6,9},{10,7,7},{9,8,7},{8,9,5},{7,10,5},
              {7,6,4},{9,9,5},{0,1,8},{3,4,4}}},
        {8,  {{7,0,5},{6,1,4},{5,2,5},{4,3,4},{3,4,5},{2,5,6},{1,6,5},{0,7,6},{10,8,5},{9,9,6},{8,10,7},
              {8,7,5},{10,10,6},{1,2,9},{3,5,4}}},
        {9,  {{8,0,7},{7,1,7},{6,2,9},{5,3,4},{4,4,6},{3,5,8},{2,6,8},{1,7,3},{0,8,5},{10,9,7},{9,10,9},
              {9,8,6},{0,0,9},{2,3,4},{4,6,5}}},
        {10, {{9,0,9},{8,1,5},{7,2,8},{6,3,4},{5,4,7},{4,5,3},{3,6,6},{2,7,9},{1,8,5},{0,9,8},{10,10,4},
              {10,9,7},{1,1,4},{3,4,5},{5,7,6}}},
        {11, {{10,0,4},{9,1,8},{8,2,5},{7,3,9},{6,4,6},{5,5,5},{4,6,9},{3,7,6},{2,8,3},{1,9,7},{0,10,6},
              {1,10,4},{2,2,5},{4,5,6},{6,8,7}}}
    };
    return r;
}


inline RawInstance make_balanced_cyclic_high_flex(
    const std::string& name,
    int n,
    const std::vector<int>& shifts
) {
    RawInstance r;
    r.name = name;
    r.num_jobs = n;
    r.num_machines = n;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    auto nmod = [&](int a) -> int {
        int v = a % n;
        return v < 0 ? v + n : v;
    };

    for (int m = 1; m <= n; ++m) {
        std::vector<RawEntry> entries;
        entries.reserve(n * (int)shifts.size());

        for (int j = 0; j < n; ++j) {
            for (int a = 0; a < (int)shifts.size(); ++a) {
                int op = nmod((m - 1) - j - shifts[a]);

                // Same modular processing-time structure as the two successful instances.
                int p = 3 + nmod(5 * op + 7 * j + 4 * m + 3 * a);

                entries.push_back({op, j, p});
            }
        }

        r.machine_pattern[m] = std::move(entries);
    }

    return r;
}

inline RawInstance make_bhf_09j09m_k06() {
    return make_balanced_cyclic_high_flex(
        "bhf_09j09m_k06",
        9,
        {0, 1, 2, 4, 6, 8}
    );
}

inline RawInstance make_bhf_10j10m_k06() {
    return make_balanced_cyclic_high_flex(
        "bhf_10j10m_k06",
        10,
        {0, 1, 2, 4, 6, 8}
    );
}

inline RawInstance make_bhf_11j11m_k07() {
    return make_balanced_cyclic_high_flex(
        "bhf_11j11m_k07",
        11,
        {0, 1, 2, 3, 5, 7, 9}
    );
}

inline RawInstance make_bhf_14j14m_k08() {
    return make_balanced_cyclic_high_flex(
        "bhf_14j14m_k08",
        14,
        {0, 1, 2, 3, 5, 7, 9, 11}
    );
}

inline RawInstance make_bhf_14j14m_k09_dense() {
    return make_balanced_cyclic_high_flex(
        "bhf_14j14m_k09_dense",
        14,
        {0, 1, 2, 3, 4, 6, 8, 10, 12}
    );
}

inline RawInstance make_bhf_15j15m_k09() {
    return make_balanced_cyclic_high_flex(
        "bhf_15j15m_k09",
        15,
        {0, 1, 2, 3, 5, 7, 9, 11, 13}
    );
}

inline RawInstance make_bhf_15j15m_k10_dense() {
    return make_balanced_cyclic_high_flex(
        "bhf_15j15m_k10_dense",
        15,
        {0, 1, 2, 3, 4, 6, 8, 10, 12, 14}
    );
}

inline RawInstance make_bhf_16j16m_k10() {
    return make_balanced_cyclic_high_flex(
        "bhf_16j16m_k10",
        16,
        {0, 1, 2, 3, 5, 7, 9, 11, 13, 15}
    );
}

inline RawInstance make_xlarge_15j15m_very_high_flex()
{
    RawInstance r;

    r.name =
        "xlarge_15jobs_15machines_very_high_flex";

    r.num_jobs     = 15;
    r.num_machines = 15;
    r.tuple_order  = TupleOrder::OP_JOB_TIME;

    /*
     * 15 jobs × 15 operations × 8 alternative machines
     * = 1800 operation-machine tuples.
     */
    const int n = 15;

    const std::vector<int> shifts = {
        0, 1, 2, 3, 5, 7, 9, 11
    };

    auto nmod = [&](int value) -> int {
        const int result = value % n;
        return result < 0
            ? result + n
            : result;
    };

    for (int m = 1; m <= n; ++m) {
        std::vector<RawEntry> entries;

        // 15 jobs × 8 alternatives = 120 entries per machine.
        entries.reserve(
            n * static_cast<int>(shifts.size())
        );

        for (int j = 0; j < n; ++j) {
            for (int a = 0;
                 a < static_cast<int>(shifts.size());
                 ++a)
            {
                const int op =
                    nmod(
                        (m - 1)
                        - j
                        - shifts[a]
                    );

                const int p =
                    3
                    + nmod(
                        5  * op
                        + 7  * j
                        + 4 * m
                        + 3  * a
                    );

                entries.push_back({
                    op,
                    j,
                    p
                });
            }
        }

        r.machine_pattern[m] =
            std::move(entries);
    }

    return r;
}

inline RawInstance make_xlarge_14j14m_very_high_flex()
{
    RawInstance r;

    r.name =
        "xlarge_14jobs_14machines_very_high_flex";

    r.num_jobs = 14;
    r.num_machines = 14;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    /*
     * 14 jobs × 14 operations × 8 alternative machines
     * = 1568 operation-machine tuples.
     */
    const int n = 14;

    const std::vector<int> shifts = {
        0, 1, 2, 3, 5, 7, 9, 11
    };

    auto nmod = [&](int value) -> int {
        const int result = value % n;
        return result < 0
            ? result + n
            : result;
    };

    for (int m = 1; m <= n; ++m) {
        std::vector<RawEntry> entries;

        // 14 jobs × 8 alternatives = 112 entries per machine.
        entries.reserve(
            n * static_cast<int>(shifts.size())
        );

        for (int j = 0; j < n; ++j) {
            for (int a = 0;
                 a < static_cast<int>(shifts.size());
                 ++a)
            {

                const int op =
                    nmod(
                        (m - 1)
                        - j
                        - shifts[a]
                    );

                // const int p =
                //     3
                //     + nmod(
                //         5  * op
                //         + 9  * j
                //         + 11 * m
                //         + 3  * a
                //     );

                const int p =
                    3
                    + nmod(
                        5  * op
                        + 9  * j
                        + 10 * m
                        + 3  * a
                    );

                // const int p =
                //     3 + nmod(
                //         5 * op
                //         + 8 * j
                //         + 5 * m
                //         + 4 * a
                //     );

                entries.push_back({
                    op,
                    j,
                    p
                });
            }
        }

        r.machine_pattern[m] =
            std::move(entries);
    }

    return r;
}


inline RawInstance make_xlarge_13j13m_very_high_flex() {
    RawInstance r;
    r.name = "xlarge_13jobs_13machines_very_high_flex";
    r.num_jobs = 13;
    r.num_machines = 13;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    // Very high flexibility: 13 jobs x 13 ops x 8 alternative machines = 1352 tuples.
    const int n = 13;
    const std::vector<int> shifts = {0, 1, 2, 3, 5, 7, 9, 11};

    auto nmod = [&](int a) -> int {
        int v = a % n;
        return v < 0 ? v + n : v;
    };

    for (int m = 1; m <= n; ++m) {
        std::vector<RawEntry> entries;
        entries.reserve(n * (int)shifts.size());
        for (int j = 0; j < n; ++j) {
            for (int a = 0; a < (int)shifts.size(); ++a) {
                int op = nmod((m - 1) - j - shifts[a]);
                int p  = 3 + nmod(5 * op + 7 * j + 4 * m + 3 * a);
                entries.push_back({op, j, p});
            }
        }
        r.machine_pattern[m] = std::move(entries);
    }

    return r;
}


inline RawInstance make_xlarge_12j12m_very_high_flex() {
    RawInstance r;
    r.name = "xlarge_12jobs_12machines_very_high_flex";
    r.num_jobs = 12;
    r.num_machines = 12;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    // Very high flexibility: 12 jobs x 12 ops x 8 alternative machines = 1152 tuples.
    // Shifts define the 8 alternative machines per (job, op) pair.
    const int n = 12;
    const std::vector<int> shifts = {0, 1, 2, 3, 5, 7, 9, 11};

    auto nmod = [&](int a) -> int {
        int v = a % n;
        return v < 0 ? v + n : v;
    };

    for (int m = 1; m <= n; ++m) {
        std::vector<RawEntry> entries;
        entries.reserve(n * (int)shifts.size());
        for (int j = 0; j < n; ++j) {
            for (int a = 0; a < (int)shifts.size(); ++a) {
                int op = nmod((m - 1) - j - shifts[a]);
                int p  = 3 + nmod(5 * op + 7 * j + 4 * m + 3 * a);
                entries.push_back({op, j, p});
            }
        }
        r.machine_pattern[m] = std::move(entries);
    }

    return r;
}


inline RawInstance make_xlarge_11j11m_very_high_flex()
{
    RawInstance r;

    r.name =
        "xlarge_11jobs_11machines_very_high_flex";

    r.num_jobs = 11;
    r.num_machines = 11;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    /*
     * Very high flexibility:
     * 11 jobs × 11 operations × 8 alternative machines
     * = 968 operation-machine tuples.
     */
    const int n = 11;

    /*
     * Eight distinct shifts modulo 11.
     *
     * Shift 11 from the 12×12 instance cannot be used,
     * because 11 mod 11 = 0, which would duplicate shift 0.
     */
    const std::vector<int> shifts = {
        0, 1, 2, 3, 5, 7, 9, 10
    };

    auto nmod = [&](int value) -> int {
        const int result = value % n;

        return result < 0
            ? result + n
            : result;
    };

    for (int m = 1; m <= n; ++m) {
        std::vector<RawEntry> entries;

        /*
         * Every machine receives:
         * 11 jobs × 8 alternatives = 88 entries.
         */
        entries.reserve(
            n * static_cast<int>(shifts.size())
        );

        for (int j = 0; j < n; ++j) {
            for (
                int a = 0;
                a < static_cast<int>(shifts.size());
                ++a
            )
            {
                /*
                 * Same cyclic eligibility rule used by
                 * the 12×12 instance.
                 */
                const int op =
                    nmod(
                        (m - 1)
                        - j
                        - shifts[a]
                    );

                /*
                 * Same processing-time formula used by
                 * the 12×12 instance.
                 *
                 * Processing times range from 3 to 13.
                 */
                const int p =
                    3 + nmod(
                        6 * op
                        + 7 * j
                        + 7 * m
                        + 3 * a
                    );

                entries.push_back({
                    op,
                    j,
                    p
                });
            }
        }

        r.machine_pattern[m] =
            std::move(entries);
    }

    return r;
}


inline RawInstance
make_xlarge_11j11m_variable_8_to_10_flex()
{
    RawInstance r;

    r.name =
        "xlarge_11jobs_11machines_8_to_10_alternatives_variable_flex";

    r.num_jobs = 11;
    r.num_machines = 11;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 11;
    constexpr int min_alternatives = 8;
    constexpr int max_alternatives = 10;

    /*
     * 121 operations in total.
     *
     * Distribution:
     * 40 operations × 8 machines
     * 40 operations × 9 machines
     * 41 operations × 10 machines
     *
     * Total tuples:
     * 40×8 + 40×9 + 41×10 = 1090.
     *
     * Average alternatives per operation:
     * 1090 / 121 ≈ 9.01.
     */

    /*
     * The first eight shifts are the same as in the
     * original 8-alternative instance.
     *
     * Shifts 4 and 6 are added for operations having
     * 9 or 10 alternatives.
     *
     * All shifts are distinct modulo 11.
     */
    const std::array<int, max_alternatives> shifts = {
        0, 1, 2, 3, 5, 7, 9, 10, 4, 6
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    /*
     * Approximately:
     * 1090 tuples / 11 machines ≈ 99 entries per machine.
     *
     * Reserving 110 avoids vector reallocations.
     */
    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(
            n * max_alternatives
        );
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            const int operation_index =
                j * n + op;

            /*
             * operation_index mod 3:
             *
             * 0 -> 10 alternatives
             * 1 ->  9 alternatives
             * 2 ->  8 alternatives
             *
             * Because there are 121 operations, this produces:
             * 41 operations with 10 alternatives,
             * 40 operations with 9 alternatives,
             * 40 operations with 8 alternatives.
             */
            const int num_alternatives =
                min_alternatives
                + nmod(
                    2 * operation_index + 2,
                    max_alternatives
                        - min_alternatives
                        + 1
                );

            for (int a = 0;
                 a < num_alternatives;
                 ++a)
            {
                /*
                 * Equivalent to the eligibility rule:
                 *
                 * op = (machine - 1 - job - shift[a]) mod n.
                 *
                 * Every operation receives 8, 9 or 10
                 * distinct eligible machines.
                 */
                const int machine_zero_based =
                    nmod(
                        op + j + shifts[a],
                        n
                    );

                const int machine =
                    machine_zero_based + 1;

                /*
                 * Same processing-time formula as the
                 * original instance.
                 *
                 * Since the result of nmod is in [0,10],
                 * processing times range from 3 to 13.
                 */
                const int p =
                    3 + nmod(
                        6 * op
                        + 7 * j
                        + 7 * machine
                        + 3 * a,
                        n
                    );

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}


inline RawInstance
make_xlarge_11j11m_very_high_flex_10alt()
{
    RawInstance r;

    r.name =
        "xlarge_11jobs_11machines_10_alternatives_cg_branching_friendly_v2";

    r.num_jobs = 11;
    r.num_machines = 11;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 11;
    constexpr int alternatives = 10;

    /*
     * Two orientations of the eligible-machine sequence.
     * Both steps visit all 11 machines.
     */
    const std::array<int, 2> machine_steps = {
        1, 2
    };

    /*
    * Profile used by ordinary operations:
    * four competitive alternatives.
    */
    const std::array<int, alternatives>
        regular_quality_penalty = {
            0, 1, 2, 3, 7,
            9, 11, 13, 15, 17
        };

    /*
    * Profile used by critical operations:
    * only two genuinely competitive alternatives.
    */
    const std::array<int, alternatives>
        critical_quality_penalty = {
            0, 1, 5, 7, 9,
            11, 13, 15, 17, 19
        };

    /*
     * Irregular permutations used to balance the excluded
     * machine exactly across all operations.
     */
    const std::array<int, n> job_code = {
        0, 4, 9, 2, 7, 1,
        10, 5, 3, 8, 6
    };

    const std::array<int, n> operation_code = {
        0, 7, 3, 10, 5, 1,
        9, 4, 8, 2, 6
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    /*
     * There are 1210 tuples:
     *
     * 11 jobs × 11 operations × 10 machines.
     *
     * Every machine receives exactly 110 tuples.
     */
    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(
            n * alternatives
        );
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            /*
            * Only 22 of the 121 operations are treated as critical.
            * This should increase the tree moderately.
            */
            const bool critical_operation =
                op == 4 || op == 7;

            /*
             * Every operation excludes exactly one machine.
             * Every machine is excluded exactly 11 times.
             */
            const int excluded_machine =
                nmod(
                    job_code[j]
                    + operation_code[op],
                    n
                );

            /*
             * Controlled machine neighbourhoods.
             *
             * Operations belonging to similar job groups
             * and processing stages receive overlapping
             * preferred machines. This produces meaningful
             * sequence conflicts without overloading one machine.
             */
            const int job_group =
                j % 3;

            const int operation_stage =
                op / 3;

            const int job_layer =
                j / 3;

            // const int machine_start =
            //     nmod(
            //         4 * job_group
            //         + 5 * operation_stage
            //         + 2 * (op % 3)
            //         + job_layer,
            //         n
            //     );


            const int machine_start =
                nmod(
                    4 * job_group
                    + 5 * operation_stage
                    + 2 * (op % 3)
                    + job_layer
                    + 3 * ((j * op) % 4),
                    n
                );

            // const int machine_step =
            //     machine_steps[
            //         nmod(
            //             j + op + operation_stage,
            //             static_cast<int>(
            //                 machine_steps.size()
            //             )
            //         )
            //     ];

            const int step_selector =
                nmod(
                    j
                    + 2 * op
                    + (j * op) % 3,
                    4
                );

            const int machine_step =
                step_selector == 0
                    ? 2
                    : 1;


            /*
             * Common offset for the complete operation.
             * Processing times remain in the range 3--15.
             */
            const int operation_offset =
                nmod(
                    j
                    + 2 * op
                    + (j * op) % 3,
                    2
                );

            /*
             * Swap only the two competitive alternatives.
             * This avoids associating the best processing
             * time permanently with the first machine.
             */

            // const int fast_swap =
            //     nmod(
            //         j
            //         + 2 * op
            //         + (j * op) % 2,
            //         2
            //     );

            const int fast_rotation =
                nmod(
                    2 * j
                    + 3 * op
                    + (j * op) % 3,
                    4
                );

            const int fast_step =
                nmod(
                    j + op + operation_stage,
                    2
                ) == 0
                    ? 1
                    : 3;

            int alternative_index = 0;

            for (int position = 0;
                 position < n;
                 ++position)
            {
                const int machine_zero_based =
                    nmod(
                        machine_start
                        + machine_step * position,
                        n
                    );

                if (machine_zero_based ==
                    excluded_machine)
                {
                    continue;
                }

                const int machine =
                    machine_zero_based + 1;

                // int quality_rank =
                //     alternative_index;

                // /*
                //  * The first two machines remain the competitive
                //  * alternatives, but their quality order is swapped.
                //  */
                // if (alternative_index < 2) {
                //     quality_rank =
                //         nmod(
                //             alternative_index
                //             + fast_swap,
                //             2
                //         );
                // }

                // int quality_rank =
                //     alternative_index;

                int quality_rank;

                if (critical_operation) {

                    /*
                    * The first two eligible machines receive ranks 0 and 1.
                    * Their order is swapped between operations.
                    */
                    const int fast_swap =
                        nmod(
                            j + op,
                            2
                        );

                    if (alternative_index < 2) {
                        quality_rank =
                            nmod(
                                alternative_index
                                + fast_swap,
                                2
                            );
                    }
                    else {
                        quality_rank =
                            alternative_index;
                    }
                }
                else {

                    /*
                    * Four competitive alternatives for ordinary operations.
                    */
                    const int fast_rotation =
                        nmod(
                            2 * j
                            + 3 * op
                            + (j * op) % 3,
                            4
                        );

                    const int fast_step =
                        nmod(
                            j + op + op / 3,
                            2
                        ) == 0
                            ? 1
                            : 3;

                    if (alternative_index < 4) {
                        quality_rank =
                            nmod(
                                fast_rotation
                                + fast_step * alternative_index,
                                4
                            );
                    }
                    else {
                        quality_rank =
                            alternative_index;
                    }
                }

                if (alternative_index < 4) {
                    quality_rank =
                        nmod(
                            fast_rotation
                            + fast_step * alternative_index,
                            4
                        );
                }

                const auto& selected_penalty =
                    critical_operation
                        ? critical_quality_penalty
                        : regular_quality_penalty;

                const int p =
                    3
                    + operation_offset
                    + selected_penalty[quality_rank];

                // const int p =
                //     3
                //     + operation_offset
                //     + quality_penalty[quality_rank];

                // r.machine_pattern[machine].push_back({
                //     op,
                //     j,
                //     p
                // });

                ++alternative_index;
            }

            if (alternative_index != alternatives) {
                throw std::runtime_error(
                    "Incorrect number of alternatives"
                );
            }
        }
    }

    return r;
}


inline RawInstance make_xlarge_11j11m_high_flex_6alt()
{
    RawInstance r;

    r.name =
        "xlarge_11jobs_11machines_6_alternatives_high_flex";

    r.num_jobs = 11;
    r.num_machines = 11;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    /*
     * 11 jobs × 11 operations × 6 alternative machines
     * = 726 operation-machine tuples.
     *
     * Relative flexibility:
     * 6 / 11 = 54.55%.
     */
    constexpr int n = 11;
    constexpr int alternatives = 6;

    /*
     * Six distinct shifts modulo 11.
     *
     * These are taken from the original eight-shift set,
     * preserving the same cyclic eligibility structure.
     */
    const std::array<int, alternatives> shifts = {
        0, 1, 2, 3, 5, 7
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    for (int m = 1; m <= n; ++m) {
        std::vector<RawEntry> entries;

        /*
         * Every machine receives:
         *
         * 11 jobs × 6 alternatives = 66 entries.
         */
        entries.reserve(
            n * alternatives
        );

        for (int j = 0; j < n; ++j) {
            for (int a = 0; a < alternatives; ++a) {

                /*
                 * Cyclic eligibility rule.
                 *
                 * For every (job, operation) pair, exactly
                 * six distinct machines are generated.
                 */
                const int op =
                    nmod(
                        (m - 1)
                        - j
                        - shifts[a],
                        n
                    );

                /*
                 * Processing times range from 3 to 13.
                 */
                const int p =
                    3
                    + nmod(
                        6 * op
                        + 7 * j
                        + 7 * m
                        + 3 * a,
                        n
                    );

                entries.push_back({
                    op,
                    j,
                    p
                });
            }
        }

        r.machine_pattern[m] =
            std::move(entries);
    }

    return r;
}


inline RawInstance
// make_xlarge_12j12m_variable_6_to_8_flex()
make_xlarge_11j11m_variable_6_to_8_flex()
{
    RawInstance r;

    r.name =
        "xlarge_12jobs_12machines_6_to_8_alternatives_bp_structured";

    r.num_jobs = 12;
    r.num_machines = 12;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 12;
    constexpr int min_alternatives = 6;
    constexpr int max_alternatives = 8;

    /*
     * 12 jobs × 12 operations = 144 operations.
     *
     * Distribution:
     *
     * 48 operations × 6 alternatives
     * 48 operations × 7 alternatives
     * 48 operations × 8 alternatives
     *
     * Total tuples:
     *
     * 48×6 + 48×7 + 48×8 = 1008.
     *
     * Average alternatives per operation = 7.
     */

    /*
     * Regular operations have four competitive alternatives.
     */
    const std::array<int, max_alternatives>
        regular_quality_penalty = {
            0, 1, 2, 3,
            5, 7, 9, 11
        };

    /*
     * Critical operations have only two strongly
     * competitive alternatives.
     *
     * Their OMA fix/forbid decisions should have a
     * greater effect on the node lower bound.
     */
    const std::array<int, max_alternatives>
        critical_quality_penalty = {
            0, 2, 5, 7,
            9, 11, 13, 15
        };

    /*
     * All these values are coprime with 12.
     * Therefore, each one generates a permutation
     * of all twelve machines.
     */
    const std::array<int, 4> machine_steps = {
        1, 5, 7, 11
    };

    /*
     * Mild machine heterogeneity.
     */
    const std::array<int, n> machine_penalty = {
        0, 1, 0, 1,
        0, 1, 0, 1,
        0, 1, 0, 1
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    /*
     * Maximum possible entries per machine:
     *
     * 12 jobs × 8 alternatives = 96.
     *
     * Actual average:
     *
     * 1008 / 12 = 84.
     */
    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(
            n * max_alternatives
        );
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            const int operation_index =
                j * n + op;

            /*
             * Exact distribution:
             *
             * operation_index mod 3:
             *
             * 0 -> 8 alternatives
             * 1 -> 7 alternatives
             * 2 -> 6 alternatives
             */
            const int num_alternatives =
                min_alternatives
                + nmod(
                    2 * operation_index + 2,
                    max_alternatives
                        - min_alternatives
                        + 1
                );

            /*
             * Four operation stages:
             *
             * stage 0: o0,o1,o2
             * stage 1: o3,o4,o5
             * stage 2: o6,o7,o8
             * stage 3: o9,o10,o11
             */
            const int operation_stage =
                op / 3;

            const int job_group =
                j % 4;

            const int job_layer =
                j / 4;

            /*
             * Jobs from related groups and operations from
             * adjacent stages have overlapping preferred
             * machine neighbourhoods.
             */
            const int machine_start =
                nmod(
                    2 * operation_stage
                    + job_group
                    + (job_layer + op) % 2,
                    n
                );

            /*
             * Change the orientation of the eligible-machine
             * sequence between operations.
             */
            const int machine_step =
                machine_steps[
                    nmod(
                        3 * j
                        + 5 * op
                        + j * op,
                        static_cast<int>(
                            machine_steps.size()
                        )
                    )
                ];

            /*
             * Operations o4 and o7 are critical.
             *
             * There are 24 critical operations:
             *
             * 12 jobs × 2 critical operations.
             */
            const bool critical_operation =
                op == 4 || op == 7;

            const auto& quality_penalty =
                critical_operation
                ? critical_quality_penalty
                : regular_quality_penalty;

            const int competitive_count =
                critical_operation
                ? 2
                : 4;

            /*
             * Base processing time in {3,4,5}.
             *
             * It varies between jobs and operations, preventing
             * all jobs from having the same critical-path length.
             */
            const int operation_base =
                3
                + nmod(
                    5 * j
                    + 7 * op
                    + 2 * j * op
                    + j * j
                    + op * op,
                    3
                );

            std::array<bool, n> machine_used{};

            for (int a = 0;
                 a < num_alternatives;
                 ++a)
            {
                /*
                 * machine_step is coprime with 12, so the generated
                 * machine indices are distinct.
                 */
                const int machine_zero_based =
                    nmod(
                        machine_start
                        + machine_step * a,
                        n
                    );

                if (machine_used[machine_zero_based]) {
                    throw std::runtime_error(
                        "Duplicate eligible machine"
                    );
                }

                machine_used[machine_zero_based] = true;

                const int machine =
                    machine_zero_based + 1;

                const bool backup =
                    a >= competitive_count;

                /*
                 * Do not perturb competitive alternatives.
                 */
                const int backup_jitter =
                    backup
                    ? nmod(
                        5 * op
                        + 7 * j
                        + 3 * machine
                        + a,
                        2
                    )
                    : 0;

                const int machine_effect =
                    backup
                    ? machine_penalty[
                        machine_zero_based
                    ]
                    : 0;

                const int p =
                    operation_base
                    + quality_penalty[a]
                    + machine_effect
                    + backup_jitter;

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}


inline RawInstance make_xlarge_11j11m_very_high_flex_9alt()
{
    RawInstance r;

    r.name =
        "xlarge_10jobs_10machines_9_alternatives_very_high_flex";

    r.num_jobs = 11;
    r.num_machines = 11;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 11;
    constexpr int alternatives = 9;

    /*
     * 10 jobs × 10 operations × 9 alternatives
     * = 900 operation-machine tuples.
     *
     * Relative flexibility:
     * 9 / 10 = 90%.
     *
     * Every operation can be processed on nine
     * out of the ten available machines.
     */
    const std::array<int, alternatives> shifts = {
        0, 1, 2, 3, 4, 5, 6, 8, 9
    };

    /*
     * One fastest alternative,
     * one competitive alternative,
     * and seven progressively slower backups.
     */
    const std::array<int, alternatives> quality_penalty = {
        0, 1, 3, 5, 7, 9, 11, 13, 15
    };

    /*
    * All these values are coprime with 9.
    * Therefore, each one generates a permutation
    * of quality ranks 0,...,8.
    */
    const std::array<int, 6> quality_steps = {
        1, 2, 4, 5, 7, 8
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    /*
     * Average entries per machine:
     * 900 tuples / 10 machines = 90 entries.
     */
    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(
            n * alternatives
        );
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            /*
             * Balanced but not completely symmetric
             * machine eligibility.
             */
            const int base_machine =
                nmod(
                    3 * op
                    + 7 * j
                    + (j * op) % 3,
                    n
                );

            /*
             * Operation-dependent offset in {0,1}.
             * It introduces moderate processing-time
             * asymmetry between operations.
             */
            const int operation_offset =
                nmod(
                    2 * j
                    + 3 * op
                    + (j * op) % 3
                    + (j + op * op) % 5,
                    2
                );

            for (int a = 0; a < alternatives; ++a) {

                /*
                 * The nine shifts are distinct modulo 10.
                 * Therefore, every operation receives exactly
                 * nine distinct eligible machines.
                 */
                const int machine_zero_based =
                    nmod(
                        base_machine + shifts[a],
                        n
                    );

                const int machine =
                    machine_zero_based + 1;

                /*
                 * Every operation receives all nine quality
                 * ranks exactly once.
                 */
                // const int quality_rank =
                //     nmod(
                //         a
                //         + 3 * j
                //         + (j * op) % 3,
                //         alternatives
                //     );

                /*
                * Irregular rotation of the quality ranks.
                */
                const int quality_rotation =
                    nmod(
                        8 * j
                        + 3 * op
                        + 5 * ((j * op) % 3)
                        + 4 * ((j * j + op * op) % 2),
                        alternatives
                    );

                /*
                * Change the order of quality ranks between operations.
                * Since every quality_step is coprime with 9,
                * all ranks from 0 to 8 appear exactly once.
                */
                const int quality_step =
                    quality_steps[
                        nmod(
                            5 * j + 9 * op,
                            static_cast<int>(quality_steps.size())
                        )
                    ];

                const int quality_rank =
                    nmod(
                        quality_rotation
                        + quality_step * a,
                        alternatives
                    );

                /*
                 * The two competitive alternatives are not
                 * perturbed. Jitter is added only to backups.
                 */
                const int backup_jitter =
                    quality_rank >= 2
                    ? nmod(
                        3 * op
                        + 4 * j
                        + 3 * machine
                        + 2 * a,
                        2
                    )
                    : 0;

                /*
                 * Processing times range approximately
                 * from 3 to 20.
                 */
                const int p =
                    3
                    + operation_offset
                    + quality_penalty[quality_rank]
                    + backup_jitter;

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}

inline RawInstance make_xlarge_10j10m_very_high_flex()
{
    RawInstance r;

    r.name =
        "xlarge_10jobs_10machines_very_high_flex";

    r.num_jobs = 10;
    r.num_machines = 10;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 10;
    constexpr int alternatives = 6;

    /*
     * 10 jobs × 10 operations × 6 alternatives
     * = 600 operation-machine tuples.
     *
     * Relative flexibility = 6/10 = 60%.
     */
    const std::array<int, alternatives> shifts = {
        0, 1, 3, 6, 8, 9
    };

    /*
     * Alternatives now have clearly separated qualities.
     *
     * Including jitter, their processing-time intervals are:
     * rank 0: 3--4
     * rank 1: 5--6
     * rank 2: 7--8
     * rank 3: 9--10
     * rank 4: 11--12
     * rank 5: 13--14
     *
     * Therefore, alternatives of different ranks cannot have
     * identical processing times.
     */
    const std::array<int, alternatives> quality_penalty = {
        0, 2, 4, 6, 8, 10
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(
            n * alternatives
        );
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            /*
             * Two overlapping job groups.
             *
             * Jobs in the same group tend to use similar machines,
             * preserving competition and sequence-branching candidates.
             */
            const int job_group = j % 2;

            /*
             * Irregular but balanced base machine.
             */
            const int base_machine =
                nmod(
                    3 * op
                    + 2 * job_group
                    + (j * op) % 3,
                    n
                );

            for (int a = 0; a < alternatives; ++a) {

                const int machine_zero_based =
                    nmod(
                        base_machine + shifts[a],
                        n
                    );

                const int machine =
                    machine_zero_based + 1;

                /*
                 * Each operation receives exactly one alternative
                 * from every quality rank.
                 */
                const int quality_rank =
                    nmod(
                        a
                        + 2 * j
                        + 3 * op
                        + (j * op) % 2,
                        alternatives
                    );

                /*
                 * Small perturbation independent of quality rank.
                 */
                const int jitter =
                    nmod(
                        3 * op
                        + 5 * j
                        + 7 * machine,
                        2
                    );

                const int p =
                    3
                    + quality_penalty[quality_rank]
                    + jitter;

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}


inline RawInstance make_xlarge_10j10m_high_flex()
{
    RawInstance r;

    r.name =
        "xlarge_10jobs_10machines_8_alternatives_cg_branching_friendly";

    r.num_jobs = 10;
    r.num_machines = 10;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 10;
    constexpr int alternatives = 8;

    /*
     * 10 jobs × 10 operations × 8 alternatives
     * = 800 operation-machine tuples.
     */
    const std::array<int, alternatives> shifts = {
        0, 1, 2, 3, 4, 6, 8, 9
    };

    /*
     * One fastest alternative,
     * one competitive alternative,
     * and six slower backup alternatives.
     */\

    // const std::array<int, alternatives> quality_penalty = {
    //     0, 1, 4, 5, 6, 7, 8, 9
    // };

    // const std::array<int, alternatives> quality_penalty = {
    //     0, 2, 5, 7, 9, 11, 13, 15
    // };

    const std::array<int, alternatives> quality_penalty = {
        0, 1, 3, 5, 7, 9, 11, 13
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;
        return result < 0
            ? result + modulus
            : result;
    };

    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(
            n * alternatives
        );
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            /*
             * Balanced but not completely symmetric eligibility.
             */
            const int base_machine =
                nmod(
                    3 * op
                    + 7 * j
                    + (j * op) % 3,
                    n
                );

            /*
             * Small operation-dependent offset in {0,1}.
             */

            // const int operation_offset =
            //     nmod(
            //         j
            //         + 2 * op
            //         + (j * op) % 2,
            //         2
            //     );

            // const int operation_offset =
            //     nmod(
            //         3 * j
            //         + 5 * op
            //         + (j * op) % 4
            //         + (j * j + op * op) % 5,
            //         3
            //     );

            const int operation_offset =
                nmod(
                    2 * j
                    + 3 * op
                    + (j * op) % 3
                    + (j + op * op) % 5,
                    2
                );

            for (int a = 0; a < alternatives; ++a) {

                const int machine_zero_based =
                    nmod(
                        base_machine + shifts[a],
                        n
                    );

                const int machine =
                    machine_zero_based + 1;

                /*
                 * Rotate the quality assignment.
                 * Every operation receives every quality rank
                 * exactly once.
                 */
                const int quality_rank =
                    nmod(
                        a
                        + 3 * j
                        + (j * op) % 3,
                        alternatives
                    );

                /*
                 * Do not perturb the two competitive alternatives.
                 * This preserves a clear primary/secondary structure.
                 */
                const int backup_jitter =
                    quality_rank >= 2
                    ? nmod(
                        5 * op
                        + 7 * j
                        + 3 * machine
                        + a,
                        2
                    )
                    : 0;

                const int p =
                    3
                    + operation_offset
                    + quality_penalty[quality_rank]
                    + backup_jitter;

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}


inline RawInstance make_xlarge_10j10m_very_high_flex_9alt()
{
    RawInstance r;

    r.name =
        "xlarge_10jobs_10machines_9_alternatives_very_high_flex";

    r.num_jobs = 10;
    r.num_machines = 10;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 10;
    constexpr int alternatives = 9;

    /*
     * 10 jobs × 10 operations × 9 alternatives
     * = 900 operation-machine tuples.
     *
     * Relative flexibility:
     * 9 / 10 = 90%.
     *
     * Every operation can be processed on nine
     * out of the ten available machines.
     */
    const std::array<int, alternatives> shifts = {
        0, 1, 2, 3, 4, 5, 6, 8, 9
    };

    /*
     * One fastest alternative,
     * one competitive alternative,
     * and seven progressively slower backups.
     */
    const std::array<int, alternatives> quality_penalty = {
        0, 1, 3, 5, 7, 9, 11, 13, 15
    };

    /*
    * All these values are coprime with 9.
    * Therefore, each one generates a permutation
    * of quality ranks 0,...,8.
    */
    const std::array<int, 6> quality_steps = {
        1, 2, 4, 5, 7, 8
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    /*
     * Average entries per machine:
     * 900 tuples / 10 machines = 90 entries.
     */
    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(
            n * alternatives
        );
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            /*
             * Balanced but not completely symmetric
             * machine eligibility.
             */
            const int base_machine =
                nmod(
                    3 * op
                    + 7 * j
                    + (j * op) % 3,
                    n
                );

            /*
             * Operation-dependent offset in {0,1}.
             * It introduces moderate processing-time
             * asymmetry between operations.
             */
            const int operation_offset =
                nmod(
                    2 * j
                    + 3 * op
                    + (j * op) % 3
                    + (j + op * op) % 5,
                    2
                );

            for (int a = 0; a < alternatives; ++a) {

                /*
                 * The nine shifts are distinct modulo 10.
                 * Therefore, every operation receives exactly
                 * nine distinct eligible machines.
                 */
                const int machine_zero_based =
                    nmod(
                        base_machine + shifts[a],
                        n
                    );

                const int machine =
                    machine_zero_based + 1;

                /*
                 * Every operation receives all nine quality
                 * ranks exactly once.
                 */
                // const int quality_rank =
                //     nmod(
                //         a
                //         + 3 * j
                //         + (j * op) % 3,
                //         alternatives
                //     );

                /*
                * Irregular rotation of the quality ranks.
                */
                const int quality_rotation =
                    nmod(
                        8 * j
                        + 3 * op
                        + 5 * ((j * op) % 3)
                        + 4 * ((j * j + op * op) % 2),
                        alternatives
                    );

                /*
                * Change the order of quality ranks between operations.
                * Since every quality_step is coprime with 9,
                * all ranks from 0 to 8 appear exactly once.
                */
                const int quality_step =
                    quality_steps[
                        nmod(
                            5 * j + 9 * op,
                            static_cast<int>(quality_steps.size())
                        )
                    ];

                const int quality_rank =
                    nmod(
                        quality_rotation
                        + quality_step * a,
                        alternatives
                    );

                /*
                 * The two competitive alternatives are not
                 * perturbed. Jitter is added only to backups.
                 */
                const int backup_jitter =
                    quality_rank >= 2
                    ? nmod(
                        5 * op
                        + 7 * j
                        + 3 * machine
                        + a,
                        2
                    )
                    : 0;

                /*
                 * Processing times range approximately
                 * from 3 to 20.
                 */
                const int p =
                    3
                    + operation_offset
                    + quality_penalty[quality_rank]
                    + backup_jitter;

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}

inline RawInstance make_xlarge_10j10m_variable_6_to_9_flex()
{
    RawInstance r;

    r.name =
        "xlarge_10jobs_10machines_6_to_9_alternatives_cg_branching_friendly";

    r.num_jobs = 10;
    r.num_machines = 10;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 10;

    /*
     * Distribution of flexibility:
     *
     * 40 operations × 6 alternatives
     * 30 operations × 7 alternatives
     * 20 operations × 8 alternatives
     * 10 operations × 9 alternatives
     *
     * Total tuples:
     * 40×6 + 30×7 + 20×8 + 10×9 = 700.
     *
     * Average alternatives per operation = 7.
     */
    // const std::array<int, 9> quality_penalty = {
    //     0, 2, 4, 6, 8, 10, 12, 14, 16
    // };

    const std::array<int, 9> quality_penalty = {
        0, 1, 5, 8, 11, 14, 17, 20, 23
    };

    /*
     * Both steps are coprime with 10.
     * Using only two steps creates controlled asymmetry.
     */
    const std::array<int, 2> machine_steps = {
        3, 7
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(70);
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            const int operation_index =
                j * n + op;

            /*
             * Every block of ten operations contains:
             *
             * 4 operations with 6 alternatives,
             * 3 operations with 7 alternatives,
             * 2 operations with 8 alternatives,
             * 1 operation with 9 alternatives.
             */
            // const int flexibility_bucket =
            //     operation_index % 10;

            const int flexibility_key =
                nmod(
                    37 * operation_index + 17,
                    100
                );

            int num_alternatives;

            if (flexibility_key < 50) {
                num_alternatives = 6;
            }
            else if (flexibility_key < 80) {
                num_alternatives = 7;
            }
            else if (flexibility_key < 95) {
                num_alternatives = 8;
            }
            else {
                num_alternatives = 9;
            }

            /*
             * Balanced cyclic base machine.
             */
            const int base_machine =
                nmod(
                    3 * op
                    + 7 * j
                    + (j * op) % 3,
                    n
                );

            /*
             * Only two possible machine-set orientations.
             */
            const int machine_step =
                machine_steps[
                    nmod(
                        j
                        + 2 * op
                        + (j * op) % 3,
                        static_cast<int>(
                            machine_steps.size()
                        )
                    )
                ];

            /*
             * Moderate processing-time asymmetry in {0,1}.
             */
            const int operation_offset =
                nmod(
                    2 * j
                    + 3 * op
                    + (j * op) % 3
                    + (j + op * op) % 5,
                    2
                );

            /*
             * Simpler quality rotation.
             */
            // const int quality_rotation =
            //     nmod(
            //         3 * j
            //         + 5 * op
            //         + (j * op) % 3,
            //         num_alternatives
            //     );

            const int quality_rotation =
                nmod(
                    2 * j
                    + 3 * op,
                    num_alternatives
                );

            /*
             * Both 1 and 5 are coprime with 6, 7, 8 and 9.
             * Therefore, every operation receives each quality
             * rank exactly once.
             */

            // const int quality_step =
            //     nmod(
            //         j
            //         + 2 * op
            //         + (j * op) % 3,
            //         2
            //     ) == 0
            //         ? 1
            //         : 5;

            for (int a = 0;
                 a < num_alternatives;
                 ++a)
            {
                const int machine_zero_based =
                    nmod(
                        base_machine
                        + machine_step * a,
                        n
                    );

                const int machine =
                    machine_zero_based + 1;

                // const int quality_rank =
                //     nmod(
                //         quality_rotation
                //         + quality_step * a,
                //         num_alternatives
                //     );

                const int quality_rank =
                    nmod(
                        quality_rotation + a,
                        num_alternatives
                    );

                /*
                 * No jitter for the two best alternatives.
                 */
                const int backup_jitter =
                    quality_rank >= 2
                    ? nmod(
                        5 * op
                        + 7 * j
                        + 3 * machine
                        + a,
                        2
                    )
                    : 0;

                const int p =
                    3
                    + operation_offset
                    + quality_penalty[quality_rank]
                    + backup_jitter;

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}


inline RawInstance
make_xlarge_10j10m_variable_6_to_8_flex()
{
    RawInstance r;

    r.name =
        "xlarge_10jobs_10machines_6_to_8_alternatives_cg_branching_friendly";

    r.num_jobs = 10;
    r.num_machines = 10;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 10;
    constexpr int max_alternatives = 8;

    /*
     * Distribution of flexibility:
     *
     * 40 operations × 6 alternatives
     * 35 operations × 7 alternatives
     * 25 operations × 8 alternatives
     *
     * Total tuples:
     *
     * 40×6 + 35×7 + 25×8
     * = 240 + 245 + 200
     * = 685.
     *
     * Average alternatives per operation:
     *
     * 685 / 100 = 6.85.
     */

    /*
     * One fastest alternative,
     * one competitive alternative,
     * and progressively slower backups.
     */
    const std::array<int, max_alternatives>
        quality_penalty = {
            0, 1, 5, 8,
            11, 14, 17, 20
        };

    /*
     * Both steps are coprime with 10.
     *
     * Therefore, the generated machines are distinct
     * for up to eight alternatives.
     */
    const std::array<int, 2> machine_steps = {
        3, 7
    };

    auto nmod = [](int value, int modulus) -> int {
        const int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    /*
     * Average entries per machine:
     *
     * 685 / 10 = 68.5.
     */
    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(70);
    }

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            const int operation_index =
                j * n + op;

            /*
             * Since 37 is coprime with 100,
             * flexibility_key is a permutation of
             * the values 0,...,99.
             *
             * This guarantees exactly:
             *
             * 40 operations with 6 alternatives,
             * 35 operations with 7 alternatives,
             * 25 operations with 8 alternatives.
             */
            const int flexibility_key =
                nmod(
                    37 * operation_index + 17,
                    100
                );

            int num_alternatives;

            if (flexibility_key < 40) {
                num_alternatives = 6;
            }
            else if (flexibility_key < 75) {
                num_alternatives = 7;
            }
            else {
                num_alternatives = 8;
            }

            /*
             * Balanced but not completely symmetric
             * base machine.
             */
            const int base_machine =
                nmod(
                    3 * op
                    + 7 * j
                    + (j * op) % 3,
                    n
                );

            /*
             * Two possible orientations of the
             * eligible-machine set.
             */
            const int machine_step =
                machine_steps[
                    nmod(
                        j
                        + 2 * op
                        + (j * op) % 3,
                        static_cast<int>(
                            machine_steps.size()
                        )
                    )
                ];

            /*
             * Operation-dependent offset in {0,1}.
             */
            const int operation_offset =
                nmod(
                    2 * j
                    + 3 * op
                    + (j * op) % 3
                    + (j + op * op) % 5,
                    2
                );

            /*
             * Rotate the quality ranks between
             * jobs and operations.
             */
            const int quality_rotation =
                nmod(
                    2 * j
                    + 3 * op,
                    num_alternatives
                );

            for (int a = 0;
                 a < num_alternatives;
                 ++a)
            {
                /*
                 * Because machine_step is either 3 or 7,
                 * all generated machines are distinct.
                 */
                const int machine_zero_based =
                    nmod(
                        base_machine
                        + machine_step * a,
                        n
                    );

                const int machine =
                    machine_zero_based + 1;

                /*
                 * Every operation receives each quality rank
                 * from 0 to num_alternatives-1 exactly once.
                 */
                const int quality_rank =
                    nmod(
                        quality_rotation + a,
                        num_alternatives
                    );

                /*
                 * Do not perturb the two best alternatives.
                 */
                const int backup_jitter =
                    quality_rank >= 2
                    ? nmod(
                        5 * op
                        + 7 * j
                        + 3 * machine
                        + a,
                        2
                    )
                    : 0;

                const int p =
                    3
                    + operation_offset
                    + quality_penalty[quality_rank]
                    + backup_jitter;

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}


inline RawInstance make_xlarge_10j10m_low_flex()
{
    RawInstance r;

    r.name =
        "xlarge_10jobs_10machines_4_to_6_alternatives_asymmetric";

    r.num_jobs = 10;
    r.num_machines = 10;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    constexpr int n = 10;
    constexpr int min_alternatives = 4;
    constexpr int max_alternatives = 6;

    /*
     * Increased processing-time levels.
     *
     * Since machine_penalty and jitter are both in {0,1},
     * the resulting intervals are:
     *
     * rank 0:  5--7
     * rank 1:  7--9
     * rank 2: 10--12
     * rank 3: 13--15
     * rank 4: 16--18
     * rank 5: 19--21
     *
     * Therefore:
     *   4-alternative operations have p in 5--15;
     *   5-alternative operations have p in 5--18;
     *   6-alternative operations have p in 5--21.
     */
    const std::array<int, max_alternatives>
        quality_penalty = {
            0, 1, 3, 5, 7, 9
        };

        // quality_penalty = {
        //     0, 2, 5, 8, 11, 14
        // };

    /*
     * Small heterogeneous machine penalty.
     */
    const std::array<int, n> machine_penalty = {
        0, 0, 1, 2, 0,
        1, 0, 2, 1, 0
    };

    /*
     * Deterministic mixing function.
     */
    auto mix = [](std::uint32_t x) -> std::uint32_t {
        x ^= x >> 16;
        x *= 0x7feb352dU;
        x ^= x >> 15;
        x *= 0x846ca68bU;
        x ^= x >> 16;

        return x;
    };

    /*
     * The average number of tuples per machine is 49.9.
     */
    for (int m = 1; m <= n; ++m) {
        r.machine_pattern[m].reserve(
            n * max_alternatives
        );
    }

    /*
     * Number of eligibility tuples already assigned
     * to each machine.
     */
    std::array<int, n> eligibility_count{};

    for (int j = 0; j < n; ++j) {
        for (int op = 0; op < n; ++op) {

            const int operation_id =
                j * n + op;

            /*
             * Permute operation IDs before assigning the
             * number of alternatives.
             *
             * Result:
             *   34 operations with 4 alternatives
             *   33 operations with 5 alternatives
             *   33 operations with 6 alternatives
             */
            const int permuted_operation_id =
                (
                    37 * operation_id
                    + 17
                ) % (n * n);

            // const int operation_alternatives =
            //     min_alternatives
            //     + (
            //         permuted_operation_id
            //         % (
            //             max_alternatives
            //             - min_alternatives
            //             + 1
            //         )
            //     );

            int operation_alternatives;

            // if (permuted_operation_id < 60) {
            //     operation_alternatives = 4;
            // }
            // else if (permuted_operation_id < 90) {
            //     operation_alternatives = 5;
            // }
            // else {
            //     operation_alternatives = 6;
            // }

            if (permuted_operation_id < 50) {
                operation_alternatives = 4;
            }
            else if (permuted_operation_id < 85) {
                operation_alternatives = 5;
            }
            else {
                operation_alternatives = 6;
            }

            /*
             * Prevent duplicate machines within the same
             * operation.
             */
            std::array<bool, n> machine_used{};

            /*
             * Build and irregularly shuffle the quality ranks.
             *
             * An operation with:
             *   4 alternatives uses ranks 0,1,2,3;
             *   5 alternatives uses ranks 0,1,2,3,4;
             *   6 alternatives uses ranks 0,1,2,3,4,5.
             */
            std::array<int, max_alternatives> quality_order = {
                0, 1, 2, 3, 4, 5
            };

            const std::uint32_t quality_seed =
                mix(
                    2654435761U
                        * static_cast<std::uint32_t>(j + 1)
                    ^
                    2246822519U
                        * static_cast<std::uint32_t>(op + 1)
                    ^
                    3266489917U
                        * static_cast<std::uint32_t>(
                            operation_alternatives
                        )
                );

            if ((quality_seed & 1U) != 0U) {
                std::swap(
                    quality_order[0],
                    quality_order[1]
                );
            }

            if (
                operation_alternatives >= 5
                && ((quality_seed >> 1U) & 1U) != 0U
            ) {
                std::swap(
                    quality_order[3],
                    quality_order[4]
                );
            }

            if (
                operation_alternatives == 6
                && ((quality_seed >> 2U) & 1U) != 0U
            ) {
                std::swap(
                    quality_order[4],
                    quality_order[5]
                );
            }

            for (
                int a = 0;
                a < operation_alternatives;
                ++a
            ) {
                /*
                 * Generate an irregular eligibility key.
                 */
                const std::uint32_t eligibility_key =
                    mix(
                        73856093U
                            * static_cast<std::uint32_t>(j + 1)
                        ^
                        19349663U
                            * static_cast<std::uint32_t>(op + 1)
                        ^
                        83492791U
                            * static_cast<std::uint32_t>(a + 1)
                    );

                int machine_zero_based = -1;
                int best_selection_score = 1'000'000'000;

                /*
                 * Select a machine by considering:
                 *
                 * 1. global eligibility balance;
                 * 2. an irregular hash-based tie-break;
                 * 3. distinctness inside the operation.
                 */
                for (
                    int candidate = 0;
                    candidate < n;
                    ++candidate
                ) {
                    if (machine_used[candidate]) {
                        continue;
                    }

                    const std::uint32_t tie_key =
                        mix(
                            eligibility_key
                            ^
                            (
                                2654435761U
                                * static_cast<std::uint32_t>(
                                    candidate + 1
                                )
                            )
                        );

                    const int tie_break =
                        static_cast<int>(
                            tie_key % 997U
                        );

                    /*
                     * The first term guarantees balanced
                     * eligibility counts. The second term
                     * breaks ties asymmetrically.
                     */

                    // const int selection_score =
                    //     1000 * eligibility_count[candidate]
                    //     + tie_break;

                    const int selection_score =
                        25 * eligibility_count[candidate]
                        + tie_break;

                    if (
                        selection_score
                        < best_selection_score
                    ) {
                        best_selection_score =
                            selection_score;

                        machine_zero_based =
                            candidate;
                    }
                }

                machine_used[machine_zero_based] = true;
                ++eligibility_count[machine_zero_based];

                const int machine =
                    machine_zero_based + 1;

                /*
                 * Quality rank assigned irregularly to the
                 * selected machine.
                 */
                const int quality_rank =
                    quality_order[a];

                /*
                 * Tuple-specific processing-time perturbation.
                 */
                const std::uint32_t processing_key =
                    mix(
                        374761393U
                            * static_cast<std::uint32_t>(j + 1)
                        ^
                        668265263U
                            * static_cast<std::uint32_t>(op + 1)
                        ^
                        2246822519U
                            * static_cast<std::uint32_t>(machine)
                        ^
                        3266489917U
                            * static_cast<std::uint32_t>(a + 1)
                    );

                const int jitter =
                    static_cast<int>(
                        processing_key % 2U
                    );

                const int p =
                    5
                    + quality_penalty[quality_rank]
                    + machine_penalty[machine_zero_based]
                    + jitter;

                r.machine_pattern[machine].push_back({
                    op,
                    j,
                    p
                });
            }
        }
    }

    return r;
}

inline RawInstance make_xlarge_9j9m_very_high_flex()
{
    RawInstance r;

    r.name =
        "xlarge_9jobs_9machines_18ops_very_high_flex";

    r.num_jobs = 9;
    r.num_machines = 9;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    // 9 jobs × 18 operations × 8 alternative machines
    // = 1296 operation-machine tuples.
    constexpr int J = 9;
    constexpr int M = 9;
    constexpr int O = 18;
    constexpr int K = 8;

    auto nmod = [](int value, int modulus) -> int {
        int result = value % modulus;

        return result < 0
            ? result + modulus
            : result;
    };

    // Initialise all machine lists.
    for (int m = 1; m <= M; ++m) {
        // Average number of tuples per machine:
        // 1296 / 9 = 144.
        r.machine_pattern[m].reserve(144);
    }

    // Generate every (job, operation) pair explicitly.
    for (int j = 0; j < J; ++j) {
        for (int o = 0; o < O; ++o) {

            // Nonlinear cyclic rotation of the eligible-machine set.
            // Every operation excludes one different machine.
            const int base_machine =
                nmod(
                    4 * j
                    + 5 * o
                    + 2 * ((j * o) % 7)
                    + (j % 3) * (o % 4),
                    M
                );

            // Base duration depends mainly on the operation and job.
            const int base_time =
                4 + nmod(
                    5 * o
                    + 7 * j
                    + 2 * j * o,
                    8
                );

            for (int a = 0; a < K; ++a) {

                // Eight distinct machines modulo 9.
                const int m =
                    1 + nmod(
                        base_machine + a,
                        M
                    );

                // Alternatives have similar processing times.
                // This creates many near-equivalent assignments.
                const int variation =
                    nmod(
                        3 * m
                        + 5 * a
                        + j
                        + o,
                        4
                    );

                // Processing times range approximately from 4 to 14.
                const int p =
                    base_time + variation;

                r.machine_pattern[m].push_back({
                    o,
                    j,
                    p
                });
            }
        }
    }

    return r;
}


inline RawInstance make_xlarge_11j10m_very_high_flex() {
    RawInstance r;

    r.name = "xlarge_11jobs_10machines_very_high_flex";
    r.num_jobs = 11;
    r.num_machines = 10;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    // Rectangular very-high-flexibility instance:
    //
    // 11 jobs × 10 operations/job × 8 alternative machines
    // = 880 operation-machine tuples.
    const int num_jobs = 11;
    const int num_machines = 10;
    const int num_operations = 10;

    // Eight distinct shifts modulo 10.
    // Consequently, every operation is eligible on exactly 8 machines.
    const std::vector<int> shifts = {
        0, 1, 2, 3, 4, 5, 7, 9
    };

    auto machine_mod = [&](int value) -> int {
        int result = value % num_machines;

        return result < 0
            ? result + num_machines
            : result;
    };

    // A different modulus is used for processing times so that jobs J0
    // and J10 do not have identical processing-time structures.
    auto time_mod = [&](int value) -> int {
        int result = value % num_jobs;

        return result < 0
            ? result + num_jobs
            : result;
    };

    for (int m = 1; m <= num_machines; ++m) {
        std::vector<RawEntry> entries;

        // Every machine receives:
        // 11 jobs × 8 alternatives = 88 entries.
        entries.reserve(
            num_jobs *
            static_cast<int>(shifts.size())
        );

        for (int j = 0; j < num_jobs; ++j) {
            for (int a = 0;
                 a < static_cast<int>(shifts.size());
                 ++a)
            {
                // Cyclic machine eligibility.
                // Operation indices range from 0 to 9.
                const int op =
                    machine_mod(
                        (m - 1)
                        - j
                        - shifts[a]
                    );

                // Processing times range from 3 to 13.
                const int p =
                    3 + time_mod(
                        5 * op
                        + 7 * j
                        + 4 * m
                        + 3 * a
                    );

                entries.push_back({
                    op,
                    j,
                    p
                });
            }
        }

        r.machine_pattern[m] =
            std::move(entries);
    }

    return r;
}


inline RawInstance
make_xlarge_14j14m_asymmetric()
{
    RawInstance r;

    r.name =
        "xlarge_14jobs_14machines_asymmetric";

    r.num_jobs = 14;
    r.num_machines = 14;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    // 14 jobs × 14 operations × 8 eligible machines
    // = 196 operations and 1568 operation-machine alternatives.
    const int n = 14;

    const std::vector<int> shifts = {
        0, 1, 2, 3, 5, 7, 9, 11
    };

    auto nmod = [&](int value) -> int {
        int result = value % n;

        return result < 0
            ? result + n
            : result;
    };

    for (int m = 1; m <= n; ++m) {
        std::vector<RawEntry> entries;

        // Every machine receives exactly:
        // 14 jobs × 8 alternatives = 112 entries.
        entries.reserve(
            n * static_cast<int>(shifts.size())
        );

        for (int j = 0; j < n; ++j) {
            for (int a = 0;
                 a < static_cast<int>(shifts.size());
                 ++a)
            {
                // Cyclic eligibility:
                // every operation is eligible on exactly eight machines.
                const int op =
                    nmod(
                        (m - 1)
                        - j
                        - shifts[a]
                    );

                // Balanced processing times with controlled ambiguity.
                //
                // 6 - 5 = 1 is coprime with 14:
                // processing times remain uniform across machines.
                //
                // 5 + 4 = 9 is coprime with 14:
                // processing times remain uniform across job operations.
                //
                // Coefficient 5*a creates more fastest-machine ties than
                // the previous balanced 14×14 variant, making routing
                // decisions harder for standalone CP-SAT.
                const int p =
                    3 +
                    nmod(
                        5 * op
                        + 6 * j
                        + 4 * m
                        + 5 * a
                    );

                entries.push_back({
                    op,
                    j,
                    p
                });
            }
        }

        r.machine_pattern[m] =
            std::move(entries);
    }

    return r;
}


inline RawInstance make_rectangular_cyclic_vhf(
    const std::string& name,
    int num_jobs,
    int num_machines,
    const std::vector<int>& shifts,
    int coef_op,
    int coef_job,
    int coef_machine,
    int coef_alt)
{
    RawInstance r;
    r.name = name;
    r.num_jobs = num_jobs;
    r.num_machines = num_machines;
    r.tuple_order = TupleOrder::OP_JOB_TIME;

    auto nmod = [num_machines](int value) -> int {
        int result = value % num_machines;
        return result < 0 ? result + num_machines : result;
    };

    for (int m = 1; m <= num_machines; ++m) {
        std::vector<RawEntry> entries;
        entries.reserve(
            num_jobs * static_cast<int>(shifts.size()));

        for (int j = 0; j < num_jobs; ++j) {
            for (int a = 0;
                 a < static_cast<int>(shifts.size());
                 ++a)
            {
                int op = nmod(
                    (m - 1)
                    - (j % num_machines)
                    - shifts[a]);

                int p = 3 + nmod(
                    coef_op * op
                    + coef_job * j
                    + coef_machine * m
                    + coef_alt * a);

                entries.push_back({op, j, p});
            }
        }

        r.machine_pattern[m] = std::move(entries);
    }

    return r;
}


inline RawInstance make_rect_10j12m_very_high_flex()
{
    return make_rectangular_cyclic_vhf(
        "rectangular_10jobs_12machines_very_high_flex",
        10,
        12,
        {0, 1, 2, 4, 5, 7, 8, 10},
        4,   // operation coefficient
        3,   // job coefficient
        8,   // machine coefficient
        11   // alternative coefficient
    );
}

inline RawInstance make_rect_13j11m_very_high_flex()
{
    return make_rectangular_cyclic_vhf(
        "rectangular_13jobs_11machines_very_high_flex",
        13,
        11,
        {0, 1, 2, 3, 5, 6, 8, 10},
        4,
        3,
        7,
        10
    );
}

inline RawInstance make_rect_14j12m_very_high_flex()
{
    return make_rectangular_cyclic_vhf(
        "rectangular_14jobs_12machines_very_high_flex",
        14,
        12,
        {0, 1, 2, 4, 5, 7, 8, 10},
        4,
        3,
        8,
        11
    );
}



// LA03-based 10 jobs, 5 machines (job_op_time order)
inline RawInstance make_large_machine_based_10j5m() {
    RawInstance r;
    r.name = "LA03_10jobs_5machines_machine_based";
    r.num_jobs = 10; r.num_machines = 5;
    r.tuple_order = TupleOrder::JOB_OP_TIME;  // (job, op, time)
    r.machine_pattern = {
        {1, {{0,2,82},{0,3,84},{1,0,21},{1,2,18},{1,3,41},{2,2,16},{2,3,52},
             {3,0,37},{3,1,54},{4,0,57},{4,1,81},{5,1,79},{5,2,89},
             {6,2,91},{6,3,20},{6,4,66},{7,2,32},{7,3,55},
             {8,1,7},{8,3,64},{9,0,40},{9,2,19}}},
        {2, {{0,0,23},{1,1,29},{1,2,18},{1,4,50},{2,1,54},{2,4,52},
             {3,3,62},{4,2,61},{4,3,68},{5,0,81},{5,2,89},{5,3,89},
             {6,0,33},{6,1,20},{6,4,66},{7,0,24},{7,1,84},
             {8,1,7},{8,3,64},{8,4,39},{9,1,83},{9,2,19}}},
        {3, {{0,1,45},{0,2,82},{1,0,21},{1,4,50},{2,0,38},{2,3,52},
             {3,2,74},{3,4,57},{4,3,68},{4,4,30},{5,3,89},
             {6,0,33},{6,1,20},{7,3,55},{8,3,64},{9,1,83},{9,3,8}}},
        {4, {{0,4,38},{1,3,41},{1,4,50},{2,1,54},{2,3,52},
             {3,0,37},{3,1,54},{3,2,74},{3,4,57},{4,3,68},
             {5,1,79},{5,4,11},{6,0,33},{7,2,32},{7,4,8},
             {8,2,54},{9,0,40},{9,4,7}}},
        {5, {{0,3,84},{1,2,18},{1,3,41},{2,0,38},{2,2,16},{2,4,52},
             {3,0,37},{3,2,74},{3,4,57},{4,0,57},
             {5,0,81},{5,2,89},{5,3,89},{5,4,11},{6,3,20},
             {7,0,24},{8,0,56},{8,2,54},{8,4,39},{9,0,40}}}
    };
    return r;
}

// ---------------------------------------------------------------------------
// Named instance registry
// ---------------------------------------------------------------------------

// enum class InstanceName {
//     SMALL_6J6M,
//     MEDIUM_7J7M,
//     LARGE_8J8M,
//     XLARGE_9J9M,
//     XLARGE_10J10M,
//     XLARGE_11J11M,
//     XLARGE_11J11M_VHF,
//     XLARGE_12J12M,
//     XLARGE_13J13M_VHF,
//     LARGE_MACHINE_BASED_10J5M,
// };

enum class InstanceName {
    SMALL_6J6M,
    MEDIUM_7J7M,
    LARGE_8J8M,
    XLARGE_9J9M,
    XLARGE_10J10M,
    XLARGE_10J10M_HF,
    XLARGE_10J10M_9ALT,
    XLARGE_10J10M_VAR_6_9,
    XLARGE_10J10M_VAR_6_8,
    XLARGE_11J11M,
    XLARGE_11J11M_6ALT,
    XLARGE_11J11M_10ALT,
    XLARGE_11J11M_VAR_8_10,
    XLARGE_11J11M_9ALT,
    XLARGE_11J11M_VAR_6_8,
    XLARGE_11J11M_VHF,
    XLARGE_12J12M,
    XLARGE_13J13M_VHF,
    XLARGE_14J14M_ASYMMETRIC,

    XLARGE_9J9M_HF,

    XLARGE_10J10M_VHF,
    // XLARGE_10J10M_VHF,

    XLARGE_10J10M_LOW_FLEX,

    // XLARGE_11J11M_VHF,
    XLARGE_11J10M_VHF,
    XLARGE_14J14M_VHF,
    XLARGE_9J9M_VHF,
    XLARGE_15J15M_VHF,

    RECT_10J12M_VHF,
    RECT_13J11M_VHF,
    RECT_14J12M_VHF,

    // New BHF benchmark instances
    // BHF_09J09M_K06,
    BHF_10J10M_K06,
    BHF_11J11M_K07,
    BHF_14J14M_K08,
    BHF_14J14M_K09_DENSE,
    BHF_15J15M_K09,
    BHF_15J15M_K10_DENSE,
    BHF_16J16M_K10,

    LARGE_MACHINE_BASED_10J5M,
};

inline Instance get_instance(InstanceName name) {
    switch (name) {
        case InstanceName::SMALL_6J6M:               return load_raw_instance(make_small_6j6m());
        case InstanceName::MEDIUM_7J7M:              return load_raw_instance(make_medium_7j7m());
        case InstanceName::LARGE_8J8M:               return load_raw_instance(make_large_8j8m());
        case InstanceName::XLARGE_9J9M:              return load_raw_instance(make_xlarge_9j9m());
        case InstanceName::XLARGE_10J10M:            return load_raw_instance(make_xlarge_10j10m());

        // case InstanceName::XLARGE_9J9M_HF:
        //     return load_raw_instance(
        //         make_xlarge_9j9m_high_flex()
        //     );

        case InstanceName::XLARGE_10J10M_VHF:
            return load_raw_instance(
                make_xlarge_10j10m_very_high_flex()
            );

        case InstanceName::XLARGE_10J10M_HF:
            return load_raw_instance(
                make_xlarge_10j10m_high_flex()
            );

        case InstanceName::XLARGE_10J10M_9ALT:
            return load_raw_instance(
                make_xlarge_10j10m_very_high_flex_9alt()
            );

        case InstanceName::XLARGE_10J10M_VAR_6_9:
            return load_raw_instance(
                make_xlarge_10j10m_variable_6_to_9_flex()
            );

        case InstanceName::XLARGE_10J10M_VAR_6_8:
            return load_raw_instance(
                make_xlarge_10j10m_variable_6_to_8_flex()
            );

        case InstanceName::XLARGE_10J10M_LOW_FLEX:
            return load_raw_instance(
                make_xlarge_10j10m_low_flex()
            );

        case InstanceName::XLARGE_11J11M_VAR_6_8:
            return load_raw_instance(
                make_xlarge_11j11m_variable_6_to_8_flex()
            );

        case InstanceName::XLARGE_11J11M_VHF:
            return load_raw_instance(
                make_xlarge_11j11m_very_high_flex()
            );


        case InstanceName::XLARGE_11J11M_6ALT:
            return load_raw_instance(
                make_xlarge_11j11m_high_flex_6alt()
            );


        case InstanceName::XLARGE_11J11M_10ALT:
            return load_raw_instance(
                make_xlarge_11j11m_very_high_flex_10alt()
            );

        case InstanceName::XLARGE_11J11M_VAR_8_10:
            return load_raw_instance(
                make_xlarge_11j11m_variable_8_to_10_flex()
            );

        case InstanceName::XLARGE_11J11M_9ALT:
            return load_raw_instance(
                make_xlarge_11j11m_very_high_flex_9alt()
            );

        case InstanceName::XLARGE_11J10M_VHF:
            return load_raw_instance(
                make_xlarge_11j10m_very_high_flex()
            );

        case InstanceName::XLARGE_14J14M_VHF:
            return load_raw_instance(
                make_xlarge_14j14m_very_high_flex()
            );

        case InstanceName::XLARGE_15J15M_VHF:
            return load_raw_instance(
                make_xlarge_15j15m_very_high_flex()
            );

        case InstanceName::XLARGE_9J9M_VHF:
            return load_raw_instance(
                make_xlarge_9j9m_very_high_flex()
            );

        case InstanceName::XLARGE_11J11M:             return load_raw_instance(make_xlarge_11j11m());
        // case InstanceName::XLARGE_11J11M_VHF:          return load_raw_instance(make_xlarge_11j11m_very_high_flex());
        case InstanceName::XLARGE_12J12M:             return load_raw_instance(make_xlarge_12j12m_very_high_flex());
        case InstanceName::XLARGE_13J13M_VHF:          return load_raw_instance(make_xlarge_13j13m_very_high_flex());
        case InstanceName::XLARGE_14J14M_ASYMMETRIC:   return load_raw_instance(make_xlarge_14j14m_asymmetric());
        case InstanceName::LARGE_MACHINE_BASED_10J5M: return load_raw_instance(make_large_machine_based_10j5m());

        case InstanceName::RECT_10J12M_VHF:
            return load_raw_instance(
                make_rect_10j12m_very_high_flex());

        case InstanceName::RECT_13J11M_VHF:
            return load_raw_instance(
                make_rect_13j11m_very_high_flex());

        case InstanceName::RECT_14J12M_VHF:
            return load_raw_instance(
                make_rect_14j12m_very_high_flex());

        // case InstanceName::BHF_09J09M_K06:
        //     return load_raw_instance(
        //         make_bhf_09j09m_k06());

        case InstanceName::BHF_10J10M_K06:
            return load_raw_instance(
                make_bhf_10j10m_k06());

        case InstanceName::BHF_11J11M_K07:
            return load_raw_instance(
                make_bhf_11j11m_k07());

        case InstanceName::BHF_14J14M_K08:
            return load_raw_instance(
                make_bhf_14j14m_k08());

        case InstanceName::BHF_14J14M_K09_DENSE:
            return load_raw_instance(
                make_bhf_14j14m_k09_dense());

        case InstanceName::BHF_15J15M_K09:
            return load_raw_instance(
                make_bhf_15j15m_k09());

        case InstanceName::BHF_15J15M_K10_DENSE:
            return load_raw_instance(
                make_bhf_15j15m_k10_dense());

        case InstanceName::BHF_16J16M_K10:
            return load_raw_instance(
                make_bhf_16j16m_k10());

        }

    return load_raw_instance(make_medium_7j7m());
}
