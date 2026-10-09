// =============================================================================
// main.cpp  –  Entry point for FJSSP Branch-and-Price solver
//
// Equivalent to running branch_and_price_machine_based.py.
// Default instance: MEDIUM_7J7M (mirrors the Python default).
// =============================================================================

#include "branch_and_price.h"
#include "instances.h"
#include "types.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <limits>

// ---------------------------------------------------------------------------
// Print instance summary
// ---------------------------------------------------------------------------
static void print_instance_info(const Instance& inst) {
    std::cout << "=== FJSSP Instance: " << inst.name << " ===\n";
    std::cout << "  Jobs:     " << inst.num_jobs     << "\n";
    std::cout << "  Machines: " << inst.num_machines  << "\n";
    std::cout << "  Ops:      ";
    int total = 0;
    for (int j = 0; j < inst.num_jobs; ++j) {
        total += inst.num_ops[j];
        if (j > 0) std::cout << ", ";
        std::cout << inst.num_ops[j];
    }
    std::cout << "  (total=" << total << ")\n";
}

// ---------------------------------------------------------------------------
// Print final schedule
// ---------------------------------------------------------------------------
static void print_schedule(const Schedule& sched, double makespan) {
    std::cout << "\n=== Best Schedule (Cmax=" << makespan << ") ===\n";
    // Group by machine
    std::map<int, std::vector<std::tuple<int,int,int,int>>> by_mach;  // m -> [(j,o,s,e)]
    for (auto& [jop, mse] : sched) {
        auto [j, o]    = jop;
        auto [m, s, e] = mse;
        by_mach[m].push_back({j, o, s, e});
    }
    for (auto& [m, ops] : by_mach) {
        std::sort(ops.begin(), ops.end(),
                  [](const auto& a, const auto& b) {
                      return std::get<2>(a) < std::get<2>(b);
                  });
        std::cout << "  M" << m << ":";
        for (auto& [j, o, s, e] : ops)
            std::cout << " J" << j << "o" << o << "[" << s << "," << e << ")";
        std::cout << "\n";
    }
}

// ---------------------------------------------------------------------------
// Read a Brandimarte / standard FJSP .fjs or .txt instance file
//
// Supported header formats (auto-detected from first line):
//   3-token: num_jobs num_machines average_flexibility  — machines are 1-based (.fjs)
//   2-token: num_jobs num_machines                      — machines are 0-based (.txt)
// ---------------------------------------------------------------------------
static Instance load_fjs_instance(const std::string& file_path)
{
    std::ifstream input(file_path);
    if (!input.is_open())
        throw std::runtime_error("Cannot open FJSP instance file: " + file_path);

    Instance inst;

    std::string header_line;
    while (std::getline(input, header_line))
        if (!header_line.empty()) break;

    std::istringstream header_ss(header_line);
    std::vector<int> header_tokens;
    int tok;
    while (header_ss >> tok) header_tokens.push_back(tok);

    if (header_tokens.size() < 2)
        throw std::runtime_error("Invalid FJS header in file: " + file_path);

    inst.num_jobs     = header_tokens[0];
    inst.num_machines = header_tokens[1];
    const bool one_based_machines = (header_tokens.size() >= 3);

    if (inst.num_jobs <= 0 || inst.num_machines <= 0)
        throw std::runtime_error("Invalid number of jobs or machines in file: " + file_path);

    const std::size_t last_sep = file_path.find_last_of("/\\");
    inst.name = (last_sep == std::string::npos) ? file_path : file_path.substr(last_sep + 1);

    inst.num_ops.resize(inst.num_jobs);
    inst.elig_mach.resize(inst.num_jobs);
    inst.proc_time.resize(inst.num_jobs);
    inst.all_machines.resize(inst.num_machines);
    std::iota(inst.all_machines.begin(), inst.all_machines.end(), 0);

    for (int j = 0; j < inst.num_jobs; ++j) {
        int num_ops = 0;
        if (!(input >> num_ops))
            throw std::runtime_error("Cannot read op count for job " + std::to_string(j));
        if (num_ops <= 0)
            throw std::runtime_error("Job " + std::to_string(j) + " has invalid op count.");

        inst.num_ops[j] = num_ops;
        inst.elig_mach[j].resize(num_ops);
        inst.proc_time[j].resize(num_ops);

        for (int o = 0; o < num_ops; ++o) {
            int num_alts = 0;
            if (!(input >> num_alts))
                throw std::runtime_error("Cannot read alt count for op (" + std::to_string(j) + "," + std::to_string(o) + ").");
            if (num_alts <= 0)
                throw std::runtime_error("Op (" + std::to_string(j) + "," + std::to_string(o) + ") has no eligible machine.");

            auto& elig = inst.elig_mach[j][o];
            auto& pt   = inst.proc_time[j][o];
            elig.reserve(num_alts);

            for (int a = 0; a < num_alts; ++a) {
                int file_mach = 0, ptime = 0;
                if (!(input >> file_mach >> ptime))
                    throw std::runtime_error("Cannot read alt " + std::to_string(a) + " for op (" + std::to_string(j) + "," + std::to_string(o) + ").");

                const int m = one_based_machines ? file_mach - 1 : file_mach;
                if (m < 0 || m >= inst.num_machines)
                    throw std::runtime_error("Invalid machine ID " + std::to_string(file_mach) + " for op (" + std::to_string(j) + "," + std::to_string(o) + ").");
                if (ptime <= 0)
                    throw std::runtime_error("Invalid processing time for op (" + std::to_string(j) + "," + std::to_string(o) + ").");
                if (pt.count(m))
                    throw std::runtime_error("Duplicate machine " + std::to_string(file_mach) + " for op (" + std::to_string(j) + "," + std::to_string(o) + ").");

                elig.push_back(m);
                pt[m] = ptime;
            }
            std::sort(elig.begin(), elig.end());
        }
    }
    return inst;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Parse instance name from command line (optional)
    InstanceName inst_name = InstanceName::MEDIUM_7J7M;
    std::string instance_path;
    // std::string rl_log_path;
    // std::string rl_policy_path;

    bool save_convergence = false;
    bool enable_node_cp_sat = true;
    std::string convergence_log_path = "convergence.csv";
    std::string variant_label = "variant";

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--small"   || arg == "-s")  inst_name = InstanceName::SMALL_6J6M;
        else if (arg == "--medium"  || arg == "-m")  inst_name = InstanceName::MEDIUM_7J7M;
        else if (arg == "--large"   || arg == "-l")  inst_name = InstanceName::LARGE_8J8M;
        else if (arg == "--xlarge"  || arg == "-xl") inst_name = InstanceName::XLARGE_9J9M;
        else if (arg == "--xxlarge"  || arg == "-xxl") inst_name = InstanceName::XLARGE_10J10M;

        else if (arg == "--11j11m-vhf") {
            inst_name = InstanceName::XLARGE_11J11M_VHF;
        }

        else if (arg == "--11j11m-6alt") {
            inst_name =
                InstanceName::XLARGE_11J11M_6ALT;
        }

        else if (arg == "--11j11m-10alt") {
            inst_name =
                InstanceName::XLARGE_11J11M_10ALT;
        }

        else if (arg == "--11j11m-var-8-10") {
            inst_name =
                InstanceName::XLARGE_11J11M_VAR_8_10;
        }

        else if (arg == "--11j11m-9alt") {
            inst_name =
                InstanceName::XLARGE_11J11M_9ALT;
        }

        else if (arg == "--11j11m-var-6-8") {
            inst_name =
                InstanceName::XLARGE_11J11M_VAR_6_8;
        }

        else if (arg == "--10j10m-vhf")
            inst_name = InstanceName::XLARGE_10J10M_VHF;

        else if (arg == "--10j10m-hf")
            inst_name = InstanceName::XLARGE_10J10M_HF;

        else if (arg == "--10j10m-9alt")
            inst_name = InstanceName::XLARGE_10J10M_9ALT;

        else if (arg == "--10j10m-var-6-9")
            inst_name =
                InstanceName::XLARGE_10J10M_VAR_6_9;

        else if (arg == "--10j10m-var-6-8")
            inst_name =
                InstanceName::XLARGE_10J10M_VAR_6_8;

        else if (arg == "--10j10m-low-flex")
            inst_name = InstanceName::XLARGE_10J10M_LOW_FLEX;

        else if (arg == "--9j9m-hf") {
            inst_name =
                InstanceName::XLARGE_9J9M_HF;
        }

        else if (arg == "--11j10m-vhf")
            inst_name = InstanceName::XLARGE_11J10M_VHF;

        else if (arg == "--14j14m-vhf")
            inst_name = InstanceName::XLARGE_14J14M_VHF;

        else if (arg == "--15j15m-vhf")
            inst_name = InstanceName::XLARGE_15J15M_VHF;

        else if (arg == "--9j9m-vhf")
            inst_name = InstanceName::XLARGE_9J9M_VHF;

        else if (arg == "--xxxlarge" || arg == "-xxxl")inst_name = InstanceName::XLARGE_11J11M;
        else if (arg == "--xxxlarge-vhf")               inst_name = InstanceName::XLARGE_11J11M_VHF;
        else if (arg == "--xxxxlarge"|| arg == "-xxxxl")inst_name = InstanceName::XLARGE_12J12M;
        else if (arg == "--xxxxlarge-vhf")              inst_name = InstanceName::XLARGE_13J13M_VHF;
        else if (arg == "--xlarge14-asym" || arg == "--14j14m-asym")   inst_name = InstanceName::XLARGE_14J14M_ASYMMETRIC;
        else if (arg == "--la03")                    inst_name = InstanceName::LARGE_MACHINE_BASED_10J5M;
        else if (arg == "--10j12m-vhf")
            inst_name = InstanceName::RECT_10J12M_VHF;

        else if (arg == "--13j11m-vhf")
            inst_name = InstanceName::RECT_13J11M_VHF;

        else if (arg == "--14j12m-vhf")
            inst_name = InstanceName::RECT_14J12M_VHF;

        // else if (arg == "--bhf-09j09m-k06")
        //     inst_name = InstanceName::BHF_09J09M_K06;

        else if (arg == "--bhf-10j10m-k06")
            inst_name = InstanceName::BHF_10J10M_K06;

        else if (arg == "--bhf-11j11m-k07")
            inst_name = InstanceName::BHF_11J11M_K07;

        else if (arg == "--11j10m-vhf")
            inst_name = InstanceName::XLARGE_11J10M_VHF;

        else if (arg == "--14j14m-vhf")
            inst_name = InstanceName::XLARGE_14J14M_VHF;

        else if (arg == "--15j15m-vhf")
            inst_name = InstanceName::XLARGE_15J15M_VHF;

        else if (arg == "--9j9m-vhf")
            inst_name = InstanceName::XLARGE_9J9M_VHF;

        else if (arg == "--xxxlarge" || arg == "-xxxl")inst_name = InstanceName::XLARGE_11J11M;
        else if (arg == "--xxxlarge-vhf")               inst_name = InstanceName::XLARGE_11J11M_VHF;
        else if (arg == "--xxxxlarge"|| arg == "-xxxxl")inst_name = InstanceName::XLARGE_12J12M;
        else if (arg == "--xxxxlarge-vhf")              inst_name = InstanceName::XLARGE_13J13M_VHF;
        else if (arg == "--xlarge14-asym" || arg == "--14j14m-asym")   inst_name = InstanceName::XLARGE_14J14M_ASYMMETRIC;
        else if (arg == "--la03")                    inst_name = InstanceName::LARGE_MACHINE_BASED_10J5M;
        else if (arg == "--10j12m-vhf")
            inst_name = InstanceName::RECT_10J12M_VHF;

        else if (arg == "--13j11m-vhf")
            inst_name = InstanceName::RECT_13J11M_VHF;

        else if (arg == "--14j12m-vhf")
            inst_name = InstanceName::RECT_14J12M_VHF;

        // else if (arg == "--bhf-09j09m-k06")
        //     inst_name = InstanceName::BHF_09J09M_K06;

        else if (arg == "--bhf-10j10m-k06")
            inst_name = InstanceName::BHF_10J10M_K06;

        else if (arg == "--bhf-11j11m-k07")
            inst_name = InstanceName::BHF_11J11M_K07;

        else if (arg == "--bhf-14j14m-k08")
            inst_name = InstanceName::BHF_14J14M_K08;

        else if (arg == "--bhf-14j14m-k09-dense")
            inst_name = InstanceName::BHF_14J14M_K09_DENSE;

        else if (arg == "--bhf-15j15m-k09")
            inst_name = InstanceName::BHF_15J15M_K09;

        else if (arg == "--bhf-15j15m-k10-dense")
            inst_name = InstanceName::BHF_15J15M_K10_DENSE;

        else if (arg == "--bhf-16j16m-k10")
            inst_name = InstanceName::BHF_16J16M_K10;

        else if ((arg == "--file" || arg == "-f") && i + 1 < argc) {
            instance_path = argv[++i];
        }
        // else if ((arg == "--rl-log") && i + 1 < argc) {
        //     rl_log_path = argv[++i];  // next token is the CSV path
        // }
        // else if ((arg == "--rl-policy") && i + 1 < argc) {
        //     rl_policy_path = argv[++i];  // trained weights file
        // }
        else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: bap [--small | --medium | --large | --xlarge | --xxlarge | --10j10m-vhf | --la03]\n";
            std::cout << "           [--file <path.fjs|path.txt>]\n";
            // std::cout << "           [--rl-log <path.csv>]\n";
            std::cout << "  --small   : 6 jobs, 6 machines\n";
            std::cout << "  --medium  : 7 jobs, 7 machines (default)\n";
            std::cout << "  --large   : 8 jobs, 8 machines\n";
            std::cout << "  --xlarge  : 9 jobs, 9 machines\n";
            std::cout << "  --xxlarge : 10 jobs, 10 machines\n";
            std::cout << "  --la03    : LA03 10 jobs, 5 machines\n";
            std::cout << "  --file    : load any .fjs (1-based) or .txt (0-based) instance\n";
            std::cout << "  --bhf-09j09m-k06       : 9 jobs, 9 machines, 6 alternatives/op\n";
            std::cout << "  --bhf-10j10m-k06       : 10 jobs, 10 machines, 6 alternatives/op\n";
            std::cout << "  --bhf-11j11m-k07       : 11 jobs, 11 machines, 7 alternatives/op\n";
            std::cout << "  --bhf-14j14m-k08       : 14 jobs, 14 machines, 8 alternatives/op\n";
            std::cout << "  --bhf-14j14m-k09-dense : 14 jobs, 14 machines, 9 alternatives/op\n";
            std::cout << "  --bhf-15j15m-k09       : 15 jobs, 15 machines, 9 alternatives/op\n";
            std::cout << "  --bhf-15j15m-k10-dense : 15 jobs, 15 machines, 10 alternatives/op\n";
            std::cout << "  --bhf-16j16m-k10       : 16 jobs, 16 machines, 10 alternatives/op\n";
            std::cout
                << "  --10j10m-vhf : 10 jobs, 10 machines, "
                << "8 alternative machines per operation\n";
            std::cout
                << "  --disable-node-cpsat : disable CP-SAT calls after CG at B&P nodes; "
                << "keep initial-UB CP-SAT enabled\n";
            // std::cout << "  --rl-log    : append RL training rows to CSV file\n";
            // std::cout << "  --rl-policy : load trained RL policy weights file\n";
            return 0;
        }
        else if (arg == "--save-convergence") {
            save_convergence = true;
        }
        else if (arg == "--conv-log" && i + 1 < argc) {
            save_convergence = true;
            convergence_log_path = argv[++i];
        }
        else if (arg == "--variant-label" && i + 1 < argc) {
            variant_label = argv[++i];
        }
        else if (arg == "--disable-node-cpsat" ||
                 arg == "--no-node-cpsat") {
            enable_node_cp_sat = false;
        }
    }

    Instance inst = instance_path.empty()
        ? get_instance(inst_name)
        : load_fjs_instance(instance_path);

    // ── Run CP-SAT first to get real UB and tighten horizon ──────────────────
    // Optimization mode: minimize Cmax within CP_SAT_INITIAL_TL seconds.
    // The resulting makespan becomes the initial UB; T = UB - 1 is the B&P horizon.
    int    horizon = compute_horizon(inst);   // pessimistic fallback
    std::cout << "[Init] Computing initial upper bound via CP-SAT (" << CP_SAT_INITIAL_TL << "s)...\n";
    std::cout << std::flush;
    BranchAndPriceNode cp_root;
    auto cp_res = cp_sat_solve(inst, cp_root, CP_SAT_INITIAL_TL);
    double init_ub = std::numeric_limits<double>::infinity();
    if (cp_res.found) {
        init_ub = cp_res.makespan;
        horizon = static_cast<int>(init_ub);  // CG horizon = UB, not UB-1
    }

    // ── Print instance header (with real UB and tight horizon) ───────────────
    print_instance_info(inst);
    if (std::isfinite(init_ub))
        std::cout << "  UB (CP-SAT): " << init_ub << "\n";
    else
        std::cout << "  UB (CP-SAT): not found\n";
    std::cout << "  Horizon:     " << horizon << "\n\n";
    std::cout << std::flush;


    std::cout << "=== Branch-and-Price ===\n";

    BranchAndPriceSolver bap(inst);
    bap.set_node_cp_sat_enabled(enable_node_cp_sat);

    if (cp_res.found) {
        bap.set_initial_ub(cp_res.makespan, cp_res.schedule);
    }

    // Convergence logging: LB/UB per processed node.
    // Disabled by default unless --save-convergence or --conv-log is used.
    std::ofstream convergence_log;

    double best_lb_curve =
        -std::numeric_limits<double>::infinity();

    auto convergence_t0 =
        std::chrono::steady_clock::now();

    if (save_convergence) {
        convergence_log.open(convergence_log_path);

        if (!convergence_log.is_open()) {
            throw std::runtime_error(
                "Cannot open convergence log file: " +
                convergence_log_path
            );
        }

        convergence_log
            << "variant,node_index,node_id,elapsed_sec,node_lb,best_lb,ub,gap_pct\n";

        int convergence_node_index = 0;

        bap.on_node_cb =
            [&, convergence_node_index](int node_id, double node_lb, double ub) mutable
            {
                ++convergence_node_index;

                const auto now =
                    std::chrono::steady_clock::now();

                const double elapsed_sec =
                    std::chrono::duration<double>(
                        now - convergence_t0
                    ).count();

                // if (std::isfinite(node_lb)) {
                //     best_lb_curve =
                //         std::max(best_lb_curve, node_lb);
                // }

                if (std::isfinite(node_lb)) {
                    best_lb_curve = node_lb;
                }

                const double best_lb_to_write =
                    std::isfinite(best_lb_curve)
                        ? best_lb_curve
                        : node_lb;

                double gap_pct = 0.0;

                if (std::isfinite(ub) &&
                    std::isfinite(best_lb_to_write) &&
                    ub > 1e-9)
                {
                    // gap_pct =
                    //     (ub - best_lb_to_write) / ub * 100.0;

                    gap_pct =
                        std::max(
                            0.0,
                            (ub - best_lb_to_write) / ub * 100.0
                        );
                }

                convergence_log
                    << variant_label << ","
                    << convergence_node_index << ","
                    << node_id << ","
                    << std::setprecision(10) << elapsed_sec << ","
                    << node_lb << ","
                    << best_lb_to_write << ","
                    << ub << ","
                    << gap_pct
                    << "\n";

                convergence_log.flush();
            };

        std::cout
            << "[Convergence] Logging enabled: "
            << convergence_log_path
            << " label=" << variant_label
            << "\n";
    }

    auto res = bap.solve();

    if (save_convergence) {
        convergence_log.close();
    }

    std::cout << "\n=== Results ===\n";
    std::cout << "  Node-level CP-SAT:       "
              << (enable_node_cp_sat ? "enabled" : "disabled") << "\n";
    std::cout << "  B&P time limit (s):      "
              << BAP_GLOBAL_TIME_LIMIT << "\n";
    std::cout << "  Cmax (makespan):        " << res.best_obj << "\n";
    std::cout << "  Best UB:                " << res.best_obj << "\n";
    std::cout << "  Best LB:                " << res.best_bound << "\n";
    {
        const double gap = (std::isfinite(res.best_obj) &&
                            std::isfinite(res.best_bound) &&
                            res.best_obj > 1e-9)
            ? (res.best_obj - res.best_bound) / res.best_obj * 100.0
            : std::numeric_limits<double>::quiet_NaN();
        std::cout << "  Gap:            "
                  << std::fixed << std::setprecision(2);
        if (std::isfinite(gap)) std::cout << gap << "%\n";
        else                    std::cout << "N/A\n";
    }
    // std::cout << "  Proven optimal: " << (res.proven_optimal ? "YES" : "NO") << "\n";
    std::cout << "  Time (s):                 " << res.elapsed_seconds << "\n";
    std::cout << "  Nodes explored:           " << res.nodes_explored << "\n";
    std::cout << "  Nodes pruned by bound:    " << res.pruned_bound << "\n";
    std::cout << "  Pruned infeasible:        " << res.pruned_infeasible << "\n";
    std::cout << "  Total Cols:               " << res.total_cols_generated << "\n";
    std::cout << "  T_pricing (wall s):       " << res.total_pricing_time << "\n";
    std::cout << "  T_master (s):             " << res.total_master_time << "\n";
    std::cout << "  Avg P/M (pricing):   "
              << std::fixed << std::setprecision(1)
              << res.avg_beam_frac * 100.0 << "%\n";

    if (!res.best_schedule.empty())
        print_schedule(res.best_schedule, res.best_obj);

    return 0;
}
