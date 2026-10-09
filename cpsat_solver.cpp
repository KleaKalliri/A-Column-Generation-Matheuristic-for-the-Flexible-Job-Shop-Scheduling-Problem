// =============================================================================
// cpsat_solver.cpp  –  Standalone CP-SAT solver for FJSSP
//
// Solves the Flexible Job Shop Scheduling Problem (minimize makespan)
// using OR-Tools CP-SAT only — no column generation, no branch-and-price.
//
// Usage:
//   cpsat [--small | --medium | --large | --xlarge | --xxlarge | --la03]
//         [--time <seconds>]
// =============================================================================

#include "instances.h"
#include "types.h"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/sat_parameters.pb.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace sat = operations_research::sat;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static void print_instance_info(const Instance& inst) {
    std::cout << "=== FJSSP Instance: " << inst.name << " ===\n";
    std::cout << "  Jobs:     " << inst.num_jobs     << "\n";
    std::cout << "  Machines: " << inst.num_machines  << "\n";
    std::cout << "  Ops/job:  ";
    int total = 0;
    for (int j = 0; j < inst.num_jobs; ++j) {
        if (j) std::cout << ", ";
        std::cout << inst.num_ops[j];
        total += inst.num_ops[j];
    }
    std::cout << "  (total=" << total << ")\n";
}

static int compute_horizon(const Instance& inst) {
    int h = 0;
    for (int j = 0; j < inst.num_jobs; ++j)
        for (int o = 0; o < inst.num_ops[j]; ++o)
            for (auto& [m, pt] : inst.proc_time[j][o])
                h += pt;
    return h;
}

static void print_schedule(const Schedule& sched, double makespan) {
    std::cout << "\n=== Solution  Cmax = " << makespan << " ===\n";
    std::map<int, std::vector<std::tuple<int,int,int,int>>> by_mach;
    for (auto& [jop, mse] : sched) {
        auto [j, o]    = jop;
        auto [m, s, e] = mse;
        by_mach[m].emplace_back(j, o, s, e);
    }
    for (auto& [m, ops] : by_mach) {
        std::sort(ops.begin(), ops.end(),
                  [](const auto& a, const auto& b) { return std::get<2>(a) < std::get<2>(b); });
        std::cout << "  M" << m << ":";
        for (auto& [j, o, s, e] : ops)
            std::cout << "  J" << j << "o" << o << "[" << s << "," << e << ")";
        std::cout << "\n";
    }
}

// ---------------------------------------------------------------------------
// CP-SAT model
// ---------------------------------------------------------------------------

struct CPSatResult {
    bool     found    = false;
    bool     optimal  = false;
    double   makespan = -1.0;
    Schedule schedule;
};

static CPSatResult solve(const Instance& inst, double time_limit_sec) {
    CPSatResult result;
    const int J = inst.num_jobs;
    const int horizon = compute_horizon(inst);

    sat::CpModelBuilder model;

    // ------------------------------------------------------------------
    // Variables
    // ------------------------------------------------------------------
    // For each (job j, op o):
    //   start[j][o]   : IntVar in [earliest_start, horizon]
    //   end[j][o]     : IntVar in [0, horizon]
    //   presence[j][o][m] : BoolVar — is op (j,o) assigned to machine m?
    //   interval[j][o][m] : optional IntervalVar on machine m

    // Precompute earliest start (critical-path lower bound per op)
    auto earliest_start = [&](int j, int o) -> int {
        int lb = 0;
        for (int k = 0; k < o; ++k) {
            int mn = std::numeric_limits<int>::max();
            for (auto& [m, pt] : inst.proc_time[j][k]) mn = std::min(mn, pt);
            lb += mn;
        }
        return lb;
    };

    using Key3 = std::tuple<int,int,int>;   // (j, o, m)
    using Key2 = std::pair<int,int>;        // (j, o)

    std::map<Key3, sat::IntervalVar> intervals;
    std::map<Key3, sat::BoolVar>     presence;
    std::map<Key2, sat::IntVar>      start_vars;
    std::map<Key2, sat::IntVar>      end_vars;

    for (int j = 0; j < J; ++j) {
        for (int o = 0; o < inst.num_ops[j]; ++o) {
            int slb = earliest_start(j, o);

            auto sv = model.NewIntVar({slb, horizon});
            auto ev = model.NewIntVar({slb, horizon});
            start_vars[{j,o}] = sv;
            end_vars[{j,o}]   = ev;

            const auto& eligs = inst.elig_mach[j][o];
            std::vector<sat::BoolVar> mach_pres;
            sat::LinearExpr pt_sum;

            for (int m : eligs) {
                int pt = inst.proc_time[j][o].at(m);
                auto pv = model.NewBoolVar();
                auto iv = model.NewOptionalFixedSizeIntervalVar(sv, (int64_t)pt, pv);
                intervals[{j,o,m}] = iv;
                presence[{j,o,m}]  = pv;
                mach_pres.push_back(pv);
                pt_sum += sat::LinearExpr::Term(pv, (int64_t)pt);
            }

            // Exactly one machine is chosen
            model.AddExactlyOne(mach_pres);

            // end = start + processing_time_of_chosen_machine
            model.AddEquality(sat::LinearExpr(ev) - sat::LinearExpr(sv), pt_sum);
        }
    }

    // ------------------------------------------------------------------
    // No-overlap constraints per machine
    // ------------------------------------------------------------------
    for (int m : inst.all_machines) {
        std::vector<sat::IntervalVar> mach_ivs;
        for (int j = 0; j < J; ++j)
            for (int o = 0; o < inst.num_ops[j]; ++o)
                if (intervals.count({j,o,m}))
                    mach_ivs.push_back(intervals.at({j,o,m}));
        if (!mach_ivs.empty())
            model.AddNoOverlap(mach_ivs);
    }

    // ------------------------------------------------------------------
    // Job precedence: op o must finish before op o+1 starts
    // ------------------------------------------------------------------
    for (int j = 0; j < J; ++j)
        for (int o = 0; o < inst.num_ops[j] - 1; ++o)
            model.AddLessOrEqual(end_vars[{j,o}], start_vars[{j,o+1}]);

    // ------------------------------------------------------------------
    // Objective: minimize makespan
    // ------------------------------------------------------------------
    sat::IntVar cmax = model.NewIntVar({0, horizon});
    for (int j = 0; j < J; ++j) {
        int last = inst.num_ops[j] - 1;
        model.AddLessOrEqual(end_vars[{j, last}], cmax);
    }
    model.Minimize(cmax);

    // ------------------------------------------------------------------
    // Solve
    // ------------------------------------------------------------------
    sat::SatParameters params;
    if (time_limit_sec > 0.0)
        params.set_max_time_in_seconds(time_limit_sec);
    params.set_log_search_progress(true);

    auto t0 = std::chrono::steady_clock::now();
    sat::CpSolverResponse response = sat::SolveWithParameters(model.Build(), params);
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::cout << "\nSolver status : ";
    switch (response.status()) {
        case sat::CpSolverStatus::OPTIMAL:   std::cout << "OPTIMAL\n";   break;
        case sat::CpSolverStatus::FEASIBLE:  std::cout << "FEASIBLE\n";  break;
        case sat::CpSolverStatus::INFEASIBLE:std::cout << "INFEASIBLE\n";break;
        default:                              std::cout << "UNKNOWN\n";   break;
    }
    std::cout << "Wall time     : " << std::fixed << std::setprecision(2) << elapsed << " s\n";

    if (response.status() == sat::CpSolverStatus::OPTIMAL ||
        response.status() == sat::CpSolverStatus::FEASIBLE)
    {
        result.found   = true;
        result.optimal = (response.status() == sat::CpSolverStatus::OPTIMAL);
        result.makespan = (double)sat::SolutionIntegerValue(response, cmax);

        for (int j = 0; j < J; ++j) {
            for (int o = 0; o < inst.num_ops[j]; ++o) {
                int64_t sv = sat::SolutionIntegerValue(response, start_vars[{j,o}]);
                for (int m : inst.elig_mach[j][o]) {
                    if (!presence.count({j,o,m})) continue;
                    if (sat::SolutionBooleanValue(response, presence[{j,o,m}])) {
                        int pt = inst.proc_time[j][o].at(m);
                        result.schedule[{j,o}] = {m, (int)sv, (int)(sv + pt)};
                        break;
                    }
                }
            }
        }
    }

    return result;
}

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
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    InstanceName inst_name = InstanceName::MEDIUM_7J7M;
    double       time_limit = 0.0;
    std::string instance_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if      (arg == "--small"   || arg == "-s")   inst_name  = InstanceName::SMALL_6J6M;
        else if (arg == "--medium"  || arg == "-m")   inst_name  = InstanceName::MEDIUM_7J7M;
        else if (arg == "--large"   || arg == "-l")   inst_name  = InstanceName::LARGE_8J8M;
        // else if (arg == "--xlarge"  || arg == "-xl")  inst_name  = InstanceName::XLARGE_9J9M;
        else if (arg == "--xxlarge"   || arg == "-xxl")   inst_name = InstanceName::XLARGE_10J10M;
        else if (arg == "--xxxlarge"  || arg == "-xxxl")  inst_name = InstanceName::XLARGE_11J11M;
        else if (arg == "--xxxlarge-vhf")                   inst_name = InstanceName::XLARGE_11J11M_VHF;
        else if (arg == "--xxxxlarge" || arg == "-xxxxl") inst_name = InstanceName::XLARGE_12J12M;
        else if (arg == "--xxxxlarge-vhf")                  inst_name = InstanceName::XLARGE_13J13M_VHF;
        else if (arg == "--xlarge14-asym" || arg == "--14j14m-asym") inst_name = InstanceName::XLARGE_14J14M_ASYMMETRIC;

        else if (arg == "--10j12m-vhf")
            inst_name = InstanceName::RECT_10J12M_VHF;

        else if (arg == "--13j11m-vhf")
            inst_name = InstanceName::RECT_13J11M_VHF;

        else if (arg == "--14j12m-vhf")
            inst_name = InstanceName::RECT_14J12M_VHF;

        // else if (arg == "--14j12m-vhf")
        //     inst_name = InstanceName::RECT_14J12M_VHF;

        else if (arg == "--10j10m-vhf")
            inst_name = InstanceName::XLARGE_10J10M_VHF;

        else if (arg == "--10j10m-hf")
            inst_name = InstanceName::XLARGE_10J10M_HF;

        else if (arg == "--10j10m-9alt")
            inst_name = InstanceName::XLARGE_10J10M_9ALT;

        else if (arg == "--10j10m-var-6-9")
            inst_name =
                InstanceName::XLARGE_10J10M_VAR_6_9;

        else if (arg == "--10j10m-low-flex")
            inst_name = InstanceName::XLARGE_10J10M_LOW_FLEX;

        else if (arg == "--11j11m-var-6-8") {
            inst_name =
                InstanceName::XLARGE_11J11M_VAR_6_8;
        }

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


        else if (arg == "--9j9m-hf") {
            inst_name =
                InstanceName::XLARGE_9J9M_HF;
        }

        else if (arg == "--14j14m-vhf")
            inst_name = InstanceName::XLARGE_14J14M_VHF;

        else if (arg == "--15j15m-vhf")
            inst_name = InstanceName::XLARGE_15J15M_VHF;

        else if (arg == "--9j9m-vhf")
            inst_name = InstanceName::XLARGE_9J9M_VHF;

        // New BHF instances
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

        else if (arg == "--la03")                     inst_name  = InstanceName::LARGE_MACHINE_BASED_10J5M;
        else if ((arg == "--time"   || arg == "-t") && i + 1 < argc)
            time_limit = std::stod(argv[++i]);
        else if (arg == "--help" || arg == "-h") {
            std::cout <<
                "Usage: cpsat [--small|--medium|--large|--xlarge|--xxlarge|--xxxlarge|--xxxxlarge|--la03]\n"
                "             [--time <seconds>]\n"
                "\n"
                "  --small      : 6 jobs, 6 machines\n"
                "  --medium     : 7 jobs, 7 machines (default)\n"
                "  --large      : 8 jobs, 8 machines\n"
                "  --xlarge     : 9 jobs, 9 machines\n"
                "  --xxlarge    : 10 jobs, 10 machines\n"
                "  --xxxlarge   : 11 jobs, 11 machines\n"
                "  --xxxxlarge      : 12 jobs, 12 machines (very high flex)\n"
                "  --xxxxlarge-vhf  : 13 jobs, 13 machines (very high flex)\n"
                "  --la03           : LA03, 10 jobs, 5 machines\n"
                "  --time <s>   : time limit in seconds (default 60)\n";
            std::cout << "  --bhf-09j09m-k06       : 9 jobs, 9 machines, 6 alternatives/op\n";
            std::cout << "  --bhf-10j10m-k06       : 10 jobs, 10 machines, 6 alternatives/op\n";
            std::cout << "  --bhf-11j11m-k07       : 11 jobs, 11 machines, 7 alternatives/op\n";
            std::cout << "  --bhf-14j14m-k08       : 14 jobs, 14 machines, 8 alternatives/op\n";
            std::cout << "  --bhf-14j14m-k09-dense : 14 jobs, 14 machines, 9 alternatives/op\n";
            std::cout << "  --bhf-15j15m-k09       : 15 jobs, 15 machines, 9 alternatives/op\n";
            std::cout << "  --bhf-15j15m-k10-dense : 15 jobs, 15 machines, 10 alternatives/op\n";
            std::cout << "  --bhf-16j16m-k10       : 16 jobs, 16 machines, 10 alternatives/op\n";
            return 0;
        }
    }

    // Instance inst = get_instance(inst_name);

    Instance inst = instance_path.empty()
        ? get_instance(inst_name)
        : load_fjs_instance(instance_path);

    print_instance_info(inst);
    std::cout << "  Horizon:  " << compute_horizon(inst) << "\n";
    if (time_limit > 0.0)
        std::cout << "  Time limit: " << time_limit << " s\n\n" << std::flush;
    else
        std::cout << "  Time limit: none (running until optimal)\n\n" << std::flush;

    CPSatResult res = solve(inst, time_limit);

    if (res.found) {
        print_schedule(res.schedule, res.makespan);
        std::cout << "\nMakespan: " << res.makespan
                  << (res.optimal ? "  (OPTIMAL)" : "  (FEASIBLE)") << "\n";
    } else {
        std::cout << "No solution found.\n";
    }

    return res.found ? 0 : 1;
}
