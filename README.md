# A Column Generation Matheuristic for the Flexible Job-Shop Scheduling Problem

This repository contains the implementation accompanying the research paper **“A Column Generation Matheuristic for the Flexible Job-Shop Scheduling Problem.”**

The project addresses the Flexible Job-Shop Scheduling Problem (FJSSP), where each operation can be processed on one of several eligible machines. The objective is to minimize the **makespan**, i.e., the completion time of the last operation.

## Method

The approach combines:

- **Column generation**, using machine patterns that represent ordered and timed sequences of operations.
- **Bounded beam-search pricing**, to generate patterns with negative reduced cost.
- **A branch-and-price-style tree search**, with operation–machine assignment and machine-sequence branching.
- **Constraint programming (CP-SAT)**, to construct and improve feasible schedules using neighborhoods derived from the master solution.

Additional components include propagation-based bounds, dual stabilization, parallel pricing, and column inheritance.

Because pricing is heuristic and CP-SAT calls are time-limited, the method **does not certify optimality**.

## Computational Results

The method was evaluated on nine synthetic asymmetric instances containing 100–225 operations. It matched or improved the best makespan obtained by standalone CP-SAT within a two-hour limit, requiring less computational time on every tested instance.

The table below compares the proposed matheuristic with standalone CP-SAT on nine synthetic asymmetric instances. The objective is to minimize makespan; lower values indicate better schedules.

| Instance | Operations | CP-SAT Makespan | Matheuristic Makespan | CP-SAT Time (min) | Matheuristic Time (min) |
|---|---:|---:|---:|---:|---:|
| I1 | 100 | 48* | **48** | 10.66 | 3.14 |
| I2 | 121 | 53 | **52** | 120.00† | 73.74 |
| I3 | 144 | 56 | **56** | 120.00† | 24.93 |
| I4 | 169 | 63 | **63** | 120.00† | 68.56 |
| I5 | 196 | 71 | **70** | 120.00† | 80.58 |
| I6 | 225 | 71 | **71** | 120.00† | 117.86 |
| I7 | 100 | 45* | **45** | 52.91 | 5.19 |
| I8 | 100 | 45 | **45** | 120.00† | 1.63 |
| I9 | 100 | 46* | **46** | 36.00 | 8.58 |

**Notes:**

- All computation times are reported in minutes and represent the running times reported in the paper.
- * indicates optimality proven by standalone CP-SAT.
- † indicates that standalone CP-SAT reached the two-hour time limit without proving optimality.
- The proposed matheuristic uses heuristic pricing and time-limited CP-SAT calls and does not itself certify optimality.

The matheuristic matched CP-SAT's makespan on seven instances and improved it on two instances (I2 and I5), requiring less computational time on every tested instance. It also reached the optimal makespan established by CP-SAT on I1, I7, and I9.

## Implementation

The implementation is written in **C++** and uses:

- **IBM ILOG CPLEX** for the restricted master problem.
- **Google OR-Tools CP-SAT** for constraint programming.

These dependencies must be installed and configured before building and running the project.
