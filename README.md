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

## Implementation

The implementation is written in **C++** and uses:

- **IBM ILOG CPLEX** for the restricted master problem.
- **Google OR-Tools CP-SAT** for constraint programming.

These dependencies must be installed and configured before building and running the project.
