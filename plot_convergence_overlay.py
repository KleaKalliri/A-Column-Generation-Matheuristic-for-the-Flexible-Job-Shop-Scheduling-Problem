import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "csv_file",
        help="Convergence CSV produced by the Branch-and-Price solver."
    )
    parser.add_argument(
        "--x",
        choices=["node", "time"],
        default="node",
        help="Use processed nodes or elapsed seconds on x-axis."
    )
    parser.add_argument(
        "--out",
        default="convergence.svg",
        help="Output image path."
    )

    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)

    required = {
        "variant",
        "node_index",
        "elapsed_sec",
        "best_lb",
        "ub",
    }

    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"{args.csv_file} is missing columns: {sorted(missing)}"
        )

    df = df.sort_values("node_index").reset_index(drop=True)

    # Make convergence curves monotone.
    df["best_lb"] = df["best_lb"].cummax()
    df["ub"] = df["ub"].cummin()

    variant = str(df["variant"].iloc[0])

    if args.x == "time":
        x = df["elapsed_sec"]
        x_label = "Time [s]"
    else:
        x = df["node_index"]
        x_label = "Processed nodes"

    plt.figure(figsize=(10, 6))

    plt.step(
        x,
        df["best_lb"],
        where="post",
        linewidth=2,
        label=f"{variant} LB"
    )

    plt.step(
        x,
        df["ub"],
        where="post",
        linewidth=2,
        linestyle="--",
        label=f"{variant} UB"
    )

    plt.xlabel(x_label)
    plt.ylabel("Cmax bound")
    plt.title("Branch-and-Price convergence")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(args.out, bbox_inches="tight", dpi=250)

    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()