#!/usr/bin/env python3

# TODO regenerate data for J1J2-line with new optimization level settings.
import sys

sys.path.insert(1, "../")
from matplotlib import pyplot as plt
import json
import matplotlib
import numpy as np

matplotlib.use("Agg")  # Use non-interactive backend

plt.rcParams.update({"text.usetex": True, "font.family": "Computer Modern"})

plt.ioff()


def plot_data(transpiler, router, linear, size, paper=False):
    if transpiler[5:9] == "line":
        target = "chain"
        markersize = 2
        if size == "large":
            i_max = 125
        elif size == "small":
            i_max = 25
    elif transpiler[5:11] == "square":
        target = "square"
        markersize = 4
        if size == "large":
            i_max = 16
        elif size == "small":
            i_max = 5

    plt.clf()
    with open(f"{transpiler}_{router}.json", "r") as f:
        lines = json.load(f)

    lines = [[eval(k), v] for k, v in lines.items() if eval(k) <= i_max]
    lines.sort()

    depth = [line[1]["depth"] for line in lines]
    swaps = [line[1]["swaps"] for line in lines]
    quantile_depth = [line[1]["quantile_depth"] for line in lines]
    quantile_swaps = [line[1]["quantile_swaps"] for line in lines]

    if target == "chain":
        wall_clock = [line[1]["wall_clock"] for line in lines]
        quantile_wall_clock = [line[1]["quantile_wall_clock"] for line in lines]
        ns = [4 * line[0] + 2 for line in lines]
    elif target == "square":
        wall_clock = [line[1]["wall_clock"] / 60 for line in lines]
        quantile_wall_clock = [line[1]["quantile_wall_clock"] / 60 for line in lines]
        ns = [(line[0] ** 2) * 4 + 3 * (2 * line[0]) + 1 for line in lines]

    if paper == False:
        fig, ax = plt.subplots(figsize=(4, 3))
    else:
        fig, ax = plt.subplots(figsize=(4, 2.5))
    fig.subplots_adjust(right=0.75)
    fig.set_facecolor("none")

    twin1 = ax.twinx()
    twin2 = ax.twinx()

    twin2.spines.right.set_position(("axes", 1.2))

    print(transpiler, router, linear, size, paper)
    print(max(ns))
    # Plot router
    (p1,) = ax.plot(ns, depth, "C0-o", label="Depth", markersize=markersize, zorder=2)
    (p2,) = twin1.plot(
        ns, swaps, "C1-o", label="Swaps", markersize=markersize, zorder=1
    )
    if target == "chain":
        label = "Wall-clock time (s)"
    elif target == "square":
        label = "Wall-clock time (m)"
    (p3,) = twin2.plot(
        ns, wall_clock, "C2-o", label=label, markersize=markersize, zorder=1
    )

    # Plot QuanTile
    (p4,) = ax.plot(
        ns,
        quantile_depth,
        "C0--",
        label="Depth",
        markersize=markersize,
        dash_capstyle="round",
        zorder=2,
    )
    (p5,) = twin1.plot(
        ns,
        quantile_swaps,
        "C1--",
        label="Swaps",
        markersize=markersize,
        dash_capstyle="round",
        zorder=1,
    )
    (p6,) = twin2.plot(
        ns,
        quantile_wall_clock,
        "C2--",
        label=label,
        markersize=markersize,
        dash_capstyle="round",
        zorder=1,
    )

    twin1.set(ylabel="Swaps")
    twin2.set(ylabel=label)
    if not paper:
        ax.set(xlabel="$N$", ylabel="Depth")
    else:
        ax.set(xlabel="Number of qudits", ylabel="Depth")

    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())
    ax.yaxis.label.set_color(p1.get_color())

    twin1.tick_params(axis="y", colors=p2.get_color())
    twin2.tick_params(axis="y", colors=p3.get_color())
    ax.tick_params(axis="y", colors=p1.get_color())

    ax.set_zorder(3)
    twin1.set_zorder(2)
    twin2.set_zorder(2)
    ax.patch.set_visible(False)

    for p in [p1, p2, p3, p4]:
        p.set_linewidth(1)

    if paper:
        p4.set_linewidth(2.5)

    if not paper:
        ax.legend(
            handles=[p1, p2, p3],
            loc="upper left",
            bbox_to_anchor=(0, 0.95),
            labels=["Depth", "Swaps", "Wall-clock time"],
        )
    else:
        ax.legend(
            handles=[p1, p2, p3],
            loc="upper left",
            labels=["Depth", "Swaps", "Wall-clock time"],
        )

    twin1.set_yscale("log")
    twin2.set_yscale("log")
    ax.set_xscale("log")
    ax.set_yscale("log")
    if paper == True:
        ax.annotate(
            "\\large\\textbf{5}",
            (3.25, 4.3),
            annotation_clip=False,
            color=p1.get_color(),
        )

    mapping = {
        "qiskitAI": "Qiskit AIRouter",
        "qiskit": "Qiskit SabreSwap",
        "cirq": "Cirq",
    }
    title = f"QuanTile vs. {mapping[router]}, J1J2-{target} on {target}"
    ax.set_title(title, x=0.7)
    if not paper:
        fname = f"/Users/tycho/Postdoc/Papers/Optimal qubit routing for periodic circuits/SM/figures/{target}_{router}_{size}.pdf"
        fname = f"/Users/tycho/Postdoc/Papers/Optimal qubit routing for periodic circuits/SM/figures/{target}_{router}_{size}.pdf"
    else:
        fname = f"/Users/tycho/Postdoc/Papers/Optimal qubit routing for periodic circuits/Paper/benchmark.pdf"
    plt.savefig(
        fname,
        bbox_inches="tight",
    )


# for router in ["qiskitAI", "qiskit", "cirq"]:
#    for size in ["large", "small"]:
#        plot_data("J1J2-line(4,1) --> line(4,1)", router, True, size)
#        plot_data("J1J2-square(2,2) --> square(2,2)", router, False, size)

# plot_data("J1J2-line(4,1) --> line(4,1)", "qiskitAI", True, "large", paper=False)

plot_data("J1J2-line(4,1) --> line(4,1)", "qiskitAI", True, "large", True)
