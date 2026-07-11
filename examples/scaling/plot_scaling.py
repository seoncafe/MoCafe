#!/usr/bin/env python3
"""Analyze and plot the MPI scaling benchmark (results.csv from run_scaling.sh).

Produces scaling_strong.pdf (speedup + parallel efficiency, equal-share vs
master-slave) and scaling_weak.pdf (weak-scaling efficiency), and prints the
numbers as a LaTeX-ready table.
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    data = {}
    with open(os.path.join(HERE, "results.csv")) as f:
        for row in csv.DictReader(f):
            key = (row["test"], row["mode"])
            data.setdefault(key, []).append((int(row["np"]), float(row["wall_seconds"])))
    return {k: np.array(sorted(v)) for k, v in data.items()}


def main():
    d = load()

    # ---- strong scaling ---------------------------------------------------
    eq = d[("strong", "eq")]
    ms = d[("strong", "ms")]
    T1 = eq[eq[:, 0] == 1, 1][0]        # equal-share single-rank reference

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.3))
    for arr, lab, c, m in ((eq, "equal share", "#1f77b4", "o"),
                           (ms, "master-slave", "#d62728", "s")):
        npr, t = arr[:, 0], arr[:, 1]
        ax1.loglog(npr, T1 / t, m + "-", color=c, label=lab)
        ax2.semilogx(npr, 100 * T1 / t / npr, m + "-", color=c, label=lab)
    nn = np.array([1, 64])
    ax1.loglog(nn, nn, "k:", lw=1, label="ideal")
    nps = ms[:, 0]
    ax1.loglog(nps, nps - 1, ":", color="#d62728", lw=1,
               label=r"$N_{\rm p}-1$ (worker limit)")
    ax2.semilogx(nps, 100 * (nps - 1) / nps, ":", color="#d62728", lw=1)
    ax1.set_xlabel(r"MPI ranks $N_{\rm p}$"); ax1.set_ylabel(r"speed-up $T_1/T_{N_{\rm p}}$")
    ax2.set_xlabel(r"MPI ranks $N_{\rm p}$"); ax2.set_ylabel("parallel efficiency [\\%]")
    ax2.set_ylim(0, 110); ax2.axhline(100, color="k", ls=":", lw=1)
    for ax in (ax1, ax2):
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "scaling_strong.pdf"))
    print("wrote scaling_strong.pdf")

    # ---- weak scaling ------------------------------------------------------
    fig2, ax = plt.subplots(figsize=(6.0, 4.2))
    for key, lab, c, m in ((("weak", "eq"), "equal share", "#1f77b4", "o"),
                           (("weak", "ms"), "master-slave", "#d62728", "s")):
        if key not in d:
            continue
        arr = d[key]; npr, t = arr[:, 0], arr[:, 1]
        Tref = d[("weak", "eq")][0, 1]           # np=1 reference
        ax.semilogx(npr, 100 * Tref / t, m + "-", color=c, label=lab)
    ax.axhline(100, color="k", ls=":", lw=1)
    ax.set_xlabel(r"MPI ranks $N_{\rm p}$ (work $\propto N_{\rm p}$)")
    ax.set_ylabel("weak-scaling efficiency [\\%]")
    ax.set_ylim(0, 115)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
    ax.legend(frameon=False, fontsize=9)
    fig2.tight_layout()
    fig2.savefig(os.path.join(HERE, "scaling_weak.pdf"))
    print("wrote scaling_weak.pdf")

    # ---- LaTeX table -------------------------------------------------------
    print("\n%--- strong scaling table rows (np & T_eq & S_eq & E_eq & T_ms & S_ms & E_ms)")
    for npv in sorted(set(eq[:, 0]) | set(ms[:, 0])):
        te = eq[eq[:, 0] == npv, 1]
        tm = ms[ms[:, 0] == npv, 1]
        se = f"{te[0]:8.1f} & {T1/te[0]:5.2f} & {100*T1/te[0]/npv:5.1f}" if te.size else "  --   &  --   &  --  "
        sm = f"{tm[0]:8.1f} & {T1/tm[0]:5.2f} & {100*T1/tm[0]/npv:5.1f}" if tm.size else "  --   &  --   &  --  "
        print(f"{int(npv):3d} & {se} & {sm} \\\\")


if __name__ == "__main__":
    main()
