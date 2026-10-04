"""Plotting for `kplex plot`."""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .model import double_sigmoid, single_sigmoid, load_curve


def plot(csv_path: str, out_png: str, fit: dict | None = None, title: str | None = None) -> str:
    x, y, ycol = load_curve(csv_path)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(x, y, s=14, color="black", zorder=3, label="k-plexity")

    if fit and "error" not in fit:
        xs = np.linspace(x.min(), x.max(), 400)
        if fit.get("model") == "double_sigmoid":
            p = [fit[k] for k in ("L1", "s1", "k0_1", "L2", "s2", "k0_2")]
            ax.plot(xs, double_sigmoid(xs, *p), color="crimson", lw=2,
                    label=f"double sigmoid (R²={fit.get('r_squared', float('nan')):.3f})")
            # decomposition
            ax.plot(xs, single_sigmoid(xs, p[0], p[1], p[2]), "--", color="#4477AA", lw=1.2,
                    label="component 1 (L1)")
            ax.plot(xs, single_sigmoid(xs, p[3], p[4], p[5]), "--", color="#EE6677", lw=1.2,
                    label="component 2 (L2)")
            ax.axhline(fit["asymptote"], color="grey", ls=":", lw=1,
                       label=f"asymptote={fit['asymptote']:.3f}")
        else:
            ax.plot(xs, single_sigmoid(xs, fit["L"], fit["s"], fit["k0"]), color="crimson", lw=2,
                    label=f"single sigmoid (R²={fit.get('r_squared', float('nan')):.3f})")

    ax.set_xlabel("k")
    ax.set_ylabel("fraction unique (U/T)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(title or Path(csv_path).stem)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    return out_png
