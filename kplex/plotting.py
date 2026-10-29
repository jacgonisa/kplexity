"""Plotting for `kplex plot`."""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .model import load_curve, first_rise, curve_from_fit


def plot(csv_path: str, out_png: str, fit: dict | None = None, title: str | None = None) -> str:
    x, y, ycol = load_curve(csv_path)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(x, y, s=14, color="black", zorder=3, label="k-plexity")

    if fit and "error" not in fit:
        xs = np.linspace(x.min(), x.max(), 400)
        if str(fit.get("model", "")).startswith("kplexity"):
            yfit = curve_from_fit(xs, fit)
            B = first_rise(xs, fit["b"], fit["k0"])
            ax.fill_between(xs, yfit, B * fit["A"], color="#EE6677", alpha=0.25, lw=0,
                            label=f"near-copy gap (S={fit['S']:.3f} at k0, divergence {100 * fit['divergence']:.1f}%)")
            ax.plot(xs, B * fit["A"], "--", color="#4477AA", lw=1.2, label="chance collisions x A")
            ax.plot(xs, yfit, color="crimson", lw=2,
                    label=f"kplexity model (R²={fit.get('r_squared', float('nan')):.4f})")
            ax.axhline(fit["A"], color="grey", ls=":", lw=1,
                       label=f"A={fit['A']:.3f} (identical 1-A={fit['identical']:.3f})")
            ax.axvline(fit["k0"], color="grey", ls=":", lw=0.8)

    ax.set_xlabel("k")
    ax.set_ylabel("fraction distinct (U/T)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(title or Path(csv_path).stem)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    return out_png
