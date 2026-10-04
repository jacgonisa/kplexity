"""Double-sigmoid model and fitting for `kplex fit`.

k-plexity curve (fraction unique vs k) is modelled as a sum of two logistics:
    y(k) = L1/(1+exp(-s1*(k-k0_1))) + L2/(1+exp(-s2*(k-k0_2)))

The six parameters describe the curve's shape:
    asymptote = L1 + L2            -> the high-k plateau (overall uniqueness)
    L2_fraction = L2/(L1+L2)       -> relative height of the second component
    k0_1, k0_2                     -> the two midpoints (where each component rises)
    s1, s2                         -> the steepness of each rise
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def single_sigmoid(x, L, s, k0):
    return L / (1.0 + np.exp(-s * (x - k0)))


def double_sigmoid(x, L1, s1, k0_1, L2, s2, k0_2):
    return (L1 / (1.0 + np.exp(-s1 * (x - k0_1)))
            + L2 / (1.0 + np.exp(-s2 * (x - k0_2))))


def load_curve(csv_path: str) -> tuple[np.ndarray, np.ndarray, str]:
    """Load (k, y). Prefer the theoretical fraction when available."""
    df = pd.read_csv(csv_path)
    ycol = "fraction_unique_theoretical" if "fraction_unique_theoretical" in df.columns else "fraction_unique"
    df = df[["k", ycol]].apply(pd.to_numeric, errors="coerce").dropna()
    df = df[(df[ycol] >= 0) & (df[ycol] <= 1)].sort_values("k")
    return df["k"].to_numpy(float), df[ycol].to_numpy(float), ycol


def _r_squared(y, yhat):
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return np.nan if ss_tot == 0 else 1.0 - ss_res / ss_tot


def _aicc(residuals, n_params, n):
    rss = float(np.sum(residuals ** 2))
    if rss <= 0 or n <= n_params + 2:
        return np.nan
    p = n_params + 1  # include sigma^2
    aic = n * np.log(rss / n) + 2 * p
    return aic + (2 * p * (p + 1)) / (n - p - 1)


def fit_double_sigmoid(x: np.ndarray, y: np.ndarray, n_starts: int = 12) -> dict:
    """Multi-start bounded fit. Enforces L1+L2<=1 region and k0_1<k0_2 (post-hoc)."""
    xmin, xmax = float(x.min()), float(x.max())
    # bounds: [L1, s1, k0_1, L2, s2, k0_2]
    lb = [0.0, 1e-3, xmin, 0.0, 1e-3, xmin]
    ub = [1.0, 5.0, xmax, 1.0, 5.0, xmax]
    rng = np.random.default_rng(0)
    best = None
    for i in range(n_starts):
        frac = 0.3 + 0.4 * (i / max(n_starts - 1, 1))
        p0 = [0.5, 0.2, xmin + frac * (xmax - xmin) * 0.5,
              0.4, 0.1, xmin + (0.4 + 0.5 * frac) * (xmax - xmin)]
        if i:  # jitter subsequent starts
            p0 = [np.clip(v * (1 + 0.3 * rng.standard_normal()), lo, hi)
                  for v, lo, hi in zip(p0, lb, ub)]
        try:
            popt, pcov = curve_fit(double_sigmoid, x, y, p0=p0, bounds=(lb, ub), maxfev=20000)
        except Exception:
            continue
        resid = y - double_sigmoid(x, *popt)
        rss = float(np.sum(resid ** 2))
        if best is None or rss < best[0]:
            best = (rss, popt, pcov)
    if best is None:
        return {"error": "fit failed", "model": "double_sigmoid"}

    _, popt, pcov = best
    L1, s1, k0_1, L2, s2, k0_2 = popt
    # identifiability: order components by k0 (component 2 = later/larger)
    if k0_1 > k0_2:
        L1, s1, k0_1, L2, s2, k0_2 = L2, s2, k0_2, L1, s1, k0_1
    yhat = double_sigmoid(x, L1, s1, k0_1, L2, s2, k0_2)
    resid = y - yhat
    se = np.sqrt(np.diag(pcov)) if pcov is not None and np.all(np.isfinite(pcov)) else [np.nan] * 6
    asymptote = float(L1 + L2)
    return {
        "model": "double_sigmoid",
        "L1": float(L1), "s1": float(s1), "k0_1": float(k0_1),
        "L2": float(L2), "s2": float(s2), "k0_2": float(k0_2),
        "L1_SE": float(se[0]), "L2_SE": float(se[3]),
        "asymptote": asymptote,
        "L2_fraction": float(L2 / asymptote) if asymptote > 0 else np.nan,
        "r_squared": _r_squared(y, yhat),
        "rmse": float(np.sqrt(np.mean(resid ** 2))),
        "aicc": _aicc(resid, 6, len(x)),
        "n_points": int(len(x)),
    }


def fit_single_sigmoid(x: np.ndarray, y: np.ndarray) -> dict:
    xmin, xmax = float(x.min()), float(x.max())
    try:
        popt, pcov = curve_fit(
            single_sigmoid, x, y,
            p0=[float(np.max(y)), 0.1, float(np.median(x))],
            bounds=([0.0, 1e-3, xmin], [1.0, 5.0, xmax]), maxfev=20000)
    except Exception as e:
        return {"error": str(e), "model": "single_sigmoid"}
    yhat = single_sigmoid(x, *popt)
    resid = y - yhat
    return {
        "model": "single_sigmoid",
        "L": float(popt[0]), "s": float(popt[1]), "k0": float(popt[2]),
        "asymptote": float(popt[0]), "r_squared": _r_squared(y, yhat),
        "rmse": float(np.sqrt(np.mean(resid ** 2))), "aicc": _aicc(resid, 3, len(x)),
        "n_points": int(len(x)),
    }
