"""Curve models and fitting for `kplex fit`.

The k-plexity curve y(k) = distinct canonical k-mers / k-mer positions is fitted by default with the
kplexity model (5 parameters):

    y(k) = B(k) * [A - S * exp(-lam * (k - k0))]
    B(k) = (1 - exp(-L)) / L,   L = 2 * exp(-b * (k - k0))

  B(k)   chance collisions: the exact curve of random sequence when b = ln 4 and k0 = log4(G)
  k0     first-rise position (about log4 of the genome size; an internal control)
  b      steepness of the first rise
  A      resolvable fraction: the plateau; 1 - A is the identical-copy fraction (copies no k <= kmax separates)
  S      near-copy share at k0: k-mers in diverged repeat copies that are still shared at k0
  lam    decay rate of that share; divergence between copies = 1 - exp(-lam)

A k-mer in a repeat copy becomes distinct once it spans a difference from the other copies, with probability
1 - (1 - divergence)^k; this gives the exponential. The gamma variant (`--model gamma`, one more parameter a)
lets the divergence vary across copies: exp(-lam (k - k0)) -> ((1 + lam k / a) / (1 + lam k0 / a))^-a.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


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


def spacing_weights(x: np.ndarray) -> np.ndarray:
    """k-span each point stands for (all 1 for a step-1 curve). For an unevenly sampled curve (a k list such as
    a size-anchored sketch), weighting by these makes the fit target the same objective as the full curve."""
    x = np.asarray(x, float)
    e = np.concatenate([[x[0]], (x[1:] + x[:-1]) / 2, [x[-1]]])
    w = np.diff(e)
    w[0] += 0.5
    w[-1] += 0.5
    return w


# ---------------------------------------------------------------------------------------------- kplexity model
KPLEX_PARAMS = ["b", "k0", "A", "R", "lam"]   # internal: y = B * A * (1 - R exp(-lam k)); S = A R exp(-lam k0)
_LO = np.array([0.05, 3, 0, 0, 1e-4])
_HI = np.array([5, 60, 1, 1, 2])
LN4 = float(np.log(4))


def first_rise(k, b, k0):
    """Chance-collision term B(k)."""
    L = 2 * np.exp(np.clip(-b * (np.asarray(k, float) - k0), -50, 50))
    return -np.expm1(-L) / L


def random_curve(k, genome_bp):
    """Expected y(k) of random sequence with genome_bp k-mer positions (no free parameters)."""
    L = 2 * genome_bp / 4.0 ** np.asarray(k, float)
    return -np.expm1(-L) / L


def _decay(k, lam, a=None):
    k = np.asarray(k, float)
    return np.exp(-lam * k) if a is None else (1 + lam * k / a) ** (-a)


def kplexity_curve(k, b, k0, A, S, lam, a=None):
    """The kplexity model in its reported form (S = near-copy share at k0)."""
    rel = _decay(k, lam, a) / _decay(k0, lam, a)
    return first_rise(k, b, k0) * (A - S * rel)


def _internal(k, p, gamma):
    b, k0, A, R, lam = p[:5]
    return first_rise(k, b, k0) * A * (1 - R * _decay(k, lam, p[5] if gamma else None))


def fit_kplexity(x: np.ndarray, y: np.ndarray, genome_bp: float | None = None, gamma: bool = False,
                 weights: np.ndarray | None = None, n_starts: int = 14, seed: int = 0) -> dict:
    """Multi-start bounded least squares fit of the kplexity model.
    genome_bp (total k-mer positions) seeds k0 at log4(G); without it k0 starts where y passes half its maximum.
    Unevenly spaced k are weighted by spacing_weights() unless weights are given."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    name = "kplexity_gamma" if gamma else "kplexity"
    uneven = len(x) > 2 and np.ptp(np.diff(x)) > 0
    if weights is None:
        weights = spacing_weights(x) if uneven else np.ones(len(x))
    sw = np.sqrt(np.asarray(weights, float))
    lo, hi = (np.append(_LO, 0.05), np.append(_HI, 200)) if gamma else (_LO, _HI)
    k0g = np.log(genome_bp) / LN4 if genome_bp else float(x[np.argmax(y > 0.5 * y.max())])
    A0 = min(float(y[-1]) + 0.02, 1.0)
    rng = np.random.default_rng(seed)
    best = None
    for i in range(n_starts):
        if i < 8:
            p0 = [[1.4, 0.8, 2.0, 0.5][i % 4], k0g, A0, [0.3, 0.1, 0.5, 0.2][(i // 2) % 4],
                  [0.02, 0.05, 0.1, 0.01, 0.2, 0.005][i % 6]] + ([[1, 3, 0.5][i % 3]] if gamma else [])
        else:
            p0 = lo + (hi - lo) * rng.uniform(0.02, 0.5, len(lo))
        p0 = np.clip(np.asarray(p0, float), lo + 1e-9, hi - 1e-9)
        try:
            r = least_squares(lambda p: (_internal(x, p, gamma) - y) * sw, p0, bounds=(lo, hi), max_nfev=5000)
        except Exception:
            continue
        if np.isfinite(r.cost) and (best is None or r.cost < best.cost):
            best = r
    if best is None:
        return {"error": "fit failed", "model": name}
    p = best.x
    b, k0, A, R, lam = (float(v) for v in p[:5])
    a = float(p[5]) if gamma else None
    S = A * R * float(_decay(k0, lam, a))
    yhat = _internal(x, p, gamma)
    resid = y - yhat
    out = {"model": name, "A": A, "identical": 1 - A, "S": S, "lam": lam,
           "divergence": float(1 - np.exp(-lam)), "k0": k0, "b": b}
    if gamma:
        out["a"] = a
    gap_end = float(S * _decay(x.max(), lam, a) / _decay(k0, lam, a))
    out.update({
        "gap_at_kmax": gap_end,
        "plateau_reached": bool(gap_end < 0.02),
        "r_squared": _r_squared(y, yhat),
        "rmse": float(np.sqrt(np.mean(resid ** 2))),
        "aicc": float(_aicc(resid, len(p), len(x))),
        "n_points": int(len(x)),
        "weighted": bool(uneven),
        "genome_bp": float(genome_bp) if genome_bp else None,
        "params_internal": dict(zip(KPLEX_PARAMS + (["a"] if gamma else []), map(float, p))),
    })
    return out


def curve_from_fit(x: np.ndarray, fit: dict) -> np.ndarray:
    """Evaluate any fitted model at k = x."""
    m = fit.get("model")
    if m in ("kplexity", "kplexity_gamma"):
        return kplexity_curve(x, fit["b"], fit["k0"], fit["A"], fit["S"], fit["lam"], fit.get("a"))
    raise ValueError(f"unknown model {m!r}")


def genome_bp_from_csv(csv_path: str) -> float | None:
    """Total k-mer positions at the smallest k (+ k - 1), from a Kplex curve CSV; None if not recorded."""
    df = pd.read_csv(csv_path)
    for col in ("total_kmers_theoretical", "total_kmers"):
        if col in df.columns:
            d = df[["k", col]].apply(pd.to_numeric, errors="coerce").dropna().sort_values("k")
            if len(d) and d[col].iloc[0] > 0:
                return float(d[col].iloc[0] + d.k.iloc[0] - 1)
    return None


def fit_curve(x: np.ndarray, y: np.ndarray, model: str = "kplexity", genome_bp: float | None = None) -> dict:
    """Dispatch: model in {kplexity, gamma}."""
    if model == "kplexity":
        return fit_kplexity(x, y, genome_bp)
    if model == "gamma":
        return fit_kplexity(x, y, genome_bp, gamma=True)
    raise ValueError(f"unknown model {model!r}")
