// k-plexity curve fitting in the browser (and in node, for tests).
// Double sigmoid  y = L1/(1+e^(-s1(k-k0_1))) + L2/(1+e^(-s2(k-k0_2)))  and a single sigmoid, fitted by
// bounded Levenberg-Marquardt from several starts — same starts, bounds and constraints as the fitter used for the
// DToL reference curves (checked against it by test/fit_check.mjs). Unevenly spaced k are weighted by k-span.

export const PARAMS = ["L1", "s1", "k0_1", "L2", "s2", "k0_2"];

const sig = (x, s, k0) => 1 / (1 + Math.exp(-s * (x - k0)));

export function doubleSigmoid(x, p) {
  return p[0] * sig(x, p[1], p[2]) + p[3] * sig(x, p[4], p[5]);
}
export function singleSigmoid(x, p) {
  return p[0] * sig(x, p[1], p[2]);
}

// analytic Jacobian rows for a sum of sigmoids; p = [L, s, k0, L, s, k0, ...]
function jacRow(x, p) {
  const row = new Array(p.length);
  for (let j = 0; j < p.length; j += 3) {
    const [L, s, k0] = [p[j], p[j + 1], p[j + 2]];
    const g = sig(x, s, k0), d = g * (1 - g);
    row[j] = g;
    row[j + 1] = L * d * (x - k0);
    row[j + 2] = -L * d * s;
  }
  return row;
}
function model(x, p) {
  let y = 0;
  for (let j = 0; j < p.length; j += 3) y += p[j] * sig(x, p[j + 1], p[j + 2]);
  return y;
}

function solve(A, b) {                       // Gaussian elimination with partial pivoting
  const n = b.length, M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    [M[c], M[piv]] = [M[piv], M[c]];
    if (Math.abs(M[c][c]) < 1e-300) return null;
    for (let r = c + 1; r < n; r++) {
      const f = M[r][c] / M[c][c];
      for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
    }
  }
  const out = new Array(n);
  for (let r = n - 1; r >= 0; r--) {
    let s = M[r][n];
    for (let k = r + 1; k < n; k++) s -= M[r][k] * out[k];
    out[r] = s / M[r][r];
  }
  return out;
}

const clip = (p, lo, hi) => p.map((v, i) => Math.min(hi[i], Math.max(lo[i], v)));

function rss(x, y, p, w) {
  let s = 0;
  for (let i = 0; i < x.length; i++) { const r = y[i] - model(x[i], p); s += w[i] * r * r; }
  return s;
}

// k-span each point stands for (all 1 for a step-1 curve); weights for sketched, unevenly spaced curves
export function spacingWeights(x) {
  const n = x.length, w = new Array(n);
  for (let i = 0; i < n; i++) {
    const left = i === 0 ? x[0] : (x[i] + x[i - 1]) / 2, right = i === n - 1 ? x[n - 1] : (x[i] + x[i + 1]) / 2;
    w[i] = right - left + (i === 0 ? 0.5 : 0) + (i === n - 1 ? 0.5 : 0);
  }
  return w;
}
const isUneven = x => x.length > 2 && x.slice(1).some((v, i) => Math.abs((v - x[i]) - (x[1] - x[0])) > 1e-9);

// bounded LM: steps are projected onto the box; lambda adapts on accept/reject
export function levenbergMarquardt(x, y, p0, lo, hi, w = null, maxIter = 2000) {
  w = w || x.map(() => 1);
  let p = clip(p0, lo, hi), f = rss(x, y, p, w), lambda = 1e-3;
  const n = p.length;
  for (let it = 0; it < maxIter; it++) {
    const JtJ = Array.from({ length: n }, () => new Array(n).fill(0)), Jtr = new Array(n).fill(0);
    for (let i = 0; i < x.length; i++) {
      const J = jacRow(x[i], p), r = y[i] - model(x[i], p);
      for (let a = 0; a < n; a++) {
        Jtr[a] += w[i] * J[a] * r;
        for (let b = a; b < n; b++) JtJ[a][b] += w[i] * J[a] * J[b];
      }
    }
    for (let a = 0; a < n; a++) for (let b = 0; b < a; b++) JtJ[a][b] = JtJ[b][a];
    let improved = false;
    for (let tries = 0; tries < 30; tries++) {
      const A = JtJ.map((row, a) => row.map((v, b) => (a === b ? v + lambda * (v + 1e-12) : v)));
      const step = solve(A, Jtr);
      if (!step) { lambda *= 10; continue; }
      const pn = clip(p.map((v, a) => v + step[a]), lo, hi), fn = rss(x, y, pn, w);
      if (fn < f) {
        const rel = (f - fn) / Math.max(f, 1e-300);
        p = pn; f = fn; lambda = Math.max(lambda / 10, 1e-12); improved = true;
        if (rel < 1e-12) return { p, rss: f };
        break;
      }
      lambda *= 10;
      if (lambda > 1e12) break;
    }
    if (!improved) break;
  }
  return { p, rss: f };
}

// deterministic PRNG (mulberry32)
function rng(seed) {
  return () => {
    seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const uni = (r, a, b) => a + (b - a) * r();

function stats(x, y, f, p, k) {
  const n = x.length;
  let rs = 0, mean = y.reduce((a, b) => a + b, 0) / n, tss = 0;
  for (let i = 0; i < n; i++) { const r = y[i] - f(x[i], p); rs += r * r; tss += (y[i] - mean) ** 2; }
  const aicc = n * Math.log(rs / n) + 2 * k + (2 * k * (k + 1)) / Math.max(n - k - 1, 1);
  return { rss: rs, R2: 1 - rs / tss, AICc: aicc, RMSE: Math.sqrt(rs / n) };
}

export function fitDoubleSigmoid(x, y, nStarts = 10, weights = undefined) {
  const w = weights === undefined ? (isUneven(x) ? spacingWeights(x) : null) : weights;
  const xmin = Math.min(...x), xmax = Math.max(...x), xr = xmax - xmin, ymax = Math.max(...y);
  const lo = [0.05, 0.005, xmin, 0.05, 0.005, xmin], hi = [0.95, 5.0, xmax, 0.95, 5.0, xmax];
  const r = rng(42);
  let best = null;
  for (let i = 0; i < nStarts; i++) {
    let p0;
    if (i === 0) {
      const i1 = y.findIndex(v => v > 0.3 * ymax), i2 = y.findIndex(v => v > 0.8 * ymax);
      p0 = [ymax * 0.55, 0.3, i1 >= 0 ? x[i1] : xmin + xr * 0.33,
            ymax * 0.40, 0.1, i2 >= 0 ? x[i2] : xmin + xr * 0.67];
    } else {
      const L1 = uni(r, 0.1, 0.7), L2 = uni(r, 0.1, Math.min(0.7, 1 - L1));
      p0 = [L1, uni(r, 0.05, 1), uni(r, xmin, xmin + 0.6 * xr), L2, uni(r, 0.05, 1), uni(r, xmin + 0.4 * xr, xmax)];
    }
    const res = levenbergMarquardt(x, y, p0, lo, hi, w);
    if (res.p[0] + res.p[3] > 1.05) continue;
    if (!best || res.rss < best.rss) best = res;
  }
  if (!best) return null;
  let p = best.p;
  if (p[2] > p[5]) p = [p[3], p[4], p[5], p[0], p[1], p[2]];          // identifiability: k0_1 < k0_2
  const out = Object.fromEntries(PARAMS.map((n, i) => [n, p[i]]));
  out.p = p;
  out.asymptote = p[0] + p[3];
  out.L2_frac = p[3] / (p[0] + p[3]);
  out.weighted = !!w;
  Object.assign(out, stats(x, y, doubleSigmoid, p, 6));
  return out;
}

export function fitSingleSigmoid(x, y, nStarts = 6) {
  const w = isUneven(x) ? spacingWeights(x) : null;
  const xmin = Math.min(...x), xmax = Math.max(...x), ymax = Math.max(...y);
  const lo = [0.01, 0.005, xmin], hi = [1.0, 5.0, xmax];
  const r = rng(7);
  let best = null;
  for (let i = 0; i < nStarts; i++) {
    const i1 = y.findIndex(v => v > 0.5 * ymax);
    const p0 = i === 0 ? [ymax, 0.3, i1 >= 0 ? x[i1] : (xmin + xmax) / 2]
                       : [uni(r, 0.3, 1), uni(r, 0.05, 1), uni(r, xmin, xmin + 0.6 * (xmax - xmin))];
    const res = levenbergMarquardt(x, y, p0, lo, hi, w);
    if (!best || res.rss < best.rss) best = res;
  }
  const p = best.p;
  return { L: p[0], s: p[1], k0: p[2], p, ...stats(x, y, singleSigmoid, p, 3) };
}
