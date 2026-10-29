// k-plexity curve fitting in the browser (and in node, for tests).
// kplexity model (default)  y = B(k) [A - S e^(-lam (k - k0))],  B = (1 - e^-L)/L,  L = 2 e^(-b (k - k0))
//   fitted as in kplex/model.py (fit_kplexity): internal form y = B A (1 - R e^(-lam k)), same bounds and 14 starts
//   (8 grid + the same 6 random starts), checked against it by test/fit_check_kplexity.mjs.
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

function rss(x, y, p, w, f = model) {
  let s = 0;
  for (let i = 0; i < x.length; i++) { const r = y[i] - f(x[i], p); s += w[i] * r * r; }
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

// numeric Jacobian row (central differences), for models without an analytic one
export function numericJac(fn, rel = 1e-6) {
  return (x, p) => p.map((v, j) => {
    const h = rel * Math.max(Math.abs(v), 1e-3), pp = [...p], pm = [...p];
    pp[j] = v + h; pm[j] = v - h;
    return (fn(x, pp) - fn(x, pm)) / (2 * h);
  });
}

// bounded LM: steps are projected onto the box; lambda adapts on accept/reject.
// opts.f(x, p) is the model and opts.jac(x, p) its Jacobian row (default: the sum of sigmoids, analytic);
// opts.freeze: parameters sitting on a bound whose gradient points out of the box are held fixed for that step
// (an active set, so a fit that ends on a bound, e.g. A = 1, converges instead of stalling on clipped steps).
export function levenbergMarquardt(x, y, p0, lo, hi, w = null, maxIter = 2000, opts = {}) {
  const fm = opts.f || model, jac = opts.jac || jacRow, freeze = !!opts.freeze;
  w = w || x.map(() => 1);
  let p = clip(p0, lo, hi), f = rss(x, y, p, w, fm), lambda = 1e-3;
  const n = p.length;
  for (let it = 0; it < maxIter; it++) {
    const JtJ = Array.from({ length: n }, () => new Array(n).fill(0)), Jtr = new Array(n).fill(0);
    for (let i = 0; i < x.length; i++) {
      const J = jac(x[i], p), r = y[i] - fm(x[i], p);
      for (let a = 0; a < n; a++) {
        Jtr[a] += w[i] * J[a] * r;
        for (let b = a; b < n; b++) JtJ[a][b] += w[i] * J[a] * J[b];
      }
    }
    for (let a = 0; a < n; a++) for (let b = 0; b < a; b++) JtJ[a][b] = JtJ[b][a];
    // free parameters: all, or (freeze) those not pinned on a bound by a gradient pointing outwards
    const free = [];
    for (let a = 0; a < n; a++) {
      const pinned = freeze && ((p[a] <= lo[a] && Jtr[a] < 0) || (p[a] >= hi[a] && Jtr[a] > 0));
      if (!pinned) free.push(a);
    }
    if (!free.length) break;
    let improved = false;
    for (let tries = 0; tries < 30; tries++) {
      const A = free.map(a => free.map(b => (a === b ? JtJ[a][b] + lambda * (JtJ[a][b] + 1e-12) : JtJ[a][b])));
      const sub = solve(A, free.map(a => Jtr[a]));
      if (!sub) { lambda *= 10; continue; }
      const step = new Array(n).fill(0);
      free.forEach((a, i) => { step[a] = sub[i]; });
      const pn = clip(p.map((v, a) => v + step[a]), lo, hi), fn = rss(x, y, pn, w, fm);
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

// ---------------------------------------------------------------------------------------------- kplexity model
export const KPLEX_PARAMS = ["A", "S", "lam", "k0", "b"];          // reported form
const LN4 = Math.log(4);
const KLO = [0.05, 3, 0, 0, 1e-4], KHI = [5, 60, 1, 1, 2];        // internal (b, k0, A, R, lam), as kplex/model.py
// numpy default_rng(0).uniform(0.02, 0.5, 5), drawn six times: the random starts of kplex/model.py fit_kplexity
const RANDOM_STARTS = [
  [0.32574160991429807, 0.14949762260665775, 0.03966729148937345, 0.02793326505369396, 0.4103697148161308],
  [0.4581226770933065, 0.31118517236824633, 0.3701583492723192, 0.28093999590340296, 0.46883476341812874],
  [0.41160970597833546, 0.02131448008167109, 0.4315540527620333, 0.03612107614662289, 0.3702346142863732],
  [0.10431469788922831, 0.43432588272794553, 0.279901385719564, 0.16386170745794468, 0.22288986617487605],
  [0.03359344214982222, 0.07965597271979069, 0.3418997190529426, 0.33065096555564005, 0.3153848535110019],
  [0.20416522604570403, 0.4986607691788213, 0.4908009626125904, 0.34906015255073347, 0.3322204526085518],
];

const riseL = (k, b, k0) => 2 * Math.exp(Math.min(50, Math.max(-50, -b * (k - k0))));
const riseB = L => (L < 1e-8 ? 1 - L / 2 : -Math.expm1(-L) / L);
// dB/dL: series for small L, where the closed form loses precision
const riseDB = L => (L < 1e-4 ? -0.5 + L / 3 - L * L / 8 : (L * Math.exp(-L) + Math.expm1(-L)) / (L * L));

// chance-collision term B(k); with b = ln 4 and k0 = log4(G) the exact curve of random sequence
export function firstRise(k, b, k0) { return riseB(riseL(k, b, k0)); }
// expected y(k) of random sequence with G k-mer positions (no free parameters)
export function randomCurve(k, G) { const L = 2 * G / 4 ** k; return riseB(L); }

// reported form: p = {A, S, lam, k0, b}  (S = near-copy share at k0)
export function kplexityCurve(k, p) {
  return firstRise(k, p.b, p.k0) * (p.A - p.S * Math.exp(-p.lam * (k - p.k0)));
}
// internal form: q = [b, k0, A, R, lam],  y = B A (1 - R e^(-lam k))
export function kpInternal(k, q) {
  return firstRise(k, q[0], q[1]) * q[2] * (1 - q[3] * Math.exp(-q[4] * k));
}
export function kpJac(k, q) {
  const [b, k0, A, R, lam] = q, e = -b * (k - k0), clipped = e > 50 || e < -50;
  const L = riseL(k, b, k0), B = riseB(L), dB = riseDB(L), E = Math.exp(-lam * k), D = 1 - R * E;
  const dBdb = clipped ? 0 : dB * -(k - k0) * L, dBdk0 = clipped ? 0 : dB * b * L;
  return [dBdb * A * D, dBdk0 * A * D, B * D, -B * A * E, B * A * R * k * E];
}

// Multi-start bounded fit of the kplexity model, as kplex/model.py fit_kplexity.
// G (k-mer positions, about the genome size) seeds k0 at log4(G); without it k0 starts where y passes half its
// maximum. weights: undefined = spacing weights for uneven k (1 otherwise), null = unweighted, or an array.
export function fitKplexity(x, y, G = null, weights = undefined) {
  const w = weights === undefined ? (isUneven(x) ? spacingWeights(x) : null) : weights;
  const ymax = Math.max(...y);
  const k0g = G ? Math.log(G) / LN4 : x[Math.max(0, y.findIndex(v => v > 0.5 * ymax))];
  const A0 = Math.min(y[y.length - 1] + 0.02, 1);
  const eps = 1e-9, lo = KLO, hi = KHI;
  let best = null;
  for (let i = 0; i < 14; i++) {
    let q0;
    if (i < 8) q0 = [[1.4, 0.8, 2.0, 0.5][i % 4], k0g, A0, [0.3, 0.1, 0.5, 0.2][(i >> 1) % 4],
                     [0.02, 0.05, 0.1, 0.01, 0.2, 0.005][i % 6]];
    else q0 = RANDOM_STARTS[i - 8].map((u, j) => lo[j] + (hi[j] - lo[j]) * u);
    q0 = q0.map((v, j) => Math.min(hi[j] - eps, Math.max(lo[j] + eps, v)));
    const opts = { f: kpInternal, jac: kpJac, freeze: true };
    let res = levenbergMarquardt(x, y, q0, lo, hi, w, 3000, opts);
    for (let polish = 0; polish < 5; polish++) {       // restart (damping reset) until the curved b-lam valley is done
      const again = levenbergMarquardt(x, y, res.p, lo, hi, w, 3000, opts);
      const gain = (res.rss - again.rss) / Math.max(res.rss, 1e-300);
      if (again.rss < res.rss) res = again;
      if (!(gain > 1e-9)) break;
    }
    if (Number.isFinite(res.rss) && (!best || res.rss < best.rss)) best = res;
  }
  if (!best) return null;
  const q = best.p, [b, k0, A, R, lam] = q;
  const S = A * R * Math.exp(-lam * k0);
  const kmax = Math.max(...x), gap = S * Math.exp(-lam * (kmax - k0));
  const out = { A, identical: 1 - A, S, lam, divergence: 1 - Math.exp(-lam), k0, b, R, q,
                gap_at_kmax: gap, plateau_reached: gap < 0.02, kmax, G: G || null, weighted: !!w };
  Object.assign(out, stats(x, y, kpInternal, q, 5));
  return out;
}
