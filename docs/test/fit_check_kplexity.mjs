// Check the browser kplexity-model fitter against the Python one (kplex/model.py fit_kplexity).
//   node docs/test/fit_check_kplexity.mjs          (--all: every reference curve instead of 30)
// 1) 30 full DToL curves from data/dtol_reference.json (parameters fitted in Python, rounded to 4 d.p.; unweighted,
//    as the reference fits were made, also for the few curves with missing k after 142)
// 2) sketched curves (about 25 k, size-anchored; spacing-weighted fit) from test/sketch_fixture_kplexity.json
// Pass: max relative difference in A, S, lam, k0, b <= 1% (floor 1e-3), or the JS fit's RSS no more than
// 1e-3 (relative) above the Python fit's — the second covers near-flat directions (e.g. A against S when the curve
// has not levelled off by the last k), where different parameter values describe the curve equally well.
import fs from "fs";
import { fitKplexity, kplexityCurve, KPLEX_PARAMS } from "../fit.js";

const here = new URL(".", import.meta.url).pathname;
const ref = JSON.parse(fs.readFileSync(here + "../data/dtol_reference.json"));
const sketch = JSON.parse(fs.readFileSync(here + "sketch_fixture_kplexity.json"));
const rss = (k, y, p) => k.reduce((s, x, i) => s + (y[i] - kplexityCurve(x, p)) ** 2, 0);

function check(name, k, y, G, pPy, weights = undefined) {
  const f = fitKplexity(k, y, G, weights);
  const rel = KPLEX_PARAMS.map(n => Math.abs(f[n] - pPy[n]) / Math.max(Math.abs(pPy[n]), 1e-3));
  const maxrel = Math.max(...rel), worst = KPLEX_PARAMS[rel.indexOf(maxrel)];
  const dRss = (rss(k, y, f) - rss(k, y, pPy)) / rss(k, y, pPy);   // unweighted RSS, both fits
  return { name, ok: maxrel <= 0.01 || dRss <= 1e-3, maxrel, worst, dRss, plateau: f.plateau_reached };
}

const nFull = process.argv.includes("--all") ? ref.species.length : 30;
const step = Math.max(1, Math.floor(ref.species.length / nFull));
const full = ref.species.filter((_, i) => i % step === 0).slice(0, nFull)
  .map(s => {                                   // some curves stop before k = 151 (null beyond)
    const ok = s.y.map(v => v !== null);
    return check(s.name, ref.k.filter((_, i) => ok[i]), s.y.filter(v => v !== null), s.Gbp, s, null);
  });
const sk = sketch.map(s => check(`${s.name} (${s.k.length} k)`, s.k, s.y, s.G, s.p));
let bad = 0;
for (const r of [...full, ...sk]) {
  if (!r.ok) bad++;
  console.log(r.ok ? "ok  " : "FAIL", r.name.padEnd(42), "max rel diff", r.maxrel.toFixed(4), `(${r.worst})`.padEnd(6),
              " dRSS", r.dRss.toExponential(1), r.plateau ? "" : " plateau not reached");
}
console.log(`${full.length + sk.length - bad}/${full.length + sk.length} agree with Python ` +
            `(${full.length} full curves, ${sk.length} sketches)`);
process.exit(bad ? 1 : 0);
