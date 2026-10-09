// Check the browser fitter against the Python fitter (kplex / paper analysis).
//   node docs/test/fit_check.mjs
// 1) 30 full DToL curves from data/dtol_reference.json (parameters fitted in Python, rounded to 4 d.p.)
// 2) 10 sketched curves (25 k, size-anchored; spacing-weighted fit) from test/sketch_fixture.json
import fs from "fs";
import { fitDoubleSigmoid, doubleSigmoid, PARAMS } from "../fit.js";

const here = new URL(".", import.meta.url).pathname;
const ref = JSON.parse(fs.readFileSync(here + "../data/dtol_reference.json"));
const sketch = JSON.parse(fs.readFileSync(here + "sketch_fixture.json"));
const rss = (k, y, p) => k.reduce((s, x, i) => s + (y[i] - doubleSigmoid(x, p)) ** 2, 0);

function check(name, k, y, pPy) {
  const f = fitDoubleSigmoid(k, y);
  const rel = f.p.map((v, i) => Math.abs(v - pPy[i]) / Math.max(Math.abs(pPy[i]), 1e-3));
  const maxrel = Math.max(...rel);
  const dRss = (rss(k, y, f.p) - rss(k, y, pPy)) / rss(k, y, pPy);   // unweighted RSS, both fits
  return { name, ok: maxrel <= 0.01 || dRss <= 1e-3, maxrel, dRss };
}

const step = Math.max(1, Math.floor(ref.species.length / 30));
const full = ref.species.filter((_, i) => i % step === 0).slice(0, 30)
  .map(s => {                                   // some curves stop before k = 151 (null beyond)
    const ok = s.y.map(v => v !== null);
    return check(s.name, ref.k.filter((_, i) => ok[i]), s.y.filter(v => v !== null), PARAMS.map(p => s[p]));
  });
const sk = sketch.map(s => check(s.name + " (25 k)", s.k, s.y, s.p));
let bad = 0;
for (const r of [...full, ...sk]) {
  if (!r.ok) bad++;
  console.log(r.ok ? "ok  " : "FAIL", r.name.padEnd(42), "max rel diff", r.maxrel.toFixed(4), " dRSS", r.dRss.toExponential(1));
}
console.log(`${full.length + sk.length - bad}/${full.length + sk.length} agree with Python`);
process.exit(bad ? 1 : 0);
