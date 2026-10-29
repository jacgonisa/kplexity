// "Reading the parameters": evidence plots (reference genomes + simulations) and the statistics quoted in the text,
// all for the kplexity model. Statistics are computed here from data/dtol_reference.json.
const pearson = (x, y) => {
  const n = x.length, mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxy / Math.sqrt(sxx * syy);
};
const median = v => { const s = [...v].filter(Number.isFinite).sort((a, b) => a - b), m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
const LN4 = Math.log(4);
const log4G = s => Math.log(s.Gbp) / LN4;
const pairs = (sp, fx, fy) => sp.map(s => [fx(s), fy(s)]).filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b));
const r = (sp, fx, fy) => { const p = pairs(sp, fx, fy); return pearson(p.map(q => q[0]), p.map(q => q[1])).toFixed(2); };
// least-squares slope through the origin
const slope0 = (sp, fx, fy) => { const p = pairs(sp, fx, fy);
  return p.reduce((a, [x, y]) => a + x * y, 0) / p.reduce((a, [x]) => a + x * x, 0); };
const te = s => (s.te_pct == null ? NaN : s.te_pct / 100);
const sat = s => (s.sat_pct == null ? NaN : s.sat_pct);
const sims = (ev, series) => ev.sims.filter(s => s.series === series);
const tandem300 = ev => sims(ev, "tandem").filter(s => /_300bp$/.test(s.name)).sort((a, b) => a.true_identical - b.true_identical);
const aging = ev => sims(ev, "aging").sort((a, b) => a.true_div - b.true_div);
const unit = s => (s.name.match(/_(\d+k?b)p?$/) || [, ""])[1];

export function fillStats(REF) {
  const sp = REF.species, ev = REF.evidence, ag = aging(ev), td = tandem300(ev);
  const agAt = d => ag.find(a => Math.abs(a.true_div - d) < 1e-6);
  const st = {
    r_k0: r(sp, log4G, s => s.k0), off_k0: median(sp.map(s => s.k0 - log4G(s))).toFixed(2),
    rand_k0: ev.random_k0_kplexity.toFixed(2), sim_log4G: (Math.log(ev.sim_G) / LN4).toFixed(2),
    rand_b: ev.random_b.toFixed(3), med_b: median(sp.map(s => s.b)).toFixed(2),
    r_G_b: r(sp, log4G, s => s.b), r_te_b: r(sp, te, s => s.b),
    sim_b: (() => { const b = ev.sims.map(s => s.b); return `${Math.min(...b).toFixed(2)}–${Math.max(...b).toFixed(2)}`; })(),
    n_te: String(pairs(sp, te, s => s.S).length), r_te_S: r(sp, te, s => s.S), slope_te_S: slope0(sp, te, s => s.S).toFixed(2),
    r_sat_S: r(sp, sat, s => s.S - slope0(sp, te, q => q.S) * te(s)),
    ag_S: [0.01, 0.1, 0.2].map(d => `${agAt(d).S.toFixed(2)} at ${Math.round(d * 100)}%`).join(", "),
    r_sat_id: r(sp, sat, s => 1 - s.A), r_te_id: r(sp, te, s => 1 - s.A), med_id: median(sp.map(s => 1 - s.A)).toFixed(3),
    td_err: Math.max(...td.map(s => Math.abs(1 - s.A - s.true_identical))).toFixed(3),
    td_range: `${Math.round(100 * td[0].true_identical)}–${Math.round(100 * td[td.length - 1].true_identical)}%`,
    ag_div: [0.01, 0.05, 0.1].map(d => `${(100 * agAt(d).divergence).toFixed(1)}% at ${Math.round(d * 100)}%`).join(", "),
    ag_div_hi: [0.15, 0.2].map(d => `${(100 * agAt(d).divergence).toFixed(0)}% at ${Math.round(d * 100)}%`).join(", "),
    med_div: (100 * median(sp.map(s => s.divergence))).toFixed(1),
    p05_div: (100 * [...sp.map(s => s.divergence)].sort((a, b) => a - b)[Math.floor(0.05 * sp.length)]).toFixed(1),
    p95_div: (100 * [...sp.map(s => s.divergence)].sort((a, b) => a - b)[Math.floor(0.95 * sp.length)]).toFixed(1),
  };
  document.querySelectorAll("[data-stat]").forEach(el => { el.textContent = st[el.dataset.stat] ?? "?"; });
}

export function renderEvidence(det, ctx) {
  const target = det.querySelector(":scope > .ev-plot");
  if (!target) return;
  const { REF, USER, COLORS, dark, bgColor } = ctx, sp = REF.species, ev = REF.evidence;
  const ink = dark() ? "#e6e6e1" : "#1f2328", grid = dark() ? "#2e343b" : "#ececE6", me = dark() ? "#ffffff" : "#000000";
  const clades = Object.keys(COLORS);
  const scatter = (fx, fy, xa = "x", ya = "y") => clades.map(c => {
    const s = sp.filter(x => x.clade === c && Number.isFinite(fx(x)) && Number.isFinite(fy(x)));
    return { x: s.map(fx), y: s.map(fy), xaxis: xa, yaxis: ya, mode: "markers", name: c, legendgroup: c,
             showlegend: true, marker: { size: 5, color: COLORS[c], opacity: 0.7 }, text: s.map(x => x.name),
             hovertemplate: "%{text}<br>%{x:.3f}, %{y:.3f}<extra></extra>" };
  });
  const byClade = (fy, xa, ya) => clades.map((c, i) => {
    const s = sp.filter(x => x.clade === c);
    return { x: s.map(() => i + (Math.random() - 0.5) * 0.5), y: s.map(fy), xaxis: xa, yaxis: ya, mode: "markers",
             name: c, legendgroup: c, showlegend: false, marker: { size: 5, color: COLORS[c], opacity: 0.7 },
             text: s.map(x => x.name), hovertemplate: "%{text}<br>%{y:.3f}<extra></extra>" };
  });
  // simulations: filled markers joined by a line (a series), open markers (other sims)
  const sim = (pts, fx, fy, xa, ya, name, line = true, symbol = "diamond") => ({
    x: pts.map(fx), y: pts.map(fy), xaxis: xa, yaxis: ya, mode: line ? "lines+markers" : "markers", name,
    line: { color: ink, width: 1.5 }, showlegend: true, text: pts.map(s => s.name),
    marker: { size: 8, color: line ? ink : bgColor(), line: { color: ink, width: 1.5 }, symbol: line ? "circle" : symbol },
    hovertemplate: "%{text}<br>%{x:.3f}, %{y:.3f}<extra>simulation</extra>" });
  const diag = (a, b, xa = "x", ya = "y", name = "equality", slope = 1) => ({ x: [a, b], y: [a * slope, b * slope],
    xaxis: xa, yaxis: ya, mode: "lines", line: { color: ink, dash: "dash", width: 1 }, name, hoverinfo: "skip" });
  const hline = (v, xa, ya, label, dash = "solid", col = me) => (v == null || !Number.isFinite(v) ? [] : [{
    type: "line", xref: `${xa} domain`, x0: 0, x1: 1, yref: ya, y0: v, y1: v, line: { color: col, width: dash === "solid" ? 2 : 1, dash },
    label: { text: label, textposition: "end", font: { size: 10, color: ink } } }]);
  const ax = (title, extra = {}) => ({ title: { text: title, font: { size: 11 } }, gridcolor: grid, zeroline: false, ...extra });
  const two = { grid: { rows: 1, columns: 2, pattern: "independent" } };
  const catAx = title => ax(title, { tickvals: clades.map((_, i) => i), ticktext: clades });
  const u = USER?.kp, ag = aging(ev), td = tandem300(ev);
  let tr = [], lay = {};
  switch (det.dataset.plot) {
    case "k0": {
      tr = scatter(log4G, s => s.k0);
      tr.push(diag(8, 17, "x", "y", "k0 = log₄G"));
      if (u && USER.G) tr.push({ x: [Math.log(USER.G) / LN4], y: [u.k0], mode: "markers", name: "your genome",
                                 marker: { size: 13, color: me, symbol: "star" } });
      lay = { xaxis: ax("log₄(genome size)"), yaxis: ax("k0") };
      break;
    }
    case "b":
      tr = [...scatter(log4G, s => s.b),
            sim(ev.sims, s => s.true_repeat + s.true_identical, s => s.b, "x2", "y2", "simulated genomes", false)];
      lay = { ...two, xaxis: ax("log₄(genome size)"), yaxis: ax("b (reference genomes)", { range: [0.3, 1.6] }),
              xaxis2: ax("repeat share of the simulated genome"), yaxis2: ax("b (simulation)", { range: [0.3, 1.6] }),
              shapes: [...hline(LN4, "x", "y", "ln 4", "dash", ink), ...hline(LN4, "x2", "y2", "ln 4", "dash", ink),
                       ...hline(u?.b, "x", "y", "yours")] };
      break;
    case "S": {
      const sl = slope0(sp, te, s => s.S), tp = sims(ev, "transposon"), tp3 = tp.filter(s => /_3kb$/.test(s.name))
        .sort((a, b) => a.true_repeat - b.true_repeat);
      tr = [...scatter(s => 100 * te(s), s => s.S),
            { x: [0, 90], y: [0, 0.9 * sl], mode: "lines", line: { color: ink, dash: "dash", width: 1 },
              name: `S = ${sl.toFixed(2)} × TE share`, hoverinfo: "skip" },
            sim(tp3, s => s.true_repeat, s => s.S, "x2", "y2", "transposons, 15% divergence, 3 kb"),
            sim(tp.filter(s => !/_3kb$/.test(s.name)), s => s.true_repeat, s => s.S, "x2", "y2",
                "transposons, 30%, 1–20 kb", false),
            sim(ag, s => s.true_repeat, s => s.S, "x2", "y2", "30% in 5-kb copies, divergence 0–40%", false, "square"),
            diag(0, 0.5, "x2", "y2")];
      lay = { ...two, xaxis: ax("transposable elements (% of genome, EDTA)"), yaxis: ax("S"),
              xaxis2: ax("repeat share of the simulated genome"), yaxis2: ax("S (simulation)", { rangemode: "tozero" }),
              shapes: hline(u?.S, "x", "y", "yours") };
      break;
    }
    case "identical": {
      const other = sims(ev, "tandem").filter(s => !/_300bp$/.test(s.name));
      tr = [sim(td, s => s.true_identical, s => 1 - s.A, "x", "y", "tandem arrays, 300-bp units"),
            { ...sim(other, s => s.true_identical, s => 1 - s.A, "x", "y", "30% tandem, other unit lengths", false),
              text: other.map(s => `${s.name} (${unit(s)} units)`) },
            diag(0, 0.52),
            ...scatter(sat, s => 1 - s.A, "x2", "y2")];
      lay = { ...two, xaxis: ax("perfect tandem arrays (share of simulated genome)"), yaxis: ax("1 − A (simulation)"),
              xaxis2: ax("satellites (% of genome, TRASH)"), yaxis2: ax("1 − A (reference genomes)"),
              shapes: hline(u ? 1 - u.A : null, "x2", "y2", "yours") };
      break;
    }
    case "divergence": {
      const tp = sims(ev, "transposon");
      tr = [sim(ag, s => s.true_div, s => s.divergence, "x", "y", "30% in 5-kb copies, divergence 0–40%"),
            sim(tp, s => s.true_div, s => s.divergence, "x", "y", "transposons, 15% divergence", false),
            diag(0, 0.42),
            ...byClade(s => 100 * s.divergence, "x2", "y2")];
      lay = { ...two, xaxis: ax("divergence between copies (simulated)"), yaxis: ax("fitted divergence"),
              xaxis2: catAx(""), yaxis2: ax("divergence, % (reference genomes)", { rangemode: "tozero" }),
              shapes: hline(u ? 100 * u.divergence : null, "x2", "y2", "yours") };
      break;
    }
    default:
      return;
  }
  Plotly.react(target, tr, {
    font: { family: "Arial, Helvetica, sans-serif", color: ink, size: 11 },
    paper_bgcolor: bgColor(), plot_bgcolor: bgColor(), margin: { l: 56, r: 12, t: 10, b: 44 },
    legend: { orientation: "h", y: -0.22, font: { size: 10 } }, hovermode: "closest", ...lay,
  }, { displayModeBar: false, responsive: true });
}
