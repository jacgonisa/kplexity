// "Reading the parameters": evidence plots (reference genomes + simulations) and the statistics quoted in the text.
const pearson = (x, y) => {
  const n = x.length, mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxy / Math.sqrt(sxx * syy);
};
const median = v => { const s = [...v].sort((a, b) => a - b), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
const log4G = s => Math.log(s.G * 1e6) / Math.log(4);
const pairs = (sp, fx, fy) => sp.map(s => [fx(s), fy(s)]).filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b));
const r = (sp, fx, fy) => { const p = pairs(sp, fx, fy); return pearson(p.map(q => q[0]), p.map(q => q[1])).toFixed(2); };

export function fillStats(REF) {
  const sp = REF.species, ev = REF.evidence, A = ev.aging;
  const at = d => A.find(a => a.div === d);
  const st = {
    r_k01: r(sp, log4G, s => s.k0_1), off_k01: median(sp.map(s => s.k0_1 - log4G(s))).toFixed(1),
    rand_k0: ev.random_k0.toFixed(1), rand_s: ev.random_s.toFixed(1),
    r_te_L2: r(sp, s => s.te_pct, s => s.L2), r_sat_L2: r(sp, s => s.sat_pct, s => s.L2),
    r_te_k02: r(sp, s => s.te_pct, s => s.k0_2), r_G_k02: r(sp, log4G, s => s.k0_2),
    r_sat_na: r(sp, s => s.sat_pct, s => 1 - s.asym), r_te_na: r(sp, s => s.te_pct, s => 1 - s.asym),
    med_s1: median(sp.map(s => s.s1)).toFixed(2), r_te_s1: r(sp, s => s.te_pct, s => s.s1),
    r_te_s2: r(sp, s => s.te_pct, s => s.s2),
    aging_k02: [1, 5, 15].map(d => `${Math.round(at(d).k0_2)} at ${d}% divergence`).join(", "),
    aging_s2: [1, 15].map(d => `${at(d).s2.toFixed(2)} at ${d}%`).join(" to "),
  };
  document.querySelectorAll("[data-stat]").forEach(el => { el.textContent = st[el.dataset.stat] ?? "?"; });
}

export function renderEvidence(det, ctx) {
  const { REF, USER, COLORS, dark, bgColor } = ctx, sp = REF.species, ev = REF.evidence;
  const ink = dark() ? "#e6e6e1" : "#1f2328", grid = dark() ? "#2e343b" : "#ececE6", me = dark() ? "#ffffff" : "#000000";
  const clades = Object.keys(COLORS);
  const scatter = (fx, fy, xa = "x", ya = "y") => clades.map(c => {
    const s = sp.filter(x => x.clade === c && Number.isFinite(fx(x)) && Number.isFinite(fy(x)));
    return { x: s.map(fx), y: s.map(fy), xaxis: xa, yaxis: ya, mode: "markers", name: c, legendgroup: c,
             showlegend: xa === "x", marker: { size: 5, color: COLORS[c], opacity: 0.7 }, text: s.map(x => x.name),
             hovertemplate: "%{text}<br>%{x:.2f}, %{y:.3f}<extra></extra>" };
  });
  // open markers: no separate second transition (k0_2 within 2 of k0_1), so the L1/L2 split is not meaningful
  const sim = (xs, ys, xa, ya, name, merged = null) => ({ x: xs, y: ys, xaxis: xa, yaxis: ya, mode: "lines+markers", name,
    line: { color: ink, width: 1.5 }, showlegend: false,
    marker: { size: 8, color: merged ? merged.map(m => (m ? bgColor() : ink)) : ink, line: { color: ink, width: 1.5 } },
    hovertemplate: "%{x}: %{y:.3f}<extra>simulation</extra>" });
  const hline = (v, xa, ya, label) => (v == null ? [] : [{ type: "line", xref: `${xa} domain`, x0: 0, x1: 1, yref: ya,
    y0: v, y1: v, line: { color: me, width: 2 }, label: { text: label, textposition: "end", font: { size: 10, color: ink } } }]);
  const ax = (title, extra = {}) => ({ title: { text: title, font: { size: 11 } }, gridcolor: grid, zeroline: false, ...extra });
  const two = { grid: { rows: 1, columns: 2, pattern: "independent" } };
  const A = ev.aging, T = ev.tandem, u = USER?.fit, merged = A.map(a => a.k0_2 - a.k0_1 < 2);
  let tr = [], lay = {};
  switch (det.dataset.plot) {
    case "k0_1": {
      tr = scatter(log4G, s => s.k0_1);
      tr.push({ x: [8, 17], y: [8, 17], mode: "lines", line: { color: ink, dash: "dash", width: 1 }, name: "k0_1 = log₄G",
                hoverinfo: "skip" });
      if (u && USER.G) tr.push({ x: [Math.log(USER.G) / Math.log(4)], y: [u.k0_1], mode: "markers", name: "your genome",
                                 marker: { size: 13, color: me, symbol: "star" } });
      lay = { xaxis: ax("log₄(genome size)"), yaxis: ax("k0_1") };
      break;
    }
    case "L2":
      tr = [...scatter(s => s.te_pct, s => s.L2),
            sim(A.map(a => a.div), A.map(a => a.L2), "x2", "y2", "L2", merged)];
      lay = { ...two, xaxis: ax("transposable elements (% of genome)"), yaxis: ax("L2"),
              xaxis2: ax("divergence between copies (%)"), yaxis2: ax("L2 (simulation, 30% repeats)", { rangemode: "tozero" }),
              shapes: hline(u?.L2, "x", "y", "yours") };
      break;
    case "k0_2":
      tr = [sim(A.map(a => a.div), A.map(a => a.k0_2), "x", "y", "k0_2", merged),
            ...clades.map((c, i) => {
              const s = sp.filter(x => x.clade === c);
              return { x: s.map(() => i + (Math.random() - 0.5) * 0.5), y: s.map(x => x.k0_2), xaxis: "x2", yaxis: "y2",
                       mode: "markers", name: c, marker: { size: 5, color: COLORS[c], opacity: 0.7 }, text: s.map(x => x.name),
                       hovertemplate: "%{text}<br>k0_2 = %{y:.1f}<extra></extra>" };
            })];
      lay = { ...two, xaxis: ax("divergence between copies (%)"), yaxis: ax("k0_2 (simulation)"),
              xaxis2: ax("", { tickvals: clades.map((_, i) => i), ticktext: clades }), yaxis2: ax("k0_2 (reference genomes)"),
              shapes: hline(u?.k0_2, "x2", "y2", "yours") };
      break;
    case "asym":
      tr = [sim(T.map(t => t.pct), T.map(t => 1 - t.asym), "x", "y"),
            { x: [0, 52], y: [0, 0.52], mode: "lines", line: { color: ink, dash: "dash", width: 1 }, hoverinfo: "skip",
              showlegend: false },
            ...scatter(s => s.sat_pct, s => 1 - s.asym, "x2", "y2").map(t => ({ ...t, showlegend: true }))];
      lay = { ...two, xaxis: ax("perfect tandem arrays (% of simulated genome)"), yaxis: ax("1 − asymptote (simulation)"),
              xaxis2: ax("satellites (% of genome)"), yaxis2: ax("1 − asymptote"),
              shapes: hline(u ? 1 - u.asymptote : null, "x2", "y2", "yours") };
      break;
    case "s":
      tr = [...scatter(s => s.te_pct, s => s.s1),
            sim(A.map(a => a.div), A.map(a => a.s2), "x2", "y2", "s2", merged)];
      lay = { ...two, xaxis: ax("transposable elements (% of genome)"), yaxis: ax("s1", { range: [0, Math.max(1.9, ev.random_s + 0.3)] }),
              xaxis2: ax("divergence between copies (%)"), yaxis2: ax("s2 (simulation)"),
              shapes: [{ type: "line", xref: "x domain", x0: 0, x1: 1, y0: ev.random_s, y1: ev.random_s,
                         line: { color: ink, dash: "dash", width: 1 },
                         label: { text: "random sequence", textposition: "end", font: { size: 10, color: ink } } },
                       ...hline(u?.s1, "x", "y", "yours")] };
      break;
  }
  Plotly.react(det.querySelector(".ev-plot"), tr, {
    font: { family: "Arial, Helvetica, sans-serif", color: ink, size: 11 },
    paper_bgcolor: bgColor(), plot_bgcolor: bgColor(), margin: { l: 56, r: 12, t: 10, b: 44 },
    legend: { orientation: "h", y: -0.22, font: { size: 10 } }, hovermode: "closest", ...lay,
  }, { displayModeBar: false, responsive: true });
}
