import { fitDoubleSigmoid, fitSingleSigmoid, doubleSigmoid, PARAMS } from "./fit.js";

const COLORS = { Viridiplantae: "#228833", Vertebrates: "#4477AA", Invertebrates: "#CC6677", Fungi: "#AA4499" };
const SHORT = { Viridiplantae: "plants", Vertebrates: "vertebr.", Invertebrates: "invertebr.", Fungi: "fungi" };
const SPECIES_COLORS = ["#EE7733", "#009988", "#33BBEE", "#EE3377", "#BBBBBB"];
const $ = id => document.getElementById(id);
const dark = () => matchMedia("(prefers-color-scheme: dark)").matches;

let REF = null;              // reference data
let USER = null;             // {name, k, y, fit, single}
const shownClades = new Set(["Viridiplantae", "Vertebrates", "Invertebrates"]);
const shownSpecies = [];

// ---------------------------------------------------------------- parsing
function parseCurve(text) {
  const lines = text.trim().split(/\r?\n/).filter(l => l.trim() && !l.startsWith("#"));
  if (lines.length < 5) throw new Error("need at least 5 rows");
  const split = l => l.trim().split(/[,\t; ]+/);
  let header = split(lines[0]).map(h => h.toLowerCase());
  const hasHeader = header.some(h => isNaN(parseFloat(h)));
  const rows = (hasHeader ? lines.slice(1) : lines).map(split);
  if (!hasHeader) header = rows[0].length === 2 ? ["k", "fraction"] : header.map((_, i) => `c${i}`);
  const col = name => header.indexOf(name);
  const iK = col("k") >= 0 ? col("k") : 0;
  let get;
  if (col("unique_kmers") >= 0 && col("total_kmers") >= 0) {
    get = r => parseFloat(r[col("unique_kmers")]) / parseFloat(r[col("total_kmers")]);
  } else if (col("fraction_unique_theoretical") >= 0) {
    get = r => parseFloat(r[col("fraction_unique_theoretical")]);
  } else if (col("fraction_unique") >= 0) {
    get = r => parseFloat(r[col("fraction_unique")]);
  } else {
    get = r => parseFloat(r[iK === 0 ? 1 : 0]);
  }
  const pts = rows.map(r => [parseFloat(r[iK]), get(r)])
                  .filter(([k, y]) => Number.isFinite(k) && Number.isFinite(y))
                  .sort((a, b) => a[0] - b[0]);
  if (pts.length < 7) throw new Error("fewer than 7 numeric rows (need k and a fraction or U and T)");
  if (pts.some(([, y]) => y < 0 || y > 1.0001)) throw new Error("fractions must lie between 0 and 1");
  return { k: pts.map(p => p[0]), y: pts.map(p => p[1]) };
}

// ---------------------------------------------------------------- user curve
function setStatus(msg, err = false) {
  const s = $("status"); s.textContent = msg; s.classList.toggle("err", err);
}

function loadText(text, name) {
  try {
    const { k, y } = parseCurve(text);
    setStatus(`Fitting ${k.length} points…`);
    setTimeout(() => {
      const fit = fitDoubleSigmoid(k, y), single = fitSingleSigmoid(k, y);
      if (!fit) { setStatus("The double-sigmoid fit failed on this curve.", true); return; }
      USER = { name, k, y, fit, single };
      const warn = (k[0] > 5 || k[k.length - 1] < 151) ? ` (k = ${k[0]}–${k[k.length - 1]}; the reference uses 5–151)` : "";
      const wnote = fit.weighted ? "; uneven k, so each point is weighted by the k-span it covers" : "";
      setStatus(`${name}: ${k.length} points, R² = ${fit.R2.toFixed(4)}${warn}${wnote}`);
      $("btn-csv").disabled = false; $("btn-png").disabled = false;
      render();
    }, 20);
  } catch (e) {
    setStatus(`Could not read the curve: ${e.message}`, true);
  }
}

// ---------------------------------------------------------------- plot
const hexA = (hex, a) => {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
};

function band(k, lo, hi, color, a, name) {
  return [
    { x: k, y: lo, mode: "lines", line: { width: 0 }, hoverinfo: "skip", showlegend: false },
    { x: k, y: hi, mode: "lines", line: { width: 0 }, fill: "tonexty", fillcolor: hexA(color, a),
      hoverinfo: "skip", showlegend: false, name },
  ];
}

function plotTraces() {
  const tr = [], K = REF.k, mode = $("bands").value;
  for (const c of Object.keys(COLORS)) {
    if (!shownClades.has(c)) continue;
    const q = REF.clades[c], col = COLORS[c];
    if (mode === "both" || mode === "90") tr.push(...band(K, q.p05, q.p95, col, 0.10, c));
    if (mode === "both" || mode === "50") tr.push(...band(K, q.p25, q.p75, col, 0.22, c));
    tr.push({ x: K, y: q.p50, mode: "lines", line: { color: col, width: 2.2 },
              name: `${c} median (n = ${q.n})`, hovertemplate: `${c} median<br>k=%{x}: %{y:.3f}<extra></extra>` });
  }
  shownSpecies.forEach((s, i) => tr.push({
    x: K, y: s.y, mode: "lines", line: { color: SPECIES_COLORS[i % SPECIES_COLORS.length], width: 1.6, dash: "dot" },
    name: s.name, hovertemplate: `${s.name}<br>k=%{x}: %{y:.3f}<extra></extra>` }));
  if (USER) {
    const f = USER.fit, kk = Array.from({ length: 400 }, (_, i) => 5 + i * 146 / 399);
    tr.push({ x: USER.k, y: USER.y, mode: "markers", marker: { size: 5, color: dark() ? "#e6e6e1" : "#1f2328" },
              name: `${USER.name} (data)`, hovertemplate: "k=%{x}: %{y:.4f}<extra></extra>" });
    tr.push({ x: kk, y: kk.map(x => doubleSigmoid(x, f.p)), mode: "lines",
              line: { color: dark() ? "#ffffff" : "#000000", width: 2.6 }, name: "double-sigmoid fit" });
    if ($("show-components").checked) {
      tr.push({ x: kk, y: kk.map(x => f.L1 / (1 + Math.exp(-f.s1 * (x - f.k0_1)))), mode: "lines",
                line: { color: "#888", width: 1.4, dash: "dash" }, name: "component 1 (L1)" });
      tr.push({ x: kk, y: kk.map(x => f.L2 / (1 + Math.exp(-f.s2 * (x - f.k0_2)))), mode: "lines",
                line: { color: "#d1495b", width: 1.6, dash: "dash" }, name: "component 2 (L2)" });
    }
  }
  return tr;
}

function layout() {
  const ink = dark() ? "#e6e6e1" : "#1f2328", grid = dark() ? "#2e343b" : "#ececE6";
  const shapes = [], ann = [];
  if (USER) {
    for (const [key, lab] of [["k0_1", "k0_1"], ["k0_2", "k0_2"]]) {
      shapes.push({ type: "line", x0: USER.fit[key], x1: USER.fit[key], y0: 0, y1: 1, yref: "paper",
                    line: { color: "#999", width: 1, dash: "dot" } });
      ann.push({ x: USER.fit[key], y: 1, yref: "paper", text: `${lab} = ${USER.fit[key].toFixed(1)}`,
                 showarrow: false, yanchor: "bottom", font: { size: 11, color: "#888" } });
    }
  }
  return {
    font: { family: "Arial, Helvetica, sans-serif", color: ink, size: 13 },
    paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
    margin: { l: 58, r: 12, t: 30, b: 50 },
    xaxis: { title: "k", range: [0, 155], gridcolor: grid, zeroline: false },
    yaxis: { title: "fraction of unique k-mers", range: [0, 1.02], gridcolor: grid, zeroline: false },
    legend: { orientation: "h", y: -0.18, font: { size: 11 } },
    shapes, annotations: ann, hovermode: "closest",
  };
}

// ---------------------------------------------------------------- results
const fmt = (v, d = 3) => (Math.abs(v) >= 100 ? v.toFixed(1) : v.toFixed(d));
const ROWS = [
  ["L1", "L1", 3], ["s1", "s1", 3], ["k0_1", "k0_1", 2], ["L2", "L2", 3], ["s2", "s2", 3], ["k0_2", "k0_2", 2],
  ["asymptote", "asym", 3], ["L2/(L1+L2)", "frac", 3],
];
const userVal = key => ({ asym: USER.fit.asymptote, frac: USER.fit.L2_frac }[key] ?? USER.fit[key]);

function percentile(vals, v) {
  let below = 0;
  for (const x of vals) if (x < v) below++;
  return Math.round(100 * below / vals.length);
}

function renderResults() {
  if (!USER) return;
  const clades = Object.keys(COLORS).filter(c => shownClades.has(c));
  const byClade = Object.fromEntries(clades.map(c => [c, REF.species.filter(s => s.clade === c)]));
  let h = `<p class="muted small">Columns after the value: percentile of your genome within each clade.</p>` +
          `<div class="table-wrap"><table><tr><th>parameter</th><th>value</th>` +
          clades.map(c => `<th title="percentile within ${c} (0 = lowest, 100 = highest)"><span class="swatch" style="display:inline-block;background:${COLORS[c]}"></span> ${SHORT[c]}</th>`).join("") + `</tr>`;
  for (const [lab, key, d] of ROWS) {
    const v = userVal(key);
    h += `<tr><td>${lab}</td><td><b>${fmt(v, d)}</b></td>` +
         clades.map(c => `<td>${percentile(byClade[c].map(s => s[key]), v)}</td>`).join("") + `</tr>`;
  }
  h += `</table></div>`;
  const dA = USER.single.AICc - USER.fit.AICc;
  h += `<p class="model">R² = ${USER.fit.R2.toFixed(4)} · RMSE = ${USER.fit.RMSE.toExponential(2)}<br>` +
       `Double vs single sigmoid: ΔAICc = ${dA.toFixed(1)} → ` +
       (dA > 10 ? "<b>double sigmoid</b> clearly preferred (a repeat-driven second transition)."
        : dA > 2 ? "double sigmoid preferred."
        : "<b>single sigmoid</b> is enough (no clear second transition, as in most bacteria).") + `</p>`;
  $("results").innerHTML = h;
}

function renderStrips() {
  const box = $("strips");
  box.innerHTML = "";
  const clades = Object.keys(COLORS).filter(c => shownClades.has(c));
  for (const [lab, key] of [["asymptote (L1+L2)", "asym"], ["L2", "L2"], ["k0_2", "k0_2"], ["k0_1", "k0_1"]]) {
    const div = document.createElement("div"); div.className = "strip"; box.appendChild(div);
    const tr = clades.map((c, i) => {
      const s = REF.species.filter(x => x.clade === c);
      return { x: s.map(x => x[key]), y: s.map(() => i + (Math.random() - 0.5) * 0.55), mode: "markers",
               marker: { size: 4, color: hexA(COLORS[c], 0.45) }, text: s.map(x => x.name),
               hovertemplate: `%{text}<br>${lab} = %{x:.3f}<extra></extra>`, showlegend: false };
    });
    const shapes = USER ? [{ type: "line", x0: userVal(key), x1: userVal(key), y0: 0, y1: 1, yref: "paper",
                             line: { color: dark() ? "#fff" : "#000", width: 2.5 } }] : [];
    const ink = dark() ? "#e6e6e1" : "#1f2328";
    Plotly.newPlot(div, tr, {
      font: { family: "Arial, Helvetica, sans-serif", color: ink, size: 11 },
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      margin: { l: 90, r: 10, t: 22, b: 26 }, title: { text: lab, font: { size: 12 }, x: 0.01, y: 0.97 },
      yaxis: { tickvals: clades.map((_, i) => i), ticktext: clades, range: [-0.6, clades.length - 0.4],
               gridcolor: "rgba(0,0,0,0)", zeroline: false },
      xaxis: { gridcolor: dark() ? "#2e343b" : "#ececE6", zeroline: false }, shapes,
    }, { displayModeBar: false, responsive: true });
  }
}

function render() {
  Plotly.react("plot", plotTraces(), layout(), { responsive: true, displaylogo: false });
  renderResults();
  renderStrips();
}

// ---------------------------------------------------------------- controls
function buildControls() {
  const box = $("clade-toggles");
  for (const c of Object.keys(COLORS)) {
    const lab = document.createElement("label"); lab.className = "toggle";
    lab.innerHTML = `<input type="checkbox" ${shownClades.has(c) ? "checked" : ""}>` +
                    `<span class="swatch" style="background:${COLORS[c]}"></span>${c} <span class="muted">(n = ${REF.clades[c].n})</span>`;
    lab.querySelector("input").addEventListener("change", e => {
      e.target.checked ? shownClades.add(c) : shownClades.delete(c); render();
    });
    box.appendChild(lab);
  }
  $("species-list").innerHTML = REF.species.map(s => `<option value="${s.name}">`).join("");
  $("ref-n").textContent = REF.species.length.toLocaleString();
  $("bands").addEventListener("change", render);
  $("show-components").addEventListener("change", render);
  const addSpecies = () => {
    const name = $("species").value.trim(), s = REF.species.find(x => x.name.toLowerCase() === name.toLowerCase());
    if (!s || shownSpecies.includes(s)) return;
    shownSpecies.push(s); $("species").value = ""; renderChips(); render();
  };
  $("btn-species").addEventListener("click", addSpecies);
  $("species").addEventListener("keydown", e => { if (e.key === "Enter") addSpecies(); });

  const drop = $("drop"), file = $("file");
  const readFile = f => f && f.text().then(t => loadText(t, f.name.replace(/\.(csv|tsv|txt)$/i, "")));
  file.addEventListener("change", () => readFile(file.files[0]));
  ["dragenter", "dragover"].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.add("over"); }));
  ["dragleave", "drop"].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.remove("over"); }));
  drop.addEventListener("drop", e => readFile(e.dataTransfer.files[0]));
  $("btn-paste").addEventListener("click", () => loadText($("paste").value, "pasted curve"));
  $("btn-example").addEventListener("click", () => loadText(REF.example.csv, "Arabidopsis thaliana (TAIR12)"));
  $("btn-csv").addEventListener("click", downloadCsv);
  $("btn-png").addEventListener("click", () => Plotly.downloadImage("plot", { format: "png", width: 1400, height: 800,
                                                                             filename: `${USER.name}_kplexity` }));
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
}

function renderChips() {
  const box = $("species-chips");
  box.innerHTML = "";
  shownSpecies.forEach((s, i) => {
    const chip = document.createElement("span"); chip.className = "chip";
    chip.innerHTML = `<span class="swatch" style="background:${SPECIES_COLORS[i % SPECIES_COLORS.length]}"></span>` +
                     `<i>${s.name}</i> <span class="muted">(${s.clade})</span><button type="button" aria-label="remove">×</button>`;
    chip.querySelector("button").addEventListener("click", () => { shownSpecies.splice(i, 1); renderChips(); render(); });
    box.appendChild(chip);
  });
}

function downloadCsv() {
  const f = USER.fit;
  const rows = [["parameter", "value"], ...PARAMS.map(p => [p, f[p]]), ["asymptote", f.asymptote],
                ["L2_over_L1_plus_L2", f.L2_frac], ["R2", f.R2], ["RMSE", f.RMSE], ["AICc_double", f.AICc],
                ["AICc_single", USER.single.AICc]];
  const blob = new Blob([rows.map(r => r.join(",")).join("\n") + "\n"], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = `${USER.name}_kplexity_fit.csv`; a.click();
  URL.revokeObjectURL(a.href);
}

// ---------------------------------------------------------------- start
fetch("data/dtol_reference.json")
  .then(r => r.json())
  .then(ref => {
    REF = ref; buildControls(); render();
    setStatus("Reference loaded. Load a curve or try the example.");
    const q = new URLSearchParams(location.search);
    if (q.has("example")) loadText(REF.example.csv, "Arabidopsis thaliana (TAIR12)");
    if (q.get("species")) q.get("species").split(",").forEach(n => {
      const s = REF.species.find(x => x.name.toLowerCase() === n.trim().toLowerCase());
      if (s && !shownSpecies.includes(s)) shownSpecies.push(s);
    });
    renderChips(); render();
  })
  .catch(e => setStatus(`Could not load the reference data: ${e.message}`, true));
