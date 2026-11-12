/* Public DNA-PAINT measurements. Rendering does not fit models or select on test outcomes. */

const dataURL = file => new URL("data/" + file, import.meta.url).href;
const protocolURL = "https://github.com/kangk1204/blink-atlas/blob/main/analysis/palette-protocol.md";
const rawProtocolURL = "https://github.com/kangk1204/blink-atlas/blob/main/analysis/palette-raw.md";
const running = new WeakMap();
const specimens = {origami: "DNA origami", npc: "Nuclear pore complexes"};
const rawMetrics = {
  darkTimeS: {label: "Dark time", unit: "s"},
  eventPhotons: {label: "Photons per binding event", unit: "photons"},
};
const roleNames = {
  target_field: "Target · full field", target_sites: "Target · docking-site picks",
  target_npc: "Target · NPC picks", target_cytoplasm: "Target · cytoplasm picks",
  reference_field: "Cy3B reference · full field", reference_npc: "Cy3B reference · NPC picks",
  reference_cytoplasm: "Cy3B reference · cytoplasm picks",
};
const finite = value => typeof value === "number" && Number.isFinite(value);
const number = (value, digits = 2) => finite(value)
  ? value.toLocaleString("en-US", {maximumFractionDigits: digits}) : "—";

function node(tag, text, className) {
  const out = document.createElement(tag);
  if (text !== undefined) out.textContent = String(text);
  if (className) out.className = className;
  return out;
}

function anchor(text, url, title) {
  const a = node("a", text);
  const destination = new URL(url, import.meta.url);
  if (!["http:", "https:"].includes(destination.protocol)) throw Error("Invalid source address");
  a.href = destination.href;
  if (title) a.title = title;
  return a;
}

function choices(selector, items, selected) {
  selector.replaceChildren(...items.map(([value, label]) => {
    const option = node("option", label);
    option.value = String(value);
    return option;
  }));
  selector.value = String(selected);
}

function sourceList(items) {
  const list = node("div", undefined, "palette-links");
  for (const [name, url] of items) list.append(anchor(name, url));
  return list;
}

function table(headers, list, label) {
  const wrapper = node("div", undefined, "palette-scroll");
  wrapper.tabIndex = 0;
  wrapper.setAttribute("role", "region");
  wrapper.setAttribute("aria-label", label);
  const element = node("table", undefined, "palette-table");
  const caption = node("caption", label, "sr");
  const head = element.createTHead().insertRow();
  for (const h of headers) {
    const th = node("th", h); th.scope = "col"; head.append(th);
  }
  const body = element.createTBody();
  for (const row of list) {
    const tr = body.insertRow();
    for (const [i, content] of row.entries()) {
      const cell = node(i ? "td" : "th");
      if (!i) cell.scope = "row";
      if (content instanceof Node) cell.append(content); else cell.textContent = content;
      tr.append(cell);
    }
  }
  element.prepend(caption);
  wrapper.append(element);
  return wrapper;
}

async function readData(file) {
  const response = await fetch(dataURL(file));
  if (!response.ok) throw new Error(file + " returned HTTP " + response.status);
  return response.json();
}

const skeleton = `
  <div class="palette-controls" role="group" aria-label="Measurements to display">
    <label>Excitation channel<select id="palette-wave"></select></label>
    <label>Specimen<select id="palette-system"><option value="origami">DNA origami</option><option value="npc">Nuclear pore complexes</option></select></label>
    <label>Measurement<select id="palette-metric"></select></label>
    <label>Dye to inspect<select id="palette-dye"></select></label>
  </div>
  <div class="palette-grid">
    <article class="palette-card palette-plot-panel" aria-labelledby="palette-plot-title">
      <div class="palette-panel-head"><h2 id="palette-plot-title"></h2><span class="palette-panel-meta" id="palette-plot-unit"></span></div>
      <p id="palette-plot-context" class="palette-description"></p>
      <div class="palette-key"><span><i class="pal-author"></i>Published mean</span><span><i class="pal-raw"></i>HDF5 mean</span></div>
      <div id="palette-series"></div>
    </article>
    <article class="palette-card" aria-labelledby="palette-repeat-title">
      <div class="palette-panel-head"><h2 id="palette-repeat-title"></h2><span class="palette-panel-meta" id="palette-repeat-tag"></span></div>
      <p id="palette-repeat-context" class="palette-description"></p>
      <div id="palette-repeats"></div><div id="palette-dye-note"></div>
      <details class="palette-details" id="palette-agreement"><summary id="palette-agreement-summary">Source agreement</summary><div id="palette-agreement-body"></div></details>
      <details class="palette-details"><summary>Files, group definitions &amp; quality checks</summary><div id="palette-raw-inspect"></div></details>
    </article>
  </div>
  <article class="palette-card palette-evaluation" aria-labelledby="palette-evaluation-title">
    <div class="palette-panel-head"><h2 id="palette-evaluation-title">Held-out NPC comparison</h2><span id="palette-gate" role="status"></span></div>
    <p class="palette-description" id="palette-eval-summary"></p>
    <div class="palette-eval-controls"><label for="palette-eval-source">Evaluation input</label><select id="palette-eval-source"><option value="raw">Deposited HDF5 reanalysis</option><option value="published">Published source tables</option></select></div>
    <div id="palette-eval-results"></div>
  </article>
  <article class="palette-card palette-source-card"><h2>Data and protocol</h2><p class="palette-description" id="palette-coverage"></p><div id="palette-sources"></div></article>
  <p class="sr" aria-live="polite" id="palette-announcement"></p>`;

function replicateValue(row, metric) {
  const value = row?.metrics?.[metric];
  return finite(value) ? value : null;
}

function principalMean(rows, metric) {
  const values = rows.filter(r => r.principal !== false && r.repeat <= 3)
    .map(r => replicateValue(r, metric)).filter(finite);
  return values.length ? values.reduce((a, b) => a + b, 0) / values.length : null;
}

function measureRows(catalog, raw, state) {
  return catalog.dyes.filter(dye => Number(dye.wavelength) === state.wave).map(dye => {
    const author = dye.systems?.[state.system]?.metrics?.[state.metric];
    const repeats = (raw?.replicates || []).filter(r => r.dye === dye.id && r.system === state.system);
    return {dye, author, repeats, rawMean: principalMean(repeats, state.metric)};
  });
}

function metricDetails(catalog, raw, system) {
  const keys = new Set();
  for (const dye of catalog.dyes) {
    for (const [key, metric] of Object.entries(dye.systems?.[system]?.metrics || {})) {
      if (finite(metric.mean) || metric.points?.length) keys.add(key);
    }
  }
  for (const r of raw?.replicates || []) {
    if (r.system === system) for (const [k, v] of Object.entries(r.metrics || {})) {
      if (finite(v) && (catalog.metrics[k] || rawMetrics[k])) keys.add(k);
    }
  }
  return [...keys].map(k => [k, (catalog.metrics[k] || rawMetrics[k]).label]);
}

function drawPlot(root, rows, metric, state, selectDye) {
  root.querySelector("#palette-plot-title").textContent = metric.label;
  root.querySelector("#palette-plot-unit").textContent = metric.unit;
  root.querySelector("#palette-plot-context").textContent =
    specimens[state.system] + " · " + state.wave + " nm · " + rows.length + " dyes. Select a dye to inspect its repeats.";
  const values = rows.flatMap(r => [r.author?.mean, r.rawMean]).filter(finite);
  const min = Math.min(0, ...values);
  const max = Math.max(0, ...values);
  const range = max - min || 1;
  const holder = root.querySelector("#palette-series");
  holder.replaceChildren();
  if (!values.length) {
    holder.append(node("p", "This measurement has no reported values for this specimen and channel.", "palette-empty"));
    return;
  }
  const header = node("div", undefined, "palette-series-head");
  for (const title of ["Dye", "Measurement", "Published", "HDF5"]) header.append(node("span", title));
  holder.append(header);
  for (const row of rows) {
    const button = node("button", undefined, "palette-series-row");
    button.type = "button";
    button.dataset.dye = row.dye.id;
    const isSelected = row.dye.id === state.dye;
    button.setAttribute("aria-pressed", String(isSelected));
    const original = row.author?.mean;
    button.setAttribute("aria-label", row.dye.name + ", published " + number(original) + ", HDF5 " + number(row.rawMean) + " " + metric.unit);
    const track = node("span", undefined, "palette-track");
    track.setAttribute("aria-hidden", "true");
    for (const [value, cls] of [[original, "pal-author"], [row.rawMean, "pal-raw"]]) {
      if (finite(value)) {
        const mark = node("i", undefined, cls);
        mark.style.left = ((value - min) / range) * 100 + "%";
        track.append(mark);
      }
    }
    button.append(node("span", row.dye.name, "palette-dye-name"), track,
      node("span", number(original), "palette-number"), node("span", number(row.rawMean), "palette-number"));
    button.addEventListener("click", () => selectDye(row.dye.id));
    holder.append(button);
  }
  const scale = node("div", undefined, "palette-scale");
  const boundaries = node("span");
  boundaries.append(node("span", number(min)), node("span", number(max)));
  scale.append(node("span"), boundaries);
  holder.append(scale);
}

function rawDetail(row, state) {
  const content = node("div", undefined, "palette-raw-details");
  content.append(node("p", state.system === "origami"
    ? "The full field supports NeNA. Docking-site picks group the binding events used for brightness, kinetics and loss."
    : "NPC and cytoplasm picks are subsets of their full field. Paired Cy3B reference fields and picks are used for the specificity ratio."));
  content.append(node("p", "A repeat is an acquisition. Row counts describe overlapping subsets and can include multiple localizations of one binding event."));
  if (!row.repeats.length) content.append(node("p", "No HDF5 record was loaded for this dye and specimen."));
  for (const record of [...row.repeats].sort((a,b) => a.repeat-b.repeat)) {
    const detail = node("details", undefined, "palette-file-details");
    detail.append(node("summary", "Repeat " + record.repeat + " · " + (record.files?.length || 0) + " files"));
    const list = node("ul", undefined, "palette-file-list");
    for (const file of record.files || []) {
      const li = node("li");
      li.append(node("strong", roleNames[file.role] || file.role),
        node("span", number(file.rows, 0) + " rows", "palette-file-count"),
        node("code", file.path));
      list.append(li);
    }
    detail.append(list);
    if (record.counts && Object.keys(record.counts).length) {
      const counts = node("dl", undefined, "palette-counts");
      for (const [key, val] of Object.entries(record.counts)) {
        if (!finite(val)) continue;
        const label = key.replace(/([a-z])([A-Z])/g, "$1 $2").replaceAll("_", " ");
        counts.append(node("dt", label), node("dd", number(val, 0)));
      }
      detail.append(counts);
    }
    for (const [role, fit] of Object.entries(record.fieldFits || {})) {
      if (finite(fit.pairCount)) detail.append(node("p", (roleNames[role] || role) +
        " · NeNA fit used " + number(fit.pairCount, 0) + " adjacent-frame pairs" +
        (finite(fit.relativeResidualL2) ? ", relative L2 residual " + number(fit.relativeResidualL2, 4) : "") + ".", "palette-qc"));
    }
    const notes = Array.isArray(record.qc) ? record.qc : [];
    for (const s of notes) detail.append(node("p", String(s), "palette-qc"));
    if (!notes.length) detail.append(node("p", "No analysis note recorded for this repeat.", "palette-qc"));
    content.append(detail);
  }
  return content;
}

function drawRepeats(root, catalog, row, state, metric) {
  const source = catalog.sources[state.system];
  root.querySelector("#palette-repeat-title").textContent = row.dye.name;
  root.querySelector("#palette-repeat-tag").textContent = state.wave + " nm";
  root.querySelector("#palette-repeat-context").textContent = metric.label + " (" + metric.unit + ") · repeat measurements";
  const points = row.author?.points || [];
  const repeatIds = [...new Set([...points, ...row.repeats].map(r => r.repeat))].sort((a,b) => a-b);
  const cellLink = address => anchor(source.sheet + "!" + address, source.local, "Open the Figure " + source.figure + " workbook");
  const lines = repeatIds.map(repeat => {
    const p = points.find(x => x.repeat === repeat);
    const rawRecord = row.repeats.find(x => x.repeat === repeat);
    return [String(repeat) + (repeat > 3 ? " · extra" : ""), number(p?.value, 4), number(replicateValue(rawRecord, state.metric), 4), p?.cell ? cellLink(p.cell) : "—"];
  });
  const body = root.querySelector("#palette-repeats");
  body.replaceChildren(table(["Repeat", "Published", "HDF5", "Source cell"], lines, row.dye.name + " " + metric.label + " by acquisition"));
  const notes = root.querySelector("#palette-dye-note");
  notes.replaceChildren();
  if (row.author?.meanCell) {
    const p = node("p", undefined, "palette-note");
    p.append("Published mean ", node("strong", number(row.author.mean)), " " + metric.unit + " · ", cellLink(row.author.meanCell));
    notes.append(p);
  }
  notes.append(node("p", "HDF5 means use available principal repeats 1–3. A dash indicates a missing value.", "palette-note"));
  if (repeatIds.some(i => i > 3)) notes.append(node("p",
    "The additional Cy3B acquisition is listed separately. The published mean can include it, depending on the metric.", "palette-note"));
  const sample = row.dye.systems?.[state.system] || {};
  const conditions = [sample.buffer ? "Buffer " + sample.buffer : null, finite(sample.concentrationPm) ? number(sample.concentrationPm) + " pM imager" : null].filter(Boolean);
  if (conditions.length) notes.append(node("p", conditions.join(" · "), "palette-note"));
  root.querySelector("#palette-raw-inspect").replaceChildren(rawDetail(row, state));
}

function drawAgreement(root, catalog, audit, state, metric) {
  const summary = root.querySelector("#palette-agreement-summary");
  const body = root.querySelector("#palette-agreement-body");
  const source = catalog.sources[state.system];
  body.replaceChildren();
  if (!catalog.dyes.some(d => d.systems?.[state.system]?.metrics?.[state.metric]?.points?.length)) {
    summary.textContent = "Source agreement · HDF5 only";
    body.append(node("p", "The authors’ workbook has no replicate cells for this measurement. The values shown come from the deposited localization tables."));
    return;
  }
  const match = audit?.metrics?.find(m => m.system === state.system && m.metric === state.metric);
  if (!match) {
    summary.textContent = "Source agreement · audit unavailable";
    body.append(node("p", "Per-repeat values and workbook cells are available above. The agreement summary for this measurement could not be loaded."));
    return;
  }
  const unit = metric.unit === "%" ? "percentage points" : metric.unit === "ratio" ? "ratio units" : metric.unit;
  const precise = v => finite(v) && v !== 0 && Math.abs(v) < .0001 ? v.toExponential(2) : number(v, 4);
  summary.textContent = "Source agreement · " + match.withinHalfGrid + "/" + match.comparable + " within half-grid";
  body.append(node("p", "Across all excitation channels for " + specimens[state.system] + ", " +
    match.comparable + " of " + match.sourceCells + " workbook repeat values have a corresponding HDF5 estimate."));
  if (finite(match.medianAbsoluteDifference)) body.append(node("p", "The median absolute difference is " +
    precise(match.medianAbsoluteDifference) + " " + unit + ". The mean signed difference (HDF5 minus published) is " +
    precise(match.meanSignedDifference) + " " + unit + "."));
  body.append(node("p", "The descriptive grid has a step of " + precise(match.gridUnit) + " " + unit +
    ". The agreement count includes absolute differences up to half a step. This grid describes numerical agreement with the workbook."));
  const worst = match.largestDifferences?.[0];
  if (worst) {
    const name = catalog.dyes.find(d => d.id === worst.dye)?.name || worst.dye;
    const p = node("p");
    p.append(node("strong", "Largest observed difference"), " · " + name + ", repeat " + worst.repeat +
      ". Published " + precise(worst.published) + " and HDF5 " + precise(worst.recomputed) + " " + metric.unit +
      " (absolute difference " + precise(Math.abs(worst.difference)) + " " + unit + "). ",
      anchor(source.sheet + "!" + worst.cell, source.local, "Open the Figure " + source.figure + " workbook"));
    body.append(p);
  }
}

function drawEvaluation(root, catalog, evaluation, which) {
  const byId = new Map(catalog.dyes.map(x => [x.id, x.name]));
  const names = ids => Array.isArray(ids) && ids.length ? ids.map(id => byId.get(id) || id).join(" + ") : "—";
  const raw = evaluation?.raw;
  const gate = root.querySelector("#palette-gate");
  const approved = raw?.status === "evaluated" && raw.gatePassed === true;
  gate.dataset.state = approved ? "passed" : "held";
  gate.textContent = approved ? "Raw-data gate passed" : "Recommendation withheld";
  const available = raw?.status === "evaluated";
  root.querySelector("#palette-eval-summary").textContent = approved
    ? "The rule met the prespecified criteria on held-out NPC repeat 3. Only pairs from the deposited HDF5 reanalysis qualify as suggestions for this DNA-PAINT panel."
    : available ? "The rule did not meet the prespecified criteria on held-out NPC repeat 3. The pairs below show the retrospective comparison."
      : "The raw-data evaluation is unavailable. The published-table replay can be inspected independently.";
  const view = root.querySelector("#palette-eval-results");
  view.replaceChildren();
  const stats = evaluation?.[which];
  if (stats?.status !== "evaluated") {
    view.append(node("p", (stats?.reasons || []).join(" ") || "Evaluation is unavailable for this input.", "palette-empty"));
    return;
  }
  const summary = node("div", undefined, "palette-eval-metrics");
  for (const [value, label] of [[number(stats.macroRegretNm) + " nm", "Rule · mean channel regret"],
    [number(stats.baselineMacroRegretNm) + " nm", "Fixed pair · mean channel regret"],
    [String(stats.improvedChannels ?? "—") + " / " + (stats.channels?.length || 0), "Channels improved"]]) {
    const item = node("div");
    item.append(node("strong", value), node("span", label));
    summary.append(item);
  }
  view.append(summary);
  const results = (stats.channels || []).map(c => [
    number(c.wavelength,0) + " nm", names(c.pair), names(c.baselinePair),
    number(c.meanNm), number(c.baselineMeanNm), number(c.regretNm), number(c.baselineRegretNm),
  ]);
  view.append(table(["Channel", "Rule pair", "Fixed pair", "Rule NeNA (nm)", "Fixed NeNA (nm)", "Rule regret (nm)", "Fixed regret (nm)"],
    results, "Held-out NPC pair scores and regrets against the best observed pair within each channel"));
  view.append(node("p", "Both policies chose two dyes using repeats 1 and 2. The rule can combine NPC and origami ranks. The fixed policy uses NPC ranks alone. Lower NeNA is better. Regret is the gap in mean NeNA to the best pair observed on NPC repeat 3.", "palette-note"));
  if (which === "published") view.append(node("p", "The published-table replay " + (stats.gatePassed ? "meets" : "falls short of") + " the same criteria. Only the deposited HDF5 evaluation determines whether suggestions are shown.", "palette-note"));
  const oracle = node("details", undefined, "palette-details");
  oracle.append(node("summary", "Best observed pairs and decision criteria"));
  for (const c of stats.channels || []) oracle.append(node("p", number(c.wavelength, 0) + " nm · " + names(c.oraclePair) + " · " + number(c.oracleMeanNm) + " nm"));
  if (finite(stats.weight)) oracle.append(node("p", "Repeat 2 selected an origami-rank weight of " + number(stats.weight) + " from 0, 0.5 and 1."));
  oracle.append(node("p", "The gate requires lower mean regret across the three equally weighted channels and lower regret in at least two of them."));
  for (const reason of (stats.reasons || [])) oracle.append(node("p", reason));
  oracle.append(anchor("Read the locked evaluation protocol", protocolURL));
  view.append(oracle);
}

async function mount(container) {
  container.classList.add("palette-workbench");
  container.setAttribute("aria-busy", "true");
  container.replaceChildren(node("p", "Loading DNA-PAINT measurements…", "palette-empty"));
  const results = await Promise.allSettled(["palette.json","palette-raw.json","palette-evaluation.json","palette-audit.json"].map(readData));
  const [catalog, raw, evaluation, audit] = results.map(result => result.status === "fulfilled" ? result.value : null);
  if (!catalog?.dyes?.length || !catalog.metrics || !catalog.sources) {
    container.replaceChildren(node("p", "The source tables could not be loaded.", "palette-empty"));
    const again = node("button", "Retry loading measurements", "button secondary");
    again.type = "button";
    again.onclick = () => { running.delete(container); renderPalette(container); };
    container.append(again);
    container.removeAttribute("aria-busy");
    return;
  }
  container.innerHTML = skeleton;
  const get = id => container.querySelector("#palette-" + id);
  const wavelengths = [...new Set(catalog.dyes.map(d => Number(d.wavelength)))].sort((a,b) => a-b);
  const state = {wave: wavelengths.includes(640) ? 640 : wavelengths[0], system: "origami", metric: "nenaNm", dye: "janeliafluor646"};
  choices(get("wave"), wavelengths.map(w => [w, w + " nm"]), state.wave);
  function update() {
    const available = metricDetails(catalog, raw, state.system);
    if (!available.some(([k]) => k === state.metric)) state.metric = available[0][0];
    choices(get("metric"), available, state.metric);
    const rows = measureRows(catalog, raw, state);
    if (!rows.some(r => r.dye.id === state.dye)) state.dye = rows[0].dye.id;
    choices(get("dye"), rows.map(r => [r.dye.id, r.dye.name]), state.dye);
    const metric = catalog.metrics[state.metric] || rawMetrics[state.metric];
    drawPlot(container, rows, metric, state, id => {
      state.dye = id;
      update();
      [...container.querySelectorAll(".palette-series-row")].find(el => el.dataset.dye === id)?.focus({preventScroll:true});
    });
    const selected = rows.find(r => r.dye.id === state.dye);
    drawRepeats(container, catalog, selected, state, metric);
    if (typeof raw?.metrics?.[state.metric] === "string") {
      get("dye-note").append(node("p", "HDF5 metric definition · " + raw.metrics[state.metric], "palette-note"));
    }
    drawAgreement(container, catalog, audit, state, metric);
    get("announcement").textContent = selected.dye.name + " · " + metric.label + " · " + specimens[state.system];
  }
  for (const key of ["wave","system","metric","dye"]) get(key).addEventListener("change", event => {
    state[key] = key === "wave" ? Number(event.target.value) : event.target.value;
    update();
  });
  update();
  const input = get("eval-source");
  input.value = evaluation?.raw?.status === "evaluated" ? "raw" :
    evaluation?.published?.status === "evaluated" ? "published" : "raw";
  input.addEventListener("change", () => drawEvaluation(container, catalog, evaluation, input.value));
  drawEvaluation(container, catalog, evaluation, input.value);
  const n = raw?.coverage;
  get("coverage").textContent = n ? number(n.completedCases,0) + " of " + number(n.expectedCases,0) +
    " deposited acquisition groups processed · " + number(n.totalFiles,0) + " localization files. Related field and pick files overlap."
    : "Workbook values are available. The deposited localization reanalysis could not be loaded.";
  if (audit) get("coverage").append(" The audit verified " + number(audit.summaryCellsReproduced,0) +
    " workbook means and standard deviations from the published repeats and traced " +
    number(audit.mainFigureHdf5,0) + " main-figure localization files in an archive of " + number(audit.hdf5Members,0) + " HDF5 files.");
  get("sources").append(sourceList([
    ["Steen et al. paper", catalog.article], ["Localization archive", catalog.record],
    ["Figure 2 workbook", catalog.sources.origami.local], ["Figure 3 workbook", catalog.sources.npc.local],
    ["Archive manifest JSON", dataURL("palette-manifest.json")],
    ["Published measurements JSON", dataURL("palette.json")],
    ["Raw-derived measurements JSON", dataURL("palette-raw.json")], ["Source agreement audit JSON", dataURL("palette-audit.json")],
    ["HDF5 reanalysis protocol", rawProtocolURL], ["Evaluation results JSON", dataURL("palette-evaluation.json")], ["Evaluation protocol", protocolURL],
  ]));
  container.removeAttribute("aria-busy");
}

export function renderPalette(container) {
  if (!running.has(container)) running.set(container, mount(container));
  return running.get(container);
}
