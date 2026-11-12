import {mountTrial} from './probe-form.js?v=unseen-1';

const metrics = {
  useful: {label: 'Useful pixels', unit: '%', delta: 'percentage points'},
  snr: {label: 'Signal-to-noise ratio', unit: 'SNR', delta: 'SNR'},
  tau: {label: 'Decorrelation time', unit: 'frames', delta: 'frames'},
};
const ids = id => document.getElementById('probe-' + id);
const fmt = (value, digits = 2) => Number(value).toLocaleString('en-US', {maximumFractionDigits: digits});
const ratio = value => value.toExponential(2);

function el(tag, text, className) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (className) node.className = className;
  return node;
}

function sourceLink(text, href) {
  const link = el('a', text);
  link.href = href;
  return link;
}

function sourceCell(data, src) {
  const link = sourceLink(src.cells, data.sources[src.workbook].file);
  link.title = src.sheet + ' · ' + src.cells;
  return link;
}

function setOptions(select, choices, selected) {
  select.replaceChildren(...choices.map(([value, label]) => new Option(label, value)));
  select.value = selected;
}

function cells(values) {
  const row = el('tr');
  for (const value of values) {
    const cell = el('td');
    cell.append(value);
    row.append(cell);
  }
  return row;
}

function structures(pair, data) {
  const draw = (dye, i) => {
    const figure = el('figure', undefined, 'probe-molecule series-' + i);
    const head = el('div', undefined, 'probe-molecule-name');
    head.append(el('i'), el('strong', dye.id));
    const picture = el('img');
    picture.src = dye.structure.image;
    picture.alt = dye.id + ', closed hydroxymethyl Si-rhodamine with ' + dye.auxochrome.toLowerCase() + ' substituents.';
    const caption = el('figcaption');
    const group = el('strong', dye.auxochrome);
    if (dye.auxochrome.includes(')')) {
      const [substituent, ring] = dye.auxochrome.split(')');
      group.replaceChildren(substituent + ')', el('wbr'), ring);
    }
    caption.append(group,
      el('span', 'Parent ' + dye.compound + ' · ligand ' + dye.ligand),
      sourceLink('View original structure ↗', data.si.url + '#page=' + dye.structure.pdfPage));
    figure.append(head, picture, caption);
    return figure;
  };
  ids('molecules').replaceChildren(...pair.map(draw));
  ids('property-head').replaceChildren(...['Property', ...pair.map(d => d.id)].map(s => el('th', s)));
  const qy = d => fmt(d.quantumYield.mean, 3);
  const labels = [
    ['Parent apparent pKa', d => d.pKaFigure1],
    ['Absorption / emission (nm)', d => d.absorptionNm + ' / ' + d.emissionNm],
    ['Fluorescence quantum yield', qy],
    ['H2B ligand ON:OFF ratio', d => ratio(d.ligandRatio.mean)],
  ];
  ids('properties').replaceChildren(...labels.map(([name, get]) => cells([name, ...pair.map(get)])));
}

function cellular(pair, ctx, metric, data) {
  const spec = metrics[metric];
  const states = ['live', 'fixed'];
  const readings = state => pair.map(d => ctx.dyes[d.id][state][metric]);
  const all = states.flatMap(readings);
  const ceiling = metric === 'useful' ? 100 : Math.max(...all.map(v => v.mean + v.sem)) * 1.04;
  const cellRows = [];
  const barGroups = [];
  for (const state of states) {
    const group = el('div', undefined, 'probe-state');
    group.append(el('h3', state === 'live' ? 'Live cells' : 'Fixed cells'));
    for (const [i, v] of readings(state).entries()) {
      const dye = pair[i];
      const row = el('div', undefined, 'probe-bar-row series-' + i);
      const track = el('div', undefined, 'probe-track');
      track.setAttribute('aria-hidden', 'true');
      const bar = el('i'), whisker = el('b');
      bar.style.width = 100 * v.mean / ceiling + '%';
      const low = Math.max(0, v.mean - v.sem), high = Math.min(ceiling, v.mean + v.sem);
      whisker.style.left = 100 * low / ceiling + '%';
      whisker.style.width = 100 * (high - low) / ceiling + '%';
      track.append(bar, whisker);
      row.append(el('span', dye.id), track, el('strong', fmt(v.mean) + (metric === 'useful' ? '%' : '')));
      group.append(row);
      cellRows.push(cells([dye.id, state, String(v.mean), String(v.sem), sourceCell(data, v.source)]));
    }
    barGroups.push(group);
  }
  ids('context-bars').replaceChildren(...barGroups);
  ids('context-title').textContent = spec.label + ' in live and fixed cells';
  ids('cell-unit').textContent = ctx.name + ' · ' + spec.unit + ' · mean ± SEM, n = ' + ctx.nCells + ' cells';
  ids('context-table').replaceChildren(...cellRows);
  ids('cell-source').textContent = ctx.fig + ' · ' + ctx.protocol + '.';
  ids('lysosome-note').hidden = ctx !== data.contexts.lysosome;
  const difference = state => {
    const [a, b] = readings(state);
    return a.mean - b.mean;
  };
  const change = difference('live') * difference('fixed') < 0;
  const sign = v => (v > 0 ? '+' : '') + fmt(v, 2);
  const outcome = ids('outcome');
  outcome.dataset.orderChange = String(change);
  outcome.replaceChildren(el('strong', change ? 'The observed order reverses.' : 'The observed order is preserved.'),
    document.createTextNode(' ' + pair[0].id + ' minus ' + pair[1].id + ' is ' + sign(difference('live')) +
      ' ' + spec.delta + ' live and ' + sign(difference('fixed')) + ' after fixation.'));
}

function drawTitration(canvas, pair, selected) {
  const width = canvas.clientWidth, height = canvas.clientHeight;
  if (!width || !height) return;
  const scale = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * scale);
  canvas.height = Math.round(height * scale);
  const g = canvas.getContext('2d');
  g.scale(scale, scale);
  const left = 36, right = 12, top = 10, bottom = 26;
  const x = v => left + (v - 1) / 6 * (width - left - right);
  const y = v => height - bottom - v * (height - bottom - top);
  g.font = '11px monospace';
  g.fillStyle = '#c6d8cf';
  g.textAlign = 'right';
  for (const tick of [0, 0.25, 0.5, 0.75, 1]) {
    g.strokeStyle = '#42645b';
    g.setLineDash([]);
    g.beginPath(); g.moveTo(left, y(tick)); g.lineTo(width - right, y(tick)); g.stroke();
    g.fillText(String(tick), left - 6, y(tick) + 4);
  }
  g.textAlign = 'center';
  for (const tick of [1, 2, 3, 4, 5, 6, 7]) g.fillText(String(tick), x(tick), height - 8);
  g.strokeStyle = '#b8cabe';
  g.setLineDash([3, 4]);
  g.beginPath(); g.moveTo(x(selected), top); g.lineTo(x(selected), height - bottom); g.stroke();
  for (const [i, dye] of pair.entries()) {
    const color = i ? '#7ed1e1' : '#c8ee92';
    g.strokeStyle = color;
    g.fillStyle = color;
    g.lineWidth = 2;
    g.setLineDash(i ? [5, 3] : []);
    g.beginPath();
    dye.titration.forEach((p, j) => j ? g.lineTo(x(p.pH), y(p.mean)) : g.moveTo(x(p.pH), y(p.mean)));
    g.stroke();
    g.setLineDash([]);
    for (const p of dye.titration) {
      const size = p.pH === selected ? 4.2 : 2.2;
      g.beginPath();
      if (i) g.rect(x(p.pH) - size, y(p.mean) - size, size * 2, size * 2);
      else g.arc(x(p.pH), y(p.mean), size, 0, 2 * Math.PI);
      g.fill();
    }
  }
}

function titrationTable(pair, data, previousPH) {
  const points = [...new Set(pair.flatMap(d => d.titration.map(p => p.pH)))].sort((a, b) => a-b);
  const ph = points.reduce((a, b) => Math.abs(a-previousPH) <= Math.abs(b-previousPH) ? a : b);
  const picker = ids('ph');
  setOptions(picker, points.map(v => [String(v), v.toFixed(2)]), String(ph));
  const records = [];
  for (const dye of pair) {
    for (const p of dye.titration) {
      const src = {workbook: '9', sheet: data.sources[9].sheet, cells: p.cells};
      records.push(cells([dye.id, String(p.pH), ...p.readings.map(String), sourceCell(data, src)]));
    }
  }
  ids('titration-table').replaceChildren(...records);
  ids('curve-legend').replaceChildren(...pair.map((d, i) => {
    const key = el('span');
    key.append(el('i', undefined, 'line key-' + (i ? 'b' : 'a')), d.id);
    return key;
  }));
}

function readingAtPH(pair) {
  const ph = Number(ids('ph').value);
  const readings = pair.map(d => {
    const point = d.titration.find(p => p.pH === ph);
    return d.id + '  ' + (point ? fmt(point.mean, 4) : 'no reading');
  });
  ids('ph-reading').textContent = readings.join('   ·   ');
  drawTitration(ids('titration'), pair, ph);
}

function drawCalibration(canvas, calibration) {
  const w = canvas.clientWidth, h = canvas.clientHeight;
  if (!w || !h) return;
  const dpi = window.devicePixelRatio || 1;
  canvas.width = Math.round(w * dpi); canvas.height = Math.round(h * dpi);
  const g = canvas.getContext('2d');
  g.scale(dpi, dpi);
  const x = p => 44 + (p - 1) / 3.6 * (w - 57);
  const y = log => h - 24 - (log + 4.15) / 1.36 * (h - 36);
  g.font = '10px monospace'; g.textAlign = 'right';
  for (const [log, label] of [[-4, '0.0001'], [-3.5, '0.00032'], [-3, '0.001']]) {
    g.strokeStyle = '#d9e3d8'; g.fillStyle = '#50665a';
    g.beginPath(); g.moveTo(45, y(log)); g.lineTo(w - 9, y(log)); g.stroke();
    g.fillText(label, 42, y(log) + 3);
  }
  g.textAlign = 'center';
  for (const tick of [1, 2, 3, 4]) g.fillText(String(tick), x(tick), h - 7);
  g.strokeStyle = '#688777'; g.lineWidth = 1.6; g.setLineDash([4,3]); g.beginPath();
  const {pKaRange, intercept, slope} = calibration.fit;
  g.moveTo(x(pKaRange[0]), y(intercept + slope * pKaRange[0]));
  g.lineTo(x(pKaRange[1]), y(intercept + slope * pKaRange[1])); g.stroke(); g.setLineDash([]);
  for (const p of calibration.pairs) {
    g.beginPath(); g.fillStyle = '#226a51'; g.arc(x(p.pKa), y(Math.log10(p.ratio)), 4, 0, Math.PI * 2); g.fill();
  }
}

export async function renderProbeDesign(root) {
  const [data, benchmark] = await Promise.all(['probe-design', 'analogue-benchmark'].map(async name => {
    const res = await fetch('data/' + name + '.json');
    if (!res.ok) throw new Error(name + ' (' + res.status + ').');
    return res.json();
  }));
  const development = benchmark.sets.find(s => s.mode === 'development');
  root.innerHTML = `
    <div class="card probe-controls">
      <label>Cellular assay<select id="probe-context" aria-label="Cellular assay"></select></label>
      <label>First dye<select id="probe-a"></select></label>
      <label>Second dye<select id="probe-b"></select></label>
      <label>SOFI readout<select id="probe-metric"></select></label>
    </div>
    <div class="probe-pair-grid">
      <article class="card probe-card"><div class="eyebrow">FREE DYE · CLOSED FORM</div><h2>What changed in the molecule?</h2>
        <div id="probe-molecules"></div>
        <div class="table-wrap"><table class="probe-properties"><thead><tr id="probe-property-head"></tr></thead><tbody id="probe-properties"></tbody></table></div>
        <p class="footnote">Spectra and quantum yield at pH 3 with 0.1% SDS. The ON:OFF ratio comes from the corresponding ligand on H2B–HaloTag in fixed COS-7 cells.</p>
      </article>
      <article class="card probe-card probe-cell-card"><div class="eyebrow">CELLULAR SOFI · HALOTAG IN U2OS</div><h2 id="probe-context-title"></h2>
        <p class="probe-sub" id="probe-cell-unit"></p><div id="probe-context-bars"></div>
        <p id="probe-outcome" role="status"></p>
        <p id="probe-cell-source" class="footnote"></p><p class="footnote" id="probe-lysosome-note" hidden>The methods identify CD63 and the figure captions identify LAMP. Results are shown at the compartment level. The two ligands had different labeling concentrations.</p>
        <details><summary>Means, errors and source cells</summary><div class="table-wrap"><table><thead><tr><th>Dye</th><th>State</th><th>Mean</th><th>SEM</th><th>Cells in sheet</th></tr></thead><tbody id="probe-context-table"></tbody></table></div><p class="footnote">These are cell-level summary statistics. Individual matched-cell differences are unavailable.</p></details>
        <p class="footnote">For H2B localization counts at the 20 nm cutoff, see <a href="#precision">precision evidence</a>.</p>
      </article>
    </div>
    <div class="probe-buffer-grid">
      <article class="plot-card probe-plot"><div class="plot-label"><i></i> FREE DYE · TWO TECHNICAL READINGS PER POINT</div>
        <div class="plot-head"><h2>Absorbance across pH</h2></div><div class="plot-unit">Normalized absorbance</div>
        <canvas id="probe-titration" role="img" aria-label="Normalized absorbance versus pH for the two selected free dyes. Exact readings are in the titration table."></canvas>
        <div class="x-axis">Buffer pH</div><div id="probe-curve-legend" class="chart-legend"></div>
        <div class="probe-ph-picker"><label for="probe-ph">Inspect measured pH</label><select id="probe-ph"></select><output id="probe-ph-reading" for="probe-ph"></output></div>
        <details><summary>Exact titration readings</summary><div class="table-wrap"><table><thead><tr><th>Dye</th><th>pH</th><th>Read 1</th><th>Read 2</th><th>Source cells</th></tr></thead><tbody id="probe-titration-table"></tbody></table></div></details>
      </article>
      <article class="card probe-card" id="probe-validation"><div class="eyebrow">FIXED H2B · RETROSPECTIVE CHECK</div><h2>Can pKa predict blinking?</h2>
        <p class="probe-sub">We held out one dye at a time and predicted its ligand’s mean ON:OFF time ratio from the free parent’s apparent pKa.</p>
        <figure class="probe-fit"><figcaption>Mean ligand ON:OFF ratio · log scale</figcaption>
          <canvas id="probe-calibration" role="img" aria-label="Six published pKa and ratio pairs, with the regression fitted to all six. Exact coordinates follow in the calibration table."></canvas>
          <div class="probe-fit-axis">Parent pKa</div>
        </figure>
        <div id="probe-benchmarks"></div><p class="footnote">Mean absolute log₁₀ error on the same held-out dyes. Lower is better.</p>
        <p id="probe-validation-copy"></p>
        <p><a class="inline-link" href="#unseen-dyes">Compare local interpolation on unseen dyes →</a></p>
        <details><summary>Six calibration pairs and fitted line</summary><p class="footnote" id="probe-equation"></p><div class="table-wrap"><table><thead><tr><th>Dye</th><th>Parent pKa</th><th>Mean ratio</th><th>SD · n = 3</th><th>Source cells</th></tr></thead><tbody id="probe-pairs"></tbody></table></div><p class="footnote">HM-SiR and JF646b have complex titrations. JF614b has only an upper limit on pKa (&lt;1). Extended Data Fig. 4j omits these three.</p></details>
        <details><summary>Inspect each held-out dye</summary><div class="table-wrap"><table><thead><tr><th>Held-out dye</th><th>pKa</th><th>Local error</th><th>Regression error</th><th>Nearest error</th><th>Median error</th></tr></thead><tbody id="probe-folds"></tbody></table></div></details>
      </article>
    </div>
    <article class="card probe-card" id="probe-trial"></article>
    <article class="card probe-card probe-sources"><h2>Measurements &amp; reproducibility</h2>
      <p>Holland, Plutkis, English et al., <em>Nature Methods</em> (2026). Parent structures are cropped from Supplementary Note 2. Spectral peaks and displayed pKa values are transcribed from Fig. 1a. ND means undetermined. The six exact calibration pairs come from Extended Data Fig. 4j.</p>
      <div id="probe-source-links"></div>
    </article>`;
  setOptions(ids('context'), Object.entries(data.contexts).map(([key, c]) => [key, key === 'er' ? 'ER · Sec61β' : c.name]), 'er');
  let pair;
  function update() {
    const ctx = data.contexts[ids('context').value];
    pair = [ids('a').value, ids('b').value].map(name => data.dyes.find(d => d.id === name));
    structures(pair, data);
    cellular(pair, ctx, ids('metric').value, data);
    titrationTable(pair, data, ids('ph').value ? Number(ids('ph').value) : 3.08);
    readingAtPH(pair);
  }
  function contextChanged() {
    const key = ids('context').value, ctx = data.contexts[key];
    const dyes = Object.keys(ctx.dyes).map(name => [name, name]);
    setOptions(ids('a'), dyes, 'JF639b');
    setOptions(ids('b'), dyes, key === 'er' ? 'JF630b' : 'JF614b');
    setOptions(ids('metric'), ctx.metrics.map(k => [k, metrics[k].label]), ctx.metrics[0]);
    update();
  }
  for (const [i, key] of ['a', 'b'].entries()) {
    ids(key).onchange = () => {
      if (ids('a').value === ids('b').value) ids(i ? 'a' : 'b').value = pair[i].id;
      update();
    };
  }
  ids('metric').onchange = () => cellular(pair, data.contexts[ids('context').value], ids('metric').value, data);
  ids('context').onchange = contextChanged;
  ids('ph').onchange = () => readingAtPH(pair);
  contextChanged();
  new ResizeObserver(() => drawTitration(ids('titration'), pair, Number(ids('ph').value))).observe(ids('titration'));
  const c = data.calibration, val = c.validation;
  new ResizeObserver(() => drawCalibration(ids('calibration'), c)).observe(ids('calibration'));
  ids('equation').textContent = 'The dashed line fits all six dyes with log₁₀(ratio) = ' + c.fit.intercept.toFixed(4) +
    ' + ' + c.fit.slope.toFixed(4) + ' × pKa. Each test fold uses only the other five dyes for all four forecasts.';
  ids('pairs').replaceChildren(...c.pairs.map(p => cells([p.id, String(p.pKa), String(p.ratio), String(p.sd), sourceCell(data, p.source)])));
  for (const [label, value] of [['Local interpolation', development.summary.mae.local], ['Nearest pKa', val.nearestMae],
    ['Regression', val.mae], ['Median', val.baselineMae]]) {
    const tile = el('div');
    tile.append(el('span', label), el('strong', value.toFixed(3)));
    ids('benchmarks').append(tile);
  }
  ids('validation-copy').textContent = val.scored + ' of ' + val.total + ' dyes fell inside the remaining calibration range. Local interpolation had the lowest average error, beating nearest pKa for ' +
    development.summary.wins + ' of ' + val.scored + ' development dyes. Both endpoint dyes were outside their folds’ ranges.';
  ids('folds').replaceChildren(...c.folds.map(f => {
    if (f.status === 'scored') {
      const local = development.rows.find(r => r.id === f.held).errors.local;
      return cells([f.held, String(f.pKa), local.toFixed(3), f.error.toFixed(3), f.nearestError.toFixed(3), f.baselineError.toFixed(3)]);
    }
    const row = cells([f.held, String(f.pKa), 'Outside training range']);
    row.lastChild.colSpan = 4; return row;
  }));
  const links = [['Article ↗', data.article], ['Supplementary Note 2 ↗', data.si.url],
    ...Object.values(data.sources).map(s => [s.sheet + ' workbook ↓', s.file]),
    ['Download source bundle (JSON)', 'data/probe-design.json'],
    ['Reproduction and assay protocol ↗', 'https://github.com/kangk1204/blink-atlas/blob/main/analysis/probe-design.md'],
    ['CC BY 4.0 ↗', 'https://creativecommons.org/licenses/by/4.0/']];
  ids('source-links').replaceChildren(...links.map(([label, href]) => sourceLink(label, href)));
  mountTrial(ids('trial'), data);
  ids('jump').disabled = false;
  ids('jump').onclick = () => {ids('trial-heading').focus({preventScroll:true}); ids('trial').scrollIntoView({block:'start'});};
}
