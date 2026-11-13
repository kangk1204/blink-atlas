import {renderRatioAudit, renderRatioEvidence} from './ratio-audit-view.js?v=ratio-2';

const $ = key => document.getElementById('unseen-' + key);
const labels = {local:'Local interpolation', nearest:'Nearest pKa', regression:'Original regression'};
const f3 = x => x == null ? '—' : x.toFixed(3);
const signed = x => (x > 0 ? '+' : '') + f3(x);
const small = x => x.toExponential(3);

function el(tag, text, cls) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (cls) node.className = cls;
  return node;
}
function link(text, href) {
  const node = el('a', text); node.href = href; return node;
}
async function load(name) {
  const res = await fetch('data/' + name + '.json?v=ratio-2');
  if (!res.ok) throw new Error('Could not read ' + name + ' (' + res.status + ').');
  return res.json();
}
function row(entries) {
  const tr = el('tr');
  for (const entry of entries) {
    const td = el('td');
    td.append(entry); tr.append(td);
  }
  return tr;
}
function tile(label, value) {
  const item = el('div');
  item.append(el('span', label), el('strong', value));
  return item;
}

const scopes = {
  primary: 'Primary assay · Predictions are fixed using all six Holland dyes. The test requires independently acquired fixed COS-7 H2B–HaloTag measurements on new dyes.',
  development: 'Development check · The rule was designed with knowledge of this publication. Each test dye is removed before all three methods are fitted.',
  'transfer-frozen': 'External transfer · All six Holland dyes are fixed as the training set. The test laboratory, chemistry and specimen change.',
  'transfer-local': 'Study-specific calibration · Each held-out dye is predicted using only the other dyes measured under this same condition. The external study supplies new training dyes.',
};

function differenceChart(set) {
  const area = $('chart');
  area.replaceChildren();
  const scored = set.rows.filter(r => r.status === 'scored');
  if (!scored.length) {
    area.append(el('p', set.rows.length ? 'No eligible paired errors in this selection.' : 'No independent matched dyes were located in the public sources reviewed.', 'unseen-empty'));
    $('axis').hidden = true; return;
  }
  const limit = Math.max(...scored.map(r => Math.abs(r.delta)), 0.001);
  for (const r of scored) {
    const line = el('div', undefined, 'unseen-bar');
    const plot = el('div', undefined, 'unseen-track'), bar = el('b');
    const width = 50 * Math.abs(r.delta) / limit;
    bar.className = r.delta >= 0 ? 'gain' : 'loss';
    bar.style.left = (r.delta >= 0 ? 50 : 50 - width) + '%';
    bar.style.width = width + '%';
    plot.setAttribute('aria-hidden', 'true'); plot.append(bar);
    line.append(el('span', r.dye), plot, el('strong', signed(r.delta)));
    area.append(line);
  }
  $('axis').hidden = false;
  $('limits').textContent = 'Symmetric scale ±' + f3(limit) + ' log₁₀ units';
}

function verificationTable(set) {
  const body = $('rows'); body.replaceChildren();
  for (const r of set.rows) {
    if (r.status !== 'scored') {
      const reason = [r.status, r.rangeStatus === 'outside training range' ? r.rangeStatus : null].filter(Boolean).join(' · ') +
        (r.range && (r.status === 'outside training range' || r.rangeStatus === 'outside training range') ? ' (' + r.range.join('–') + ')' : '');
      const tr = row([r.dye, String(r.pKa), small(r.observed), reason]);
      tr.lastChild.colSpan = 5; body.append(tr); continue;
    }
    const info = el('details'), summary = el('summary', 'References');
    info.append(summary, el('p','Lower ' + r.lower.join(', ') + ' · upper ' + r.upper.join(', ') + ' · upper weight ' + f3(r.upperWeight)),
      el('p','Nearest ' + r.nearestDyes.join(', ')), el('p','Training identities: ' + r.trainingIds.join(', ')));
    body.append(row([r.dye, String(r.pKa), small(r.observed), ...Object.keys(labels).map(k => small(r.predictions[k])),
      signed(r.delta), info]));
  }
  if (!set.rows.length) {
    const tr = row(['No independent matching measurements.']);
    tr.firstChild.colSpan = 8; body.append(tr);
  }
}

function show(set, meta, matchedSource) {
  const s = set.summary, base = meta.studies[set.study];
  const override = base.conditions?.[set.condition];
  const study = override ? {...base, ...override, sourceUrl: base.conjugateSourceUrl} : base;
  $('scope').textContent = scopes[set.mode];
  $('measure').textContent = study.endpoint || 'Mean ligand ON:OFF time ratio in fixed COS-7 H2B–HaloTag.';
  $('assay-info').textContent = (study.chemistry || '') + ' ' + (study.transferNote || '');
  $('source').replaceChildren(document.createTextNode(study.source + ' '), link('Open source ↗', study.sourceUrl || study.url));
  if (set.mode === 'primary' && matchedSource) $('source').append(' · ', link('Independent measurements ↗', matchedSource));
  $('coverage').textContent = s.scored + ' of ' + s.total + ' dyes scored';
  $('metrics').replaceChildren(...Object.entries(labels).map(([m,label]) => tile(label, f3(s.mae?.[m]))));
  if (s.scored) {
    $('contrast').textContent = 'Nearest-pKa error minus local-model error averaged ' + signed(s.delta) +
      '. Local interpolation had ' + s.wins + ' wins, ' + s.ties + ' ties and ' + s.losses + ' losses.';
    const interval = s.pairedBootstrap95;
    let uncertainty = 'A single scored dye cannot estimate between-dye variability.';
    if (interval) {
      uncertainty = 'Paired dye bootstrap, 95% interval ' + signed(interval[0]) + ' to ' + signed(interval[1]) +
        '. This small-panel interval is descriptive.';
    } else if (s.omitOneRange) {
      uncertainty = 'Omitting one scored dye at a time from the average gives ' + signed(s.omitOneRange[0]) +
        ' to ' + signed(s.omitOneRange[1]) + '. The leave-one-dye-out fits share training data.';
    }
    $('uncertainty').textContent = uncertainty;
    const comparison = s.localVsRegression;
    $('original').textContent = 'Against the original regression, local interpolation had ' + comparison.wins +
      ' wins, ' + comparison.ties + ' ties and ' + comparison.losses + ' losses. Its mean improvement was ' +
      signed(comparison.delta) + (comparison.pairedBootstrap95 ? ' (paired 95% interval ' +
        comparison.pairedBootstrap95.map(signed).join(' to ') + ').' : '.');
  } else {
    $('contrast').textContent = set.mode === 'primary'
      ? 'A gain in the fixed H2B assay remains unconfirmed.'
      : 'A within-range target and a comparable endpoint are required to score a prediction.';
    $('uncertainty').textContent = ''; $('original').textContent = '';
  }
  differenceChart(set); verificationTable(set);
}

function showConjugates(sets, source) {
  const partners = ['IgG', 'HaloTag7', 'Nanobody'];
  const measurements = partners.map(condition => sets.find(s => s.study === 'remmel' && s.condition === condition &&
    s.mode === 'transfer-frozen').rows.find(r => r.id === 'remmel-4'));
  $('conjugates').replaceChildren(...measurements.map((r, i) => row([partners[i], String(r.pKa), small(r.observed), small(r.predictions.local)])));
  const fold = measurements[0].observed / measurements[1].observed;
  $('conjugate-effect').textContent = 'The measured ON:OFF ratio on IgG was ' + fold.toFixed(1) +
    '× the ratio on purified HaloTag7. The pKa-only model gives all three conjugates the same prediction. Both the attachment chemistry and protein environment can contribute.';
  $('conjugate-source').href = source;
}

function showTitrations(analysis, source) {
  const conditions = analysis.byPower, scored = conditions.filter(c => c.mae);
  const wins = scored.filter(c => c.mae.bound.local < c.mae.free.local).length;
  $('titrations').replaceChildren(...conditions.map(c => row([
    String(c.power), c.scored + '/' + c.total, f3(c.mae?.free.local), f3(c.mae?.bound.local),
    f3(c.mae?.free.nearest), f3(c.mae?.bound.nearest)])));
  $('titration-effect').textContent = 'Replacing free-parent pKa with HaloTag-bound pKa lowered local-model error at ' +
    wins + ' of ' + scored.length + ' laser powers. The same two held-out dyes, SRhB and STMR, are scored at every power.';
  $('titration-source').href = source;
}

export async function renderAnalogueBenchmark(root) {
  const [result, meta, challenge, ratios, evidence] = await Promise.all([load('analogue-benchmark'), load('analogue-studies'),
    load('h2b-cal664sb-lock'), load('ratio-sensitivity'), load('ratio-evidence')]);
  root.innerHTML = `
    <article class="card unseen-headline"><div><div class="eyebrow">PRIMARY ASSAY · FIXED H2B–HALOTAG</div>
      <h2 id="unseen-gate"></h2><p id="unseen-counts"></p></div>
      <a class="button secondary" href="https://github.com/kangk1204/blink-atlas/blob/main/analysis/analogue-benchmark.md">Protocol &amp; development plan ↗</a></article>
    <article class="card unseen-audit"><div class="eyebrow">REGISTERED FORECAST · AWAITING MEASUREMENT</div>
      <h2>Cal664sb on H2B–HaloTag</h2>
      <p id="unseen-challenge-input"></p>
      <div class="table-wrap"><table><thead><tr><th>Local interpolation</th><th>Nearest pKa</th><th>Original regression</th></tr></thead>
      <tbody id="unseen-challenge-predictions"></tbody></table></div><p id="unseen-challenge-boundary"></p>
      <div class="unseen-links"><a id="unseen-registration">Dated registration ↗</a><a href="data/h2b-cal664sb-lock.json" download>Download locked forecasts</a>
      <a id="unseen-challenge-figure">Source titration ↗</a><a href="https://github.com/kangk1204/blink-atlas/tree/main/validation/h2b_blind">Run an independent panel ↗</a></div></article>
    <article class="card unseen-audit" id="unseen-ratio-audit"></article>
    <div class="card unseen-selectors"><label>Test specimen<select id="unseen-assay"></select></label>
      <label>Training data<select id="unseen-training"></select></label><p id="unseen-scope"></p></div>
    <div class="unseen-grid">
      <article class="plot-card unseen-chart-card"><div class="plot-label"><i></i> HELD-OUT DYES</div><h2>Which prediction was closer?</h2>
        <p class="unseen-chart-sub">Absolute log₁₀ error of nearest pKa minus local interpolation</p>
        <div id="unseen-chart" role="group" aria-label="Paired error differences by dye"></div>
        <div id="unseen-axis"><div><span>← Nearest pKa closer</span><span>Interpolation closer →</span></div><small id="unseen-limits"></small></div>
      </article>
      <article class="card unseen-result" aria-live="polite"><p class="eyebrow" id="unseen-coverage"></p><h2>Prediction error</h2>
        <div id="unseen-metrics"></div><p class="footnote">Mean absolute log₁₀ error on the same dyes. Lower is better.</p>
        <p id="unseen-contrast"></p><p id="unseen-uncertainty"></p><p id="unseen-original"></p></article>
    </div>
    <article class="card unseen-audit"><h2>Predictions and observations</h2><p id="unseen-measure"></p>
      <div class="table-wrap"><table><thead><tr><th>Test dye</th><th>pKa</th><th>Observed</th><th>Local</th><th>Nearest</th><th>Regression</th><th>Error gain</th><th>Training</th></tr></thead><tbody id="unseen-rows"></tbody></table></div>
      <p id="unseen-assay-info"></p><p id="unseen-source"></p></article>
    <article class="card unseen-audit"><h2>One dye core, three binding partners</h2>
      <p>Remmel and colleagues measured three conjugates of compound 4. All share the free parent’s pKa₂.</p>
      <div class="table-wrap"><table><thead><tr><th>Binding partner</th><th>Parent pKa₂</th><th>Measured ON:OFF</th><th>Frozen local prediction</th></tr></thead><tbody id="unseen-conjugates"></tbody></table></div>
      <p id="unseen-conjugate-effect"></p><a id="unseen-conjugate-source" class="inline-link">Remmel et al., Supplementary Table S1 ↗</a></article>
    <article class="card unseen-audit"><div class="eyebrow">EXPLORATORY · PURIFIED HALOTAG</div>
      <h2>Does measuring the bound dye help?</h2>
      <p>Zheng and colleagues measured both free and HaloTag-bound titrations. We changed the pKa input and repeated the same within-study holdouts. Errors below are for occupancy odds over a 60 s recording.</p>
      <div class="table-wrap"><table><thead><tr><th>Power (W/cm²)</th><th>Dyes scored</th><th>Local · free</th><th>Local · bound</th><th>Nearest · free</th><th>Nearest · bound</th></tr></thead><tbody id="unseen-titrations"></tbody></table></div>
      <p id="unseen-titration-effect"></p><a id="unseen-titration-source" class="inline-link">Zheng et al., supporting tables S1 and S2 ↗</a></article>
    <article class="card unseen-audit" id="unseen-ratio-evidence"></article>
    <article class="card unseen-audit"><h2>Evidence and exclusions</h2><p>The candidate interpolates the log ratios of two training dyes that bracket the measured pKa. The comparator copies the closest training dye. Both abstain outside the calibration range.</p>
      <p id="unseen-screen-date"></p><details><summary>Inspect screened sources and exclusions</summary><ul id="unseen-screen"></ul></details>
      <div class="unseen-links"><a href="data/analogue-observations.csv" download>Download source measurements (CSV)</a>
        <a href="data/analogue-benchmark.json" download>Download all predictions (JSON)</a><a id="unseen-freeze">Frozen evaluation protocol ↗</a>
        <a href="#probe-design">Lock a new H2B prediction →</a></div></article>`;
  renderRatioAudit($('ratio-audit'), ratios, evidence);
  renderRatioEvidence($('ratio-evidence'), evidence);
  const primary = result.sets.find(s => s.mode === 'primary');
  const candidate = challenge.dyes[0];
  $('challenge-input').textContent = 'Tyson et al. measured a free-parent pKcycl of ' + candidate.pKa +
    ' in sodium phosphate. These conditional mean ON:OFF predictions use that value. A new titration in the Holland buffer and a matched cellular acquisition are still needed for the primary test.';
  $('challenge-predictions').replaceChildren(row(Object.keys(labels).map(k => small(candidate.predictions[k]))));
  $('challenge-boundary').textContent = 'At this input, local interpolation beats both baselines if the observed mean ON:OFF exceeds ' +
    small(challenge.equalErrorBoundaries[1].ratio) + '. A single dye tests this forecast for that dye.';
  $('challenge-figure').href = challenge.source.figure;
  $('registration').href = 'https://github.com/kangk1204/blink-atlas/blob/' + challenge.registrationRevision + '/validation/h2b_blind/registration.md';
  $('gate').textContent = primary.summary.scored ? result.primaryVerdict : 'Independent H2B validation is still open';
  $('counts').textContent = result.uniqueExternalDyes + ' additional parent dyes, ' + result.externalConditions + ' dye–condition measurements from three public studies. ' +
    primary.summary.scored + ' new dyes scored under the matched H2B protocol.';
  const groups = new Map();
  for (const set of result.sets) {
    const group = set.study === 'holland' ? 'holland' : set.id.replace(/-(frozen|local)$/, '');
    if (!groups.has(group)) groups.set(group, []);
    groups.get(group).push(set);
  }
  for (const [key, sets] of groups) {
    const title = key === 'holland' ? 'Holland · fixed H2B' : sets[0].label.split(' · ').slice(0,-1).join(' · ');
    $('assay').add(new Option(title, key));
  }
  function selectAssay() {
    const options = groups.get($('assay').value);
    const names = {primary:'Independent matched dyes · frozen fit',development:'Development dyes · leave one out',
      'transfer-frozen':'Frozen Holland fit','transfer-local':'Other dyes in this assay'};
    $('training').replaceChildren(...options.map(s => new Option(names[s.mode], s.id)));
    show(options[0], meta, result.matchedSource);
  }
  $('assay').value = 'remmel-igg'; $('assay').onchange = selectAssay;
  $('training').onchange = () => show(result.sets.find(s => s.id === $('training').value), meta, result.matchedSource);
  selectAssay();
  showConjugates(result.sets, meta.studies.remmel.conjugateSourceUrl);
  showTitrations(result.exploratoryTitrations, meta.studies.zheng.sourceUrl);
  $('screen-date').textContent = 'Literature and repository screening date ' + meta.searchDate +
    '. A different specimen, an already-known dye or an incompatible response cannot confirm the primary target.';
  for (const s of meta.screened) {
    const li = el('li'); li.append(link(s.label, s.url), document.createTextNode(' — ' + s.reason)); $('screen').append(li);
  }
  $('freeze').href = 'https://github.com/kangk1204/blink-atlas/blob/' + result.frozenRevision + '/analysis/analogue-benchmark.md';
}
