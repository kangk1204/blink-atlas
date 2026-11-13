const dataDir = 'data/structure-pk/';
const transferDir = 'data/structure-transfer/';
const repo = 'https://github.com/kangk1204/blink-atlas/blob/main/';
const fixed = value => value.toFixed(2);
const names = {nearest: 'Nearest structure', ridge: 'Kernel ridge'};
const near = (x, y) => Number.isFinite(x) && Math.abs(x - y) < 1e-7;

function node(tag, value, className) {
  const el = document.createElement(tag);
  if (value !== undefined) el.textContent = value;
  if (className) el.className = className;
  return el;
}

function link(label, href, download = false) {
  const a = node('a', label);
  a.href = href;
  if (download) a.download = '';
  return a;
}

function links(items) {
  const group = node('div', undefined, 'pk-links');
  group.append(...items);
  return group;
}

function detail(title) {
  const el = node('details');
  el.append(node('summary', title));
  return el;
}

function table(headings, values) {
  const wrap = node('div', undefined, 'table-wrap');
  const t = node('table');
  const header = t.createTHead().insertRow();
  headings.forEach(h => {const th = node('th', h); th.scope = 'col'; header.append(th);});
  const body = t.createTBody();
  values.forEach(row => {
    const tr = body.insertRow();
    row.forEach(cell => tr.insertCell().textContent = cell);
  });
  wrap.append(t);
  return wrap;
}

function candidate(c, i, updated) {
  const box = node('section', undefined, 'pk-candidate');
  const title = i === 0 ? 'Fluoromethyl · CH₂F' : 'Difluoromethyl · CHF₂';
  const graphic = node('img');
  graphic.src = dataDir + c.dye + '.svg';
  graphic.alt = 'Neutral closed Si-rhodamine with ' + c.substitution;
  graphic.width = 520; graphic.height = 265;
  const metrics = node('dl', undefined, 'pk-estimates');
  for (const k of ['nearest', 'ridge']) {
    const pair = node('div');
    pair.append(node('dt', names[k] + (k === 'ridge' ? ' (original)' : '')), node('dd', fixed(c[k])));
    metrics.append(pair);
  }
  const next = node('div');
  next.append(node('dt', 'Core and donor regression'), node('dd', fixed(updated.result.predictions.coreDonor)));
  metrics.append(next);
  box.append(node('h3', title), graphic, metrics,
    node('p', 'Closest measured dye: ' + c.neighbors.map(n => n.dye + ' (' + fixed(n.measured) + ')').join(', ') +
      '. Tanimoto ' + c.similarity.toFixed(3) + '.'),
    links([link('SDF ↓', dataDir + c.dye + '.sdf', true), link('Original forecast ↓', dataDir + c.dye + '.json', true),
      link('2D structure ↓', graphic.src, true), link('Azetidine building block ↗', 'https://pubchem.ncbi.nlm.nih.gov/compound/' + [66596287,55299483][i]),
      link('Matched synthesis precedent ↗', 'https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM1_ESM.pdf#page=44'),
      link('Updated predictions ↓', transferDir + updated.dye + '.json', true)]));
  return box;
}

function benchmark(b) {
  const view = detail('Original 26-parent benchmark and quantum comparison');
  view.id = 'pk-evaluation';
  view.append(node('p', 'Retrospective mean absolute error, in pH units. Each test study and all its parent structures were removed from training. Only exact, single-transition absorbance midpoints were scored. Each dye contributes once.'),
    table(['Held-out study', 'Scored / parents', 'Training dyes', 'Nearest MAE', 'Ridge MAE', 'Shared parents'],
      b.folds.map(f => [f.study, f.eligible + ' / ' + f.total, f.train.length, f.mae.nearest.toFixed(3),
        f.mae.ridge.toFixed(3), f.overlap])));
  for (const f of b.folds) {
    const d = detail(f.study + ' · individual errors and sensitivity');
    d.append(table(['Free dye', 'Measured', 'Nearest', 'Abs. error', 'Ridge', 'Abs. error', 'Nearest training dye(s)', 'Tanimoto'],
      f.rows.map(r => [r.dye, fixed(r.measured), fixed(r.nearest), fixed(r.nearestError), fixed(r.ridge),
        fixed(r.ridgeError), r.neighbors.map(n => n.dye).join(', '), r.similarity.toFixed(3)])),
      node('p', 'MAE after removing one test dye from the study average. Training and predictions stay fixed.'),
      table(['Omitted dye', 'Nearest MAE', 'Ridge MAE'], f.leaveOneOut.map(r => [r.omit, r.nearest.toFixed(3), r.ridge.toFixed(3)])));
    const reversals = f.leaveOneOut.filter(r => (r.ridge < r.nearest) !== (f.mae.ridge < f.mae.nearest));
    if (reversals.length) view.append(node('p', f.study + ' changes the deployment decision. Omitting ' +
      reversals.map(r => r.omit).join(', ') + ' reverses the MAE ranking. Both studies had to favor ridge after every single-dye omission.'));
    view.append(d);
  }
  const match = b.quantumSummary.primaryMatched;
  const qm = detail('Compare the published quantum calculations');
  qm.append(node('p', 'On ' + match.n + ' Table 1 parents with a matching original single-transition measurement, mean absolute errors are ' +
    match.nearest.toFixed(3) + ' (nearest), ' + match.ridge.toFixed(3) + ' (ridge) and ' + match.quantum.toFixed(3) +
    ' (published quantum calculation). Each structure model excludes the entire test parent. The quantum protocol was developed on some of these parents.'),
    node('p', 'All 17 original table rows are retained. HMRG and AMRG use the authors’ revised five-water results. The download includes the initial three-water values and identifies the development rows.'),
    links([link('17-row comparison ↓', dataDir + 'quantum.csv', true), link('Tachibana Table 1 ↗', 'https://www.nature.com/articles/s42004-020-0326-x/tables/1')]));
  view.append(qm, links([link('Per-dye errors ↓', dataDir + 'holdouts.csv', true),
    link('Folds and models ↓', dataDir + 'benchmark.json', true)]));
  return view;
}

// Verify numerical predictions against the released model using the CLI-supplied fingerprint.
// RDKit assigns the chemical identity on the user's machine.
function verifyImport(r, b) {
  if (!r || r.schema !== 'blink-atlas/free-parent-pk/1' || r.modelId !== b.modelId ||
      r.endpoint !== 'Free-parent pKcycl' || r.preferred !== b.preferred || r.trainingParents !== b.trainingParents ||
      typeof r.dye !== 'string' || !r.dye.length || r.dye.length > 256 ||
      !/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(r.inchiKey) || typeof r.smiles !== 'string' ||
      !r.smiles.length || r.smiles.length > 12000) throw Error('This file does not match the released structure model.');
  const bits = r.fingerprint;
  if (!Array.isArray(bits) || !bits.length || bits.length > 2048 ||
      bits.some((n, i) => !Number.isInteger(n) || n < 0 || n >= 2048 || (i > 0 && n <= bits[i-1]))) throw Error('Invalid molecular fingerprint.');
  const v = b.verification, set = new Set(bits);
  const sim = v.parents.map(p => {
    const shared = p.bits.filter(n => set.has(n)).length;
    return shared / (set.size + p.bits.length - shared);
  });
  const max = Math.max(...sim), closest = v.parents.filter((_, i) => near(sim[i], max));
  const ref = v.parents.find(p => p.inchiKey === r.inchiKey);
  const obs = v.observed[r.inchiKey] || [];
  const status = ref ? 'known exact' : obs.length ? 'observed complex or limited' : 'unmeasured in corpus';
  const publicStructure = b.candidates.find(c => c.inchiKey === r.inchiKey);
  if (r.status !== status || (ref && !near(sim[v.parents.indexOf(ref)], 1)) ||
      (publicStructure && (r.smiles !== publicStructure.smiles || String(bits) !== String(publicStructure.fingerprint))) ||
      !near(r.similarity, max) || !Array.isArray(r.neighbors) || r.neighbors.length !== closest.length ||
      !closest.every((p, i) => r.neighbors[i]?.inchiKey === p.inchiKey)) throw Error('Structure identity or nearest-neighbor check failed.');
  const limited = status === 'observed complex or limited';
  const nearest = closest.reduce((s, p) => s + p.measured, 0) / closest.length;
  const avg = sim.reduce((a, n) => a + n, 0) / sim.length;
  const ridge = sim.reduce((s, n, i) => s + (n - avg - v.columnMeans[i] + v.grandMean) * v.weights[i], v.mean);
  if ((ref ? !near(r.knownMeasured, ref.measured) : r.knownMeasured !== null) ||
      (limited ? [r.nearest, r.ridge, r.modelGap].some(n => n !== null) :
        (!near(r.nearest, nearest) || !near(r.ridge, ridge) || !near(r.modelGap, Math.abs(ridge - nearest))))) {
    throw Error('The numerical results disagree with the released model.');
  }
  return {status, obs, closest, limited};
}

function transferSection(b) {
  const block = detail('Why did the original model miss the Holland dyes?');
  block.id = 'pk-transfer';
  const hol = b.folds[0], audit = b.audit;
  block.append(node('p', 'The expanded reference set contains ' + b.trainingParents + ' unique measured parents. All six held-out Holland dyes have a silicon bridge. The original nearest neighbors were the oxygen-bridged HMTMR and HMRB. None of the ' +
    audit.nTrain + ' original training dyes contained an azetidine. ' + audit.belowTrain + ' test midpoints fell below the training minimum of ' + fixed(audit.trainRange[0]) + '.'),
    node('p', 'Adding the measured azetidinyl Si-rhodamine HM-DS655 changes all six nearest neighbors. Their common estimate is 5.30. Holland MAE falls from ' +
      audit.priorMAE.nearest.toFixed(3) + ' to ' + hol.summary.metrics.nearest.mae.toFixed(3) +
      ' pH units. Removing HM-DS655 from the expanded training set gives ' + b.hollandWithoutChi.summary.metrics.nearest.mae.toFixed(3) + ' again. The six-way tie leaves their order unresolved.'),
    node('h3', 'Unseen-study performance'),
    node('p', 'MAE in pH units for the same dyes and training parents. The new ridge regression uses the bridge, spiro ring and amino-donor environment. Kernel ridge uses the whole-molecule fingerprint. Each held-out parent is removed across studies.'),
    table(['Held-out panel', 'Scored', 'Nearest', 'Fingerprint ridge', 'Core and donor', 'Training mean'],
      [...b.folds, b.newParents].map(f => [f.study === 'New parents' ? 'New parents · three studies excluded together' : f.study,
        f.rows.length, ...['nearest','kernel','coreDonor','trainMean'].map(k => f.summary ? f.summary.metrics[k].mae.toFixed(3) : '—')])));
  const n = b.newParents;
  block.append(node('p', 'The combined check holds out Sakabe, Uno 2018 and Chi 2020 together. It scores the three new eligible identities, from Sakabe and Chi. They overlap the other analyses. The training-mean baseline also outperforms the new regression on this small panel.'),
    table(['New parent', 'Measured', 'Nearest', 'Fingerprint ridge', 'Core and donor'],
      n.rows.map(r => [r.dye, fixed(r.measured), ...['nearest','kernel','coreDonor'].map(k => fixed(r.predictions[k]))])),
    node('p', 'Nearest structure remains the default. The new regression misses the whole-study target on both Holland and Tyson. With Holland siblings available, holding out just one Holland dye gives MAEs of ' +
      b.withinHolland.summary.metrics.coreDonor.mae.toFixed(3) + ' for the new regression and ' +
      b.withinHolland.summary.metrics.nearest.mae.toFixed(3) + ' for nearest. The ranking changes when individual test dyes are omitted.'),
    node('p', 'Uno reports a midpoint of 7.2 for HMbisTFER and both optical channels. We could not verify which channel supplied that fit. Allowing it in a separate sensitivity check changes the combined panel to four parents. HM-DS531 is permanently charged and outside the neutral-parent domain.'),
    links([link('Per-dye errors ↓', transferDir + 'per-dye.csv', true), link('Added assays ↓', transferDir + 'assays.csv', true),
      link('Folds, sensitivity and models ↓', transferDir + 'benchmark.json', true),
      link('Methods and reproduction ↗', repo + 'analysis/structure-transfer.md')]));
  return block;
}

function checkTransfer(r, b) {
  if (r.modelId !== b.modelId || r.preferred !== b.preferred || typeof r.dye !== 'string' ||
      !r.dye.length || r.dye.length > 256 || !/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(r.inchiKey) ||
      typeof r.smiles !== 'string' || !r.smiles.length) throw Error('This file does not match the released transfer model.');
  const fp = r.fingerprint;
  if (!Array.isArray(fp) || !fp.length || fp.length > 2048 ||
      fp.some((n,i) => !Number.isInteger(n) || n < 0 || n >= 2048 || (i && n <= fp[i-1]))) throw Error('Invalid molecular fingerprint.');
  const set = new Set(fp), v = b.verification;
  const sims = v.bits.map(bits => {
    const same = bits.filter(i => set.has(i)).length;
    return same / (bits.length + fp.length - same);
  });
  const best = Math.max(...sims), avg = sims.reduce((a,n) => a+n, 0) / sims.length;
  const ties = b.training.filter((p,i) => near(sims[i], best));
  const known = b.training.findIndex(p => p.inchiKey === r.inchiKey), obs = b.observed[r.inchiKey] || [];
  const status = known >= 0 ? 'known exact' : obs.length ? 'observed complex or limited' : 'unmeasured in corpus';
  if (r.status !== status || (known >= 0 ? (!near(r.knownMeasured,b.training[known].measured) || !near(sims[known],1)) :
      r.knownMeasured !== null)) throw Error('Parent identity or measurement check failed.');
  if (status === 'observed complex or limited') {
    if (r.result !== null) throw Error('A single-midpoint result is withheld for this parent.');
    return r.dye + ' · ' + obs.map(x => x.status + ' (' + x.reported + ')').join(', ') + '. A single midpoint is withheld.';
  }
  const z = r.result;
  if (!z || !z.predictions || !z.descriptors || !near(z.similarity, best) || !Array.isArray(z.neighbors) ||
      ties.length !== z.neighbors.length || !ties.every((p,i) => z.neighbors[i]?.inchiKey === p.inchiKey)) {
    throw Error('Nearest structure check failed.');
  }
  const standard = b.candidates.find(c => c.inchiKey === r.inchiKey);
  if (standard && (r.smiles !== standard.smiles || String(fp) !== String(standard.fingerprint) ||
      b.fullModel.names.some(n => !near(z.descriptors[n],standard.result.descriptors[n])))) throw Error('Candidate structure check failed.');
  const linear = m => m.names.reduce((acc,n,i) => acc + (z.descriptors[n] - m.center[i]) * m.weights[i] / m.scale[i], m.intercept);
  const predicted = {
    nearest: ties.reduce((a,p) => a+p.measured,0) / ties.length,
    kernel: sims.reduce((a,s,i) => a+(s-avg-v.means[i]+v.grandMean)*v.weights[i], v.intercept),
    coreDonor: linear(b.fullModel), withoutWithdrawal: linear(v.ablated), trainMean: v.intercept,
  };
  if (!b.fullModel.names.every(n => Number.isFinite(z.descriptors[n])) ||
      !Object.keys(predicted).every(k => near(z.predictions[k], predicted[k]))) throw Error('The transfer forecasts disagree with the released model.');
  return r.dye + ' · Nearest ' + fixed(predicted.nearest) + ' · Fingerprint ridge ' + fixed(predicted.kernel) +
    ' · Core and donor ' + fixed(predicted.coreDonor) + ' · Closest measured: ' + ties.map(p => p.dye).join(', ') +
    (known < 0 ? '. No eligible measurement in this corpus.' : '. Known training parent, measured ' + fixed(r.knownMeasured) + '.');
}

function importer(b, transfer) {
  const block = detail('Open a result from your own structure');
  block.id = 'pk-import';
  const label = node('label', 'Local prediction JSON ');
  const input = node('input'); input.type = 'file'; input.accept = '.json,application/json'; input.id = 'pk-file';
  label.append(input);
  const status = node('p'); status.id = 'pk-file-status'; status.setAttribute('role', 'status');
  block.append(node('p', 'Run the local CLI with a closed free-parent SMILES or single-record SDF. Open its JSON here to check the forecasts against the published weights. The file stays in your browser.'),
    node('code', 'python -m analysis.structure_transfer predict --sdf dye.sdf --output result.json'),
    links([link('CLI setup and SMILES example ↗', repo + 'analysis/structure-transfer.md'),
      link('Sample JSON ↓', transferDir + 'HM-SiR-CH2F.json', true)]), label, status);
  let generation = 0;
  input.onchange = async () => {
    const run = ++generation, file = input.files[0];
    status.replaceChildren();
    if (!file) return;
    try {
      if (file.size > 100000) throw Error('The file exceeds 100 kB.');
      const r = JSON.parse(await file.text());
      if (r?.schema === 'blink-atlas/transfer-query/1') {
        const message = checkTransfer(r, transfer);
        if (run === generation) status.textContent = message;
        return;
      }
      const checked = verifyImport(r, b);
      if (run !== generation) return;
      const msg = checked.limited ? 'Published curve or limit: ' + checked.obs.map(x => x.dye + ' ' + x.reported + ' (' + x.status + ')').join(', ') +
        '. A single midpoint is withheld.' : 'Nearest ' + fixed(r.nearest) + ' · Ridge ' + fixed(r.ridge) +
        ' · Closest measured: ' + checked.closest.map(x => x.dye).join(', ') +
        (r.knownMeasured === null ? '. No eligible measurement in this corpus.' : '. Known training parent, measured ' + fixed(r.knownMeasured) + '.');
      status.textContent = r.dye + ' · ' + msg;
    } catch (e) {if (run === generation) status.textContent = 'Could not open file. ' + e.message;}
  };
  return block;
}

export async function mountStructure(root) {
  try {
    const [b, transfer] = await Promise.all([dataDir, transferDir].map(async d => {
      const res = await fetch(d + 'benchmark.json');
      if (!res.ok) throw Error('HTTP ' + res.status);
      return res.json();
    }));
    const [a, c, ctrl] = b.candidates;
    const grid = node('div', undefined, 'pk-grid');
    grid.append(candidate(a, 0, transfer.candidates[0]), candidate(c, 1, transfer.candidates[1]));
    const parts = node('div', undefined, 'pk-control');
    parts.append(node('strong', 'Matched control · JF626b (CF₃)'),
      node('p', 'Holland measured ' + fixed(ctrl.knownMeasured) + ' for this free dye. Its value is included in the full training set.'),
      links([link('Control SDF ↓', dataDir + 'JF626b.sdf', true), link('Control JSON ↓', dataDir + 'JF626b.json', true),
        link('Structure ↗', dataDir + 'JF626b.svg'), link('Measurement and synthesis, compound 8 ↗', 'https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM1_ESM.pdf#page=44')]));
    const assay = detail('Matched titration and synthesis precedent');
    assay.append(node('p', 'Run both free dyes and JF626b together in the Holland buffers on each of three independent preparation days. Measure absorbance across pH and compare within-day differences. Keep the full traces for complex transitions.'),
      links([link('Titration protocol ↗', repo + 'validation/structure_pk/titration.md'),
        link('Ditriflate coupling, Holland SI ↗', 'https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM1_ESM.pdf#page=41'),
        link('Azetidine substitution, Grimm et al. ↗', 'https://doi.org/10.1038/nmeth.4403'),
        link('Full-parent identity searches ↗', repo + 'validation/structure_pk/identity-search.json')]));
    const info = detail('Original source structures, measurements and exclusions');
    info.append(node('p', b.uniqueParents + ' unique parents are represented, including an open noncyclizing control. ' +
      b.trainingParents + ' have eligible original absorbance measurements. Separate fields preserve secondary quotations, bounds, complex curves and alternate fluorescence readouts.'),
      table(['Study', 'Parents represented', 'Eligible original parents'], Object.entries(b.counts).map(([s, x]) => [s, x.represented, x.exactOriginal])),
      links([link('All measurements and exclusions ↓', dataDir + 'measurements.csv', true), link('Parent structures ↓', dataDir + 'parents.sdf', true),
        link('Source files and hashes ↗', repo + 'tools/structure_sources.py'), link('Frozen design ↗', repo + 'validation/structure_pk/registration.json')]));
    const title = node('h2', 'Which free dye should we titrate?'); title.id = 'pk-title';
    root.replaceChildren(title, node('p', 'Two virtual Si-rhodamines differ by one fluorine on each azetidine arm. Their free-dye pKcycl values await measurement.'),
      grid, parts, node('p', 'The default nearest-structure model gives both candidates the same midpoint. The original fingerprint ridge separates them by ' +
        Math.abs(a.ridge - c.ridge).toFixed(2) + ' pH units. The new regression has substantial whole-study error. Their order awaits titration. H2B registration requires a measured free-dye pKa.'),
      transferSection(transfer), benchmark(b), assay, importer(b, transfer), info);
  } catch (e) {root.append(node('p', 'The structure benchmark could not load. ' + e.message));}
}
