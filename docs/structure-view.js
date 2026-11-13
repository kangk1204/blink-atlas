const dataDir = 'data/structure-pk/';
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

function candidate(c, i) {
  const box = node('section', undefined, 'pk-candidate');
  const title = i === 0 ? 'Fluoromethyl · CH₂F' : 'Difluoromethyl · CHF₂';
  const graphic = node('img');
  graphic.src = dataDir + c.dye + '.svg';
  graphic.alt = 'Neutral closed Si-rhodamine with ' + c.substitution;
  graphic.width = 520; graphic.height = 265;
  const metrics = node('dl', undefined, 'pk-estimates');
  for (const k of ['nearest', 'ridge']) {
    const pair = node('div');
    pair.append(node('dt', names[k]), node('dd', fixed(c[k])));
    metrics.append(pair);
  }
  box.append(node('h3', title), graphic, metrics,
    node('p', 'Closest measured dye: ' + c.neighbors.map(n => n.dye + ' (' + fixed(n.measured) + ')').join(', ') +
      '. Tanimoto ' + c.similarity.toFixed(3) + '. The estimates differ by ' + fixed(c.modelGap) + ' pH units.'),
    links([link('SDF ↓', dataDir + c.dye + '.sdf', true), link('Prediction JSON ↓', dataDir + c.dye + '.json', true),
      link('2D structure ↓', graphic.src, true), link('Azetidine building block ↗', 'https://pubchem.ncbi.nlm.nih.gov/compound/' + [66596287,55299483][i]),
      link('Matched synthesis precedent ↗', 'https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM1_ESM.pdf#page=44')]));
  return box;
}

function benchmark(b) {
  const view = detail('How well do the models predict an unseen study?');
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

function importer(b) {
  const block = detail('Open a result from your own structure');
  block.id = 'pk-import';
  const label = node('label', 'Local prediction JSON ');
  const input = node('input'); input.type = 'file'; input.accept = '.json,application/json'; input.id = 'pk-file';
  label.append(input);
  const status = node('p'); status.id = 'pk-file-status'; status.setAttribute('role', 'status');
  block.append(node('p', 'Run the local CLI with a closed free-parent SMILES or single-record SDF. Open its JSON here to check both predictions against the released model. The file stays in your browser.'),
    node('code', 'python -m analysis.structure_pk predict --sdf dye.sdf --output result.json'),
    links([link('CLI setup and SMILES example ↗', repo + 'analysis/structure-pk.md'),
      link('Sample JSON ↓', dataDir + 'HM-SiR-CH2F.json', true)]), label, status);
  let generation = 0;
  input.onchange = async () => {
    const run = ++generation, file = input.files[0];
    status.replaceChildren();
    if (!file) return;
    try {
      if (file.size > 100000) throw Error('The file exceeds 100 kB.');
      const r = JSON.parse(await file.text());
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
    const res = await fetch(dataDir + 'benchmark.json');
    if (!res.ok) throw Error('HTTP ' + res.status);
    const b = await res.json(), [a, c, ctrl] = b.candidates;
    const grid = node('div', undefined, 'pk-grid');
    grid.append(candidate(a, 0), candidate(c, 1));
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
    const info = detail('Source structures, measurements and exclusions');
    info.append(node('p', b.uniqueParents + ' unique parents are represented, including an open noncyclizing control. ' +
      b.trainingParents + ' have eligible original absorbance measurements. Separate fields preserve secondary quotations, bounds, complex curves and alternate fluorescence readouts.'),
      table(['Study', 'Parents represented', 'Eligible original parents'], Object.entries(b.counts).map(([s, x]) => [s, x.represented, x.exactOriginal])),
      links([link('All measurements and exclusions ↓', dataDir + 'measurements.csv', true), link('Parent structures ↓', dataDir + 'parents.sdf', true),
        link('Source files and hashes ↗', repo + 'tools/structure_sources.py'), link('Frozen design ↗', repo + 'validation/structure_pk/registration.json')]));
    const title = node('h2', 'Which free dye should we titrate?'); title.id = 'pk-title';
    root.replaceChildren(title, node('p', 'Two virtual Si-rhodamines differ by one fluorine on each azetidine arm. Their free-dye pKcycl values await measurement.'),
      grid, parts, node('p', 'The default nearest-structure model gives the candidates the same midpoint. Ridge separates them by ' +
        Math.abs(a.ridge - c.ridge).toFixed(2) + ' pH units. Their order is unresolved. H2B registration requires a measured free-dye pKa.'),
      benchmark(b), assay, importer(b), info);
  } catch (e) {root.append(node('p', 'The structure benchmark could not load. ' + e.message));}
}
