// Prospective predictions and measurement scoring, shared by the UI and round-trip tests.
const LEGACY = 'blink-atlas/probe-trial/1';
const SCHEMA = 'blink-atlas/probe-trial/2';
const RULE = 'Three independent preparations. Average the three mean ON:OFF time ratios arithmetically. Compare absolute log10 errors to the frozen predictions. Lower is better.';
export const methods = {regression: 'pKa regression', nearest: 'Nearest measured pKa', median: 'Training median',
  local: 'Local log interpolation'};

function measuredNumber(input, label) {
  if (typeof input === 'number' && Number.isFinite(input)) return input;
  if (typeof input === 'string' && /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(input.trim())) {
    const number = Number(input);
    if (Number.isFinite(number)) return number;
  }
  throw new Error('Enter a finite numeric ' + label + '.');
}

export function predict(calibration, pKaInput, scope, schema = SCHEMA) {
  if (scope !== calibration.scope) {
    throw new Error('A prediction needs a monotonic free-parent titration from this scaffold and the matched fixed COS-7 H2B–HaloTag assay.');
  }
  const pKa = measuredNumber(pKaInput, 'pKa');
  const fit = calibration.fit;
  const [low, high] = fit.pKaRange;
  if (pKa < low || pKa > high) throw new Error('Outside calibration. The measured pKa range is ' + low + ' to ' + high + '.');
  const sorted = [...calibration.pairs].sort((a, b) => a.pKa - b.pKa);
  const byDistance = row => Math.round(Math.abs(row.pKa - pKa) * 1e12);
  const nearest = sorted.reduce((a, b) => byDistance(a) <= byDistance(b) ? a : b);
  const lower = sorted.filter(r => r.pKa <= pKa).at(-1);
  const upper = sorted.find(r => r.pKa >= pKa);
  const forecast = {pKa, references: [...new Set([lower.id, upper.id])], nearestDye: nearest.id,
    ratios: {regression: 10 ** (fit.intercept + fit.slope * pKa), nearest: nearest.ratio, median: 10 ** fit.medianLog}};
  if (schema === LEGACY) return forecast;
  const at = x => sorted.filter(r => r.pKa === x);
  const averageLog = rows => rows.reduce((sum, r) => sum + Math.log10(r.ratio), 0) / rows.length;
  const lowerRows = at(lower.pKa), upperRows = at(upper.pKa);
  const weight = lower.pKa === upper.pKa ? 0 : (pKa - lower.pKa) / (upper.pKa - lower.pKa);
  forecast.references = [...new Set([...lowerRows, ...upperRows].map(r => r.id))];
  if (at(nearest.pKa).length > 1) {
    forecast.nearestDye = at(nearest.pKa).map(r => r.id).join(', ');
    forecast.ratios.nearest = 10 ** averageLog(at(nearest.pKa));
  }
  forecast.ratios.local = 10 ** ((1 - weight) * averageLog(lowerRows) + weight * averageLog(upperRows));
  return forecast;
}

function makePlan(data, labelInput, pKa, scope, createdUTC, schema = SCHEMA) {
  const label = typeof labelInput === 'string' ? labelInput.trim() : '';
  if (!label || label.length > 80) throw new Error('Give the new compound a trial ID of 1–80 characters.');
  const key = v => v.toLowerCase().replace(/[\s–_-]/g, '');
  if (data.dyes.some(d => [d.id, d.ligand].some(name => key(name) === key(label)))) {
    throw new Error('This dye already has published measurements. Use a new analogue ID for a prospective trial.');
  }
  const date = new Date(createdUTC);
  if (!Number.isFinite(date.getTime()) || date.toISOString() !== createdUTC) throw new Error('Invalid card timestamp.');
  const c = data.calibration;
  const forecast = predict(c, pKa, scope, schema);
  const version = schema === LEGACY ? {} : {candidate: 'bracket-logratio-v1'};
  return {schema, model: c.id, ...version, sourceSha256: c.sourceSha256, article: data.article,
    createdUTC, compound: label, scope, target: c.target, ...forecast, scoringRule: RULE};
}

export function lockTrial(data, name, pKa, scope, createdUTC = new Date().toISOString()) {
  return {plan: makePlan(data, name, pKa, scope, createdUTC), measurement: null};
}

export function scoreTrial(plan, entries) {
  if (!Array.isArray(entries) || entries.length !== 3) throw new Error('Enter one ON:OFF ratio from each of three independent preparations.');
  const ratios = entries.map(value => measuredNumber(value, 'ratio for each preparation'));
  if (ratios.some(r => r <= 0)) throw new Error('Every measured time ratio must be positive.');
  const max = Math.max(...ratios);
  const observed = max * (ratios.reduce((sum, r) => sum + r / max, 0) / ratios.length);
  const error = {};
  for (const [method, predicted] of Object.entries(plan.ratios)) {
    error[method] = Math.abs(Math.log10(observed) - Math.log10(predicted));
  }
  const best = Math.min(...Object.values(error));
  return {ratios, observed, error, closest: Object.keys(error).filter(m => Math.abs(error[m] - best) <= 1e-12)};
}

function sameSnapshot(saved, expected) {
  if (typeof expected === 'number') {
    return Number.isFinite(saved) && Math.abs(saved - expected) <= 1e-12 * Math.max(Number.MIN_VALUE, Math.abs(expected));
  }
  if (expected === null || typeof expected !== 'object') return saved === expected;
  if (saved === null || typeof saved !== 'object' || Array.isArray(saved) !== Array.isArray(expected)) return false;
  const keys = Object.keys(expected);
  return keys.length === Object.keys(saved).length && keys.every(k => Object.hasOwn(saved, k) && sameSnapshot(saved[k], expected[k]));
}

export function restoreTrial(data, record) {
  if (![SCHEMA, LEGACY].includes(record?.plan?.schema) || !Object.hasOwn(record, 'measurement')) throw new Error('Choose a Blink Atlas probe trial JSON.');
  const p = record.plan;
  const trial = {plan: makePlan(data, p.compound, p.pKa, p.scope, p.createdUTC, p.schema), measurement: null};
  if (!sameSnapshot(p, trial.plan)) throw new Error('The saved plan disagrees with the source, model or predictions for this release.');
  if (record.measurement !== null) {
    trial.measurement = scoreTrial(trial.plan, record.measurement?.ratios);
    if (!sameSnapshot(record.measurement, trial.measurement)) throw new Error('The saved score disagrees with its three measurements.');
  }
  return trial;
}
