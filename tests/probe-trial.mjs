import assert from 'node:assert/strict';
import {test} from 'node:test';
import {readFile} from 'node:fs/promises';
import {lockTrial, predict, restoreTrial, scoreTrial} from '../docs/probe-trial.js';

const data = JSON.parse(await readFile(new URL('../docs/data/probe-design.json', import.meta.url), 'utf8'));
const c = data.calibration;
const benchmark = JSON.parse(await readFile(new URL('../docs/data/analogue-benchmark.json', import.meta.url), 'utf8'));
const create = () => lockTrial(data, 'blind-01', '3.20', c.scope, '2026-11-12T01:02:03.000Z');

test('locked predictions, endpoints and a deterministic nearest reference at a midpoint', () => {
  const p = create().plan;
  assert.equal(p.schema, 'blink-atlas/probe-trial/2');
  assert.equal(p.candidate, benchmark.method);
  assert.equal(p.nearestDye, 'JF635b');
  assert.deepEqual(p.references, ['JF635b','JF639b']);
  assert.ok(Math.abs(p.ratios.regression - .0004297281986991795) < 1e-14);
  assert.equal(p.ratios.nearest, .000483915);
  assert.equal(predict(c, 3.46, c.scope).nearestDye, 'JF635b');
  assert.deepEqual(predict(c, 1.22, c.scope).references, ['JF626b']);
  assert.deepEqual(predict(c, 4.48, c.scope).references, ['JFX650b']);
});

test('browser forecasts match the CPU benchmark on every eligible frozen-test row', () => {
  for (const set of benchmark.sets.filter(s => s.mode === 'transfer-frozen')) {
    for (const row of set.rows.filter(r => r.status === 'scored')) {
      const forecast = predict(c, row.pKa, c.scope);
      for (const key of ['nearest', 'regression', 'local']) {
        assert.ok(Math.abs(forecast.ratios[key] / row.predictions[key] - 1) < 1e-10, set.id + ' ' + key);
      }
      assert.deepEqual(forecast.references, [...new Set([...row.lower, ...row.upper])]);
    }
  }
});

test('prior-release three-method card and score round-trip exactly', () => {
  // Snapshot produced by the published bde6ae9 implementation.
  const trial = {plan: {schema:'blink-atlas/probe-trial/1', model:'holland-h2b-pka-log-ratio-v1',
    sourceSha256:'1e26a58d558f7d96518aca0264b0efed6e5574b67ba5f18899965802fedf5031',
    article:'https://doi.org/10.1038/s41592-026-03062-5', createdUTC:'2026-11-12T01:02:03.000Z',
    compound:'archived-01', scope:'hydroxymethyl-sir-monotonic-fixed-cos7-h2b',
    target:'Mean ON:OFF time ratio of the HaloTag ligand in fixed COS-7 H2B–HaloTag.', pKa:3.2,
    references:['JF635b','JF639b'], nearestDye:'JF635b',
    ratios:{regression:.0004297281986962477, nearest:.000483915, median:.0003744094454760885},
    scoringRule:'Three independent preparations. Average the three mean ON:OFF time ratios arithmetically. Compare absolute log10 errors to the frozen predictions. Lower is better.'},
  measurement:{ratios:[.0004,.0005,.0006], observed:.0005,
    error:{regression:.06577615138201898, nearest:.014200920114844084, median:.12562320789601866}, closest:['nearest']}};
  assert.deepEqual(restoreTrial(data, trial), trial);
  assert.deepEqual(scoreTrial(trial.plan, trial.measurement.ratios), trial.measurement);
});

test('unsupported assays and unmeasured covariates abstain', () => {
  for (const x of ['', '  ', NaN, true, null, Infinity, '0x3', '1.219', '4.481']) {
    assert.throws(() => predict(c, x, c.scope));
  }
  assert.throws(() => predict(c, 3, 'complex'));
  assert.throws(() => predict(c, 3, 'outside'));
  for (const id of ['', 'JF639b', 'jf-639B', '4HTL']) {
    assert.throws(() => lockTrial(data, id, 3, c.scope));
  }
});

test('three independent ratios determine the target and support a tie', () => {
  const p = create().plan, r = p.ratios.nearest;
  const m = scoreTrial(p, [r/2, r, 1.5*r]);
  assert.ok(Math.abs(m.observed - r) < 1e-15);
  assert.equal(m.error.nearest, 0);
  assert.deepEqual(m.closest, ['nearest']);
  const tie = Math.sqrt(p.ratios.nearest * p.ratios.regression);
  assert.deepEqual(scoreTrial(p, [tie,tie,tie]).closest, ['regression','nearest']);
  for (const values of [[], [1,2], [1,2,3,4], [1,'',2], [1,0,2], [1,-1,1]]) {
    assert.throws(() => scoreTrial(p, values));
  }
});

test('saved plans and measured cards round-trip; altered results and models are rejected', () => {
  const trial = create();
  assert.deepEqual(restoreTrial(data, structuredClone(trial)), trial);
  trial.measurement = scoreTrial(trial.plan, [.0004,.0005,.0006]);
  assert.deepEqual(restoreTrial(data, structuredClone(trial)), trial);
  const edits = [
    v => {v.plan.model = 'other';}, v => {v.plan.sourceSha256 = 'x';},
    v => {v.plan.ratios.regression *= 2;}, v => {v.plan.pKa = 2;},
    v => {v.plan.nearestDye = 'fake';}, v => {v.measurement.ratios[0] *= 2;},
    v => {v.measurement.error.median = 3;}, v => {v.measurement = {};},
  ];
  for (const edit of edits) {
    const copy = structuredClone(trial); edit(copy);
    assert.throws(() => restoreTrial(data, copy));
  }
});
