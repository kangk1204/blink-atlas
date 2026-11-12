import assert from 'node:assert/strict';
import {test} from 'node:test';
import {readFile} from 'node:fs/promises';
import {lockTrial, predict, restoreTrial, scoreTrial} from '../docs/probe-trial.js';

const data = JSON.parse(await readFile(new URL('../docs/data/probe-design.json', import.meta.url), 'utf8'));
const c = data.calibration;
const create = () => lockTrial(data, 'blind-01', '3.20', c.scope, '2026-11-12T01:02:03.000Z');

test('locked predictions, endpoints and a deterministic nearest reference at a midpoint', () => {
  const p = create().plan;
  assert.equal(p.nearestDye, 'JF635b');
  assert.deepEqual(p.references, ['JF635b','JF639b']);
  assert.ok(Math.abs(p.ratios.regression - .0004297281986991795) < 1e-14);
  assert.equal(p.ratios.nearest, .000483915);
  assert.equal(predict(c, 3.46, c.scope).nearestDye, 'JF635b');
  assert.deepEqual(predict(c, 1.22, c.scope).references, ['JF626b']);
  assert.deepEqual(predict(c, 4.48, c.scope).references, ['JFX650b']);
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
