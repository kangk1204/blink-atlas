import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { compareFields, validateComparison, validateField } from "../docs/precision-compare.js";

function equal(actual,expected) {
  if (typeof expected==="number") {
    assert.ok(Math.abs(actual-expected) <= 1e-9*Math.max(1,Math.abs(expected)));
  } else if (Array.isArray(expected)) {
    assert.equal(actual.length,expected.length);
    expected.forEach((v,i)=>equal(actual[i],v));
  } else if (expected && typeof expected==="object") {
    assert.deepEqual(Object.keys(actual).sort(),Object.keys(expected).sort());
    Object.keys(expected).forEach(k=>equal(actual[k],expected[k]));
  } else assert.equal(actual,expected);
}
const cases=JSON.parse(readFileSync(0,"utf8"));
for (const {fields,dyes,expected} of cases) {
  const result=compareFields(fields,dyes);
  equal(validateComparison(result),expected);
}
const field=structuredClone(cases[0].fields[0]);
field.denominators.elapsed_s*=2;
assert.throws(()=>validateField(field),/denominators/);
const comparison=structuredClone(cases[0].expected);
comparison.preparations[0].denominator.qualified=9999;
assert.throws(()=>validateComparison(comparison),/counts|rates/);
