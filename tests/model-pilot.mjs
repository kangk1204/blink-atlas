import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import { evaluate, probabilities } from "../docs/model-core.js";

const read = path => JSON.parse(fs.readFileSync(new URL(path, import.meta.url), "utf8"));
const dyes = read("../docs/data/dyes.json"), model = read("../docs/data/model.json");
const settings = { densityMin: 10, densityMax: 50000 };

// A matched synthetic pilot, used solely as a known-input software test.
function matched(rho = 10000) {
  const m = dyes.dyes.find(d => d.id === "JF630b").metrics;
  const p = probabilities(m.onMs.mean, m.offMs.mean, m.photonsPerMs.mean, 25, 100);
  const mu = rho * p.qualified * Math.exp(-rho * Math.PI * 0.25**2 * p.any) * 0.0625;
  const frames = 200, tilesPerFrame = 400, total = frames * tilesPerFrame;
  const empty = Math.round(total * Math.exp(-mu));
  const single = Math.round(total * mu * Math.exp(-mu));
  const multiple = total - empty - single, locs = single + 2 * multiple;
  return {
    schema: "smlm-pilot/1",
    acquisition: {
      sample_state: "fixed", cell_line: "COS-7", target: "H2B", tag: "HaloTag",
      reference_dye: "JF630b", exposure_ms: 25, frame_interval_ms: 27, excitation_nm: 647,
      width_nm: 5000, height_nm: 5000, total_frames: 5000
    },
    pilot: {
      start_frame: 0, min_photons: 100, frames, area_um2: 25, localizations: locs,
      tiles: { side_nm: 250, area_um2: 0.0625, tiles_per_frame: tilesPerFrame,
               observations: total, empty, single, multiple, localizations: locs }
    }
  };
}

test("a matched pilot conditions the broad density prior with a single tile likelihood", () => {
  const pilot = matched(), updated = evaluate(dyes, model, settings, pilot);
  assert.equal(updated.pilot.status, "applied");
  const [lo, hi] = updated.density;
  assert.ok(lo > settings.densityMin && lo <= 10000 && hi >= 10000 && hi < settings.densityMax);
  assert.equal(updated.pilot.effectiveTileObservations, 400);
  // More detections within the multiple class add no second copy of its evidence.
  pilot.pilot.localizations += 500;
  pilot.pilot.tiles.localizations += 500;
  assert.deepEqual(evaluate(dyes, model, settings, pilot).pilot, updated.pilot);
});

test("matched calibration uses only the pilot and retains identical measurements for all dyes", () => {
  const pilot = matched(1000), base = JSON.stringify(dyes);
  const before = evaluate(dyes, model, settings, pilot);
  pilot.evaluation = { blocks: [{ observed: 900000 }, { observed: 0 }] };
  assert.deepEqual(before, evaluate(dyes, model, settings, pilot));
  assert.equal(JSON.stringify(dyes), base);
});

test("live public case, different photon cutoffs and inconsistent data withhold recommendations", () => {
  const external = read("../docs/data/pilot.json");
  const wrongThreshold = matched(), wrongUnits = matched(), empty = matched(), incompatible = matched();
  wrongThreshold.pilot.min_photons = 0;
  wrongUnits.pilot.tiles.area_um2 = 62500;
  empty.pilot.tiles.empty = empty.pilot.tiles.observations;
  empty.pilot.tiles.single = empty.pilot.tiles.multiple = empty.pilot.tiles.localizations = empty.pilot.localizations = 0;
  const t = incompatible.pilot.tiles;
  t.empty = t.single = 0; t.multiple = t.observations;
  t.localizations = incompatible.pilot.localizations = 2 * t.observations;
  for (const data of [external, wrongThreshold, wrongUnits, empty, incompatible]) {
    const output = evaluate(dyes, model, settings, data);
    assert.equal(output.pilot.status, "rejected");
    assert.deepEqual(output.recommendation, []);
    assert.ok(output.pilot.reason.length > 0);
  }
});

test("invalid numeric settings and unsupported acquisition have distinct responses", () => {
  for (const s of [{ minPhotons: NaN }, { exposureMs: 40 }, { densityMin: -1 }, {densityMax: 0}]) {
    assert.equal(evaluate(dyes, model, s).status, "invalid");
  }
  for (const s of [{cellLine:"U2OS"}, {excitationNm: 640}, {frameIntervalMs: 100}]) {
    assert.equal(evaluate(dyes, model, s).status, "out_of_scope");
  }
});
