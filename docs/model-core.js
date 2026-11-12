/** Finite-exposure event-yield model. All state and source records are passed in. */
const unitKey = value => String(value ?? "").toLowerCase().replace(/[^a-z0-9]/g, "");
const positive = value => Number.isFinite(value) && value > 0;
const count = value => Number.isSafeInteger(value) && value >= 0;

function emptyResult(settings) {
  const density = [settings.densityMin, settings.densityMax];
  return {
    status: "invalid", recommendation: [], contenders: [], ranking: [], density,
    centerDensity: density.every(positive) ? Math.sqrt(density[0] * density[1]) : null,
    reasons: [], nextMeasurement: "",
    pilot: { status: "none", reason: "Publication measurements only.", interval: null }
  };
}

function settingsError(s, config) {
  if (!["all","validation"].includes(s.candidatePanel)) return "Choose the full dye panel or the three-ligand validation panel.";
  const numbers = [s.densityMin, s.densityMax, s.minPhotons, s.exposureMs, s.frameIntervalMs, s.excitationNm];
  if (!numbers.every(positive)) return "A finite positive number is required for each acquisition setting.";
  const [lower, upper] = config.limits.density;
  if (s.densityMin < lower || s.densityMin > s.densityMax || s.densityMax > upper)
    return `Enter an ordered density range between ${lower} and ${upper} labels/µm².`;
  const [minSignal, maxSignal] = config.limits.minPhotons;
  if (s.minPhotons < minSignal || s.minPhotons > maxSignal)
    return `Choose a photon threshold between ${minSignal} and ${maxSignal}.`;
  if (s.exposureMs > s.frameIntervalMs) return "Exposure must fit within the camera frame interval.";
  return null;
}

function inMeasuredDomain(s, defaults) {
  const keys = ["sampleState", "cellLine", "target", "tag"];
  return keys.every(k => unitKey(s[k]) === unitKey(defaults[k])) &&
    Number.isFinite(s.excitationNm) && Math.abs(s.excitationNm - defaults.excitationNm) <= 1 &&
    Number.isFinite(s.exposureMs) && Math.abs(s.exposureMs - defaults.exposureMs) <= 0.2 &&
    Number.isFinite(s.frameIntervalMs) && Math.abs(s.frameIntervalMs - defaults.frameIntervalMs) <= 0.2;
}

export function probabilities(on, off, photonsPerMs, exposure, threshold) {
  const p = on / (on + off);
  const firstBy = time => p - (1 - p) * Math.expm1(-time / off);
  const needed = threshold / photonsPerMs;
  return {
    any: firstBy(exposure),
    qualified: needed > exposure ? 0 : Math.exp(-needed / on) * firstBy(exposure - needed)
  };
}

function central(dye, s) {
  const { onMs, offMs, photonsPerMs } = dye.metrics;
  return probabilities(onMs.mean, offMs.mean, photonsPerMs.mean, s.exposureMs, s.minPhotons);
}

function sensitivity(dye, config, s) {
  const m = dye.metrics;
  return config.offTimeScales.map(scale => {
    const values = [];
    for (const on of m.onMs.replicates)
      for (const off of m.offMs.replicates)
        for (const photons of m.photonsPerMs.replicates)
          values.push(probabilities(on, off * scale, photons, s.exposureMs, s.minPhotons));
    return values;
  });
}

function area(config) {
  return Math.PI * (config.collisionRadiusNm / 1000) ** 2;
}

function eligibleDyes(dyeData, config, s) {
  return s.candidatePanel === "validation" ? dyeData.dyes.filter(d => config.validationDyes.includes(d.id)) : dyeData.dyes;
}

function rate(rho, p, disc, cycle) {
  return rho * p.qualified * Math.exp(-rho * disc * p.any) * 1000 / cycle;
}

function normRange(grid, accum, probability) {
  let sum = 0;
  for (let i = 0; i < accum.length; i++) {
    sum += accum[i];
    if (sum >= probability) return grid[i];
  }
  return grid[grid.length - 1];
}

function mixtureLogLikelihood(logs) {
  const largest = Math.max(...logs);
  return largest + Math.log(logs.reduce((v, x) => v + Math.exp(x - largest), 0) / logs.length);
}

function occupancyTerms(mean) {
  const p0 = Math.exp(-mean), p1 = mean * p0;
  const pMany = mean < 0.0001 ? mean * mean * (0.5 - mean / 3 + mean * mean / 8)
                               : -Math.expm1(-mean) - p1;
  return [p0, p1, pMany];
}

function validatePilot(input, s, defaults, dyes) {
  if (!input || input.schema !== "smlm-pilot/1") return "Select a smlm-pilot/1 localization summary.";
  const a = input.acquisition || {}, p = input.pilot || {}, t = p.tiles;
  const mapped = {
    sampleState: a.sample_state, cellLine: a.cell_line, target: a.target, tag: a.tag,
    excitationNm: a.excitation_nm, exposureMs: a.exposure_ms, frameIntervalMs: a.frame_interval_ms
  };
  if (!inMeasuredDomain(mapped, defaults))
    return "The reference pilot must use fixed COS-7, H2B–HaloTag, 647 nm, 25 ms exposure and a 27 ms camera cycle.";
  if (!dyes.some(d => unitKey(d.id) === unitKey(a.reference_dye))) return "Supply a measured reference dye.";
  if (p.min_photons !== s.minPhotons) return "Use the same photon threshold in the pilot importer and candidate comparison.";
  if (!positive(a.width_nm) || !positive(a.height_nm) || !count(a.total_frames) ||
      !count(p.start_frame) || !count(p.frames) || p.frames < 1 || p.start_frame + p.frames > a.total_frames)
    return "Pilot field dimensions or acquisition frame counts are missing or invalid.";
  if (!t || ![t.empty, t.single, t.multiple, t.localizations, t.observations, t.tiles_per_frame, p.localizations].every(count) ||
      !positive(t.side_nm) || !positive(t.area_um2) ||
      Math.abs(t.area_um2 - (t.side_nm / 1000) ** 2) > t.area_um2 * 1e-5 ||
      !positive(p.area_um2) || Math.abs(p.area_um2 - a.width_nm * a.height_nm / 1e6) > p.area_um2 * 1e-5 ||
      t.tiles_per_frame !== Math.floor(a.width_nm / t.side_nm + 1e-10) * Math.floor(a.height_nm / t.side_nm + 1e-10) ||
      t.observations < 1 || t.observations !== t.empty + t.single + t.multiple ||
      t.observations !== t.tiles_per_frame * p.frames || t.localizations > p.localizations ||
      t.localizations < t.single + 2 * t.multiple || (t.multiple === 0 && t.localizations !== t.single))
    return "The pilot needs consistent counts from complete equal-area square tiles, including empty frames.";
  if (t.single + t.multiple === 0) return "No qualifying localizations were detected in the complete pilot tiles.";
  return null;
}

/** Log-uniform grid prior; a single joint tile histogram likelihood avoids double-counting. */
function updateDensity(input, s, config, dye) {
  const t = input.pilot.tiles, obs = [t.empty, t.single, t.multiple].map(n => n / t.observations);
  const cfg = config.pilotConfig;
  const effective = Math.min(t.tiles_per_frame, cfg.maxEffectiveTiles);
  const saturated = effective * obs.reduce((ll, p) => ll + (p ? p * Math.log(p) : 0), 0);
  const reference = sensitivity(dye, config, s).flat();
  const disc = area(config);
  const steps = s.densityMin === s.densityMax ? 1 : cfg.gridPoints;
  const grid = [], logs = [], priors = [];
  let maxLL = -Infinity;
  for (let i = 0; i < steps; i++) {
    const f = steps === 1 ? 0 : i / (steps - 1);
    const rho = i === steps - 1 ? s.densityMax : s.densityMin * (s.densityMax / s.densityMin) ** f;
    grid.push(rho);
    const scenarioLL = reference.map(p => {
      const mu = rho * p.qualified * Math.exp(-rho * disc * p.any) * t.area_um2;
      const probabilities = occupancyTerms(mu);
      return effective * obs.reduce((ll, observed, k) =>
        ll + observed * Math.log(Math.max(1e-300, probabilities[k])), 0);
    });
    maxLL = Math.max(maxLL, ...scenarioLL);
    // Trapezoid weights in log density. All kinetic scenarios get equal mass.
    const priorWeight = steps > 1 && (i === 0 || i === steps - 1) ? 0.5 : 1;
    priors.push(priorWeight);
    logs.push(mixtureLogLikelihood(scenarioLL) + Math.log(priorWeight));
  }
  const deviance = Math.max(0, 2 * (saturated - maxLL));
  if (deviance > cfg.maxDeviance)
    return { status: "rejected", interval: null, deviance, reason:
      "Observed tile counts fall outside the reference sensitivity scenarios. Check density, calibration, background and field selection." };
  const max = Math.max(...logs), raw = logs.map(ll => Math.exp(ll - max));
  const total = raw.reduce((a, b) => a + b), posterior = raw.map(v => v / total);
  const priorSum = priors.reduce((x, y) => x + y);
  const informationNats = posterior.reduce((kl, p, i) =>
    kl + (p > 0 ? p * Math.log(p * priorSum / priors[i]) : 0), 0);
  const tail = (1 - cfg.credibleMass) / 2;
  // A weak pilot must not narrow the original band merely by trimming prior tails.
  const informative = informationNats >= cfg.minInformationNats;
  const interval = informative ? [normRange(grid, posterior, tail), normRange(grid, posterior, 1 - tail)]
                               : [s.densityMin, s.densityMax];
  return {
    status: "applied", interval, referenceDye: dye.id, effectiveTileObservations: effective,
    informationNats, deviance, credibleMass: cfg.credibleMass,
    reason: informative ? "Matched pilot updated the effective density range (90% conditional posterior interval)."
                        : "The pilot supplies little additional density information; the full input range is retained."
  };
}

function dyeBounds(dye, config, s, lo, hi, center, disc) {
  const nominal = central(dye, s);
  let lower = Infinity, upper = 0;
  const endpoints = sensitivity(dye, config, s).map(scenarios => {
    const ext = [[Infinity, 0], [Infinity, 0]];
    for (const p of scenarios) {
      const atLo = rate(lo, p, disc, s.frameIntervalMs), atHi = rate(hi, p, disc, s.frameIntervalMs);
      ext[0][0] = Math.min(ext[0][0], atLo);
      ext[0][1] = Math.max(ext[0][1], atLo);
      ext[1][0] = Math.min(ext[1][0], atHi);
      ext[1][1] = Math.max(ext[1][1], atHi);
      lower = Math.min(lower, atLo, atHi);
      const peak = Math.max(lo, Math.min(hi, 1 / (disc * p.any)));
      upper = Math.max(upper, rate(peak, p, disc, s.frameIntervalMs));
    }
    return ext;
  });
  return {
    endpoints, row: {
      id: dye.id, name: dye.name, rate: rate(center, nominal, disc, s.frameIntervalMs),
      low: lower, high: upper, detectProbability: nominal.qualified, anyOnProbability: nominal.any,
      overlapProbability: -Math.expm1(-center * disc * nominal.any)
    }
  };
}

export function evaluate(dyeData, config, overrides = {}, pilot = null) {
  const s = { ...config.defaults, ...overrides };
  const result = emptyResult(s);
  const error = settingsError(s, config);
  if (error) {
    result.reasons = [error];
    result.nextMeasurement = "Check the input and units.";
    return result;
  }
  if (!inMeasuredDomain(s, config.defaults)) {
    result.status = "out_of_scope";
    result.reasons = ["The measured reference covers fixed COS-7 H2B–HaloTag at 647 nm, 25 ms exposure and a 27 ms frame interval."];
    result.nextMeasurement = "Measure on/off durations and photon rates in this acquisition.";
    return result;
  }
  if (pilot !== null) {
    const invalid = validatePilot(pilot, s, config.defaults, dyeData.dyes);
    result.pilot = invalid ? { status: "rejected", reason: invalid, interval: null } :
      updateDensity(pilot, s, config, dyeData.dyes.find(d => unitKey(d.id) === unitKey(pilot.acquisition.reference_dye)));
    if (result.pilot.status === "rejected") {
      result.status = "out_of_scope";
      result.reasons = [result.pilot.reason];
      result.nextMeasurement = "Use a matched reference pilot or return to the publication-only scenario.";
      return result;
    }
    result.density = result.pilot.interval;
    result.centerDensity = Math.sqrt(result.density[0] * result.density[1]);
  }
  const bounds = new Map(), disc = area(config);
  for (const dye of eligibleDyes(dyeData, config, s)) {
    const { endpoints, row } = dyeBounds(dye, config, s, ...result.density, result.centerDensity, disc);
    bounds.set(dye.id, endpoints);
    result.ranking.push(row);
  }
  result.ranking.sort((a, b) => b.rate - a.rate);
  function dominates(a, b) {
    // Log rate ratios are affine in density for fixed kinetic settings.
    return bounds.get(a).every((edges, j) =>
      edges.every((limits, k) => limits[0] > bounds.get(b)[j][k][1] + 1e-12));
  }
  const rows = result.ranking, pair = rows.slice(0, 2).map(row => row.id);
  result.contenders = rows.filter(a => a.high > 0 && rows.filter(b => a.id !== b.id && dominates(b.id, a.id)).length < 2)
                          .map(a => a.id);
  if (rows.filter(row => row.high > 0).length < 2) {
    result.status = "low_signal";
    result.reasons = ["Fewer than two dyes can reach the required expected photon count in one exposure."];
    result.nextMeasurement = "Measure the photon distribution and camera background in a sparse field before setting the detection threshold.";
  } else if (rows[1].low > 0 && pair.every(id => rows.filter(row => !pair.includes(row.id)).every(row => dominates(id, row.id)))) {
    result.status = "stable";
    result.recommendation = pair;
    result.reasons = ["The same unordered pair leads across the density interval, measured replicate combinations and two shared dark-time assumptions."];
    result.nextMeasurement = "Test these two ligands in matched wells at equal acquisition times.";
  } else {
    result.status = "uncertain";
    result.reasons = ["The top two can change across the measurement and density sensitivity checks."];
    result.nextMeasurement = pilot === null ?
      "Use a sparse reference field at the same exposure and photon threshold to refine effective density." :
      "Collect matched dark-time and photon measurements for the competing ligands.";
  }
  return result;
}

/** Nominal (publication-mean) yield curves, independent of a pilot's density interval. */
export function curve(dyeData, config, overrides = {}, densities = []) {
  const s = { ...config.defaults, ...overrides }, disc = area(config);
  if (settingsError(s, config) || !inMeasuredDomain(s, config.defaults)) return [];
  const ps = eligibleDyes(dyeData, config, s).map(d => [d.id, central(d, s)]);
  return densities.filter(rho => Number.isFinite(rho) && rho >= config.limits.density[0] && rho <= config.limits.density[1])
    .map(density => ({
      density, rates: Object.fromEntries(ps.map(([id, p]) => [id, rate(density, p, disc, s.frameIntervalMs)]))
    }));
}
