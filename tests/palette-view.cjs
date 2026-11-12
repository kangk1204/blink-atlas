// The published catalog is real. Test-only fixtures are served from memory on
// loopback. Pass "msedge --release" to additionally verify the release JSON.
const assert = require("node:assert/strict");
const fs = require("node:fs/promises");
const http = require("node:http");
const path = require("node:path");
const {chromium} = require("playwright");

const docs = path.join(__dirname, "..", "docs");
const metric = "nenaNm";
const specimen = "origami";
const dyeId = "janeliafluor646";
const html = `<!DOCTYPE html><html lang="en"><meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<link rel="stylesheet" href="style.css"><link rel="stylesheet" href="palette-view.css">
<main class="container"><h1>DNA-PAINT test fixture</h1><div id="mount"></div></main>
<script type="module">
import {renderPalette} from './palette-view.js';
const target=document.querySelector('#mount'); window.ready=renderPalette(target);
window.reused=renderPalette(target)===window.ready;
</script></html>`;

function rawFixture(catalog) {
  const replicates = [];
  for (const dye of catalog.dyes) {
    for (const [system, measurement] of Object.entries(dye.systems)) {
      for (const {repeat, value} of measurement.metrics[metric].points) {
        replicates.push({
          system, dye: dye.id, repeat, principal: repeat <= 3,
          metrics: {nenaNm: repeat === 4 ? 90 : value + 0.0246, darkTimeS: 5, eventPhotons: 950},
          files: [{role: "target_field", path: "TEST_ONLY/" + dye.id + "/n" + repeat + "/field.hdf5", rows: 321}],
          counts: {centerFrames: 11}, fieldFits: {target_field: {pairCount: 80, relativeResidualL2: 0.0321}},
          qc: ["Test-only quality note."],
        });
      }
    }
  }
  return {replicates, metrics: {darkTimeS: "test-only dark-time definition"}, coverage: {
    expectedCases: replicates.length, completedCases: replicates.length, totalFiles: replicates.length,
  }};
}

function evaluationFixture(catalog, passed) {
  function result(ok) {
    const channels = [488, 560, 640].map(wavelength => {
      const ids = catalog.dyes.filter(d => d.wavelength === wavelength).map(d => d.id);
      return {wavelength, pair: ids.slice(0, 2), baselinePair: ids.slice(1, 3),
        oraclePair: ids.slice(-2), meanNm: 5.1, baselineMeanNm: ok ? 6.3 : 4.4,
        oracleMeanNm: 4, regretNm: 1.1, baselineRegretNm: ok ? 2.3 : 0.4};
    });
    return {status: "evaluated", gatePassed: ok, weight: .5, improvedChannels: ok ? 3 : 0,
      macroRegretNm: 1.1, baselineMacroRegretNm: ok ? 2.3 : .4, channels,
      reasons: ok ? [] : ["Test-only failed gate."]};
  }
  return {raw: result(passed), published: result(!passed)};
}

async function main() {
  const catalog = JSON.parse(await fs.readFile(path.join(docs, "data", "palette.json"), "utf8"));
  const raw = rawFixture(catalog);
  let mode = "held";
  const server = http.createServer(async (req, res) => {
    const u = req.url.split("?")[0];
    function respond(status, type, body) {
      res.writeHead(status, {"content-type": type, "cache-control": "no-store"});
      res.end(body);
    }
    if (u === "/") return respond(200, "text/html; charset=utf-8", html);
    if (mode === "release" && ["/data/palette-raw.json", "/data/palette-evaluation.json", "/data/palette-audit.json"].includes(u)) {
      try {return respond(200, "application/json", await fs.readFile(path.join(docs, u.slice(1))));}
      catch {return respond(404, "text/plain", "Missing release asset");}
    }
    if (u === "/data/palette-raw.json") return mode === "missing"
      ? respond(404, "text/plain", "Missing test raw") : respond(200, "application/json", JSON.stringify(raw));
    if (u === "/data/palette-audit.json") {
      if (mode === "missing") return respond(404, "text/plain", "Missing test audit");
      return respond(200, "application/json", JSON.stringify({
        summaryCellsReproduced: 123, hdf5Members: 222, mainFigureHdf5: 210,
        metrics: [{system: specimen, metric, sourceCells: 54, comparable: 53, withinHalfGrid: 23, gridUnit: .01,
          medianAbsoluteDifference: .006, meanSignedDifference: -.003, largestDifferences: [{
            dye: dyeId, repeat: 2, cell: "R8", published: 4.34, recomputed: 4.5, difference: .16,
          }], missing: []}],
      }));
    }
    if (u === "/data/palette-evaluation.json") {
      const evaluation = evaluationFixture(catalog, mode === "passed");
      if (mode === "missing") evaluation.raw = {status: "unavailable", gatePassed: false, channels: [], reasons: ["Missing test group."]};
      return respond(200, "application/json", JSON.stringify(evaluation));
    }
    const allowed = {
      "/palette-view.js": "text/javascript", "/palette-view.css": "text/css", "/style.css": "text/css",
      "/data/palette.json": "application/json", "/assets/dm-sans.woff2": "font/woff2",
      "/assets/plex-mono.woff2": "font/woff2",
    };
    if (!allowed[u]) return respond(404, "text/plain", "Missing test URL");
    try {respond(200, allowed[u], await fs.readFile(path.join(docs, u.slice(1))));}
    catch {respond(500, "text/plain", "Could not read local asset");}
  });
  await new Promise(resolve => server.listen(0, "127.0.0.1", resolve));
  let browser;
  try {
    const channel = process.argv.slice(2).find(arg => !arg.startsWith("--"));
    browser = await chromium.launch({headless: true, ...(channel ? {channel} : {})});
    const page = await browser.newPage({viewport: {width: 1366, height: 900}});
    const errors = [];
    page.on("pageerror", e => errors.push(e.message));
    const load = async () => {
      await page.goto("http://127.0.0.1:" + server.address().port + "/");
      await page.locator(".palette-series-row").first().waitFor();
      await page.evaluate(() => window.ready);
    };
    await load();
    assert.equal(await page.evaluate(() => window.reused), true);
    assert.equal(await page.locator("#palette-wave").inputValue(), "640");
    assert.equal(await page.locator("#palette-dye").inputValue(), dyeId);
    assert.equal(await page.locator(".palette-series-row").count(), 10);
    const first = page.locator("#palette-repeats tbody tr").first();
    const expected = catalog.dyes.find(d => d.id === dyeId).systems[specimen].metrics[metric].points[0];
    assert.equal(await first.locator("td").nth(0).innerText(), String(expected.value));
    assert.equal(await first.locator("td").nth(1).innerText(), (expected.value + .0246).toFixed(4));
    assert.equal(await first.locator("a").innerText(), "Sheet1!" + expected.cell);
    assert.match(await first.locator("a").getAttribute("href"), /data\/steen-figure2\.xlsx$/);
    assert.match(await page.locator("#palette-agreement-summary").innerText(), /23\/53/);
    await page.locator("#palette-agreement-summary").click();
    const agreement = page.locator("#palette-agreement-body");
    assert.match(await agreement.innerText(), /53 of 54 workbook repeat values/);
    assert.match(await agreement.innerText(), /median absolute difference is 0.006 nm/);
    assert.match(await agreement.innerText(), /grid has a step of 0.01 nm/);
    assert.match(await agreement.innerText(), /JF646, repeat 2/);
    assert.equal(await agreement.locator("a").innerText(), "Sheet1!R8");
    assert.match(await page.locator("#palette-coverage").innerText(), /123 workbook means/);
    await page.locator("summary", {hasText: "Files, group definitions"}).click();
    await page.locator("#palette-raw-inspect details > summary").first().click();
    assert.equal(await page.locator("#palette-raw-inspect code").first().isVisible(), true);
    assert.match(await page.locator("#palette-raw-inspect").innerText(), /80 adjacent-frame pairs/);
    assert.match(await page.locator("#palette-raw-inspect").innerText(), /Test-only quality note/);
    await page.locator("#palette-metric").selectOption("darkTimeS");
    assert.equal(await page.locator("#palette-plot-title").innerText(), "Dark time");
    assert.match(await page.locator("#palette-agreement-summary").innerText(), /HDF5 only/);
    assert.equal(await first.locator("td").nth(0).innerText(), "—");
    assert.equal(await first.locator("td").nth(1).innerText(), "5");
    assert.match(await page.locator("#palette-dye-note").innerText(), /test-only dark-time definition/);
    await page.locator("#palette-metric").selectOption("eventPhotons");
    assert.equal(await page.locator("#palette-plot-unit").innerText(), "photons");
    await page.locator("#palette-wave").selectOption("560");
    assert.equal(await page.locator(".palette-series-row").count(), 4);
    await page.locator("#palette-system").selectOption("npc");
    await page.locator("#palette-dye").selectOption("cy3b");
    await page.locator("#palette-metric").selectOption(metric);
    assert.equal(await page.locator("#palette-repeats tbody tr").count(), 4);
    assert.match(await page.locator("#palette-repeats tbody tr").last().innerText(), /extra/);
    assert.equal(await page.locator(".palette-series-row[data-dye=cy3b] .palette-number").last().innerText(), "4.71");
    await page.locator(".palette-series-row[data-dye=cf568]").focus();
    await page.keyboard.press("Enter");
    assert.equal(await page.locator("#palette-dye").inputValue(), "cf568");
    assert.equal(await page.locator(".palette-series-row[aria-pressed=true]").count(), 1);
    assert.match(await page.locator("#palette-gate").innerText(), /withheld/i);
    assert.equal(await page.locator("#palette-eval-results tbody tr").count(), 3);
    await page.locator("#palette-eval-source").selectOption("published");
    assert.match(await page.locator("#palette-eval-results").innerText(), /replay meets/);
    assert.match(await page.locator("#palette-gate").innerText(), /withheld/i);
    const urls = await page.locator("#palette-sources a").evaluateAll(as => as.map(a => a.href));
    assert(urls.some(u => u.endsWith("/analysis/palette-protocol.md")));
    assert(urls.some(u => u.endsWith("/data/palette-manifest.json")));
    assert(urls.some(u => u.endsWith("/data/palette-audit.json")));
    assert(urls.some(u => u.endsWith("/analysis/palette-raw.md")));
    for (const width of [320, 390, 850, 1366]) {
      await page.setViewportSize({width, height: 850});
      assert(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth), "Page overflow at " + width);
    }
    if (process.env.BLINK_PALETTE_SCREEN) {
      await load();
      await page.screenshot({path: process.env.BLINK_PALETTE_SCREEN, fullPage: true});
    }
    mode = "passed";
    await load();
    assert.equal(await page.locator("#palette-gate").getAttribute("data-state"), "passed");
    mode = "missing";
    await load();
    assert.equal(await page.locator("#palette-eval-source").inputValue(), "published");
    assert.match(await page.locator("#palette-gate").innerText(), /withheld/i);
    assert.match(await page.locator("#palette-coverage").innerText(), /could not be loaded/);
    assert.match(await page.locator("#palette-agreement-summary").innerText(), /audit unavailable/);
    assert.equal(await page.locator("#palette-metric option[value=darkTimeS]").count(), 0);
    await page.locator("#palette-eval-source").selectOption("raw");
    assert.match(await page.locator("#palette-eval-results").innerText(), /Missing test group/);
    if (process.argv.includes("--release")) {
      const read = async file => JSON.parse(await fs.readFile(path.join(docs, "data", file), "utf8"));
      const [measurements, audit, evaluation] = await Promise.all([
        read("palette-raw.json"), read("palette-audit.json"), read("palette-evaluation.json"),
      ]);
      mode = "release";
      await load();
      const expectedRow = measurements.replicates.find(r => r.system === specimen && r.dye === dyeId && r.repeat === 1);
      const displayed = Number((await first.locator("td").nth(1).innerText()).replaceAll(",", ""));
      assert(Math.abs(displayed - expectedRow.metrics[metric]) < .000051);
      const agreementRow = audit.metrics.find(r => r.system === specimen && r.metric === metric);
      assert((await page.locator("#palette-agreement-summary").innerText()).includes(agreementRow.withinHalfGrid + "/" + agreementRow.comparable));
      assert.equal(await page.locator("#palette-gate").getAttribute("data-state"), evaluation.raw.gatePassed ? "passed" : "held");
      await page.locator("#palette-wave").selectOption("488");
      await page.locator("#palette-system").selectOption("npc");
      await page.locator("#palette-metric").selectOption("specificity");
      assert.equal(await page.locator("#palette-repeats tbody tr").count(), 3);
      console.log("Release assets matched the rendered measurement, source agreement and recommendation gate.");
    }
    assert.deepEqual(errors, []);
    console.log("DNA-PAINT view passed desktop, mobile, source trace, holdout and unavailable-input fixture checks.");
  } finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
}
main().catch(e => {console.error(e); process.exitCode = 1;});
