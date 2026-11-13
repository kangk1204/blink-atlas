// Exercise the anonymous public page, browser-side model check and downloadable structures.
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const path = require('node:path');
const {chromium} = require('playwright');

const base = process.argv[2] || 'http://127.0.0.1:8187/';
const directory = path.resolve(__dirname, '..', 'docs', 'data', 'structure-pk');
(async () => {
  const browser = await chromium.launch({headless: true, channel: process.argv[3] || 'chrome'});
  try {
    const session = await browser.newContext({viewport: {width: 1380, height: 920}, acceptDownloads: true});
    const page = await session.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push(e.message));
    await page.goto(new URL('#probe-design', base).href);
    await page.locator('.pk-candidate').first().waitFor();
    assert.equal(await page.locator('.pk-candidate').count(), 2);
    const work = page.locator('#pk-workbench');
    const bundle = JSON.parse(await fs.readFile(path.join(directory, 'benchmark.json')));
    for (const c of bundle.candidates) {
      for (const ext of ['sdf', 'svg', 'json']) {
        const filename = c.dye + '.' + ext;
        const a = work.locator('a[href$="' + filename + '"]').first();
        const url = new URL(await a.getAttribute('href'), base).href;
        const response = await session.request.get(url);
        assert.ok(response.ok(), url);
        assert.deepEqual(await response.body(), await fs.readFile(path.join(directory, filename)));
      }
    }
    await work.locator('#pk-evaluation > summary').click();
    assert.match(await work.locator('#pk-evaluation').innerText(), /Yale676sb reverses/);
    const transferDir = path.join(directory, '..', 'structure-transfer');
    const next = JSON.parse(await fs.readFile(path.join(transferDir, 'benchmark.json')));
    await work.locator('#pk-transfer > summary').click();
    const readout = await work.locator('#pk-transfer').innerText();
    for (const f of [next.folds[0], next.folds[1], next.newParents]) {
      for (const m of ['nearest', 'kernel', 'coreDonor', 'trainMean']) assert.ok(readout.includes(f.summary.metrics[m].mae.toFixed(3)));
    }
    for (const filename of await fs.readdir(transferDir)) {
      const res = await session.request.get(new URL('data/structure-transfer/' + filename, base).href);
      assert.ok(res.ok(), filename);
      assert.deepEqual(await res.body(), await fs.readFile(path.join(transferDir,filename)));
    }
    await work.locator('#pk-import > summary').click();
    async function upload(value, pattern) {
      await work.locator('#pk-file').setInputFiles({name: 'roundtrip.json', mimeType: 'application/json',
        buffer: Buffer.from(JSON.stringify(value))});
      await page.waitForFunction(([s]) => new RegExp(s).test(document.getElementById('pk-file-status').textContent), [pattern.source]);
    }
    for (const c of bundle.candidates) {
      const data = JSON.parse(await fs.readFile(path.join(directory, c.dye + '.json')));
      await upload(data, new RegExp(c.dye + ' · Nearest ' + c.nearest.toFixed(2)));
    }
    const candidate = JSON.parse(await fs.readFile(path.join(directory, 'HM-SiR-CH2F.json')));
    await upload({...candidate, ridge: candidate.ridge + 0.5}, /disagree with the released model/);
    await upload({...candidate, modelId: 'old-release'}, /does not match/);
    const forecast = next.candidates[0];
    await upload(forecast, new RegExp('Core and donor ' + forecast.result.predictions.coreDonor.toFixed(2)));
    const altered = JSON.parse(JSON.stringify(forecast));
    altered.result.predictions.coreDonor += 0.5;
    await upload(altered, /forecasts disagree/);
    if (process.argv[6]) {
      await upload(JSON.parse(await fs.readFile(process.argv[6])), /single midpoint is withheld/);
    }
    if (process.argv[5]) {
      const complex = JSON.parse(await fs.readFile(process.argv[5]));
      await upload(complex, /single midpoint is withheld/);
    }
    await work.locator('#pk-evaluation > summary').click();
    for (const width of [320, 390, 760, 1000, 1380]) {
      await page.setViewportSize({width, height: 900});
      assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'Horizontal overflow at ' + width);
    }
    if (process.argv[4]) {
      await fs.mkdir(process.argv[4], {recursive:true});
      await work.scrollIntoViewIfNeeded();
      await page.screenshot({path:path.join(process.argv[4], 'structure-workbench.png')});
    }
    assert.deepEqual(errors, []);
    console.log('Verified structure downloads, original and transfer forecasts, study tables, JSON import integrity and mobile layout at ' + base);
  } finally { await browser.close(); }
})().catch(e => {console.error(e); process.exitCode = 1;});
