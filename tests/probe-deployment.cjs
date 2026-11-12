// Exercise the researcher's workflow in a fresh browser at localhost or on Pages.
const assert = require('node:assert/strict');
const {createHash} = require('node:crypto');
const fs = require('node:fs/promises');
const path = require('node:path');
const {chromium} = require('playwright');

const base = process.argv[2] || 'http://127.0.0.1:8187/';
const channel = process.argv[3] || 'chrome';
const screenshotDir = process.argv[4];
const repo = path.join(__dirname, '..', 'docs');
const hash = body => createHash('sha256').update(body).digest('hex');

async function verifyUnseen(page, bytes) {
  const result = JSON.parse(await bytes('data/analogue-benchmark.json'));
  const workbook = await bytes('data/analogue-observations.csv');
  const studies = await bytes('data/analogue-studies.json');
  assert.equal(hash(workbook), result.sourceHashes.observations);
  assert.equal(hash(studies), result.sourceHashes.studies);
  const item = key => page.locator('#unseen-' + key);
  await page.setViewportSize({width:1366,height:900});
  await page.locator('#view-probe-design .title-actions a[href="#unseen-dyes"]').click();
  await item('chart').locator('.unseen-bar').first().waitFor();
  assert.match(await item('counts').innerText(), /11 additional parent dyes, 25/);
  assert.match(await item('gate').innerText(), /validation is still open/);
  assert.match(await item('conjugate-effect').innerText(), /16.5×/);
  assert.equal(await item('conjugates').locator('tr').count(), 3);
  assert.match(await item('titration-effect').innerText(), /at 1 of 4 laser powers/);
  assert.equal(await item('titrations').locator('tr').count(), 4);
  assert.equal(await item('training').inputValue(), 'remmel-igg-frozen');
  assert.equal(await page.locator('[data-nav=probe-design]').getAttribute('aria-current'), 'location');
  if (screenshotDir) await page.screenshot({path:path.join(screenshotDir, 'unseen-desktop.png')});
  for (const set of result.sets) {
    const group = set.study === 'holland' ? 'holland' : set.id.replace(/-(frozen|local)$/, '');
    await item('assay').selectOption(group);
    await item('training').selectOption(set.id);
    const s = set.summary;
    assert.equal(await item('coverage').textContent(), s.scored + ' of ' + s.total + ' dyes scored', set.id);
    assert.equal(await item('chart').locator('.unseen-bar').count(), s.scored);
    const table = item('rows').locator('tr');
    assert.equal(await table.count(), Math.max(1, s.total));
    for (const [i, r] of set.rows.entries()) {
      const cells = table.nth(i).locator('td');
      assert.equal(await cells.nth(0).innerText(), r.dye);
      if (r.status === 'scored') {
        assert.equal(await cells.nth(3).innerText(), r.predictions.local.toExponential(3));
        await cells.last().locator('summary').click();
        assert.ok((await cells.last().innerText()).includes(r.trainingIds.join(', ')));
      } else {
        assert.ok((await cells.last().innerText()).includes(r.status));
      }
    }
    assert.deepEqual(await item('metrics').locator('strong').allTextContents(),
      ['local','nearest','regression'].map(k => s.mae ? s.mae[k].toFixed(3) : '—'));
  }
  await item('assay').selectOption('remmel-igg');
  for (const width of [320,390,760,850,851,940,1060,1191]) {
    await page.setViewportSize({width,height:844});
    assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'Unseen-dye viewport ' + width);
  }
  if (screenshotDir) {
    await page.setViewportSize({width:390,height:844});
    await page.evaluate(() => scrollTo(0,0));
    await page.screenshot({path:path.join(screenshotDir, 'unseen-mobile.png')});
  }
  await page.goto(base + '#unseen-dyes');
  await page.reload();
  await item('chart').locator('.unseen-bar').first().waitFor();
  assert.ok((await item('freeze').getAttribute('href')).includes(result.frozenRevision));
  await page.locator('#view-unseen-dyes .page-title a[href="#probe-design"]').click();
}

(async () => {
  const browser = await chromium.launch({headless:true, channel});
  try {
    const ctx = await browser.newContext({viewport:{width:1366,height:900},acceptDownloads:true});
    const page = await ctx.newPage(), errors = [];
    page.on('pageerror', err => errors.push(err.message));
    page.on('response', r => { if (r.url().startsWith(base) && r.status() >= 400) errors.push(r.url()); });
    const bytes = async name => {
      const res = await ctx.request.get(new URL(name, base).href);
      assert.equal(res.status(), 200, name);
      const buf = await res.body();
      assert.equal(hash(buf), hash(await fs.readFile(path.join(repo, name))), name);
      return buf;
    };
    const data = JSON.parse(await bytes('data/probe-design.json'));
    for (const source of Object.values(data.sources)) assert.equal(hash(await bytes(source.file)), source.sha256);
    for (const dye of data.dyes) await bytes(dye.structure.image);
    await page.goto(base + '#probe-design');
    const item = name => page.locator('#probe-' + name);
    await item('outcome').waitFor();
    assert.equal(await item('outcome').getAttribute('data-order-change'), 'true');
    assert.match(await item('outcome').innerText(), /\+12\.14 percentage points live and -36\.88/);
    assert.equal(await item('a').locator('option').count(), 9);
    assert.match(await item('benchmarks').innerText(), /0\.216/);
    if (screenshotDir) {
      await fs.mkdir(screenshotDir, {recursive:true});
      await page.evaluate(() => document.fonts.ready);
      await page.screenshot({path:path.join(screenshotDir, 'workbench.png')});
    }
    for (const d of data.dyes) {
      await item('a').selectOption(d.id);
      assert.notEqual(await item('a').inputValue(), await item('b').inputValue());
      assert.equal(await item('molecules').locator('img').first().getAttribute('src'), d.structure.image);
      for (const m of data.contexts.er.metrics) {
        await item('metric').selectOption(m);
        const rows = item('context-table').locator('tr');
        for (const [index, state] of ['live','fixed'].entries()) {
          assert.equal(Number(await rows.nth(index*2).locator('td').nth(2).textContent()), data.contexts.er.dyes[d.id][state][m].mean);
        }
      }
    }
    await item('context').selectOption('lysosome');
    assert.equal(await item('a').locator('option').count(), 2);
    assert.equal(await item('metric').locator('option').count(), 2);
    assert.match(await item('lysosome-note').innerText(), /different labeling concentrations/);
    assert.equal(await item('outcome').getAttribute('data-order-change'), 'true');
    await item('ph').selectOption('1.14');
    assert.match(await item('ph-reading').innerText(), /JF639b\s+no reading/);
    await item('jump').click();
    assert.equal(await item('trial-heading').evaluate(e => e === document.activeElement), true);
    await item('compound').fill('Blinded sample 01');
    await item('pka').fill('0.9');
    await item('scope').selectOption(data.calibration.scope);
    await item('lock').click();
    assert.match(await item('trial-status').innerText(), /Outside calibration/);
    assert.equal(await item('score-form').isVisible(), false);
    await item('pka').fill('3.20');
    await item('scope').selectOption('outside');
    assert.match(await item('preview').innerText(), /matched fixed COS-7/);
    await item('scope').selectOption(data.calibration.scope);
    const download = async button => {
      const [d] = await Promise.all([page.waitForEvent('download'), button.click()]);
      return JSON.parse(await fs.readFile(await d.path(), 'utf8'));
    };
    const card = await download(item('lock'));
    assert.equal(card.plan.schema, 'blink-atlas/probe-trial/2');
    assert.equal(card.plan.candidate, 'bracket-logratio-v1');
    assert.ok(card.plan.ratios.local > 0);
    assert.equal(card.plan.nearestDye,'JF635b');
    assert.equal(card.measurement, null);
    assert.equal(await item('pka').isDisabled(),true);
    const importCard = async record => {
      await item('import-details').evaluate(el => {el.open = true;});
      await item('import').setInputFiles({name:'test.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(record))});
    };
    await page.reload();
    await item('compound').waitFor();
    await importCard(card);
    await item('score-form').waitFor();
    await item('score').click();
    assert.match(await item('score-result').innerText(), /finite numeric/);
    for (const [i,v] of [.0004,.0005,.0006].entries()) await item('rep-' + (i+1)).fill(String(v));
    await item('score').click();
    assert.match(await item('score-result').innerText(), /Nearest measured pKa/);
    const scored = await download(item('save-score'));
    assert.deepEqual(scored.measurement.closest, ['nearest']);
    assert.equal(scored.measurement.ratios.length, 3);
    assert.deepEqual((await download(item('save-plan'))).plan, card.plan);
    await item('clear').click();
    assert.equal(await item('pka').inputValue(), '');
    await importCard(scored);
    await item('save-score').waitFor();
    const bad = structuredClone(scored);
    bad.plan.ratios.nearest = 9;
    await importCard(bad);
    await page.waitForFunction(() => document.getElementById('probe-trial-status').textContent.includes('Could not open'));
    assert.equal(await item('save-score').isVisible(),true);
    await item('clear').click();
    for (const width of [320,390,760,850,851,940,941,1021,1051,1191]) {
      await page.setViewportSize({width,height:844});
      assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'Viewport ' + width);
    }
    if (screenshotDir) {
      await page.setViewportSize({width:390,height:844});
      await page.evaluate(() => scrollTo(0,0));
      await page.screenshot({path:path.join(screenshotDir, 'mobile.png')});
    }
    await page.locator('[data-nav=compare]').click();
    await page.locator('#candidate-rows tr').first().waitFor();
    await page.locator('[data-nav=probe-design]').click();
    await item('outcome').waitFor();
    await verifyUnseen(page, bytes);
    await item('outcome').waitFor();
    assert.deepEqual(errors, []);
    console.log('Verified workbench sources, prospective trial, all unseen-dye evaluations, measured conjugates, responsive layouts and public source hashes at ' + base);
  } finally { await browser.close(); }
})().catch(e => {console.error(e);process.exitCode=1;});
