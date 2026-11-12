// Verify the real mounted Pages route, provenance files and anonymous access.
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const path = require('node:path');
const {chromium} = require('playwright');

const base = process.argv[2] || 'http://127.0.0.1:8109/';
const channel = process.argv[3] || 'chrome';
const fmt = v => Number.isFinite(v) ? v.toLocaleString('en-US',{maximumFractionDigits:2}) : '—';

(async()=>{
  const browser=await chromium.launch({headless:true,channel});
  try {
    const ctx=await browser.newContext({viewport:{width:1366,height:900}}), page=await ctx.newPage();
    const errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    page.on('response',r=>{if(r.url().startsWith(base) && r.status()>=400) errors.push(r.status()+' '+r.url());});
    const objects={};
    for (const name of ['palette.json','palette-raw.json','palette-manifest.json','palette-audit.json','palette-evaluation.json',
                        'steen-figure2.xlsx','steen-figure3.xlsx']) {
      const r=await ctx.request.get(new URL('data/'+name,base).href);
      assert.equal(r.status(),200,name);
      const content=await r.body();
      assert.equal(Buffer.compare(content,await fs.readFile(path.join(__dirname,'..','docs','data',name))),0,name);
      if(name.endsWith('json')) objects[name]=JSON.parse(content);
    }
    const catalog=objects['palette.json'], raw=objects['palette-raw.json'], evaluation=objects['palette-evaluation.json'];
    await page.goto(base+'#dna-paint');
    await page.locator('#palette-gate').waitFor();
    assert.equal(await page.locator('#palette-gate').getAttribute('data-state'),evaluation.raw.gatePassed?'passed':'held');
    assert.match(await page.locator('#palette-coverage').innerText(),/109 of 109/);
    for (const specimen of ['origami','npc']) {
      await page.locator('#palette-system').selectOption(specimen);
      await page.locator('#palette-metric').selectOption('nenaNm');
      for (const wavelength of [488,560,640]) {
        await page.locator('#palette-wave').selectOption(String(wavelength));
        const panel=catalog.dyes.filter(d=>d.wavelength===wavelength);
        assert.equal(await page.locator('.palette-series-row').count(),panel.length);
        for (const dye of panel) {
          const source=dye.systems[specimen].metrics.nenaNm;
          const avg=raw.dyes.find(d=>d.dye===dye.id && d.system===specimen).metrics.nenaNm.mean;
          const row=page.locator('.palette-series-row[data-dye="'+dye.id+'"] .palette-number');
          assert.deepEqual(await row.allTextContents(),[fmt(source.mean),fmt(avg)],specimen+'/'+dye.id);
        }
      }
    }
    await page.locator('#palette-dye').selectOption('janeliafluor646');
    await page.locator('#palette-metric').selectOption('centerPhotons');
    await page.locator('#palette-agreement-summary').click();
    assert.match(await page.locator('#palette-agreement-body').innerText(),/workbook repeat values/);
    await page.locator('#palette-eval-source').selectOption('published');
    assert.equal(await page.locator('#palette-gate').getAttribute('data-state'),evaluation.raw.gatePassed?'passed':'held');
    for (const width of [320,390,760,851,941,1051]) {
      await page.setViewportSize({width,height:844});
      assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),'Viewport '+width);
    }
    await page.locator('[data-nav=precision]').click();
    await page.locator('#fixed-comparison').waitFor();
    await page.locator('[data-nav=dna-paint]').click();
    await page.locator('#palette-dye').waitFor();
    assert.deepEqual(errors,[]);
    console.log('Verified deployed file bytes, 36 dye/specimen NeNA comparisons, source audit, gate, navigation and responsive views.');
  } finally { await browser.close(); }
})().catch(e=>{console.error(e);process.exitCode=1;});
