// Run against a localhost preview or the deployed GitHub Pages URL.
// Usage: node tests/site-smoke.cjs URL [chrome]
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const {chromium} = require('playwright');
const url = process.argv[2] || 'http://localhost:8000/';
const channel = process.argv[3];

(async () => {
  const browser=await chromium.launch({headless:true,...(channel?{channel}:{})});
  try {
    const context=await browser.newContext({viewport:{width:1366,height:900},acceptDownloads:true});
    const page=await context.newPage(), errors=[];
    page.on('pageerror',error=>errors.push(error.message));
    page.on('response',response=>{
      if(response.url().startsWith(url) && response.status()>=400)errors.push(response.status()+' '+response.url());
    });
    const status=()=>page.locator('#decision').getAttribute('data-status');
    await page.goto(url);await page.locator('#candidate-rows tr').first().waitFor();
    assert.equal(await page.locator('#candidate-rows tr').count(),9);
    assert.equal(await status(),'uncertain');
    await page.locator('[data-preset=dense]').click();assert.equal(await status(),'stable');
    assert.equal(await page.locator('#decision-title').innerText(),'JF630b + JF629b');
    await page.locator('#chart-mode').selectOption('all');
    assert.equal(await page.locator('#other-legend span').count(),7);
    const [download]=await Promise.all([page.waitForEvent('download'),page.locator('#export').click()]);
    const saved=JSON.parse(await fs.readFile(await download.path(),'utf8'));
    assert.deepEqual(saved.result.recommendation,['JF630b','JF629b']);
    await page.locator('#threshold').fill('2200');assert.equal(await status(),'low_signal');
    await page.locator('#density-min').fill('0');assert.equal(await status(),'invalid');
    await page.locator('#reset').click();
    await page.locator('#condition').selectOption('live');assert.equal(await status(),'out_of_scope');
    assert.equal(await page.locator('#candidate-rows .error-row').count(),1);
    const importJSON=async value=>page.locator('#import-file').setInputFiles({
      name:'case.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(value))});
    await importJSON(saved);
    await page.waitForFunction(()=>document.getElementById('decision').dataset.status==='stable');
    assert.equal(await page.locator('#chart-mode').inputValue(),'all');
    const pilot=await (await context.request.get(new URL('data/pilot.json',url).href)).json();
    await importJSON(pilot);await page.waitForFunction(()=>document.getElementById('decision').dataset.status==='out_of_scope');
    await page.locator('#remove-pilot').click();assert.equal(await status(),'stable');
    await page.locator('#import-file').setInputFiles({name:'broken.json',mimeType:'application/json',buffer:Buffer.from('{')});
    await page.waitForFunction(()=>document.getElementById('import-status').textContent.includes('Import failed'));
    await page.locator('[data-nav=pilot]').click();await page.locator('#view-pilot').waitFor();
    assert.equal(await page.locator('#forecast-rows tr').count(),8);
    assert.match(await page.locator('#pilot-metrics').innerText(),/456.32/);
    await page.locator('.media-controls button').click();
    assert.match(await page.locator('#pilot-media img').getAttribute('src'),/preview.webp/);
    await page.locator('.media-controls button').click();
    await page.locator('[data-nav=methods]').click();await page.getByRole('button',{name:'JF635b',exact:true}).click();
    assert.match(await page.locator('#trace').innerText(),/JF635b/);
    const pdf=await context.request.get(new URL('media/validation-proposal.pdf',url).href);
    assert.equal(pdf.status(),200);assert.equal((await pdf.body()).subarray(0,5).toString(),'%PDF-');
    await page.locator('[data-nav=compare]').click();await page.locator('.title-actions [data-watch]').click();
    await page.waitForFunction(()=>document.querySelector('#guide-video').readyState>=2);
    const duration=await page.locator('#guide-video').evaluate(v=>v.duration);
    assert.ok(Math.abs(duration-120)<0.15,'Video duration: '+duration);
    await page.keyboard.press('Escape');assert.equal(await page.locator('#guide').isVisible(),false);
    await page.setViewportSize({width:390,height:844});
    for(const id of ['compare','pilot','methods']){
      await page.locator('[data-nav='+id+']').click();
      await page.waitForFunction(()=>document.documentElement.scrollWidth<=window.innerWidth);
    }
    assert.deepEqual(errors,[]);
    console.log('Verified comparison, abstentions, import/export, S21, source rows, PDF, 120 s video and mobile layouts at '+url);
    await context.close();
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
