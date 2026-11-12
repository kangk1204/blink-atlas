// Run against a preview or a freshly deployed GitHub Pages root.
const assert = require('node:assert/strict');
const {chromium} = require('playwright');
const url = process.argv[2] || 'http://127.0.0.1:8109/';

(async () => {
  const browser = await chromium.launch({headless:true, channel:process.argv[3] || 'chrome'});
  try {
    const ctx = await browser.newContext({viewport:{width:1320,height:860}}), page=await ctx.newPage(), errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    await page.goto(url+'#precision');
    await page.locator('#epfl-calibration').waitFor();
    const content=page.locator('#precision-content');
    assert.match(await content.innerText(),/88.51–89.48%/);
    assert.match(await content.innerText(),/110,694 of 214,794/);
    assert.match(await content.innerText(),/7.37%/);
    assert.match(await content.innerText(),/12 of 20/);
    await page.getByText('Inspect individual counts, yields and held-out errors').click();
    await page.locator('#qc-filter').selectOption('HaloTag');
    assert.equal(await page.getByRole('table',{name:'Individual quality checks'}).locator('tbody tr').count(),3);
    const data=await (await ctx.request.get(new URL('data/precision.json',url).href)).json();
    const json=(name,value)=>({name,mimeType:'application/json',buffer:Buffer.from(JSON.stringify(value))});
    const input=page.locator('#precision-file'), result=page.locator('#precision-import-result');
    await input.setInputFiles(json('observed.json',data.s21));
    await result.getByText(/214,794/).first().waitFor();
    assert.match(await result.innerText(),/27.3659/);
    await input.setInputFiles([json('one.json',data.s21),json('two.json',data.auxiliary.find(f=>f.acquisition.reference_dye==='JF635b'))]);
    await result.getByRole('heading',{name:'Matched comparison unavailable'}).waitFor();
    assert.match(await result.innerText(),/fixed COS-7/i);
    await input.setInputFiles(json('invalid.json',{}));
    await page.getByText(/Import failed/).waitFor();
    await page.setViewportSize({width:390,height:844});
    await page.waitForFunction(()=>document.documentElement.scrollWidth<=window.innerWidth);
    await page.locator('[data-nav=compare]').click();
    assert.equal(await page.locator('#candidate-rows tr').count(),9);
    await page.locator('[data-nav=precision]').click();
    assert.equal(await page.locator('#fixed-comparison').isVisible(),true);
    assert.deepEqual(errors,[]);
    console.log('Precision cutoff, calibration, archive filter, portable import, abstention and mobile routes passed.');
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
