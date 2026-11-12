// Capture the actual site in Chromium, including input changes and the public movie.
// Usage: node tools/capture-tour.cjs URL OUTPUT.webm [chrome]
const {chromium} = require('playwright');
const fs = require('node:fs');
const path = require('node:path');
const {execFileSync} = require('node:child_process');
const [url, output, channel] = process.argv.slice(2);
if (!url || !output) throw Error('Usage: node tools/capture-tour.cjs URL OUTPUT.webm [chrome]');

(async () => {
  const browser = await chromium.launch({headless:true, ...(channel ? {channel} : {})});
  try {
    const dir=path.dirname(output);fs.mkdirSync(dir,{recursive:true});
    const proposalImage=path.join(dir,'walkthrough-proposal');
    execFileSync('pdftoppm',['-png','-singlefile','-scale-to','1350',
      path.join(__dirname,'../docs/media/validation-proposal.pdf'),proposalImage],{stdio:'ignore'});
    const proposal=fs.readFileSync(proposalImage+'.png').toString('base64');
    const context=await browser.newContext({viewport:{width:1280,height:800},deviceScaleFactor:1,
      recordVideo:{dir,size:{width:1280,height:800}}});
    const page=await context.newPage();
    await page.goto(url);await page.locator('#candidate-rows tr').first().waitFor();
    await page.evaluate(()=>document.fonts.ready);
    await page.addStyleTag({content:
      'html{scroll-behavior:smooth}#record-pointer{position:fixed;z-index:99999;pointer-events:none;width:17px;height:17px;'+
      'border:2px solid #183c30;outline:2px solid #deffbe;border-radius:100%;box-shadow:0 0 0 5px #96c97844;transition:left .25s,top .25s;left:50%;top:40%}'});
    await page.evaluate(()=>{const cursor=document.createElement('div');cursor.id='record-pointer';document.body.append(cursor);
      window.addEventListener('pointermove',e=>{cursor.style.left=e.clientX+'px';cursor.style.top=e.clientY+'px';});});
    const start=performance.now();
    const at=async(s)=>page.waitForTimeout(Math.max(0,s*1000-(performance.now()-start)));
    async function focus(selector) {
      await page.locator(selector).first().evaluate(el=>el.scrollIntoView({behavior:'smooth',block:'center'}));
      await page.waitForTimeout(650);await page.locator(selector).first().hover();
    }
    async function route(name) {await page.evaluate(v=>{location.hash=v;},name);await page.waitForTimeout(600);}
    await at(5);await page.locator('#yield-chart').hover({position:{x:410,y:90}});
    await at(11);await page.locator('#plot-message').evaluate(()=>window.scrollTo({top:0,behavior:'smooth'}));
    await at(14);await focus('.reference-card');
    await at(22);await page.locator('#pair-rows').hover();
    await at(30);await focus('[data-preset=broad]');
    await at(35);await page.locator('[data-preset=dense]').click();
    await at(38);await focus('#decision');
    await at(44);await focus('#threshold');
    await at(47);await page.locator('#threshold').fill('2200');await focus('#decision');
    await at(53);await page.locator('#condition').selectOption('live');await focus('#decision');
    await at(58);await page.locator('#reset').click();
    await at(60);await route('pilot');await page.locator('.media-controls button').click();
    await at(68);await page.locator('#acquisition-rows').hover();
    await at(76);await focus('.holdout');await page.locator('#pilot-metrics').hover();
    await at(85);await page.locator('#forecast-chart').hover({position:{x:850,y:150}});
    await at(90);await page.locator('.holdout summary').click();
    await at(93);await route('methods');await page.getByRole('button',{name:'JF630b',exact:true}).click();await focus('#trace');
    await at(100);await focus('#epfl-results');
    await at(107);await focus('.evidence');
    await at(108);await page.getByRole('link',{name:/Joint validation proposal/}).hover();
    await at(109);
    await page.setContent('<html lang="en"><title>Joint validation proposal</title><body style="margin:0;background:#26342e;text-align:center">'+
      '<img alt="One-page joint validation proposal" style="width:830px;max-width:100%;display:block;margin:0 auto" src="data:image/png;base64,'+
      proposal+'"></body></html>');
    await at(114);await page.evaluate(()=>window.scrollTo({top:document.body.scrollHeight,behavior:'smooth'}));
    await at(120.8);
    const video=page.video();await context.close();await video.saveAs(output);
    if(path.resolve(await video.path())!==path.resolve(output))await video.delete();
  } finally {await browser.close();}
})().catch(err=>{console.error(err);process.exitCode=1;});
