import {curve, evaluate} from './model-core.js';
import {drawForecast, drawYield, dyeColor} from './charts.js';

const $ = id => document.getElementById(id);
const fmt = (number, decimals=1) => Number(number).toLocaleString('en-US', {maximumFractionDigits:decimals});
const pair = ['JF630b','JF635b'];
const titles = {compare:'Fluorophore comparison',pilot:'Pilot recording',methods:'Methods & data'};
const store = {page:null,config:null,dyes:null,s21:null,epfl:null,gpu:null,frames:null,settings:{},pilot:null,
  result:null,points:[],hover:null,chartBounds:null,importNotice:'',traceDye:'JF630b'};

function element(tag, value='', className='') {
  const node = document.createElement(tag);
  if (value !== '') node.textContent = value;
  if (className) node.className = className;
  return node;
}
function numericCell(value) { return element('td', value, 'num'); }
function tableRow(cells) { const row=element('tr'); row.append(...cells); return row; }
function anchor(text, href, className='') { const link=element('a',text,className); link.href=href; return link; }
function paragraphs(target, sentences) { target.replaceChildren(...sentences.map(s=>element('p',s))); }
function dataLink(text, href) { return anchor(text,href,'inline-link'); }
async function readJSON(name) {
  const response = await fetch('data/' + name + '.json');
  if (!response.ok) throw new Error('Could not load ' + name + ' (' + response.status + ').');
  return response.json();
}
function saveJSON(name, value) {
  const blob=new Blob([JSON.stringify(value,null,2)+'\n'],{type:'application/json'});
  const url=URL.createObjectURL(blob), link=anchor('',url);
  link.download=name; document.body.append(link); link.click(); link.remove();
  setTimeout(()=>URL.revokeObjectURL(url),1000);
}
function navigate(scroll) {
  const name=location.hash.slice(1);
  store.page=Object.hasOwn(titles,name) ? name : 'compare';
  for (const key of Object.keys(titles)) $('view-'+key).hidden = key !== store.page;
  document.querySelectorAll('[data-nav]').forEach(link=>{
    if (link.dataset.nav === store.page) link.setAttribute('aria-current','page'); else link.removeAttribute('aria-current');
  });
  document.title = 'Blink Atlas · '+titles[store.page];
  if(scroll) window.scrollTo(0,0);
  requestAnimationFrame(redrawVisible);
}

function setFields() {
  const s=store.settings;
  $('density-min').value=s.densityMin; $('density-max').value=s.densityMax;
  $('threshold').value=s.minPhotons; $('condition').value=s.sampleState==='fixed' ? 'fixed' : 'live';
  $('candidate-panel').value=s.candidatePanel;
  $('threshold').max=store.config.limits.minPhotons[1];
  for (const id of ['density-min','density-max']) {
    $(id).min=store.config.limits.density[0]; $(id).max=store.config.limits.density[1];
  }
  const fixed = s.sampleState==='fixed';
  $('protocol-note').textContent=fmt(s.excitationNm)+' nm excitation · '+fmt(s.exposureMs)+' ms exposure · '+fmt(s.frameIntervalMs)+' ms cycle';
  $('protocol-note').title=fixed ? 'Reference protocol from Holland et al.' : 'External S21 acquisition';
}
function configureSlider() {
  const s=store.settings, limits=store.config.limits.density, slider=$('density-slider');
  const valid=s.densityMin>=limits[0] && s.densityMax<=limits[1] && s.densityMin<=s.densityMax;
  slider.disabled=!valid;
  if(valid) {
    const half=Math.log10(s.densityMax/s.densityMin)/2;
    slider.min=Math.log10(limits[0])+half; slider.max=Math.log10(limits[1])-half;
    slider.value=(Math.log10(s.densityMin)+Math.log10(s.densityMax))/2;
    slider.setAttribute('aria-valuetext',fmt(s.densityMin)+' to '+fmt(s.densityMax)+' labels per square micrometre');
  }
}
function updateComparison() {
  store.result=evaluate(store.dyes,store.config,store.settings,store.pilot);
  const r=store.result, s=store.settings;
  const d=$('decision'); d.dataset.status=r.status;
  const labels={stable:'Two candidates to test',uncertain:'Several candidates remain plausible',low_signal:'Photon threshold leaves too few candidates',
    out_of_scope:'Recommendation withheld',invalid:'Check the acquisition settings'};
  $('decision-title').textContent=r.status==='stable' ? r.recommendation.join(' + ') : labels[r.status];
  const panelSize=s.candidatePanel==='validation'?' · VALIDATION PANEL':'';
  $('decision-label').textContent=(r.status==='stable' ? 'MODEL-SELECTED PAIR' : 'CANDIDATE SELECTION')+panelSize;
  const contenderText=r.status==='uncertain' && r.contenders.length ? ' Plausible choices include '+r.contenders.join(', ')+'.' : '';
  $('decision-copy').textContent=r.reasons.join(' ')+contenderText;
  $('decision-next').textContent=r.nextMeasurement;
  const priorText=$('prior-comparison');
  priorText.hidden=r.pilot.status!=='applied';
  if(!priorText.hidden){
    const prior=evaluate(store.dyes,store.config,s);
    const describe=record=>record.status==='stable'?record.recommendation.join(' + '):'pair unresolved';
    priorText.textContent='Same settings, publication only: '+describe(prior)+' at '+prior.density.map(v=>fmt(v)).join('–')+
      ' labels/µm². With pilot: '+describe(r)+' at '+r.density.map(v=>fmt(v)).join('–')+' labels/µm².';
  }
  $('remove-pilot').hidden=!store.pilot;
  $('import-status').textContent=store.importNotice+(store.pilot && !store.importNotice ? r.pilot.reason : '');
  if(r.pilot.status==='applied') {
    $('import-status').textContent=r.pilot.reason+' Effective range: '+r.density.map(x=>fmt(x)).join('–')+' labels/µm².';
  }
  const tbody=$('candidate-rows'); tbody.replaceChildren();
  for(const row of r.ranking) {
    const dye=element('td'), dot=element('i','','dye-dot');
    dot.style.backgroundColor=dyeColor[row.id] || '#829e8a'; dot.style.border='1px solid #0002';
    dye.append(dot,document.createTextNode(row.name));
    tbody.append(tableRow([dye,numericCell(fmt(row.rate,2)),numericCell(fmt(row.low,2)+' – '+fmt(row.high,2)),
      numericCell(fmt(row.overlapProbability*100,1)+'%')]));
  }
  if(!r.ranking.length) {
    const td=element('td','No yield estimates for this input.','error-row'); td.colSpan=4; tbody.append(tableRow([td]));
  }
  $('center-density').textContent=r.ranking.length ? '@ '+fmt(r.centerDensity,1)+' labels/µm²' : '';
  document.querySelectorAll('[data-preset]').forEach(button=>{
    const [lo,hi]=button.dataset.preset==='broad' ? [100,1000] : [10000,10500];
    button.setAttribute('aria-pressed',String(s.densityMin===lo && s.densityMax===hi));
  });
  const minPower=s.densityMin>0 && s.densityMin<10 ? 0 : 1, maxPower=5;
  const densities=Array.from({length:181},(_,i)=>10**(minPower+(maxPower-minPower)*i/180));
  const study=s.candidatePanel==='validation';
  $('chart-mode').options[1].textContent=study?'All three dyes':'All nine dyes';
  $('other-legend').querySelectorAll('[data-dye]').forEach(item=>{item.hidden=study && !store.config.validationDyes.includes(item.dataset.dye);});
  store.points=r.ranking.length ? curve(store.dyes,store.config,s,densities) : [];
  store.hover=null; $('chart-tooltip').hidden=true;
  const message=$('plot-message'); message.hidden=!!store.points.length;
  message.textContent=message.hidden?'':r.reasons.join(' ');
  configureSlider(); redrawVisible();
}
function redrawVisible() {
  if(store.page==='compare' && store.result) {
    const all=$('chart-mode').value==='all';
    $('other-legend').hidden=!all;
    if(store.points.length) {
      store.chartBounds=drawYield($('yield-chart'),store.points,store.result,all,store.hover);
    } else {
      const canvas=$('yield-chart'); canvas.getContext('2d').clearRect(0,0,canvas.width,canvas.height);
      store.chartBounds=null;
    }
  }
  if(store.page==='pilot' && store.s21) drawForecast($('forecast-chart'),store.s21.evaluation.blocks);
}
function evidenceTable(id) {
  const dye=store.dyes.dyes.find(d=>d.id===id);
  if(!dye)return;
  store.traceDye=dye.id;
  const trace=$('trace'); trace.replaceChildren(element('h3',dye.name+' · '+dye.ligand+' · three cellular imaging experiments'));
  for(const [label,key] of [['ON duration','onMs'],['OFF duration','offMs'],['Photon emission rate','photonsPerMs'],['Blinks per molecule','blinksPerMolecule'],['ON:OFF ratio','onOffRatio']]) {
    const m=dye.metrics[key];
    const cells=m.sheet+' · '+m.cells+(m.medianCells?' / '+m.medianCells:'');
    const summary=m.replicates ? m.replicates.map(v=>fmt(v,3)).join(', ') : 'mean '+m.mean.toPrecision(4)+' · SD '+m.sd.toPrecision(3);
    trace.append(element('p',label+' ('+m.unit+')   '+summary+'   |   '+cells));
  }
}
function renderReference() {
  const legend=$('other-legend');
  for(const d of store.dyes.dyes.filter(d=>!pair.includes(d.id))){
    const item=element('span'), swatch=element('i','','line');item.dataset.dye=d.id;swatch.style.color=dyeColor[d.id];
    item.append(swatch,document.createTextNode(d.name));legend.append(item);
  }
  const refs=store.dyes.dyes.filter(d=>pair.includes(d.id));
  const ordered=pair.map(id=>refs.find(d=>d.id===id));
  const group=element('div','','ref-framegrid');
  for(const id of pair) {
    const figure=element('figure'), img=new Image(); img.src=store.frames.crops[id].image;
    img.alt='Published fixed COS-7 H2B–HaloTag field labeled with '+id;
    figure.append(img,element('figcaption',id)); group.append(figure);
  }
  const figure=element('figure','','ref-fig'), caption=element('figcaption');
  caption.append('First-frame crops · author display intensities. ',anchor('Holland et al., Video 1 (CC BY 4.0) ↗',store.frames.source));
  figure.append(group,caption); $('ref-images').replaceChildren(figure);
  for(const [label,key,scale,dp] of [['ON time (ms)','onMs',1,2],['OFF time (s)','offMs',1000,1],['Photons per ms','photonsPerMs',1,1]]) {
    $('pair-rows').append(tableRow([element('td',label),...ordered.map(d=>numericCell(fmt(d.metrics[key].mean/scale,dp)))]));
  }
}
function renderMethods() {
  for(const d of store.dyes.dyes) {
    const btn=element('button',d.name,'trace-toggle'); btn.type='button';
    btn.addEventListener('click',()=>{evidenceTable(d.id);$('trace').scrollIntoView({block:'nearest'});});
    const cell=element('td'); cell.append(btn);
    const m=d.metrics;
    $('measured-rows').append(tableRow([cell,numericCell(fmt(m.onMs.mean,2)),numericCell(fmt(m.onMs.median,2)),
      numericCell(fmt(m.offMs.mean/1000,2)),numericCell(fmt(m.offMs.median/1000,2)),numericCell(fmt(m.photonsPerMs.mean,2))]));
  }
  evidenceTable(store.traceDye);
  paragraphs($('method-description'),[
    'Each candidate must emit enough expected photons in one ON segment and have no active neighbor within 250 nm during the exposure. The model assumes stationary, independent blinking and locally uniform labels. Shot noise, long-run bleaching and camera background can alter practical detectability.',
    'We cross the three measured replicates for each input and shorten reported OFF times by a factor of two in a sensitivity scenario. A recommendation requires the same unordered pair to outrank every alternative throughout the selected density interval.',
    'A matched reference pilot updates a log-uniform density prior using counts of empty, single- and multiple-detection tiles. One capped effective field limits repeated counting of correlated frames. A mismatched recording or an incompatible occupancy pattern withholds the recommendation.'
  ]);
  const ep=store.epfl.at_demo_radius;
  paragraphs($('epfl-results'),[
    'Across '+fmt(store.epfl.frames)+' simulated tubulin frames, '+fmt(store.epfl.emitters)+' emitter coordinates were available. At 250 nm, '+fmt(ep.overlap_fraction*100,2)+'% had a neighbor; the whole-field uniform-density calculation gave '+fmt(ep.uniform_poisson_overlap_fraction*100,2)+'%.',
    'Spatial-bin and exhaustive distance checks agreed for every emitter. The large uniform-model discrepancy shows why a clustered specimen needs local density measurements.'
  ]);
  if(store.gpu) $('epfl-results').append(element('p','An independent CUDA exposure and neighbor simulation also checked the analytical probabilities for JF630b and JF635b.'),dataLink('Download the GPU verification ↓','data/gpu-check.json'));
}
function renderPilot() {
  const {acquisition:a,pilot:p,evaluation:e,media:m}=store.s21;
  const img=new Image(); img.src=m.poster; img.alt='Fixed-contrast crop of the live H2B field with cyan rings at author-fitted localizations.';
  $('pilot-media').replaceChildren(img);
  const controls=element('div','','media-controls'), play=element('button','Play 40-frame clip ▷'); play.type='button';
  let playing=false;
  play.addEventListener('click',()=>{playing=!playing;img.src=playing?m.preview:m.poster;play.textContent=playing?'Show first frame ▪':'Play 40-frame clip ▷';play.setAttribute('aria-pressed',String(playing));});
  controls.append(play,anchor('Open full field ↗',m.field));
  $('pilot-media').after(controls);
  $('media-caption').textContent='40 frames · 224 × 224 px crop · 2 µm scale bar. Rings show author-fitted positions. Walterspiel et al., CC BY 4.0.';
  const meta=[['Sample','Live '+a.cell_line+' · '+a.target],['Tag / reference dye',a.tag+' / '+a.reference_dye],['Excitation',fmt(a.excitation_nm)+' nm'],['Exposure / cycle',fmt(a.exposure_ms,3)+' / '+fmt(a.frame_interval_ms,3)+' ms'],['Camera field',a.width_px+' × '+a.height_px+' px'],['Pixel size',fmt(a.pixel_size_nm)+' nm'],['Recorded length',fmt(a.total_frames)+' frames']];
  for(const [key,value] of meta) {const row=element('div');row.append(element('dt',key),element('dd',value));$('acquisition-rows').append(row);}
  for(const [label,value] of [['Pilot localizations',fmt(p.localizations)],['Trend MAE / block',fmt(e.model_mae,2)],['Baseline MAE / block',fmt(e.baseline_mae,2)]]) {
    const box=element('div');box.append(element('small',label),element('strong',value));$('pilot-metrics').append(box);
  }
  const relative=100*(1-e.model_mae/e.baseline_mae);
  $('pilot-verdict').textContent=e.model_improved ? 'The damped trend reduced mean absolute block error by '+fmt(relative,1)+'% in this recording. It still overestimated most later blocks as detection counts declined.' : 'The constant pilot mean performed better on this recording. The trend is retained for quality review.';
  for(const b of e.blocks) $('forecast-rows').append(tableRow([element('td',(b.start_frame+1)+'–'+(b.start_frame+b.frames)),
    ...[b.observed,b.prediction,b.baseline,Math.abs(b.observed-b.prediction),Math.abs(b.observed-b.baseline)].map(n=>numericCell(fmt(n,2)))]));
}

function showHover(index) {
  if(!store.chartBounds || !store.points.length) return;
  store.hover=Math.max(0,Math.min(store.points.length-1,index));
  const p=store.points[store.hover], tip=$('chart-tooltip');
  const series=$('chart-mode').value==='all'?Object.keys(p.rates).sort((a,b)=>p.rates[b]-p.rates[a]):pair;
  tip.textContent=fmt(p.density)+' labels/µm²\n'+series.map(id=>id+'   '+fmt(p.rates[id],2)).join('\n');
  tip.hidden=false;
  const {left,right}=store.chartBounds, box=$('yield-container');
  const xx=left+(right-left)*store.hover/(store.points.length-1);
  tip.style.left=Math.max(0,Math.min(box.clientWidth-tip.offsetWidth,xx+9))+'px';tip.style.top='1px';
  redrawVisible();
}
function hookEvents() {
  window.addEventListener('hashchange',()=>navigate(true));
  new ResizeObserver(redrawVisible).observe($('yield-container'));
  new ResizeObserver(redrawVisible).observe($('forecast-chart').parentElement);
  document.fonts.ready.then(redrawVisible);
  for(const [id,key] of [['density-min','densityMin'],['density-max','densityMax'],['threshold','minPhotons']]) {
    $(id).addEventListener('input',()=>{store.settings[key]=$(id).valueAsNumber;updateComparison();});
  }
  $('density-slider').addEventListener('input',()=>{
    const half=Math.log10(store.settings.densityMax/store.settings.densityMin)/2, center=Number($('density-slider').value);
    const limits=store.config.limits.density;
    store.settings.densityMin=Math.max(limits[0],Math.round(10**(center-half)));
    store.settings.densityMax=Math.min(limits[1],Math.round(10**(center+half)));
    setFields();updateComparison();
  });
  for(const button of document.querySelectorAll('[data-preset]'))button.addEventListener('click',()=>{
    const high=button.dataset.preset==='dense';store.settings.densityMin=high?10000:100;store.settings.densityMax=high?10500:1000;setFields();updateComparison();
  });
  $('condition').addEventListener('change',()=>{
    const {defaults}=store.config, a=store.s21.acquisition;
    const keep={densityMin:store.settings.densityMin,densityMax:store.settings.densityMax,minPhotons:store.settings.minPhotons,candidatePanel:store.settings.candidatePanel};
    store.settings=$('condition').value==='fixed'?{...defaults,...keep}:{...defaults,...keep,sampleState:'live',cellLine:a.cell_line,
      exposureMs:a.exposure_ms,excitationNm:a.excitation_nm,frameIntervalMs:a.frame_interval_ms};
    setFields(); updateComparison();
  });
  $('reset').addEventListener('click',()=>{store.settings={...store.config.defaults};store.pilot=null;store.importNotice='';$('chart-mode').value='reference';setFields();updateComparison();});
  $('chart-mode').addEventListener('change',redrawVisible);
  $('candidate-panel').addEventListener('change',()=>{store.settings.candidatePanel=$('candidate-panel').value;updateComparison();});
  const canvas=$('yield-chart');
  canvas.addEventListener('pointermove',e=>{if(!store.chartBounds)return;const {left,right}=store.chartBounds;showHover(Math.round((e.clientX-canvas.getBoundingClientRect().left-left)/(right-left)*(store.points.length-1)));});
  function closeTip(){store.hover=null;$('chart-tooltip').hidden=true;redrawVisible();}
  canvas.addEventListener('pointerleave',closeTip);canvas.addEventListener('blur',closeTip);
  canvas.addEventListener('keydown',e=>{if(['ArrowLeft','ArrowRight'].includes(e.key)){e.preventDefault();showHover((store.hover??90)+(e.key==='ArrowRight'?3:-3));}if(e.key==='Escape')closeTip();});
  $('export').addEventListener('click',()=>saveJSON('blink-atlas-analysis.json',{schema:'blink-atlas-analysis/1',model:store.config.schema,
    measurements:store.dyes.schema,measurementDoi:'10.1038/s41592-026-03062-5',settings:store.settings,pilot:store.pilot,
    chartMode:$('chart-mode').value,result:store.result}));
  $('export-pilot').addEventListener('click',()=>saveJSON('s21-pilot.json',store.s21));
  for(const link of document.querySelectorAll('[data-watch]'))link.addEventListener('click',event=>{
    event.preventDefault(); const video=$('guide-video');
    if(!video.getAttribute('src'))video.src=link.href;
    $('guide').showModal();video.play().catch(()=>{});
  });
  $('close-guide').addEventListener('click',()=>$('guide').close());
  $('guide').addEventListener('close',()=>$('guide-video').pause());
  $('import').addEventListener('click',()=>$('import-file').click());
  $('remove-pilot').addEventListener('click',()=>{store.pilot=null;store.importNotice='';updateComparison();});
  $('import-file').addEventListener('change',async()=>{
    const file=$('import-file').files[0];if(!file)return;
    try{
      if(file.size>8*1024*1024)throw new Error('Choose a summary JSON under 8 MB.');
      const value=JSON.parse(await file.text());
      if(!value || typeof value!=='object')throw new Error('Expected a pilot summary or an exported analysis.');
      if(value.schema==='smlm-pilot/1') {
        store.pilot=value;store.importNotice='';
      } else if(value.schema==='blink-atlas-analysis/1' && value.model===store.config.schema && value.measurements===store.dyes.schema && value.settings) {
        const next={};for(const key of Object.keys(store.config.defaults)){
          if(key==='candidatePanel' && value.settings[key]===undefined){next[key]='all';continue;}
          if(typeof value.settings[key]!==typeof store.config.defaults[key])throw new Error('Saved acquisition has an invalid '+key+' value.');
          next[key]=value.settings[key];
        }
        if(value.pilot!=null && value.pilot.schema!=='smlm-pilot/1')throw new Error('Invalid embedded pilot.');
        store.settings=next;store.pilot=value.pilot??null;store.importNotice='Analysis restored and recomputed.';
        $('chart-mode').value=value.chartMode==='all'?'all':'reference';setFields();
      } else {throw new Error('Choose a compatible smlm-pilot/1 JSON or Blink Atlas export.');}
      updateComparison();
    }catch(error){store.importNotice='Import failed. '+error.message;$('import-status').textContent=store.importNotice;}
    finally{$('import-file').value='';}
  });
}
async function start() {
  const [dyes,config,s21,epfl,frames,gpu]=await Promise.all(['dyes','model','pilot','epfl','reference-frames','gpu-check'].map(readJSON));
  Object.assign(store,{dyes,config,s21,epfl,frames,gpu,settings:{...config.defaults}});
  renderReference();renderMethods();renderPilot();setFields();hookEvents();
  $('loading').hidden=true;navigate(false);updateComparison();
}
start().catch(error=>{$('loading').hidden=false;$('loading').textContent='The data could not be opened. Reload the site or download the analysis from Source code. '+error.message;});
