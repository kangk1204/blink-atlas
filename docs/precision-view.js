import {compareFields, validateComparison, validateField} from './precision-compare.js';

const format = (n, places=2) => Number.isFinite(n) ? n.toLocaleString('en-US', {maximumFractionDigits:places}) : 'Unavailable';
const pct = n => Number.isFinite(n) ? format(n*100)+'%' : 'Unavailable';
const number = n => format(n, 0);
const repo = 'https://github.com/kangk1204/blink-atlas/blob/main/';

function node(tag, value, cls) {
  const el = document.createElement(tag);
  if (value !== undefined) el.textContent = value;
  if (cls) el.className = cls;
  return el;
}
function link(label, href) {
  const el = node('a', label, 'inline-link');
  el.href = href;
  return el;
}
function links(...items) {
  const bar = node('div', undefined, 'ev-links');
  bar.append(...items);
  return bar;
}
function saveButton(label, filename, value) {
  const button=node('button',label,'text-button');
  button.type='button';
  button.onclick=()=>{
    const blob=new Blob([JSON.stringify(value,null,2)+'\n'],{type:'application/json'});
    const url=URL.createObjectURL(blob), a=node('a');
    a.href=url; a.download=filename; document.body.append(a); a.click(); a.remove();
    setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  return button;
}
function card(kicker, title) {
  const el = node('article', undefined, 'card ev-card');
  el.append(node('div', kicker, 'eyebrow'), node('h2', title));
  return el;
}
function table(headers, rows, label) {
  const wrap = node('div', undefined, 'table-wrap'), t = node('table'), head = t.createTHead().insertRow();
  t.setAttribute('aria-label', label);
  for (const h of headers) { const th = node('th', h); th.scope = 'col'; head.append(th); }
  const body = t.createTBody();
  for (const row of rows) {
    const tr = body.insertRow();
    for (const v of row) tr.append(node('td', v));
  }
  wrap.append(t);
  return wrap;
}
function bar(fraction, boundary=0) {
  const el = node('div', undefined, 'ev-bar');
  el.setAttribute('aria-hidden','true');
  for (const [name, val] of [['passed',fraction],['boundary',boundary]]) {
    const segment = node('span', undefined, name);
    segment.style.width = 100*val+'%';
    el.append(segment);
  }
  return el;
}
function publicResults(evidence) {
  const h = evidence.holland, f = evidence.s21, model = node('article', undefined, 'card ev-status');
  model.id = 'fixed-comparison';
  model.append(node('h2','Low-density JF639b versus JF630b remains unresolved'),
    node('p','The public sources examined contain no matched, density-calibrated fixed COS-7 fields for this pair. The two measurements below use different samples and uncertainty definitions. An observed yield comparison requires calibrated field areas, equal acquisition times and the same uncertainty definition.'),
    links(link('Inspect the source search ↗',repo+'validation/h2b_precision/README.md'),link('Open modeled yields →','#compare')));
  const columns = node('div', undefined, 'ev-columns');
  const fixed = card('FIXED COS-7 · HOLLAND ET AL.', 'JF630b histogram');
  const [low, high] = h.qualifying_fraction_bounds;
  fixed.append(node('strong',format(100*low)+'–'+pct(high),'ev-number'),
    node('p','of the displayed scalar precision histogram meets the 20 nm cutoff.'),
    bar(low, high-low),
    node('p',number(h.qualifying_count_bounds[0])+' counts fall entirely below the threshold. The striped interval represents '+number(h.threshold_bin_count)+' counts in the '+h.threshold_bin_nm.join('–')+' nm bin. The denominator is '+number(h.histogram_total)+' counts across '+h.bins+' bins. Field area and acquisition time are unavailable.','footnote'),
    links(link('Extended Data Fig. 5b workbook ↓',h.source),link('Imaging methods ↗','https://doi.org/10.1038/s41592-026-03062-5')));
  const live = card('LIVE U2OS · WALTERSPIEL ET AL.', 'JF630b before illumination');
  live.append(node('strong',pct(f.result.fraction),'ev-number'),
    node('p',number(f.counts.qualified)+' of '+number(f.counts.roi_valid_window)+' localizations have positive fitted uncertainties of at most 20 nm in both x and y.'),
    bar(f.result.fraction),
    node('p','The full '+number(f.denominators.frames)+'-frame field covers '+format(f.denominators.area_um2)+' µm² in an estimated '+format(f.denominators.elapsed_s)+' s. Yield was '+format(f.result.per_um2)+' localizations/µm², or '+format(f.result.per_um2_s,5)+'/µm²/s. The pixel calibration is '+number(f.acquisition.pixel_size_nm)+' nm.','footnote'),
    links(link('Open S21 source ↗',evidence.source.url),link('View the recording →','#pilot'),
      saveButton('↓ S21 field JSON','s21-precision-field.json',f)));
  columns.append(fixed, live);
  return [model, columns];
}
function calibration(epfl) {
  const el = card('EPFL · SIMULATED TUBULIN WITH KNOWN COORDINATES', 'Does fitted uncertainty track position error?');
  el.id = 'epfl-calibration';
  el.append(node('p','Each camera image is fitted with Picasso and matched one-to-one to emitters in the same frame within 150 nm. Recall uses all true emitter appearances. Calibration uses only matched fits.'));
  const data = epfl.datasets;
  el.append(table(['Measurement', ...data.map(d=>d.id==='sparse'?'Sparse sequence':'Dense sequence')],[
    ['Sensor frames', ...data.map(d=>number(d.frames))],
    ['True appearances / accepted fits', ...data.map(d=>number(d.truth_count)+' / '+number(d.detection_count))],
    ['Recall / false discovery rate', ...data.map(d=>pct(d.recall)+' / '+pct(d.fdr))],
    ['Median actual radial error', ...data.map(d=>format(d.median_radial_error_nm)+' nm')],
    ['True positions inside nominal 95% ellipse', ...data.map(d=>pct(d.coverage95))],
    ['Actual radial error ≤20 nm among matched fits passing both axes', ...data.map(d=>pct(d.precision20.fraction_actual_radial_le20_given_matched_pass))]
  ],'Calibration against position ground truth'));
  const detail = node('details');
  detail.append(node('summary','Inspect crowding and chance matches'));
  detail.append(table(['Sequence','Recall, neighbor ≤250 nm','Recall, no such neighbor','Wrong-frame matches','Same-frame matches'],data.map(d=>[
    d.id, pct(d.ground_truth_crowding.neighbor_within_250nm.recall),
    pct(d.ground_truth_crowding.no_neighbor_within_250nm.recall),
    number(d.wrong_frame_control.matches),number(d.matches)]),'Sensitivity to crowding and chance matching'));
  el.append(detail,node('p','A nearby emitter can shift the fitted center while the fitter still reports a small uncertainty. Chance matches are frequent in the dense sequence. These simulated coordinates test localization performance under crowding.','footnote'),
    links(link('Source sequence · sparse ↗',data[0].source_url),link('Dense ↗',data[1].source_url),
      link('Calibration JSON ↓','data/epfl-precision.json'),link('Fitting and matching protocol ↗',repo+'analysis/epfl-precision.md')));
  return el;
}
function archiveQC(evidence) {
  const el = card('LIVE-CELL ARCHIVE', 'Recording-level quality checks');
  const entries = [evidence.s21,...evidence.auxiliary], manifest = evidence.manifest;
  const excluded = manifest.filter(m=>m.status==='quarantined');
  const better = entries.filter(f=>f.qc_forecast?.model_improved).length;
  el.append(node('p',entries.length+' of '+manifest.length+' localization tables have unambiguous recording metadata. The trend forecast had lower count error than the first-200-frame mean in '+better+' of '+entries.length+' recordings. Before, during and recovery segments can come from the same experiment. No across-recording independence is assumed.'));
  const filter = node('select'); filter.id='qc-filter';
  const tags = [...new Set(entries.map(f=>f.acquisition.tag))];
  for (const value of ['All tags',...tags]) filter.add(new Option(value,value));
  const label = node('label','Show recording tag'); label.append(filter);
  const controls=node('div',undefined,'ev-controls'); controls.append(label);
  const detail=node('details'); detail.append(node('summary','Inspect individual counts, yields and held-out errors'));
  const results=node('div');
  function paint(){
    const rows=entries.filter(f=>filter.value==='All tags' || f.acquisition.tag===filter.value).map(f=>{
      const a=f.acquisition, q=f.qc_forecast;
      return [f.source.figure, a.field_id, a.tag, a.reference_dye, a.illumination_state,
        number(f.counts.qualified)+' / '+number(f.counts.roi_valid_window),pct(f.result.fraction),
        format(f.result.per_um2_s,5),format(q.model_mae),format(q.baseline_mae)];
    });
    results.replaceChildren(table(['Figure','Recording ID','Tag','Dye','Illumination','Passing / in ROI','Passing','Per µm²/s','Trend MAE','Constant MAE'],rows,'Individual quality checks'));
  }
  paint(); filter.addEventListener('change',paint);
  detail.append(controls,results,node('p','Forecast errors are absolute count errors averaged over the next eight 100-frame blocks. A missing time calibration leaves the rate unavailable.','footnote'));
  for (const m of excluded) detail.append(node('p','Excluded from the field summaries: Fig. '+m.figure+' · '+m.recording+'. '+(m.reason || 'The archive and recording labels disagree on the tag.'),'footnote'));
  el.append(detail,links(link('Archive and original recordings ↗',evidence.source.url)));
  return el;
}

function displayField(f) {
  const el = node('div');
  el.append(node('h3',(f.acquisition.reference_dye || 'Unknown dye')+' · '+(f.acquisition.field_id || 'Field summary')),
    node('p',f.definition.name+' · '+f.definition.threshold_nm+' nm cutoff · '+f.result.status),
    table(['In ROI','Valid precision','Missing','Invalid','Qualifying','Fraction','Per µm²','Per µm²/s'],[[
      number(f.counts.roi_valid_window),number(f.counts.valid_precision),number(f.counts.missing_precision),
      number(f.counts.invalid_precision),number(f.counts.qualified),pct(f.result.fraction),
      format(f.result.per_um2,4),format(f.result.per_um2_s,5)]],'Imported field precision'),
    node('p','Area '+format(f.denominators.area_um2)+' µm² · '+number(f.denominators.frames)+' frames · '+format(f.denominators.elapsed_s,4)+' s.','footnote'));
  return el;
}
function displayComparison(c) {
  const el = node('div');
  if(c.status!=='matched') {
    const reasons=node('ul'); for(const reason of c.reasons) reasons.append(node('li',reason));
    el.append(node('h3','Matched comparison unavailable'),reasons); return el;
  }
  const [a,b] = c.dyes;
  el.append(node('h3',a+' minus '+b),
    node('p','Mean paired difference '+format(c.mean_paired_difference_per_um2_s,5)+' localizations/µm²/s across '+number(c.preparations.length)+' independent preparation ID(s).'),
    table(['Preparation',a+' fields',b+' fields',a+' per µm²/s',b+' per µm²/s','Difference per µm²/s','Ratio'],c.preparations.map(p=>[
      p.prep_id,number(p.numerator.n_fields),number(p.denominator.n_fields),format(p.numerator.per_um2_s,5),
      format(p.denominator.per_um2_s,5),format(p.difference_per_um2_s,5),format(p.ratio_per_um2_s,3)]),'Observed differences within preparation'));
  el.append(node('p',c.model_domain==='reference_acquisition' ? 'Acquisition settings match the H2B reference protocol.' : 'Acquisition settings fall outside the modeled H2B reference protocol.','footnote'));
  return el;
}
function importer() {
  const el = card('YOUR LOCAL ANALYSIS', 'Compare matched fields');
  const file=node('input'); file.type='file'; file.multiple=true; file.accept='.json,application/json'; file.hidden=true; file.id='precision-file';
  const open=node('button','↑ Open field or comparison JSON','button secondary'); open.type='button'; open.onclick=()=>file.click();
  const clear=node('button','Clear','text-button'); clear.type='button'; clear.hidden=true;
  const status=node('p',undefined,'ev-status-message'); status.setAttribute('role','status');
  const result=node('div',undefined,'ev-report'); result.id='precision-import-result';
  const options=node('div',undefined,'ev-controls'), actions=node('div',undefined,'ev-controls');
  actions.append(open,file,clear);
  const reset=()=>{result.replaceChildren();options.replaceChildren();status.textContent='';clear.hidden=true;};
  clear.onclick=reset;
  let selection=0;
  file.onchange=async()=>{
    const chosen=Array.from(file.files), run=++selection; file.value='';
    if(!chosen.length) return;
    reset();
    try {
      const objects=await Promise.all(chosen.map(async f=>JSON.parse(await f.text())));
      if(run!==selection) return;
      if(objects.length===1 && objects[0].schema==='smlm-precision-comparison/1'){
        validateComparison(objects[0]); result.append(displayComparison(objects[0]));
      } else {
        objects.forEach(validateField);
        if(objects.length===1) result.append(displayField(objects[0]));
        else {
          const dyes=[...new Set(objects.map(f=>f.acquisition.reference_dye).filter(Boolean))];
          const selectDye=(label, preference)=>{
            const select=node('select');
            for(const dye of dyes) select.add(new Option(dye,dye));
            if(dyes.includes(preference)) select.value=preference;
            const wrap=node('label',label); wrap.append(select); options.append(wrap); return select;
          };
          const a=selectDye('First dye',dyes.includes('JF639b')?'JF639b':dyes[0]);
          const b=selectDye('Second dye',dyes.includes('JF630b') && a.value!=='JF630b'?'JF630b':dyes.find(x=>x!==a.value));
          const paint=()=>{
            if(!a.value || !b.value || a.value===b.value){result.replaceChildren(node('p','Select two distinct dyes from the loaded fields.'));return;}
            const fields=objects.filter(f=>[a.value,b.value].includes(f.acquisition.reference_dye));
            const comparison=compareFields(fields,[a.value,b.value]);
            result.replaceChildren(displayComparison(comparison),
              saveButton('↓ Export comparison JSON','precision-comparison.json',comparison));
          };
          a.onchange=paint; b.onchange=paint; paint();
        }
      }
      clear.hidden=false; status.textContent=number(objects.length)+' local file(s) opened.';
    } catch(error) {reset(); status.textContent='Import failed. '+error.message;}
  };
  el.append(node('p','Generate summaries locally from a Picasso table with metadata, or a calibrated CSV. Open a single field to inspect its denominator or select multiple files for a within-preparation comparison. Files stay in this browser.'),
    actions,options,status,result,links(link('CLI and required metadata ↗',repo+'analysis/precision-input.md')));
  return el;
}
export async function renderPrecision(container) {
  const read=async name=>{
    const r=await fetch('data/'+name+'.json');
    if(!r.ok) throw Error(name+' returned HTTP '+r.status);
    return r.json();
  };
  const [evidence,epfl]=await Promise.all([read('precision'),read('epfl-precision')]);
  container.replaceChildren(...publicResults(evidence),importer(),calibration(epfl),archiveQC(evidence));
}
