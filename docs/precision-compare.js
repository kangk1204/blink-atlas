// Portable summaries from analysis.precision. All counts remain browser-local.
const FIELD = "smlm-precision-field/1";
const COMPARISON = "smlm-precision-comparison/1";
const AXES = "both lateral 1-sigma uncertainties";
const STR = ["sample_state", "cell_line", "target", "tag", "illumination_state",
  "density_calibration_id", "acquisition_id", "analysis_id", "roi_protocol_id", "density_origin", "density_unit"];
const NUM = ["density_labels_per_um2", "excitation_nm", "exposure_ms", "frame_interval_ms", "start_frame", "window_frames"];
const number = n => typeof n === "number" && Number.isFinite(n);
const count = n => Number.isSafeInteger(n) && n >= 0;
const positive = n => number(n) && n > 0;
const object = o => o !== null && typeof o === "object" && !Array.isArray(o);
const label = s => typeof s === "string" && s.trim().length > 0 && s.trim().length <= 200;
const norm = s => typeof s === "string" ? s.trim().toLowerCase() : "";
const near = (a,b) => b === null ? a === null : number(a) && number(b) &&
  Math.abs(a-b) <= Math.max(1e-10,1e-8*Math.max(Math.abs(a),Math.abs(b)));
const sameSetting = (a,b) => number(a) && number(b) &&
  Math.abs(a-b) <= Math.max(1e-9,1e-9*Math.max(Math.abs(a),Math.abs(b)));
function require(ok, message) { if (!ok) throw new Error(message); }
function checkDefinition(d) {
  require(object(d) && d.threshold_nm === 20 && d.unit === "nm" &&
    ((d.kind === "axes" && d.name === AXES) || (d.kind === "scalar" && label(d.name))),
  "Unsupported precision definition.");
}
function sameDefinition(a,b) { return ["threshold_nm","unit","name","kind"].every(k=>a[k]===b[k]); }
function duration(a) {
  if (a.window_frames === 1) return a.exposure_ms/1000;
  return positive(a.frame_interval_ms) ? ((a.window_frames-1)*a.frame_interval_ms+a.exposure_ms)/1000 : null;
}
function inAssay(a) {
  return norm(a.sample_state)==="fixed" && norm(a.cell_line).replaceAll("-","")==="cos7" &&
    norm(a.target)==="h2b" && norm(a.tag)==="halotag";
}
function modelDomain(a) {
  return Math.abs(a.excitation_nm-647)<=1 && Math.abs(a.exposure_ms-25)<=.2 && Math.abs(a.frame_interval_ms-27)<=.2
    ? "reference_acquisition" : "outside_reference_acquisition";
}

export function validateField(f) {
  require(object(f) && f.schema===FIELD, "Expected a smlm-precision-field/1 summary.");
  const {acquisition:a,definition:d,counts:c,denominators:v,result:r,source:src} = f;
  require([a,d,c,v,r,src].every(object), "Incomplete field summary.");
  require(["csv","picasso"].includes(src.format) && typeof src.table_sha256==="string" &&
    /^[0-9a-f]{64}$/.test(src.table_sha256) && typeof src.precision_columns_present==="boolean", "Invalid localization source.");
  checkDefinition(d);
  for (const k of ["width_nm","height_nm","camera_width_nm","camera_height_nm","exposure_ms"])
    require(positive(a[k]), "Invalid acquisition "+k+".");
  for (const k of ["frame_interval_ms","excitation_nm","pixel_size_nm","density_labels_per_um2"])
    require(a[k]===null || positive(a[k]), "Invalid acquisition "+k+".");
  for (const k of ["sample_state","cell_line","target","tag","reference_dye","prep_id","field_id","illumination_state",
    "density_calibration_id","density_origin","acquisition_id","analysis_id","roi_protocol_id","timing_method"])
    require(a[k]===null || label(a[k]), "Invalid acquisition label "+k+".");
  require(a.density_unit==="labels/um2", "Invalid density unit.");
  if (src.format==="picasso") require(positive(a.pixel_size_nm), "Missing Picasso pixel calibration.");
  require([a.total_frames,a.start_frame,a.window_frames].every(count) && a.window_frames>0 &&
    a.start_frame+a.window_frames<=a.total_frames,"Invalid acquisition frame window.");
  require(Array.isArray(a.roi_nm) && a.roi_nm.length===4 && a.roi_nm.every(number),"Invalid physical ROI.");
  const [left,top,right,bottom]=a.roi_nm;
  require(left>=0 && right>left && right<=a.camera_width_nm && top>=0 && bottom>top && bottom<=a.camera_height_nm &&
    near(a.width_nm,right-left) && near(a.height_nm,bottom-top),"Inconsistent physical ROI.");
  require(["source","outside_window","outside_roi","invalid_geometry","roi_valid_window","missing_precision","invalid_precision",
    "valid_precision"].every(k=>count(c[k])), "Invalid localization counts.");
  require(c.source===c.outside_window+c.outside_roi+c.invalid_geometry+c.roi_valid_window &&
    c.roi_valid_window===c.missing_precision+c.invalid_precision+c.valid_precision,"Inconsistent localization totals.");
  const n=a.window_frames, area=a.width_nm*a.height_nm/1e6, elapsed=duration(a);
  require(v.frames===n && count(v.empty_frames) && n-v.empty_frames>=0 && n-v.empty_frames<=Math.min(n,c.roi_valid_window) &&
    near(v.area_um2,area) && near(v.elapsed_s,elapsed), "Invalid area/time or empty-frame denominators.");
  let status="measured";
  if (!src.precision_columns_present || (c.roi_valid_window && !c.valid_precision)) status="unavailable";
  else if (c.invalid_geometry || c.missing_precision || c.invalid_precision) status="partial";
  const q=c.qualified;
  require(r.status===status && (status==="unavailable" ? q===null : count(q) && q<=c.valid_precision),
    "Invalid precision coverage or qualified counts.");
  require(near(r.fraction, q!==null && c.valid_precision ? q/c.valid_precision : null) &&
    near(r.per_um2,q!==null ? q/area : null) && near(r.per_um2_s,q!==null && elapsed ? q/(area*elapsed) : null),
    "Invalid precision rates.");
  return f;
}

function groupFields(fields) {
  let area=0, areaTime=0, qualified=0;
  for (const f of fields) {
    area+=f.denominators.area_um2;
    areaTime+=f.denominators.area_um2*f.denominators.elapsed_s;
    qualified+=f.counts.qualified;
  }
  return {n_fields:fields.length,qualified,area_um2:area,area_time_um2_s:areaTime,per_um2:qualified/area,per_um2_s:qualified/areaTime};
}

export function compareFields(fields, dyes) {
  require(Array.isArray(dyes) && dyes.length===2 && dyes.every(label) && norm(dyes[0])!==norm(dyes[1]), "Choose two distinct dyes (numerator first).");
  require(Array.isArray(fields), "Supply an array of field summaries.");
  const result={schema:COMPARISON,dyes:[...dyes],status:"unavailable",reasons:[],preparations:[],n_fields:fields.length,definition:null};
  if (!fields.length) { result.reasons=["No field summaries supplied."]; return result; }
  fields.forEach(validateField);
  const reference=fields[0], reasons=new Set(), seen=new Set(), groups=new Map();
  for (const f of fields) {
    const a=f.acquisition, prep=a.prep_id;
    const dye=dyes.find(x=>norm(x)===norm(a.reference_dye));
    if (dye===undefined) reasons.add("Unexpected dye in the supplied fields.");
    if (f.result.status!=="measured") reasons.add("Each field needs complete precision and geometry (a verified zero is allowed).");
    if (!prep || !a.field_id) reasons.add("Both prep_id and field_id are required for independent pairing.");
    if (!inAssay(a)) reasons.add("The paired endpoint requires fixed COS-7 H2B–HaloTag fields.");
    const key=JSON.stringify([prep,dye??null,a.field_id]);
    if (seen.has(key)) reasons.add("Duplicate prep/dye/field identity.");
    seen.add(key);
    if (prep && dye!==undefined) {
      if (!groups.has(prep)) groups.set(prep,new Map());
      const g=groups.get(prep);
      if (!g.has(dye)) g.set(dye,[]);
      g.get(dye).push(f);
    }
    for (const k of STR) {
      if (!label(a[k]) || ["unknown","unavailable"].includes(norm(a[k]))) reasons.add("Missing matched metadata: "+k+".");
      else if (norm(reference.acquisition[k])!==norm(a[k])) reasons.add("Incompatible "+k+".");
    }
    if (norm(a.density_origin)!=="measured" || a.density_unit!=="labels/um2")
      reasons.add("Actual independently calibrated labeling density in labels/um2 is required.");
    for (const k of NUM) {
      if (!number(a[k]) || (k==="start_frame" ? a[k]<0 : a[k]<=0)) reasons.add("Missing matched metadata: "+k+".");
      else if (!sameSetting(a[k],reference.acquisition[k])) reasons.add("Incompatible "+k+".");
    }
    if (!sameDefinition(f.definition,reference.definition)) reasons.add("Precision definition, units or threshold differ.");
    if (!positive(f.denominators.area_um2) || !positive(f.denominators.elapsed_s)) reasons.add("A calibrated ROI area and acquisition duration are required.");
  }
  for (const [prep,g] of groups) if (g.size!==2) reasons.add("Preparation "+prep+" lacks one of the two dyes.");
  if (reasons.size) { result.reasons=[...reasons].sort(); return result; }
  result.status="matched";
  result.definition=reference.definition;
  result.matched_context=Object.fromEntries([...STR,...NUM].map(k=>[k,reference.acquisition[k]]));
  result.model_domain=modelDomain(reference.acquisition);
  for (const prep of [...groups.keys()].sort()) {
    const g=groups.get(prep);
    const [n,d]=dyes.map(dye=>groupFields(g.get(dye)));
    result.preparations.push({prep_id:prep,numerator:n,denominator:d,difference_per_um2:n.per_um2-d.per_um2,
      difference_per_um2_s:n.per_um2_s-d.per_um2_s,ratio_per_um2_s:d.per_um2_s>0 ? n.per_um2_s/d.per_um2_s : null});
  }
  result.mean_paired_difference_per_um2_s=result.preparations.reduce((s,p)=>s+p.difference_per_um2_s,0)/result.preparations.length;
  return result;
}

export function validateComparison(c) {
  require(object(c) && c.schema===COMPARISON && Array.isArray(c.dyes) && c.dyes.length===2 && c.dyes.every(label) &&
    norm(c.dyes[0])!==norm(c.dyes[1]) && count(c.n_fields) && Array.isArray(c.preparations) && Array.isArray(c.reasons),
    "Invalid comparison summary.");
  if (c.status==="unavailable") {
    require(c.preparations.length===0 && c.reasons.length>0 && c.reasons.every(s=>typeof s==="string" && s.length>0) &&
      c.definition===null && c.mean_paired_difference_per_um2_s==null,"Unavailable comparisons cannot contain an effect estimate.");
    return c;
  }
  require(c.status==="matched" && c.reasons.length===0 && c.preparations.length>0,"Invalid comparison state.");
  checkDefinition(c.definition);
  const m=c.matched_context;
  require(object(m) && STR.every(k=>label(m[k]) && !["unknown","unavailable"].includes(norm(m[k]))) && inAssay(m) &&
    norm(m.density_origin)==="measured" && m.density_unit==="labels/um2" && NUM.every(k=>k==="start_frame" ? count(m[k]) : positive(m[k])) &&
    count(m.window_frames) && c.model_domain===modelDomain(m),"Invalid matched acquisition context.");
  const seen=new Set(); let nFields=0,sum=0;
  for (const p of c.preparations) {
    require(object(p) && label(p.prep_id) && !seen.has(p.prep_id),"Invalid or duplicated preparation.");
    seen.add(p.prep_id);
    for (const g of [p.numerator,p.denominator]) {
      require(object(g) && count(g.n_fields) && g.n_fields>0 && count(g.qualified) && positive(g.area_um2) && positive(g.area_time_um2_s) &&
        near(g.area_time_um2_s,g.area_um2*duration(m)) && near(g.per_um2,g.qualified/g.area_um2) && near(g.per_um2_s,g.qualified/g.area_time_um2_s),
        "Invalid preparation-level counts or rates.");
      nFields+=g.n_fields;
    }
    const n=p.numerator,d=p.denominator;
    require(near(p.difference_per_um2,n.per_um2-d.per_um2) && near(p.difference_per_um2_s,n.per_um2_s-d.per_um2_s) &&
      near(p.ratio_per_um2_s,d.per_um2_s>0 ? n.per_um2_s/d.per_um2_s : null),"Invalid preparation-level contrast.");
    sum+=p.difference_per_um2_s;
  }
  require(nFields===c.n_fields && near(c.mean_paired_difference_per_um2_s,sum/c.preparations.length),
    "Invalid cross-preparation mean or field count.");
  return c;
}
