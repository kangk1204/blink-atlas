const models = {local:'Local interpolation', nearest:'Nearest pKa', regression:'Original regression'};
const science = value => value.toExponential(3);
const range = (values, format) => values.map(format).join(' – ');
const f4 = value => value.toFixed(4);
const percent = n => (n >= 0 ? '+' : '') + n.toFixed(3) + '%';

function tableRows(body, rows) {
  body.replaceChildren(...rows.map(entries => {
    const tr = document.createElement('tr');
    for (const value of entries) {
      const td = document.createElement('td');
      td.textContent = value;
      tr.append(td);
    }
    return tr;
  }));
}

export function renderRatioAudit(root, data, evidence) {
  root.innerHTML = '<div class="eyebrow">HOLLAND · EXTENDED DATA FIGURE 4</div>' +
    '<h2>Does the way we average ON and OFF times change the result?</h2>' +
    '<p id="ratio-intro"></p><div class="ratio-controls"><label>Time-ratio calculation<select id="ratio-rule"></select></label>' +
    '<label>Source dye<select id="ratio-dye"></select></label></div>' +
    '<div aria-live="polite"><div class="table-wrap"><table><thead><tr><th>Model</th><th>Held-out mean absolute log₁₀ error</th><th>Cal664sb forecast</th></tr></thead>' +
    '<tbody id="ratio-models"></tbody></table></div><p id="ratio-boundaries"></p><p id="ratio-rank"></p></div>' +
    '<details><summary>Inspect the source dye and its six possible pairings</summary>' +
    '<p id="ratio-cells"></p><p id="ratio-reported"></p><div class="table-wrap"><table><thead><tr><th>OFF order for ON columns B, C, D</th>' +
    '<th>Mean of three ON/OFF ratios</th><th>Difference from reported</th></tr></thead><tbody id="ratio-pairs"></tbody></table></div></details>' +
    '<details><summary>How were the uncertainty bars drawn?</summary><p id="ratio-sd"></p>' +
    '<a id="ratio-figure" class="inline-link">Inspect the original figure ↗</a></details>' +
    '<div class="unseen-links"><a href="data/ratio-audit.csv" download>9-dye source audit (CSV)</a>' +
    '<a href="data/ratio-sensitivity.json" download>Model results (JSON)</a>' +
    '<a href="data/ratio-scenarios.csv.gz" download>All pairing scenarios (CSV.GZ)</a>' +
    '<a href="https://github.com/kangk1204/blink-atlas/blob/main/analysis/ratio-audit.md">Methods &amp; reproduction ↗</a></div>';
  const get = name => root.querySelector('#ratio-' + name);
  const scan = data.pairings, n = scan.scenarioCount.toLocaleString('en-US');
  const matched = data.audit.filter(r => Math.abs(r.ratio_of_means_diff_pct) <= .1).length;
  const exception = data.audit.find(r => r.dye === 'JF626b');
  get('intro').textContent = 'For ' + matched + ' of ' + data.audit.length + ' dyes, mean ON divided by mean OFF reproduces the published ratio within 0.1%. ' +
    exception.dye + ' differs by ' + exception.ratio_of_means_diff_pct.toFixed(2) + '%. The ON and OFF recordings used different labeling concentrations. The workbook gives no preparation IDs linking their columns.';
  const rule = get('rule'), dye = get('dye');
  for (const [name, label] of Object.entries({
    exhaustive:'All possible pairings · sensitivity range', j_reported:'Published ratios · registration inputs',
    ratio_of_means:'Mean ON divided by mean OFF', paired_BCD:'Mean of three ratios · same-column pairing'
  })) rule.add(new Option(label, name));
  for (const r of data.audit) dye.add(new Option(r.dye, r.dye));
  function showModel() {
    const all = rule.value === 'exhaustive', set = data.estimands[rule.value];
    tableRows(get('models'), Object.entries(models).map(([key, label]) => [label,
      all ? range(scan.overall.errors[key], f4) : f4(set.metrics.mae[key]),
      all ? range(scan.calForecastRanges[key], science) : science(set.cal[key])]));
    const bounds = all ? scan.calOrder.flatMap(o => o.adjacentBoundaries) : set.calBoundaries;
    get('boundaries').textContent = 'Cal664sb equal-error thresholds at pKa ' + data.registered.pKa + ' · ' +
      bounds.map(v => v.models.map(m => models[m]).join(' / ') + ' ' +
        (all ? range(v.range, science) : science(v.ratio))).join(' · ') +
      '. Above the upper boundary local interpolation is closest. The registered thresholds and forecasts are available above.';
    const wins = scan.folds.filter(f => f.comparisons?.local_vs_nearest.firstWins === scan.scenarioCount);
    get('rank').textContent = all ?
      'All ' + n + ' scenarios rank local interpolation first by average error across ' + scan.scoredDyes +
      ' in-range holdouts. Only ' + wins.map(f => f.dye).join(', ') + ' favors it over nearest pKa on its individual fold. ' +
      'The two end-point dyes fall outside training support when held out. The ranges describe sensitivity to pairing.' :
      set.metrics.scored + ' of ' + set.metrics.total + ' development dyes are in range after holding each dye out. Local interpolation beats nearest pKa on ' +
      set.metrics.wins + ' individual fold and loses on ' + set.metrics.losses + '.';
  }
  function showDye() {
    const r = data.audit.find(d => d.dye === dye.value);
    get('cells').textContent = 'ON (ms) ' + ['B','C','D'].map(c => r['on_'+c+'_ms']).join(', ') + ' · ' + r.on_cells_ms +
      '. OFF (ms) ' + ['B','C','D'].map(c => r['off_'+c+'_ms']).join(', ') + ' · ' + r.off_cells_ms + '.';
    get('reported').textContent = 'Published ' + science(r.reported) + ' · ' + r.reported_cells_ratio_sd_n +
      '. Ratio of means ' + science(r.ratio_of_means) + ' (' + percent(r.ratio_of_means_diff_pct) + ').';
    tableRows(get('pairs'), scan.pairingLabels.map(p => [p, science(r['paired_'+p]), percent(r['paired_'+p+'_diff_pct'])]));
  }
  const u = evidence.uncertaintyAudit;
  get('sd').textContent = 'Panels 4h–i label their bars SD. The same ' + u.sharedSdCells + ' SD values appear in worksheet 4j, whose figure caption specifies SEM. ' +
    'For JFX650b and JF639b the visible half-bars agree with SD divided by √3 within the figure’s pixel resolution. The model uses the means.';
  get('figure').href = data.input.figure;
  rule.onchange = showModel;
  dye.onchange = () => {
    showDye();
    get('cells').closest('details').open = true;
  };
  dye.value = exception.dye;
  showModel(); showDye();
}

export function renderRatioEvidence(root, evidence) {
  root.innerHTML = '<h2>Which additional published dyes could be scored?</h2><p id="ratio-coverage"></p>' +
    '<p id="ratio-examples"></p><details><summary>Inspect each parent dye and its eligibility</summary>' +
    '<div class="table-wrap"><table><thead><tr><th>Study</th><th>Parent</th><th>Free-dye input</th><th>Measurement and score status</th>' +
    '</tr></thead><tbody id="ratio-screen-rows"></tbody></table></div></details>' +
    '<div class="unseen-links"><a href="data/ratio-evidence.json" download>Per-dye evidence and scores (JSON)</a>' +
    '<a id="ratio-uno">Uno et al. · Table 1 ↗</a><a id="ratio-chen">Chen et al. · supporting figures ↗</a></div>';
  root.querySelector('#ratio-coverage').textContent = evidence.countsByStudy['Uno 2014'] + ' parent dyes from Uno et al. and ' +
    evidence.countsByStudy['Chen 2024'] + ' from Chen et al. were screened. ' + evidence.matchedScored +
    ' had an independent matched H2B measurement. ' + evidence.frozenTransferScored + ' supplied a scorable transfer pair for the frozen models.';
  const ama = evidence.records.find(r => r.dye === 'AMAcRG');
  root.querySelector('#ratio-examples').textContent = 'AMAcRG is in range at pKcycl ' + ama.pKa +
    '. The flash experiment reports an open-state relaxation time of ' + (ama.flashOpenLifetimeSeconds*1e6).toFixed(2) +
    ' µs. An independent mean OFF dwell is missing. Chen’s SMI-NIR1–Biotin panel contributes ' +
    evidence.chenByCondition.length + ' duty measurements on one parent. Its free-dye optical midpoint is above the calibration range.';
  tableRows(root.querySelector('#ratio-screen-rows'), evidence.records.map(r => [r.study, r.dye,
    r.pKa ?? r.pKaBound ?? ('~' + r.opticalMidpointApproxPH.toFixed(1) + ' optical midpoint'), r.reason]));
  root.querySelector('#ratio-uno').href = evidence.articles['Uno 2014'].table;
  root.querySelector('#ratio-chen').href = evidence.articles['Chen 2024'].deposit;
}
