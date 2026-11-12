import {lockTrial, methods, predict, restoreTrial, scoreTrial} from './probe-trial.js';

const id = name => document.getElementById('probe-' + name);
const format = value => value.toExponential(3);

function download(trial, includeMeasurements) {
  const contents = {plan: trial.plan, measurement: includeMeasurements ? trial.measurement : null};
  const blob = new Blob([JSON.stringify(contents, null, 2) + '\n'], {type:'application/json'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'blink-atlas-' + (includeMeasurements ? 'scored' : 'plan') + '-' +
    trial.plan.createdUTC.slice(0, 19).replace(/:/g, '-') + '.json';
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function mountTrial(root, data) {
  root.innerHTML = `
    <div class="eyebrow">PROSPECTIVE EXPERIMENT</div><h2 id="probe-trial-heading" tabindex="-1">Test a new analogue</h2>
    <p class="probe-sub">Give the compound a blinded ID and enter the apparent pKa measured on its free parent in buffer. Save the predictions before imaging the corresponding HaloTag ligand.</p>
    <div class="probe-trial-grid">
      <div>
        <form id="probe-plan-form" novalidate>
          <fieldset id="probe-plan-inputs"><legend class="sr">New compound and assay</legend>
            <div class="probe-two-fields"><label>Blinded compound ID<input id="probe-compound" type="text" maxlength="80" autocomplete="off" placeholder="New analogue ID"></label>
              <label>Measured parent pKa<input id="probe-pka" type="number" step="any" inputmode="decimal" placeholder="e.g. 3.20"></label></div>
            <label>Chemistry and planned assay<select id="probe-scope"><option value="">Select the experimental conditions</option></select></label>
          </fieldset>
          <div id="probe-preview" role="status"></div>
          <div class="probe-actions"><button class="button primary" type="submit" id="probe-lock">Lock &amp; download predictions</button>
            <button class="button secondary" type="button" id="probe-save-plan" hidden>Download locked plan</button>
            <button class="text-button" type="button" id="probe-clear" hidden>Start a new card</button></div>
        </form>
        <details class="probe-import" id="probe-import-details"><summary>Open a saved trial card</summary>
          <label>Local prediction or scored JSON<input id="probe-import" type="file" accept=".json,application/json"></label>
        </details>
        <p id="probe-trial-status" role="status"></p>
      </div>
      <div class="probe-readout">
        <h3>Matched H2B measurement</h3>
        <p>Use fixed COS-7 H2B–HaloTag under the published subsaturating labeling protocol. Acquire ON times at 2 ms exposure and OFF times at 25 ms exposure. Include the reference ligands printed on the locked card.</p>
        <p>Enter one mean ON:OFF time ratio per independent preparation. The card scores the three frozen predictions against the arithmetic mean of the three ratios. OFF times in this assay are upper-bound estimates.</p>
        <form id="probe-score-form" hidden novalidate>
          <fieldset id="probe-observations"><legend>Independent preparation ratios</legend><div class="probe-measurements">
            <label>Prep 1<input id="probe-rep-1" type="number" step="any" inputmode="decimal"></label>
            <label>Prep 2<input id="probe-rep-2" type="number" step="any" inputmode="decimal"></label>
            <label>Prep 3<input id="probe-rep-3" type="number" step="any" inputmode="decimal"></label>
          </div></fieldset>
          <button id="probe-score" class="button secondary" type="submit">Score measurements</button>
          <p id="probe-score-result" role="status"></p>
          <button id="probe-save-score" class="button primary" type="button" hidden>Download scored card</button>
        </form>
      </div>
    </div>`;
  const c = data.calibration;
  id('scope').append(new Option('Monotonic hydroxymethyl SiR · fixed COS-7 H2B–HaloTag', c.scope),
    new Option('Nonmonotonic or unmeasured titration', 'complex'), new Option('Different scaffold, tag, target or protocol', 'outside'));
  let trial = null;
  function preview() {
    const target = id('preview');
    target.replaceChildren();
    if (!trial && !id('pka').value && !id('scope').value) return;
    let prediction;
    try {
      prediction = trial?.plan || predict(c, id('pka').value, id('scope').value);
    } catch (err) {
      target.textContent = err.message;
      return;
    }
    const caption = document.createElement('strong');
    caption.textContent = trial ? 'Locked ON:OFF predictions' : 'Predicted mean ON:OFF time ratios';
    const list = document.createElement('dl');
    for (const [key, label] of Object.entries(methods)) {
      const term = document.createElement('dt'), value = document.createElement('dd');
      term.textContent = label + (key === 'nearest' ? ' · ' + prediction.nearestDye : '');
      value.textContent = format(prediction.ratios[key]);
      list.append(term, value);
    }
    const refs = document.createElement('p');
    refs.textContent = 'The reference ligands at the bounding pKa values are ' + prediction.references.join(' and ') + '.';
    target.append(caption, list, refs);
  }
  const inputs = [id('rep-1'), id('rep-2'), id('rep-3')];
  function renderState() {
    const locked = Boolean(trial), measured = Boolean(trial?.measurement);
    id('plan-inputs').disabled = locked;
    id('lock').hidden = locked;
    id('save-plan').hidden = !locked;
    id('clear').hidden = !locked;
    id('score-form').hidden = !locked;
    id('score').hidden = measured;
    id('observations').disabled = measured;
    id('save-score').hidden = !measured;
    if (locked) {
      id('compound').value = trial.plan.compound;
      id('pka').value = trial.plan.pKa;
      id('scope').value = trial.plan.scope;
    }
    inputs.forEach((input, i) => {input.value = measured ? trial.measurement.ratios[i] : '';});
    preview();
    const text = id('score-result');
    text.textContent = '';
    if (measured) {
      const m = trial.measurement;
      const errors = Object.keys(methods).map(k => methods[k] + ' ' + m.error[k].toFixed(3));
      text.textContent = 'Observed mean ' + format(m.observed) + '. Absolute log₁₀ errors — ' + errors.join(' · ') +
        '. Closest for this compound — ' + m.closest.map(k => methods[k]).join(' and ') + '.';
    }
  }
  id('pka').oninput = id('scope').onchange = () => {id('trial-status').textContent = ''; preview();};
  id('compound').oninput = () => {id('trial-status').textContent = '';};
  id('plan-form').onsubmit = event => {
    event.preventDefault();
    if (trial) return;
    try {
      trial = lockTrial(data, id('compound').value, id('pka').value, id('scope').value);
      renderState();
      download(trial, false);
      id('trial-status').textContent = 'Predictions locked for ' + trial.plan.compound + '. Keep the downloaded card before collecting measurements.';
    } catch (err) {id('trial-status').textContent = err.message;}
  };
  id('score-form').onsubmit = event => {
    event.preventDefault();
    if (!trial || trial.measurement) return;
    try {
      trial.measurement = scoreTrial(trial.plan, inputs.map(input => input.value));
      renderState();
      id('trial-status').textContent = '';
    } catch (err) {id('score-result').textContent = err.message;}
  };
  id('save-plan').onclick = () => download(trial, false);
  id('save-score').onclick = () => download(trial, true);
  id('clear').onclick = () => {
    trial = null;
    id('plan-form').reset();
    id('trial-status').textContent = '';
    renderState();
    id('compound').focus();
  };
  id('import').onchange = async event => {
    const picker = event.target, file = picker.files[0];
    if (!file) return;
    try {
      const restored = restoreTrial(data, JSON.parse(await file.text()));
      trial = restored;
      renderState();
      id('trial-status').textContent = 'Restored ' + trial.plan.compound + ' and verified the stored predictions against this release.';
    } catch (err) {
      id('trial-status').textContent = 'Could not open card. ' + err.message;
    } finally {picker.value = '';}
  };
}
