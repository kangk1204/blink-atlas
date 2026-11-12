const MONO = '"IBM Plex Mono", monospace';
export const dyeColor = {
  JF630b: '#c4ef87', JF635b: '#73c6dc', 'HM-SiR': '#d5a88f',
  JFX650b: '#afa1e5', JF646b: '#e5b86e', JF639b: '#f1a9b6',
  JF629b: '#75c5ad', JF626b: '#a9ba81', JF614b: '#a9b6ce'
};
const format = new Intl.NumberFormat('en-US', {maximumFractionDigits: 1});
function compact(n) { return n >= 1000 ? format.format(n / 1000) + 'k' : format.format(n); }

function surface(canvas) {
  const {width, height} = canvas.getBoundingClientRect();
  if (!width || !height) return null;
  const ratio = Math.min(window.devicePixelRatio || 1, 3);
  canvas.width = Math.round(width * ratio);
  canvas.height = Math.round(height * ratio);
  const ctx = canvas.getContext('2d');
  ctx.scale(ratio, ratio);
  ctx.font = '10px ' + MONO;
  ctx.lineJoin = 'round';
  return {ctx, width, height};
}
function line(ctx, points, color, width = 1, dash = []) {
  ctx.beginPath();
  for (let i = 0; i < points.length; i++) {
    if (i === 0) ctx.moveTo(...points[i]); else ctx.lineTo(...points[i]);
  }
  ctx.strokeStyle = color;
  ctx.lineWidth = width;
  ctx.setLineDash(dash);
  ctx.stroke();
  ctx.setLineDash([]);
}
function dot(ctx, x, y, color, radius) {
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.arc(x, y, radius, 0, Math.PI * 2);
  ctx.fill();
}
function niceStep(rough) {
  if (!(rough > 0)) return 1;
  const scale = 10 ** Math.floor(Math.log10(rough));
  return ([1, 2, 2.5, 5, 10].find(n => n >= rough / scale) || 10) * scale;
}

export function drawYield(canvas, points, result, showAll, hoverIndex) {
  const setup = surface(canvas);
  if (!setup || !points.length) return null;
  const {ctx, width, height} = setup;
  const left = 34, right = width - 16, top = 19, bottom = height - 23;
  const lo = Math.log10(points[0].density), hi = Math.log10(points.at(-1).density);
  const x = density => left + (Math.log10(density) - lo) / (hi - lo) * (right - left);
  const reference = ['JF635b', 'JF630b'];
  const ids = showAll ? Object.keys(points[0].rates).filter(k => !reference.includes(k)).concat(reference) : reference;
  const maxRate = Math.max(...points.flatMap(p => ids.map(id => p.rates[id])));
  const step = niceStep(maxRate / 4);
  const maxY = Math.max(step * Math.ceil(maxRate / step), step);
  const y = n => bottom - (n / maxY) * (bottom - top);
  const bx = result.density.map(x);
  ctx.fillStyle = '#b6d89518';
  ctx.fillRect(bx[0], top, Math.max(1.5, bx[1] - bx[0]), bottom - top);
  for (let i = 0; i <= Math.round(maxY / step); i++) {
    const value = i * step, yy = y(value);
    line(ctx, [[left, yy], [right, yy]], '#8db49a37', 1, [2, 4]);
    ctx.fillStyle = '#bdd1c0'; ctx.textAlign = 'right'; ctx.fillText(compact(value), left - 7, yy + 3);
  }
  for (let p = Math.ceil(lo - 1e-8); p <= Math.floor(hi + 1e-8); p++) {
    const xx = x(10 ** p);
    line(ctx, [[xx, top], [xx, bottom]], '#799e8a36', 1, [2, 4]);
    ctx.fillStyle = '#bdd1c0'; ctx.textAlign = 'center'; ctx.fillText(compact(10 ** p), xx, bottom + 15);
  }
  for (const xx of bx) line(ctx, [[xx, top], [xx, bottom]], '#accd9878', 1, [3, 3]);
  if (bx[1] - bx[0] > 75) {
    ctx.fillStyle='#c5d8bd'; ctx.textAlign='center'; ctx.font='9px '+MONO;
    ctx.fillText('INPUT RANGE', (bx[0]+bx[1])/2, top-7);
  }
  for (const id of ids) {
    const focus = reference.includes(id);
    line(ctx, points.map(p => [x(p.density), y(p.rates[id])]), dyeColor[id], focus ? 2.9 : 1.3);
  }
  if (hoverIndex !== null && points[hoverIndex]) {
    const point = points[hoverIndex], xx = x(point.density);
    line(ctx, [[xx, top], [xx, bottom]], '#e9f4e5aa', 1, [3, 2]);
    for (const id of reference) dot(ctx, xx, y(point.rates[id]), dyeColor[id], 4);
  }
  return {left, right};
}

export function drawForecast(canvas, blocks) {
  const setup = surface(canvas);
  if (!setup) return;
  const {ctx, width, height} = setup;
  const left = 53, right = width - 16, top = 21, bottom = height - 34;
  const data = blocks.flatMap(b => [b.observed, b.prediction, b.baseline]);
  const step = niceStep((Math.max(...data) - Math.min(...data)) / 4);
  const low = Math.floor(Math.min(...data) / step) * step, high = Math.max(low + step, Math.ceil(Math.max(...data) / step) * step);
  const x = i => left + i / (blocks.length - 1) * (right - left);
  const y = value => bottom - (value - low) / (high - low) * (bottom - top);
  for (let n = low; n <= high; n += step) {
    line(ctx, [[left,y(n)], [right,y(n)]], '#d7dfd6', 1, [2,3]);
    ctx.fillStyle = '#526658'; ctx.textAlign='right'; ctx.fillText(format.format(n), left - 7, y(n) + 3);
  }
  for (const [key,color,dash] of [['baseline','#6e7e94',[6,4]],['prediction','#986a25',[]],['observed','#1b714d',[]]]) {
    line(ctx, blocks.map((b,i) => [x(i),y(b[key])]),color,key === 'observed' ? 3 : 2,dash);
    if (key !== 'baseline') blocks.forEach((b,i) => dot(ctx,x(i),y(b[key]),color,3));
  }
  blocks.forEach((b,i) => {
    if (width < 500 && i % 2 === 0 && i !== 0) return;
    ctx.fillStyle='#50645a';ctx.textAlign='center';ctx.fillText(String(b.start_frame + b.frames),x(i),bottom+15);
  });
  ctx.textAlign='left';ctx.fillText('Localizations per block',left,11);
  ctx.textAlign='center';ctx.fillText('Last frame in block', (left+right)/2, height-1);
}
