// Measures how compressible per-organism genomes are: distinct genomes across all births,
// and how many conditions differ between child and parent.
// usage: node analysis/measure-genomes.js . 10000 7
const fs = require('fs'), path = require('path'), vm = require('vm');
const repo = path.resolve(process.argv[2]), ticks = +process.argv[3] || 10000, seed = +process.argv[4] || 7;
const stub = new Proxy({}, {get: () => () => ({width: 0})});
const sb = {console, TextEncoder, TextDecoder, document: {getElementById: () => ({textContent: ''}), createElement: () => ({getContext: () => stub}), body: {appendChild() {}}}, gameEngine: {ctx: stub}};
vm.createContext(sb);
let src = ['params.js', 'util.js', 'organismlog.js', 'hexca.js'].map(f => fs.readFileSync(path.join(repo, f), 'utf8')).join('\n;\n');
// Hook: after each birth, report (child cell, parent cell). The parent's genome is unchanged by reproduction.
const hook = 'birthTick[targetKey]  = this.tick;';
if (!src.includes(hook)) throw new Error('birth hook point not found in hexca.js');
src = src.replace(hook, hook + ' if (this.onBirth) this.onBirth(targetKey, i);');
src = `Math.random = (function () { let s = ${seed} ^ 0xA5A5A5A5; return function () { s = (s + 0x6D2B79F5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; })();\n` + src;
src = `(function () { const Math = globalThis.Math, document = globalThis.document; ${src}\n globalThis.h = new HexCA(); })();`;
vm.runInContext(src, sb);

const h = sb.h, nc = h.ncond, G = h.genomes;
const key = i => Buffer.from(G.buffer, G.byteOffset + i * nc, nc).toString('latin1');
const seen = new Map();   // genome -> number of organisms carrying it
let births = 0, diffSum = 0, zeroDiff = 0;
const diffHist = new Array(16).fill(0);
for (let i = 0; i < h.N; i++) if (h.color[i] !== -1) { const k = key(i); seen.set(k, (seen.get(k) || 0) + 1); }
const founders = [...seen.values()].reduce((a, b) => a + b, 0);
h.onBirth = (child, parent) => {
    births++;
    const k = key(child);
    seen.set(k, (seen.get(k) || 0) + 1);
    let d = 0;
    for (let c = 0; c < nc; c++) if (G[child * nc + c] !== G[parent * nc + c]) d++;
    diffSum += d; if (d === 0) zeroDiff++; diffHist[Math.min(15, d)]++;
};
const checkpoints = [];
for (let t = 1; t <= ticks; t++) {
    h.update();
    if (t % (ticks / 5) === 0) checkpoints.push({tick: t, organisms: founders + births, distinct: seen.size});
}
const counts = [...seen.values()].sort((a, b) => b - a);
console.log(JSON.stringify({
    ticks, seed, founders, births, organisms: founders + births, distinctGenomes: seen.size,
    distinctShare: +(seen.size / (founders + births)).toFixed(4),
    organismsInTop100Genomes: counts.slice(0, 100).reduce((a, b) => a + b, 0),
    largestGenomeGroup: counts[0], singletonGenomes: counts.filter(c => c === 1).length,
    meanChangesPerBirth: +(diffSum / births).toFixed(2), identicalToParent: +(zeroDiff / births).toFixed(4),
    changesHistogram: diffHist, checkpoints,
}, null, 1));
