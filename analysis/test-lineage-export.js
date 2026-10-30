// Checks the lineage + genome export end to end: runs the simulation, exports through
// hexca.exportLineageGenomes, parses the file, and verifies that every living cell's
// genome rebuilt from parent deltas matches the genome actually in the simulation, and
// that parents, founders, births and deaths agree with the family-tree JSON export.
// usage: node analysis/test-lineage-export.js
const fs = require('fs'), path = require('path'), vm = require('vm');
const repo = path.join(__dirname, '..');
const OrganismLog = require(path.join(repo, 'organismlog.js'));

function simulate(params, ticks, seed) {
    const stub = new Proxy({}, {get: () => () => ({width: 0})});
    const blobs = [];
    const sb = {
        console, TextEncoder, TextDecoder,
        document: {getElementById: () => ({textContent: ''}), createElement: () => ({getContext: () => stub, click() {}, remove() {}}), body: {appendChild() {}}},
        Blob: class { constructor(parts) { blobs.push(parts); } },
        URL: {createObjectURL: () => '', revokeObjectURL() {}},
        gameEngine: {ctx: stub},
    };
    vm.createContext(sb);
    let src = ['params.js', 'util.js', 'organismlog.js', 'hexca.js'].map(f => fs.readFileSync(path.join(repo, f), 'utf8')).join('\n;\n');
    src = `Math.random = (function () { let s = ${seed} ^ 0xA5A5A5A5; return function () { s = (s + 0x6D2B79F5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; })();\n` + src;
    src = `(function () { const Math = globalThis.Math, document = globalThis.document; ${src}
        Object.assign(PARAMETERS, ${JSON.stringify(params)});
        globalThis.h = new HexCA(); })();`;
    vm.runInContext(src, sb);
    const h = sb.h;
    for (let t = 0; t < ticks; t++) h.update();
    blobs.length = 0;
    h.exportLineageGenomes(null);
    const buf = blobs[0][0];
    blobs.length = 0;
    h.exportGenealogy();
    const json = JSON.parse(blobs[0].join(''));
    return {h, buf, json};
}

const cases = [
    {name: 'default sync, 3 colours', params: {}, ticks: 1500},
    {name: 'async updates', params: {asyncUpdate: true}, ticks: 800},
    {name: 'history rules', params: {historyRules: true}, ticks: 800},
    {name: '4 colours', params: {n: 4}, ticks: 800},
];
let failures = 0;
for (const c of cases) {
    const {h, buf, json} = simulate(c.params, c.ticks, 11);
    const d = OrganismLog.parse(buf);
    const nc = h.ncond, g = new Int8Array(nc);
    let living = 0, genomeMismatch = 0;
    for (let i = 0; i < h.N; i++) {
        if (h.color[i] === -1) continue;
        living++;
        d.genomeOf(h.organismId[i], g);
        for (let k = 0; k < nc; k++) if (g[k] !== h.genomes[i * nc + k]) { genomeMismatch++; break; }
    }
    let rowMismatch = 0;
    for (const o of json.organisms) {
        const death = o.deathTick === null ? -1 : o.deathTick;
        if (d.parent[o.id] !== o.parentId || d.founder[o.id] !== o.founderId || d.birth[o.id] !== o.birthTick || d.death[o.id] !== death) rowMismatch++;
    }
    const aliveInLog = d.header.liveAtExport.length;
    const ok = genomeMismatch === 0 && rowMismatch === 0 && living === aliveInLog && d.header.organisms === json.organisms.length;
    if (!ok) failures++;
    console.log(`${ok ? 'PASS' : 'FAIL'}  ${c.name}: ${d.header.organisms.toLocaleString()} organisms, ${living} living genomes rebuilt, ` +
        `${genomeMismatch} genome mismatches, ${rowMismatch} row mismatches, ${(buf.byteLength / 1e6).toFixed(1)} MB ` +
        `(${(d.header.deltas / d.header.organisms).toFixed(1)} changes per organism)`);
}
process.exit(failures ? 1 : 0);
