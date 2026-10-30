// Headless runner for Genome-Phase-Transition.
// Loads the repo's real simulation/stat code into a VM context with a seeded RNG,
// steps it without drawing, and writes per-sample metrics plus the genealogy export.
//
// usage: node analysis/run.js --seed 1 --ticks 200000 --out runs/s1 [--genealogy] [--lineageGenomes]

const fs = require('fs');
const path = require('path');
const vm = require('vm');
const runStorage = require('./storage/local.js');

const args = {};
for (let i = 2; i < process.argv.length; i++) {
    const a = process.argv[i];
    if (a.startsWith('--')) {
        const next = process.argv[i + 1];
        if (next === undefined || next.startsWith('--')) args[a.slice(2)] = true;
        else { args[a.slice(2)] = next; i++; }
    }
}
const repo = path.resolve(args.repo || path.join(__dirname, '..'));
const seed = parseInt(args.seed || '1');
const ticks = parseInt(args.ticks || '50000');
const out = path.resolve(args.out || path.join(__dirname, 'runs', 's' + seed));
const overrides = args.params ? JSON.parse(args.params) : {};
fs.mkdirSync(out, {recursive: true});

// ── Browser stubs ────────────────────────────────────────────────────────────
const stubCtx = new Proxy({}, {get: () => () => ({width: 0})});
const blobs = [];
const sandbox = {
    console, TextEncoder, TextDecoder,
    document: {
        getElementById: () => ({textContent: '', value: '', checked: false}),
        createElement: () => ({getContext: () => stubCtx, click() {}, remove() {}, style: {}}),
        body: {appendChild() {}},
    },
    // Keeps each export's parts as written: strings for JSON, ArrayBuffers for binary.
    Blob: class { constructor(parts) { blobs.push(parts); } },
    URL: {createObjectURL: () => '', revokeObjectURL() {}},
    gameEngine: {ctx: stubCtx},
};
vm.createContext(sandbox);

const files = ['params.js', 'util.js', 'graph.js', 'histogram.js', 'views.js', 'rulespace.js', 'stats.js', 'organismlog.js', 'hexca.js'];
// Seeded RNG (mulberry32), defined in the same script as the sim so V8 can inline it.
// State lives in a closure: global vars in a VM context go through a slow interceptor.
let src = `Math.random = (function () { let s = ${seed >>> 0} ^ 0xA5A5A5A5;
    return function () { s = (s + 0x6D2B79F5) | 0;
        let t = Math.imul(s ^ (s >>> 15), 1 | s);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; })();\n`;
src += files.map(f => fs.readFileSync(path.join(repo, f), 'utf8')).join('\n;\n');
src += `
;Object.assign(PARAMETERS, ${JSON.stringify(overrides)});
globalThis.__hexca = new HexCA();
globalThis.__stats = new Stats(globalThis.__hexca);
globalThis.__PARAMETERS = PARAMETERS;`;
// One function scope makes the sim's top-level names locals instead of slow VM globals.
src = `(function () { const Math = globalThis.Math, document = globalThis.document, gameEngine = globalThis.gameEngine;\n${src}\n})();`;
vm.runInContext(src, sandbox, {filename: 'sim.js'});

const hexca = sandbox.__hexca, stats = sandbox.__stats, P = sandbox.__PARAMETERS;
const {N, ncond, n} = hexca;
const period = P.reportingPeriod;
const SYM = '.' + 'RGBYCM'.slice(0, n);   // '.' = no rule, else output colour

function genomeKey(i) {
    return Buffer.from(hexca.genomes.buffer, hexca.genomes.byteOffset + i * ncond, ncond).toString('latin1');
}
function genomeString(key) {
    let s = '';
    for (let c = 0; c < ncond; c++) { const v = key.charCodeAt(c); s += v === 255 ? '.' : SYM[v + 1]; }
    return s;
}

// The organism log grows by ~150 records per tick, so long runs can keep genealogy only
// for the first --historyTicks ticks. Founder (ancestorId) counts are unaffected.
const historyTicks = args.historyTicks !== undefined ? parseInt(args.historyTicks) : Infinity;
let historyCapped = false;

const samples = [];
const quietThreshold = Number(args.quietThreshold || 10);
const stableSamples = Number(args.stableSamples || 5);
let quietStreak = 0;
let stopReason = 'max-ticks';
const storedRun = runStorage.createRun(out);
function saveSamples() {
    fs.writeFileSync(path.join(out, 'samples.json'), JSON.stringify({seed, ticks: hexca.tick, params: P, conditions: hexca.conditions, samples, dominantGenomes}));
}
const dominantGenomes = {};   // genome string -> first tick seen as dominant
const t0 = Date.now();
for (let t = 1; t <= ticks; t++) {
    if (!historyCapped && hexca.tick >= historyTicks) {
        hexca.log.enabled = false;
        historyCapped = true;
    }
    hexca.update();
    stats.update();
    if (hexca.tick % period !== 0) continue;

    // Extra per-sample metrics the in-browser stats don't compute.
    const genomeCounts = new Map(), founderCounts = new Map();
    let alive = 0;
    for (let i = 0; i < N; i++) {
        if (hexca.color[i] === -1) continue;
        alive++;
        const k = genomeKey(i);
        genomeCounts.set(k, (genomeCounts.get(k) || 0) + 1);
        const f = hexca.ancestorId[i];
        founderCounts.set(f, (founderCounts.get(f) || 0) + 1);
    }
    const gSorted = [...genomeCounts.entries()].sort((a, b) => b[1] - a[1]);
    const genomeSnapshot = gSorted.map(([key, count]) => [genomeString(key), count]);
    const last = a => a.length ? a[a.length - 1] : null;
    const sampleBirths = last(stats.birthSeries) || 0;
    const sampleRuleDeaths = last(stats.ruleDeathSeries) || 0;
    const sampleRandomDeaths = last(stats.randomDeathSeries) || 0;
    const sampleTotalDeaths = sampleRuleDeaths + sampleRandomDeaths;
    fs.appendFileSync(snapshotPath, JSON.stringify({
        tick: hexca.tick, population: alive, distinctGenomes: genomeSnapshot.length,
        genomes: genomeSnapshot
    }) + '\\n');
    const fSorted = [...founderCounts.entries()].sort((a, b) => b[1] - a[1]);
    let simpson = 0;
    for (const [, c] of genomeCounts) simpson += (c / alive) ** 2;
    const domStr = gSorted.length ? genomeString(gSorted[0][0]) : '';
    if (domStr && !(domStr in dominantGenomes)) dominantGenomes[domStr] = hexca.tick;

    const rs = stats.rs;

    // Nearly every cell's genome is unique, so describe the population by its consensus:
    // per condition, the most common output (from rs.P, rows of [absent, colour0..]).
    const nOut = n + 1;
    let consensus = '', agree = 0;
    const Prow = [];
    for (let c = 0; c < ncond; c++) {
        let best = 0;
        for (let o = 1; o < nOut; o++) if (rs.P[c * nOut + o] > rs.P[c * nOut + best]) best = o;
        consensus += SYM[best];
        agree += rs.P[c * nOut + best];
        for (let o = 0; o < nOut; o++) Prow.push(Math.round(rs.P[c * nOut + o] * 1000));
    }
    // Mean Hamming distance of living cells from the consensus.
    let hamSum = 0;
    for (let i = 0; i < N; i++) {
        if (hexca.color[i] === -1) continue;
        const base = i * ncond;
        for (let c = 0; c < ncond; c++) if (SYM[hexca.genomes[base + c] + 1] !== consensus[c]) hamSum++;
    }
    const enc = last(rs.encFreq);
    samples.push({
        tick: hexca.tick,
        pop: alive / N,
        colors: stats.colorSeries.map(s => last(s)),
        meanGenome: last(stats.genomeSeries),
        births: last(stats.birthSeries),
        ruleDeaths: sampleRuleDeaths,
        randomDeaths: sampleRandomDeaths,
        totalDeaths: sampleTotalDeaths,
        cluster: last(stats.clusterSeries),
        meanAge: last(stats.meanAgeSeries),
        velocity: last(rs.velocitySeries),
        pi: last(rs.piSeries),
        entropy: last(rs.entropySeries),
        fixed: last(rs.fixedSeries),
        selfLoop: last(rs.selfLoopSeries),
        transMat: rs.transMat.length ? Array.from(last(rs.transMat), v => +v.toFixed(4)) : null,
        distinctGenomes: genomeCounts.size,
        domGenomeShare: gSorted.length ? gSorted[0][1] / alive : 0,
        top5GenomeShare: gSorted.slice(0, 5).reduce((s, e) => s + e[1], 0) / (alive || 1),
        effGenomes: simpson > 0 ? 1 / simpson : 0,
        domGenomeSize: domStr ? domStr.replace(/\./g, '').length : 0,
        domGenome: domStr,
        foundersAlive: founderCounts.size,
        topFounder: fSorted.length ? fSorted[0][0] : -1,
        topFounderShare: fSorted.length ? fSorted[0][1] / alive : 0,
        top3Founders: fSorted.slice(0, 3).map(([f, c]) => [f, +(c / alive).toFixed(4)]),
        organismsRecorded: hexca.log.length,
        consensus,
        consensusSize: consensus.replace(/\./g, '').length,
        consensusAgreement: agree / ncond,
        meanDistFromConsensus: alive ? hamSum / alive : 0,
        P: Prow,                                              // per-mille, ncond × (n+1)
        enc: enc ? Array.from(enc, v => +v.toFixed(5)) : null, // encounters per cell-tick
    });
    if (alive === 0) { stopReason = 'extinction'; console.log(`seed ${seed}: extinct at tick ${hexca.tick}`); break; }
    // Checkpoint so a stopped run still leaves data behind.
    if (samples.length % 10 === 0) saveSamples();
    // Stop only after births and all deaths stay below the threshold for several samples.
    quietStreak = sampleBirths < quietThreshold && sampleTotalDeaths < quietThreshold ? quietStreak + 1 : 0;
    if (quietStreak >= stableSamples) {
        stopReason = 'births-and-deaths-below-threshold';
        console.log(`seed ${seed}: births and deaths below ${quietThreshold} for ${stableSamples} samples; stopping at tick ${hexca.tick}`);
        break;
    }
    if (hexca.tick % (period * 20) === 0) {
        const s = samples[samples.length - 1];
        console.log(`seed ${seed} t=${hexca.tick} pop=${s.pop.toFixed(2)} genome=${s.meanGenome.toFixed(1)} distinct=${s.distinctGenomes} dom=${s.domGenomeShare.toFixed(2)} founders=${s.foundersAlive} ${((Date.now() - t0) / 1000).toFixed(0)}s`);
    }
}

saveSamples();
runStorage.completeRun(storedRun, {
    format: 'hex-ca-run-v1', seed, ticks: hexca.tick, stopReason,
    quietThreshold, stableSamples, params: P, conditions: hexca.conditions,
    snapshotFile: 'genome-snapshots.jsonl', sampleCount: samples.length
});

// Real genealogy export path (hexca.exportGenealogy), captured through the Blob stub.
if (args.genealogy) {
    blobs.length = 0;
    hexca.exportGenealogy();
    fs.writeFileSync(path.join(out, 'genealogy.json'), blobs[0].join(''));
}
// Lineage + genome export (hexca.exportLineageGenomes), the file lineage-explorer.html reads.
if (args.lineageGenomes) {
    blobs.length = 0;
    hexca.exportLineageGenomes(stats);
    fs.writeFileSync(path.join(out, 'lineage-genomes.hxlg'), Buffer.from(blobs[0][0]));
}
// Compact lineage table straight from the organism log (works even when the full export is huge).
const fl = new Map(), L = hexca.log;
for (let id = 0; id < L.length; id++) {
    const f = L.founder[id];
    let r = fl.get(f);
    if (!r) fl.set(f, r = {founder: f, members: 0, alive: 0, lastBirth: 0, lastDeath: 0});
    r.members++;
    if (L.death[id] === -1) r.alive++; else r.lastDeath = Math.max(r.lastDeath, L.death[id]);
    r.lastBirth = Math.max(r.lastBirth, L.birth[id]);
}
fs.writeFileSync(path.join(out, 'lineages.json'), JSON.stringify([...fl.values()]));
console.log(`seed ${seed} done: ${hexca.tick} ticks, ${hexca.log.length} organisms, ${((Date.now() - t0) / 1000).toFixed(0)}s`);
