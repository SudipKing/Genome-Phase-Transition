// Do independent runs end in the same genetic place, up to relabeling the colours?
// Unlike symmetry.js (which compares one consensus genome per run), this compares the
// actual final genome populations from genome-snapshots.jsonl(.gz).
//
// Each population is summarised by its rule profile: for every condition, the share of
// organisms whose genome has no rule ('.') or outputs R, G or B there. Two profiles are
// compared per condition by total-variation distance, weighted by how much either
// population actually uses that condition (share with a rule), so the many conditions
// neither population has a rule for don't inflate the score. similarity = 1 - weighted TV.
// Every pair is scored raw and under all colour permutations (n! = 6 for three colours);
// the permutation is applied to the condition's self colour, its neighbour counts and the
// output, so a rule maps to the matching rule of the relabeled world.
//
// usage: node analysis/endpoints.js analysis/runs/e1 analysis/runs/e2 ... > endpoints.json
//        [--present 0.08] [--cluster 0.75]

const fs = require('fs');
const path = require('path');
const {readSnapshots} = require('./storage/local.js');

const args = {}, dirs = [];
for (let i = 2; i < process.argv.length; i++) {
    const a = process.argv[i];
    if (a.startsWith('--')) args[a.slice(2)] = process.argv[++i];
    else dirs.push(a);
}
const PRESENT = Number(args.present || 0.08);   // same colour-presence cut as analyze.js
const CLUSTER = Number(args.cluster || 0.75);   // similarity needed to join a cluster
if (dirs.length < 2) { console.error('usage: node endpoints.js <run dir> <run dir> ...'); process.exit(1); }

const SYM = 'RGBYCM';

function load(dir) {
    const run = JSON.parse(fs.readFileSync(path.join(dir, 'run.json'), 'utf8'));
    const samples = JSON.parse(fs.readFileSync(path.join(dir, 'samples.json'), 'utf8')).samples;
    const snaps = readSnapshots(dir);
    return {name: path.basename(dir), seed: run.seed, run, last: samples[samples.length - 1], snaps};
}
const runs = dirs.map(load);
const conds = runs[0].run.conditions;            // [self, count colour 0, ..., count colour n-1]
const n = runs[0].run.params.n, ncond = conds.length, nOut = n + 1;
if (runs[0].run.params.historyRules) { console.error('endpoints.js does not support historyRules conditions'); process.exit(1); }
for (const r of runs) if (r.run.conditions.length !== ncond || r.run.params.n !== n)
    { console.error(`${r.name}: different condition space`); process.exit(1); }

const condIndex = new Map(conds.map((c, i) => [c.join(','), i]));
const outIndex = ch => ch === '.' ? 0 : SYM.indexOf(ch) + 1;

// Rule profile: Float64Array(ncond × (n+1)), each condition's row sums to 1.
function profile(snap) {
    const P = new Float64Array(ncond * nOut);
    let total = 0;
    for (const [g, count] of snap.genomes) {
        total += count;
        for (let c = 0; c < ncond; c++) P[c * nOut + outIndex(g[c])] += count;
    }
    for (let i = 0; i < P.length; i++) P[i] /= total || 1;
    return P;
}

function permutations(k) {
    if (k === 1) return [[0]];
    const out = [];
    for (const p of permutations(k - 1))
        for (let i = 0; i <= p.length; i++) out.push([...p.slice(0, i), k - 1, ...p.slice(i)]);
    return out;
}
const perms = permutations(n).sort((a, b) => a.join('').localeCompare(b.join('')));
const permName = p => p.map(k => SYM[k]).join('');   // "BGR": R→B, G→G, B→R

// Condition map under colour permutation p (old colour k → p[k]).
const condMaps = perms.map(p => conds.map(([s, ...cnt]) => {
    const m = new Array(n).fill(0);
    for (let k = 0; k < n; k++) m[p[k]] = cnt[k];
    return condIndex.get([p[s], ...m].join(','));
}));
function permuteProfile(P, pi) {
    const p = perms[pi], cm = condMaps[pi], Q = new Float64Array(P.length);
    for (let c = 0; c < ncond; c++) {
        const j = cm[c];
        Q[j * nOut] = P[c * nOut];
        for (let k = 0; k < n; k++) Q[j * nOut + 1 + p[k]] = P[c * nOut + 1 + k];
    }
    return Q;
}
function similarity(A, B) {
    let num = 0, den = 0;
    for (let c = 0; c < ncond; c++) {
        const w = Math.max(1 - A[c * nOut], 1 - B[c * nOut]);
        if (w <= 0) continue;
        let tv = 0;
        for (let o = 0; o < nOut; o++) tv += Math.abs(A[c * nOut + o] - B[c * nOut + o]);
        num += w * tv / 2; den += w;
    }
    return den ? 1 - num / den : 1;
}
// Best similarity of A to B over every relabeling of A.
function bestMatch(A, B) {
    let best = {sim: -1, perm: ''};
    for (let pi = 0; pi < perms.length; pi++) {
        const s = similarity(permuteProfile(A, pi), B);
        if (s > best.sim) best = {sim: s, perm: permName(perms[pi])};
    }
    return best;
}

// Endpoint "box": colours present above PRESENT, most common first.
function box(sample) {
    if (sample.pop < 0.01) return {label: 'extinct', order: []};
    const order = sample.colors.map((f, k) => [f, k]).filter(x => x[0] > PRESENT).sort((a, b) => b[0] - a[0]).map(x => x[1]);
    return {label: order.map(k => SYM[k]).join('>') + (sample.pop >= 0.95 ? ' full' : ''), order};
}
// Relabel so the dominant colour is R, the next G, and so on; absent colours fill the rest.
function canonicalPerm(order) {
    const p = new Array(n), used = new Set();
    order.forEach((k, rank) => { p[k] = rank; used.add(k); });
    let next = order.length;
    for (let k = 0; k < n; k++) if (!used.has(k)) p[k] = next++;
    return perms.findIndex(q => q.every((v, i) => v === p[i]));
}

const finals = runs.map(r => {
    const snap = r.snaps[r.snaps.length - 1];
    const b = box(r.last);
    const P = profile(snap);
    return {name: r.name, seed: r.seed, tick: snap.tick, stopReason: r.run.stopReason,
        pop: +r.last.pop.toFixed(3), colors: r.last.colors.map(c => +c.toFixed(3)),
        births: +r.last.births.toFixed(1), totalDeaths: +r.last.totalDeaths.toFixed(1),
        distinctGenomes: snap.distinctGenomes, box: b.label, order: b.order,
        P, canon: permuteProfile(P, canonicalPerm(b.order)),
        genomes: new Set(snap.genomes.map(g => g[0]))};
});

// Pairwise scores.
const pairs = [], N = finals.length;
const S = Array.from({length: N}, () => new Array(N).fill(1));
for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
    const a = finals[i], b = finals[j];
    const raw = similarity(a.P, b.P), best = bestMatch(a.P, b.P);
    let shared = 0;
    for (const g of a.genomes) if (b.genomes.has(g)) shared++;
    S[i][j] = S[j][i] = best.sim;
    pairs.push({a: a.name, b: b.name, boxA: a.box, boxB: b.box,
        raw: +raw.toFixed(3), relabeled: +best.sim.toFixed(3), perm: best.perm,
        canonical: +similarity(a.canon, b.canon).toFixed(3), sharedGenomes: shared});
}

// Average-linkage clustering on the relabeled similarity.
let clusters = finals.map((_, i) => [i]);
const link = (A, B) => { let s = 0; for (const i of A) for (const j of B) s += S[i][j]; return s / (A.length * B.length); };
for (;;) {
    let best = null;
    for (let x = 0; x < clusters.length; x++) for (let y = x + 1; y < clusters.length; y++) {
        const s = link(clusters[x], clusters[y]);
        if (s >= CLUSTER && (!best || s > best.s)) best = {x, y, s};
    }
    if (!best) break;
    clusters[best.x] = clusters[best.x].concat(clusters[best.y]);
    clusters.splice(best.y, 1);
}
const clusterOut = clusters.map(c => ({runs: c.map(i => finals[i].name), boxes: c.map(i => finals[i].box),
    cohesion: c.length > 1 ? +link(c, c).toFixed(3) : null}));

// When did each run reach its endpoint? First sample from which every later sample stays
// within 0.9 similarity of the final population.
const convergence = runs.map((r, i) => {
    const sims = r.snaps.map(s => similarity(profile(s), finals[i].P));
    let k = sims.length - 1;
    while (k > 0 && sims[k - 1] >= 0.9) k--;
    return {name: r.name, convergedAt: r.snaps[k].tick, trajectory: sims.map((s, j) => [r.snaps[j].tick, +s.toFixed(3)])};
});

// Rules most runs share once every endpoint is put in canonical colours (dominant = R).
function ruleText(c, out) {
    const [s, ...cnt] = conds[c];
    const nb = cnt.map((v, k) => v ? v + SYM[k] : '').filter(Boolean);
    const empty = 6 - cnt.reduce((x, y) => x + y, 0);
    if (empty) nb.push(empty + ' empty');
    return `${SYM[s]} with ${nb.join(' ')} → ${out}`;
}
const shared = [];
for (let c = 0; c < ncond; c++) for (let o = 1; o < nOut; o++) {
    const shares = finals.map(f => f.canon[c * nOut + o]);
    const carriers = shares.filter(v => v >= 0.5).length;
    if (carriers >= 2) shared.push({cond: c, out: SYM[o - 1], runs: carriers,
        meanShare: +(shares.reduce((x, y) => x + y, 0) / N).toFixed(3), text: ruleText(c, SYM[o - 1])});
}
shared.sort((x, y) => y.runs - x.runs || y.meanShare - x.meanShare);

const result = {
    method: {present: PRESENT, cluster: CLUSTER, similarity: '1 - rule-use-weighted total variation of per-condition output shares'},
    runs: finals.map(({P, canon, genomes, order, ...f}) => f),
    pairs, clusters: clusterOut, convergence, sharedRules: shared,
};
process.stdout.write(JSON.stringify(result, null, 1));

// Readable summary on stderr.
const pad = (s, w) => String(s).padEnd(w);
console.error('\nEndpoints');
for (const f of result.runs) console.error(`  ${pad(f.name, 6)} ${pad(f.box, 12)} pop ${f.pop}  colours ${f.colors.join('/')}  distinct ${f.distinctGenomes}  ${f.stopReason} @ ${f.tick}`);
console.error('\nPairs                 raw   relabeled  perm  canonical  shared genomes');
for (const p of pairs) console.error(`  ${pad(p.a + '–' + p.b, 18)} ${pad(p.raw, 6)} ${pad(p.relabeled, 10)} ${pad(p.perm, 5)} ${pad(p.canonical, 10)} ${p.sharedGenomes}`);
console.error(`\nClusters (relabeled similarity ≥ ${CLUSTER})`);
for (const c of clusterOut) console.error(`  [${c.runs.join(', ')}]  ${c.boxes.join(' | ')}${c.cohesion !== null ? '  cohesion ' + c.cohesion : ''}`);
console.error('\nConverged (within 0.9 of final from)');
for (const c of convergence) console.error(`  ${pad(c.name, 6)} tick ${c.convergedAt}`);
console.error(`\nRules held by ≥50% of organisms in ≥2 runs (canonical colours): ${shared.length}`);
for (const r of shared.slice(0, 12)) console.error(`  ${r.runs} runs  ${r.text}`);
