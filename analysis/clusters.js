// Genome cluster analysis across runs: who is similar to whom, without lineage data.
// From each run's final snapshot it draws organisms at random (weighted by how many carry
// each genome, seeded so it is repeatable), relabels each run so its surviving colours take
// the first labels in the order that best matches the other runs (colours.js), and measures
// Hamming distance over all conditions.
// k-medoids clusters the pooled genomes for k = 2..kmax and keeps the k with the best
// silhouette. The report says which runs fall into which cluster, the core rules nearly every
// run shares, how far apart runs are (mean distance between and within runs, raw and
// relabeled) and which rules define each cluster.
//
// usage: node analysis/clusters.js analysis/runs/c1 analysis/runs/c2 ... > clusters.json
//        [--sample 300] [--kmax 8] [--seed 1] [--present 0.08] [--tick T] [--noAlign]

const fs = require('fs');
const path = require('path');
const {readSnapshots} = require('./storage/local.js');
const {colourSpace, endpointBox} = require('./colours.js');

const args = {}, dirs = [];
for (let i = 2; i < process.argv.length; i++) {
    const a = process.argv[i];
    if (a.startsWith('--')) args[a.slice(2)] = process.argv[++i];
    else dirs.push(a);
}
const SAMPLE = parseInt(args.sample || '300');
const KMAX = parseInt(args.kmax || '8');
const PRESENT = Number(args.present || 0.08);
const TICK = args.tick ? parseInt(args.tick) : null;   // compare at this tick instead of each run's end
if (dirs.length < 2) { console.error('usage: node clusters.js <run dir> <run dir> ...'); process.exit(1); }

// Seeded RNG (mulberry32) so a sample is reproducible.
let rs = (parseInt(args.seed || '1') >>> 0) ^ 0x9E3779B9;
function rand() {
    rs = (rs + 0x6D2B79F5) | 0;
    let t = Math.imul(rs ^ (rs >>> 15), 1 | rs);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
}

const runs = dirs.map(dir => {
    const run = JSON.parse(fs.readFileSync(path.join(dir, 'run.json'), 'utf8'));
    const samples = JSON.parse(fs.readFileSync(path.join(dir, 'samples.json'), 'utf8')).samples;
    const snaps = readSnapshots(dir);
    const snap = TICK === null ? snaps[snaps.length - 1] : snaps.find(s => s.tick === TICK);
    if (!snap) { console.error(`${dir}: no snapshot at tick ${TICK}`); process.exit(1); }
    const sample = samples.find(s => s.tick === snap.tick);
    return {name: path.basename(dir), seed: run.seed, run, snap, sample};
});
const conds = runs[0].run.conditions, n = runs[0].run.params.n;
if (runs[0].run.params.historyRules) { console.error('clusters.js does not support historyRules conditions'); process.exit(1); }
for (const r of runs) if (r.run.conditions.length !== conds.length || r.run.params.n !== n)
    { console.error(`${r.name}: different condition space`); process.exit(1); }
const cs = colourSpace(conds, n), {ncond, SYM} = cs;

// Draw organisms without replacement.
const raw = [], owner = [];
for (let ri = 0; ri < runs.length; ri++) {
    const r = runs[ri];
    const box = endpointBox(r.sample, PRESENT);
    r.box = box.label;
    r.order = box.order;
    const pool = [];
    r.snap.genomes.forEach(([, count], gi) => { for (let c = 0; c < count; c++) pool.push(gi); });
    const k = Math.min(SAMPLE, pool.length);
    for (let i = 0; i < k; i++) {
        const j = i + Math.floor(rand() * (pool.length - i));
        [pool[i], pool[j]] = [pool[j], pool[i]];
        raw.push(cs.encode(r.snap.genomes[pool[i]][0]));
        owner.push(ri);
    }
    r.sampled = k;
}

// Relabel each run so its surviving colours take the first labels (R, G, …) and extinct
// ones the rest. Survivors often hold nearly equal shares, so ranking them by share is a
// coin flip; instead try every order of the survivors and keep the one whose rule profile
// is closest to the other runs', repeating until no run changes. --noAlign ranks by share.
const nOut = n + 1;
function profileOf(ri, pi) {
    const P = new Float64Array(ncond * nOut);
    for (let i = 0; i < raw.length; i++) if (owner[i] === ri) {
        const g = cs.relabelGenome(raw[i], pi);
        for (let c = 0; c < ncond; c++) P[c * nOut + g[c]]++;
    }
    for (let i = 0; i < P.length; i++) P[i] /= runs[ri].sampled;
    return P;
}
function profileDist(A, B) { let d = 0; for (let i = 0; i < A.length; i++) d += Math.abs(A[i] - B[i]); return d; }
for (const r of runs) {
    const m = r.order.length, present = new Set(r.order);
    r.candidates = cs.perms.map((p, pi) => pi).filter(pi => cs.perms[pi].every((v, k) => present.has(k) ? v < m : v >= m));
    r.pi = cs.canonicalPerm(r.order);
}
if (!args.noAlign) {
    const prof = runs.map((r, ri) => profileOf(ri, r.pi));
    for (let pass = 0; pass < 10; pass++) {
        let changed = false;
        for (let ri = 0; ri < runs.length; ri++) {
            const r = runs[ri];
            let best = r.pi, bd = Infinity, bp = prof[ri];
            for (const pi of r.candidates) {
                const P = profileOf(ri, pi);
                let d = 0;
                for (let o = 0; o < runs.length; o++) if (o !== ri) d += profileDist(P, prof[o]);
                if (d < bd - 1e-9) { bd = d; best = pi; bp = P; }
            }
            if (best !== r.pi) { r.pi = best; prof[ri] = bp; changed = true; }
        }
        if (!changed) break;
    }
}
const canon = raw.map((g, i) => cs.relabelGenome(g, runs[owner[i]].pi));
for (const r of runs) r.perm = cs.permName(r.pi);
const N = canon.length, R = runs.length;

function hamming(a, b) { let d = 0; for (let c = 0; c < ncond; c++) if (a[c] !== b[c]) d++; return d; }
function matrix(G) {
    const D = new Uint16Array(N * N);
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) D[i * N + j] = D[j * N + i] = hamming(G[i], G[j]);
    return D;
}
const Dc = matrix(canon), Dr = matrix(raw);

// Mean distance between (and within, on the diagonal) runs.
function runDistances(D) {
    const sum = Array.from({length: R}, () => new Array(R).fill(0)), cnt = Array.from({length: R}, () => new Array(R).fill(0));
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
        const a = owner[i], b = owner[j];
        sum[a][b] += D[i * N + j]; cnt[a][b]++;
        if (a !== b) { sum[b][a] += D[i * N + j]; cnt[b][a]++; }
    }
    return sum.map((row, a) => row.map((s, b) => cnt[a][b] ? +(s / cnt[a][b]).toFixed(1) : null));
}

// k-medoids (alternating assignment / medoid update) with k-medoids++ starts.
function kmedoids(D, k, restarts = 6) {
    let best = null;
    for (let rep = 0; rep < restarts; rep++) {
        const med = [Math.floor(rand() * N)];
        const near = new Float64Array(N).fill(Infinity);
        while (med.length < k) {
            let tot = 0;
            for (let i = 0; i < N; i++) { near[i] = Math.min(near[i], D[i * N + med[med.length - 1]]); tot += near[i] ** 2; }
            let x = rand() * tot, pick = 0;
            for (let i = 0; i < N; i++) { x -= near[i] ** 2; if (x <= 0) { pick = i; break; } }
            med.push(pick);
        }
        const lab = new Int32Array(N);
        for (let iter = 0; iter < 50; iter++) {
            for (let i = 0; i < N; i++) {
                let bd = Infinity;
                for (let c = 0; c < k; c++) { const d = D[i * N + med[c]]; if (d < bd) { bd = d; lab[i] = c; } }
            }
            let changed = false;
            for (let c = 0; c < k; c++) {
                const mem = [];
                for (let i = 0; i < N; i++) if (lab[i] === c) mem.push(i);
                if (!mem.length) continue;
                let bm = med[c], bs = Infinity;
                for (const i of mem) { let s = 0; for (const j of mem) s += D[i * N + j]; if (s < bs) { bs = s; bm = i; } }
                if (bm !== med[c]) { med[c] = bm; changed = true; }
            }
            if (!changed) break;
        }
        let cost = 0;
        for (let i = 0; i < N; i++) cost += D[i * N + med[lab[i]]];
        if (!best || cost < best.cost) best = {cost, med: med.slice(), lab: Int32Array.from(lab)};
    }
    return best;
}
function silhouette(D, lab, k) {
    const size = new Array(k).fill(0);
    for (let i = 0; i < N; i++) size[lab[i]]++;
    let total = 0;
    const sums = new Float64Array(k);
    for (let i = 0; i < N; i++) {
        sums.fill(0);
        for (let j = 0; j < N; j++) if (j !== i) sums[lab[j]] += D[i * N + j];
        const own = lab[i];
        if (size[own] <= 1) continue;
        const a = sums[own] / (size[own] - 1);
        let b = Infinity;
        for (let c = 0; c < k; c++) if (c !== own && size[c]) b = Math.min(b, sums[c] / size[c]);
        total += (b - a) / Math.max(a, b);
    }
    return total / N;
}

const trials = [];
for (let k = 2; k <= Math.min(KMAX, N - 1); k++) {
    const km = kmedoids(Dc, k);
    trials.push({k, km, sil: silhouette(Dc, km.lab, k)});
}
const chosen = trials.reduce((a, b) => b.sil > a.sil ? b : a);
const {k, km} = chosen;

// Order clusters by size.
const order = [...Array(k).keys()].sort((a, b) => km.lab.filter(l => l === b).length - km.lab.filter(l => l === a).length);
const rank = new Map(order.map((c, i) => [c, i]));
const clusters = order.map(c => {
    const mem = [];
    for (let i = 0; i < N; i++) if (km.lab[i] === c) mem.push(i);
    const perRun = runs.map((r, ri) => +(mem.filter(i => owner[i] === ri).length / r.sampled).toFixed(3));
    let spread = 0;
    for (const i of mem) spread += Dc[i * N + km.med[c]];
    return {cluster: String.fromCharCode(65 + rank.get(c)), size: mem.length,
        meanDistToMedoid: +(spread / mem.length).toFixed(1),
        runs: Object.fromEntries(runs.map((r, ri) => [r.name, perRun[ri]]).filter(([, v]) => v > 0)),
        medoid: cs.decode(canon[km.med[c]]), medoidRules: canon[km.med[c]].filter(v => v).length, mem};
});

// Defining rules: carried by ≥80% of a cluster and ≤30% of everyone else.
for (const cl of clusters) {
    const inSet = new Set(cl.mem), others = N - cl.mem.length;
    const rules = [];
    for (let c = 0; c < ncond; c++) for (let o = 1; o <= n; o++) {
        let inside = 0, outside = 0;
        for (let i = 0; i < N; i++) if (canon[i][c] === o) { if (inSet.has(i)) inside++; else outside++; }
        const fin = inside / cl.mem.length, fout = others ? outside / others : 0;
        if (fin >= 0.8 && fout <= 0.3) rules.push({text: cs.ruleText(c, SYM[o - 1]), inCluster: +fin.toFixed(2), elsewhere: +fout.toFixed(2)});
    }
    rules.sort((a, b) => (b.inCluster - b.elsewhere) - (a.inCluster - a.elsewhere));
    cl.definingRules = rules;
    delete cl.mem;
}

// Core rulebook: rules most organisms carry in nearly every run, after relabeling.
const coreRules = [];
for (let c = 0; c < ncond; c++) for (let o = 1; o <= n; o++) {
    const shares = runs.map(() => 0);
    for (let i = 0; i < N; i++) if (canon[i][c] === o) shares[owner[i]]++;
    shares.forEach((v, ri) => { shares[ri] = v / runs[ri].sampled; });
    const held = shares.filter(v => v >= 0.5).length;
    if (held >= Math.max(2, R - 1)) coreRules.push({text: cs.ruleText(c, SYM[o - 1]), runs: held,
        meanShare: +(shares.reduce((x, y) => x + y, 0) / R).toFixed(2)});
}
coreRules.sort((x, y) => y.runs - x.runs || y.meanShare - x.meanShare);

// Runs whose sampled organisms mostly share a cluster belong together.
const runSummary = runs.map((r, ri) => {
    const shares = clusters.map(cl => [cl.cluster, cl.runs[r.name] || 0]).sort((a, b) => b[1] - a[1]);
    return {name: r.name, seed: r.seed, tick: r.snap.tick, stopReason: r.run.stopReason, box: r.box, relabel: r.perm,
        sampled: r.sampled, distinctGenomes: r.snap.distinctGenomes, mainCluster: shares[0][0], mainShare: shares[0][1]};
});
const runDistCanon = runDistances(Dc), runDistRaw = runDistances(Dr);

const result = {
    method: {sample: SAMPLE, distance: 'Hamming over all ' + ncond + ' conditions (no rule counts as a value)',
        relabel: args.noAlign ? 'each run relabeled so its colours, most common first, become R, G, B' : 'surviving colours take the first labels, ordered to best match the other runs; extinct colours take the rest', kmax: KMAX, seed: parseInt(args.seed || '1')},
    k, silhouette: trials.map(t => ({k: t.k, silhouette: +t.sil.toFixed(3)})),
    runs: runSummary, clusters, coreRules,
    runDistance: {names: runs.map(r => r.name), relabeled: runDistCanon, raw: runDistRaw},
};
process.stdout.write(JSON.stringify(result, null, 1));

const pad = (s, w) => String(s).padEnd(w);
console.error(`\n${N} genomes from ${R} runs; best k = ${k} (silhouette ${trials.map(t => t.k + ':' + t.sil.toFixed(2)).join(' ')})`);
console.error('\nRuns');
for (const r of runSummary) console.error(`  ${pad(r.name, 5)} ${pad(r.box, 12)} tick ${pad(r.tick, 7)} ${pad(r.stopReason, 34)} cluster ${r.mainCluster} (${(r.mainShare * 100).toFixed(0)}%)`);
console.error('\nClusters');
for (const cl of clusters) console.error(`  ${cl.cluster}  ${pad(cl.size, 5)} spread ${pad(cl.meanDistToMedoid, 5)} rules ${pad(cl.medoidRules, 4)} ${Object.entries(cl.runs).map(([k, v]) => k + ' ' + (v * 100).toFixed(0) + '%').join(', ')}`);
console.error('\nMean Hamming distance between runs (relabeled; diagonal = within run)');
console.error('        ' + runs.map(r => pad(r.name, 6)).join(''));
runDistCanon.forEach((row, a) => console.error('  ' + pad(runs[a].name, 6) + row.map(v => pad(v, 6)).join('')));
console.error(`
Core rules (held by ≥50% of organisms in at least ${Math.max(2, R - 1)} of ${R} runs): ${coreRules.length}`);
for (const r of coreRules.slice(0, 15)) console.error(`  ${r.runs}/${R}  ${r.text}   mean ${r.meanShare}`);
for (const cl of clusters) if (cl.definingRules.length) {
    console.error(`\nCluster ${cl.cluster} defining rules (${cl.definingRules.length})`);
    for (const r of cl.definingRules.slice(0, 6)) console.error(`  ${r.text}   ${r.inCluster} vs ${r.elsewhere}`);
}
