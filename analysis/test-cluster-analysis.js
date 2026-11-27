// Checks for the cluster analysis core (GA) in genome-cluster-explorer.html, <script id="cluster-analysis-core">:
// distances, rule conservation, compactness, lineage relations, rank statistics and cluster comparison.
// Every expected number below is worked out by hand in the comments next to it.
// usage: node analysis/test-cluster-analysis.js
const assert = require('assert');
const fs = require('fs');
const path = require('path');
const gc = require('./genome-clusters.js');

function loadGA() {
    const html = fs.readFileSync(path.join(__dirname, '..', 'genome-cluster-explorer.html'), 'utf8');
    const m = html.match(/<script id="cluster-analysis-core">([\s\S]*?)<\/script>/);
    if (!m) throw new Error('cluster-analysis-core block not found');
    return new Function('module', m[1] + '\nreturn GA;')({});
}
const GA = loadGA();

let rs = 424242;
const rand = () => { rs = (rs * 1103515245 + 12345) >>> 0; return rs / 4294967296; };
const OUT = '.RGB';
const pad = s => s + '.'.repeat(252 - s.length);
const randomGenome = () => Array.from({length: 252}, () => OUT[Math.floor(rand() * 4)]).join('');
function mutate(g, k) {   // k random changes: toggle presence or change colour
    const a = [...g];
    for (let i = 0; i < k; i++) {
        const c = Math.floor(rand() * a.length), r = rand();
        a[c] = r < 0.4 ? (a[c] === '.' ? 'RGB'[Math.floor(rand() * 3)] : '.') : 'RGB'[Math.floor(rand() * 3)];
    }
    return a.join('');
}
const close = (a, b, msg) => assert.ok(Math.abs(a - b) < 1e-9, `${msg}: ${a} != ${b}`);
function noNaN(v, where = 'value') {   // walks objects/arrays; fails on any NaN
    if (typeof v === 'number') assert.ok(!Number.isNaN(v), `NaN at ${where}`);
    else if (v && typeof v === 'object') for (const k of Object.keys(v)) noNaN(v[k], where + '.' + k);
}
const E = (g, count, ids = null, pat = 0) => ({g, count, ids, pat});

let passed = 0;
const ok = (name) => { passed++; console.log('ok  ' + name); };

// ================================================================ 1. distances
{
    // a: rules at 0 (R), 2 (G), 4 (B); b: rules at 0 (G), 2 (G), 5 (R).
    // presence: 4 (a only) and 5 (b only) → 2. shared: 0, 2; output differs at 0 (R vs G) → 1. fine = 3.
    const a = pad('R.G.B'), b = pad('G.G..R');
    assert.deepStrictEqual(GA.genomeDistance(a, b), {presence: 2, output: 1, fine: 3, shared: 2});
    assert.deepStrictEqual(GA.pairDist(GA.encode(a, 3), GA.encode(b, 3), 3), {presence: 2, output: 1, fine: 3, shared: 2});
    // identical genomes: all zero, shared = 3 rules
    assert.deepStrictEqual(GA.pairDist(GA.encode(a, 3), GA.encode(a, 3), 3), {presence: 0, output: 0, fine: 0, shared: 3});
    // empty genome vs a: presence = a's 3 rules, nothing shared
    const z = '.'.repeat(252);
    assert.deepStrictEqual(GA.pairDist(GA.encode(z, 3), GA.encode(a, 3), 3), {presence: 3, output: 0, fine: 3, shared: 0});
    // word boundaries: conditions 31, 32 and 251 (last). c: R@31, G@32, B@251; d: G@31, G@32, nothing @251
    // presence 1 (251); shared 31, 32; output 1 (31) → fine 2
    const c = [...z], d = [...z]; c[31] = 'R'; c[32] = 'G'; c[251] = 'B'; d[31] = 'G'; d[32] = 'G';
    assert.deepStrictEqual(GA.pairDist(GA.encode(c.join(''), 3), GA.encode(d.join(''), 3), 3), {presence: 1, output: 1, fine: 2, shared: 2});
    ok('distances: hand examples (incl. word boundaries 31/32/251)');

    for (let i = 0; i < 500; i++) {
        const g = randomGenome(), h = rand() < 0.5 ? mutate(g, Math.floor(rand() * 40)) : randomGenome();
        const want = GA.genomeDistance(g, h), got = GA.pairDist(GA.encode(g, 3), GA.encode(h, 3), 3);
        assert.deepStrictEqual(got, want);
        assert.strictEqual(want.presence, gc.coarseDistance(g, h));
        assert.strictEqual(want.presence, gc.bitDistance(gc.presenceBits(g), gc.presenceBits(h)));
        assert.strictEqual(want.fine, gc.fineDistance(g, h));
        const cmp = gc.compareGenomes(g, h);
        assert.strictEqual(want.output, cmp.differentOutputs);
        assert.strictEqual(want.shared, cmp.sharedRules);
    }
    ok('distances: bit version == string version == genome-clusters.js coarse/fine on 500 random pairs');
}

// ================================================================ 2. conservation
{
    // N = 4 organisms, 8 conditions, 3 colours.
    //   idx:  0 1 2 3 4 5 6 7
    //   e1:   R R R . G B . .   ×2
    //   e2:   R R G . . B . R   ×1
    //   e3:   R . B . . G G G   ×1
    const C = GA.conservation([E('RRR.GB..', 2), E('RRG..B.R', 1), E('R.B..GGG', 1)], 8, 3);
    const want = [
        // x  carriers outs     top topCount tie   agreement unanimous level
        [0, 4, [4, 0, 0], 0, 4, false, 1, true, 'all'],
        [1, 3, [3, 0, 0], 0, 3, false, 1, true, 'most'],
        [2, 4, [2, 1, 1], 0, 2, false, 0.5, false, 'all'],     // present in all, outputs differ
        [3, 0, [0, 0, 0], -1, 0, false, null, false, 'absent'],
        [4, 2, [0, 2, 0], 1, 2, false, 1, true, 'half'],
        [5, 4, [0, 1, 3], 2, 3, false, 0.75, false, 'all'],
        [6, 1, [0, 1, 0], 1, 1, false, 1, true, 'minority'],
        [7, 2, [1, 1, 0], 0, 1, true, 0.5, false, 'half'],     // R vs G tie
    ];
    assert.strictEqual(C.N, 4);
    for (const [x, carriers, outs, top, topCount, tie, agreement, unanimous, level] of want) {
        const r = C.rows[x];
        assert.deepStrictEqual({x: r.x, carriers: r.carriers, outs: r.outs, top: r.top, topCount: r.topCount, tie: r.tie,
            agreement: r.agreement, unanimous: r.unanimous, level: r.level},
            {x, carriers, outs, top, topCount, tie, agreement, unanimous, level}, `row ${x}`);
        assert.strictEqual(r.prevalence, carriers / 4);
    }
    ok('conservation: carriers/prevalence/outs/agreement/unanimous/top/tie per condition (N=4)');

    const s = C.summary;
    // levels: absent x3; minority x6; half x4,x7; most x1; all x0,x2,x5
    assert.deepStrictEqual(s.levels, {absent: 1, minority: 1, half: 2, most: 1, all: 3});
    // atLeast 100: x0,x2,x5; 75: + x1; 50: + x4,x7; 25: + x6
    assert.deepStrictEqual(s.atLeast, {100: 3, 75: 4, 50: 6, 25: 7});
    assert.strictEqual(s.presentAny, 7);
    assert.strictEqual(s.unanimous, 4);              // x0 x1 x4 x6
    assert.strictEqual(s.partialDisagreement, 3);    // x2 x5 x7 (x2 and x5 are 'all' yet disagree)
    assert.deepStrictEqual(s.unanimousAt, {100: 1, 75: 2, 50: 3});   // {x0}, {x0,x1}, {x0,x1,x4}
    assert.strictEqual(s.fullyConserved, 1);         // only x0 is in all members with one output
    // presenceAgreement = mean over present x of max(c, N-c)/N = (1 + .75 + 1 + .5 + 1 + .75 + .5) / 7 = 5.5/7
    close(s.presenceAgreement, 5.5 / 7, 'presenceAgreement');
    // outputAgreement = sum topCount / sum carriers = (4+3+2+2+3+1+1) / (4+3+4+2+4+1+2) = 16/20
    close(s.outputAgreement, 0.8, 'outputAgreement');
    ok('conservation: summary counts, presenceAgreement 5.5/7, outputAgreement 0.8');
    // presence is not conflated with output: x2 is 'all' and counts as partial disagreement
    assert.strictEqual(C.rows[2].level, 'all');
    assert.ok(!C.rows[2].unanimous);
    ok('conservation: a rule in every member with differing outputs is level all AND partialDisagreement');
}
{
    // Level boundaries. N = 4: organism i has a rule at x iff i < x, so carriers(x) = x.
    //   0 absent, 1 = 25% minority, 2 = 50% half, 3 = 75% most, 4 = 100% all
    const gen = (i, n) => Array.from({length: n + 1}, (_, x) => i < x ? 'R' : '.').join('');
    const C4 = GA.conservation([0, 1, 2, 3].map(i => E(gen(i, 4), 1)), 5, 3);
    assert.deepStrictEqual(C4.rows.map(r => r.level), ['absent', 'minority', 'half', 'most', 'all']);
    assert.deepStrictEqual(C4.summary.atLeast, {100: 1, 75: 2, 50: 3, 25: 4});
    // N = 8: carriers 0..8 → absent, min(1/8), min(2/8=25%), min(3/8), half(4/8), half(5/8), most(6/8), most(7/8), all
    const C8 = GA.conservation([0, 1, 2, 3, 4, 5, 6, 7].map(i => E(gen(i, 8), 1)), 9, 3);
    assert.deepStrictEqual(C8.rows.map(r => r.level), ['absent', 'minority', 'minority', 'minority', 'half', 'half', 'most', 'most', 'all']);
    // atLeast 100: {8}; 75: {6,7,8}; 50: {4..8}; 25: {2..8}
    assert.deepStrictEqual(C8.summary.atLeast, {100: 1, 75: 3, 50: 5, 25: 7});
    // same N = 8 as four entries of count 2 (entry j stands for organisms 2j, 2j+1): carriers(x) = 2*ceil(x/2)
    const C8w = GA.conservation([0, 1, 2, 3].map(j => E(gen(2 * j, 8), 2)), 9, 3);
    assert.deepStrictEqual(C8w.rows.map(r => r.carriers), [0, 2, 2, 4, 4, 6, 6, 8, 8]);
    assert.deepStrictEqual(C8w.rows.map(r => r.level), ['absent', 'minority', 'minority', 'half', 'half', 'most', 'most', 'all', 'all']);
    assert.deepStrictEqual([1, 2, 3, 4].map(c => GA.levelOf(c, 4)), ['minority', 'half', 'most', 'all']);
    assert.deepStrictEqual([0, 2, 4, 6, 8].map(c => GA.levelOf(c, 8)), ['absent', 'minority', 'half', 'most', 'all']);
    ok('conservation: level boundaries exactly at 25/50/75/100% for N=4 and N=8 (weighted and unweighted)');
}
{
    // single organism 'R.G': both rules 'all', unanimous; presenceAgreement 1, outputAgreement 1
    const C = GA.conservation([E('R.G', 1)], 3, 3);
    assert.deepStrictEqual(C.rows.map(r => r.level), ['all', 'absent', 'all']);
    assert.strictEqual(C.summary.fullyConserved, 2);
    assert.strictEqual(C.summary.presenceAgreement, 1);
    assert.strictEqual(C.summary.outputAgreement, 1);
    ok('conservation: single-organism cluster');
    // identical genomes: one entry ×5 == entries ×2 and ×3 of the same genome
    const A = GA.conservation([E('RG.B', 5)], 4, 3), B = GA.conservation([E('RG.B', 2), E('RG.B', 3)], 4, 3);
    assert.deepStrictEqual(A.rows, B.rows);
    assert.deepStrictEqual(A.summary, B.summary);
    assert.strictEqual(A.summary.fullyConserved, 3);
    assert.strictEqual(A.summary.partialDisagreement, 0);
    ok('conservation: cluster of identical genomes (split or not) is fully conserved');
    // distinct genomes, same presence: 'RG.' and 'GG.' → x0 all, outs [1,1,0] tie, agreement .5; x1 unanimous
    // presenceAgreement (1+1)/2 = 1; outputAgreement (1+2)/(2+2) = .75
    const D = GA.conservation([E('RG.', 1), E('GG.', 1)], 3, 3);
    assert.deepStrictEqual([D.rows[0].level, D.rows[0].tie, D.rows[0].agreement, D.rows[0].unanimous], ['all', true, 0.5, false]);
    assert.strictEqual(D.summary.presenceAgreement, 1);
    close(D.summary.outputAgreement, 0.75, 'outputAgreement');
    assert.deepStrictEqual([D.summary.atLeast[100], D.summary.fullyConserved, D.summary.partialDisagreement], [2, 1, 1]);
    ok('conservation: distinct genomes with identical presence: presence agrees fully, output does not');
    // N = 0
    const Z = GA.conservation([], 4, 3);
    assert.strictEqual(Z.N, 0);
    assert.ok(Z.rows.every(r => r.level === 'absent' && r.prevalence === 0 && r.agreement === null));
    assert.strictEqual(Z.summary.presenceAgreement, null);
    assert.strictEqual(Z.summary.outputAgreement, null);
    noNaN(Z, 'conservation([])');
    ok('conservation: empty entries (N=0) give absent rows and no NaN');
}

// ================================================================ 3. compactness vs brute force over organisms
function bruteCompact(entries, repG) {
    const orgs = [];
    for (const e of entries) for (let k = 0; k < e.count; k++) orgs.push(e.g);
    const N = orgs.length, out = {N};
    if (repG) {
        const d = orgs.map(g => GA.genomeDistance(repG, g));
        out.toRep = {meanPresence: d.reduce((s, x) => s + x.presence, 0) / N, maxPresence: Math.max(...d.map(x => x.presence)),
            meanFine: d.reduce((s, x) => s + x.fine, 0) / N, maxFine: Math.max(...d.map(x => x.fine))};
    }
    const st = () => ({sum: 0, min: Infinity, max: -Infinity}), P = st(), O = st(), F = st(), hist = [];
    let pairs = 0, differ = 0;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
        const d = GA.genomeDistance(orgs[i], orgs[j]); pairs++;
        if (orgs[i] !== orgs[j]) differ++;
        for (const [a, v] of [[P, d.presence], [O, d.output], [F, d.fine]]) { a.sum += v; a.min = Math.min(a.min, v); a.max = Math.max(a.max, v); }
        hist[d.presence] = (hist[d.presence] || 0) + 1;
    }
    out.pairs = pairs;
    out.simpson = pairs ? differ / pairs : 0;
    const fin = a => ({mean: a.sum / pairs, min: a.min, max: a.max});
    out.pairwise = pairs ? {presence: fin(P), output: fin(O), fine: fin(F), weight: pairs, presenceHistogram: Array.from(hist, v => v || 0)} : null;
    return out;
}
function checkCompact(entries, repG, label) {
    const got = GA.compactness(entries, 3, repG), want = bruteCompact(entries, repG);
    assert.strictEqual(got.organisms, want.N, label);
    assert.strictEqual(got.pairs, want.pairs, label);
    assert.strictEqual(got.sampled, false, label);
    close(got.simpson, want.simpson, label + ' simpson');
    if (repG) for (const k of Object.keys(want.toRep)) close(got.toRep[k], want.toRep[k], `${label} toRep.${k}`);
    if (!want.pairwise) { assert.strictEqual(got.pairwise, null, label); return got; }
    assert.strictEqual(got.pairwise.weight, want.pairwise.weight, label);
    assert.deepStrictEqual(got.pairwise.presenceHistogram, want.pairwise.presenceHistogram, label + ' hist');
    for (const m of ['presence', 'output', 'fine']) {
        close(got.pairwise[m].mean, want.pairwise[m].mean, `${label} ${m}.mean`);
        assert.strictEqual(got.pairwise[m].min, want.pairwise[m].min, `${label} ${m}.min`);
        assert.strictEqual(got.pairwise[m].max, want.pairwise[m].max, `${label} ${m}.max`);
    }
    return got;
}
{
    // Hand case: A = 'RRR..' ×2, B = 'R...GB' ... keep it small: A = 'RRR' ×2, B = 'G..' ×1 (ncond 3)
    //   presence A–B = 2 (x1, x2), output A–B = 1 (x0 R vs G), fine 3. Pairs: A–A ×1 (0), A–B ×2.
    //   presence mean 4/3, min 0 (only because A has count 2), max 2; hist [1, 0, 2]; simpson 1 - 2/6 = 2/3
    //   rep = A: toRep mean presence = (0+0+2)/3 = 2/3, max 2; fine mean 3/3 = 1, max 3
    const got = checkCompact([E('RRR', 2), E('G..', 1)], 'RRR', 'hand');
    assert.deepStrictEqual(got.pairwise.presenceHistogram, [1, 0, 2]);
    assert.strictEqual(got.pairwise.presence.min, 0);
    close(got.pairwise.presence.mean, 4 / 3, 'mean');
    close(got.simpson, 2 / 3, 'simpson');
    close(got.toRep.meanPresence, 2 / 3, 'toRep'); assert.strictEqual(got.toRep.maxFine, 3);
    // without the duplicate the min is 2, not 0
    assert.strictEqual(checkCompact([E('RRR', 1), E('G..', 1)], 'RRR', 'no dup').pairwise.presence.min, 2);
    ok('compactness: hand case (min 0 only via an entry of count 2), simpson 2/3, toRep 2/3');

    // N = 1: no pairs, simpson 0, distance to itself 0
    const one = checkCompact([E('RG.', 1)], 'RG.', 'N=1');
    assert.strictEqual(one.pairwise, null); assert.strictEqual(one.simpson, 0); assert.strictEqual(one.toRep.maxFine, 0);
    // all identical, ×4: 6 pairs at 0, simpson 0
    const same = checkCompact([E('RG.', 4)], 'RG.', 'identical');
    assert.deepStrictEqual(same.pairwise.presenceHistogram, [6]); assert.strictEqual(same.simpson, 0);
    ok('compactness: N=1 and all-identical clusters');

    // random clusters: a few families, counts 1..4
    for (let trial = 0; trial < 60; trial++) {
        const base = randomGenome(), n = 1 + Math.floor(rand() * 10), entries = [], seen = new Set();
        let g = base;
        while (entries.length < n) {
            g = mutate(rand() < 0.5 ? base : g, 1 + Math.floor(rand() * 6));
            if (seen.has(g)) continue; seen.add(g);
            entries.push(E(g, 1 + Math.floor(rand() * 4), null, Math.floor(rand() * 3)));
        }
        checkCompact(entries, entries[0].g, `trial ${trial}`);
    }
    ok('compactness: 60 random clusters equal brute-force organism expansion (all pairwise stats, histogram, simpson, toRep)');

    // sampled path: 40 entries, maxEntries 7
    const entries = Array.from({length: 40}, (_, i) => E(mutate(randomGenome(), 0), 1 + (i % 3)));
    const s = GA.compactness(entries, 3, entries[0].g, {maxEntries: 7});
    assert.strictEqual(s.sampled, true);
    assert.ok(s.pairwise && s.pairwise.weight > 0);
    noNaN(s, 'sampled');
    // N = 0
    const z = GA.compactness([], 3, null);
    assert.strictEqual(z.organisms, 0); noNaN(z, 'compactness([])');
    ok('compactness: sampled path flags sampled and stays finite; N=0 does not throw');
}

// ================================================================ 4. lineage on a hand-built tree
//   founder 0 (born 0) ─ 2 (10) ─ 4 (20) ─ 7 (30)
//                      │        └ 5 (22)
//                      └ 3 (12) ─ 6 (25)
//   founder 1 (born 0) ─ 8 (15) ─ 9 (28)
//   deaths: 0@20 1@30 2@30 3@40 4@50 8@28, others alive (-1). Export tick 100.
const L = {
    parent: Int32Array.from([-1, -1, 0, 0, 2, 2, 3, 4, 1, 8]),
    founder: Int32Array.from([0, 1, 0, 0, 0, 0, 0, 0, 1, 1]),
    birth: Int32Array.from([0, 0, 10, 12, 20, 22, 25, 30, 15, 28]),
    death: Int32Array.from([20, 30, 30, 40, 50, -1, -1, -1, 28, -1]),
    tick: 100,
};
{
    assert.deepStrictEqual(GA.chainOf(L, 7), [0, 2, 4, 7]);
    assert.deepStrictEqual(GA.chainOf(L, 9), [1, 8, 9]);
    assert.deepStrictEqual(GA.chainOf(L, 0), [0]);
    assert.strictEqual(GA.commonPrefix([0, 2, 4, 7], [0, 2, 5]), 1);
    assert.strictEqual(GA.commonPrefix([0, 2, 4, 7], [0, 3, 6]), 0);
    assert.strictEqual(GA.commonPrefix([0, 2, 4, 7], [1, 8, 9]), -1);
    assert.strictEqual(GA.commonPrefix([0, 2], [0, 2, 4, 7]), 1);       // ancestor
    assert.strictEqual(GA.commonPrefix([0, 2, 4, 7], [0, 2, 4, 7]), 3); // self
    assert.strictEqual(GA.commonPrefix([], [0]), -1);
    const rel = (a, b) => GA.pairRelation(L, GA.chainOf(L, a), GA.chainOf(L, b));
    // 7,5: MRCA 2 (born 10): 7→4→2 = 2 gens, 5→2 = 1; treeDist 3, age 90
    assert.deepStrictEqual(rel(7, 5), {sameFounder: true, mrca: 2, genA: 2, genB: 1, treeDist: 3, ageTicks: 90});
    assert.deepStrictEqual(rel(7, 6), {sameFounder: true, mrca: 0, genA: 3, genB: 2, treeDist: 5, ageTicks: 100});
    assert.deepStrictEqual(rel(4, 7), {sameFounder: true, mrca: 4, genA: 0, genB: 1, treeDist: 1, ageTicks: 80});
    assert.deepStrictEqual(rel(7, 9), {sameFounder: false, mrca: null, genA: null, genB: null, treeDist: null, ageTicks: null});
    ok('lineage: chainOf, commonPrefix, pairRelation on a hand-built tree');
}
{
    // Cluster: E0 'RR..' ids [7,5]; E1 'RRR.' [6]; E2 'R...' [4]; E3 '..GG' [9]  → member order 7,5,6,4,9
    // presence: E0–E0 0, E0–E1 1, E0–E2 1, E0–E3 4, E1–E2 2, E1–E3 3, E2–E3 3
    // pairs (pres, tree):  7-5 (0, 3: mrca 2, 2/1)   7-6 (1, 5)   7-4 (1, 1 parent-child)   7-9 (4, -)
    //                      5-6 (1, 4: 2/2)   5-4 (1, 2 siblings)   5-9 (4, -)   6-4 (2, 4: 2/2)   6-9 (3, -)   4-9 (3, -)
    const entries = [E('RR..', 2, [7, 5]), E('RRR.', 1, [6]), E('R...', 1, [4]), E('..GG', 1, [9])];
    const R = GA.lineageRelations(L, entries, 3, {closeDistance: 1, farGenerations: 4, nearGenerations: 3, farDistance: 1});
    assert.strictEqual(R.available, true);
    assert.deepStrictEqual([R.organisms, R.considered, R.sampled, R.missing], [5, 5, false, 0]);
    assert.deepStrictEqual([R.pairs, R.sameFounder, R.parentChild, R.siblings, R.identicalGenomePairs], [10, 6, 1, 1, 1]);
    close(R.sameFounderShare, 0.6, 'sameFounderShare');
    // max(genA, genB): 7-5 2, 7-6 3, 7-4 1, 5-6 2, 5-4 1, 6-4 2 → ≤1: 2, ≤2: 5, ≤4 and up: 6
    assert.deepStrictEqual(R.sharedAncestorWithin.map(s => s.pairs), [2, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6]);
    assert.strictEqual(R.clusterMRCA, null);   // two founders
    assert.deepStrictEqual(R.founders, [{founder: 0, count: 4, share: 0.8}, {founder: 1, count: 1, share: 0.2}]);
    // closeUnrelated (presence ≤ 1, different founder or treeDist ≥ 4): 7-6 (td 5), 5-6 (td 4); sorted td desc
    assert.deepStrictEqual(R.closeUnrelated.map(p => [p.a, p.b, p.treeDist]), [[7, 6, 5], [5, 6, 4]]);
    assert.strictEqual(R.closeUnrelatedTotal, 2);
    // farRelated (same founder, treeDist ≤ 3, presence ≥ 1): 7-4 (td 1), 5-4 (td 2)
    assert.deepStrictEqual(R.farRelated.map(p => [p.a, p.b, p.treeDist]), [[7, 4, 1], [5, 4, 2]]);
    // nearest genome vs nearest relative, per organism:
    //   7: nearest genome 5 (0), its tree dist 3; nearest relative 4 (1) → miss
    //   5: nearest genome 7 (td 3); nearest relative 4 (td 2) → miss
    //   6: nearest genomes 7, 5 (1), best td 4; nearest relative td 4 → hit
    //   4: nearest genomes 7, 5 (1), best td 1; nearest relative 7 (1) → hit
    //   9: no relative → miss                                         → 2/5
    close(R.nearestGenomeIsNearestRelative, 0.4, 'nearestGenomeIsNearestRelative');
    // by presence distance: 0: 1 pair, age 90, tree 3; 1: 4 pairs, ages 100 80 100 90 → median 95, tree 5 1 4 2 → 3;
    //   2: 1 pair (age 100, tree 4); 3: 2 pairs, both different founders; 4: 2 pairs, different founders
    assert.deepStrictEqual(R.byPresenceDistance, [
        {presence: 0, pairs: 1, differentFounder: 0, medianAgeTicks: 90, medianTreeDist: 3},
        {presence: 1, pairs: 4, differentFounder: 0, medianAgeTicks: 95, medianTreeDist: 3},
        {presence: 2, pairs: 1, differentFounder: 0, medianAgeTicks: 100, medianTreeDist: 4},
        {presence: 3, pairs: 2, differentFounder: 2, medianAgeTicks: null, medianTreeDist: null},
        {presence: 4, pairs: 2, differentFounder: 2, medianAgeTicks: null, medianTreeDist: null}]);
    // Spearman presence vs tree over same-founder pairs: x [0,1,1,1,1,2], y [3,5,1,4,2,4]
    //   ranks x [0,2.5,2.5,2.5,2.5,5], y [2,5,0,3.5,1,3.5]; sxy 3.75, sxx 12.5, syy 17 → 3.75/sqrt(212.5)
    close(R.spearmanPresenceVsTree, 3.75 / Math.sqrt(212.5), 'spearman');
    assert.deepStrictEqual(R.organismsUsed, [7, 5, 6, 4, 9]);
    ok('lineageRelations: pair counts, founders, sharedAncestorWithin, close/far lists, nearest-relative 0.4, spearman');

    // single-founder cluster 7,5,6,4: MRCA 0 (age 100); deepest member 7 is 3 generations below it
    const R1 = GA.lineageRelations(L, entries.slice(0, 3), 3);
    assert.deepStrictEqual(R1.clusterMRCA, {id: 0, ageTicks: 100, maxGenerations: 3});
    // 7, 5: MRCA 2 (born 10, age 90), 7 is 2 generations below
    assert.deepStrictEqual(GA.lineageRelations(L, [E('RR..', 2, [7, 5])], 3).clusterMRCA, {id: 2, ageTicks: 90, maxGenerations: 2});
    // one member: its own MRCA, 0 generations
    assert.deepStrictEqual(GA.lineageRelations(L, [E('RR..', 1, [7])], 3).clusterMRCA, {id: 7, ageTicks: 70, maxGenerations: 0});
    // ancestor and descendant in one cluster: MRCA is the ancestor itself
    assert.deepStrictEqual(GA.lineageRelations(L, [E('R...', 1, [2]), E('RR..', 1, [7])], 3).clusterMRCA, {id: 2, ageTicks: 90, maxGenerations: 2});
    ok('lineageRelations: clusterMRCA (id, age, max generations), incl. single member and ancestor-in-cluster');

    // missing data
    assert.strictEqual(GA.lineageRelations(null, entries, 3).available, false);
    const noIds = GA.lineageRelations(L, [E('RR..', 2, [7, 5]), E('RRR.', 3, null)], 3);
    assert.strictEqual(noIds.available, false); assert.ok(/^3 organisms/.test(noIds.reason));
    const part = GA.lineageRelations(L, [E('RR..', 3, [7, 99, -1]), E('R...', 1, [4])], 3);
    assert.deepStrictEqual([part.available, part.organisms, part.missing, part.pairs], [true, 2, 2, 1]);
    assert.strictEqual(GA.lineageRelations(L, [E('RR..', 2, [99, 100])], 3).available, false);
    ok('lineageRelations: L null, entries without ids, ids out of range');

    // sampling: maxOrganisms 3 of 5
    const Rs = GA.lineageRelations(L, entries, 3, {maxOrganisms: 3});
    assert.deepStrictEqual([Rs.sampled, Rs.considered, Rs.organisms, Rs.pairs], [true, 3, 5, 3]);
    ok('lineageRelations: maxOrganisms samples and flags sampled');
}
{
    // DFS order: 4 [0,2,4], 7 [0,2,4,7], 5 [0,2,5], 3 [0,3], 6 [0,3,6], 9 [1,8,9]
    assert.deepStrictEqual(GA.lineageOrder(L, [9, 6, 5, 7, 4, 3]), [4, 7, 5, 3, 6, 9]);
    const byId = nodes => Object.fromEntries(nodes.map(n => [n.id, [n.parent, n.gens, n.member, n.founder]]));
    // members 7, 5, 6: keep founder 0, branch point 2; 4 and 3 are straight runs and disappear
    //   0 (root, 0)  2 (← 0, 1)  7 (← 2, 2 via 4)  5 (← 2, 1)  6 (← 0, 2 via 3)
    assert.deepStrictEqual(byId(GA.compressedTree(L, [7, 5, 6])), {
        0: [-1, 0, false, 0], 2: [0, 1, false, 0], 7: [2, 2, true, 0], 5: [2, 1, true, 0], 6: [0, 2, true, 0]});
    // member that is an ancestor of another member: 4 kept (← 0, 2 via 2), 7 (← 4, 1)
    assert.deepStrictEqual(byId(GA.compressedTree(L, [4, 7])), {0: [-1, 0, false, 0], 4: [0, 2, true, 0], 7: [4, 1, true, 0]});
    // two founders: 7 (← 0, 3), 9 (← 1, 2)
    assert.deepStrictEqual(byId(GA.compressedTree(L, [7, 9])), {0: [-1, 0, false, 0], 7: [0, 3, true, 0], 1: [-1, 0, false, 1], 9: [1, 2, true, 1]});
    ok('lineageOrder DFS order; compressedTree keeps members, branch points, founders with correct gens');

    // aliveAt(T) = population at the START of tick T (export convention): born before T, death ≥ T or none.
    //   T=28: 9 born at 28 → not yet; 8 died at 28 → still alive.  → 1 2 3 4 5 6 8
    //   T=29: 9 born at 28 → alive; 8 died at 28 → gone.  → 1 2 3 4 5 6 9
    assert.deepStrictEqual(GA.aliveAt(L, 28), [1, 2, 3, 4, 5, 6, 8]);
    assert.deepStrictEqual(GA.aliveAt(L, 29), [1, 2, 3, 4, 5, 6, 9]);
    //   T=31: 1, 2 died at 30 → gone; 7 born at 30 → alive  → 3 4 5 6 7 9
    assert.deepStrictEqual(GA.aliveAt(L, 31), [3, 4, 5, 6, 7, 9]);
    assert.deepStrictEqual(GA.aliveAt(L, 1), [0, 1]);
    assert.deepStrictEqual(GA.aliveAt(L, 0), []);   // nobody is born before tick 0
    assert.deepStrictEqual(GA.aliveAt(L, 101), [5, 6, 7, 9]);
    ok('aliveAt boundary semantics (start of tick T: born at T not yet alive, died at T still alive)');
}

// ================================================================ 5. spearman, median
{
    assert.strictEqual(GA.spearman([1, 2, 3], [10, 20, 30]), 1);
    assert.strictEqual(GA.spearman([1, 2, 3, 4], [9, 7, 5, 1]), -1);
    // ties: x [1,2,2,3] → ranks [0,1.5,1.5,3]; y [1,3,2,4] → [0,2,1,3]; sxy 4.5, sxx 4.5, syy 5 → sqrt(0.9)
    close(GA.spearman([1, 2, 2, 3], [1, 3, 2, 4]), Math.sqrt(0.9), 'tied spearman');
    assert.strictEqual(GA.spearman([1, 2], [2, 1]), null);
    assert.strictEqual(GA.spearman([1, 1, 1], [1, 2, 3]), null);   // no variance
    const a = [4, 1, 3, 2];
    assert.strictEqual(GA.median([3, 1, 2]), 2);
    assert.strictEqual(GA.median(a), 2.5);
    assert.deepStrictEqual(a, [4, 1, 3, 2]);   // input untouched
    assert.strictEqual(GA.median([10, 9]), 9.5);   // numeric, not string, sort
    assert.strictEqual(GA.median([]), null);
    ok('spearman (±1, ties = sqrt(0.9), n<3 null, constant null); median odd/even/empty');
}

// ================================================================ 6. compareConservation
{
    // A: 'RRR.' ×2, 'R.G.' ×2 (N=4). prevalence [1, .5, 1, 0]; tops [R, R, R(tie R/G), -]
    // B: 'GR.B' ×1, 'GRBB' ×1 (N=2). prevalence [1, 1, .5, 1]; tops [G, R, B, B]
    const A = GA.conservation([E('RRR.', 2), E('R.G.', 2)], 4, 3), B = GA.conservation([E('GR.B', 1), E('GRBB', 1)], 4, 3);
    const c = GA.compareConservation(A, B);
    // x0: dPrev 0, tv 1 (R vs G), outScore 1·min(1,1) = 1 → score 1, output
    // x1: dPrev -.5, tv 0 → score .5, presence
    // x2: dPrev .5, tv = (|.5-0|+|.5-0|+|0-1|)/2 = 1, outScore 1·.5 = .5 → score .5, presence (tie goes to presence)
    // x3: dPrev -1, tv null → score 1, presence
    // x2: A's majority output is a tie (R/G), so it is neither "differs" nor "same" (outputTie true)
    assert.deepStrictEqual(c.rows.map(r => [r.dPrev, r.tv, r.score, r.kind, r.outputDiffers, r.outputTie]), [
        [0, 1, 1, 'output', true, false], [-0.5, 0, 0.5, 'presence', false, false], [0.5, 1, 0.5, 'presence', false, true], [-1, null, 1, 'presence', false, false]]);
    // level 1 and .75: A {0,2}, B {0,1,3} → shared {0} (tops differ), onlyA {2}, onlyB {1,3}, jaccard 1/4
    // level .5 and .25: A {0,1,2}, B {0,1,2,3} → shared 3 (same top only x1), onlyB {3}, jaccard 3/4
    assert.deepStrictEqual(c.sets[1], {shared: 1, sharedSameOutput: 0, onlyA: 1, onlyB: 2, jaccard: 0.25});
    assert.deepStrictEqual(c.sets[0.75], {shared: 1, sharedSameOutput: 0, onlyA: 1, onlyB: 2, jaccard: 0.25});
    assert.deepStrictEqual(c.sets[0.5], {shared: 3, sharedSameOutput: 1, onlyA: 0, onlyB: 1, jaccard: 0.75});
    assert.deepStrictEqual(c.sets[0.25], {shared: 3, sharedSameOutput: 1, onlyA: 0, onlyB: 1, jaccard: 0.75});
    // prevalenceSimilarity = 1 - (0 + .5 + .5 + 1)/4 = .5
    assert.strictEqual(c.prevalenceSimilarity, 0.5);
    // distinguishing: score desc, stable: x0 (1), x3 (1), x1 (.5), x2 (.5)
    assert.deepStrictEqual(c.distinguishing.map(r => r.x), [0, 3, 1, 2]);
    // a cluster compared with itself: nothing distinguishes, similarity 1
    const self = GA.compareConservation(A, A);
    assert.strictEqual(self.distinguishing.length, 0); assert.strictEqual(self.prevalenceSimilarity, 1);
    ok('compareConservation: rows, tv, score/kind, sets per level, prevalenceSimilarity 0.5, distinguishing order');
}

// ================================================================ 7. populationProfile, compatible
{
    // organisms: 'RR.' ×2, 'GR.' ×1, 'R..' ×1 → N 4; x0 outs [3,1,0], x1 [3,0,0], x2 none
    // meanRules (4+3)/4 = 1.75; outputShare R 6/7, G 1/7; prevalence 1 → decile 9, .75 → decile 7
    const pats = [{exact: [['RR.', 2], ['GR.', 1]]}, {exact: [['R..', 1]]}];
    const p = GA.populationProfile(pats, 3, 3);
    assert.deepStrictEqual([p.survivors, p.uniqueGenomes, p.uniquePatterns, p.meanRules], [4, 3, 2, 1.75]);
    close(p.outputShare[0], 6 / 7, 'R share'); close(p.outputShare[1], 1 / 7, 'G share'); assert.strictEqual(p.outputShare[2], 0);
    assert.deepStrictEqual(p.prevalenceHistogram, [0, 0, 0, 0, 0, 0, 0, 1, 0, 1]);
    const ep = GA.populationProfile([], 3, 3);
    assert.strictEqual(ep.meanRules, null); noNaN(ep, 'empty profile');
    ok('populationProfile: survivors, meanRules 1.75, output shares 6/7 1/7, prevalence deciles; empty population');

    const run = (colours, conditions) => ({colours, conditions});
    const cond = [[0, 1], [1, 2], [2, 0]];
    assert.deepStrictEqual(GA.compatible(run(3, cond), run(3, cond.map(c => c.slice()))), {ok: true});
    assert.strictEqual(GA.compatible(run(3, cond), run(4, cond)).ok, false);
    assert.ok(/colour/.test(GA.compatible(run(3, cond), run(4, cond)).reason));
    assert.strictEqual(GA.compatible(run(3, cond), run(3, cond.slice(0, 2))).ok, false);
    assert.deepStrictEqual(GA.compatible(run(3, cond), run(3, [[0, 1], [2, 1], [2, 0]])), {ok: false, reason: 'condition 1 differs'});
    ok('compatible: same, different colours, different condition count, different condition list');
}

// ================================================================ 8. robustness
{
    // A cluster larger than the engine's argument limit (~120k) must not break clusterMRCA. Star tree: founder 0
    // with 200,000 children; MRCA 0, 1 generation.
    const n = 200001, Lb = {parent: new Int32Array(n).fill(0), founder: new Int32Array(n), birth: new Int32Array(n).fill(5), death: new Int32Array(n).fill(-1), tick: 10};
    Lb.parent[0] = -1; Lb.birth[0] = 0;
    const ids = Array.from({length: n - 1}, (_, i) => i + 1);
    const R = GA.lineageRelations(Lb, [E('R', n - 1, ids)], 3, {maxOrganisms: 20});
    assert.deepStrictEqual(R.clusterMRCA, {id: 0, ageTicks: 10, maxGenerations: 1});
    ok('lineageRelations: 200,000-member cluster (beyond the spread-argument limit) gives MRCA 0, 1 generation');
}

console.log(`\n${passed} checks passed`);
