// Checks for genome-clusters.js: the bit distance matches the plain definition, and clusters
// cut from the spanning tree match brute-force connected components at every threshold.
// usage: node analysis/test-genome-clusters.js
const assert = require('assert');
const gc = require('./genome-clusters.js');

let rs = 12345;
const rand = () => { rs = (rs * 1103515245 + 12345) >>> 0; return rs / 4294967296; };
const OUT = '.RGB';
const randomGenome = () => Array.from({length: 252}, () => OUT[Math.floor(rand() * 4)]).join('');
function mutate(g, k) {
    const a = [...g];
    for (let i = 0; i < k; i++) { const c = Math.floor(rand() * 252); a[c] = a[c] === '.' ? 'RGB'[Math.floor(rand() * 3)] : '.'; }
    return a.join('');
}
// Same partition, ignoring label names.
function samePartition(x, y) {
    const m = new Map(), n = new Map();
    for (let i = 0; i < x.length; i++) {
        if (m.has(x[i]) && m.get(x[i]) !== y[i]) return false;
        if (n.has(y[i]) && n.get(y[i]) !== x[i]) return false;
        m.set(x[i], y[i]); n.set(y[i], x[i]);
    }
    return true;
}
let passed = 0;
const ok = (name) => { passed++; console.log('ok  ' + name); };

// 1. Distances
{
    const a = 'R.G.B' + '.'.repeat(247), b = 'G.G..' + 'R' + '.'.repeat(246);
    // presence: a has 0,2,4; b has 0,2,5 → differs at 4 and 5
    assert.strictEqual(gc.coarseDistance(a, b), 2);
    assert.strictEqual(gc.bitDistance(gc.presenceBits(a), gc.presenceBits(b)), 2);
    assert.strictEqual(gc.fineDistance(a, b), 3);   // + colour change at 0
    const cmp = gc.compareGenomes(a, b);
    assert.deepStrictEqual([cmp.sharedRules, cmp.onlyInA, cmp.onlyInB, cmp.differentOutputs, cmp.coarse, cmp.fine], [2, 1, 1, 1, 2, 3]);
    for (let i = 0; i < 500; i++) {
        const g = randomGenome(), h = mutate(g, Math.floor(rand() * 30));
        assert.strictEqual(gc.bitDistance(gc.presenceBits(g), gc.presenceBits(h)), gc.coarseDistance(g, h));
    }
    ok('bit-packed coarse distance equals the per-condition definition');
}

// 2. Colour changes don't separate genomes in the coarse measure
{
    const g = randomGenome(), recoloured = g.replace(/R/g, 'G');
    const pats = gc.groupPatterns([[g, 3], [recoloured, 2]]);
    assert.strictEqual(pats.length, 1);
    assert.strictEqual(pats[0].count, 5);
    assert.strictEqual(pats[0].exact.length, 2);
    ok('genomes differing only in output colours share one presence pattern');
}

// 3. MST cut == brute-force connected components, on families of mutants
for (let trial = 0; trial < 5; trial++) {
    const genomes = [];
    for (let f = 0; f < 4; f++) {
        const founder = randomGenome();
        let g = founder;
        for (let k = 0; k < 40; k++) { g = mutate(rand() < 0.5 ? founder : g, 1 + Math.floor(rand() * 3)); genomes.push([g, 1 + Math.floor(rand() * 5)]); }
    }
    const pats = gc.groupPatterns(genomes), edges = gc.spanningTree(pats);
    assert.strictEqual(edges.length, pats.length - 1);
    for (let t = 0; t <= 40; t++) assert.ok(samePartition(gc.clusterAt(pats, edges, t), gc.bruteClusters(pats, t)), `trial ${trial} t ${t}`);
    const total = pats.reduce((s, p) => s + p.count, 0);
    const curve = gc.mergeCurve(pats, edges, total);
    for (const row of curve) assert.strictEqual(row.clusters, new Set(gc.clusterAt(pats, edges, row.t)).size);
    assert.strictEqual(curve[curve.length - 1].clusters, 1);
}
ok('spanning-tree clusters equal brute-force connected components at t = 0..40');
ok('merge curve cluster counts match');

// 4. Chris's example: one founder + 1-rule mutants → one cluster at t = 1, separate at t = 0
{
    const founder = randomGenome(), genomes = [[founder, 50]];
    const used = new Set();
    for (let i = 0; i < 10; i++) {
        let c; do c = Math.floor(rand() * 252); while (used.has(c)); used.add(c);
        const a = [...founder]; a[c] = a[c] === '.' ? 'B' : '.'; genomes.push([a.join(''), 1]);
    }
    const pats = gc.groupPatterns(genomes), edges = gc.spanningTree(pats);
    assert.strictEqual(new Set(gc.clusterAt(pats, edges, 0)).size, 11);
    assert.strictEqual(new Set(gc.clusterAt(pats, edges, 1)).size, 1);
    const [cl] = gc.summarize(pats, gc.clusterAt(pats, edges, 1), 60, 252);
    assert.strictEqual(cl.organisms, 60);
    assert.strictEqual(cl.representative, founder);
    assert.strictEqual(cl.maxDistToRep, 1);
    ok('founder plus single-rule mutants: 11 clusters at t=0, 1 cluster at t=1, founder is representative');
}

// 5. Two unrelated families stay apart at mutant-scale thresholds
{
    const A = randomGenome(), B = randomGenome(), genomes = [];
    for (let i = 0; i < 30; i++) { genomes.push([mutate(A, 2), 1]); genomes.push([mutate(B, 2), 1]); }
    const pats = gc.groupPatterns(genomes), edges = gc.spanningTree(pats);
    const cl = gc.summarize(pats, gc.clusterAt(pats, edges, 6), 60, 252);
    assert.strictEqual(cl.length, 2);
    assert.deepStrictEqual(cl.map(c => c.organisms), [30, 30]);
    ok('two unrelated families give two clusters of 30 at t=6');
}
// 6. Representative clustering never chains: every member within t of its representative
{
    const A = randomGenome(), genomes = [];
    let g = A;
    for (let i = 0; i < 60; i++) { g = mutate(g, 1); genomes.push([g, 1 + (i % 3)]); }   // a chain drifting away
    const pats = gc.groupPatterns(genomes), edges = gc.spanningTree(pats);
    for (const t of [1, 2, 4, 8]) {
        const labels = gc.representativeClusters(pats, t);
        labels.forEach((r, i) => assert.ok(gc.bitDistance(pats[r].bits, pats[i].bits) <= t));
        assert.ok(new Set(labels).size >= new Set(gc.clusterAt(pats, edges, t)).size);
    }
    ok('representative clusters keep every member within t, never fewer clusters than connected components');
}
console.log(`\n${passed} checks passed`);
