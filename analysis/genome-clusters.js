// Coarse genome clustering of a surviving population, at a threshold the researcher chooses.
//
// Measure: each genome becomes 252 bits, rule present / absent; output colours are ignored.
// The difference between two genomes is the number of conditions where one has a rule and the
// other does not (Hamming distance on presence).
//
// Clustering: connected components (single linkage). Two genomes are linked when they differ in
// at most t rules, and links chain: if A–B = 1 and B–C = 1, then A, B and C share a cluster at
// t = 1 even though A–C = 2. Each cluster reports its spread and an upper bound on its diameter
// so chaining stays visible.
//
// Every threshold at once: organisms are grouped into unique presence patterns (with counts),
// and the minimum spanning tree over the patterns is built once. The clusters at t are the
// MST with every edge longer than t removed, so the explorer's slider never recomputes distances.
//
// usage: node analysis/genome-clusters.js <run dir> [<run dir> ...] > genome-clusters.json
//        [--tick T[,T...]]   checkpoints to analyze (default: each run's final snapshot)
//        [--every N]         every snapshot whose tick is a multiple of N, plus the final one
//        [--ticks all]       every snapshot (large output: ~1.5 MB per checkpoint)
//        [--thresholds 0,1,2,3,4,6,8,12,16]   thresholds summarized in the JSON and the report
//        [--top 8]           clusters listed per threshold in the report
//        [--present 0.08]    share a colour needs to count as surviving in the phase label
//        [--slim]            leave out per-pattern genomes (no explorer; report only)

const fs = require('fs');
const path = require('path');
const zlib = require('zlib');

const SYM = 'RGBYCM';
const WORDS = 8;   // 8 × 32 bits ≥ 252 conditions

// ---------------------------------------------------------------- measurement

// Genome string ('.' = no rule) → presence bits.
function presenceBits(g) {
    const b = new Uint32Array(Math.ceil(g.length / 32));
    for (let c = 0; c < g.length; c++) if (g[c] !== '.') b[c >>> 5] |= 1 << (c & 31);
    return b;
}
function popcount(x) {
    x -= (x >>> 1) & 0x55555555;
    x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
    return (((x + (x >>> 4)) & 0x0F0F0F0F) * 0x01010101) >>> 24;
}
function bitDistance(a, b) {
    let d = 0;
    for (let w = 0; w < a.length; w++) d += popcount(a[w] ^ b[w]);
    return d;
}
// Reference definition, used by the tests: count conditions where presence differs.
function coarseDistance(g1, g2) {
    let d = 0;
    for (let c = 0; c < g1.length; c++) if ((g1[c] === '.') !== (g2[c] === '.')) d++;
    return d;
}
// Fine distance (deferred in the analysis, provided for spot checks): any difference in output.
function fineDistance(g1, g2) {
    let d = 0;
    for (let c = 0; c < g1.length; c++) if (g1[c] !== g2[c]) d++;
    return d;
}
// Full comparison of two genomes, as described in the meeting notes.
function compareGenomes(g1, g2) {
    let shared = 0, onlyA = 0, onlyB = 0, colourDiff = 0;
    for (let c = 0; c < g1.length; c++) {
        const a = g1[c] !== '.', b = g2[c] !== '.';
        if (a && b) { shared++; if (g1[c] !== g2[c]) colourDiff++; }
        else if (a) onlyA++;
        else if (b) onlyB++;
    }
    const coarse = onlyA + onlyB, n = g1.length;
    return {identical: g1 === g2, sharedRules: shared, onlyInA: onlyA, onlyInB: onlyB, differentOutputs: colourDiff,
        coarse, fine: coarse + colourDiff, coarseSimilarity: 1 - coarse / n, fineSimilarity: 1 - (coarse + colourDiff) / n};
}

// [[genome, count], ...] → unique presence patterns, most organisms first.
function groupPatterns(genomes) {
    const byKey = new Map();
    for (const [g, count] of genomes) {
        const key = g.replace(/[^.]/g, '1').replace(/\./g, '0');
        let p = byKey.get(key);
        if (!p) byKey.set(key, p = {key, bits: presenceBits(g), count: 0, rules: 0, exact: []});
        p.count += count;
        p.exact.push([g, count]);
    }
    const patterns = [...byKey.values()];
    for (const p of patterns) {
        p.rules = [...p.key].filter(ch => ch === '1').length;
        p.exact.sort((a, b) => b[1] - a[1]);
    }
    patterns.sort((a, b) => b.count - a.count || (a.key < b.key ? -1 : 1));
    return patterns;
}

// ---------------------------------------------------------------- clustering

// Prim's MST over all pattern pairs, O(U²). Returns edges [i, j, w] sorted by w.
function spanningTree(patterns) {
    const U = patterns.length, edges = [];
    if (U < 2) return edges;
    const inTree = new Uint8Array(U), best = new Int32Array(U).fill(1 << 30), from = new Int32Array(U).fill(-1);
    let cur = 0;
    inTree[0] = 1;
    for (let added = 1; added < U; added++) {
        const cb = patterns[cur].bits;
        let next = -1, nd = 1 << 30;
        for (let j = 0; j < U; j++) {
            if (inTree[j]) continue;
            const d = bitDistance(cb, patterns[j].bits);
            if (d < best[j]) { best[j] = d; from[j] = cur; }
            if (best[j] < nd) { nd = best[j]; next = j; }
        }
        inTree[next] = 1;
        edges.push([from[next], next, nd]);
        cur = next;
    }
    edges.sort((a, b) => a[2] - b[2] || a[0] - b[0] || a[1] - b[1]);
    return edges;
}

function unionFind(n) {
    const parent = new Int32Array(n).map((_, i) => i);
    const find = x => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
    return {parent, find, union(a, b) { a = find(a); b = find(b); if (a === b) return false; parent[b] = a; return true; }};
}

// Cluster label per pattern at threshold t.
function clusterAt(patterns, edges, t) {
    const uf = unionFind(patterns.length);
    for (const [i, j, w] of edges) { if (w > t) break; uf.union(i, j); }
    return patterns.map((_, i) => uf.find(i));
}

// Representative clustering (no chaining): patterns, most organisms first, join the first
// representative within t rules; a pattern with none becomes a new representative. Every
// member is then within t of its representative, so a cluster's diameter is at most 2t.
function representativeClusters(patterns, t) {
    const reps = [], labels = new Array(patterns.length);
    for (let i = 0; i < patterns.length; i++) {
        let home = -1;
        for (const r of reps) if (bitDistance(patterns[r].bits, patterns[i].bits) <= t) { home = r; break; }
        if (home < 0) { reps.push(i); home = i; }
        labels[i] = home;
    }
    return labels;
}

// Brute force connected components, for the tests.
function bruteClusters(patterns, t) {
    const uf = unionFind(patterns.length);
    for (let i = 0; i < patterns.length; i++) for (let j = i + 1; j < patterns.length; j++)
        if (bitDistance(patterns[i].bits, patterns[j].bits) <= t) uf.union(i, j);
    return patterns.map((_, i) => uf.find(i));
}

// Clusters vs threshold: count, largest share, clusters holding ≥1% of the population.
function mergeCurve(patterns, edges, total) {
    const U = patterns.length, uf = unionFind(U), size = patterns.map(p => p.count);
    const maxT = edges.length ? edges[edges.length - 1][2] : 0;
    let clusters = U, largest = Math.max(0, ...size), e = 0;
    const big = () => { let k = 0; for (let i = 0; i < U; i++) if (uf.parent[i] === i && size[i] >= 0.01 * total) k++; return k; };
    const curve = [];
    for (let t = 0; t <= maxT; t++) {
        while (e < edges.length && edges[e][2] <= t) {
            const a = uf.find(edges[e][0]), b = uf.find(edges[e][1]);
            if (a !== b) { uf.parent[b] = a; size[a] += size[b]; largest = Math.max(largest, size[a]); clusters--; }
            e++;
        }
        curve.push({t, clusters, largestShare: +(largest / total).toFixed(4), clustersOver1pct: big()});
    }
    return curve;
}

// Summaries of every cluster at one threshold, largest first.
function summarize(patterns, labels, total, ncond) {
    const groups = new Map();
    labels.forEach((l, i) => { if (!groups.has(l)) groups.set(l, []); groups.get(l).push(i); });
    const out = [...groups.values()].map(mem => {
        mem.sort((a, b) => patterns[b].count - patterns[a].count || a - b);
        const rep = patterns[mem[0]];
        let organisms = 0, exact = 0, spread = 0, far = 0;
        const freq = new Float64Array(ncond);
        for (const i of mem) {
            const p = patterns[i];
            organisms += p.count; exact += p.exact.length;
            const d = bitDistance(rep.bits, p.bits);
            spread += d * p.count; far = Math.max(far, d);
            for (let c = 0; c < ncond; c++) if (p.key[c] === '1') freq[c] += p.count;
        }
        let consensus = '';
        for (let c = 0; c < ncond; c++) { freq[c] /= organisms; consensus += freq[c] >= 0.5 ? '1' : '0'; }
        return {organisms, share: +(organisms / total).toFixed(4), patterns: mem.length, exactGenomes: exact,
            representative: rep.exact[0][0], representativeRules: rep.rules,
            consensusRules: [...consensus].filter(ch => ch === '1').length,
            meanDistToRep: +(spread / organisms).toFixed(2), maxDistToRep: far, diameterAtMost: 2 * far,
            memberPatterns: mem};
    });
    out.sort((a, b) => b.organisms - a.organisms);
    return out;
}

// ---------------------------------------------------------------- run directories

function readSnapshots(dir) {
    const gz = path.join(dir, 'genome-snapshots.jsonl.gz'), plain = path.join(dir, 'genome-snapshots.jsonl');
    const text = fs.existsSync(gz) ? zlib.gunzipSync(fs.readFileSync(gz)).toString('utf8') : fs.readFileSync(plain, 'utf8');
    return text.split('\n').filter(l => l.trim()).map(l => JSON.parse(l));
}

// Phase label in the style of endpointBox (colours.js): colours above `present`, most common first.
function phaseLabel(sample, present) {
    if (!sample) return '';
    if (sample.pop < 0.01) return 'extinct';
    const order = sample.colors.map((f, k) => [f, k]).filter(x => x[0] > present).sort((a, b) => b[0] - a[0]).map(x => x[1]);
    return order.map(k => SYM[k]).join('>') + (sample.pop >= 0.95 ? ' full' : '');
}

function analyzeCheckpoint(snap, ncond, thresholds, slim) {
    const patterns = groupPatterns(snap.genomes);
    const total = patterns.reduce((s, p) => s + p.count, 0);
    const edges = spanningTree(patterns);
    const curve = mergeCurve(patterns, edges, total);
    const byThreshold = {};
    const strip = cl => cl.slice(0, 25).map(({memberPatterns, ...c}) => c);
    for (const t of thresholds) {
        const cc = summarize(patterns, clusterAt(patterns, edges, t), total, ncond);
        const rep = summarize(patterns, representativeClusters(patterns, t), total, ncond);
        byThreshold[t] = {clusters: cc.length, top: strip(cc),
            representative: {clusters: rep.length, top: strip(rep)}};
    }
    const cp = {tick: snap.tick, survivors: total, uniqueGenomes: snap.genomes.length, uniquePatterns: patterns.length,
        nearestNeighbour: nearestHistogram(patterns, edges), mergeCurve: curve, byThreshold, mst: edges};
    if (!slim) cp.patterns = patterns.map(p => p.exact);
    return cp;
}

// How far each pattern's nearest other pattern is, from the MST (the MST holds every pattern's
// nearest neighbour edge, so the minimum incident MST edge is the nearest-neighbour distance).
function nearestHistogram(patterns, edges) {
    const nn = new Array(patterns.length).fill(Infinity);
    for (const [i, j, w] of edges) { nn[i] = Math.min(nn[i], w); nn[j] = Math.min(nn[j], w); }
    const hist = {};
    nn.forEach((d, i) => { if (d !== Infinity) hist[d] = (hist[d] || 0) + patterns[i].count; });
    return hist;
}

function main() {
    const args = {}, dirs = [];
    for (let i = 2; i < process.argv.length; i++) {
        const a = process.argv[i];
        if (a === '--slim') args.slim = true;
        else if (a.startsWith('--')) args[a.slice(2)] = process.argv[++i];
        else dirs.push(a);
    }
    if (!dirs.length) { console.error('usage: node analysis/genome-clusters.js <run dir> [...] > genome-clusters.json'); process.exit(1); }
    const thresholds = (args.thresholds || '0,1,2,3,4,6,8,12,16').split(',').map(Number);
    const TOP = parseInt(args.top || '8'), PRESENT = Number(args.present || 0.08);

    const runs = [];
    for (const dir of dirs) {
        const run = JSON.parse(fs.readFileSync(path.join(dir, 'run.json'), 'utf8'));
        if (run.params && run.params.historyRules) { console.error(`${dir}: historyRules conditions are not supported`); process.exit(1); }
        let samples = [];
        try { samples = JSON.parse(fs.readFileSync(path.join(dir, 'samples.json'), 'utf8')).samples; } catch (e) { /* optional */ }
        const snaps = readSnapshots(dir), last = snaps[snaps.length - 1];
        let chosen;
        if (args.ticks === 'all') chosen = snaps;
        else if (args.every) { const N = parseInt(args.every); chosen = snaps.filter(s => s.tick % N === 0 || s === last); }
        else if (args.tick) {
            const want = args.tick.split(',').map(Number);
            chosen = want.map(t => snaps.find(s => s.tick === t) || (console.error(`${dir}: no snapshot at tick ${t}`), process.exit(1)));
        } else chosen = [last];
        const ncond = run.conditions.length;
        const checkpoints = chosen.map(snap => {
            const t0 = Date.now();
            const cp = analyzeCheckpoint(snap, ncond, thresholds, args.slim);
            cp.phase = phaseLabel(samples.find(s => s.tick === snap.tick), PRESENT);
            cp.final = snap === last;
            console.error(`${path.basename(dir)} tick ${snap.tick}: ${cp.survivors} survivors, ${cp.uniqueGenomes} genomes, ${cp.uniquePatterns} presence patterns (${Date.now() - t0} ms)`);
            return cp;
        });
        runs.push({name: path.basename(path.resolve(dir)), seed: run.seed, stopReason: run.stopReason, ticks: run.ticks,
            colours: run.params.n, conditions: run.conditions, checkpoints});
    }

    const result = {
        format: 'genome-clusters/1',
        method: {distance: 'coarse: Hamming distance over rule presence (rule present / absent), output colours ignored',
            clustering: 'connected components (single linkage) on unique presence patterns; patterns within t rules are linked and links chain',
            thresholds, present: PRESENT},
        runs,
    };
    process.stdout.write(JSON.stringify(result));

    // Readable report, in the shape agreed with Chris.
    const pct = x => (x * 100).toFixed(1) + '%';
    for (const r of runs) for (const cp of r.checkpoints) {
        console.error(`\n=== ${r.name} (seed ${r.seed}) tick ${cp.tick}${cp.final ? ' [final]' : ''}  phase ${cp.phase || '?'}`);
        console.error(`Survivors: ${cp.survivors}   unique genomes: ${cp.uniqueGenomes}   unique rule-presence patterns: ${cp.uniquePatterns}`);
        const nn = Object.entries(cp.nearestNeighbour).sort((a, b) => a[0] - b[0]).slice(0, 8)
            .map(([d, c]) => `${d}:${pct(c / cp.survivors)}`).join(' ');
        console.error(`Organisms by distance to nearest other pattern: ${nn}`);
        console.error('Threshold  clusters  largest  clusters>=1%');
        for (const row of cp.mergeCurve.filter(x => x.t <= Math.max(...thresholds)))
            console.error(`  ${String(row.t).padStart(4)}  ${String(row.clusters).padStart(9)}  ${pct(row.largestShare).padStart(7)}  ${String(row.clustersOver1pct).padStart(12)}`);
        const line = c => `    ${String(c.organisms).padStart(6)} organisms ${pct(c.share).padStart(6)}  patterns ${String(c.patterns).padStart(5)}  rules ${c.representativeRules}  spread ${c.meanDistToRep}  max ${c.maxDistToRep}`;
        for (const t of thresholds) {
            const b = cp.byThreshold[t];
            console.error(`\n  t = ${t}: ${b.clusters} connected-component clusters`);
            for (const c of b.top.slice(0, TOP)) console.error(line(c));
            console.error(`  t = ${t}: ${b.representative.clusters} representative clusters (no chaining)`);
            for (const c of b.representative.top.slice(0, TOP)) console.error(line(c));
        }
    }
}

module.exports = {presenceBits, bitDistance, coarseDistance, fineDistance, compareGenomes, groupPatterns, spanningTree,
    clusterAt, representativeClusters, bruteClusters, mergeCurve, summarize, analyzeCheckpoint};
if (require.main === module) main();
