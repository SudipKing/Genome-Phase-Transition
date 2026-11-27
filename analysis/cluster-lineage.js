// Do genome clusters follow the family tree? Links the coarse genome clusters of the organisms alive in a
// .hxlg export (as in genome-cluster-explorer.html) to their lineage, which the same file records.
//
// For the living population it reports:
//   - survivors, unique genomes and presence patterns, rules per genome, output colour shares
//   - surviving founders, mean generations, and whether the whole population has a common ancestor
//   - nearest-neighbour presence distances, and the most recent common ancestor (MRCA) of each organism and
//     its nearest genome, against random pairs (age in ticks, generations, share with different founders)
//   - median MRCA age as a function of presence distance: how far genome distance tracks ancestry
//   - output-colour differences between nearest neighbours that the coarse (presence-only) distance hides
//   - the connected-cluster merge curve (clusters and largest share at each t) and the t where the largest
//     cluster first reaches 1, 5, 10, 25, 50, 90 and 99% of the population
//   - the clusters at --t (connected) and the top representative clusters at several t, each with its
//     founders, the age and depth of its members' common ancestor, and pairwise MRCA ages
//
// Distances and clustering match genome-clusters.js and the explorer: 252 presence bits, Hamming distance,
// connected components from the minimum spanning tree, representative clusters over patterns ordered by
// organism count then pattern key. Cluster numbers therefore match the explorer's table.
//
// usage: node analysis/cluster-lineage.js <file.hxlg> > cluster-lineage.json
//        [--t 4]          threshold for the connected-cluster table
//        [--top 12]       clusters listed at --t
//        [--pairs 4000]   pairs sampled per presence distance for the MRCA-by-distance table
//        [--random 20000] random pairs sampled for the baseline
//        [--seed 1]       seed for the sampling, so reruns give the same numbers
//
// Run on hex-ca-lineage-genomes-tick-35400.hxlg (2,502 survivors) it takes about 10 s.

const fs = require('fs');

const args = process.argv.slice(2);
const opt = (name, def) => { const i = args.indexOf('--' + name); return i >= 0 ? +args[i + 1] : def; };
const file = args.find(a => !a.startsWith('--') && args[args.indexOf(a) - 1]?.startsWith('--') !== true);
if (!file) { console.error('usage: node analysis/cluster-lineage.js <file.hxlg> [--t 4] [--top 12] [--pairs 4000] [--random 20000] [--seed 1]'); process.exit(1); }
const T_TABLE = opt('t', 4), TOP = opt('top', 12), PAIRS = opt('pairs', 4000), RANDOM = opt('random', 20000);

// Seeded PRNG (mulberry32), so sampled figures are reproducible.
let seed = opt('seed', 1) >>> 0;
const rand = () => { seed = (seed + 0x6D2B79F5) >>> 0; let x = seed; x = Math.imul(x ^ (x >>> 15), x | 1); x ^= x + Math.imul(x ^ (x >>> 7), x | 61); return ((x ^ (x >>> 14)) >>> 0) / 4294967296; };

// ---------------------------------------------------------------- read the .hxlg (same format as organismlog.js)
const buf = fs.readFileSync(file);
const ab = buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength);
const u8 = new Uint8Array(ab);
if (u8[0] !== 72 || u8[1] !== 88 || u8[2] !== 76 || u8[3] !== 71) throw new Error('Not a lineage + genomes file (missing HXLG marker).');
const headerLen = new DataView(ab).getUint32(8, true);
const H = JSON.parse(new TextDecoder().decode(u8.subarray(12, 12 + headerLen)));
if (H.useHistory) throw new Error('This run used history rules; this analysis does not support those conditions.');
const start = (12 + headerLen + 3) & ~3, types = {Int32Array, Uint32Array, Uint16Array, Int8Array};
const D = {};
for (const s of H.sections) D[s.name] = new types[s.type](ab, start + s.offset, s.length);
const T = H.exportedAtTick, NC = H.ncond, SYM = '.' + 'RGBYCM'.slice(0, H.n), WORDS = Math.ceil(NC / 32);

// ---------------------------------------------------------------- living genomes and their ancestor chains
function genomeOf(id) {   // replay rule changes from the founder down to `id`
  const g = new Int8Array(NC).fill(-1), chain = [];
  for (let a = id; a >= 0; a = D.parent[a]) chain.push(a);
  for (let k = chain.length - 1; k >= 0; k--) {
    const a = chain[k];
    for (let d = D.deltaStart[a]; d < D.deltaStart[a + 1]; d++) g[D.deltaCond[d]] = D.deltaVal[d];
  }
  return {g, chain};
}
function popcount(x) { x -= (x >>> 1) & 0x55555555; x = (x & 0x33333333) + ((x >>> 2) & 0x33333333); return (((x + (x >>> 4)) & 0x0F0F0F0F) * 0x01010101) >>> 24; }

const raw = H.liveAtExport.filter(id => id >= 0 && id < H.organisms);
const rawData = raw.map(id => {
  const {g, chain} = genomeOf(id);
  let s = ''; const bits = new Uint32Array(WORDS);
  for (let c = 0; c < NC; c++) { s += SYM[g[c] + 1]; if (g[c] >= 0) bits[c >>> 5] |= 1 << (c & 31); }
  return {id, g, chain, s, bits, key: s.replace(/[^.]/g, '1').replace(/\./g, '0')};
});
// Order organisms as the explorer orders presence patterns (most organisms first, then pattern key), so
// representative clusters and cluster numbering come out the same.
const patCount = new Map();
for (const o of rawData) patCount.set(o.key, (patCount.get(o.key) || 0) + 1);
rawData.sort((a, b) => patCount.get(b.key) - patCount.get(a.key) || (a.key < b.key ? -1 : a.key > b.key ? 1 : 0));
const N = rawData.length;
const live = rawData.map(o => o.id), G = rawData.map(o => o.g), CH = rawData.map(o => o.chain), B = rawData.map(o => o.bits);

const pres = (i, j) => { let d = 0; for (let w = 0; w < WORDS; w++) d += popcount(B[i][w] ^ B[j][w]); return d; };
const fine = (i, j) => { let d = 0; for (let c = 0; c < NC; c++) if (G[i][c] !== G[j][c]) d++; return d; };
const outDiff = (i, j) => { let d = 0; for (let c = 0; c < NC; c++) if (G[i][c] >= 0 && G[j][c] >= 0 && G[i][c] !== G[j][c]) d++; return d; };
const summarize = arr => {
  if (!arr.length) return {n: 0};
  const a = arr.slice().sort((x, y) => x - y), q = p => a[Math.floor(p * (a.length - 1))];
  return {n: a.length, median: q(0.5), p10: q(0.1), p90: q(0.9), mean: Math.round(a.reduce((s, v) => s + v, 0) / a.length)};
};

// Pair MRCA: walk j's chain until it meets i's. No meeting point = different founders.
const depthMaps = CH.map(c => { const m = new Map(); c.forEach((a, k) => m.set(a, k)); return m; });
function mrca(i, j) {
  const m = depthMaps[i], c = CH[j];
  for (let k = 0; k < c.length; k++) if (m.has(c[k])) return {id: c[k], age: T - D.birth[c[k]], gens: m.get(c[k]) + k};
  return null;
}
// Cluster MRCA, as the explorer computes it: the oldest point where members' chains meet the first member's.
function clusterAncestor(mem) {
  const founders = new Set(mem.map(i => D.founder[live[i]]));
  if (mem.length < 2 || founders.size > 1) return {founders: founders.size, ancestorAge: null, generations: null};
  let deepest = 0; const hits = [];
  for (let k = 1; k < mem.length; k++) {
    const c = CH[mem[k]]; let s = 0;
    while (s < c.length && !depthMaps[mem[0]].has(c[s])) s++;
    const at = depthMaps[mem[0]].get(c[s]);
    hits.push(s - at); if (at > deepest) deepest = at;
  }
  const anc = CH[mem[0]][deepest];
  return {founders: 1, ancestorAge: T - D.birth[anc], generations: Math.max(deepest, ...hits.map(h => h + deepest))};
}
function clusterLineage(mem, maxPairs = 300) {
  const sub = mem.length > maxPairs ? mem.slice(0, maxPairs) : mem, ages = [];
  for (let a = 0; a < sub.length; a++) for (let b = a + 1; b < sub.length; b++) { const m = mrca(sub[a], sub[b]); if (m) ages.push(m.age); }
  return {size: mem.length, rulesInRep: G[mem[0]].reduce((s, v) => s + (v >= 0), 0), ...clusterAncestor(mem),
    pairMRCAageTicks: summarize(ages), pairsSampledFrom: sub.length};
}

const out = {file: file.split(/[\\/]/).pop(), tick: T, seed: opt('seed', 1)};

// ---------------------------------------------------------------- population
out.survivors = N;
out.uniqueGenomes = new Set(rawData.map(o => o.s)).size;
out.uniquePatterns = patCount.size;
const rules = G.map(g => g.reduce((s, v) => s + (v >= 0), 0));
out.rulesPerGenome = {mean: +(rules.reduce((s, v) => s + v, 0) / N).toFixed(1), min: Math.min(...rules), max: Math.max(...rules)};
const outCol = new Array(H.n).fill(0); for (const g of G) for (const v of g) if (v >= 0) outCol[v]++;
const outTot = outCol.reduce((s, v) => s + v, 0);
out.outputColourShare = Object.fromEntries(outCol.map((v, k) => [SYM[k + 1], +(v / outTot).toFixed(3)]));
if (H.series && H.series.colors) { const k = H.series.population.length - 1; out.livingCellColourShare = Object.fromEntries(H.series.colors.map((c, i) => [SYM[i + 1], +c[k].toFixed(3)])); }

// ---------------------------------------------------------------- lineage basics
const founders = new Map(); for (const id of live) founders.set(D.founder[id], (founders.get(D.founder[id]) || 0) + 1);
out.survivingFounders = founders.size;
out.largestFounderShare = +(Math.max(...founders.values()) / N).toFixed(3);
out.generationsFromFounderMean = Math.round(CH.reduce((s, c) => s + c.length - 1, 0) / N);
{
  let common = new Set(CH[0]);
  for (let i = 1; i < N && common.size; i++) { const s = depthMaps[i]; common = new Set([...common].filter(a => s.has(a))); }
  const best = CH[0].find(a => common.has(a));
  out.populationMRCA = best === undefined ? null : {birthTick: D.birth[best], ageTicks: T - D.birth[best]};
}

// ---------------------------------------------------------------- nearest neighbours, distance bins, MST
const nn = new Int32Array(N).fill(-1), nnd = new Int32Array(N).fill(1 << 30), nnTies = new Int32Array(N);
const bins = new Map();   // presence distance -> reservoir sample of pairs
const binSeen = new Map();
for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
  const d = pres(i, j);
  for (const [x, y] of [[i, j], [j, i]]) {   // ties broken uniformly (reservoir over equal distances)
    if (d < nnd[x]) { nnd[x] = d; nn[x] = y; nnTies[x] = 1; }
    else if (d === nnd[x] && rand() * ++nnTies[x] < 1) nn[x] = y;
  }
  const seen = (binSeen.get(d) || 0) + 1; binSeen.set(d, seen);
  let l = bins.get(d); if (!l) bins.set(d, l = []);
  if (l.length < PAIRS) l.push([i, j]); else { const r = Math.floor(rand() * seen); if (r < PAIRS) l[r] = [i, j]; }
}
const hist = {}; for (const d of nnd) hist[d] = (hist[d] || 0) + 1;
out.nearestNeighbour = {
  histogram: hist,
  within4Pct: +(100 * [...nnd].filter(d => d <= 4).length / N).toFixed(1),
  identicalPattern: [...nnd].filter(d => d === 0).length,
};

const nnM = []; let nnNone = 0;
for (let i = 0; i < N; i++) { const m = mrca(i, nn[i]); if (m) nnM.push(m); else nnNone++; }
out.nearestNeighbour.mrcaAgeTicks = summarize(nnM.map(m => m.age));
out.nearestNeighbour.mrcaGenerations = summarize(nnM.map(m => m.gens));
out.nearestNeighbour.differentFounder = nnNone;
const od = []; for (let i = 0; i < N; i++) od.push(outDiff(i, nn[i]));
out.nearestNeighbour.hiddenOutputDiffs = summarize(od);
out.nearestNeighbour.shareWithHiddenOutputDiff = +(od.filter(x => x > 0).length / N).toFixed(2);

const rnd = []; let rndNone = 0, rndTried = 0;
while (rndTried < RANDOM) {
  const i = Math.floor(rand() * N), j = Math.floor(rand() * N); if (i === j) continue; rndTried++;
  const m = mrca(i, j); if (!m) { rndNone++; continue; }
  rnd.push({age: m.age, gens: m.gens, pres: pres(i, j), fine: fine(i, j)});
}
out.randomPairs = {
  sampled: RANDOM, differentFounder: rndNone,
  mrcaAgeTicks: summarize(rnd.map(r => r.age)), mrcaGenerations: summarize(rnd.map(r => r.gens)),
  presenceDistance: summarize(rnd.map(r => r.pres)), fineDistance: summarize(rnd.map(r => r.fine)),
};

out.mrcaAgeByPresenceDistance = {};
for (const d of [...bins.keys()].sort((a, b) => a - b)) {
  const l = bins.get(d), ms = l.map(([i, j]) => mrca(i, j)), ages = ms.filter(Boolean).map(m => m.age);
  out.mrcaAgeByPresenceDistance[d] = {pairs: binSeen.get(d), sampled: l.length, medianAgeTicks: summarize(ages).median ?? null,
    differentFounder: ms.filter(m => !m).length};
}

// Minimum spanning tree (Prim's, O(N²)); connected clusters at t = MST minus edges longer than t.
const mst = [];
{
  const inTree = new Uint8Array(N), best = new Int32Array(N).fill(1 << 30), from = new Int32Array(N).fill(-1);
  let cur = 0; inTree[0] = 1;
  for (let added = 1; added < N; added++) {
    let next = -1, nd = 1 << 30;
    for (let j = 0; j < N; j++) {
      if (inTree[j]) continue;
      const d = pres(cur, j); if (d < best[j]) { best[j] = d; from[j] = cur; }
      if (best[j] < nd) { nd = best[j]; next = j; }
    }
    inTree[next] = 1; mst.push([from[next], next, nd]); cur = next;
  }
  mst.sort((a, b) => a[2] - b[2]);
}
function ccAt(t) {
  const p = Int32Array.from({length: N}, (_, i) => i);
  const f = x => { while (p[x] !== x) { p[x] = p[p[x]]; x = p[x]; } return x; };
  for (const [i, j, w] of mst) { if (w > t) break; const a = f(i), b = f(j); if (a !== b) p[b] = a; }
  const g = new Map(); for (let i = 0; i < N; i++) { const r = f(i); if (!g.has(r)) g.set(r, []); g.get(r).push(i); }
  return [...g.values()].sort((a, b) => b.length - a.length);   // stable: ties keep explorer order
}
const maxT = mst.length ? mst[mst.length - 1][2] : 0;
out.mergeCurve = [];
for (let t = 0; t <= maxT; t++) { const c = ccAt(t); out.mergeCurve.push({t, clusters: c.length, largestPct: +(100 * c[0].length / N).toFixed(1)}); }
out.largestClusterReaches = Object.fromEntries([1, 5, 10, 25, 50, 90, 99].map(p => [p + '%', (out.mergeCurve.find(r => r.largestPct >= p) || {}).t ?? null]));

// ---------------------------------------------------------------- clusters with lineage
const cc = ccAt(T_TABLE), multi = cc.filter(c => c.length >= 2);
const within = []; for (const c of multi) for (let a = 0; a < c.length; a++) for (let b = a + 1; b < c.length; b++) { const m = mrca(c[a], c[b]); if (m) within.push(m.age); }
out.connected = {t: T_TABLE, clusters: cc.length, singletons: cc.length - multi.length, largest: cc[0].length,
  singleFounderClusters: multi.filter(c => clusterAncestor(c).founders === 1).length, multiOrganismClusters: multi.length,
  pairMRCAageTicksAllClusters: summarize(within), top: cc.slice(0, TOP).map(c => clusterLineage(c))};
// the giant cluster's ancestry as t crosses the horizon
out.largestClusterAncestry = [8, 16, 24, 28, 30, 32, 33, 34, 36, 40].filter(t => t <= maxT).map(t => {
  const c = ccAt(t)[0], a = clusterAncestor(c);
  return {t, size: c.length, founders: a.founders, ancestorAgeTicks: a.ancestorAge, generations: a.generations};
});

function repClusters(t) {   // explorer's representative mode, over organisms in pattern order
  const reps = [], lab = new Int32Array(N);
  for (let i = 0; i < N; i++) { let h = -1; for (const r of reps) if (pres(r, i) <= t) { h = r; break; } if (h < 0) { reps.push(i); h = i; } lab[i] = h; }
  const g = new Map(); for (let i = 0; i < N; i++) { if (!g.has(lab[i])) g.set(lab[i], []); g.get(lab[i]).push(i); }
  return [...g.values()].sort((a, b) => b.length - a.length);
}
out.representative = [4, 8, 12, 16, 20, 24].map(t => {
  const r = repClusters(t);
  return {t, clusters: r.length, largestPct: +(100 * r[0].length / N).toFixed(2), top: r.slice(0, 3).map(c => clusterLineage(c))};
});

console.log(JSON.stringify(out, null, 1));
const nnA = out.nearestNeighbour.mrcaAgeTicks, rA = out.randomPairs.mrcaAgeTicks;
console.error(`${out.file}: ${N} survivors, ${out.uniqueGenomes} genomes, ${out.survivingFounders} surviving founders.\n` +
  `Nearest genome shares an ancestor a median ${nnA.median} ticks back (random pair: ${rA.median}; ${rndNone} of ${RANDOM} random pairs share none).\n` +
  `At t = ${T_TABLE}: ${cc.length} clusters, ${out.connected.singleFounderClusters} of ${multi.length} multi-organism clusters descend from one founder.`);
