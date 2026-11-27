// Real-data verification of the cluster analysis core (GA) in genome-cluster-explorer.html. Not a unit test:
// it needs data files, and each part is skipped with a note when its file is absent.
//
//   1. genome-clusters.json: clusters at t = 4 and 16 rebuilt as the explorer does (MST cut); conservation
//      carriers/outputs equal an independent recount from genome strings; compactness toRep mean/max equal the
//      explorer's spread / far for the top 10 clusters; timings on the largest cluster.
//   2. a .hxlg export: aliveAt(exportedAtTick) == liveAtExport; cluster MRCAs at t = 4 against the explorer's
//      clusterLineage and a brute-force ancestor intersection; 2,000 random pair MRCAs against a naive walk;
//      lineageRelations timing on the largest cluster.
//   3. an earlier .hxlg of the same run: is it the same simulation? If so, the population at its tick rebuilt
//      from the later file (aliveAt + delta replay) is compared organism by organism.
//
// usage: node analysis/verify-cluster-analysis.js [--json analysis/runs/genome-clusters.json]
//          [--hxlg <late.hxlg>] [--early <early.hxlg>]
const fs = require('fs');
const path = require('path');
const gc = require('./genome-clusters.js');

const html = fs.readFileSync(path.join(__dirname, '..', 'genome-cluster-explorer.html'), 'utf8');
const GA = new Function('module', html.match(/<script id="cluster-analysis-core">([\s\S]*?)<\/script>/)[1] + '\nreturn GA;')({});

const args = process.argv.slice(2), opt = (k, d) => { const i = args.indexOf('--' + k); return i >= 0 ? args[i + 1] : d; };
const DL = path.join(process.env.USERPROFILE || process.env.HOME || '', 'Downloads');
const JSONF = opt('json', path.join(__dirname, 'runs', 'genome-clusters.json'));
const LATE = opt('hxlg', path.join(DL, 'hex-ca-lineage-genomes-tick-35400.hxlg'));
const EARLY = opt('early', path.join(DL, 'hex-ca-lineage-genomes-tick-5100.hxlg'));

let seed = 7;
const rand = () => { seed = (seed + 0x6D2B79F5) >>> 0; let x = seed; x = Math.imul(x ^ (x >>> 15), x | 1); x ^= x + Math.imul(x ^ (x >>> 7), x | 61); return ((x ^ (x >>> 14)) >>> 0) / 4294967296; };
let fails = 0, checks = 0;
const check = (cond, msg) => { checks++; if (!cond) { fails++; console.log('FAIL ' + msg); } };
const ms = f => { const t = process.hrtime.bigint(); const r = f(); return [r, Number(process.hrtime.bigint() - t) / 1e6]; };

// ---------------------------------------------------------------- explorer copies (main <script>)
function clustersAt(pats, mst, t) {   // labelsCC + clustersAt (connected), same as the explorer
    const lab = gc.clusterAt(pats, mst, t), groups = new Map();
    lab.forEach((l, i) => { if (!groups.has(l)) groups.set(l, []); groups.get(l).push(i); });
    return [...groups.values()].map(mem => {
        mem.sort((a, b) => pats[b].count - pats[a].count || a - b);
        const rep = pats[mem[0]];
        let organisms = 0, spread = 0, far = 0;
        for (const i of mem) { const d = gc.bitDistance(rep.bits, pats[i].bits); organisms += pats[i].count; spread += d * pats[i].count; if (d > far) far = d; }
        return {mem, rep: mem[0], organisms, spread: spread / organisms, far};
    }).sort((a, b) => b.organisms - a.organisms);
}
function modelFromPatterns(patterns) {   // buildModel
    return patterns.map(exact => ({exact, count: exact.reduce((s, e) => s + e[1], 0), bits: gc.presenceBits(exact[0][0])}));
}

// ================================================================ 1. genome-clusters.json
if (!fs.existsSync(JSONF)) console.log(`skip 1: ${JSONF} not found`);
else {
    const J = JSON.parse(fs.readFileSync(JSONF, 'utf8'));
    for (const run of J.runs.filter(r => ['c1', 'e1'].includes(r.name))) {
        const cp = run.checkpoints[run.checkpoints.length - 1], ncond = run.conditions.length, ncol = run.colours;
        const pats = modelFromPatterns(cp.patterns);
        for (const t of [4, 16]) {
            const cl = clustersAt(pats, cp.mst, t);
            let recountOK = 0, spreadOK = 0;
            for (const c of cl.slice(0, 10)) {
                const entries = GA.entriesFromPatterns(pats, c.mem);
                const cons = GA.conservation(entries, ncond, ncol);
                // independent recount from genome strings
                const car = new Array(ncond).fill(0), outs = Array.from({length: ncond}, () => new Array(ncol).fill(0));
                let N = 0;
                for (const i of c.mem) for (const [g, k] of cp.patterns[i]) {
                    N += k;
                    for (let x = 0; x < ncond; x++) if (g[x] !== '.') { car[x] += k; outs[x]['RGBYCM'.indexOf(g[x])] += k; }
                }
                let good = N === c.organisms && cons.N === N;
                for (let x = 0; x < ncond; x++) {
                    const r = cons.rows[x];
                    if (r.carriers !== car[x] || r.prevalence !== car[x] / N || r.outs.join() !== outs[x].join()) good = false;
                }
                check(good, `${run.name} t=${t} cluster of ${c.organisms}: conservation recount`);
                recountOK += good;
                const comp = GA.compactness(entries, ncol, pats[c.rep].exact[0][0]);
                const sOK = comp.toRep.meanPresence === c.spread && comp.toRep.maxPresence === c.far;
                check(sOK, `${run.name} t=${t}: toRep ${comp.toRep.meanPresence}/${comp.toRep.maxPresence} vs spread ${c.spread}/far ${c.far}`);
                spreadOK += sOK;
            }
            const big = cl[0], entries = GA.entriesFromPatterns(pats, big.mem);
            const [, tc] = ms(() => GA.conservation(entries, ncond, ncol));
            const [comp, tk] = ms(() => GA.compactness(entries, ncol, pats[big.rep].exact[0][0]));
            console.log(`1. ${run.name} tick ${cp.tick} t=${t}: ${cl.length} clusters; top 10: recount ${recountOK}/10, spread/far ${spreadOK}/10 exact; ` +
                `largest ${big.organisms} organisms, ${entries.length} genomes: conservation ${tc.toFixed(1)} ms, compactness ${tk.toFixed(0)} ms` +
                ` (sampled ${comp.sampled}, mean presence ${comp.pairwise ? comp.pairwise.presence.mean.toFixed(2) : '-'})`);
        }
    }
}

// ================================================================ .hxlg reader (same as parseHxlg)
function readHxlg(file) {
    const buf = fs.readFileSync(file), ab = buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength), u8 = new Uint8Array(ab);
    if (String.fromCharCode(...u8.subarray(0, 4)) !== 'HXLG') throw new Error(file + ': no HXLG marker');
    const hl = new DataView(ab).getUint32(8, true), H = JSON.parse(new TextDecoder().decode(u8.subarray(12, 12 + hl)));
    const start = (12 + hl + 3) & ~3, types = {Int32Array, Uint32Array, Uint16Array, Int8Array}, D = {H};
    for (const s of H.sections) D[s.name] = new types[s.type](ab, start + s.offset, s.length);
    D.genome = id => {
        const g = new Int8Array(H.ncond).fill(-1), ch = [];
        for (let a = id; a >= 0; a = D.parent[a]) ch.push(a);
        for (let k = ch.length - 1; k >= 0; k--) for (let d = D.deltaStart[ch[k]]; d < D.deltaStart[ch[k] + 1]; d++) g[D.deltaCond[d]] = D.deltaVal[d];
        let s = ''; const sym = '.' + 'RGBYCM'.slice(0, H.n); for (let c = 0; c < H.ncond; c++) s += sym[g[c] + 1];
        return s;
    };
    D.L = {parent: D.parent, founder: D.founder, birth: D.birth, death: D.death, tick: H.exportedAtTick};
    return D;
}

// ================================================================ 2. lineage on the late .hxlg
let late = null;
if (!fs.existsSync(LATE)) console.log(`skip 2: ${LATE} not found`);
else {
    let rt; [late, rt] = ms(() => readHxlg(LATE));
    const H = late.H, L = late.L, T = H.exportedAtTick;
    console.log(`2. ${path.basename(LATE)}: ${H.organisms.toLocaleString()} organisms ever, read in ${(rt / 1000).toFixed(1)} s`);
    const live = H.liveAtExport.filter(id => id >= 0 && id < H.organisms);
    const [alive, ta] = ms(() => GA.aliveAt(L, T));
    const A = new Set(alive), B = new Set(live);
    const onlyAlive = alive.filter(id => !B.has(id)), onlyLive = live.filter(id => !A.has(id));
    check(!onlyAlive.length && !onlyLive.length, 'aliveAt == liveAtExport');
    console.log(`   aliveAt(${T}): ${alive.length} vs liveAtExport ${live.length} (${H.liveAtExport.length - live.length} out of range dropped); ` +
        `only in aliveAt ${onlyAlive.length}, only in liveAtExport ${onlyLive.length} (${ta.toFixed(0)} ms)` + (onlyAlive.length ? ` e.g. ${onlyAlive.slice(0, 5)}` : ''));

    // clusters at t = 4 from the living genomes, as fromHxlg does
    const byG = new Map();
    for (const id of live) { const s = late.genome(id), e = byG.get(s); if (e) { e.count++; e.ids.push(id); } else byG.set(s, {count: 1, ids: [id]}); }
    const byKey = new Map();
    for (const [g, e] of byG) {
        const key = g.replace(/[^.]/g, '1').replace(/\./g, '0');
        let p = byKey.get(key); if (!p) byKey.set(key, p = {key, bits: gc.presenceBits(g), count: 0, exact: []});
        p.count += e.count; p.exact.push([g, e.count, e.ids]);
    }
    const pats = [...byKey.values()];
    for (const p of pats) p.exact.sort((a, b) => b[1] - a[1]);
    pats.sort((a, b) => b.count - a.count || (a.key < b.key ? -1 : 1));
    const cl = clustersAt(pats, gc.spanningTree(pats), 4);

    // the explorer's clusterLineage, verbatim logic
    function explorerLineage(ids) {
        const founders = new Set(ids.map(id => L.founder[id]));
        if (founders.size > 1 || ids.length < 2) return {founders: founders.size, ancestor: null, ageTicks: null, generations: 0};
        const chain = [], pos = new Map();
        for (let a = ids[0]; a >= 0; a = L.parent[a]) { pos.set(a, chain.length); chain.push(a); }
        let deepest = 0; const hits = [];
        for (let k = 1; k < ids.length; k++) {
            let a = ids[k], steps = 0;
            while (a >= 0 && !pos.has(a)) { a = L.parent[a]; steps++; }
            if (a < 0) return {founders: founders.size, ancestor: null, ageTicks: null, generations: 0};
            const at = pos.get(a); hits.push(steps - at); if (at > deepest) deepest = at;
        }
        const generations = Math.max(deepest, ...hits.map(h => h + deepest));
        return {founders: 1, ancestor: chain[deepest], ageTicks: T - L.birth[chain[deepest]], generations};
    }
    // brute force: ancestors (inclusive) common to every member; the MRCA is the deepest of them
    function bruteMRCA(ids) {
        const depth = id => { let d = 0; for (let a = L.parent[id]; a >= 0; a = L.parent[a]) d++; return d; };
        let common = new Set(); for (let a = ids[0]; a >= 0; a = L.parent[a]) common.add(a);
        for (const id of ids.slice(1)) { const s = new Set(); for (let a = id; a >= 0; a = L.parent[a]) if (common.has(a)) s.add(a); common = s; }
        if (!common.size) return null;
        let best = -1, bd = -1; for (const a of common) { const d = depth(a); if (d > bd) { bd = d; best = a; } }
        return {id: best, ageTicks: T - L.birth[best], maxGenerations: Math.max(...ids.map(depth)) - bd};
    }
    let multi = 0, agreeExp = 0, agreeBrute = 0, withMRCA = 0, single = 0;
    for (const c of cl) {
        const entries = GA.entriesFromPatterns(pats, c.mem), ids = entries.flatMap(e => e.ids);
        const R = GA.lineageRelations(L, entries, H.n, {maxOrganisms: 2});
        const br = bruteMRCA(ids);
        const bOK = JSON.stringify(R.clusterMRCA) === JSON.stringify(br);
        check(bOK, `cluster of ${ids.length}: GA ${JSON.stringify(R.clusterMRCA)} vs brute ${JSON.stringify(br)}`);
        agreeBrute += bOK;
        if (ids.length < 2) { single++; continue; }
        multi++;
        const ex = explorerLineage(ids), g = R.clusterMRCA;
        const eOK = ex.ancestor === null ? g === null : (g && g.id === ex.ancestor && g.ageTicks === ex.ageTicks && g.maxGenerations === ex.generations);
        check(eOK, `cluster of ${ids.length}: GA ${JSON.stringify(g)} vs explorer ${JSON.stringify(ex)}`);
        agreeExp += !!eOK; if (g) withMRCA++;
    }
    console.log(`   t=4: ${cl.length} clusters; clusterMRCA == brute force for ${agreeBrute}/${cl.length}; ` +
        `== explorer clusterLineage for ${agreeExp}/${multi} clusters with ≥2 organisms (${withMRCA} have a common ancestor; ${single} singletons, where the explorer reports none and GA reports the organism itself)`);

    // 2,000 random pairs: binary-search prefix vs naive ancestor-set walk
    let pairOK = 0;
    for (let k = 0; k < 2000; k++) {
        const a = live[Math.floor(rand() * live.length)], b = live[Math.floor(rand() * live.length)];
        const r = GA.pairRelation(L, GA.chainOf(L, a), GA.chainOf(L, b));
        const anc = new Map(); let d = 0; for (let x = a; x >= 0; x = L.parent[x]) anc.set(x, d++);
        let m = -1, gb = 0; for (let x = b; x >= 0; x = L.parent[x], gb++) if (anc.has(x)) { m = x; break; }
        const want = m < 0 ? {mrca: null} : {mrca: m, genA: anc.get(m), genB: gb, ageTicks: T - L.birth[m]};
        const good = m < 0 ? r.mrca === null : (r.mrca === m && r.genA === want.genA && r.genB === want.genB && r.ageTicks === want.ageTicks);
        check(good, `pair ${a},${b}: ${JSON.stringify(r)} vs ${JSON.stringify(want)}`); pairOK += good;
    }
    console.log(`   random pair MRCAs: ${pairOK}/2000 equal the naive ancestor walk`);
    const bigE = GA.entriesFromPatterns(pats, cl[0].mem);
    const [R, tl] = ms(() => GA.lineageRelations(L, bigE, H.n, {maxOrganisms: 300}));
    console.log(`   lineageRelations on the largest t=4 cluster (${cl[0].organisms} organisms, ${R.considered} used): ${tl.toFixed(0)} ms; ` +
        `MRCA ${JSON.stringify(R.clusterMRCA)}, same founder ${(R.sameFounderShare * 100).toFixed(1)}%, spearman(presence, tree) ${R.spearmanPresenceVsTree?.toFixed(3)}`);
    // t = 4 clusters are tiny here, so also time a large cluster (t = 16) and the whole living population
    const mst = gc.spanningTree(pats), c16 = clustersAt(pats, mst, 16)[0];
    for (const [label, mem] of [['largest t=16 cluster', c16.mem], ['whole population', pats.map((_, i) => i)]]) {
        const en = GA.entriesFromPatterns(pats, mem), org = en.reduce((s, e) => s + e.count, 0);
        const [Rx, tx] = ms(() => GA.lineageRelations(L, en, H.n, {maxOrganisms: 300}));
        const [, tcx] = ms(() => GA.conservation(en, H.ncond, H.n));
        const [, tkx] = ms(() => GA.compactness(en, H.n, en[0].g));
        const bx = bruteMRCA(en.flatMap(e => e.ids));
        check(JSON.stringify(Rx.clusterMRCA) === JSON.stringify(bx), label + ' MRCA');
        console.log(`   ${label} (${org} organisms): lineageRelations(maxOrganisms 300, sampled ${Rx.sampled}) ${tx.toFixed(0)} ms, ` +
            `conservation ${tcx.toFixed(0)} ms, compactness ${tkx.toFixed(0)} ms; MRCA ${JSON.stringify(Rx.clusterMRCA)} == brute ${JSON.stringify(Rx.clusterMRCA) === JSON.stringify(bx)}`);
    }
}

// ================================================================ 3. temporal: early file vs the late one
if (!late || !fs.existsSync(EARLY)) console.log(`skip 3: ${!late ? 'late file missing' : EARLY + ' not found'}`);
else {
    const E = readHxlg(EARLY), T0 = E.H.exportedAtTick, n = Math.min(E.H.organisms, late.H.organisms);
    let pDiff = 0, fDiff = 0, bDiff = 0, dDiff = 0, firstDiff = -1;
    for (let id = 0; id < n; id++) {
        const p = E.parent[id] !== late.parent[id], f = E.founder[id] !== late.founder[id], b = E.birth[id] !== late.birth[id];
        pDiff += p; fDiff += f; bDiff += b;
        if ((p || f || b) && firstDiff < 0) firstDiff = id;
        if (E.death[id] >= 0 && E.death[id] !== late.death[id]) dDiff++;   // a death recorded early must persist
    }
    const same = !pDiff && !fDiff && !bDiff;
    console.log(`3. ${path.basename(EARLY)} (tick ${T0}, ${E.H.organisms.toLocaleString()} organisms) vs late file, common ids 0..${n - 1}: ` +
        `parent differs ${pDiff}, founder ${fDiff}, birth ${bDiff}, recorded deaths ${dDiff}` + (firstDiff >= 0 ? `, first difference at id ${firstDiff}` : ''));
    if (!same) console.log('   NOT the same simulation (or ids were renumbered): no temporal reconstruction attempted.');
    else {
        check(dDiff === 0, 'early recorded deaths persist in the late file');
        const recon = GA.aliveAt(late.L, T0), truth = E.H.liveAtExport.filter(id => id >= 0 && id < E.H.organisms);
        const S = new Set(truth), R = new Set(recon);
        const onlyR = recon.filter(id => !S.has(id)), onlyT = truth.filter(id => !R.has(id));
        check(!onlyR.length && !onlyT.length, 'alive set at early tick reconstructed');
        let gBad = 0, gN = 0;
        for (const id of truth) if (R.has(id)) { gN++; if (late.genome(id) !== E.genome(id)) gBad++; }
        check(gBad === 0, 'rebuilt genomes equal');
        console.log(`   same simulation. aliveAt(late, ${T0}) = ${recon.length} vs early liveAtExport ${truth.length}: ` +
            `only reconstructed ${onlyR.length}, only in early file ${onlyT.length}; genomes rebuilt from both files differ for ${gBad}/${gN}` +
            (onlyR.length ? `; e.g. ${onlyR.slice(0, 3).map(id => `${id} born ${late.birth[id]} died ${late.death[id]} (early death ${E.death[id]})`).join(', ')}` : ''));
    }
}
console.log(`\n${checks - fails}/${checks} checks passed`);
if (fails) process.exitCode = 1;
