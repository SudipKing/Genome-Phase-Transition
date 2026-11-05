// In-depth cross-run analysis of a batch of run.js outputs.
//
// usage: node analysis/deep.js <batch dir> > deep.json
//        node analysis/deep.js <run dir> <run dir> ... > deep.json
//
// For a batch directory written by batch.js, runs are grouped by their parameter set
// (batch.json); otherwise every run is in set "runs". For every run it writes:
//   - phases (phases.js) and a phase class per phase: "3-colour churn", "2-colour churn",
//     "2-colour full", ...; fillTick = start of the final full-grid phase, if any
//   - per-phase summaries: population, genome size, the consensus rules split by number of
//     empty neighbours, where cells meet conditions and where rule deaths happen, genome
//     diversity and clade sweeps
//   - event timing relative to the fill: when the third colour is lost, when the consensus
//     becomes mostly rules for fully surrounded cells, when genome size has fallen halfway,
//     when rule deaths collapse, the biggest rule-space velocity spike and the last clade
//     sweep (one clade holding at least half the cells at the end of a clade window)
//   - compact time series for the explorer (explorer.js)
// Across runs:
//   - per parameter set: how many fill, when, extinctions, endpoints and phase sequences
//   - per phase class (pooled over runs): the same summaries
//   - genome similarity of phase consensus genomes (best of the six colour relabelings),
//     within a phase class across runs vs across phase classes, to test whether phases
//     have characteristic genomes
//   - fill-aligned series: each metric at ticks relative to the fill

const fs = require('fs');
const path = require('path');
const {segment} = require('./phases.js');
const {colourSpace, endpointBox} = require('./colours.js');

const inputs = process.argv.slice(2);
if (!inputs.length) { console.error('usage: node deep.js <batch dir> | <run dir> ...'); process.exit(1); }

let runDirs = [], setOf = new Map(), setOrder = [], setParams = {};
if (inputs.length === 1 && fs.existsSync(path.join(inputs[0], 'batch.json'))) {
    const b = JSON.parse(fs.readFileSync(path.join(inputs[0], 'batch.json'), 'utf8'));
    for (const s of b.plan.sets) { setOrder.push(s.name); setParams[s.name] = s.params || {}; }
    for (const r of b.runs) {
        const dir = path.join(inputs[0], r.dir);
        if (!fs.existsSync(path.join(dir, 'samples.json'))) continue;
        runDirs.push(dir); setOf.set(dir, r.set);
    }
} else {
    runDirs = inputs; setOrder = ['runs']; setParams.runs = {};
    for (const d of inputs) setOf.set(d, 'runs');
}

const mean = a => a.length ? a.reduce((s, v) => s + v, 0) / a.length : null;
const median = a => { if (!a.length) return null; const s = [...a].sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
const quart = a => { if (!a.length) return null; const s = [...a].sort((x, y) => x - y); const q = p => s[Math.min(s.length - 1, Math.floor(p * (s.length - 1) + 0.5))]; return [q(0.25), q(0.5), q(0.75)]; };
const r3 = v => v == null ? null : +v.toFixed(3);
const sum = a => a.reduce((s, v) => s + v, 0);

function phaseClass(lab) {
    if (lab === 'extinct') return 'extinct';
    const [cols, full] = lab.split(' ');
    return `${cols.length}-colour ${full ? 'full' : 'churn'}`;
}

let cs = null, conds = null, nDeadOf = null;
const runs = [];
for (const dir of runDirs) {
    const d = JSON.parse(fs.readFileSync(path.join(dir, 'samples.json'), 'utf8'));
    const meta = fs.existsSync(path.join(dir, 'run.json')) ? JSON.parse(fs.readFileSync(path.join(dir, 'run.json'), 'utf8')) : {};
    if (!cs) {
        conds = d.conditions; cs = colourSpace(conds, d.params.n);
        nDeadOf = conds.map(c => 6 - sum(c.slice(1)));
    }
    const S = d.samples;
    if (!S.length || S[0].consByDead === undefined) { console.error(`${dir}: no samples or missing new metrics, skipped`); continue; }
    runs.push({name: path.basename(dir), set: setOf.get(dir), seed: d.seed, ticks: d.ticks, stopReason: meta.stopReason || null, S, cladeEvery: meta.cladeEvery || 5000});
}
if (!runs.length) { console.error('no runs'); process.exit(1); }
const ncond = cs.ncond;

// Share of consensus rules for fully surrounded cells (0 empty neighbours).
const surroundShare = s => { const t = sum(s.consByDead); return t ? s.consByDead[0] / t : 0; };
const shareOf = arr => { const t = sum(arr); return arr.map(v => t ? v / t : 0); };

// Similarity of two consensus genomes (Uint8 arrays): rules agreeing / conditions either has a rule for.
function sim(a, b) {
    let u = 0, same = 0;
    for (let c = 0; c < ncond; c++) {
        if (!a[c] && !b[c]) continue;
        u++; if (a[c] === b[c]) same++;
    }
    return u ? same / u : 1;
}
function bestSim(a, b) {
    let best = 0;
    for (let pi = 0; pi < cs.perms.length; pi++) best = Math.max(best, sim(cs.relabelGenome(a, pi), b));
    return best;
}
// Majority output per condition over several consensus strings.
function modeGenome(strs) {
    const out = new Uint8Array(ncond);
    for (let c = 0; c < ncond; c++) {
        const cnt = [0, 0, 0, 0, 0, 0, 0];
        for (const g of strs) cnt[g[c] === '.' ? 0 : cs.SYM.indexOf(g[c]) + 1]++;
        let best = 0; for (let o = 1; o < cnt.length; o++) if (cnt[o] > cnt[best]) best = o;
        out[c] = best;
    }
    return out;
}

function summarise(ss) {
    const settled = ss.slice(Math.floor(ss.length / 2));
    const m = f => r3(mean(settled.map(f)));
    const vecMean = f => { const k = f(settled[0]).length; const acc = new Array(k).fill(0); for (const s of settled) f(s).forEach((v, i) => acc[i] += v / settled.length); return acc.map(r3); };
    return {
        pop: m(s => s.pop), meanGenome: m(s => s.meanGenome), consensusSize: m(s => s.consensusSize),
        genomeSizeQ: vecMean(s => s.genomeSizeQ),
        consByDead: vecMean(s => shareOf(s.consByDead)), encByDead: vecMean(s => shareOf(s.encByDead)),
        deathByDead: vecMean(s => s.deathByDead),
        surroundShare: m(surroundShare), selfLoop: m(s => s.selfLoop),
        births: m(s => s.births), ruleDeaths: m(s => s.ruleDeaths), randomDeaths: m(s => s.randomDeaths),
        shannonGenomes: m(s => s.shannonGenomes), distinctGenomes: m(s => s.distinctGenomes),
        consensusAgreement: m(s => s.consensusAgreement), meanDistFromConsensus: m(s => s.meanDistFromConsensus),
        velocity: r3(mean(ss.map(s => s.velocity || 0)) * 1000) / 1000,
        cluster: m(s => s.cluster || 0),
    };
}

// First index at/after `from` (stepping by dir) where pred holds for `hold` consecutive samples.
function firstSustained(S, from, pred, hold = 3, dir = 1) {
    for (let i = from; i >= 0 && i < S.length; i += dir) {
        let ok = true;
        for (let k = 0; k < hold && ok; k++) { const j = i + k * dir; if (j < 0 || j >= S.length || !pred(S[j], j)) ok = false; }
        if (ok) return i;
    }
    return -1;
}

const out = {conditions: conds, nDead: nDeadOf, sets: [], runs: [], phaseClasses: [], similarity: null, aligned: null, events: null};

for (const r of runs) {
    const S = r.S, ticks = S.map(s => s.tick);
    const phases = segment(S).map(p => {
        const ss = S.slice(p.start, p.end + 1);
        return {label: p.label, cls: phaseClass(p.label), start: p.start, end: p.end, startTick: ss[0].tick, endTick: ss[ss.length - 1].tick,
            summary: summarise(ss), consensus: modeGenome(ss.slice(Math.floor(ss.length / 2)).map(s => s.consensus))};
    });
    const last = phases[phases.length - 1];
    const filled = last.cls.endsWith('full');
    const fillIdx = filled ? last.start : -1;
    const fillTick = filled ? last.startTick : null;
    const box = endpointBox(S[S.length - 1], 0.08);
    const canon = cs.canonicalPerm(box.order);

    // Clade windows: the sample at the end of each window holds that window's clade shares.
    const windows = [];
    for (let i = 0; i < S.length; i++) {
        const nextStart = i + 1 < S.length ? S[i + 1].cladeWindowStart : null;
        if (nextStart !== S[i].cladeWindowStart || i === S.length - 1)
            windows.push({start: S[i].cladeWindowStart, end: S[i].tick, topShare: r3(S[i].topCladeShare), clades: S[i].cladesAlive, cls: phaseClass(phases.find(p => i >= p.start && i <= p.end).label)});
    }

    // Event timing, relative to the fill.
    let events = null;
    if (filled) {
        const evt = {};
        const relTick = i => i < 0 ? null : S[i].tick - fillTick;
        // Lost the third colour: last time 3 colours were present before the fill.
        let lastThree = -1;
        for (let i = 0; i < fillIdx; i++) if (S[i].colors.filter(c => c > 0.08).length >= 3) lastThree = i;
        evt.thirdColourLost = lastThree >= 0 ? relTick(Math.min(lastThree + 1, S.length - 1)) : null;
        // Consensus becomes mostly (>50%) rules for fully surrounded cells.
        const iSur = firstSustained(S, Math.max(0, lastThree), s => surroundShare(s) > 0.5);
        evt.surroundMajority = relTick(iSur);
        // Genome size halfway from its two-colour-churn level to its final level.
        const pre = median(S.slice(Math.max(0, lastThree + 1), Math.max(lastThree + 2, fillIdx - 2)).map(s => s.meanGenome));
        const post = median(S.slice(-5).map(s => s.meanGenome));
        if (pre != null && pre > post) { const half = (pre + post) / 2; evt.genomeHalfway = relTick(firstSustained(S, Math.max(0, lastThree), s => s.meanGenome < half)); }
        // Rule deaths fall below 10% of their two-colour-churn level.
        const preD = median(S.slice(Math.max(0, lastThree + 1), Math.max(lastThree + 2, fillIdx - 2)).map(s => s.ruleDeaths));
        if (preD) evt.ruleDeathsCollapse = relTick(firstSustained(S, Math.max(0, lastThree), s => s.ruleDeaths < 0.1 * preD));
        // Population first reaches 80%.
        evt.pop80 = relTick(firstSustained(S, 0, s => s.pop >= 0.8));
        // Biggest velocity spike within 30k ticks of the fill.
        let vi = -1;
        for (let i = 0; i < S.length; i++) if (Math.abs(S[i].tick - fillTick) <= 30000 && (vi < 0 || (S[i].velocity || 0) > (S[vi].velocity || 0))) vi = i;
        evt.velocityPeak = relTick(vi);
        // Clade sweeps (top clade ≥ 50% at window end); the one nearest the fill.
        const sweeps = windows.filter(w => w.topShare >= 0.5);
        const near = sweeps.length ? sweeps.reduce((a, b) => Math.abs(b.end - fillTick) < Math.abs(a.end - fillTick) ? b : a) : null;
        evt.nearestSweep = near ? near.end - fillTick : null;
        evt.lastFounder = relTick(S.findIndex(s => s.foundersAlive === 1));
        events = evt;
    }

    // Compact series for the explorer.
    const col = f => S.map(f);
    const series = {
        tick: ticks,
        pop: col(s => r3(s.pop)), colors: [0, 1, 2].map(k => col(s => r3(s.colors[k]))),
        meanGenome: col(s => +s.meanGenome.toFixed(1)), genomeSizeQ: col(s => s.genomeSizeQ), consensusSize: col(s => s.consensusSize),
        consByDead: col(s => s.consByDead), encByDead: col(s => shareOf(s.encByDead).map(r3)), deathByDead: col(s => s.deathByDead),
        surroundShare: col(s => r3(surroundShare(s))),
        births: col(s => +(s.births || 0).toFixed(1)), ruleDeaths: col(s => +(s.ruleDeaths || 0).toFixed(1)), randomDeaths: col(s => +(s.randomDeaths || 0).toFixed(1)),
        shannonGenomes: col(s => s.shannonGenomes), distinctGenomes: col(s => s.distinctGenomes), effGenomes: col(s => +(s.effGenomes || 0).toFixed(1)),
        velocity: col(s => +((s.velocity || 0)).toFixed(5)), entropy: col(s => r3(s.entropy || 0)), pi: col(s => r3(s.pi || 0)),
        selfLoop: col(s => r3(s.selfLoop || 0)), cluster: col(s => r3(s.cluster || 0)),
        consensusAgreement: col(s => r3(s.consensusAgreement)), meanDistFromConsensus: col(s => +s.meanDistFromConsensus.toFixed(1)),
        cladesAlive: col(s => s.cladesAlive), topCladeShare: col(s => r3(s.topCladeShare)), foundersAlive: col(s => s.foundersAlive),
    };
    // Consensus strips, raw and relabeled so the run's dominant final colour is R.
    const strips = S.map(s => s.consensus);
    const stripsCanon = strips.map(g => cs.decode(cs.relabelGenome(cs.encode(g), canon)));

    r.phases = phases; r.filled = filled; r.fillTick = fillTick;
    out.runs.push({
        name: r.name, set: r.set, seed: r.seed, ticks: r.ticks, stopReason: r.stopReason,
        filled, fillTick, endpoint: box.label, canonPerm: cs.permName(canon),
        extinct: S[S.length - 1].pop < 0.01,
        phaseSequence: phases.map(p => p.label).join(' → '),
        classSequence: phases.map(p => p.cls).filter((c, i, a) => i === 0 || a[i - 1] !== c).join(' → '),
        phases: phases.map(({consensus, ...p}) => ({...p, consensus: cs.decode(consensus)})),
        windows, events, series, strips, stripsCanon,
    });
}

// ── Per parameter set ────────────────────────────────────────────────────────
for (const name of setOrder) {
    const rs = out.runs.filter(r => r.set === name);
    if (!rs.length) continue;
    const fills = rs.filter(r => r.filled).map(r => r.fillTick);
    const counts = key => { const m = {}; for (const r of rs) m[r[key]] = (m[r[key]] || 0) + 1; return Object.entries(m).sort((a, b) => b[1] - a[1]); };
    const durOf = cls => rs.map(r => sum(r.phases.filter(p => p.cls === cls).map(p => p.endTick - p.startTick + 1000)));
    const finalPh = rs.map(r => r.phases[r.phases.length - 1]);
    out.sets.push({
        name, params: setParams[name], runs: rs.length, filled: fills.length, extinct: rs.filter(r => r.extinct).length,
        fillTicks: fills, fillQuartiles: quart(fills),
        threeColourTicks: quart(durOf('3-colour churn')), twoColourChurnTicks: quart(durOf('2-colour churn')),
        finalGenome: quart(finalPh.map(p => p.summary.meanGenome)),
        finalSurround: quart(finalPh.map(p => p.summary.surroundShare)),
        endpoints: counts('endpoint'), classSequences: counts('classSequence'),
    });
}

// ── Per phase class, pooled over the first set's runs (the baseline) ─────────────────
const baseSet = out.sets[0].name;
const baseRuns = runs.filter(r => r.set === baseSet);
const classes = ['3-colour churn', '2-colour churn', '2-colour full', '1-colour churn', '1-colour full', '3-colour full'];
for (const cls of classes) {
    const ph = baseRuns.flatMap(r => r.phases.filter(p => p.cls === cls).map(p => ({...p, run: r.name})));
    if (!ph.length) continue;
    const pick = f => r3(mean(ph.map(p => f(p.summary))));
    const vec = f => { const k = f(ph[0].summary).length; const a = new Array(k).fill(0); for (const p of ph) f(p.summary).forEach((v, i) => a[i] += v / ph.length); return a.map(r3); };
    const wins = out.runs.filter(r => r.set === baseSet).flatMap(r => r.windows.filter(w => w.cls === cls));
    out.phaseClasses.push({
        cls, phases: ph.length, runs: new Set(ph.map(p => p.run)).size,
        durationTicks: quart(ph.map(p => p.endTick - p.startTick + 1000)),
        pop: pick(s => s.pop), meanGenome: pick(s => s.meanGenome), consensusSize: pick(s => s.consensusSize),
        genomeSizeQ: vec(s => s.genomeSizeQ),
        consByDead: vec(s => s.consByDead), encByDead: vec(s => s.encByDead), deathByDead: vec(s => s.deathByDead),
        surroundShare: pick(s => s.surroundShare), selfLoop: pick(s => s.selfLoop),
        births: pick(s => s.births), ruleDeaths: pick(s => s.ruleDeaths), shannonGenomes: pick(s => s.shannonGenomes),
        consensusAgreement: pick(s => s.consensusAgreement), meanDistFromConsensus: pick(s => s.meanDistFromConsensus),
        cluster: pick(s => s.cluster),
        cladeWindows: wins.length, sweepShare: wins.length ? r3(wins.filter(w => w.topShare >= 0.5).length / wins.length) : null,
        topCladeShare: quart(wins.map(w => w.topShare)),
    });
}

// ── Genome similarity: same phase class across runs vs different classes ───────────
{
    const items = baseRuns.flatMap(r => r.phases.filter(p => p.cls !== 'extinct' && p.end - p.start + 1 >= 5).map(p => ({run: r.name, cls: p.cls, g: p.consensus})));
    const cls = [...new Set(items.map(i => i.cls))].filter(c => classes.includes(c)).sort((a, b) => classes.indexOf(a) - classes.indexOf(b));
    const M = cls.map(() => cls.map(() => []));
    const withinRunCross = [];
    for (let i = 0; i < items.length; i++)
        for (let j = i + 1; j < items.length; j++) {
            const a = items[i], b = items[j];
            const ai = cls.indexOf(a.cls), bi = cls.indexOf(b.cls);
            if (ai < 0 || bi < 0) continue;
            const s = bestSim(a.g, b.g);
            if (a.run === b.run) { if (a.cls !== b.cls) withinRunCross.push({run: a.run, a: a.cls, b: b.cls, sim: r3(s)}); continue; }
            M[ai][bi].push(s); if (ai !== bi) M[bi][ai].push(s);
        }
    out.similarity = {
        classes: cls,
        matrix: M.map(row => row.map(v => v.length ? {mean: r3(mean(v)), q: quart(v).map(r3), n: v.length} : null)),
        withinRunAcrossPhases: withinRunCross,
        note: 'Phase consensus genomes compared over conditions either genome has a rule for, best of the six colour relabelings. Off-diagonal and diagonal cells use pairs from different runs only.',
    };
    // Shared rules of each phase class: rules held by the relabeled phase consensus in at least 80% of runs.
    out.similarity.coreRules = cls.map(c => {
        const gs = baseRuns.map(r => { const p = r.phases.filter(x => x.cls === c).pop(); if (!p) return null; const box = endpointBox(r.S[p.end], 0.08); return cs.relabelGenome(p.consensus, cs.canonicalPerm(box.order)); }).filter(Boolean);
        const rules = [];
        for (let k = 0; k < ncond; k++) {
            const cnt = {}; for (const g of gs) if (g[k]) cnt[g[k]] = (cnt[g[k]] || 0) + 1;
            for (const [o, v] of Object.entries(cnt)) if (v >= 0.8 * gs.length) rules.push({cond: k, nDead: nDeadOf[k], out: cs.SYM[o - 1], share: r3(v / gs.length), text: cs.ruleText(k, cs.SYM[o - 1])});
        }
        return {cls: c, runs: gs.length, rules};
    });
}

// ── Event timing aggregated and fill-aligned series ──────────────────────────────
{
    const filled = out.runs.filter(r => r.set === baseSet && r.filled);
    const keys = ['thirdColourLost', 'pop80', 'surroundMajority', 'genomeHalfway', 'ruleDeathsCollapse', 'velocityPeak', 'nearestSweep', 'lastFounder'];
    out.events = {runs: filled.length, keys: keys.map(k => {
        const v = filled.map(r => r.events[k]).filter(x => x != null);
        return {key: k, n: v.length, q: quart(v), before: v.filter(x => x < 0).length, values: filled.map(r => ({run: r.name, v: r.events[k]}))};
    })};
    const offsets = []; for (let o = -60000; o <= 20000; o += 1000) offsets.push(o);
    const metrics = ['pop', 'meanGenome', 'consensusSize', 'surroundShare', 'ruleDeaths', 'births', 'shannonGenomes', 'velocity', 'topCladeShare', 'selfLoop', 'consensusAgreement', 'cluster'];
    const aligned = {offsets, sets: {}};
    for (const set of out.sets) {
        const rs = out.runs.filter(r => r.set === set.name && r.filled);
        aligned.sets[set.name] = {};
        for (const m of metrics) {
            const perRun = rs.map(r => offsets.map(o => { const i = r.series.tick.indexOf(r.fillTick + o); return i < 0 ? null : r.series[m][i]; }));
            const med = offsets.map((_, k) => { const v = perRun.map(a => a[k]).filter(x => x != null); return v.length >= 3 ? +median(v).toFixed(4) : null; });
            aligned.sets[set.name][m] = {median: med, runs: rs.map(r => r.name), perRun};
        }
    }
    out.aligned = aligned;
}

process.stdout.write(JSON.stringify(out));
const fmt = q => q ? q.map(v => v == null ? '–' : Math.round(v)).join('/') : '–';
for (const s of out.sets) console.error(`${s.name.padEnd(12)} runs ${s.runs} filled ${s.filled} extinct ${s.extinct} fill q1/med/q3 ${fmt(s.fillQuartiles)}  endpoints ${s.endpoints.map(e => e.join(':')).join(' ')}`);
for (const p of out.phaseClasses) console.error(`${p.cls.padEnd(15)} n=${p.phases} pop ${p.pop} genome ${p.meanGenome} surround ${p.surroundShare} sweeps ${p.sweepShare}`);
if (out.similarity) console.error('similarity', out.similarity.classes.join(' | '), JSON.stringify(out.similarity.matrix.map(r => r.map(c => c && c.mean))));
for (const e of out.events.keys) console.error(`${e.key.padEnd(20)} n=${e.n} q ${fmt(e.q)} before fill ${e.before}`);
