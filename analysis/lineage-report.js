// Lineage + genome analysis of one or more .hxlg exports.
// 1. Phases (from the stats series) vs extinctions of major clades and drops in TMRCA.
// 2. Selection vs drift: are the rule changes that founded successful clades enriched in
//    conditions cells actually encounter, compared with all rule changes at birth?
// 3. Lineage distance vs genome distance for living pairs in each phase.
// usage: node analysis/lineage-report.js run1.hxlg [run2.hxlg ...] > report.json
const fs = require('fs'), path = require('path');
const OrganismLog = require(path.join(__dirname, '..', 'organismlog.js'));
const A = require(path.join(__dirname, '..', 'lineage-analysis.js'));

function load(file) {
    const b = fs.readFileSync(file);
    return OrganismLog.parse(b.buffer.slice(b.byteOffset, b.byteOffset + b.length));
}

// Same phase rule as analyze.js and lineage-explorer.html.
function phasesOf(s) {
    const T = s.population.length, names = 'RGBYPW';
    const lab = k => s.population[k] < 0.01 ? 'extinct' : s.colors.map((c, i) => c[k] > 0.08 ? names[i] : '').join('') + (s.population[k] >= 0.95 ? ' full' : '');
    let ph = [];
    for (let k = 0; k < T; k++) { const l = lab(k); if (ph.length && ph[ph.length - 1].label === l) ph[ph.length - 1].end = k; else ph.push({label: l, start: k, end: k}); }
    let changed = true;
    while (changed && ph.length > 1) {
        changed = false;
        const k = ph.findIndex(p => p.end - p.start + 1 < 5);
        if (k >= 0) { const into = k > 0 ? k - 1 : k + 1; ph[into].start = Math.min(ph[into].start, ph[k].start); ph[into].end = Math.max(ph[into].end, ph[k].end); ph.splice(k, 1); changed = true; }
        for (let j = ph.length - 1; j > 0; j--) if (ph[j].label === ph[j - 1].label) { ph[j - 1].end = ph[j].end; ph.splice(j, 1); changed = true; }
    }
    return ph.map(p => ({label: p.label, startTick: (p.start + 1) * s.reportingPeriod, endTick: (p.end + 1) * s.reportingPeriod}));
}

const encounteredAt = (s, t) => {
    const k = Math.min(s.encounterFreq.length - 1, Math.max(0, Math.round(t / s.reportingPeriod) - 1));
    return s.encounterFreq[k].map(v => v > 0.001);
};

function analyseRun(file) {
    const d = load(file), h = d.header, s = h.series;
    const phases = phasesOf(s);
    const changes = phases.slice(1).map(p => p.startTick);
    const times = A.sampleTicks(d);
    const an = A.analyzeClades(d, times, {threshold: 0.05});
    const P = s.reportingPeriod, window = 5 * P;
    const near = t => changes.reduce((best, c) => best === null || Math.abs(c - t) < Math.abs(best - t) ? c : best, null);

    // Extinctions of clades that held >= 10%.
    const extinctions = an.clades.filter(c => c.extinctTick !== null && c.peak >= 0.1)
        .map(c => ({clade: c.id, peak: +c.peak.toFixed(3), extinctTick: c.extinctTick, nearestPhaseChange: near(c.extinctTick)}));
    extinctions.forEach(e => { e.gap = e.nearestPhaseChange === null ? null : e.extinctTick - e.nearestPhaseChange; });

    // Extinction events: nested clades on one branch die together, so group by tick.
    const byTick = new Map();
    for (const e of extinctions) {
        const g = byTick.get(e.extinctTick) || {tick: e.extinctTick, clades: 0, largestPeak: 0, nearestPhaseChange: e.nearestPhaseChange, gap: e.gap};
        g.clades++; g.largestPeak = Math.max(g.largestPeak, e.peak);
        byTick.set(e.extinctTick, g);
    }
    const extinctionEvents = [...byTick.values()].sort((a, b) => a.tick - b.tick);

    // Final-phase rulebook: majority output among living cells at export, in the
    // conditions the final phase actually encounters.
    const encFinal = s.encounterFreq[s.encounterFreq.length - 1].map(v => v > 0.001);
    const usedFinal = encFinal.map((u, c) => u ? c : -1).filter(c => c >= 0);
    const tally = usedFinal.map(() => new Map());
    for (const id of A.sample(h.liveAtExport, 400, 5)) {
        const g = d.genomeOf(id);
        usedFinal.forEach((c, i) => tally[i].set(g[c], (tally[i].get(g[c]) || 0) + 1));
    }
    const finalRule = usedFinal.map((c, i) => [...tally[i].entries()].sort((a, b) => b[1] - a[1])[0][0]);
    const finalMatch = g => usedFinal.filter((c, i) => g[c] === finalRule[i]).length;

    // Big TMRCA drops: the common ancestor of everyone jumps forward by >= 25% of elapsed
    // time. For each, compare the sweeping lineage's founder (the new MRCA) with organisms
    // alive when it was born: how many final-phase rules did each already carry?
    const sweeps = [];
    for (let k = 1; k < times.length; k++) {
        const a = an.tmrca[k - 1], b = an.tmrca[k];
        if (!(a !== null && b !== null && a - b >= 0.25 * times[k] && a - b >= 2000)) continue;
        const founder = an.mrca[k], born = d.birth[founder];
        const rivals = A.sample(A.livingAt(d, born), 300, 7).filter(id => id !== founder).map(id => finalMatch(d.genomeOf(id)));
        const mine = finalMatch(d.genomeOf(founder));
        sweeps.push({tick: times[k], tmrcaBefore: a, tmrcaAfter: b, nearestPhaseChange: near(times[k]),
            sweepFounder: founder, founderBorn: born, finalRulesUsed: usedFinal.length,
            founderMatchesFinal: mine, contemporariesMeanMatch: +(rivals.reduce((x, y) => x + y, 0) / rivals.length).toFixed(2),
            contemporariesAtLeastAsGood: rivals.filter(x => x >= mine).length, contemporariesSampled: rivals.length});
    }

    // Chance baseline for "near a phase change": share of the run within the window.
    const span = h.exportedAtTick;
    let covered = 0;
    for (let t = 0; t <= span; t += P) if (changes.some(c => Math.abs(c - t) <= window)) covered += P;
    const baseline = covered / span;

    // Selection test. Background: rule changes at birth for a random sample of births.
    const births = [];
    for (let id = 0; id < h.organisms; id++) if (d.parent[id] >= 0) births.push(id);
    const bg = {changes: 0, encountered: 0};
    for (const id of A.sample(births, 200000, 99)) {
        const enc = encounteredAt(s, d.birth[id]);
        for (let x = d.deltaStart[id]; x < d.deltaStart[id + 1]; x++) { bg.changes++; if (enc[d.deltaCond[x]]) bg.encountered++; }
    }
    // Successful clades: reached >= 50% of living cells. Count the changes between each
    // clade's founder and its parent clade's founder, judged at the clade's birth tick.
    const winners = {clades: 0, changes: 0, encountered: 0, list: []};
    for (const c of an.clades) {
        if (c.peak < 0.5 || c.parentClade < 0) continue;
        const muts = A.mutationsBetween(d, an.clades[c.parentClade].id, c.id);
        const enc = encounteredAt(s, c.birth);
        const e = muts.filter(m => enc[m.cond]).length;
        winners.clades++; winners.changes += muts.length; winners.encountered += e;
        winners.list.push({clade: c.id, birth: c.birth, peak: +c.peak.toFixed(2), changes: muts.length, encountered: e});
    }

    // Lineage vs genome distance in the last stretch of each phase.
    const distance = phases.map(p => {
        const t = p.endTick;
        const enc = encounteredAt(s, t);
        const {pairs} = A.pairDistances(d, t, {sampleSize: 80, encountered: enc});
        const ps = pairs.filter(x => x.coalescence !== null);
        if (ps.length < 20) return {phase: p.label, tick: t, pairs: ps.length};
        const fit = (key) => {          // least-squares slope per 1,000 ticks and mean
            const mx = ps.reduce((a, x) => a + x.coalescence, 0) / ps.length, my = ps.reduce((a, x) => a + x[key], 0) / ps.length;
            let sxy = 0, sxx = 0;
            for (const x of ps) { sxy += (x.coalescence - mx) * (x[key] - my); sxx += (x.coalescence - mx) ** 2; }
            return {mean: +my.toFixed(1), slopePer1k: +(sxx ? sxy / sxx * 1000 : 0).toFixed(2)};
        };
        return {phase: p.label, tick: t, pairs: ps.length, meanCoalescence: Math.round(ps.reduce((a, x) => a + x.coalescence, 0) / ps.length),
            encounteredConditions: enc.filter(Boolean).length, all: fit('genomeDist'), encountered: fit('encounteredDist')};
    });

    return {
        file: path.basename(path.dirname(file)), ticks: h.exportedAtTick, organisms: h.organisms,
        bytes: fs.statSync(file).size, changesPerOrganism: +(h.deltas / h.organisms).toFixed(2),
        phases, majorClades: an.clades.length,
        oneFounderLeftAt: times[an.tmrca.findIndex(v => v !== null)] ?? null,
        extinctionEvents, eventsNearChange: extinctionEvents.filter(e => e.gap !== null && Math.abs(e.gap) <= window).length,
        largestEvent: extinctionEvents.reduce((m, e) => !m || e.clades > m.clades ? e : m, null),
        sweeps, windowTicks: window, chanceNearChange: +baseline.toFixed(3),
        selection: {
            background: {changes: bg.changes, encounteredShare: +(bg.encountered / bg.changes).toFixed(3)},
            winners: {clades: winners.clades, changes: winners.changes, encounteredShare: winners.changes ? +(winners.encountered / winners.changes).toFixed(3) : null, list: winners.list},
        },
        distance,
        tmrca: an.tmrca.map((v, k) => [times[k], v]),
    };
}

const out = process.argv.slice(2).map(f => { const t = Date.now(); const r = analyseRun(f); console.error(r.file, 'analysed in', ((Date.now() - t) / 1000).toFixed(1), 's'); return r; });
process.stdout.write(JSON.stringify(out, null, 1));
