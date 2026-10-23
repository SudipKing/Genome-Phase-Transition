// Cross-run analysis of runs/*/samples.json produced by run.js.
// Segments each run into phases by which cell colours are present, then compares
// the dominant genome of each phase within and across runs.
//
// usage: node analyze.js runs/explore1 runs/r2 ... > report-data.json

const fs = require('fs');
const path = require('path');

const COLOR_NAMES = 'RGB';
const PRESENT = 0.08;     // a colour "is present" in a phase above this share of living cells
const MIN_PHASE = 5;      // samples; shorter blips merge into the neighbouring phase

function load(dir) {
    const d = JSON.parse(fs.readFileSync(path.join(dir, 'samples.json'), 'utf8'));
    d.name = path.basename(dir);
    return d;
}

function label(s) {
    if (s.pop < 0.01) return 'extinct';
    const colours = s.colors.map((f, i) => f > PRESENT ? COLOR_NAMES[i] : '').join('') || '?';
    // Same colours can be a churning population or a filled grid; those are different phases.
    return colours + (s.pop >= 0.95 ? ' full' : '');
}

function segment(samples) {
    let raw = [];
    for (let i = 0; i < samples.length; i++) {
        const l = label(samples[i]);
        if (raw.length && raw[raw.length - 1].label === l) raw[raw.length - 1].end = i;
        else raw.push({label: l, start: i, end: i});
    }
    // Merge short blips into the previous phase until stable.
    let changed = true;
    while (changed && raw.length > 1) {
        changed = false;
        for (let k = 0; k < raw.length; k++) {
            const p = raw[k];
            if (p.end - p.start + 1 >= MIN_PHASE) continue;
            const into = k > 0 ? k - 1 : k + 1;
            raw[into].start = Math.min(raw[into].start, p.start);
            raw[into].end = Math.max(raw[into].end, p.end);
            raw.splice(k, 1);
            changed = true;
            break;
        }
        // Re-join neighbours that now share a label.
        for (let k = raw.length - 1; k > 0; k--) {
            if (raw[k].label === raw[k - 1].label) { raw[k - 1].end = raw[k].end; raw.splice(k, 1); changed = true; }
        }
    }
    return raw;
}

const mean = a => a.reduce((s, v) => s + v, 0) / (a.length || 1);
const hamming = (a, b) => { let d = 0; for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) d++; return d; };

function mode(arr) {
    const m = new Map();
    for (const v of arr) m.set(v, (m.get(v) || 0) + 1);
    return [...m.entries()].sort((a, b) => b[1] - a[1])[0];
}

function summarisePhase(run, ph) {
    const ss = run.samples.slice(ph.start, ph.end + 1);
    // Use the back half of the phase so the summary reflects its settled state.
    const settled = ss.slice(Math.floor(ss.length / 2));
    const [domGenome, domCount] = mode(settled.map(s => s.consensus));
    // Mean encounter rate per condition over the settled samples: which rules actually fire.
    const enc = new Array(settled[0].enc ? settled[0].enc.length : 0).fill(0);
    for (const s of settled) if (s.enc) s.enc.forEach((v, c) => enc[c] += v / settled.length);
    return {
        label: ph.label,
        startTick: ss[0].tick, endTick: ss[ss.length - 1].tick,
        samples: ss.length,
        pop: mean(settled.map(s => s.pop)),
        colors: [0, 1, 2].map(c => mean(settled.map(s => s.colors[c]))),
        meanGenome: mean(settled.map(s => s.meanGenome)),
        distinctGenomes: mean(settled.map(s => s.distinctGenomes)),
        effGenomes: mean(settled.map(s => s.effGenomes)),
        domGenomeShare: mean(settled.map(s => s.domGenomeShare)),
        foundersAlive: settled[settled.length - 1].foundersAlive,
        ruleDeaths: mean(settled.map(s => s.ruleDeaths)),
        births: mean(settled.map(s => s.births)),
        velocity: mean(ss.map(s => s.velocity)),
        domGenome, domGenomeStability: domCount / settled.length,
        selfLoop: mean(settled.map(s => s.selfLoop)),
        consensusAgreement: mean(settled.map(s => s.consensusAgreement)),
        meanDistFromConsensus: mean(settled.map(s => s.meanDistFromConsensus)),
        enc: enc.map(v => +v.toFixed(5)),
    };
}

const runs = process.argv.slice(2).map(load);
const out = {runs: [], conditions: runs[0].conditions, crossRun: {}};

for (const run of runs) {
    const phases = segment(run.samples).map(ph => summarisePhase(run, ph));
    // Velocity spikes: samples more than 4 robust SDs above the run's median velocity.
    const v = run.samples.map(s => s.velocity).filter(x => x > 0);
    const sorted = [...v].sort((a, b) => a - b);
    const med = sorted[Math.floor(sorted.length / 2)] || 0;
    const mad = [...v.map(x => Math.abs(x - med))].sort((a, b) => a - b)[Math.floor(v.length / 2)] || 0;
    const spikes = run.samples.filter(s => s.velocity > med + 4 * 1.4826 * mad).map(s => s.tick);
    const lastFounderTick = (run.samples.find(s => s.foundersAlive === 1) || {}).tick ?? null;
    out.runs.push({
        name: run.name, seed: run.seed, ticks: run.ticks,
        phases, spikes, lastFounderTick,
        phaseSequence: phases.map(p => p.label).join(' → '),
        series: {
            tick: run.samples.map(s => s.tick),
            pop: run.samples.map(s => +s.pop.toFixed(4)),
            colors: [0, 1, 2].map(c => run.samples.map(s => +s.colors[c].toFixed(4))),
            meanGenome: run.samples.map(s => +s.meanGenome.toFixed(2)),
            distinctGenomes: run.samples.map(s => s.distinctGenomes),
            effGenomes: run.samples.map(s => +s.effGenomes.toFixed(2)),
            domGenomeShare: run.samples.map(s => +s.domGenomeShare.toFixed(4)),
            foundersAlive: run.samples.map(s => s.foundersAlive),
            velocity: run.samples.map(s => +(s.velocity || 0).toFixed(5)),
            ruleDeaths: run.samples.map(s => +(s.ruleDeaths || 0).toFixed(2)),
        },
        // Dominant genome per sample, for the genome-strip view.
        domGenomes: run.samples.map(s => s.consensus),
        consensusAgreement: run.samples.map(s => +s.consensusAgreement.toFixed(3)),
        meanDistFromConsensus: run.samples.map(s => +s.meanDistFromConsensus.toFixed(1)),
        consensusSize: run.samples.map(s => s.consensusSize),
    });
}

// Compare the settled (final) phase dominant genomes across runs, grouped by phase label.
const finals = out.runs.map(r => ({run: r.name, ...r.phases[r.phases.length - 1]}));
const pairs = [];
for (let i = 0; i < finals.length; i++)
    for (let j = i + 1; j < finals.length; j++) {
        const a = finals[i], b = finals[j];
        if (!a.domGenome || !b.domGenome) continue;
        // Restrict to conditions either genome has a rule for: the "used" part of rule space.
        let union = 0, diff = 0;
        for (let c = 0; c < a.domGenome.length; c++) {
            if (a.domGenome[c] === '.' && b.domGenome[c] === '.') continue;
            union++;
            if (a.domGenome[c] !== b.domGenome[c]) diff++;
        }
        // And to conditions cells actually meet (>0.1% of cell-ticks) in either run: rules that matter.
        let active = 0, diffActive = 0;
        for (let c = 0; c < a.domGenome.length; c++) {
            if (!((a.enc[c] || 0) > 0.001 || (b.enc[c] || 0) > 0.001)) continue;
            active++;
            if (a.domGenome[c] !== b.domGenome[c]) diffActive++;
        }
        pairs.push({a: a.run, b: b.run, labelA: a.label, labelB: b.label, sameLabel: a.label === b.label,
            hamming: hamming(a.domGenome, b.domGenome), unionRules: union, diffInUnion: diff,
            activeConds: active, diffActive});
    }
out.crossRun.finalPairs = pairs;

// Rules shared by every final dominant genome with the same phase label: the candidate signature.
const byLabel = {};
for (const f of finals) (byLabel[f.label] = byLabel[f.label] || []).push(f);
out.crossRun.signatures = Object.entries(byLabel).map(([lab, fs_]) => {
    const g = fs_.map(f => f.domGenome).filter(Boolean);
    const shared = [];
    if (g.length) for (let c = 0; c < g[0].length; c++)
        if (g[0][c] !== '.' && g.every(x => x[c] === g[0][c])) shared.push(c);
    return {label: lab, runs: fs_.map(f => f.run), sizes: g.map(x => x.replace(/\./g, '').length), sharedRules: shared};
});

process.stdout.write(JSON.stringify(out));
