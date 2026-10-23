// Are final-phase consensus genomes the same rulebook up to relabeling the colours?
// The model treats colours symmetrically, so a red-blue world and a red-green world can
// run identical rules with green and blue swapped. For each pair of runs, try all 6 colour
// permutations and report the best agreement over the rules either genome contains.
// usage: node symmetry.js rd.json > symmetry.json

const fs = require('fs');
const rd = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const conds = rd.conditions;                 // [self, nR, nG, nB]
const key = c => c.join(',');
const index = new Map(conds.map((c, i) => [key(c), i]));
const SYM = 'RGB';
const perms = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];

// Apply colour map p (old colour k -> p[k]) to a genome string.
function permute(g, p) {
    const out = new Array(g.length).fill('.');
    for (let i = 0; i < conds.length; i++) {
        if (g[i] === '.') continue;
        const [s, ...n] = conds[i];
        const m = [0, 0, 0];
        for (let k = 0; k < 3; k++) m[p[k]] = n[k];
        const j = index.get(key([p[s], ...m]));
        out[j] = SYM[p[SYM.indexOf(g[i])]];
    }
    return out.join('');
}
function compare(a, b) {
    let union = 0, same = 0;
    for (let i = 0; i < a.length; i++) {
        if (a[i] === '.' && b[i] === '.') continue;
        union++;
        if (a[i] === b[i]) same++;
    }
    return {union, same};
}
// Readable form of a rule, e.g. "R with 3R 3B → B".
function describe(i, out) {
    const [s, r, g, b] = conds[i];
    const nb = [[r, 'R'], [g, 'G'], [b, 'B']].filter(x => x[0]).map(x => x[0] + x[1]);
    const dead = 6 - r - g - b;
    if (dead) nb.push(dead + ' empty');
    return SYM[s] + ' with ' + nb.join(' ') + ' → ' + out;
}

const finals = rd.runs.map(r => ({name: r.name, seed: r.seed, label: r.phases[r.phases.length - 1].label, g: r.domGenomes[r.domGenomes.length - 1]}));
const pairs = [];
for (let i = 0; i < finals.length; i++)
    for (let j = i + 1; j < finals.length; j++) {
        const a = finals[i], b = finals[j];
        const raw = compare(a.g, b.g);
        let best = null;
        for (const p of perms) {
            const c = compare(permute(a.g, p), b.g);
            if (!best || c.same / c.union > best.same / best.union) best = {...c, perm: p.map(k => SYM[k]).join('')};
        }
        pairs.push({a: a.seed, b: b.seed, labelA: a.label, labelB: b.label, raw, best});
    }

// Canonical form: relabel every final genome onto red-blue (the two colours present → R, B)
// and count how many runs carry each rule.
function toRB(f) {
    const present = [...f.label.replace(' full', '')].map(ch => SYM.indexOf(ch));
    const absent = [0, 1, 2].find(k => !present.includes(k));
    const p = [0, 0, 0];
    p[present[0]] = 0; p[present[1]] = 2; p[absent] = 1;   // first present → R, second → B, absent → G
    return permute(f.g, p);
}
const canon = finals.map(f => ({seed: f.seed, label: f.label, g: toRB(f)}));
const ruleCounts = new Map();
for (const c of canon) for (let i = 0; i < c.g.length; i++) if (c.g[i] !== '.') {
    const k = i + c.g[i];
    ruleCounts.set(k, (ruleCounts.get(k) || 0) + 1);
}
const rules = [...ruleCounts.entries()].map(([k, n]) => ({cond: +k.slice(0, -1), out: k.slice(-1), runs: n}))
    .sort((x, y) => y.runs - x.runs || x.cond - y.cond)
    .map(r => ({...r, text: describe(r.cond, r.out)}));

process.stdout.write(JSON.stringify({pairs, canonical: canon, rules}, null, 1));
