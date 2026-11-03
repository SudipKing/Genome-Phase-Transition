// Colour relabeling shared by the cross-run analyses.
// A relabeling p maps old colour k to p[k]. It moves a condition's own colour, its
// neighbour counts and the rule's output together, so a rule maps to the matching rule of
// the relabeled world. Genomes here are arrays of output indices: 0 = no rule, 1.. = colour.

const SYM = 'RGBYCM';

function permutations(k) {
    if (k === 1) return [[0]];
    const out = [];
    for (const p of permutations(k - 1))
        for (let i = 0; i <= p.length; i++) out.push([...p.slice(0, i), k - 1, ...p.slice(i)]);
    return out;
}

// conds: run.json conditions, [self, count colour 0, ..., count colour n-1] (no history rules).
function colourSpace(conds, n) {
    const ncond = conds.length, nOut = n + 1;
    const condIndex = new Map(conds.map((c, i) => [c.join(','), i]));
    const perms = permutations(n).sort((a, b) => a.join('').localeCompare(b.join('')));
    const permName = pi => perms[pi].map(k => SYM[k]).join('');   // "BGR": R→B, G→G, B→R
    const condMaps = perms.map(p => conds.map(([s, ...cnt]) => {
        const m = new Array(n).fill(0);
        for (let k = 0; k < n; k++) m[p[k]] = cnt[k];
        return condIndex.get([p[s], ...m].join(','));
    }));

    // Genome string ('.RGB…') → Uint8Array of output indices.
    function encode(g) {
        const a = new Uint8Array(ncond);
        for (let c = 0; c < ncond; c++) a[c] = g[c] === '.' ? 0 : SYM.indexOf(g[c]) + 1;
        return a;
    }
    function decode(a) {
        let s = '';
        for (let c = 0; c < ncond; c++) s += a[c] ? SYM[a[c] - 1] : '.';
        return s;
    }
    function relabelGenome(a, pi) {
        const p = perms[pi], cm = condMaps[pi], b = new Uint8Array(ncond);
        for (let c = 0; c < ncond; c++) b[cm[c]] = a[c] ? p[a[c] - 1] + 1 : 0;
        return b;
    }
    // Profile: Float64Array(ncond × nOut) of output shares per condition.
    function relabelProfile(P, pi) {
        const p = perms[pi], cm = condMaps[pi], Q = new Float64Array(P.length);
        for (let c = 0; c < ncond; c++) {
            const j = cm[c];
            Q[j * nOut] = P[c * nOut];
            for (let k = 0; k < n; k++) Q[j * nOut + 1 + p[k]] = P[c * nOut + 1 + k];
        }
        return Q;
    }
    // Relabeling that makes the colours in `order` (most common first) R, G, B, …;
    // absent colours fill the remaining labels in their original order.
    function canonicalPerm(order) {
        const p = new Array(n), used = new Set();
        order.forEach((k, rank) => { p[k] = rank; used.add(k); });
        let next = order.length;
        for (let k = 0; k < n; k++) if (!used.has(k)) p[k] = next++;
        return perms.findIndex(q => q.every((v, i) => v === p[i]));
    }
    // Readable rule, e.g. "R with 2R 1G 3 empty → G".
    function ruleText(c, out) {
        const [s, ...cnt] = conds[c];
        const nb = cnt.map((v, k) => v ? v + SYM[k] : '').filter(Boolean);
        const empty = 6 - cnt.reduce((x, y) => x + y, 0);
        if (empty) nb.push(empty + ' empty');
        return `${SYM[s]} with ${nb.join(' ')} → ${out}`;
    }
    return {SYM, ncond, nOut, n, perms, permName, encode, decode, relabelGenome, relabelProfile, canonicalPerm, ruleText};
}

// Endpoint "box": colours above `present` of living cells, most common first.
function endpointBox(sample, present) {
    if (sample.pop < 0.01) return {label: 'extinct', order: []};
    const order = sample.colors.map((f, k) => [f, k]).filter(x => x[0] > present).sort((a, b) => b[0] - a[0]).map(x => x[1]);
    return {label: order.map(k => SYM[k]).join('>') + (sample.pop >= 0.95 ? ' full' : ''), order};
}

module.exports = {SYM, colourSpace, endpointBox};
