// ─── Lineage + genome analysis ────────────────────────────────────────────────
// Pure functions over a parsed lineage-genomes file (OrganismLog.parse). Used by
// lineage-explorer.html and by Node scripts in analysis/.
//
// Organism ids are assigned in birth order, so a parent always has a smaller id than
// its children and birth[] is non-decreasing. Several functions rely on that.

const LineageAnalysis = (() => {

    // Highest organism id born at or before tick t (binary search over birth[]).
    function lastIdBornBy(d, t) {
        let lo = 0, hi = d.header.organisms - 1, ans = -1;
        while (lo <= hi) {
            const mid = (lo + hi) >> 1;
            if (d.birth[mid] <= t) { ans = mid; lo = mid + 1; } else hi = mid - 1;
        }
        return ans;
    }

    const isAlive = (d, id, t) => d.birth[id] <= t && (d.death[id] === -1 || d.death[id] > t);

    // For one tick t: cnt[v] = organisms alive at t descended from v (v included),
    // and maxChild[v] = the largest such count among v's children.
    // Children have larger ids than parents, so one descending pass propagates counts.
    function descendantCounts(d, t, cnt, maxChild) {
        const top = lastIdBornBy(d, t);
        cnt.fill(0, 0, top + 1);
        if (maxChild) maxChild.fill(0, 0, top + 1);
        let living = 0;
        for (let id = top; id >= 0; id--) {
            const alive = d.death[id] === -1 || d.death[id] > t ? 1 : 0;
            living += alive;
            const c = cnt[id] + alive;
            cnt[id] = c;
            if (c === 0) continue;
            const p = d.parent[id];
            if (p >= 0) {
                cnt[p] += c;
                if (maxChild && c > maxChild[p]) maxChild[p] = c;
            }
        }
        return {top, living};
    }

    // Sample ticks: the stats series ticks when present, else evenly spaced.
    function sampleTicks(d, count) {
        const end = d.header.exportedAtTick;
        const s = d.header.series;
        if (s && s.population && s.population.length >= 10 && !count)
            return s.population.map((_, i) => (i + 1) * s.reportingPeriod).filter(t => t <= end);
        const n = count || 200;
        return Array.from({length: n}, (_, i) => Math.round((i + 1) * end / n));
    }

    // Major clades and the time to the most recent common ancestor (TMRCA).
    // A clade is the set of living descendants of one organism. It is "major" if at some
    // sample it holds >= threshold of the living population and is a branch point there
    // (no single child carries all of its living descendants).
    function analyzeClades(d, times, opts = {}) {
        const threshold = opts.threshold ?? 0.05;
        const N = d.header.organisms;
        const cnt = new Int32Array(N), maxChild = new Int32Array(N);
        const living = [], tmrca = [], mrca = [];
        const major = new Set();

        // Pass 1: find major clades and the MRCA at each sample.
        for (const t of times) {
            const {top, living: L} = descendantCounts(d, t, cnt, maxChild);
            living.push(L);
            // MRCA = youngest (highest id) organism that is an ancestor of everyone alive.
            let m = -1;
            for (let id = top; id >= 0; id--) if (cnt[id] === L) { m = id; break; }
            mrca.push(m);
            tmrca.push(m >= 0 && L > 0 ? t - d.birth[m] : null);
            if (L === 0) continue;
            const need = Math.max(2, Math.ceil(threshold * L));
            for (let id = 0; id <= top; id++)
                if (cnt[id] >= need && maxChild[id] < cnt[id]) major.add(id);
        }

        // Pass 2: each major clade's count at every sample.
        const ids = [...major].sort((a, b) => a - b);
        const counts = ids.map(() => new Int32Array(times.length));
        times.forEach((t, k) => {
            descendantCounts(d, t, cnt, null);
            ids.forEach((id, j) => { counts[j][k] = cnt[id]; });
        });

        // Nest clades: parent clade = nearest ancestor that is itself a major clade.
        const index = new Map(ids.map((id, j) => [id, j]));
        const clades = ids.map((id, j) => {
            let p = d.parent[id];
            while (p >= 0 && !index.has(p)) p = d.parent[p];
            const freq = Array.from(counts[j], (c, k) => living[k] ? c / living[k] : 0);
            let peak = 0, peakK = 0;
            freq.forEach((f, k) => { if (f > peak) { peak = f; peakK = k; } });
            // Extinct: had living descendants at some sample and none at the last one.
            const lastK = times.length - 1;
            let extinctTick = null;
            if (counts[j][lastK] === 0 && peak > 0) {
                for (let k = peakK; k <= lastK; k++) if (counts[j][k] === 0) { extinctTick = times[k]; break; }
            }
            return {id, parentClade: p >= 0 ? index.get(p) : -1, founder: d.founder[id], birth: d.birth[id],
                count: counts[j], freq, peak, peakTick: times[peakK], extinctTick};
        });
        clades.forEach((c, j) => { c.children = []; });
        clades.forEach((c, j) => { if (c.parentClade >= 0) clades[c.parentClade].children.push(j); });
        return {times, living, tmrca, mrca, clades, threshold};
    }

    // Lay out a Muller plot: for each clade and sample, [bottom, top] in 0..1 of the living
    // population. Children are centred inside their parent's band (standard Muller layout).
    function mullerLayout(an) {
        const {times, clades} = an;
        const T = times.length;
        const bands = clades.map(() => ({lo: new Float32Array(T), hi: new Float32Array(T)}));
        const roots = clades.map((c, j) => j).filter(j => clades[j].parentClade < 0)
            .sort((a, b) => clades[a].birth - clades[b].birth || a - b);
        for (let k = 0; k < T; k++) {
            // Top-level clades stack from 0 (the rest, minor lineages, sits above them);
            // centring them would shift every band whenever their total changes.
            const place = (list, start, span, centre = true) => {
                const total = list.reduce((s, j) => s + clades[j].freq[k], 0);
                let y = start + (centre ? Math.max(0, span - total) / 2 : 0);
                for (const j of list) {
                    const f = clades[j].freq[k];
                    bands[j].lo[k] = y; bands[j].hi[k] = y + f;
                    const kids = clades[j].children.slice().sort((a, b) => clades[a].birth - clades[b].birth || a - b);
                    if (kids.length) place(kids, y, f);
                    y += f;
                }
            };
            place(roots, 0, 1, false);
        }
        return bands;
    }

    // Most recent common ancestor of two organisms (ids shrink toward the root).
    function mrcaOf(d, a, b) {
        while (a !== b) {
            if (a < 0 || b < 0) return -1;
            if (a > b) a = d.parent[a]; else b = d.parent[b];
        }
        return a;
    }

    function livingAt(d, t) {
        const out = [];
        const top = lastIdBornBy(d, t);
        for (let id = 0; id <= top; id++) if (d.death[id] === -1 || d.death[id] > t) out.push(id);
        return out;
    }

    // Deterministic sample of k items (mulberry32 with a fixed seed so views are stable).
    function sample(arr, k, seed = 12345) {
        let s = seed;
        const rnd = () => { s = (s + 0x6D2B79F5) | 0; let x = Math.imul(s ^ (s >>> 15), 1 | s); x = (x + Math.imul(x ^ (x >>> 7), 61 | x)) ^ x; return ((x ^ (x >>> 14)) >>> 0) / 4294967296; };
        const a = arr.slice();
        for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; }
        return a.slice(0, Math.min(k, a.length));
    }

    // Lineage distance vs genome distance for pairs of organisms alive at tick t.
    // Lineage distance = ticks back to their common ancestor. Genome distance = number of
    // conditions where their rules differ, split by whether cells actually encounter that
    // condition (encountered[c] true) at this time.
    function pairDistances(d, t, opts = {}) {
        const k = opts.sampleSize ?? 120;
        const enc = opts.encountered || null;
        const ids = sample(livingAt(d, t), k, opts.seed);
        const ncond = d.header.ncond;
        const genomes = ids.map(id => d.genomeOf(id));
        const pairs = [];
        for (let i = 0; i < ids.length; i++)
            for (let j = i + 1; j < ids.length; j++) {
                const m = mrcaOf(d, ids[i], ids[j]);
                let all = 0, encDiff = 0;
                const ga = genomes[i], gb = genomes[j];
                for (let c = 0; c < ncond; c++) if (ga[c] !== gb[c]) { all++; if (enc && enc[c]) encDiff++; }
                pairs.push({a: ids[i], b: ids[j], coalescence: m >= 0 ? t - d.birth[m] : null, genomeDist: all, encounteredDist: encDiff});
            }
        return {ids, pairs};
    }

    // Tree of a sample of organisms alive at tick t: the union of their ancestor paths,
    // with chains of single-child ancestors collapsed. Returns nodes with x = birth tick
    // and y = leaf order, ready to draw as a time-scaled tree.
    function livingTree(d, t, opts = {}) {
        const leaves = sample(livingAt(d, t), opts.sampleSize ?? 300, opts.seed);
        const childCount = new Map();
        const inTree = new Set();
        for (const leaf of leaves) {
            for (let a = leaf; a >= 0; a = d.parent[a]) {
                const fresh = !inTree.has(a);
                inTree.add(a);
                const p = d.parent[a];
                if (p >= 0 && fresh) childCount.set(p, (childCount.get(p) || 0) + 1);
                if (!fresh) break;
            }
        }
        const leafSet = new Set(leaves);
        // Keep leaves, branch points and founders; skip single-child links.
        const keep = id => leafSet.has(id) || (childCount.get(id) || 0) > 1 || d.parent[id] < 0;
        const nodes = new Map();
        for (const id of inTree) if (keep(id)) {
            let p = d.parent[id];
            while (p >= 0 && !keep(p)) p = d.parent[p];
            nodes.set(id, {id, parent: p, birth: d.birth[id], children: [], leaf: leafSet.has(id)});
        }
        const roots = [];
        for (const n of nodes.values()) (n.parent >= 0 ? nodes.get(n.parent).children : roots).push(n.id);
        // Order leaves depth-first, children by birth, so related leaves sit together.
        let y = 0;
        const visit = id => {
            const n = nodes.get(id);
            n.children.sort((a, b) => nodes.get(a).birth - nodes.get(b).birth || a - b);
            if (!n.children.length) { n.y = y++; return n.y; }
            const ys = n.children.map(visit);
            n.y = (Math.min(...ys) + Math.max(...ys)) / 2;
            return n.y;
        };
        roots.sort((a, b) => a - b).forEach(visit);
        return {t, leaves, nodes, roots, leafCount: y};
    }

    // Rule changes along the ancestry from organism `from` (exclusive) down to `to`
    // (inclusive). `from` must be an ancestor of `to`, or -1 for "from an empty genome".
    function mutationsBetween(d, from, to) {
        const chain = [];
        for (let a = to; a !== from && a >= 0; a = d.parent[a]) chain.push(a);
        const before = from >= 0 ? d.genomeOf(from) : new Int8Array(d.header.ncond).fill(-1);
        const g = before.slice();
        const changed = new Map();
        for (let k = chain.length - 1; k >= 0; k--) {
            const a = chain[k];
            for (let x = d.deltaStart[a]; x < d.deltaStart[a + 1]; x++) {
                const c = d.deltaCond[x];
                g[c] = d.deltaVal[x];
                changed.set(c, true);
            }
        }
        const out = [];
        for (const c of changed.keys()) if (g[c] !== before[c]) out.push({cond: c, from: before[c], to: g[c]});
        return out.sort((a, b) => a.cond - b.cond);
    }

    // Readable rule, e.g. "red with 3 red, 3 blue → blue". Assumes no history rules.
    function describeRule(conditions, cond, out, names) {
        const c = conditions[cond];
        const offset = c.length - names.length;             // 1 normally, 2 with history rules
        const parts = [];
        for (let k = 0; k < names.length; k++) if (c[offset + k]) parts.push(c[offset + k] + ' ' + names[k]);
        const dead = 6 - c.slice(offset).reduce((a, b) => a + b, 0);
        if (dead) parts.push(dead + ' empty');
        return names[c[0]] + ' with ' + parts.join(', ') + ' → ' + (out < 0 ? 'no rule (dies)' : names[out]);
    }

    return {lastIdBornBy, isAlive, descendantCounts, sampleTicks, analyzeClades, mullerLayout, mrcaOf,
        livingAt, sample, pairDistances, livingTree, mutationsBetween, describeRule};
})();

if (typeof module !== 'undefined') module.exports = LineageAnalysis;
