// Draft paper figure: do runs with the same endpoint share genomes?
// Reads the JSON from endpoints.js and clusters.js and writes one self-contained HTML page:
// (A) each run's final colour mix beside the genetic clusters its organisms fall into,
// (B) every pair's similarity raw vs relabeled, (C) mean genome distance between runs after
// relabeling, (D) the core rules nearly every run shares, (E) the rules that define each cluster.
//
// usage: node analysis/endpoint-figure.js endpoints.json clusters.json figure.html

const fs = require('fs');
const [epFile, clFile, outFile] = process.argv.slice(2);
if (!outFile) { console.error('usage: node endpoint-figure.js endpoints.json clusters.json figure.html'); process.exit(1); }
const ep = JSON.parse(fs.readFileSync(epFile, 'utf8'));
const cl = JSON.parse(fs.readFileSync(clFile, 'utf8'));

const esc = s => String(s).replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c]));
const CELL = ['var(--red)', 'var(--green)', 'var(--blue)'];
const CL = ['var(--c1)', 'var(--c2)', 'var(--c3)', 'var(--c4)', 'var(--c5)', 'var(--c6)', 'var(--c7)', 'var(--c8)'];
const clIndex = new Map(cl.clusters.map((c, i) => [c.cluster, i]));
const epRun = new Map(ep.runs.map(r => [r.name, r]));

// Runs ordered by main cluster, then by how strongly they belong to it.
const runs = cl.runs.slice().sort((a, b) => clIndex.get(a.mainCluster) - clIndex.get(b.mainCluster) || b.mainShare - a.mainShare);

// (A) Rows of two stacked bars.
const W = 840, rowH = 26, top = 34, labelW = 150, barW = 250, gap = 40;
let a = `<svg viewBox="0 0 ${W} ${top + runs.length * rowH + 8}" role="img" aria-label="Endpoint colours and genetic clusters per run">`;
a += `<text x="${labelW}" y="20" class="hd">Final colour mix</text><text x="${labelW + barW + gap}" y="20" class="hd">Genetic cluster of sampled organisms</text>`;
runs.forEach((r, i) => {
    const y = top + i * rowH, e = epRun.get(r.name);
    a += `<text x="0" y="${y + 15}" class="lbl">${esc(r.name)} <tspan class="sub">${esc(r.box)}</tspan></text>`;
    let x = labelW;
    (e ? e.colors : []).forEach((f, k) => {
        const w = f * barW;
        if (w > 0.3) a += `<rect x="${x.toFixed(1)}" y="${y + 3}" width="${w.toFixed(1)}" height="16" fill="${CELL[k]}"><title>${'RGB'[k]} ${(f * 100).toFixed(1)}%</title></rect>`;
        x += w;
    });
    x = labelW + barW + gap;
    cl.clusters.forEach((c, ci) => {
        const f = c.runs[r.name] || 0, w = f * barW;
        if (w > 0.3) a += `<rect x="${x.toFixed(1)}" y="${y + 3}" width="${w.toFixed(1)}" height="16" fill="${CL[ci % CL.length]}"><title>Cluster ${c.cluster}: ${(f * 100).toFixed(0)}%</title></rect>`;
        if (w > 22) a += `<text x="${(x + w / 2).toFixed(1)}" y="${y + 15}" class="in">${c.cluster}</text>`;
        x += w;
    });
    const stop = r.stopReason === 'births-and-deaths-below-threshold' ? 'settled' : r.stopReason;
    a += `<text x="${labelW + 2 * barW + gap + 10}" y="${y + 15}" class="sub">${esc(stop)} @ ${(r.tick / 1000).toFixed(0)}k</text>`;
});
a += '</svg>';

// (B) Heatmap of mean relabeled Hamming distance.
const names = cl.runDistance.names, M = cl.runDistance.relabeled;
const idx = runs.map(r => names.indexOf(r.name));
const vals = M.flat().filter(v => v !== null);
const lo = Math.min(...vals), hi = Math.max(...vals);
const cs = 34, off = 60;
let b = `<svg viewBox="0 0 ${off + idx.length * cs + 10} ${off + idx.length * cs + 10}" role="img" aria-label="Mean genome distance between runs" class="heat">`;
idx.forEach((ri, i) => {
    b += `<text x="${off + i * cs + cs / 2}" y="${off - 8}" class="lbl" text-anchor="middle">${esc(names[ri])}</text>`;
    b += `<text x="${off - 8}" y="${off + i * cs + cs / 2 + 4}" class="lbl" text-anchor="end">${esc(names[ri])}</text>`;
    idx.forEach((rj, j) => {
        const v = M[ri][rj], t = (v - lo) / ((hi - lo) || 1);
        b += `<rect x="${off + j * cs}" y="${off + i * cs}" width="${cs - 2}" height="${cs - 2}" rx="3" fill="var(--heat)" fill-opacity="${(1 - t * 0.88).toFixed(2)}"><title>${esc(names[ri])} vs ${esc(names[rj])}: ${v}</title></rect>`;
        b += `<text x="${off + j * cs + cs / 2 - 1}" y="${off + i * cs + cs / 2 + 3}" class="${t < 0.45 ? 'cellv on' : 'cellv'}" text-anchor="middle">${Math.round(v)}</text>`;
    });
});
b += '</svg>';

// (P) Pair similarity, raw vs relabeled, split by whether the runs kept the same colours.
const colourSet = name => { const e = epRun.get(name); return e ? e.colors.map((f, k) => f > ep.method.present ? 'RGB'[k] : '').join('') : ''; };
const groups = [['Same surviving colours', ep.pairs.filter(q => colourSet(q.a) === colourSet(q.b))],
                ['Different surviving colours', ep.pairs.filter(q => colourSet(q.a) !== colourSet(q.b))]];
const PW = 760, px0 = 190, px1 = PW - 20, rowP = 46, toX = v => px0 + v * (px1 - px0);
let pp = `<svg viewBox="0 0 ${PW} ${30 + groups.length * rowP + 36}" role="img" aria-label="Similarity of every pair of runs, raw and relabeled">`;
for (let t = 0; t <= 10; t += 2) {
    const x = toX(t / 10);
    pp += `<line x1="${x}" y1="18" x2="${x}" y2="${18 + groups.length * rowP}" stroke="var(--line)"/><text x="${x}" y="${34 + groups.length * rowP}" class="sub" text-anchor="middle">${(t / 10).toFixed(1)}</text>`;
}
groups.forEach(([label, list], gi) => {
    const y = 30 + gi * rowP;
    pp += `<text x="0" y="${y + 5}" class="lbl">${label} <tspan class="sub">(${list.length})</tspan></text>`;
    list.forEach((q, j) => {
        const jit = ((j * 37) % 17 - 8) * 0.9;
        pp += `<circle cx="${toX(q.raw).toFixed(1)}" cy="${(y + jit).toFixed(1)}" r="4" fill="none" stroke="var(--mute)"><title>${esc(q.a)}–${esc(q.b)} raw ${q.raw}</title></circle>`;
        pp += `<circle cx="${toX(q.relabeled).toFixed(1)}" cy="${(y + jit).toFixed(1)}" r="4" fill="var(--heat)" fill-opacity="0.75"><title>${esc(q.a)}–${esc(q.b)} relabeled ${q.relabeled}</title></circle>`;
    });
});
pp += `<text x="${px0}" y="${52 + groups.length * rowP}" class="sub">○ raw colours   ● best colour relabeling   (1 = identical rule populations)</text></svg>`;

// Core rulebook shared by nearly every run.
const core = (cl.coreRules || []).slice(0, 14);
const coreHtml = core.length ? '<ul>' + core.map(r => `<li><code>${esc(r.text)}</code> <span class="meta">${r.runs}/${cl.runs.length} runs · ${(r.meanShare * 100).toFixed(0)}% of organisms on average</span></li>`).join('') + '</ul>'
    : '<p class="meta">No rule is carried by most organisms in nearly every run.</p>';
const swaps = ep.runs.filter(r => r.survivorSwap !== null && r.survivorSwap !== undefined).map(r => r.survivorSwap);
const swapText = swaps.length ? `Swapping a run's two surviving colours leaves its rule population ${Math.min(...swaps).toFixed(2)}–${Math.max(...swaps).toFixed(2)} similar to itself.` : '';

// (C) Defining rules per cluster.
let c = '';
cl.clusters.forEach((k, ci) => {
    const from = Object.entries(k.runs).filter(([, v]) => v >= 0.1).map(([n, v]) => `${esc(n)} ${(v * 100).toFixed(0)}%`).join(', ');
    c += `<div class="cl"><h3><span class="sw" style="background:${CL[ci % CL.length]}"></span>Cluster ${k.cluster}<span class="meta">${k.size} genomes · spread ${k.meanDistToMedoid} · ${from}</span></h3>`;
    c += k.definingRules.length ? '<ul>' + k.definingRules.slice(0, 6).map(r => `<li><code>${esc(r.text)}</code> <span class="meta">${(r.inCluster * 100).toFixed(0)}% here vs ${(r.elsewhere * 100).toFixed(0)}% elsewhere</span></li>`).join('') + '</ul>'
        : '<p class="meta">No rule separates this cluster cleanly.</p>';
    c += '</div>';
});

const settled = cl.runs.filter(r => r.stopReason === 'births-and-deaths-below-threshold').length;
const html = `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Genome Endpoints</title><style>
:root{--bg:#fbfaf7;--fg:#1d1d1f;--mute:#6b6b70;--line:#e3e1dc;--red:#d64541;--green:#3c9a4f;--blue:#3a6fd8;--heat:#3b3f8f;
--c1:#7b61c9;--c2:#d18a1f;--c3:#2a9d8f;--c4:#c2457a;--c5:#5a6f86;--c6:#8a9a2b;--c7:#b0603a;--c8:#4f8fbf}
@media (prefers-color-scheme:dark){:root{--bg:#151517;--fg:#ececee;--mute:#9a9aa2;--line:#2c2c31;--heat:#a9b0ff}}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:860px;margin:0 auto;padding:28px 16px 48px}h1{font-size:22px;margin:0 0 4px}h2{font-size:15px;margin:28px 0 8px}
h3{font-size:14px;margin:0 0 6px;display:flex;align-items:center;gap:8px;flex-wrap:wrap}.meta,.sub{color:var(--mute);font-size:12px;font-weight:400}
p.lead{color:var(--mute);margin:0 0 8px}svg{width:100%;height:auto;display:block}.heat{max-width:${off + idx.length * cs + 10}px}
svg text{fill:var(--fg);font-size:12px}svg .hd{font-weight:600}svg .sub{fill:var(--mute);font-size:11px}svg .in{fill:#fff;font-size:11px;font-weight:600;text-anchor:middle}
svg .cellv{font-size:10px;fill:var(--fg)}svg .cellv.on{fill:var(--bg)}.sw{width:12px;height:12px;border-radius:3px;display:inline-block}
.cl{border-top:1px solid var(--line);padding:10px 0}ul{margin:0;padding-left:18px}code{font-size:12.5px}
</style></head><body><main>
<h1>Do runs with the same endpoint share genomes?</h1>
<p class="lead">${cl.runs.length} independent runs (${settled} reached the births-and-deaths stop${settled < cl.runs.length ? ', the rest hit the tick limit' : ''}). ${cl.method.sample} organisms sampled from each final population, colours relabeled (${esc(cl.method.relabel)}), then clustered by Hamming distance over all 252 conditions (best k = ${cl.k}, silhouette ${cl.silhouette.find(s => s.k === cl.k).silhouette}; near 0 means no real split). ${swapText}</p>
<h2>A. Endpoint colours and genetic clusters</h2>${a}
<h2>B. Similarity of every pair of runs, before and after relabeling colours</h2>${pp}
<h2>C. Mean genome distance between runs (relabeled; diagonal = spread within a run)</h2>${b}
<h2>D. Core rulebook (relabeled; R and G are the surviving colours)</h2>${coreHtml}
<h2>E. Rules that define each cluster</h2>${c}
</main></body></html>`;
fs.writeFileSync(outFile, html);
console.error(`wrote ${outFile}`);
