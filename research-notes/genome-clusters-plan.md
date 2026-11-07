# Coarse genome clustering — plan (week of Nov 6, 2026)

Source: the two most recent chats in the ChatGPT project "Hex Phase Transition Research - Chris Mariott"
("Genome Rule Analysis" and "Weekly Project Reflection"), plus the current `Genome-Phase-Transition` repo.

## What Chris asked for

The question for this quarter, asked for each of the four phases: **do organisms that reach the same phase
share the same rules, or can different rule sets reach the same visual end state?**

This week's stopping point is a working tool, not a result:

1. **Measure.** Turn each genome into 252 bits (rule present or absent; output colour ignored). The
   difference between two genomes is the number of positions where one has a rule and the other doesn't
   (Hamming distance). Chris chose this "coarse" measure first: it needs no weighting decisions, while the
   "fine" one (does red→green count as 1 or ½?) does.
2. **Cluster.** Two genomes join the same cluster when they are within *t* differences. At *t* = 0 only
   identical rule sets group together. As *t* rises, mutants (a child differs from its parent in about 4
   rules) fold into their parent's cluster. Don't decide ahead of time how many clusters there should be.
3. **Slider.** Sweep *t* and watch clusters merge. If everything becomes one cluster at *t* ≈ 4, the
   population is one genome plus its mutants. If 4 clusters are still separate at *t* = 4, they are
   genuinely different organisms, each with its own mutants.
4. **Report** survivors, unique genomes, number of clusters, and organisms and % of the population per
   cluster.
5. **Show the genome of each cluster** (secondary): which rules it has, how consistently its members
   carry them, and what colours those rules output.

Next week's checkpoint: look at the clusters together and decide whether to keep going in this direction.

Later (this quarter, then next): run this at each phase, compare repeated runs ("is it the same dominant
genome again?"), then trace how one phase's genomes become the next phase's (the story for March).

## How it differs from the existing `clusters.js`

`clusters.js` (Nov 3) samples 300 organisms per run, distinguishes output colours, relabels colours
across runs, and uses k-medoids with k picked by silhouette. That answers a different question, and it
fixes the number of clusters by score. Chris asked for: every survivor, rule presence only, a threshold
the researcher controls, and no assumed cluster count. The new tool sits beside `clusters.js` and leaves
it unchanged.

## Design decisions

| Decision | Choice | Why |
|---|---|---|
| Distance | Hamming over 252 presence bits | Chris's coarse measure; no weights |
| Clustering method | **Connected components** (single linkage): genomes within *t* are linked, and links chain | This is the method ChatGPT recommended documenting. It shows clusters merging as *t* rises, which is what the slider is for. Chaining (A–B = 1, B–C = 1 puts A and C together even at distance 2) is reported through each cluster's spread and diameter, so it can't hide. |
| Scale | Group organisms into unique presence patterns with counts, compare patterns, then weight by count | ~6,000 survivors → ~5,700 patterns → ~16M pairs. Bit-packed popcount takes about a second. |
| Every threshold at once | Build the minimum spanning tree once. Clusters at *t* = MST with edges > *t* cut out. | The slider updates instantly in the browser, and the merge curve (clusters vs *t*) comes for free. |
| Original genomes | Kept: each pattern keeps its exact genomes and counts | Output colours stay available for the cluster viewer and for the fine measure later |
| Checkpoints | `--tick`, `--every`, or `all` | Phase-by-phase analysis needs ticks other than the last one |

## Deliverables

1. `analysis/genome-clusters.js`: Node CLI. Run dir(s) → `genome-clusters.json` + text report.
   - Per checkpoint: survivors, unique exact genomes, unique presence patterns, MST edges, merge curve,
     cluster summaries at the requested thresholds (default 0,1,2,3,4,6,8,12,16).
2. `genome-cluster-explorer.html` (repo root, like `lineage-explorer.html`): open the JSON and get
   - the threshold slider, with a live count of clusters, a size distribution, and the merge curve
   - a cluster table: organisms, % of population, unique genomes, rules in representative, spread, diameter bound
   - a cluster genome viewer: rule-presence frequency for all 252 conditions (3 self-colours × 84
     neighbourhoods), the most common output colour per rule, and the representative genome
   - a checkpoint and run picker
3. `analysis/test-genome-clusters.js`: checks the MST-cut clusters against brute-force connected
   components on random genomes, and checks the distance function.
4. A README section.

## Done when

- `node analysis/test-genome-clusters.js` passes.
- Running it on the real runs (c1–c10, e1–e5) produces a report in a few seconds per run.
- The explorer loads that JSON, and the slider changes the cluster count and distribution.

## Status (Nov 6): built and checked

- `analysis/genome-clusters.js`, `analysis/test-genome-clusters.js` (7 checks pass), `genome-cluster-explorer.html`,
  `analysis/README-genome-clusters.md` (to paste into `analysis/README.md`).
- Ran on the final snapshots of c1, c2, c5, e1 and e3: 3–6 s per 6,000-organism run. Output was written outside the repo; regenerate into `analysis/runs/` (gitignored).

### First look, for discussion with Chris (observations, not conclusions)

- Almost every survivor has its own genome: c1 has 6,000 organisms, 5,710 exact genomes and 5,697 presence
  patterns. Counting identical genomes alone tells us nothing, so the similarity measure is needed.
- 77% of organisms have another genome within 4 rules, so the population is mostly families of mutants.
- **Connected (chaining) clusters show no plateau.** Clusters merge steadily into one giant component
  (c1: 3.5% of the population at t=4, 45% at t=8, 97.5% at t=16), and its members end up as far as 35 rules
  apart. Chaining hides structure here, which is why the explorer also offers the representative method.
- **Representative clusters (no chaining) do show structure** at t=16 in c1: a 19-rule group (44.5%) and a
  27-rule group (8.4%). Both carry the same core of "full-grid" rules (e.g. G with 2G 4B, B with 3G 3B).
  Worth looking at with Chris: same core, plus extra rules?

## Deferred (on purpose)

- Fine distance (output colours, weights): add it only if coarse turns out to be too coarse.
- Matching clusters across repeated runs (needs colour relabeling, as in `colours.js`) is the next
  step after Chris reviews the clusters.
- Tracking clusters across phases over time (next quarter).
