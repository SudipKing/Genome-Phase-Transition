# Analysis tools

Headless batch runs of the simulation and the phase/genome analysis behind the progress report.

## Run simulations

`run.js` loads the simulation's own source files (`hexca.js`, `stats.js`, `rulespace.js`, ...) into Node with a seeded random number generator, so a seed reproduces a run exactly. Every reporting period (1,000 ticks) it records the existing stats plus distinct-genome counts, the consensus genome (majority output per condition), the rule-frequency matrix, condition encounter rates and surviving founder lineages. It also writes the complete living genome distribution to `genome-snapshots.jsonl`: each row is one sample tick and contains each unique genome once with its organism count.

```
node analysis/run.js --seed 1 --ticks 50000 --out analysis/runs/r1 --historyTicks 0
```

Options:

- `--historyTicks N` keeps genealogy records only for the first N ticks (they grow by ~150 per tick).
- `--genealogy` writes the same family-tree JSON as the in-browser export.
- `--lineageGenomes` writes `lineage-genomes.hxlg`, the same file as the simulator's **Export lineage + genomes** button, for `lineage-explorer.html`.
- `--ticks N` sets the hard run limit (default 50,000).
- `--quietThreshold N` sets the births and total deaths threshold (default 10).
- `--stableSamples N` requires both births and total deaths to remain below the threshold for this many consecutive samples (default 5).
- `genome-snapshots.jsonl.gz` stores the exact genome population distribution at every reporting sample. Each genome is stored once per sample with its population count; it does not store one full genome per organism. Each line is its own gzip member, so the file stays append-only and a partial run is still readable (about 3 MB for a 50,000-tick run instead of about 30 MB). `--plainSnapshots` writes uncompressed `genome-snapshots.jsonl` instead; `readSnapshots(dir)` in `storage/local.js` reads either.
- `run.json` records the seed, parameters, stop reason and snapshot filename.
- `--params '{"mutationRate":0.02}'` overrides parameters.

The runner saves `samples.json` every 10 samples and writes genome snapshots incrementally, so partial runs retain their data.

## Analyze

```
node analysis/analyze.js analysis/runs/r1 analysis/runs/r2 ... > analysis/rd.json
node analysis/symmetry.js analysis/rd.json > analysis/symmetry.json
node analysis/build-report.js analysis/rd.json analysis/findings.json analysis/report.html
```

- `analyze.js` splits each run into phases (which colours are above 8% of living cells, and whether the grid is at least 95% full) and summarizes each phase.
- `symmetry.js` compares final-phase consensus genomes across runs under all six colour relabelings.
- `build-report.js` renders `template.html` with the data and the text in `findings.json`.

### Common endpoints

```
node analysis/endpoints.js analysis/runs/e1 analysis/runs/e2 ... > analysis/runs/endpoints.json
```

`endpoints.js` asks whether independent runs end in the same genetic place, using the actual final genome populations from the snapshots rather than one consensus genome. Each population becomes a rule profile (per condition, the share of organisms with no rule or each output colour). Two profiles are compared by total-variation distance per condition, weighted by how much either population uses that condition, so similarity is 1 for identical populations. Every pair is scored raw and under all colour relabelings (the permutation moves the condition's own colour, its neighbour counts and the output together). It also reports each run's endpoint box (colours above 8%, most common first), how symmetric the endpoint is (`survivorSwap`: similarity to itself with its two main colours swapped), exact genomes shared between runs, average-linkage clusters (`--cluster`, default 0.75), when each run came within 0.9 of its final population, and the rules most runs share once every endpoint is relabeled so its dominant colour is R. A readable summary goes to stderr.

### Genome clusters

```
node analysis/clusters.js analysis/runs/c1 analysis/runs/c2 ... > analysis/runs/clusters.json
```

`clusters.js` is the cluster analysis on individual genomes, with no lineage data. From each run's final snapshot (or `--tick T`) it draws `--sample` organisms at random (default 300, weighted by genome counts, seeded with `--seed`), relabels each run so its surviving colours take the first labels (R, G) and extinct colours the rest (`colours.js`), and measures Hamming distance over all 252 conditions. Surviving colours often hold nearly equal shares, so instead of ranking them by share it tries every order of a run's survivors and keeps the one whose rule profile is closest to the other runs' (`--noAlign` ranks by share instead). k-medoids clusters the pooled genomes for k = 2..`--kmax` and keeps the k with the best silhouette. It reports each run's main cluster, each cluster's medoid genome and the runs it draws from, the rules that define a cluster (carried by at least 80% of it and at most 30% of the rest), and the mean distance between runs (raw and relabeled; the diagonal is the spread within a run). It also lists the core rules: those carried by at least half the organisms in all but at most one run.

`colours.js` holds the colour relabeling shared by `endpoints.js` and `clusters.js`.

### Figure

```
node analysis/endpoint-figure.js analysis/runs/endpoints.json analysis/runs/clusters.json analysis/runs/figure.html
```

Writes one self-contained HTML page from the two analyses: each run's final colours beside its genetic clusters, every pair's similarity raw vs relabeled (split by whether the runs kept the same surviving colours), the run-by-run distance heatmap, the core rulebook and each cluster's defining rules.
`endpoints.js` asks whether independent runs end in the same genetic place, using the actual final genome populations from the snapshots rather than one consensus genome. Each population becomes a rule profile (per condition, the share of organisms with no rule or each output colour). Two profiles are compared by total-variation distance per condition, weighted by how much either population uses that condition, so similarity is 1 for identical populations. Every pair is scored raw and under all colour relabelings (the permutation moves the condition's own colour, its neighbour counts and the output together). It also reports each run's endpoint box (colours above 8%, most common first), exact genomes shared between runs, average-linkage clusters (`--cluster`, default 0.75), when each run came within 0.9 of its final population, and the rules most runs share once every endpoint is relabeled so its dominant colour is R. A readable summary goes to stderr.

## Lineage + genomes

The simulator keeps every organism in `OrganismLog` (`organismlog.js`): parent, founder, birth and death tick, and its genome stored as the rules that differ from its parent. A genome is rebuilt by replaying those changes from the founder down.

Why changes instead of a genome → organism IDs map: `measure-genomes.js` shows that 98% of organisms in a 10k-tick run (1.8M organisms) carry a genome no other organism shares, while a child differs from its parent in only ~4 of 252 rules. Grouping identical genomes would save about 2%; storing changes is about 30× smaller than storing whole genomes.

```
node analysis/measure-genomes.js . 10000 7
node analysis/test-lineage-export.js
node analysis/lineage-report.js analysis/runs/r6/lineage-genomes.hxlg ... > analysis/lineage-report.json
```

- `test-lineage-export.js` runs the simulation in four settings (sync, async, history rules, 4 colours), exports, and checks that every living genome rebuilt from the file matches the simulation exactly and that every row matches the family-tree JSON.
- `lineage-report.js` finds major clades, extinction events and sweeps (sharp drops in the time to the most recent common ancestor), and lines them up with phases. For each sweep it checks how many of the final phase's rules the sweeping lineage's founder already carried, compared with organisms alive at the same time. It also tests whether the rule changes behind successful clades are enriched in conditions cells actually meet, and fits genome distance against lineage distance in each phase.
- `lineage-explorer.html` (repo root) shows the same analysis interactively: a Muller plot of clades coloured by genome, the TMRCA curve, a zoomable tree of the living population, a lineage-vs-genome distance plot and an extinction table.
