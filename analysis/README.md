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
- `genome-snapshots.jsonl` stores the exact genome population distribution at every reporting sample. Each genome is stored once per sample with its population count; it does not store one full genome per organism.
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
