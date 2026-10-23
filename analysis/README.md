# Analysis tools

Headless batch runs of the simulation and the phase/genome analysis behind the progress report.

## Run simulations

`run.js` loads the simulation's own source files (`hexca.js`, `stats.js`, `rulespace.js`, ...) into Node with a seeded random number generator, so a seed reproduces a run exactly. Every reporting period (1,000 ticks) it records the existing stats plus distinct-genome counts, the consensus genome (majority output per condition), the rule-frequency matrix, condition encounter rates and surviving founder lineages.

```
node analysis/run.js --seed 1 --ticks 400000 --out analysis/runs/r1 --historyTicks 0
```

Options:

- `--historyTicks N` keeps genealogy records only for the first N ticks (they grow by ~150 per tick).
- `--genealogy` writes the same family-tree JSON as the in-browser export.
- `--settle N` stops once the grid has stayed full for N samples (default 40).
- `--params '{"mutationRate":0.02}'` overrides parameters.

The runner saves `samples.json` every 10 samples.

## Analyze

```
node analysis/analyze.js analysis/runs/r1 analysis/runs/r2 ... > analysis/rd.json
node analysis/symmetry.js analysis/rd.json > analysis/symmetry.json
node analysis/build-report.js analysis/rd.json analysis/findings.json analysis/report.html
```

- `analyze.js` splits each run into phases (which colours are above 8% of living cells, and whether the grid is at least 95% full) and summarizes each phase.
- `symmetry.js` compares final-phase consensus genomes across runs under all six colour relabelings.
- `build-report.js` renders `template.html` with the data and the text in `findings.json`.
