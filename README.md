# Genome Phase Transition Research

Combined from Chris Marriott's current Hex CA implementation with research focused founder lineage tracking.

## Views
- Population: original cell state colors and population statistics.
- Genotype (G): colors by a projection of the current genome. Mutations can change the color; it is not an ancestry identifier or evidence of similar behavior.
- Founder lineage (L): every initially living cell receives a unique founder ID. Descendants inherit it even when their genomes mutate.

Founder IDs reset with each simulation reset and should not be compared as IDs across independent runs.

## Research direction
The current priority is tracking which founder lineages persist, expand, or disappear. Later, repeated runs can be compared and lineage survival related to observed phase transitions and genome characteristics. Colors support visual inspection but are not quantitative evidence by themselves.

Keys 1–4 switch statistics pages, G toggles genotype coloring, L toggles founder-lineage coloring, O cycles raster ordering, M marks a reference, and E exports rule-space JSON.


## Individual and parent tracking
- Each initially living organism receives a unique individual ID and has no parent (parent ID `-1`).
- Each offspring receives a new individual ID and records the individual ID of the organism that reproduced it as its parent ID.
- Founder lineage ID remains separate and is inherited from the original founder, even when the offspring's genome mutates.
- Individual IDs are unique within one simulation run and reset when the simulation is reset. Parent IDs are used to reconstruct ancestry and compute lineage statistics. Export the genealogy JSON from the simulation, then open `genealogy.html` in a browser and load that JSON to view the family tree and download a lineage summary CSV.


## Genealogy analysis
1. Run the simulation and click **Export family tree JSON**.
2. Open [genealogy.html](./genealogy.html) in a browser.
3. Select the exported JSON file.
4. Inspect the run overview, compare founder lineages, and select a founder to view its family tree.
5. Use **Download lineage summary CSV** to save the per-founder statistics for later analysis.

The genealogy explorer processes the selected file locally in the browser. It reports recorded descendants, living and dead members at export, maximum generation depth, birth tick range, and whether each founder lineage is extinct at export. These are descriptive ancestry measures, not proof that a lineage caused a phase transition.

## Lineage + genome analysis
1. Run the simulation and click **Export lineage + genomes**. This saves a compact `.hxlg` file with every organism's parent, founder, birth and death tick, and genome (stored as the rules that changed from its parent), plus the stats time series.
2. Open [lineage-explorer.html](./lineage-explorer.html) in a browser and select the file.
3. Use the timeline (population by colour, a Muller plot of clades coloured by genome, the time to the most recent common ancestor, births and deaths) to line up phase changes with extinctions and sweeps. Click to inspect a clade's genome and the rule changes it carries, or to set the tick for the zoomable family tree and the lineage-vs-genome distance plot.

A 45,000-tick run (about 5 million organisms) exports to about 160 MB and analyses in about 5 seconds. See [analysis/README.md](./analysis/README.md) for the file format and the batch analysis scripts.
