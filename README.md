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
- Individual IDs are unique within one simulation run and reset when the simulation is reset. Parent IDs can later be used to reconstruct ancestry; a full family tree visualization is a later step.
