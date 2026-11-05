# Project Log: Genome-Phase-Transition

How this repo came to be: what I set out to build, what I asked the AI agents for, what I decided, and what I dropped along the way. I've kept it in date order. Quotes from me are copied from my prompts, typos and all. Quotes from Chris come from the meeting transcripts I pasted into my ChatGPT project, "Hex Phase Transition Research - Chris Mariott".

- **Course:** TCSS 499 directed research, 2 credits, autumn 2026, with Prof. Chris Marriott (UW Tacoma). Logged for TCSS 490 as Artifact 4.
- **Repo:** [SudipKing/Genome-Phase-Transition](https://github.com/SudipKing/Genome-Phase-Transition). In my prompts I call it "GPT".
- **Built on:** Chris's [algorithm0r/hex-het-CA](https://github.com/algorithm0r/hex-het-CA) and my earlier fork [SudipKing/hex-het-CA-sudip](https://github.com/SudipKing/hex-het-CA-sudip).
- **Agents used:**
  - **ChatGPT, with its GitHub connector, from Sept 28 to Oct 23.** I used it for planning, to read the transcripts of my meetings with Chris, and for the commits on Oct 12 and Oct 17.
  - **Claude Code, from Oct 23 on.** I used it for the batch runs and analysis, the large-export fix, the genome deltas and the lineage explorer.
- **Authorship:** every commit is authored as me, with no co-author lines. That was my explicit instruction (see Oct 23).

---

## 1. Research question

> **Do distinct phases of the heterogeneous hex cellular automaton correspond to characteristic genomes / transition rules?**

**The model.**
- The world is a toroidal hex grid. Each cell has a colour, a genome, energy and an age.
- With 3 colours there are 252 conditions: a cell's own colour plus how many neighbours it has of each colour, counting dead ones.
- A genome maps each condition either to an output colour or to nothing. If a cell meets a condition it has no rule for, it dies.
- Cells earn energy by changing colour, and reproduce with mutation when they have enough.
- Larger genomes cost more to reproduce, so selection favours small ones.
- When Chris ran the simulation, it went through visible phase transitions. In the last one, cells stopped dying ("it's, like, immortal") and genome size fell sharply.

**Chris's hypothesis,** from his Sept 29 email: "there are certain kinds of rules (which correspond to genomes) that characterize the phases. The easiest way for us to find this is probably to do a bunch of runs and see what genomes appear in each phase and look for a pattern."

**What he wants out of it.** He wants to explain the transitions, not just show that they happen: "it would be cool for us just to go to a conference and say, say, Look, this happened. But it's way cooler for us to go and say, And this is why."

---

## 2. Brainstorm and getting started (Sept 28 – Oct 1)

### Sept 28: first meeting with Chris
Before the meeting I asked ChatGPT how my interests lined up with Chris's:

> "I have a meeting today with one of my professors Chris Marriott … I want to do some research opportunities … how can I talk to him so I can do research with him"

**The projects he showed me.** In the meeting Chris walked me through his "simulation zoo":
- hex-pipe organisms that feed each other through pipes
- the heterogeneous hex CA
- Hex-Turmites
- a traffic-safety simulation
- Plant World
- foraging agents
- game trails
- "emergent boundaries" animats

**Why the hex CA.** I chose it because it had a clear, open research question. Chris described it this way:

> "what is it about these rules, these simple rules that I put in there, that allows it to have these different phases that it goes through to get to the final state."

**What the work would be.** Chris framed it as "experimental architecture": software around the simulation that runs it under many settings, gathers data, and analyses it. In his words: "If we can get a very pretty graph that's like phase transition here, phase transition here … then we have a paper."

**Timing.** He suggested a light 2-credit autumn quarter, possibly continuing in winter and spring. I agreed, because I'm also taking two 5-credit classes, one of which is his TCSS 490.

### Sept 29 – Oct 1: proposal
Chris agreed by email and sent the repo. He suggested I use Claude to understand "the way the genes select the rules and how the rules work", then write a one-page proposal and a 10-week plan with one goal per week.

**First draft.** I didn't have Claude Code yet, so I asked ChatGPT: "So first review the repo, i dont have claude but i have u bro. then lets do the proposal".

**Decisions while revising:**
- "put 8 in simple terms, also not doing a presentation, just a research paper is the final". The deliverable is a research paper only, with no presentation.
- "make this sound simpler, way too advanced for someone new to research yknow". I rewrote the project description in plainer language.

---

## 3. Timeline

### Oct 9: understanding the code, second meeting
**Repository report.** I had a detailed report written on how the code works: genomes, rules, reproduction, mutation, stats, and the defects that limit what we can conclude from runs. That is `repository-report.txt` ("Prepared: October 9, 2026"). I committed it to this repo later, on Oct 17.

**What Chris said in the meeting** (transcript pasted into the "Summarize Genome Evolution" chat):
- **Genomes, not just genome size.** The genome-size graph already shows phases, but "it's going to be helpful for us to look at what the genome is … Oh, here's why it got smaller. It got rid of all of these rules."
- **Count and rank genomes.** Track how many distinct genomes are present, rank them by frequency, and find the most popular genome and its rules.
- **Priority one is genome colouring.** Chris showed Hex-Turmites, which gives every distinct genome its own colour: "maybe priority one should be to go to this repo, the Hex Termites repo, find the code that assigns colors … I think it might be the only thing we need."
- **A 252-pixel strip per genome.** He suggested drawing each genome as a line of 252 pixels, coloured by each rule's output, then clustering and sorting the population that way. He made this next week's goal: "this week's goal is still color assigning."
- **Stay at 3 colours.** Six, five and four colours died out or barely survived in his test, so 3 colours (252 conditions) became our setting.
- **Compare across runs.** "when we get to the ones that are always red-blue or always red-green, are they also always the same genomes … that's one of my main research questions here."
- **Use Claude to learn the codebase.** "just interrogate Claude. Ask Claude questions about the code base".

**My plan as I told it to Chris:** "I'll just run the simulation, record the genomes present in each phase, and then try to identify whatever the dominant genomes are. Just keep repeating the simulation and then … compare the genomes associated with those, like, similar phases."

### Oct 12: genome colouring, in my fork
I asked:

> "Based on me and chris previous converstion, what is the goal this week for this project?"
>
> "i want you to go find the genotype coloring code in Termites and add those changes in Genome, then ill verfiy it by running the simulation. (let me know where the code is tho at the end)"

**What changed in `SudipKing/hex-het-CA-sudip`:**
- I ported Hex-Turmites' genome-to-colour projection into `hexca.js`.
- Press **G** to toggle between cell colours and genome colours.
- I ran it to check that it worked, then asked "it does, how wil this help identify the behavior". The answer: the colours show which genomes spread or disappear, but similar colours don't prove similar behaviour.

### Oct 15–16: Chris's new visualizers, third meeting
**Chris's update.** On Oct 15 he pushed "Add rule-space visualizers and genotype coloring" to his repo: `rulespace.js` and `views.js`, with four stats pages (population, rule space, rasters, PCA trajectory). Before our call I asked ChatGPT:

> "Chris pushed some new changes to his version … and here is mine after the color changer … explain what the new chnages are. after that ill start my call with chris"

**My point in the meeting.** My genome colouring is "more just for like a visual identifier". The RGB projection makes distinct genomes look alike, so similar colours don't mean similar genomes.

**What Chris asked for:**
1. **Merge both codebases** "this week yeah so that we have both".
2. **Founder lineage colouring.** "every organism when they get born they get assigned a unique ID … we're going to track every organism back to one of the first hundred or however many it is … and a color that based on which when they came from".
3. **Parent tracking (stretch).** Record each organism's actual parent so we can build the full tree. He would "save [building the tree] for next week".

**His caveat about lineages:** "the only problem with tracking lineages is it doesn't track how mutations accumulate". He also stressed that a lineage can stay the same while its genome mutates, so later we should compare lineage and genome side by side.

**His longer-term idea:** "if every time we run this simulation certain lineages or genomes are always dying off and other ones are not that's the data we eventually [want]".

### Oct 17: a new repo, the merge, IDs, genealogy export and explorer (ChatGPT + GitHub connector)
**Decision: a new repo.** Rather than working in Chris's repo or my old fork, I started this one:

> "Since the 1st current goal is to merge both versions, i made a new repo called Genome-Phase-Transition. Since Chris doesnt wanna do any work on his repo, merge his one and the one i forked, combine both making sure his and my codes dont conflict, commit it to GPT"
>
> "ok do that. u can also read the past chats to figure out what chris wants and plan according to that."

**Merge with founder-lineage colouring (`c6b83a5`).**
- The agent found that Chris's latest version already had genotype colouring, plus the new rule-space, raster and trajectory tools. So the merge kept his code as the base and added my research features on top.
- New: **founder lineage colouring (L)**. Every cell alive at the start gets a founder ID, and its descendants inherit it even when their genomes mutate.
- `spec.md`, `research-notes/CLAUDE.md` and Chris's visualizers came in with the merge.

**Report (`14082d3`).** I committed the Oct 9 repository report as `repository-report.txt`.

**ID and parent tracking (`0dbc56c`, `51517c7`).**
- I tested both colourings, then asked: "both work, 1st part is done, what does Chris want after this" and "Okay start that".
- Each organism now has an `organismId`, a `parentId` (-1 for founders) and an `ancestorId` (its founder).
- I asked: "Run simulations so we're on the same page, use subagent if u need, Yes everything is working."
- Result: 12 headless runs of 200 ticks (sync, async, history rules, 4 colours) produced 9,252 births and 0 tracking errors.

**Permanent genealogy and family-tree export (`d1ad833`, `55082e0`, `eeacbe6`).**
- Prompt: "lets start the reconstructing".
- Every organism gets a permanent record with birth and death ticks that survives its death. An **Export family tree JSON** button saves the records.

**Genealogy explorer (`7c9d87e`, `2cd5e21`).**
- Prompt: "the export works, start the next part".
- `genealogy.html` loads the export and shows each founder's family tree by generation, with per-founder stats (descendants, alive or dead, depth, extinct or not) and a CSV download.
- The README notes that these are descriptive measures, not proof that a lineage caused a phase change.

### Oct 23 (morning): moving to Claude Code, first real data
**Catching Claude Code up.** I moved the project into Claude Code and pointed it at the history: "go thru all the chats so youre udate don the project". It read all 9 chats in the ChatGPT project through the browser, after I signed in myself. It noted that ChatGPT's own status summary was out of date: it still listed family-tree reconstruction as "next", although that work was done.

**First runs.** I asked: "run genealogy and some simulations and analyze it".

**What was built:**
- **Headless batch runner (`c3712ed`, folder `analysis/`).**
  - `run.js` loads the repo's own simulation files into Node with a seeded RNG, so a seed reproduces a run exactly.
  - Every 1,000 ticks it records the existing stats, distinct-genome counts, the consensus genome, the rule-frequency matrix, how often each condition is met, and the surviving founders.
  - `analyze.js` splits each run into phases. `symmetry.js` compares final genomes across runs, allowing the colours to be relabelled. `build-report.js` with `template.html` and `findings.json` renders a progress report.
- **Large-export fix (`7f701b5`).** `genealogy.html` crashed on exports beyond about 1k ticks. It passed arrays of hundreds of thousands of values into `Math.max` and `Math.min`, which overflowed the call stack. Using loops instead, a 5,000-tick export (740k organisms, 149 MB) now loads in about 7 seconds.

**Follow-up prompts:** "so what does all this mean", then "is it done? write a report on it so i can tell chris wats been done so far". The result was a private "Genome–Phase Progress Report" artifact (also saved as `report.html` in my OneDrive analysis folder).

**Committing and pushing:**
- "yes commit both, makke sure only my name is on the repo, i should be the sole author". Claude set this repo's local git identity to my GitHub noreply address and left my global config unchanged.
- "Push fix/genealogy-large-exports to GitHub (origin). Don't merge into main."
- I reviewed it and merged it myself as **PR #1** (`f2818d5`, Oct 23).

**Findings from 7 seeded runs** (seeds 1–7, default parameters, 100×60 grid, about 1,200 founders):
- **Same three phases every run.** Three-colour churn (about 40% alive, about 66 rules) → two-colour churn (about 50 rules) → a full two-colour grid (100% alive, about 21 rules), reached after 30k–161k ticks. Three runs ended red-blue, two red-green and two green-blue.
- **The final phase has a genome signature.** After relabelling the colours, any two runs' final consensus genomes match about 70% (range 43–92%). Without relabelling, runs ending in different colour pairs shared 0%.
  - 4 rules appear in all 7 runs and 5 more in 6 of 7. All are rules for fully surrounded cells, and 88% say "flip to the other colour". The "immortal" final state is a grid of blinkers.
- **The last transition is a switch in which situations the genome covers.** In the phase just before the grid fills, 4–27% of consensus rules are for fully surrounded cells. Once it fills, 85–96% are.
- **No single dominant genome.** About 5,600 of 6,000 cells carry a unique genome, and the most common one covers under 0.5% of cells. I switched to the *consensus genome* (the majority output for each of the 252 conditions) instead.
- **Founder lineage can't separate the phases.** In every run one founder out of about 1,200 takes over *before* the final transition. This is Chris's "survival of the first".
- **The family tree is too big for whole runs.** It grows by about 150 records per tick, so a full run would be several GB.

### Oct 23 (afternoon): fourth meeting, linking genomes to lineage
I recorded this meeting in the "Research status update" chat. Chris asked for:
1. **Export genomes alongside the lineage.** The genealogy JSON was already about 450 MB ("please don't make too many of these JSONs"). He suggested inverting the map, from ID → genome to genome → list of IDs, so each genome is stored once. He gave a caveat himself: "If every genome is unique, then this wouldn't make any change."
2. **Cluster two ways.** By lineage (how far back two organisms share an ancestor) and by genome. He wants to see whether there is selection or just drift.
3. **Make the tree readable.** "Claude made a first stab, but it was kind of awful … find some good graph visualization techniques … so it's not a big tangled mess".
4. **Look for extinctions at phase changes.** "it's very likely that those phase transitions might also be characterized by extinction events … That's a story, and that's almost a paper right there." For now we want correlation; causation comes later, and he suspects a feedback loop.
5. **Join the files in the dashboard.** The files can stay separate, and the dashboard joins them as needed.

### Oct 30: genome deltas, lineage explorer, lineage analysis (Claude Code)
I pointed Claude at the new meeting notes:

> "their is a new chat in the gpt poject file, read that and start working on the goal"

**I measured before building Chris's idea.**
- `analysis/measure-genomes.js` showed that 98% of 1.8M organisms in a 10k-tick run have a genome no other organism shares. The inverted map would save only about 2%, the case Chris had warned about.
- A child differs from its parent in only about 4 of 252 rules.

**OrganismLog genome deltas (`e89c770`).**
- `organismlog.js` stores, for each organism, only the rules that changed from its parent, in growable typed arrays. Any genome is rebuilt by replaying the changes from the founder down.
- A new **Export lineage + genomes** button writes a compact `.hxlg` file, about 30× smaller than storing whole genomes. A run with about 5M organisms is 160 MB *with* every genome; the old genealogy JSON for that run was 450 MB *without* genomes.
- `test-lineage-export.js` rebuilds every living genome from the file and matches the simulation exactly in 4 settings (sync, async, history rules, 4 colours). Seeded runs stay identical, so recording doesn't change the simulation.

**Lineage explorer (`b1a6174`).** In answer to "make the tree readable", `lineage-explorer.html` uses techniques from evolutionary biology:
- a **Muller plot** of the major clades, coloured by genome, so both of Chris's groupings appear in one picture
- a **time-to-most-recent-common-ancestor (TMRCA)** curve
- births and deaths over time
- a zoomable, pannable tree of the living population
- a lineage-distance vs genome-distance plot
- an extinction-events table lined up with the phase changes

`analysis/lineage-report.js` runs the same analysis in batch.

**Pushing and merging:**
- "Push branch". I merged **PR #2** myself (`50c1776`, Oct 30).
- "its commited now, explain what u did. then whats the analysis after the simulations"
- "can you put the project here so i can open the html files without github everytime". The clone and the analysis outputs now live in my OneDrive "Phase Transition research" folder. I kept a full seed-6 `.hxlg` there as an example and deleted about 1 GB of other run files, since any run can be regenerated from its seed.

**Findings from 5 long runs** (seeds 2, 4, 5, 6, 7, run through the final transition):

| Seed | Grid fills | One clade takes over (TMRCA collapse) | Final-phase rules its founder already had | Contemporaries (avg) |
|---|---|---|---|---|
| 6 | 30k | 29k | 7 of 14 | 1.8 |
| 4 | 51k | 48k | 7 of 12 | 1.5 |
| 2 | 46k | 49k | 4 of 11 | 2.3 |
| 7 | 37k | 43k | 7 of 11 | 2.6 |
| 5 | 71k | not detected | — | — |

- **A sweep comes with the final phase.** In 4 of 5 runs a single clade replaced everyone within about 6k ticks of the grid filling, sometimes just before and sometimes just after. That fits Chris's feedback-loop idea better than a one-way cause.
- **The winners were prepared early.** In 3 runs the sweeping clade's founder already carried about 7 of the final-phase rules, and at most 4 of 300 sampled contemporaries matched it.
- **Mostly drift, with episodic selection.** The mutations that founded successful clades are not concentrated in the rules cells actually use. Once the grid is full, though, those rules are strongly conserved: about 2–3.5 differences between cells whose common ancestor lived about 20k ticks ago.
- **Extinctions.** The biggest extinction event lands at the final fill in 4 of 5 runs. Smaller extinctions cluster around phase changes only slightly more than chance would predict.
- **Caveat.** This is 5 runs at default settings, and these are correlations only.

### Oct 30 (later): genome population snapshots (branch `research/genome-snapshots-local`, not merged)
These 10 commits were made after the Claude Code session above ended, and I have no Claude Code transcript for them. They change the batch runner to:
- write every sample's full living-genome distribution to `genome-snapshots.jsonl`, storing each unique genome once with its count, through a small append-only storage adapter (`analysis/storage/local.js`) meant to be swappable for a database later
- also record random deaths per reporting interval (`stats.js`)
- replace "stop once the grid has stayed full for `--settle` samples" with "stop once births and total deaths stay below `--quietThreshold` for `--stableSamples` samples", and lower the default tick cap from 100k to 50k

Three of the commits are follow-up fixes for the JSON Lines newline delimiter, two of them empty. The branch is still on origin and has not been merged into `main`.

### Nov 5: 60 more runs, deeper data collection, Phase Explorer (Claude Code)
Prompt: "run more simulation. Improve the graphs and data collection visualizers so more indepth analysis can be done. Then write a report based on the simulations."

**Data collection.**
- `hexca.js` now counts rule deaths per condition and carries a resettable clade ID. Neither draws random numbers, so seeds reproduce the same runs.
- `run.js` records rules, encounters and rule deaths split by empty neighbours, genome-size quantiles, the effective number of genomes and clade sweeps per 5k window. It can stop 15k ticks after the grid fills (`--settleFull`).
- `batch.js` runs a plan of seeds and parameter sets in parallel.

**Analysis and visualizers.** `deep.js` sorts phases into classes, times events relative to the fill, compares phase genomes across runs and aligns metrics on the fill. `explorer.js` builds the Phase Explorer page from it.

**Runs.** 24 default seeds (201–224), plus 6 seeds (301–306) each at mutation rate 0.005 and 0.02, random death 0.0005 and 0.002, and k 0.05 and 0.2.

**Findings:**
- **All 24 default runs** filled and ended with two colours (median fill 57k ticks). 17 went three-colour churn → two-colour churn → two-colour full.
- **Phases have characteristic genomes.** The same phase in different runs matches 48–61% after relabelling. Different phases match 0–28%, and a run's full-grid genome matches its own two-colour-churn genome only 7%.
- **The final genome shift follows the fill.** Rule deaths collapse at +1k ticks and genome size is halfway down at +2k: once every cell is surrounded, unused rules are pruned at each birth. Before the fill, rules for surrounded cells rise to about a third of the consensus, and a clade sweep is 3.5× as likely in the window containing the fill.
- **Mutation rate controls whether the grid fills.** At 0.02 no run filled in 200k ticks; at 0.005 runs filled at a median of 14k.

---

## 4. Findings that steered the next step

| Finding | What it changed |
|---|---|
| Similar genome colours don't mean similar genomes (Oct 16) | Added founder lineage colouring and ID tracking instead of relying on colour alone |
| Founder lineage is decided long before the final phase (Oct 23) | Lineage colouring alone can't separate the phases; proposed resetting founder IDs at each phase change |
| Nearly every genome is unique (Oct 23, measured again Oct 30) | Switched from "most popular genome" to the consensus genome; replaced the genome → IDs map with parent deltas |
| The final phase's rules match about 70% across runs after colour relabelling | Early support for Chris's hypothesis; the next test is seeding a grid with the about-14-rule blinker genome |
| The genealogy grows by about 150 records per tick | Full-run trees in JSON were dropped in favour of the compact `.hxlg` format |
| A clade sweep comes within about 6k ticks of the grid filling | Points toward a feedback-loop explanation; needs more seeds and direct tests |

## 5. Decisions

- **3 colours** (252 conditions), because 4 or more colours die out or barely survive.
- **Research paper only** as the final deliverable, with no presentation. The project is 2 credits in autumn and may continue in winter and spring.
- **A new repo (Genome-Phase-Transition),** keeping Chris's code as the base and adding my features on top, without touching his repo.
- **Lineage and genome are tracked separately** (`ancestorId` vs genome), because a lineage can persist while its genome mutates.
- **Consensus genome** instead of exact "dominant genome" matching.
- **Compare across runs with colour relabelling,** because red-blue, red-green and green-blue endings turned out to be the same rulebook under different colour names.
- **Measure before building:** I tested Chris's inverted-map idea against real data and chose delta encoding instead.
- **I'm the sole author:** commits use my GitHub identity, set in this repo only, with no AI co-author lines.
- **Branches and PRs, merged by me:** the agents push feature branches, and I review and merge (PR #1 on Oct 23, PR #2 on Oct 30).
- **Keep the bulky run data out of git** (`analysis/.gitignore`). Any run can be regenerated from its seed.

## 6. Abandoned or changed ideas

- **The other projects in Chris's zoo** (traffic, Plant World, emergent boundaries, Hex-Turmites itself), in favour of the hex CA.
- **A presentation** as a deliverable.
- **Working directly in Chris's repo** or my old fork.
- **Finding "the most popular genome",** because exact genomes are almost all unique.
- **Chris's genome → IDs inverted map,** because it would save only about 2% when 98% of genomes are unique. I replaced it with parent-delta encoding.
- **Whole-run family-tree JSON exports,** because they would be several GB. The genealogy explorer still works for short windows.
- **Founder lineage colouring as a phase signal,** because one founder wins before the phases I care about.
- **The assumption of about 100 founders.** Every initially living cell is a founder, which is about 1,200.
- **Stopping runs once the grid has stayed full.** On the unmerged snapshot branch it is replaced by stopping when births and deaths stay low.

## 7. Next steps

1. Verify and decide on the `research/genome-snapshots-local` branch, then merge it or drop it.
2. Run 20–30 seeds, then sweep mutation rate, `k` and random death.
3. Reset lineage IDs at each phase change, to see which sub-clade drives the next transition.
4. Test the signature directly: seed a grid with the relabelled blinker genome, then knock out single rules.
5. Characterize the middle two-colour-churn phase's genomes across runs.
6. Bring the 252-rule consensus strip into the simulator as a live stats page.
7. Write the paper.
