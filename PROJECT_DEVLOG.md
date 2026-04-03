# cpmd-interface-screen Project Dev Log

## Project identity

- Repo: https://github.com/junminchen/cpmd-interface-screen
- Focus: constant-potential molecular dynamics for interface-aware electrolyte screening
- Core direction:

`molecule / formulation -> constant-potential interface MD -> interfacial descriptors -> screening / supervised learning`

## Why this project matters

The key scientific goal is not just to simulate one electrolyte system.

The larger goal is to build an **interface-descriptor-driven electrolyte screening framework** under constant-potential MD, so that candidate solvents / additives / anions can be compared using mechanistically meaningful and computationally tractable descriptors.

This project should ideally connect three layers:

1. **molecular descriptors**
   - charge distribution
   - dipole / polarizability
   - solvation-related properties
2. **interfacial descriptors under constant potential**
   - adsorption / depletion
   - interfacial enrichment
   - orientation response to applied potential
   - anion/cation redistribution near the electrode
3. **performance hypotheses / labels**
   - CE
   - capacity retention
   - interfacial stability
   - SEI / CEI formation tendency

## Current repository status

Already present in the repo:

- descriptor extraction workflow
- descriptor aggregation workflow
- standard constant-potential MD protocol
- campaign template
- literature CE data ingestion / curation

Current descriptor families already mentioned in repo docs include:

- interfacial enrichment
- first-layer contact probability and residence time
- Li+ first-shell replacement fractions
- interfacial capacitance from electrode charge logs

This means the repo already has a good basis for moving toward descriptor-to-target modeling.

## New brainstorming triggered by the TMSB / dual-descriptor paper

A recent paper on **dual-descriptor-guided design of electric field-sensitive solubilizing additives** is especially relevant to this repo.

The most useful lesson is not the specific additive itself, but the workflow logic:

- first use **low-cost molecular descriptors** to pre-screen candidates
- then use **high-cost interface simulations** to validate mechanism

For this project, that suggests the following development principle:

`molecular pre-screening -> constant-potential interface simulation -> interface descriptor extraction -> mechanism-guided ranking`

### Key conceptual upgrade

The project should move beyond ordinary bulk descriptors.

The real opportunity is to define **field-response interfacial descriptors** under constant potential, for example:

- how a molecule reorients near the electrode under different voltages
- whether an additive is attracted to or repelled from the interface
- whether the additive promotes interfacial anion enrichment
- how the first-layer composition changes with applied potential
- how strongly these quantities change with voltage

This is likely more valuable than using only static molecular features.

## Working hypothesis for the project

A useful screening framework may require two coupled descriptor blocks:

### 1. Bulk-side descriptors

Examples:

- coordination number
- salt dissociation tendency
- local solvation structure
- diffusion-related proxies
- bulk anion ratio

### 2. Interface-side descriptors

Examples:

- interfacial number density profile `rho(z)`
- first-layer enrichment factor
- orientation order parameter
- residence time near electrode / interfacial Li+
- adsorption/depletion tendency
- interfacial anion ratio
- voltage-response slopes such as `d descriptor / dV`

The most important scientific step is to understand how bulk-side descriptors and interface-side descriptors couple together.

## Immediate scientific opportunity

The TMSB-style mechanism suggests a benchmark question that fits this repo well:

- can constant-potential MD reproduce additive depletion from a positively polarized interface?
- can it quantify voltage-dependent NO3- or other anion enrichment near the interface?
- can the repo extract descriptors that are more informative than simple Mulliken charge?

If yes, this repo becomes more than a simulation toolkit.
It becomes a platform for **descriptor-driven electrolyte design under electrochemical boundary conditions**.

## Concrete next steps

### Phase A — first descriptor expansion

Build a first explicit list of constant-potential interfacial descriptors to standardize across all cases.

Priority candidates:

- interfacial density profile `rho(z)`
- first-layer enrichment factor
- orientation order parameter
- additive / anion residence time
- Li+-anion coordination at interface
- bulk-reference coordination metrics
- voltage-response slopes such as:
  - `d_interfacial_anion_ratio_dV`
  - `d_enrichment_dV`
  - `d_orientation_dV`

### Phase B — minimal benchmark system

Do not start from the most complicated realistic system.

Instead, define a minimal benchmark:

- one electrode model
- one base electrolyte
- one or two additives
- two or three voltage points

The first goal is **descriptor discriminability**, not full electrochemical realism.

### Phase C — descriptor mapping logic

Formalize the mapping:

`molecular descriptor -> interfacial response descriptor -> electrochemical hypothesis`

Example chain:

- Mulliken charge / dipole / polarizability
- adsorption vs depletion tendency
- anion redistribution near electrode
- likely SEI / CEI precursor environment
- expected CE / stability trend

### Phase D — data/model integration

After descriptors are stable:

- align formulation ids
- assemble descriptor table
- join with CE / retention labels
- train simple tabular baselines first
- only consider neural models later if data volume supports it

## Development priorities right now

The immediate priority is **not** fancy ML architecture.

The immediate priority is:

1. descriptor definition stability
2. protocol consistency across voltages and formulations
3. interface-specific anion / additive metrics
4. reproducible extraction pipeline
5. clean label alignment
6. baseline interpretable modeling

## Notes for future updates

This file should be updated whenever one of the following happens:

- a new descriptor family is added
- a benchmark system is fixed
- a protocol decision is made
- a literature paper changes project direction
- the descriptor-to-label modeling path is revised
