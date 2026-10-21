# Physaria SSR Project

This repository contains the microsatellite data, analysis scripts, manuscript-facing reports, and rendered figures for the *Physaria globosa* SSR project.

## Directory layout

- `data/`
  Source tables, genotype matrices, metadata crosswalks, map layers, and bundled Natural Earth layers.
- `script/`
  Executable R workflows grouped by analysis type.
- `manuscript/`
  R Markdown reports and rendered HTML companions grouped to match the script structure.
- `figure/`
  Rendered figure outputs written by the analysis scripts.
- `results/`
  Intermediate software outputs used by downstream plotting workflows.
- `tools/`
  Third-party helper software kept with the project.
- `photo/`
  Project imagery and reference photos.

## Script organization

- `script/data_prep/`
  Metadata building and crosswalk generation.
- `script/maps/`
  Map generation and hydrography/overlay workflows.
- `script/ordination/`
  PCA, DAPC, and genetic-distance tree workflows.
- `script/structure/`
  STRUCTURE and structureHarvester visualization workflows.
- `script/tables/`
  Table-derived plotting workflows.
- `script/utils/`
  Shared helper functions, including Wen-label utilities.

## Report organization

- `manuscript/maps/`
  Map-oriented R Markdown reports and their rendered HTML files.
- `manuscript/ordination/`
  Ordination/tree reports and HTML outputs.
- `manuscript/structure/`
  STRUCTURE and structureHarvester reports and HTML outputs.
- `manuscript/tables/`
  Table-summary reports and HTML outputs.

## Naming conventions

- Canonical population codes remain in the `Population` field where they are part of the original table definitions.
- `Wen_label` is used as the standardized display label where appropriate.
- Outgroup labels are normalized to the same display style, for example `AR_01`, `MO_06`, and `TX_01`.

## Working expectation

When a script and a report represent the same workflow, they now live in matching thematic folders and resolve paths from the project root instead of relying on the old flat layout.
