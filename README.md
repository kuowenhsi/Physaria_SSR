# Physaria globosa SSR dataset

This repository contains microsatellite (SSR) genotype data for *Physaria globosa*, extracted manuscript tables, and mapping scripts.

## Project structure

- `data/`: raw genotype files, coordinate table, spreadsheet, and manuscript tables exported to CSV
- `script/`: R script and R Markdown workflow for the Table 1 map
- `figure/`: rendered figure outputs
- `manuscript/`: manuscript source files and rendered HTML report

## Data files

### `data/PG_All.csv`

Primary genotype matrix in CSV format.

- Total samples: 352
- Total loci: 18
- Header structure:
  - Row 1 contains summary metadata. The leading values appear to be `18` loci and `352` samples, followed by additional numeric summaries.
  - Row 2 contains species/site metadata beginning with `Physaria globosa`.
  - Row 3 contains the column labels used for the genotype table.
- Data rows begin on row 4.
- Column layout:
  - Column 1: `sample`
  - Column 2: `site`
  - Columns 3-38: diploid allele calls, stored as two columns per locus
- Loci, in order:
  - `88`, `13`, `37`, `79`, `8`, `24`, `4`, `29`, `23`, `33a`, `33b`, `2`, `30`, `122`, `12`, `82`, `26`, `39`
- Missing genotype calls are encoded as `0`.

Example data row:

```text
Pg346_02,346,136,137,238,250,...
```

### `data/PG_All_structure.txt`

Whitespace-delimited version of the same SSR matrix prepared for use with STRUCTURE-style population genetic workflows.

- Total samples: 352
- Total loci: 18
- Line 1 lists the locus names in genotype order.
- Each subsequent line contains:
  - Sample identifier, typically derived from the `.fsa` source filename
  - A numeric population/group code
  - A constant `0` field
  - Two allele values for each of the 18 loci
- Missing genotype calls are encoded as `-9`.

Example data row:

```text
883779813_Pg346_02.fsa 1 0 136 137 238 250 ...
```

### `data/pglo.coordinates.txt`

Tab-delimited site coordinate table.

- Rows: 14
- Columns:
  - Site ID
  - Latitude
  - Longitude

This file does not cover every site label present in `data/PG_All.csv`. Site labels present in the genotype matrix but absent from the coordinate table are:

- `348`
- `353`
- `355`
- `p346`

Additional data products in `data/`:

- `coordsandadmixprpo.pglo.xlsx`
- `Table_1.csv` through `Table_6.csv`

## Notes on identifiers and consistency

- The genotype files represent the same 352 individuals and 18 loci, but use different missing-data codes:
  - `data/PG_All.csv`: `0`
  - `data/PG_All_structure.txt`: `-9`
- The `site` field in `data/PG_All.csv` is mostly numeric, but one sample (`Pg346_01`) is labeled `p346` while the rest of that group uses `346`.
- The second metadata row of `data/PG_All.csv` includes a site label `253`, while the sample rows include `353`. That discrepancy should be checked against the original data source before automated downstream use.

## Suggested downstream use

- Use `data/PG_All.csv` for general tabular processing in R, Python, or spreadsheets.
- Use `data/PG_All_structure.txt` for STRUCTURE or other software expecting a whitespace-delimited SSR matrix.
- Use `data/pglo.coordinates.txt` only after confirming how missing site IDs and the `p346` label should map to sampling localities.
- Use `script/plot_Table_1_map.R` for the publication-style map and `script/plot_Table_1_map.Rmd` for the stepwise HTML workflow.
