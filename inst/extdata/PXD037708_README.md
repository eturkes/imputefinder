# PXD037708 vignette fixture

`PXD037708_brain_EV_log2.tsv.gz` contains the processed protein-group matrix
used in the ImputeFinder real-data example. The matrix contains 6,107 protein
groups and 61 human Alzheimer’s disease brain extracellular-vesicle fraction
samples. `PXD037708_brain_EV_samples.csv` records the donor and density fraction
for each sample.

The fixture was generated from the `PROTEINS` sheet in
`20221024_TAU_EV64_PRIDE.xlsx`. Positive intensities were transformed with
log2. Blank, non-finite, and non-positive entries were stored as missing values.
Samples `A_F8`, `E_F8`, and `F_F8` were excluded to match the source analysis.
The vignette's eight manual cutoffs also match that analysis at commit
`2dfb54176e82ef049084793586c9d15e4f2cdbdd`.

Source: https://www.ebi.ac.uk/pride/archive/projects/PXD037708

Analysis source: https://github.com/duff-lab-team/AD-EV-characterisation

Workbook SHA-256:
`a018b3a20904f07f49928014056d0c005dab14b933e9d4b6cceadd28764cc34b`

Regenerate the files from the repository root:

```bash
python data-raw/build-pxd037708-example.py
```

The regeneration script requires Python and openpyxl.
