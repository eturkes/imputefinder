# PXD020717 Tau-phosphosite fixture

`PXD020717_tau_phosphosites.csv` is a compact extract from the independent
human brain phosphoproteomics study archived under PRIDE accession PXD020717.
It contains three control and three Alzheimer's disease soluble BA46 samples
and seven MAPT phosphosites.

The values reproduce MaxQuant `Intensity` columns from
`txt/Phospho (STY)Sites.txt`. Zero entries are stored as missing values. Positive
entries are log2-transformed and rounded to three decimals.

Run `data-raw/build-pxd020717-example.py` from the package root to regenerate the
fixture. The script uses HTTP range requests, so it reads only the ZIP directory
and the 165 kB compressed source member rather than downloading the 934 MB
archive. It verifies the archive size, source-member size and CRC, selected
P10636-8 positions, and the generated fixture SHA256
`eb65d2e543da1ea97010f76f7f2d3c0e234001a312adb429d8c71113b6c55e1d`.

The samples, sites, and vignette cutoff are selected to demonstrate
classification behavior. They are not a disease-effect analysis.
