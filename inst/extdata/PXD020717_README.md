# PXD020717 Tau-phosphosite fixture

`PXD020717_tau_phosphosites.csv` contains seven MAPT phosphosites from the
PXD020717 human brain study in PRIDE. It covers three control and three
Alzheimer's disease soluble BA46 samples.

The values reproduce MaxQuant `Intensity` columns from
`txt/Phospho (STY)Sites.txt`. Zero entries are stored as missing values. Positive
entries are log2-transformed and rounded to three decimals.

Run `data-raw/build-pxd020717-example.py` from the package root to regenerate the
fixture. HTTP range requests read only the ZIP directory and 165 kB source
member, not the full 934 MB archive. Checks cover archive size, source-member
size and CRC, P10636-8 positions, and fixture SHA256
`eb65d2e543da1ea97010f76f7f2d3c0e234001a312adb429d8c71113b6c55e1d`.

The samples, sites, and cutoff demonstrate classification, not disease effects.
