## Proteome-wide structure-based accessibility analysis of ligandable and detectable cysteines in chemoproteomic datasets
Code and analysis scripts for [White et al.](LINK TO BE ADDED) - describing integrated structural analysis of public cysteine profiling datasets with AlphaFold-predicted amino acid accessibility. 

### Amino acid accessibility prediction
Individual amino acid accessibility predictions (and IDR prediction) were adapted from [Bludau et al](https://doi.org/10.1371/journal.pbio.3001636) and are re-distributed in accordance with the Apache 2.0 License. This script allows the calculation of prediction-aware part-sphere exposure (pPSE) values from AlphaFold-predicted proteomes from [model organism proteomes](https://alphafold.ebi.ac.uk/download) as well as IDR prediction and secondary structure annotation.

### Input data
All cysteine profiling data used in this study was downloaded from the respective publications as Supplementary Information files and reformatted for comparison. Reformatted data containing all cysteine residues detected across public cysteine profiling datasets can be generated directly from Supplementary files with the 'reformat_figure2.R' script. Additionally, fragment screening datasets can be reformatted and combined with 'reformat_figure3.R'. 

### Figure generation
All figures from the corresponding [pre-print](LINK TO BE ADDED) can be generated from source data using the 'figures.Rmd' document. 

### Acknowledgements
pPSE calculation code is adapted from [Bludau et al](https://doi.org/10.1371/journal.pbio.3001636)
