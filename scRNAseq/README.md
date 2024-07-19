# scRNAseq

## Inputs

The matrices were generated as in [Mayran et al. 2023](https://doi.org/10.1101/2023.11.22.568291) and the WT sample comes from this analysis. The code to generate matrices is available [here](https://github.com/MayranA/allScriptsFromMayranEtAl2023/blob/f766ff5a58fbe115e85ab17363bf4f66ee569ba0/scRNAseq/fastq_to_matrices).

Both samples have been registered into a [csv file](./metadata.csv).

Once matrices have been generated in galaxy they were downloaded and organized like this:

```bash
.
└── GEX
    ├── 144h.WT.rep1
    │   ├── barcodes.tsv
    │   ├── genes.tsv
    │   └── matrix.mtx
    └── 144h.Deli9-13_het.rep1
        ├── barcodes.tsv
        ├── genes.tsv
        └── matrix.mtx
```

## R configuration

The scripts have been launched on a RStudio server.

The session information is available [here](./sessionInfo.txt).

## Pipeline

The pipeline has been copied and slightly modified from https://github.com/MayranA/allScriptsFromMayranEtAl2023/tree/44924b98a177decd43cf149da9be078f8143639b/scRNAseq/matrices_to_plots

[Step1](./Step1.Seurat.Demultiplexing.Analysis.R) has been run on all available samples. This R script generates a RDS for each sample. If needed it demultiplexes CellPlex experiments (here not applicable).

[Step2](./Step2.Seurat.Analysis.and.Merging.R) has been to generate a merged Seurat object with the wt and the mutant sample.

[Step3](./Step3.qmd) contains all code used to generate figures relative to scRNAseq.

Custom functions have been collected into a [single file](./scRNAseqFunctions.R).

## baredSC

baredSC was run on SCITAS (EPFL cluster). The steps are described [here](./baredSC/run_all_baredSC.sh).

The convergence was checked manually and the plots were run locally using [this quarto](./plot_2d_2024.qmd).
