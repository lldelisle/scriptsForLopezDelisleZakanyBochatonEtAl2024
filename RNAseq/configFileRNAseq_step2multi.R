### Required for all steps ###
RNAseqFunctionPath <- "/home/ldelisle/softwares//rnaseq_rscripts/RNAseqFunctions.R"
# This file should be a tabulated file with at least one column called
# 'sample'.
samplesPlan <- "/home/ldelisle/softwares/scriptsForLopezDelisleZakanyBochatonEtAl2024/RNAseq/samplesplan_Gastruloids.txt"

#### STEP 2 - DESEQ 2 ANALYSIS ####
# Required
tableWithCounts <- "/scratch/ldelisle/HoxBstudy/RNAseq//mergedTables/AllHTSeqCounts_subset_onlypc.txt"
# Specify here the name of the column which contains the gene IDs (they need to
# be unique).
geneIDColCounts <- "Ens_ID"
# For the DESeq2 analysis you need to specify a factor on which you want to do
# the analysis: This needs to be a name of a column of the samplesPlan file.
# Here you want to do multiple analyses at once.
# You describe all what you want into the variable 'all.anayses':
# list(factorToStudy = list(loopingVariable = list(subsetting)))
all.analyses <- list(
  "genoclone" = list(
      "genoclone" = list(
        "time" = "120"
      ),
      "genoclone" = list(
        "genoclone" = c("delBdeli9-13", "delBins5-10del5-10"),
        "time" = "120"
      ),
      "genoclone" = list(
        "genoclone" = c("delBdeli9-13", "wt"),
        "time" = "120"
      )
  )
)
# First, using Line as factor, for each Line which is not the reference (here there is only one but there could be more), using only samples where Tissue is DFL.
# Then, using Tissue as factor, for each Tissue which is not the reference (here only one), using only samples where Line is Wt.

# all.analyses <- list(
#   "Line" = list(
#     "Line" = list("Tissue" = "Limbs"),
#     "Line" = list("Tissue" = "Trunk")
#   )
# )
# In the given example, it will use the factor 'Line' and do:
# 1. For each Line which is not the reference (first one in samples plan), comparison to the reference using only samples where Tissue is Limbs.
# 2. For each Line which is not the reference (first one in samples plan), comparison to the reference using only samples where Tissue is Trunk.


# Optional
# This can be table from cufflinks or cuffdiff or Biomart to annotate genes.
# You will need to choose a file with at least one column with the Ensembl Gene IDs.
# tableWithAnnotations <- "~/rnaseq_rscripts/outputs/mergedTables/AllCufflinks_Simplified_norm.txt"
# Specify here the name of the column which contains the gene IDs (it must
# match with the content of the geneID from the table with counts).
# geneIDColInAnnotations <- "gene_id"
# You can also provide a gtf:
gtfFile <- "/scratch/ldelisle/HoxBstudy/RNAseq/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf"

# A directory with outputs
pathForDESeq2 <- "/scratch/ldelisle/HoxBstudy/RNAseq/DESeq2_multi_onlypc/"
# Set the threshold for significance here
log2FC.threshold <- log2(1.5)

# For Get_summary_some_genes.R
my_gene_ids <- c("ENSMUSG00000049604", "ENSMUSG00000020875", "ENSMUSG00000018973")
