# Set working directory (github directory)
setwd("~/Documents/mygit/scriptsForLopezDelisleZakanyBochatonEtAl2024/")
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggforce")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("scales")
# Get good fontfamily
safelyLoadAPackageInCRANorBioconductor("extrafont")
# Only once
# font_import()
# loadfonts(device = "all")
fontFamily <- "DejaVu Sans"

# Parameters
samplesPlan <- "RNAseq/samplesplan_embryos.txt"
tableWithNormalizedExpression <- "RNAseq/outputs/AllCufflinks_Simplified_E.txt.gz"
gtfFile <- "publishedData/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz"
tableWithDESeq2Summary <- "RNAseq/outputs/summary_E.txt.gz"
outputFolder <- "RNAseq/outputs/"
my.pretty.names <- c("wt", "Del(i9-13)", "Del(i9-13):hd", "Hoxb13hd")
my.pretty.names.formatted <- c("'wt'", "'Del(i9-13)'", "'Del(i9-13):hd'", "italic('Hoxb13'^{'hd'})")
names(my.pretty.names) <- c("wt", "deli9-13", "deli9-13hd", "b13hd")
names(my.pretty.names.formatted) <- c("wt", "deli9-13", "deli9-13hd", "b13hd")
fixedColors <- list(pretty.geno = c("red", "green4", "grey50", "cadetblue"))
names(fixedColors$pretty.geno) <- my.pretty.names
fixedColors$pretty.geno.formatted <- fixedColors$pretty.geno
names(fixedColors$pretty.geno.formatted) <- my.pretty.names.formatted
# Create outputFolder
dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)

# Read samples plan
samplesPlanDF <- read.delim(samplesPlan, check.names = FALSE)
samplesPlanDF$pretty.geno <- my.pretty.names[samplesPlanDF$genotype]
samplesPlanDF$pretty.geno <- factor(
  samplesPlanDF$pretty.geno,
  levels = my.pretty.names
)
samplesPlanDF$pretty.geno.formatted <- my.pretty.names.formatted[samplesPlanDF$genotype]
samplesPlanDF$pretty.geno.formatted <- factor(
  samplesPlanDF$pretty.geno.formatted,
  levels = my.pretty.names.formatted
)
samplesPlanDF$stage <- factor(
  samplesPlanDF$stage,
  levels = unique(samplesPlanDF$stage)
)
samplesPlanDF$group.label <- paste0(samplesPlanDF$stage, "\n", samplesPlanDF$pretty.geno)
rownames(samplesPlanDF) <- samplesPlanDF$sample
samplesPlanDF$pretty.sample <- paste0(samplesPlanDF$stage, "_", samplesPlanDF$pretty.geno, "_", samplesPlanDF$replicate)
samplesPlanDF$pretty.sample.norep <- paste0(samplesPlanDF$stage, "_", samplesPlanDF$pretty.geno)

# Read summary table
summaryDF <- read.delim(tableWithDESeq2Summary, check.names = FALSE)

# Read FPKM table
expressionDF <- read.delim(tableWithNormalizedExpression, check.names = FALSE)
metaCols <- which(sapply(colnames(expressionDF), function(cn) {
  class(expressionDF[, cn]) != "numeric"
}))
# Remove FPKM_ from column names:
colnames(expressionDF) <- gsub("^FPKM_", "", colnames(expressionDF))
# Asign rownames
rownames(expressionDF) <- expressionDF$gene_id
# We will only consider protein coding genes:
# Read gtf
gtf <- readGFF(gtfFile)
proteinCodingGenes <- unique(gtf$gene_id[gtf$gene_biotype == "protein_coding"])
# We exclude chromosomes X,Y,M
autosomalGenes <- unique(gtf$gene_id[! gtf$seqid %in% c("chrX", "chrY", "chrM")])
# Subset expression with selected genes:
expressionDF <- expressionDF[intersect(proteinCodingGenes, autosomalGenes), ]
samplesToPlot <- intersect(colnames(expressionDF), samplesPlanDF$sample)
# Get only values (no meta)
data <- expressionDF[, samplesToPlot]
# Also subset samples plan if needed
samplesPlanDF <- samplesPlanDF[samplesToPlot, ]
# Exclude genes with no expression
sumperline <- apply(data, 1, sum)
nonZdata <- data[sumperline != 0, ]
# Log transform the data
ldata <- log2(nonZdata + 1)
rldata <- ldata
# Keep only the genes with more variation:
rldata <- ldata[
  order(
    apply(ldata, 1, var),
    decreasing = TRUE
  )[seq_len(min(nrow(ldata), 500))],
]
# Compute PCA
sample.pca <- prcomp(t(rldata), center = TRUE, scale. = FALSE)
# Add PC values to samples plan
new.df <- data.frame(
  samplesPlanDF, sample.pca$x[samplesToPlot, ],
  check.names = FALSE
)
var <- round((sample.pca$sdev)^2 / sum(sample.pca$sdev^2) * 100)
# Plot PC1 PC2
g <- ggplot(new.df, aes(PC1, PC2)) +
  geom_point(aes(color = pretty.geno.formatted, shape = stage), size = 1.5) +
  xlab(paste0("PC1: ", var[1], "% variance")) +
  ylab(paste0("PC2: ", var[2], "% variance")) +
  scale_color_manual("", values = fixedColors$pretty.geno.formatted, label = parse_format()) +
  theme_pubr(legend = "right") +
  geom_mark_ellipse(
    aes(
      color = pretty.geno.formatted,
      fill = NULL,
      label = group.label,
      group = group.label
    ),
    label.fontsize = 6,
    label.fontface = "plain",
    label.margin = margin(0, 0, 0, 0, "mm"),
    label.buffer = unit(0, "mm"),
    con.type = "none",
    con.border = "none",
    label.family = fontFamily,
    label.hjust = 1
  ) +
  theme(
    text = element_text(family = fontFamily, size = 6),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 6)
  )
ggsave("RNAseq/outputs/PCA12_E_ellipses.pdf", g, width = 4, height = 3.5)

# Plot Hoxb13 expression
temp.df <- cbind(
  samplesPlanDF,
  t(ldata[expressionDF$gene_id[expressionDF$gene_short_name == "Hoxb13"], ])
)
colnames(temp.df)[ncol(temp.df)] <- "Hoxb13"

# Retrieve adjusted p-value
comparisons <- unique(
  subset(
    samplesPlanDF, genotype != "wt",
    select = c("pretty.geno.formatted", "genotype", "stage")
  )
)
colnames(comparisons)[1] <- "group1"
comparisons$group2 <- my.pretty.names.formatted['wt']
comparisons$group2.genotype <- "wt"
comparisons <- rbind(
  comparisons,
  subset(
    comparisons,
    stage == "E10.5" & genotype == "deli9-13hd"
  )
)
comparisons$group2[nrow(comparisons)] <- my.pretty.names.formatted["deli9-13"]
comparisons$group2.genotype[nrow(comparisons)] <- "deli9-13"
comparisons$colnames <- with(
  comparisons,
  paste0("genotype_", stage, "_", genotype, "vs", group2.genotype, "_padj")
)
comparisons$colnames[nrow(comparisons)] <-
  grep(
    gsub("genotype_", "", comparisons$colnames[nrow(comparisons)]),
    colnames(summaryDF), value = TRUE
  )
comparisons$p.adj <- t(summaryDF[which(summaryDF$gene_name == "Hoxb13"), comparisons$colnames])

comparisons$p.signif <- cut(
  comparisons$p.adj,
  breaks = c(0, 1e-4, 1e-3, 1e-2, 0.05, 1),
  labels = c("****", "***", "**", "*", "ns")
)
comparisons$p.signif[is.na(comparisons$p.signif)] <- "ns"


g <- ggboxplot(
  temp.df,
  x = "pretty.geno.formatted",
  y = "Hoxb13",
  color = "pretty.geno.formatted",
  palette = fixedColors$pretty.geno.formatted,
  add = "jitter"
) +
  facet_grid(
    . ~ stage,
    scale = "free_x",
    space = "free_x"
  ) +
  stat_pvalue_manual(
    data = comparisons,
    label = "p.signif",
    y.position = c(4.5, 2.8, 3.1, 3.4, 2.8),
    family = fontFamily,
    size = 2
  ) +
  labs(color = "", y = expression(atop(italic("Hoxb13")~"expression","["~log[2](1+FPKM)~"]"))) +
  scale_x_discrete(labels = parse_format()) +
  theme(
    text = element_text(family = fontFamily, size = 6),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), # Put the genotypes vertical
    axis.title.x = element_blank()
  )

ggsave("RNAseq/outputs/Hoxb13_E.pdf", g, width = 3, height = 3.5)

# Plot all Hox genes in a heatmap
hox.genes <- intersect(
  paste0("Hox", rep(letters[1:4], each = 13), rep(1:13, 4)),
  expressionDF$gene_short_name
)

hox.genes.id <-
  expressionDF$gene_id[match(
    hox.genes,
    expressionDF$gene_short_name
  )]

mat <- ldata[
  hox.genes.id,
  samplesToPlot
]

rownames(mat) <- hox.genes

annot.samples <- subset(
  samplesPlanDF,
  select = c("pretty.geno", "stage")
)
colnames(annot.samples) <- c(" ", "stage")

# From https://www.biostars.org/p/400381/
hox.genes.formatted <- lapply(
  hox.genes,
  function(x) bquote(italic(.(x))))

pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annot.samples,
  annotation_colors = list(' ' = fixedColors$pretty.geno),
  file = "RNAseq/outputs/Hox_heatmap_E.pdf",
  labels_row = as.expression(hox.genes.formatted),
  labels_col = samplesPlanDF[colnames(mat), "pretty.sample.norep"],
  fontfamily = fontFamily,
  annotation_names_col = FALSE,
  width = 6,
  height = 6.5,
  fontsize = 6
)
