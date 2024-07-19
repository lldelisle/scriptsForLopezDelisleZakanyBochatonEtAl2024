# Set working directory (github directory)
setwd("~/Documents/mygit/scriptsForLopezDelisleZakanyBochatonEtAl2024/")
if (!require(devtools)) {
    install.packages("devtools")
    library(devtools)
}
install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("reshape2")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
# To modify the p-value font size:
safelyLoadAPackageInCRANorBioconductor("gginnards")
safelyLoadAPackageInCRANorBioconductor("ggrepel")
safelyLoadAPackageInCRANorBioconductor("ggrastr")
safelyLoadAPackageInCRANorBioconductor("ggh4x")
# Get good fontfamily
safelyLoadAPackageInCRANorBioconductor("extrafont")
# Only once
# font_import()
# loadfonts(device = "all")

# Parameters
samplesPlan <- "RNAseq/samplesplan_Gastruloids.txt"
tableWithNormalizedExpression <- "RNAseq/outputs/AllCufflinks_Simplified.txt.gz"
gtfFile <- "publishedData/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz"
tableWithDESeq2Summary <- "RNAseq/outputs/summary_selected.txt"
outputFolder <- "RNAseq/outputs/"
fixedColors <- list(genoclone = c("red", "green4", "#cd4c25", "#90e70f"))
names(fixedColors$genoclone) <- c("wt", "delBdeli9-13", "delBins5-10", "delBins5-10del5-10")
fontFamily <- "DejaVu Sans"

# Create outputFolder
dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)

# Read samples plan
samplesPlanDF <- read.delim(samplesPlan, check.names = FALSE)
samplesPlanDF$pretty.geno <- factor(
    samplesPlanDF$pretty.geno,
    levels = unique(samplesPlanDF$pretty.geno)
)
rownames(samplesPlanDF) <- samplesPlanDF$sample

# Read summary table
summaryDF <- read.delim(tableWithDESeq2Summary, check.names = FALSE)
summaryDF$subsetting <- gsub(
    "delBdeli9-13ordelBins5-10del5-10_", "",
    gsub("delBdeli9-13orwt_", "", summaryDF$subsetting)
)
summaryDF$time <- as.numeric(gsub("_$", "", summaryDF$subsetting))
summaryDF$stars <- cut(
    summaryDF$padj,
    breaks = c(0, 1e-3, 1e-2, 0.05, 1),
    labels = c("***", "**", "*", "ns")
)
summaryDF$stars[is.na(summaryDF$padj)] <- "ns"
summaryDF$starsNA <- cut(
    summaryDF$pvalue,
    breaks = c(0, 1e-3, 1e-2, 0.05, 1),
    labels = c("***", "**", "*", "ns")
)
summaryDF$starsNA[is.na(summaryDF$pvalue)] <- "ns"

# Read FPKM table
expressionDF <- read.delim(tableWithNormalizedExpression, check.names = FALSE)
metaCols <- which(sapply(colnames(expressionDF), function(cn) {
    class(expressionDF[, cn]) != "numeric"
}))
# Remove FPKM_ from column names:
colnames(expressionDF) <- gsub("^FPKM_", "", colnames(expressionDF))
# Asign rownames
rownames(expressionDF) <- expressionDF$gene_id
samplesToPlot <- intersect(samplesPlanDF$sample, colnames(expressionDF))
# Get only values (no meta)
data <- expressionDF[, samplesToPlot]
# Also subset samples plan if needed
samplesPlanDF <- samplesPlanDF[samplesToPlot, ]

samplesPlanDF <- samplesPlanDF %>%
    group_by(genoclone, time) %>%
    mutate(n = n()) %>%
    as.data.frame()

ldata <- log2(data + 1)
# Plot Hoxb1-9-13 expression

my.gene.names <- c("Hoxb1", "Hoxb9", "Hoxb13")
names(my.gene.names) <- expressionDF$gene_id[match(my.gene.names, expressionDF$gene_short_name)]

ldata.selected <- melt(
    t(ldata[names(my.gene.names), ])
)
colnames(ldata.selected) <- c("sample", "Ens_ID", "log2FPKMp1")


temp.df <- merge(
    samplesPlanDF,
    ldata.selected
)

temp.df$gene_name <- paste0("italic('", my.gene.names[temp.df$Ens_ID], "')")
temp.df$stage <- paste0("'", temp.df$time, "h'")
# First, only plot wt, delBdeli9-13, delBins5-10 at 120h and 144h
my.genoclones <- c("wt", "delBdeli9-13", "delBins5-10")
selected.expression.df <- subset(
    temp.df, genoclone %in% my.genoclones & time %in% c(120, 144)
)
genoclone.pretty <- unique(selected.expression.df[, c("genoclone", "pretty.geno")])
genoclone.pretty.conv <- genoclone.pretty[
    match(my.genoclones, genoclone.pretty$genoclone), "pretty.geno"
]
names(genoclone.pretty.conv) <- my.genoclones
selected.expression.df$pretty.geno <- factor(
    selected.expression.df$pretty.geno,
    levels = genoclone.pretty.conv
)
selected.de <- subset(
    summaryDF, value %in% my.genoclones & ref.value %in% my.genoclones
)
selected.de$group1 <- factor(
    genoclone.pretty.conv[selected.de$value],
    levels = genoclone.pretty.conv
)
selected.de$group2 <- factor(
    genoclone.pretty.conv[selected.de$ref.value],
    levels = genoclone.pretty.conv
)
selected.de$gene_name <- paste0("italic('", selected.de$gene_name, "')")
selected.de$stage <- paste0("'", selected.de$time, "h'")

my.fixed.colors <- fixedColors$genoclone
names(my.fixed.colors) <- genoclone.pretty.conv[names(my.fixed.colors)]

# Get the good order for gene names
selected.expression.df$gene_name <-
    factor(
        selected.expression.df$gene_name,
        levels = paste0("italic('", my.gene.names, "')")
    )
selected.de$gene_name <- factor(
    selected.de$gene_name,
    levels = levels(selected.expression.df$gene_name)
)

g <- ggplot(
    selected.expression.df,
    aes(
        x = pretty.geno,
        y = log2FPKMp1,
        color = pretty.geno
    )
) +
    geom_jitter(width = 0.2, height = 0, size = 1) +
    geom_boxplot(
        data = subset(selected.expression.df, n >= 3),
        fill = NA
    ) +
    expand_limits(y = 0) +
    facet_nested(
        . ~ gene_name + stage,
        scale = "free_y",
        label = label_parsed
    ) +
    stat_pvalue_manual(
        data = selected.de,
        label = "starsNA",
        y.position = c(6.5, 7, 6.5, 8, 7.5, 7.5, 6.5, 6, 6),
        family = fontFamily
    ) +
    labs(
        color = "",
        y = expression(atop("gene"~"expression","["~log[2](1+FPKM)~"]"))) +
    theme_pubr() +
    theme(
        text = element_text(family = fontFamily, size = 6),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), # Put the genotypes vertical
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(3, "mm")
    ) +
    scale_color_manual("", values = my.fixed.colors)
# Decrease the size of the p-value
# g$layers[[which_layers(g, "GeomSignif")]]$aes_params$textsize <- 3
g$layers[[4]]$aes_params$label.size <- 2
set.seed(4)
ggsave("RNAseq/outputs/FPKM_1.pdf", g, width = 5, height = 3)


# Second, only plot delBdeli9-13, delBins5-10, delBins5-10del5-10 at 120h only Hoxb13
my.genoclones <- c("delBdeli9-13", "delBins5-10", "delBins5-10del5-10")
selected.expression.df <- subset(
    temp.df, genoclone %in% my.genoclones & time %in% c(120) &
    Ens_ID %in% names(my.gene.names[my.gene.names == "Hoxb13"])
)
genoclone.pretty <- unique(selected.expression.df[, c("genoclone", "pretty.geno")])
genoclone.pretty.conv <- genoclone.pretty[
    match(my.genoclones, genoclone.pretty$genoclone), "pretty.geno"
]
names(genoclone.pretty.conv) <- my.genoclones
selected.expression.df$pretty.geno <- factor(
    selected.expression.df$pretty.geno,
    levels = genoclone.pretty.conv
)
selected.de <- subset(
    summaryDF, value %in% my.genoclones & ref.value %in% my.genoclones &
    Ens_ID %in% selected.expression.df$Ens_ID
)
selected.de$group1 <- factor(
    genoclone.pretty.conv[selected.de$value],
    levels = genoclone.pretty.conv
)
selected.de$group2 <- factor(
    genoclone.pretty.conv[selected.de$ref.value],
    levels = genoclone.pretty.conv
)

my.fixed.colors <- fixedColors$genoclone
names(my.fixed.colors) <- genoclone.pretty.conv[names(my.fixed.colors)]
g <- ggplot(
    selected.expression.df,
    aes(
        x = pretty.geno,
        y = log2FPKMp1,
        color = pretty.geno
    )
) +
    geom_jitter(width = 0.2, height = 0) +
    geom_boxplot(
        data = subset(selected.expression.df, n >= 3),
        fill = NA
    ) +
    expand_limits(y = 0) +
    stat_pvalue_manual(
        data = selected.de,
        label = "starsNA",
        y.position = c(6, 6.5, 6),
        family = fontFamily,
        size = 2
    ) +
    labs(
        color = "",
        y = expression(atop(italic("Hoxb13")~"expression","["~log[2](1+FPKM)~"]"))) +
    theme_pubr() +
    theme(
        text = element_text(family = fontFamily, size = 6),
        legend.position = "right",
        legend.text = element_text(size = 6),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), # Put the genotypes vertical
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    scale_color_manual("", values = my.fixed.colors) + 
    guides(colour = guide_legend(nrow = 3))
set.seed(6)
ggsave("RNAseq/outputs/FPKM_2.pdf", g, width = 4.5, height = 3)

# Third, only those with pretty.del only Hoxb13
my.genoclones <- c("delBins5-10", "delBins5-10del8-10", "delBins5-10del7-10", "delBdeli9-13")
selected.expression.df <- subset(
    temp.df, genoclone %in% my.genoclones &
    Ens_ID %in% names(my.gene.names[my.gene.names == "Hoxb13"]) &
    time == 120
)
del.new.pretty <- unique(selected.expression.df[, c("genoclone", "del.new.name")])
del.new.pretty.conv <- del.new.pretty[
    match(my.genoclones, del.new.pretty$genoclone), "del.new.name"
]
names(del.new.pretty.conv) <- my.genoclones
selected.expression.df$del.new.name <- factor(
    selected.expression.df$del.new.name,
    levels = del.new.pretty.conv
)
selected.de <- subset(
    summaryDF, value %in% my.genoclones & ref.value %in% my.genoclones &
    Ens_ID %in% selected.expression.df$Ens_ID
)
selected.de$group1 <- factor(
    del.new.pretty.conv[selected.de$value],
    levels = del.new.pretty.conv
)
selected.de$group2 <- factor(
    del.new.pretty.conv[selected.de$ref.value],
    levels = del.new.pretty.conv
)

# my.fixed.colors <- fixedColors$genoclone
# names(my.fixed.colors) <- genoclone.pretty.conv[names(my.fixed.colors)]
g <- ggplot(
    selected.expression.df,
    aes(
        x = del.new.name,
        y = log2FPKMp1,
        color = del.new.name
    )
) +
    geom_jitter(width = 0.05, height = 0) +
    geom_boxplot(
        data = subset(selected.expression.df, n >= 3),
        fill = NA
    ) +
    geom_text(
        data = subset(selected.de, ref.value == "delBins5-10"),
        aes(label = starsNA, x = group1),
        y = 5.5, color = "black", family = fontFamily,
        size = 2
    ) +
    ylim(0, 5.5) +
    # facet_grid(. ~ time, scale = "free_x", space = "free_x",
    #            labeller = labeller(time = function(x){paste0(x, "h")})) +
    labs(
        color = "",
        y = expression(atop(italic("Hoxb13")~"expression","["~log[2](1+FPKM)~"]"))) +
    theme_pubr() +
    theme(
        text = element_text(family = fontFamily, size = 6),
        axis.text.x = element_text(
            angle = 20, vjust = 0.5, hjust = 0.5
        ), # Put the genotypes vertical
        axis.title.x = element_blank(),
        legend.position = "none"
    ) # +
    # scale_color_manual("", values = my.fixed.colors)
set.seed(6)
ggsave(
    "RNAseq/outputs/FPKM_3.pdf", g,
    width = 7.5, height = 5,
    unit = "cm"
)

# Fourth, inversion effect
my.genoclones <- c("delBins5-10", "delBinv5-10clone1", "delBinv5-10clone2", "delBins5-10del7-10", "delBinv5-10del7-10")
selected.expression.df <- subset(
    temp.df, genoclone %in% my.genoclones &
        Ens_ID %in% names(my.gene.names[my.gene.names == "Hoxb13"])
)
# Only get conditions where both orientations are present
unique.cond <- unique(selected.expression.df[, c("time", "del.new.name", "orientation")])
unique.cond.dup <- unique.cond[duplicated(unique.cond[, c("time", "del.new.name")]), c("time", "del.new.name")]

selected.expression.df <- merge(selected.expression.df, unique.cond.dup)
del.new.pretty <- unique(selected.expression.df[, c("genoclone", "del.new.name")])
del.new.pretty.conv <- del.new.pretty[
    match(my.genoclones, del.new.pretty$genoclone), "del.new.name"
]
names(del.new.pretty.conv) <- my.genoclones
selected.de <- subset(
    summaryDF, value %in% my.genoclones & ref.value %in% my.genoclones &
    Ens_ID %in% selected.expression.df$Ens_ID
)
selected.de$del.new.name1 <- factor(
    del.new.pretty.conv[selected.de$value],
    levels = unique(del.new.pretty.conv)
)
selected.de$del.new.name2 <- factor(
    del.new.pretty.conv[selected.de$ref.value],
    levels = unique(del.new.pretty.conv)
)
selected.de <- subset(
    selected.de,
    del.new.name1 == del.new.name2
)
selected.de$del.new.name <- selected.de$del.new.name1
orientation.new <- unique(selected.expression.df[, c("genoclone", "orientation.new")])
orientation.new.conv <- orientation.new[
    match(my.genoclones, orientation.new$genoclone), "orientation.new"
]
names(orientation.new.conv) <- my.genoclones
selected.de$group1 <- factor(
    orientation.new.conv[selected.de$value],
    levels = unique(orientation.new.conv)
)
selected.de$group2 <- factor(
    orientation.new.conv[selected.de$ref],
    levels = unique(orientation.new.conv)
)
selected.expression.df$orientation.new <- factor(
    selected.expression.df$orientation.new,
    levels = unique(orientation.new.conv)
)
selected.expression.df$orientation <- factor(
    selected.expression.df$orientation,
    levels = unique(orientation.new.conv)
)
# my.fixed.colors <- fixedColors$genoclone
# names(my.fixed.colors) <- genoclone.pretty.conv[names(my.fixed.colors)]
g <- ggplot(
    selected.expression.df,
    aes(
        x = orientation.new,
        y = log2FPKMp1,
        color = orientation,
        group = genoclone
    )
) +
    geom_jitter(width = 0.2, height = 0) +
    geom_boxplot(
        data = subset(selected.expression.df, n >= 3),
        fill = NA
    ) +
    stat_pvalue_manual(
        data = selected.de,
        label = "starsNA",
        y.position = c(5, 4.5),
        family = fontFamily,
        size = 2
    ) +
    ylim(0, 5.5) +
    facet_nested(. ~ time + del.new.name, scale = "free_x", space = "free_x",
               labeller = labeller(time = function(x){paste0(x, "h")})) +
    labs(
        color = "",
        y = expression(atop(italic("Hoxb13")~"expression","["~log[2](1+FPKM)~"]"))) +
    theme_pubr() +
    theme(
        text = element_text(family = fontFamily, size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), # Put the genotypes vertical
        axis.title.x = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 6)
    # legend.box.spacing = unit(1, "pt")
    ) # +
    # scale_color_manual("", values = my.fixed.colors)
set.seed(4)
ggsave("RNAseq/outputs/FPKM_4.pdf", g, width = 3.8, height = 3)


# Plot volcano delBins5-10 vs delBdeli9-13
full.de <- read.delim("RNAseq/outputs/summary_long.txt.gz", check.names = FALSE)
my.genoclones <- c("delBins5-10", "delBdeli9-13")
current.de <- subset(full.de, value %in% my.genoclones & ref.value %in% my.genoclones)
log2FC.threshold <- log2(1.5)
signif.genes.df <- subset(
    current.de,
    abs(log2FoldChange) > log2FC.threshold & !is.na(padj) & padj < 0.05
)
g <- ggplot(current.de, aes(x = -log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_point(
        data = signif.genes.df,
        color = "blue"
    ) +
    geom_hline(yintercept = -log10(0.05), lty = 2, linewidth = .2) +
    geom_vline(xintercept = log2FC.threshold * c(-1, 1), lty = 2, linewidth = .2) +
    geom_label_repel(
        data = signif.genes.df,
        aes(label = gene_name),
        family = fontFamily
    ) +
    theme_pubr() +
    xlab("log2FC") +
    ggtitle(paste0(unique(current.de$ref.value), "vs", unique(current.de$value))) +
    theme(text = element_text(family = fontFamily))
ggsave("RNAseq/outputs/volcano_1.pdf", rasterize(g, layers = "Point", dpi = 250), width = 5, height = 5)

# Plot volcano delBdeli9-13 vs wt
my.genoclones <- c("wt", "delBdeli9-13")
current.de <- subset(full.de, value %in% my.genoclones & ref.value %in% my.genoclones)
log2FC.threshold <- log2(1.5)
signif.genes.df <- subset(
    current.de,
    abs(log2FoldChange) > log2FC.threshold & !is.na(padj) & padj < 0.05
)
g <- ggplot(current.de, aes(x = -log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_point(
        data = signif.genes.df,
        color = "blue"
    ) +
    geom_hline(yintercept = -log10(0.05), lty = 2, linewidth = .2) +
    geom_vline(xintercept = log2FC.threshold * c(-1, 1), lty = 2, linewidth = .2) +
    geom_label_repel(
        data = signif.genes.df,
        aes(label = gene_name),
        family = fontFamily
    ) +
    theme_pubr() +
    xlab("log2FC") +
    ggtitle(paste0(unique(current.de$ref.value), "vs", unique(current.de$value))) +
    theme(text = element_text(family = fontFamily))
ggsave("RNAseq/outputs/volcano_2.pdf", rasterize(g, layers = "Point", dpi = 250), width = 5, height = 5)
