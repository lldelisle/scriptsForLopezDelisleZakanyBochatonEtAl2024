# Set working directory (github directory)
setwd("~/Documents/mygit/scriptsForLopezDelisleZakanyBochatonEtAl2024/")
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("eulerr")
safelyLoadAPackageInCRANorBioconductor("grid")
# Get good fontfamily
safelyLoadAPackageInCRANorBioconductor("extrafont")
# Only once
# font_import()
# loadfonts(device = "all")
fontFamily <- "DejaVu Sans"
# Parameters
tableWithDESeq2Summary <- "RNAseq/outputs/summary_E.txt.gz"
outputFolder <- "RNAseq/outputs/"

# Create outputFolder
dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)

# Read summary table
summaryDF <- read.delim(tableWithDESeq2Summary, check.names = FALSE)

# Identify comparisons:
comparison.names <- gsub(
  "_l2fc",
  "",
  grep("_l2fc$", colnames(summaryDF), value = TRUE)
)
pretty.names <- gsub("genotype_", "", comparison.names)
pretty.names <- gsub("deli9-13ordeli9-13hd_", "", pretty.names)
pretty.names <- gsub("deli9-13", "Del(i9-13)", pretty.names)
pretty.names <- gsub(")hd", "):hd", pretty.names)

signifs.list.up <- lapply(
  comparison.names,
  function(comparison.name) {
    which(
      summaryDF[, paste0(comparison.name, "_signif")] &
        summaryDF[, paste0(comparison.name, "_l2fc")] > 0
    )
  }
)
names(signifs.list.up) <- gsub("vs", ">", pretty.names)
signifs.list.down <- lapply(
  comparison.names,
  function(comparison.name) {
    which(
      summaryDF[, paste0(comparison.name, "_signif")] &
        summaryDF[, paste0(comparison.name, "_l2fc")] < 0
    )
  }
)
names(signifs.list.down) <- gsub("vs", "<", pretty.names)

signifs.list <- c(signifs.list.down, signifs.list.up)

full.venn.labels <- c(
 "E9.5_Del(i9-13)>wt",
 "E10.5_Del(i9-13)>wt",
 "E10.5_Del(i9-13):hd<Del(i9-13)",
 "E9.5_Del(i9-13)<wt",
 "E10.5_Del(i9-13)<wt",
 "E10.5_Del(i9-13):hd>Del(i9-13)"
)
pdf("RNAseq/outputs/euler_E.pdf", width = 6, height = 3, title = "euler_E")
p <- plot(
  euler(subsetByNamesOrIndices(signifs.list, full.venn.labels)),
  quantities = list(fontfamily = fontFamily, fontsize = 6),
  fills = "transparent",
  labels = FALSE,
)
print(p)
xlims <- p$data$xlim
ylims <- p$data$ylim
for (my.id in 1:6) {
  my.data <- p$data$ellipse[my.id, ]
  my.label <- gsub("_", "\n", rownames(my.data))
  my.x <- my.data$h
  if (my.data$k > 0) {
    my.y <- my.data$k + my.data$a / 1.3 + 0.3
  } else {
    my.y <- my.data$k - my.data$a / 1.3 - 0.3
  }
  # if (my.x > 0) {
  #   my.y <- my.y - 0.3
  # }
  grid.text(
    my.label,
    x = (my.x - xlims[1]) / (xlims[2] - xlims[1]),
    y = (my.y - ylims[1]) / (ylims[2] - ylims[1]),
    gp = gpar(col = "black", fontsize = 6, fontfamily = fontFamily)
  )
}
dev.off()

summaryDF$gene_name[Reduce(intersect, subsetByNamesOrIndices(signifs.list, full.venn.labels[1:3]))]
# [1] "Chl1"     "Baiap2l1"

summaryDF$gene_name[Reduce(intersect, subsetByNamesOrIndices(signifs.list, full.venn.labels[1:2]))]
# [1] "Chl1"     "Baiap2l1" "Hoxb13"

summaryDF$gene_name[Reduce(intersect, subsetByNamesOrIndices(signifs.list, full.venn.labels[4:6]))]
# [1] "Rbpj"

summaryDF$gene_name[Reduce(intersect, subsetByNamesOrIndices(signifs.list, full.venn.labels[4:5]))]
# [1] "Rbpj"  "Tor3a"
