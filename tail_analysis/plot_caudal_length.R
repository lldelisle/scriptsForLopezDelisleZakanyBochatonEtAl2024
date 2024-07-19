# Set working directory (github directory)
setwd("~/Documents/mygit/scriptsForLopezDelisleZakanyBochatonEtAl2024/")
# Install the necessary packages:
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("scales")
# Get good fontfamily
safelyLoadAPackageInCRANorBioconductor("extrafont")
# Only once
# font_import()
# loadfonts(device = "all")
fontFamily <- "DejaVu Sans"

my.geno.order <- c(
  "r434 r434 d13hd d13hd", "r434 r434 d13 d13", "res res",
  "b13hd b13hd", "w453 w453 d13hd d13hd",
  "w453 w453 d13 d13", "wt wt",
  "gof gof", "gof gof D D", "gof gof D(I) D(I)"
)

my.labeller <- c(
    "r434 r434 d13hd d13hd" = "italic('HoxB'^{'Del(i9-13):Hoxb13hd'}~':Hoxd13'^{'hd'})",
    "r434 r434 d13 d13" = "italic('HoxB'^{'Del(i9-13):Hoxb13hd'}~':Hoxd13'^{'+/+'})",
    "res res" = "italic('HoxB'^{'Del(i9-13):Hoxb13hd'})",
    "b13hd b13hd" = "italic('Hoxb13'^{'hd'})",
    "res wt" = "italic('HoxB'^{'Del(i9-13):Hoxb13hd/+'})",
    "b13hd wt" = "italic('Hoxb13'^{'hd/+'})",
    "w453 w453 d13hd d13hd" = "italic('Hoxd13'^{'hd'})",
    "w453 w453 d13 d13" = "italic('Hoxd13'^{'+/+'})",
    "wt wt" = "wt",
    "gof wt" = "italic('HoxB'^{'Del(i9-13)/+'})",
    "gof gof" = "italic('HoxB'^{'Del(i9-13)'})",
    "gof gof D D" = "italic('HoxB'^{'Del(i9-13)'}~':HoxD'^{'+/+'})",
    "gof gof D(I) D(I)" = "italic('HoxB'^{'Del(i9-13)'}~':HoxD'^{'Del(10-12)'})"
)

my.colors <- c(
  "r434 r434 d13hd d13hd" = "#030327",
  "res res" = "royalblue4",
  "b13hd b13hd" = "#76dae3",
  "w453 w453 d13hd d13hd" = "#38666b",
  "wt wt" = "red",
  "gof gof" = "green4",
  "gof gof D(I) D(I)" = "green"
)

my.colors.pretty <- my.colors
names(my.colors.pretty) <- my.labeller[names(my.colors)]
full.df <- read.delim("tail_analysis/all_caudal_measurements_mCT.txt")
# I use the geno order:
full.df$geno <- factor(full.df$geno, levels = my.geno.order)
full.df$pretty.geno <- factor(my.labeller[as.character(full.df$geno)], levels = my.labeller)
# I also use the simplified geno
full.df$pooled.geno <- factor(full.df$pooled.geno, levels = intersect(my.geno.order, full.df$pooled.geno))
full.df$pretty.pooled.geno <- my.labeller[as.character(full.df$pooled.geno)]
full.df$pretty.pooled.geno <- factor(full.df$pretty.pooled.geno, levels = intersect(my.labeller, full.df$pretty.pooled.geno))

ggboxplot(
  full.df,
  x = "caudal.id", y = "length",
  color = "pretty.pooled.geno",
  outlier.size = 0.1
) +
  theme(
    text = element_text(family = fontFamily, size = 6),
    legend.position.inside = c(0.25, 0.25),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 6)
  ) +
  xlab("Position of caudal vertebrae") +
  ylab("Vertebrae length [mm]") +
  scale_colour_manual(
    "",
    values = my.colors.pretty,
    labels = parse_format()
  ) +
  guides(
    colour = guide_legend(position = "inside")
  )

ggsave(
  "tail_analysis/caudal_length_boxplot.pdf",
  width = 17.5, height = 8.5, unit = "cm"
)
