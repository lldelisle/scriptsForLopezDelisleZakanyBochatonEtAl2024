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
# To modify the p-value font size:
safelyLoadAPackageInCRANorBioconductor("gginnards")
# Get good fontfamily
safelyLoadAPackageInCRANorBioconductor("extrafont")
# Only once
# font_import()
# loadfonts(device = "all")
fontFamily <- "DejaVu Sans"

my.geno.order <- c("r434 r434 d13hd d13hd", "r434 r434 d13 d13", "res res",
                   "b13hd b13hd", "res wt", "b13hd wt", "w453 w453 d13hd d13hd",
                   "w453 w453 d13 d13", "wt wt", "gof wt",
                   "gof gof", "gof gof D D", "gof gof D(I) D(I)")

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
full.df <- read.delim("tail_analysis/all_vertebrae_summary.txt")
# I use the geno order:
full.df$geno <- factor(full.df$geno, levels = my.geno.order)
table(full.df$geno, full.df$technique, exclude = NULL)
full.df$pretty.geno <- factor(my.labeller[as.character(full.df$geno)], levels = my.labeller)
table(full.df$pretty.geno, full.df$technique, exclude = NULL)
# I also use the simplified geno
full.df$pooled.geno <- factor(full.df$pooled.geno, levels = intersect(my.geno.order, full.df$pooled.geno))
table(full.df$pooled.geno, full.df$technique, exclude = NULL)
full.df$pretty.pooled.geno <- my.labeller[as.character(full.df$pooled.geno)]
full.df$pretty.pooled.geno <- factor(full.df$pretty.pooled.geno, levels = intersect(my.labeller, full.df$pretty.pooled.geno))
print(with(full.df, table(pretty.pooled.geno, technique, exclude = NULL)))
#                                                       technique
# pretty.pooled.geno                                     alizarin mCT
#   italic('HoxB'^{'Del(i9-13):hd'}~':Hoxd13'^{'hd'})           0   4
#   italic('HoxB'^{'Del(i9-13):hd'})                           11   6
#   italic('Hoxb13'^{'hd'})                                     3   8
#   italic('HoxB'^{'Del(i9-13):hd/+'})                         14   0
#   italic('Hoxb13'^{'hd/+'})                                   4   0
#   italic('Hoxd13'^{'hd'})                                     0   4
#   wt                                                         28  12
#   italic('HoxB'^{'Del(i9-13)/+'})                            14   0
#   italic('HoxB'^{'Del(i9-13)'})                              17   5
#   italic('HoxB'^{'Del(i9-13)'}~':HoxD'^{'Del(10-12)'})        0   3
my.comparisons.p <- list(c("gof gof", "gof gof D(I) D(I)"),
                         c("wt wt", "gof gof"),
                         c("wt wt", "gof wt"),
                         c("res res", "r434 r434 d13hd d13hd"),
                         c("b13hd b13hd", "b13hd wt"),
                         c("b13hd b13hd", "wt wt"),
                         c("wt wt", "b13hd wt"),
                         c("wt wt", "res wt"))

my.comparisons.pretty <- lapply(my.comparisons.p, function(v) {
    my.labeller[v]
})
my.technique.pretty <- c('alizarin' = "'Alizarin redS'", 'mCT' = "mu*'CT'")
full.df$pretty.technique <- my.technique.pretty[full.df$technique]
g <- ggboxplot(
  full.df[, c("pretty.pooled.geno", "Ca", "pretty.technique")],
  x = "pretty.pooled.geno",
  y = "Ca",
  add = "jitter",
  add.params = list(size = 1),
  color = "pretty.technique",
  # size = 0.5,
  legend = "right",
  orientation = "horizontal"
) +
  stat_compare_means(
    method = "t.test", comparisons = my.comparisons.pretty,
    aes(label = ..p.signif..),
    label.y = c(35, 35, 34, 35, 34, 36, 34, 35),
    family = fontFamily
  ) + # Use stars instead of p-values
  theme(
    text = element_text(family = fontFamily, size = 6),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey"),
  ) +
  ylab("Number of caudal vertebrae") +
  scale_y_continuous(limits = c(22, 37), breaks = seq(22, 34, 2)) +
  scale_x_discrete(labels = parse_format()) +
  scale_color_discrete("Technique", labels = parse_format())

# Decrease the size of the p-value
g$layers[[which_layers(g, "GeomSignif")]]$aes_params$textsize <- 2
# Change font Family
g$layers[[which_layers(g, "GeomSignif")]]$aes_params$family <- fontFamily

ggsave(
  "tail_analysis/Nb_Caudal.pdf", g,
  width = 17.5, height = 7, unit = "cm"
)

g1 <- ggplot(full.df, aes(x = pretty.pooled.geno, fill = Ta)) +
  geom_bar() +
  theme_pubr(legend = "right") +
  labs(fill = "Anomalies in\nthoracic vertebrae", y = "Number of animals") +
  theme(
    text = element_text(family = fontFamily, size = 6),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 6)
  )
g2 <- ggplot(full.df, aes(x = pretty.pooled.geno, fill = paste0("L", L))) +
  geom_bar() +
  theme_pubr(legend = "right") +
  labs(fill = "Number of\nlumbar vertebrae", y = "Number of animals") +
  theme(
    text = element_text(family = fontFamily, size = 6),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.title.x = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 6)
  ) +
  scale_x_discrete(labels = parse_format())

# At this point, they have the same axis, but the axis lengths are unequal, so ...

# build the plots
p1.common.x <- ggplot_gtable(ggplot_build(g1))
p2.common.x <- ggplot_gtable(ggplot_build(g2))

# copy the plot height from p2 to p1
p1.common.x$widths <- p2.common.x$widths

# Put them next to the other:
ggarrange(p1.common.x, p2.common.x, ncol = 1, heights = c(1.2, 2), labels = c("A", "B"))

ggsave("tail_analysis/Homeotic.pdf", width = 4, height = 4.5)

full.df$formula <- with(full.df, paste0("C", C, ", T", T, ", L", L, ", S", S, ", C", Ca))
table(subset(full.df, pretty.pooled.geno == "wt")$formula)
# C7, T13, L5, S3, C29 C7, T13, L5, S4, C25 C7, T13, L5, S4, C28 
#                    1                    1                   12 
# C7, T13, L5, S4, C29 C7, T13, L5, S4, C30 C7, T13, L5, S4, C32 
#                   15                    2                    1 
# C7, T13, L6, S3, C29 C7, T13, L6, S4, C27 C7, T13, L6, S4, C28 
#                    1                    2                    4 
# C7, T13, L6, S4, C30 
#                    1 
table(subset(full.df, pooled.geno == "gof gof")$formula)
# C7, T12, L4, S4, C25 C7, T12, L4, S5, C25 C7, T12, L5, S4, C24 
#                    1                    2                    3 
# C7, T12, L5, S4, C25 C7, T12, L5, S4, C26 C7, T12, L5, S4, C27 
#                    8                    6                    2 