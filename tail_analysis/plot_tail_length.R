# Set working directory (github directory)
setwd("~/Documents/mygit/scriptsForLopezDelisleZakanyBochatonEtAl2024/")
# Install the necessary packages:
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
# To reshape a df:
safelyLoadAPackageInCRANorBioconductor("reshape2")
# To make pretty plots:
safelyLoadAPackageInCRANorBioconductor("ggplot2")
# To have p-values on plots:
safelyLoadAPackageInCRANorBioconductor("ggpubr")
# To modify the p-value font size:
safelyLoadAPackageInCRANorBioconductor("gginnards")
# To display equations of regression:
safelyLoadAPackageInCRANorBioconductor("ggpmisc")
# To get exponent in legend:
safelyLoadAPackageInCRANorBioconductor("scales")
# Get good fontfamily
safelyLoadAPackageInCRANorBioconductor("extrafont")
# Only once
# font_import()
# loadfonts(device = "all")
fontFamily <- "DejaVu Sans"


my.new.order <- c(
  "gof/gof", "gof/wt", "wt/wt",
  "r405/wt", "r405/r405", "r434/wt", "r434/r434",
  "m438;w453", "m438", "m424;w453", "m424"
)
my.pretty.names <- c("Del(i9-13)", "wt", "Del(i9-13):hd1", "Del(i9-13):hd2", "italic('Hoxb13'^{'hd1'})", "italic('Hoxb13'^{'hd2'})")
names(my.pretty.names) <- c("gof", "wt", "r405", "r434", "m438", "m424")
my.line.colors <- c("green4", "red", "grey50", "royalblue4", "cadetblue", "cadetblue")
names(my.line.colors) <- my.pretty.names


# Data for tail length was sumarized into 2 data tables:
df <- read.delim("tail_analysis/tail_length_per_geno.txt")
# List all comparisons to display
# Including hetero vs homo
my_comparisons_full <- c(
  lapply(intersect(setdiff(my.new.order, c("wt/wt")), df$geno), function(s) {
    c("wt/wt", s)
  }),
  list(
    c("r405/wt", "r405/r405"), c("r434/wt", "r434/r434"),
    c("gof/wt", "gof/gof")
  )
)
# I need to check the p-values
sapply(my_comparisons_full, function(v) {
  t.test(df$length[df$geno == v[1]], df$length[df$geno == v[2]])$p.value
})
# Splitting by line and status
df$line <- sapply(as.character(df$geno), function(s) {
  strsplit(s, "\\/")[[1]][1]
})
# Use pretty.names
df$line <- my.pretty.names[df$line]
df$status <- "hom"
df$status[grepl("wt", as.character(df$geno))] <- "het"
df$status[df$geno == "wt/wt"] <- "hom"
df$line <- factor(df$line, levels = names(my.line.colors))
# Try to get similar data from hoxb13hd
df2 <- read.delim("tail_analysis/tail_length_per_id_age_line_13hd.txt", check.names = FALSE)
# Get the last measure for each animal if above 60 days
df2 <- df2[order(-df2$age), ]
df2last <- df2[!duplicated(df2$ID) & df2$age >= 60, ]

my_comparisons_full2b <- c(
  lapply(intersect(setdiff(my.new.order, c("wt")), df2last$geno), function(s) {
    c("wt", s)
  }),
  list(c("m424;w453", "m424"), c("m438;w453", "m438"))
)
# I need to check the p-values
sapply(my_comparisons_full2b, function(v) {
  t.test(df2last$length[df2last$geno == v[1]], df2last$length[df2last$geno == v[2]])$p.value
})
# Splitting by line and status
df2last$line <- sapply(as.character(df2last$geno), function(s) {
  strsplit(s, ";")[[1]][1]
})
# Use pretty.names
df2last$line <- my.pretty.names[df2last$line]
df2last$status <- "hom"
df2last$status[grepl("w453", as.character(df2last$geno))] <- "het"
df2last$line <- factor(df2last$line, levels = names(my.line.colors))

# Put all together
# I mix wt together
df.both <- rbind(df, df2last[, colnames(df)])
df.both[df.both$geno == "wt", "geno"] <- "wt/wt"

my_comparisons_full2b_comb <- c(
  lapply(intersect(setdiff(my.new.order, c("wt")), df2last$geno), function(s) {
    c("wt/wt", s)
  }),
  list(c("m438;w453", "m438"), c("m424;w453", "m424"))
)
my_comparisons_both <- c(my_comparisons_full, my_comparisons_full2b_comb)
my.line.colors.parsable <- my.line.colors
names(my.line.colors.parsable)[grepl("Del\\(", names(my.line.colors.parsable))] <-
  paste0("'", names(my.line.colors.parsable)[grepl("Del\\(", names(my.line.colors.parsable))], "'")
df.both$line <- as.character(df.both$line)
df.both$line[grepl("Del\\(", df.both$line)] <-
  paste0("'", df.both$line[grepl("Del\\(", df.both$line)], "'")
print(with(df.both, table(line, status)))
#                           status
# line                       het hom
#   'Del(i9-13)'              22   5
#   'Del(i9-13):hd1'          17  14
#   'Del(i9-13):hd2'          28  13
#   italic('Hoxb13'^{'hd1'})  20   9
#   italic('Hoxb13'^{'hd2'})  22  11
#   wt                         0  60

line.order <- unique(df.both$line[match(my.new.order, df.both$geno)])
df.both$line <- factor(df.both$line, levels = line.order)

line.status <- unique(df.both[, c("geno", "status")])$status
names(line.status) <- unique(df.both[, c("geno", "status")])$geno

new.g <- ggboxplot(df.both,
  x = "geno", y = "length",
  color = "line", alpha = "status",
  outlier.shape = NA, # This is to prevent the boxplot to show outlayers (we will see them in the jitter)
  order = my.new.order
  # linewidth = 0.2 # This does not work...
) +
  geom_jitter(aes(color = line, alpha = status), size = 0.4, stroke = 0) + # The jitter need to be separately because inside ggboxplot the alpha is not used
  stat_compare_means(
    method = "t.test",
    comparisons = my_comparisons_both,
    label.y = c(
      c(
        100, 97, 103, 106, # This is manual arrangement
        109, 112, 103, 109, 97
      ),
      c(
        115, 118, 121, 124, # This is manual arrangement
        115, 121
      )
    ),
    label = "p.signif",
    family = fontFamily, # I don't know why it does not work
    tip.length = 0.02,
    vjust = 0
  ) + # I decreased the size of the bracket
  scale_color_manual("", values = my.line.colors.parsable, labels = parse_format()) + # Use my colors
  scale_alpha_discrete("", range = c(0.5, 1)) + # use exactly alpha of 0.5 and 1
  guides(alpha = "none") +
  ylab("Tail length [mm]") +
  scale_y_continuous(limits = c(70, 127), breaks = c(70, 80, 90, 100, 110)) +
  scale_x_discrete(labels = line.status[my.new.order]) +
  theme(
    text = element_text(family = fontFamily, size = 6),
    line = element_line(linewidth = 0.2),
    axis.line = element_line(linewidth = 0.2),
    axis.text.x = element_text(margin = margin(b = 20, t = -7), size = 5),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right", # Put the legend on right
    panel.grid.major.y = element_line(colour = "grey"),
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(margin = margin(l = 2)),
    legend.box.spacing = unit(1, "pt")
  ) # Add some major lines y axis

# Decrease the size of the p-value
new.g$layers[[which_layers(new.g, "GeomSignif")]]$aes_params$textsize <- 1.5
# Set fontFamily
new.g$layers[[which_layers(new.g, "GeomSignif")]]$aes_params$family <- fontFamily
# Set linewidth of boxplot
new.g$layers[[1]]$aes_params$linewidth <- 0.2
set.seed(2)
ggsave("tail_analysis/FigTailLength_allboxplots.pdf",
  new.g,
  width = 7.5, height = 4.5,
  unit = "cm"
)

# Second figure = trend lines
# Get only homozygous
df2 <- subset(df2, !grepl("w453", geno))
# pretty line:
df2$line <- my.pretty.names[df2$geno]
df2$line <- factor(df2$line, levels = my.pretty.names)
print(table(df2$line))
#   Del(i9-13)                       wt             Del(i9-13)R1
#            0                      175                        0
# Del(i9-13)R2 italic('Hoxb13'^{'hd1'}) italic('Hoxb13'^{'hd2'})
#            0                      138                      132
print(table(unique(df2[, c("ID", "line")])$line))
#   Del(i9-13)                       wt             Del(i9-13)R1
#            0                       15                        0
# Del(i9-13)R2 italic('Hoxb13'^{'hd1'}) italic('Hoxb13'^{'hd2'})
#            0                        9                       11

# Split the growth phase in 3
df2$part <- 1
df2$part[df2$age >= 56] <- 2
df2$part[df2$age >= 84] <- 3
df2$part <- as.factor(df2$part)
df2$color <- my.line.colors[df2$line]
# Plot
my.formula <- y ~ x
g2 <- ggplot(df2, aes(x = age, y = length, color = color, linetype = part)) +
  geom_point(aes(shape = line)) +
  scale_color_manual("", values = unique(df2$color)) +
  scale_shape_discrete("", labels = parse_format()) +
  stat_poly_eq(
    formula = my.formula, # Put the equations
    aes(
      label = paste0("atop(", ..eq.label.., ",", ..rr.label.., ")") # Use atop to have one on top of the other
    ),
    parse = TRUE,
    # label.x = rep(c(0.05, 0.3, 0.9), 2), # Specify x as proportion
    # label.y = (c(105, 110, 115, 65, 70, 75) - 60) / 60, # Specify y as proportion
    label.x = rep(c(0.1, 0.5, 0.9), 2), # Specify x as proportion
    label.y = (c(105, 110, 117, 65, 70, 75) - 60) / 60, # Specify y as proportion
    family = fontFamily,
    size = 2
  ) + # Change the font size
  geom_smooth(method = "lm", se = FALSE, formula = my.formula) +
  ylim(60, 120) +
  guides(
    linetype = "none", color = "none",
    shape = guide_legend(
      override.aes = list(color = my.line.colors[sort(unique(df2$line))])
    )
  ) + # Remove the linetype and color from the legend
  theme_pubr(legend = "right") + # Use the same theme as figure 2A
  theme(
    text = element_text(family = fontFamily, size = 6),
    panel.grid.major.y = element_line(colour = "grey"), # Add the horizontal lines
    axis.title.x.bottom = element_text(margin = margin(b = 20, t = -20)),
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(margin = margin(l = 2), size = 6),
    legend.box.spacing = unit(1, "pt")
  ) +
  ylab("Tail length [mm]") +
  xlab("Age after birth [days]") +
  scale_x_continuous(breaks = c(28, 56, 84, 112))

ggsave("tail_analysis/FigTailLength_time_course.pdf", g2, width = 17.5, height = 8, unit = "cm")
