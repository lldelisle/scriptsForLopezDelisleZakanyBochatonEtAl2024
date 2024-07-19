# Set working directory (github directory)
setwd("~/Documents/mygit/scriptsForLopezDelisleZakanyBochatonEtAl2024/")
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("biomaRt")
safelyLoadAPackageInCRANorBioconductor("taxizedb")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("scales")
# Get good fontfamily
safelyLoadAPackageInCRANorBioconductor("extrafont")
# Only once
# font_import()
# loadfonts(device = "all")
fontFamily <- "DejaVu Sans"

# All species (not mouse strains)
selected.ranks <- c("Hominidae", "Primates", "Muridae", "Eutheria", "Metatheria", "Aves", "Crocodylia", "Lepidosauria", "Testudines", "Latimeria chalumnae", "Actinopterygii", "Vertebrata")

pretty.names <- c("Hominidae", "Other Primates", "Muridae", "Other Placentalia", "Marsupialia", "Aves", "Reptilia", "Reptilia", "Reptilia", "Coelacanth", "Actinopterygii", "Other Vertebrata")
names(pretty.names) <- selected.ranks

if (!file.exists("b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13_diff_plot.txt")) {
  # Ensembl Release 111
  ensembl <- useMart("ensembl")
  all_datasets <- listDatasets(ensembl)
  if (!file.exists("b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13_diff.txt")) {
    if (!file.exists("b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13.txt")) {
      all.results <- as.data.frame(t(rep(NA, 8)))
      colnames(all.results) <- c(
        "ensembl_transcript_id", "chromosome_name", "exon_chrom_start",
        "exon_chrom_end", "external_gene_name", "strand", "rank",
        "mmusculus_homolog_associated_gene_name"
      )
      all.results$dataset <- NA
      # I needed to run twice the for loop because I got an error:
      # Error in curl::curl_fetch_memory(url, handle = handle) :
      #   Timeout was reached: [www.ensembl.org:443] Operation timed out after 300001 milliseconds with 8169 bytes received
      # I also got:
      # Ensembl site unresponsive, trying asia mirror
      # Error in curl::curl_fetch_memory(url, handle = handle) :
      #   SSL certificate problem: unable to get local issuer certificate

      while (!all(all_datasets$dataset %in% all.results$dataset)) {
        for (dataset in setdiff(all_datasets$dataset, all.results$dataset)) {
          print(dataset)
          # Get info for the current dataset
          current_ensembl <- tryCatch(useDataset(dataset, ensembl),
            error = function(e) {
              NULL
            }
          )
          if (is.null(current_ensembl)) {
            print("Failed to get the ensembl database")
            next
          }
          # First way = using external_gene_name
          results1 <- tryCatch(
            getBM(
              attributes = c(
                "ensembl_transcript_id", "chromosome_name",
                "exon_chrom_start", "exon_chrom_end",
                "external_gene_name", "strand", "rank"
              ),
              filters = "external_gene_name",
              values = c("Hoxb1", "Hoxb9", "Hoxb13", "HOXB1", "HOXB9", "HOXB13", "hoxb1a", "hoxb9a", "hoxb13a"),
              mart = current_ensembl
            ),
            error = function(e) {
              NULL
            }
          )
          if (is.null(results1)) {
            print("Failed to get the genes")
            next
          }
          results1$mmusculus_homolog_associated_gene_name <- gsub("HOXB", "Hoxb", gsub("a", "", results1$external_gene_name), ignore.case = T)
          # Second way = using homolog with mouse
          if ("mmusculus_homolog_associated_gene_name" %in% listAttributes(current_ensembl)$name) {
            cat("Getting orthologues...\n")
            all_homologous <- tryCatch(
              getBM(
                attributes = c(
                  "ensembl_gene_id", "ensembl_transcript_id",
                  "mmusculus_homolog_associated_gene_name"
                ),
                filters = "with_mmusculus_homolog",
                values = T,
                mart = current_ensembl
              ),
              error = function(e) {
                NULL
              }
            )
            if (is.null(all_homologous)) {
              print("Failed to get homologous with mouse")
              next
            }
            ids <- unique(all_homologous$ensembl_gene_id[all_homologous$mmusculus_homolog_associated_gene_name %in% c("Hoxb1", "Hoxb9", "Hoxb13")])
            if (length(ids) > 0) {
              cat("Getting exon info...\n")
              results2 <- tryCatch(
                getBM(
                  attributes = c(
                    "ensembl_transcript_id", "chromosome_name",
                    "exon_chrom_start", "exon_chrom_end",
                    "external_gene_name", "strand", "rank"
                  ),
                  filters = "ensembl_gene_id",
                  values = ids,
                  mart = current_ensembl
                ),
                error = function(e) {
                  NULL
                }
              )
              if (is.null(results1)) {
                print("Failed")
                next
              }
              results2 <- merge(results2, all_homologous[, c("ensembl_transcript_id", "mmusculus_homolog_associated_gene_name")])
              # Put the second approach within first results
              results1 <- unique(rbind(results1, results2))
            }
          }
          # Add this result to 'all.results'.
          if (nrow(results1) == 0) {
            results1 <- as.data.frame(t(rep(NA, 8)))
            colnames(results1) <- c(
              "ensembl_transcript_id", "chromosome_name", "exon_chrom_start",
              "exon_chrom_end", "external_gene_name", "strand", "rank",
              "mmusculus_homolog_associated_gene_name"
            )
          }
          results1$dataset <- dataset
          print(unique(results1$mmusculus_homolog_associated_gene_name))
          all.results <- rbind(all.results, results1)
        }
      }
      # Add info on the datasets:
      all.results <- merge(all.results, all_datasets, all = T)
      all.results <- unique(all.results)
      all.results <- all.results[order(all.results$dataset, all.results$external_gene_name, all.results$ensembl_transcript_id, all.results$rank), ]
      write.table(all.results,
        "b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13.txt",
        sep = "\t", row.names = F, quote = F
      )
    } else {
      all.results <- read.delim("b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13.txt")
    }
    summary_table <- NULL
    sink("b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13_diff.log")
    for (dataset in all_datasets$dataset) {
      current.result <- all.results[all.results$dataset %in% dataset, ]
      for (genes.to.check in c("Hoxb1_Hoxb9", "Hoxb9_Hoxb13")) {
        both.genes <- strsplit(genes.to.check, "_")[[1]]
        current.result.both.genes <- subset(current.result, mmusculus_homolog_associated_gene_name %in% both.genes)
        current.genes <- unique(current.result.both.genes[, c("chromosome_name", "external_gene_name")])
        if (nrow(current.genes) == 0) {
          results <- c(dataset, genes.to.check, NA, "no_homologous_gene_found")
          summary_table <- rbind(summary_table, results)
          next
        }
        if (nrow(current.genes) == 2 && length(unique(current.genes$chromosome_name)) != 1) {
          if (length(unique(current.genes$external_gene_name)) == 1) {
            results <- c(dataset, genes.to.check, NA, "only_one_gene_found")
          } else {
            results <- c(dataset, genes.to.check, NA, "two_genes_on_two_different_chr")
          }
          summary_table <- rbind(summary_table, results)
          next
        } else if (nrow(current.genes) == 1) {
          if (is.na(current.genes$chromosome_name)) {
            results <- c(dataset, genes.to.check, NA, "no_homologous_gene_found")
          } else {
            results <- c(dataset, genes.to.check, NA, "only_one_gene_found")
          }
          summary_table <- rbind(summary_table, results)
          next
        }
        if (nrow(current.genes) > 2) {
          chrom_nb <- table(current.genes$chromosome_name)
          probable_chrom <- names(chrom_nb[which.max(chrom_nb)])
          comment <- paste0("ignored_", paste(setdiff(names(chrom_nb), probable_chrom), collapse = "_and_"), "_genes__")
          cat("Ignoring genes in ", setdiff(names(chrom_nb), probable_chrom), " for ", dataset, "\n")
          print(current.genes)
          current.result.both.genes <- subset(current.result.both.genes, chromosome_name %in% probable_chrom)
        } else {
          comment <- ""
        }
        first_exons <- current.result.both.genes[current.result.both.genes$rank %in% c(1), ]
        if (length(unique(first_exons$strand)) == 1) {
          if (unique(first_exons$strand) == "-1") {
            my.col <- "exon_chrom_start"
          } else {
            my.col <- "exon_chrom_end"
          }
          extremities <- unique(first_exons[, c("external_gene_name", my.col)])
          if (nrow(extremities) == 2) {
            results <- c(dataset, genes.to.check, abs(extremities[1, my.col] - extremities[2, my.col]), comment)
          } else {
            cat("Multiple ", my.col, " in ", dataset, " will be averaged.\n")
            print(extremities)
            extremities <- aggregate(first_exons[, my.col], by = list(gene = first_exons$external_gene_name), FUN = mean)
            results <- c(dataset, genes.to.check, abs(extremities$x[1] - extremities$x[2]), paste0(comment, "multiple_first_exon_this_is_average"))
          }
        } else {
          cat("Multiple strand in ", dataset, "\n")
          print(first_exons)
          strand_table <- table(unique(first_exons[, c("external_gene_name", "strand")]))
          if (all(strand_table < 2)) {
            results <- c(dataset, genes.to.check, NA, "different_strand")
          } else {
            print(strand_table)
            results <- c(dataset, genes.to.check, NA, "different_strand_to_do")
          }
        }
        summary_table <- rbind(summary_table, results)
      }
    }
    sink()
    summary_table <- as.data.frame(summary_table)
    colnames(summary_table) <- c("dataset", "genes_considered", "diff_end_both_exons", "comment")
    summary_table <- merge(all_datasets, summary_table)
    write.table(
      summary_table,
      file = "b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13_diff.txt",
      sep = "\t", row.names = F
    )
  } else {
    summary_table <- read.delim("b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13_diff.txt")
  }

  # Only keep species with both distances
  summary_table_plot <- subset(summary_table, !is.na(diff_end_both_exons))
  summary_table_plot <- subset(
    summary_table_plot,
    version %in% summary_table_plot$version[duplicated(summary_table_plot$version)]
  )
  summary_table_plot$species <- sapply(strsplit(summary_table_plot$description, " genes"), head, 1)
  table(table(summary_table_plot$species))
  # Get class from species:
  my.species <- unique(summary_table_plot$species)
  # The list of species has been downloaded from
  # https://www.ensembl.org/info/about/species.html
  species.df <- read.csv("b9b13_evolution/Species.csv")
  colnames(species.df)[1] <- "species"
  # Fix a short cut:
  species.df$species[species.df$species == "P. kingsleyae"] <- "Paramormyrops kingsleyae"

  all(my.species %in% species.df$species)

  species.df <- subset(species.df, species %in% my.species)

  classification1 <- classification(unique(species.df$Taxon.ID))
  classification1bis <- do.call(rbind, classification1)
  colnames(classification1bis) <- paste0("rank_", colnames(classification1bis))
  classification1bis$Taxon.ID <- sapply(strsplit(rownames(classification1bis), "\\."), head, n = 1)
  classification1bis <- merge(classification1bis, species.df)

  all(my.species %in% classification1bis$species)

  all(selected.ranks %in% classification1bis$rank_name)
  classification1bis.selected <- subset(classification1bis, rank_name %in% selected.ranks)
  classification1bis.selected$rank_name <- factor(
    classification1bis.selected$rank_name,
    levels = selected.ranks
  )
  classification1bis.selected <- classification1bis.selected[order(classification1bis.selected$rank_name), ]

  classification1bis.selected <- classification1bis.selected[!duplicated(classification1bis.selected$species), ]

  classification1bis.selected$species[classification1bis.selected$rank_name == "Vertebrata"]

  summary_table_plot <- merge(summary_table_plot, classification1bis.selected[, c("species", "rank_name")], all.x = TRUE)
  table(summary_table_plot$rank_name, exclude = NULL)

  summary_table_plot$diff_end_both_exons <- as.numeric(
    summary_table_plot$diff_end_both_exons
  )

  summary_table_plot$distance <- summary_table_plot$diff_end_both_exons / 1000

  summary_table_plot$pretty.names <- pretty.names[summary_table_plot$rank_name]

  write.table(
    summary_table_plot,
    file = "b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13_diff_plot.txt",
    sep = "\t", row.names = F
  )
} else {
  summary_table_plot <- read.delim("b9b13_evolution/all_Hoxb1_Hoxb9_Hoxb13_diff_plot.txt")
}


summary_table_plot$pretty.names <- factor(
  summary_table_plot$pretty.names,
  levels = unique(pretty.names)
)

# Add a horizontal line corresponding to the deleted allele
deletion.df <- read.delim("annotations/mice_deletion.bed", header = FALSE)
deleted_d <-
  summary_table_plot$distance[summary_table_plot$species == "Mouse" &
    summary_table_plot$genes_considered == "Hoxb9_Hoxb13"] -
  (deletion.df$V3 - deletion.df$V2) / 1000

my.colors = c("Hoxb1_Hoxb9" = "grey", "Hoxb9_Hoxb13" = "black")
my.labels = c("Hoxb1_Hoxb9" = "italic('Hoxb1')*' to '*italic('Hoxb9')", "Hoxb9_Hoxb13" = "italic('Hoxb9')*' to '*italic('Hoxb13')")

summary_table_plot$genes_considered_pretty <- my.labels[summary_table_plot$genes_considered]
names(my.colors) <- my.labels[names(my.colors)]


g <- ggboxplot(
  summary_table_plot,
  x = "pretty.names",
  color = "genes_considered_pretty",
  y = "distance",
  outlier.shape = NA,
  fill = "pretty.names"
) +
  geom_point(aes(col = genes_considered, group = genes_considered), size = 0.5, position = position_jitterdodge()) +
  scale_color_manual(values = my.colors,
    labels = parse_format()) +
  ylab("Distance between genes [kb]") +
  labs(color = "", fill = "") +
  theme(
    legend.position = "right",
    text = element_text(family = fontFamily, size = 6),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm")
  ) +
  scale_y_continuous(expand = c(0, 0.01),  limits = c(0, NA)) +
  geom_hline(yintercept = deleted_d, lty = 2)

set.seed(1)
ggsave("b9b13_evolution/Distance.pdf", g, width = 7, height = 3)

# Which Muridae:
unique(summary_table_plot$species[summary_table_plot$rank_name == "Muridae"])
unique(summary_table_plot$dataset[summary_table_plot$rank_name == "Muridae"])
unique(summary_table_plot$version[summary_table_plot$rank_name == "Muridae"])
# [1] "Algerian mouse" "Mouse"          "Rat"            "Ryukyu mouse"  
# [5] "Shrew mouse"    "Steppe mouse"  
# [1] "mspretus_gene_ensembl"    "mmusculus_gene_ensembl"  
# [3] "rnorvegicus_gene_ensembl" "mcaroli_gene_ensembl"    
# [5] "mpahari_gene_ensembl"     "mspicilegus_gene_ensembl"
# [1] "SPRET_EiJ_v1"    "GRCm39"          "mRatBN7.2"       "CAROLI_EIJ_v1.1"
# [5] "PAHARI_EIJ_v1.1" "MUSP714" 

# Which species into each boxplot
temp.df <- unique(summary_table_plot[, c("species", "pretty.names")])
nrow(temp.df)        
# [1] 157
temp.df.split <- split(temp.df$species, temp.df$pretty.names)
sapply(names(temp.df.split), function(name) {
  cat("In", name, ":\n")
  cat(temp.df.split[[name]], sep = ", ")
  cat("\n")
})
# In Hominidae :
# Bonobo, Chimpanzee, Gorilla, Human, Sumatran orangutan
# In Other Primates :
# Black snub-nosed monkey, Bolivian squirrel monkey, Bushbaby, Coquerel's sifaka, Crab-eating macaque, Drill, Golden snub-nosed monkey, Greater bamboo lemur, Ma's night monkey, Macaque, Mouse Lemur, Olive baboon, Panamanian white-faced capuchin, Pig-tailed macaque, Sooty mangabey, Vervet-AGM, White-tufted-ear marmoset
# In Muridae :
# Algerian mouse, Mouse, Rat, Ryukyu mouse, Shrew mouse, Steppe mouse
# In Other Placentalia :
# American bison, American mink, Arabian camel, Arctic ground squirrel, Beluga whale, Blue whale, Cat, Chacoan peccary, Chinese hamster CHOK1GS, Cow, Degu, Dingo, Dog, Domestic yak, Donkey, Elephant, Eurasian red squirrel, Ferret, Giant panda, Goat, Golden Hamster, Greater horseshoe bat, Guinea Pig, Horse, Hybrid - Bos Indicus, Kangaroo rat, Leopard, Lesser Egyptian jerboa, Lion, Long-tailed chinchilla, Microbat, Naked mole-rat female, Narwhal, Northern American deer mouse, Pig, Pig - Bamei, Pig - Jinhua, Pig - Landrace, Pig - Pietrain, Pig - Tibetan, Pig USMARC, Prairie vole, Rabbit, Red fox, Sheep, Siberian musk deer, Sperm whale, Squirrel, Upper Galilee mountains blind mole rat, Vaquita, Wild yak, Yarkand deer
# In Marsupialia :
# Common wombat, Koala, Opossum, Tasmanian devil
# In Aves :
# Chicken, Duck, Golden eagle, Japanese quail, Kakapo, Pink-footed goose
# In Reptilia :
# Abingdon island giant tortoise, Argentine black and white tegu, Australian saltwater crocodile, Chinese softshell turtle, Eastern brown snake, Goodes thornscrub tortoise, Indian cobra, Mainland tiger snake, Painted turtle, Three-toed box turtle
# In Coelacanth :
# Coelacanth
# In Actinopterygii :
# Asian bonytongue, Atlantic cod, Atlantic herring, Barramundi perch, Bicolor damselfish, Brown trout, Burton's mouthbrooder, Channel bull blenny, Channel catfish, Chinese medaka, Chinook salmon, Climbing perch, Clown anemonefish, Coho salmon, Common carp, Denticle herring, Eastern happy, Electric eel, Fugu, Gilthead seabream, Goldfish, Greater amberjack, Guppy, Huchen, Indian medaka, Japanese medaka HdrR, Javanese ricefish, Large yellow croaker, Lumpfish, Lyretail cichlid, Mangrove rivulus, Mexican tetra, Midas cichlid, Nile tilapia, Northern pike, Orange clownfish, Paramormyrops kingsleyae, Pike-perch, Pinecone soldierfish, Platyfish, Rainbow trout, Red-bellied piranha, Reedfish, Sheepshead minnow, Siamese fighting fish, Spiny chromis, Spotted gar, Tetraodon, Tiger tail seahorse, Turbot, Yellowtail amberjack, Zebra mbuna, Zebrafish, Zig-zag eel
# In Other Vertebrata :
# Elephant shark, Platypus

annot_df <- unique(summary_table_plot[
  ,
  setdiff(colnames(summary_table_plot), c("genes_considered", "diff_end_both_exons", "comment", "distance"))
])
ratio_plot_df <- summary_table_plot %>%
  group_by(species) %>%
  summarize(ratio = distance[genes_considered == "Hoxb9_Hoxb13"] / distance[genes_considered == "Hoxb1_Hoxb9"]) %>%
  as.data.frame()

ratio_plot_df <- merge(ratio_plot_df, annot_df)

g <- ggboxplot(
  ratio_plot_df,
  x = "pretty.names",
  y = "ratio",
  outlier.shape = NA
) +
  # geom_point(aes(col = genes_considered, group = genes_considered, shape = species), size = 2, position = position_jitterdodge(), show.legend = FALSE) +
  geom_point(size = 0.5) +
  scale_shape_manual(values = rep(c(letters, LETTERS), 4)) +
  ylab("Distance between b9 and b13 / Distance between b9 and b1") +
  labs(color = "") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  ) +
  expand_limits(y = 0)

ggsave("b9b13_evolution/Ratio.pdf", g, width = 4, height = 6)
