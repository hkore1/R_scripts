setwd("~/Desktop/treestuff/6_plotting_Phylogenies")

library(ape)
library(ggplot2)
library(dplyr)
library(treeio)
library(ggtree)
library(glue)


###############################################
### First specify some necessary parameters ###
###############################################
# Outgroups and ingroup root (for RT)
outgroups_file <- "IQTREE_maxlikelihood/1_individual_tree_plots/outgroup.txt"
ingroup_root_file <- "IQTREE_maxlikelihood/1_individual_tree_plots/RT_ingroup_root.txt"

# Spreadsheet with sample information
df <- read.csv(file = "_Eriostemon_Group_TC_samples.csv")

# Name for output file
output_name = "IQTREE_maxlikelihood/1_individual_tree_plots/Unedited_Plot_phylogenetics_RT_concordance"

#########################
### Import tree files ###
#########################
# Stats
statfile = read.delim("1_final_input_trees/concordance_trees/phylogenetics_RT_concord.cf.stat", header = T, comment.char="#")

# Tree
tree = treeio::read.iqtree("1_final_input_trees/concordance_trees/phylogenetics_RT_concord.cf.branch")


###################################
### Insert stats into tree file ###
###################################
# Note: code here is a bit weird to accommodate for default mismatch in node IDs...

## Gene Concordance
#   gCF: Gene concordance factor (=gCF_N/gN %)
tree@data$gCF <- c(NA, head(statfile$gCF[match(tree@data$node, statfile$ID)], -1))

#   gDF1: Gene discordance factor for NNI-1 branch (=gDF1_N/gN %)
tree@data$gDF1 <- c(NA, head(statfile$gDF1[match(tree@data$node, statfile$ID)], -1))

#   gDF2: Gene discordance factor for NNI-2 branch (=gDF2_N/gN %)
tree@data$gDF2 <- c(NA, head(statfile$gDF2[match(tree@data$node, statfile$ID)], -1))

#   gDFP: Gene discordance factor due to polyphyly (=gDFP_N/gN %)
tree@data$gDFP <- c(NA, head(statfile$gDFP[match(tree@data$node, statfile$ID)], -1))

## Site Concordance
#   sCF: Site concordance factor averaged over 100 quartets (=sCF_N/sN %)
tree@data$sCF <- c(NA, head(statfile$sCF[match(tree@data$node, statfile$ID)], -1))


#########################
### Rename tip labels ###
#########################
# Read outgroups and rename them to match df, check if tree is RT and if so, root on specified ingroups
if (grepl("RT", output_name) == TRUE){
  outgroups <- scan(ingroup_root_file, character(), quote = "")
  outgroups <- df$raw_reads_filename[match(outgroups, df$hp2out_filename)]
} else {
  outgroups <- scan(outgroups_file, character(), quote = "")
  outgroups <- df$raw_reads_filename[match(outgroups, df$hp2out_filename)]
}

# Rename tree tips
tree@phylo$tip.label <- df$raw_reads_filename[match(tree@phylo$tip.label, df$hp2out_filename)]


##########################
### Reroot on outgroup ###
##########################
tree <- treeio::root(tree, outgroup = outgroups)


#########################################
### Create gene concordance dataframe ###
#########################################
gene_concord <- data.frame(node = tree@data$node,
                        gCF = tree@data$gCF,
                        gDF1 = tree@data$gDF1,
                        gDF2 = tree@data$gDF2,
                        gDFP = tree@data$gDFP)

pies <- nodepie(gene_concord,
                cols=2:5,
                color = c("#000000","#5D5D5D","#B7B7B7","#FFFFFF"),
                outline.color = "black",
                outline.size = 0.3)


##################################
### Create tip label dataframe ###
##################################
# Get herbarium voucher numbers
herb_voucher_number <- df$voucher_herbarium_catalog_number[match(tree@phylo$tip.label, df$raw_reads_filename)]
herb_voucher_number <- gsub("NSW","NSW ", herb_voucher_number)
herb_voucher_number <- gsub("MEL|MEL ","MEL ", herb_voucher_number)
herb_voucher_number <- gsub("MELUD|MEL UD","MELUD ", herb_voucher_number)
herb_voucher_number <- strsplit(herb_voucher_number," ")
herb_voucher_number <- plyr::ldply(herb_voucher_number, rbind)

# Get taxon names
scientific_name <- df$scientific_name[match(tree@phylo$tip.label, df$raw_reads_filename)]

# Bind into dataframe
TaxonDF <- strsplit(scientific_name," ")
TaxonDF <- plyr::ldply(TaxonDF, rbind)
TaxonDF <- cbind(row.names = tree@phylo$tip.label, TaxonDF, herb_voucher_number)
colnames(TaxonDF) <- c("Genus","Species","HerbCode","VoucherNumber")

# New column with plotmath label stylisation
TaxonDF <- dplyr::mutate(TaxonDF, 
                         lab = glue("italic({Genus})~italic({Species})~~({HerbCode}~'{VoucherNumber}')")) 

# Write to tree labels
tree@phylo$tip.label <- TaxonDF$lab


############
### Plot ###
############

## For bottom-up ladderized tree with branch lengths
# p <- ggtree(tree, ladderize = TRUE) +
#   xlim(0, max(tree@phylo$edge.length)*3) +
#   geom_text(aes(label=round(sCF,1), x=branch, fontface=1), nudge_y=0.3, nudge_x = -max(tree@phylo$edge.length)/15, size = 2) +
#   geom_tiplab(parse = T) +
#   geom_treescale(x = 0, y = 0) +
#   coord_cartesian(clip = "off")
# 
# inset(p, pies, width = max(tree@phylo$edge.length), height = max(tree@phylo$edge.length), x = "branch")
# 
# pdf(file = paste0(output_name, ".pdf"), width = 12, height = 16)
# 
# inset(p, pies,
#       width = max(tree@phylo$edge.length)*1.5,
#       height = max(tree@phylo$edge.length)*1.5, x = "branch")
# 
# dev.off()

## For top-down ladderized tree with branch lengths
p <- ggtree(tree, ladderize = TRUE) +
  scale_y_reverse() +
  xlim(0, max(tree@phylo$edge.length)*3) +
  geom_text(aes(label=round(sCF,1), x=branch, fontface=1), nudge_y=0.3, nudge_x = -max(tree@phylo$edge.length)/15, size = 2) +
  geom_text(aes(label=round(gCF,1), x=branch, fontface=2), nudge_y=-0.3, nudge_x = -max(tree@phylo$edge.length)/15, size = 2) +
  geom_tiplab(parse = T) +
  geom_treescale(x = 0, y = 0) +
  coord_cartesian(clip = "off")

pdf(file = paste0(output_name, ".pdf"), width = 12, height = 16)

inset(p, pies,
      width = max(tree@phylo$edge.length)*1.5,
      height = max(tree@phylo$edge.length)*1.5, x = "branch",
      reverse_y = TRUE)

dev.off()





