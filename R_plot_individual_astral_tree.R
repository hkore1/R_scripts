setwd("~/Desktop/treestuff/6_plotting_Phylogenies")

library(ape)
library(ggplot2)
library(treeio)
library(ggtree)
library(ggimage)
library(plyr)
library(glue)

############################################
### Read in tree and set output filename ###
############################################
astral_tree <- read.astral("1_final_input_trees/astral_trees/phylogenetics_RT_collapsed_astral.tre")
output_name = "Unedited_Plot_phylogenetics_RT_collapsed_astral"


#########################
### Rename tip labels ###
#########################

# Read spreadsheet
df <- read.csv(file = "_Eriostemon_Group_TC_samples.csv")

# Read outgroups and rename them to match df, check if tree is RT and if so, root on specified ingroups
if (grepl("RT", output_name) == TRUE){
  outgroups <- scan("RT_ingroup_root.txt", character(), quote = "")
  outgroups <- df$raw_reads_filename[match(outgroups, df$hp2out_filename)]
} else {
  outgroups <- scan("outgroup.txt", character(), quote = "")
  outgroups <- df$raw_reads_filename[match(outgroups, df$hp2out_filename)]
}

# Rename tree tips
astral_tree@phylo$tip.label <- df$raw_reads_filename[match(astral_tree@phylo$tip.label, df$hp2out_filename)]


################################
### Set NA length edges to 0 ###
################################
astral_tree@phylo$edge.length[is.na(astral_tree@phylo$edge.length) == TRUE] <- 0


##########################
### Reroot on outgroup ###
##########################
astral_tree <- treeio::root(astral_tree, outgroup = outgroups)


########################################
### Create quartet support dataframe ###
########################################
q_support <- data.frame(node = astral_tree@data$node,
                        q1 = astral_tree@data$q1,
                        q2 = astral_tree@data$q2,
                        q3 = astral_tree@data$q3)

pies <- nodepie(q_support,
                cols=2:4,
                color = c("black","grey","white"),
                outline.color = "black",
                outline.size = 0.3)


##################################
### Create tip label dataframe ###
##################################
# Get herbarium voucher numbers
herb_voucher_number <- df$voucher_herbarium_catalog_number[match(astral_tree@phylo$tip.label, df$raw_reads_filename)]
herb_voucher_number <- gsub("NSW","NSW ", herb_voucher_number)
herb_voucher_number <- gsub("MEL|MEL ","MEL ", herb_voucher_number)
herb_voucher_number <- gsub("MELUD|MEL UD","MELUD ", herb_voucher_number)
herb_voucher_number <- strsplit(herb_voucher_number," ")
herb_voucher_number <- plyr::ldply(herb_voucher_number, rbind)

# Get taxon names
scientific_name <- df$scientific_name[match(astral_tree@phylo$tip.label, df$raw_reads_filename)]

# Bind into dataframe
TaxonDF <- strsplit(scientific_name," ")
TaxonDF <- plyr::ldply(TaxonDF, rbind)
TaxonDF <- cbind(row.names = astral_tree@phylo$tip.label, TaxonDF, herb_voucher_number)
colnames(TaxonDF) <- c("Genus","Species","HerbCode","VoucherNumber")

# New column with plotmath label stylisation
TaxonDF <- dplyr::mutate(TaxonDF, 
                    lab = glue("italic({Genus})~italic({Species})~~({HerbCode}~'{VoucherNumber}')")) 

# Write to tree labels
astral_tree@phylo$tip.label <- TaxonDF$lab


############
### Plot ###
############
p <- ggtree(astral_tree, ladderize = TRUE, branch.length = "none") + 
  scale_y_reverse() +
  xlim(0, 25) +
  geom_text(aes(label=round(pp1,2), x=branch, fontface=2), nudge_y=1.25) +
  geom_label(aes(label=round(q1,2), x=branch, fontface=1), nudge_y=-1.15, fill='black', colour = "white",  size = 2) +
  geom_tiplab(parse = T) +
  coord_cartesian(clip = "off")

# Save plot
pdf(file = paste0(output_name,".pdf"), width = 12, height = 14)

inset(p, pies, width = 0.05, height = 0.05, reverse_y = TRUE, x = "branch")

dev.off()






