#!/usr/bin/env Rscript

#rm(list = ls())
### phylogenetic analysis with "ape" package with methylation data coverage and DSS diff files
# Load the libraries
#library(annotatr)
library(ggplot2)
#library(GenomicRanges)
library(dplyr)
library(stringr)
#library(RColorBrewer)
#library(viridis)
library(patchwork)
library(ape)
library(ggtree)
#library(phylotools)
#library(rtracklayer)
#library(ChIPseeker)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # TxDb for the human genome hg19
#library(org.Hs.eg.db)  # Annotation database for human genes
#library(clusterProfiler)
#library(enrichplot)
#library(dbplyr)

## load coverage files
cov_fol<- list.files(path= "/home/shsoudi/methylation/full_cohort_analysis/ape_phylogenetics/coverage/coverage_with_PosID", pattern = "*.cov", full.names = TRUE)    #### Path to the folder containing .cov files

## dss diff files pairwise comparisons
diff_fol<-list.files(path = "/home/shsoudi/methylation/full_cohort_analysis/ape_phylogenetics/diff_dss_cpg/dss_outputs_pairwise_top_CpG_regions", pattern = "*.txt", full.names = TRUE)

## dss diff files timepoint effect
#diff_fol<-list.files(path = "/home/shsoudi/methylation/full_cohort_analysis/ape_phylogenetics/diff_dss_cpg/DSS_outputs_timepoint_effect_TopCpG_regions_per_subject", pattern = "*.txt", full.names = TRUE)


#combinations<-c("N|T1","N|T2","N|T3","T1|T2","T1|T3","T2|T3") ## paired_combinations
combinations<-c("N|T3","T1|T2","T1|T3","T2|T3")
#combinations<-c("N|T1|T2|T3","N|T1|T2","N|T1|T3","T1|T2|T3") ## all_combinations

for (cc in 1:length(combinations)) {

    focal_comb<-combinations[cc]
    pattern <- gsub("\\|", "_", focal_comb)
    output_dir <- paste("/home/shsoudi/methylation/full_cohort_analysis/ape_phylogenetics/outdir/pairwise_top_CpG_regions/",pattern,"/",sep = "")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    matching_files <- grep(pattern, diff_fol, value = TRUE)
    diff_file<-read.table(file = matching_files, header = TRUE)[,c("chr","pos","fdr")]  ## fdrs in pairs
    colnames(diff_file)<-c("chr","pos","fdr")


    diff_file_Topsig <- diff_file %>% 
    filter(fdr <= 0.055) %>% 
    arrange(fdr) %>%        # Sort by fdr in ascending order
    #slice_head(n = 5000) %>% 
    #dplyr::select(chr,pos, mu1,mu2,diff,fdr) %>% 
    mutate(pos_id= paste(chr,pos , sep = "_"))

    sig_posID<-diff_file_Topsig$pos_id

   ### process coverage files

    focal_files <- grep(focal_comb, cov_fol, value = TRUE)

    # Function to filter and extract methylation percentage from a coverage file
    coverage_list <- lapply(focal_files, function(file) {
       read.table(file, header = FALSE, stringsAsFactors = FALSE)
    })

    #names(coverage_list) <- basename(cov_fol)
    clean_names <- gsub("(_merged\\.cov|_updated\\.cov)$", "", basename(focal_files))
    names(coverage_list)<-clean_names


   process_coverage_file <- function(coverage_file, sig_posID) {
        coverage_file %>%
        filter(V7 %in% sig_posID) %>%  # Filter rows matching significant CpGs, V7 is posiID: 1_2334, 1_3445
        select(pos_id = V7, rmethylation = V4)  # V4 is methylation percentage
    }

    # Process all files and combine results
    combined_data <- Reduce(
        function(x, y) full_join(x, y, by = "pos_id"),
        lapply(coverage_list, process_coverage_file, sig_posID = sig_posID)
    )


    colnames(combined_data) <- c("pos_id", clean_names)
    top_cpgs_NA10 <- combined_data[rowSums(is.na(combined_data)) <=5 , ]
    write.table(top_cpgs_NA10, file = paste(output_dir,"shared_top_cpgs_with_NA5_",pattern,".table", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    transposed_mat<-t(as.matrix(top_cpgs_NA10[, -1]))
    distance_matrix<-dist(transposed_mat)
    #distance_matrix <- dist(t(as.matrix(top_cpgs[, -1])))  

    hclust_tree <- hclust(distance_matrix, method = "average")  ###check other methods as well

     # Convert hclust object to a phylo object (for compatibility with ape)
    upgma_tree <- as.phylo(hclust_tree)  ### Save this tree
    saveRDS(upgma_tree, file = paste(output_dir,"phylo_tree_shared_top_cpgs_NA_threshold_5_",pattern,".rds", sep = ""))


    #### with NA removed ####
    top_cpgs_0NA <- combined_data[complete.cases(combined_data), ] 
    write.table(top_cpgs_0NA, file = paste(output_dir,"shared_top_cpgs_removed_NA_",pattern,".table", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    transposed_mat_0NA<-t(as.matrix(top_cpgs_0NA[, -1]))
    distance_matrix_0NA<-dist(transposed_mat_0NA)
    #distance_matrix <- dist(t(as.matrix(top_cpgs[, -1])))  

    hclust_tree_0NA <- hclust(distance_matrix_0NA, method = "average")  ###check other methods as well

     # Convert hclust object to a phylo object (for compatibility with ape)
    upgma_tree_0NA <- as.phylo(hclust_tree_0NA)  ### Save this tree
    saveRDS(upgma_tree_0NA, file = paste(output_dir,"phylo_tree_shared_top_cpgs_NA_removed_",pattern,".rds", sep = ""))

}
    
    
#### read the tree file
upgma_tree <- readRDS("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/ape_phylogenetics/top_CpG_regions/pairs/T1_T2/phylo_tree_shared_top_cpgs_NA_removed_T1_T2.rds")


# Plot the UPGMA tree
## types: "phylogram" (the default), "cladogram", "fan", "unrooted", "radial", "tidy",
### font: an integer specifying the type of font for the labels: 1 (plain text), 2 (bold), 3 (italic, the default), or 4 (bold italic).

tip_colors <- ifelse(grepl("T2", upgma_tree$tip.label), "lightslateblue", "indianred3")  # Define tip colors based on their labels


#tip_colors <- ifelse(grepl("T1", upgma_tree$tip.label), "lightslateblue",
#              ifelse(grepl("T2", upgma_tree$tip.label), "indianred3", "seagreen"))
# Define tip colors based on their labels


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/ape_phylogenetics/top_CpG_regions/pairs/T1_T2/Plot_phylo_tree_shared_top_cpgs_NA5_T1-T2.pdf", height = 14, width = 12)
plot_tree<-plot(upgma_tree, tip.col = tip_colors, cex = 1,edge.width = 5,no.margin = FALSE,font = 2)  # Adjust cex for better readability
print(plot_tree)
dev.off()

#plot(upgma_tree, type = "cladogram", main = "UPGMA Tree", edge.width = 2) 
#plot(upgma_tree, type = "radial", main = "UPGMA Tree", edge.width = 2) 
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/ape_phylogenetics/top_CpG_regions/pairs/T1_T2/Unrooted_Plot_phylo_tree_shared_top_cpgs_with_NA5_T1-T2.pdf", height = 22, width = 22)
unroot<-plot(upgma_tree, 
     type = "u",              # Unrooted tree layout
     main = "UPGMA Tree",     # Title
     edge.width = 5,          # Thicker branches
     cex = 1.5,               # Reduce font size of labels
     no.margin = TRUE,
     tip.col = tip_colors)
print(unroot)
dev.off()


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/ape_phylogenetics/top_CpG_regions/pairs/T1_T2/Fan_Plot_phylo_tree_shared_top_cpgs_with_NA5_T1-T2.pdf", height = 22, width = 22)
fan<-plot(upgma_tree, 
     type = "f",              
     #main = "UPGMA Tree",     
     edge.width = 7,
     no.margin = FALSE,
     tip.col = tip_colors,       
     cex = 3,
     font = 2)              
print(fan)
dev.off()


######

#num_clades <- length(upgma_tree$tip.label) # You can adjust this based on clades
#clade_colors <- brewer.pal(min(num_clades, 12), "Set3") # Adjust color scheme as needed

# Identify clades and assign colors
#clade_groups <- cutree(hclust(dist(cophenetic(upgma_tree))), k = 5) # Adjust 'k' as needed
#tip_colors <- clade_colors[clade_groups]

# Plot the fan tree with colors per clade
#fan <- plot(upgma_tree, 
#     type = "fan",              
#     edge.width = 7,
#     no.margin = FALSE,
#     tip.col = tip_colors,       
#     cex = 3,
#     font = 2)  

#print(fan)

# Convert the tree to a ggtree-compatible object
#gg_tree <- ggtree(as.phylo(nj_tree)) +
#  geom_tiplab() +
#  theme_tree2()

################
##### END ######
################
  




#for file in *.cov; do
#  # Exclude files with '_merged.cov' in the name
#  if [[ "$file" != *_merged.cov ]]; then
#    # Add the new column (column 1 and 2 concatenated with '_') to a new file
#    awk '{print $0, $1"_"$2}' "$file" > "coverage_with_PosID/${file%.cov}_updated.cov"
#  fi
#done