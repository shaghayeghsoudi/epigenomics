#### make alluvial plots
# Load libraries
library(ggplot2)
library(ggalluvial)
library(dplyr)
#library(alluvial)
library(GenomicRanges)
library(annotatr)
library(patchwork)

################################################
### load DSS output files for paired combinations

path<- ("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/DSS_outputs_pair_top_CpG_regions/")
files<-list.files(path="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/DSS_outputs_pair_top_CpG_regions", pattern = "^all_methylation_CpG_", full.names = TRUE)

P1<-read.table(file= paste(path,"all_methylation_CpG_T1_T2.txt", sep = ""), header= TRUE)
P2<-read.table(file= paste(path,"all_methylation_CpG_T1_T3.txt", sep = ""), header= TRUE)


ChrNames <- c(1:22,"X","Y")

P1 <- P1 %>% 
    mutate(id= paste(chr,pos,sep = "_")) %>% 
    mutate( methyl_change= ifelse(fdr <=0.059 & diff <0 , "hyper",   ### the second part of the pair is considered
                           ifelse(fdr <=0.059 & diff >0 , "hypo",
                           ifelse(fdr > 0.059 , "no change","NA")))) %>% 
    filter(chr %in% ChrNames)                       


P2 <- P2 %>% 
    mutate(id= paste(chr,pos,sep = "_")) %>% 
    mutate( methyl_change= ifelse(fdr <=0.059 & diff <0 , "hyper",
                           ifelse(fdr <=0.059 & diff >0 , "hypo",
                           ifelse(fdr > 0.059 , "no change","NA")))) %>% 
    filter(chr %in% ChrNames)                       


shared_values <- intersect(P1[,"id"], P2[, "id"])    

P1_shared<- P1 %>% 
    filter(id %in%shared_values) %>% 
    dplyr::select(chr,pos,mu1,mu2,diff,fdr,id, methyl_change)

P2_shared<- P2 %>% 
    filter(id %in%shared_values) %>% 
    dplyr::select(chr,pos,mu1,mu2,diff,fdr,id, methyl_change)

both <- P1_shared %>%      #### both is a shred files across all genomic location (for the function)
    full_join(P2_shared, by = "id") %>% 
    mutate(pair_change = paste(methyl_change.x,methyl_change.y,sep = "_")) %>% 
    mutate( alluvial= ifelse(pair_change ==  "hypo_hypo", "ConsistHypo",
                           ifelse(pair_change == "hypo_no change" , "LoseHypo",
                           ifelse(pair_change == "hypo_hyper" , "SwitchHypoToHyper",
                           ifelse(pair_change == "hyper_hyper" , "ConsistHyper",
                           ifelse(pair_change ==  "hyper_no change", "LoseHyper",
                           ifelse(pair_change == "hyper_hypo" , "SwitchHyperToHypo",
                           ifelse(pair_change ==  "no change_hypo", "GainHypo",
                           ifelse(pair_change == "no change_hyper" , "GainHyper",
                           ifelse(pair_change ==  "no change_no change", "NoChange","NA")
                           )))))))))


freq_alluvial<-data.frame("freq"=table(both$alluvial))
freq_alluvial<-freq_alluvial %>% 
      mutate(status= ifelse(freq.Var1 ==  "ConsistHypo","hypo_hypo",
                           ifelse(freq.Var1 ==  "LoseHypo","hypo_no change",
                           ifelse(freq.Var1 =="SwitchHypoToHyper", "hypo_hyper" , 
                           ifelse(freq.Var1 == "ConsistHyper","hyper_hyper" , 
                           ifelse(freq.Var1 == "LoseHyper", "hyper_no change", 
                           ifelse(freq.Var1 ==  "SwitchHyperToHypo","hyper_hypo" ,
                           ifelse(freq.Var1 ==  "GainHypo","no change_hypo", 
                           ifelse(freq.Var1 == "GainHyper","no change_hyper" , 
                           ifelse(freq.Var1 ==  "NoChange","no change_no change","NA")
                           ))))))))) %>% 
      mutate(primary = gsub("_.*$","",status)) %>% 
      mutate(radiation = gsub(".*_","",status)) %>% 
      mutate(genomic_location ="all" ) 

    

# Define colors
colors <- c(
     "ConsistHyper" = "indianred4",
     "ConsistHypo" = "blue4",
     "GainHyper" = "seagreen4",
     "GainHypo" = "seagreen1",
     "LoseHyper" = "indianred1",
     "LoseHypo" = "blue",
     "NoChange" = "olivedrab",
     "SwitchHyperToHypo" = "orangered",
     "SwitchHypoToHyper" = "cornflowerblue"
)

# **Reorder factor levels**
#data$primary <- factor(data$primary, levels = c("hypo", "hyper", "no change"))
#data$relapse <- factor(data$relapse, levels = c("hypo", "hyper", "no change"))

 # Create Alluvial Plot
 #alluvial(
 #    freq_alluvial[, 4:5],  # Select Primary and Relapse columns
 #    freq = freq_alluvial$freq.Freq,  # Frequency column
 #    col = colors, # Color based on status
 #   border = "black",       # Border color
 #   alpha = 0.6,           # Transparency
 #    cex = 0.7,             # Text size
 #    axis_labels = c("Primary", "Relapse") # Axis labels
 #)


# Create an alluvial plot
colnames(freq_alluvial)[1]<-"Methylation_dynamics"
freq_alluvial$category<-"all"
all_cps<-ggplot(freq_alluvial, aes(axis1 = primary, axis2 = radiation, y = freq.Freq)) +
  geom_alluvium(aes(fill = Methylation_dynamics), width = 4/12) +
  geom_stratum(width = 5/12, fill = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("T1 vs. T2", "T1 vs. T3")) +
  scale_fill_manual(values = colors) +
  labs(title = "Alluvial Plot of Status Changes", y = "Frequency") +
  theme_minimal()  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=17),axis.title=element_text(size=19))


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/alluvial_all_and_CGIs_T1T2_vs_T1T3.pdf", height = 8, width= 8)
print(all_cps)
dev.off()

write.table(freq_alluvial, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/DSS_outputs_pair_top_CpG_regions/out_tables/output_alluvial_table_top_CpG_All_T1T2_T1T3.table", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


##############################################################
#### annotate CpGs and make alluvial plot fot CpG Islands ####
##############################################################

P1_shared<- P1 %>% 
    filter(id %in%shared_values) %>% 
    dplyr::select(chr,pos,mu1,mu2,diff,fdr,id, methyl_change)

shared_values_data<-data.frame("id"=shared_values) 
shared_values_data<-shared_values_data %>% 
    mutate("chr"=gsub("_.*$","",id), "pos" = gsub(".*_","",id)) %>% 
    mutate("chr" = paste("chr", chr,sep = "")) %>% 
    mutate("pos" = as.numeric(pos))


## convert DML file into Granges
granges <- GRanges(
    seqnames = shared_values_data$chr,
    IRanges(start = shared_values_data$pos, end = shared_values_data$pos))

## annotate to find CpG Islands
annotations <- build_annotations(genome = "hg19", annotations = c(
  "hg19_cpg_islands",
  "hg19_cpg_shores",   # CpG Shores (up to 2kb from islands)
  "hg19_cpg_shelves"    # CpG Shelves (2–4kb from islands)
  # "hg19_intergenic" ## NOTE: annotatr does not support intergenic
))

dm_annotated <- data.frame(annotate_regions(
       regions = granges,
       annotations = annotations,
       ignore.strand = TRUE
    ))

dm_annotated_df <- as.data.frame(dm_annotated)[,c("seqnames", "start", "annot.type")]


### optimized plotting #####
# Create a Function for Filtering and Data Preparation
prepare_data <- function(annot_type, dm_annotated_df, both) {   ### annotation types (shore, shelve, island)
  # Filter data for the specific annotation type
    dm_annotated <- dm_annotated_df %>%
    mutate(chr = gsub("chr", "", seqnames), id = paste(chr, start, sep = "_")) %>%
    filter(annot.type == annot_type)
  
    both_filtered <- both[both$id %in% dm_annotated$id, ]
  
   # Create the frequency table and process data
   freq_alluvial <- data.frame("freq" = table(both_filtered$alluvial))
   freq_alluvial <- freq_alluvial %>%
   mutate(status = ifelse(freq.Var1 == "ConsistHypo", "hypo_hypo",
                   ifelse(freq.Var1 == "LoseHypo", "hypo_no change",
                   ifelse(freq.Var1 == "SwitchHypoToHyper", "hypo_hyper",
                   ifelse(freq.Var1 == "ConsistHyper", "hyper_hyper",
                   ifelse(freq.Var1 == "LoseHyper", "hyper_no change",
                   ifelse(freq.Var1 == "SwitchHyperToHypo", "hyper_hypo",
                   ifelse(freq.Var1 == "GainHypo", "no change_hypo",
                   ifelse(freq.Var1 == "GainHyper", "no change_hyper",
                   ifelse(freq.Var1 == "NoChange", "no change_no change", "NA")
                   ))))))))) %>%
    mutate(primary = gsub("_.*$", "", status)) %>%
    mutate(relapse = gsub(".*_", "", status))
  
   colnames(freq_alluvial)[1] <- "Methylation_dynamics"
  
   return(freq_alluvial)
}


# Create a Function for Alluvial Plotting
create_alluvial_plot <- function(freq_alluvial, title, colors) {
  ggplot(freq_alluvial, aes(axis1 = primary, axis2 = relapse, y = freq.Freq)) +
    geom_alluvium(aes(fill = Methylation_dynamics), width = 4/12) +
    geom_stratum(width = 5/12, fill = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("T1 vs. T2", "T1 vs. T3")) +
    scale_fill_manual(values = colors) +
    labs(title = title, y = "Frequency") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 17), axis.title = element_text(size = 19))
}


# Define Groups and Automate Filtering and Plotting
# Define annotation types and their titles
annot_types <- c("hg19_cpg_shelves", "hg19_cpg_shores", "hg19_cpg_islands")
titles <- c("Alluvial Plot of Status Changes (Shelves)",
            "Alluvial Plot of Status Changes (Shores)",
            "Alluvial Plot of Status Changes (Islands)")

# Create a list to store plots
alluvial_plots <- list()

# Loop through annotation types
for (i in seq_along(annot_types)) {
  # Prepare data for the current annotation type
  freq_alluvial <- prepare_data(annot_types[i], dm_annotated_df, both)
  
  # Create the alluvial plot
  alluvial_plots[[i]] <- create_alluvial_plot(freq_alluvial, titles[i], colors)
}


# Combine all plots vertically with a shared legend
combined_plot <- wrap_plots(alluvial_plots) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/alluvial_plots_annotatted_genome_locations_T1T2_vs_T1T3.pdf", combined_plot, width = 15, height = 18)


### create a function for output table

# Create a Function for Filtering and Data Preparation for output table
prepare_data <- function(annot_type, dm_annotated_df, both) {   ### annotation types (shore, shelve, island)
  # Filter data for the specific annotation type
    dm_annotated <- dm_annotated_df %>%
    mutate(chr = gsub("chr", "", seqnames), id = paste(chr, start, sep = "_")) %>%
    filter(annot.type == annot_type)
  
    both_filtered <- both[both$id %in% dm_annotated$id, ]
  
   # Create the frequency table and process data
   freq_alluvial <- data.frame("freq" = table(both_filtered$alluvial))
   freq_alluvial <- freq_alluvial %>%
   mutate(status = ifelse(freq.Var1 == "ConsistHypo", "hypo_hypo",
                   ifelse(freq.Var1 == "LoseHypo", "hypo_no change",
                   ifelse(freq.Var1 == "SwitchHypoToHyper", "hypo_hyper",
                   ifelse(freq.Var1 == "ConsistHyper", "hyper_hyper",
                   ifelse(freq.Var1 == "LoseHyper", "hyper_no change",
                   ifelse(freq.Var1 == "SwitchHyperToHypo", "hyper_hypo",
                   ifelse(freq.Var1 == "GainHypo", "no change_hypo",
                   ifelse(freq.Var1 == "GainHyper", "no change_hyper",
                   ifelse(freq.Var1 == "NoChange", "no change_no change", "NA")
                   ))))))))) %>%
    mutate(primary = gsub("_.*$", "", status)) %>%
    mutate(relapse = gsub(".*_", "", status)) %>% 
    mutate(genomic_location =annot_types[i] ) 

  
   colnames(freq_alluvial)[1] <- "Methylation_dynamics"
  
   return(freq_alluvial)
}


our_res<-NULL
for (i in seq_along(annot_types)) {
  # Prepare data for the current annotation type
  freq_alluvial <- prepare_data(annot_types[i], dm_annotated_df, both)
  our_res<-rbind(freq_alluvial,our_res)

}


write.table(our_res, file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/DSS_outputs_pair_top_CpG_regions/out_tables/output_alluvial_table_top_CpG_annotated_T1T2_T1T3.table", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)




###############################
############ END ##############
###############################
#### Focus on CpG Islands ####
#dm_annotated_island<-dm_annotated_df %>% 
#     mutate(chr = gsub("chr","",seqnames),id = paste (chr,start,sep = "_")) %>% 
#     filter(annot.type =="hg19_cpg_islands")
#
#
#both_island<-both[both$id %in% dm_annotated_island$id,]
#
#
#freq_alluvial_island<-data.frame("freq"=table(both_island$alluvial))
#freq_alluvial_island<-freq_alluvial_island %>%  
#      mutate(status= ifelse(freq.Var1 ==  "ConsistHypo","hypo_hypo",
#                           ifelse(freq.Var1 ==  "LoseHypo","hypo_no change",
#                           ifelse(freq.Var1 =="SwitchHypoToHyper", "hypo_hyper" , 
#                           ifelse(freq.Var1 == "ConsistHyper","hyper_hyper" , 
#                           ifelse(freq.Var1 == "LoseHyper", "hyper_no change", 
#                           ifelse(freq.Var1 ==  "SwitchHyperToHypo","hyper_hypo" ,
#                           ifelse(freq.Var1 ==  "GainHypo","no change_hypo", 
#                           ifelse(freq.Var1 == "GainHyper","no change_hyper" , 
#                           ifelse(freq.Var1 ==  "NoChange","no change_no change","NA")
#                           ))))))))) %>% 
#      mutate(primary = gsub("_.*$","",status)) %>% 
#     mutate(radiation = gsub(".*_","",status))                  
#    
#
#
#colors <- c(
#     "ConsistHyper" = "indianred4",
#     "ConsistHypo" = "blue4",
#     "GainHyper" = "seagreen4",
#     "GainHypo" = "seagreen1",
#     "LoseHyper" = "indianred1",
#     "LoseHypo" = "blue",
#     "NoChange" = "olivedrab",
#     "SwitchHyperToHypo" = "orangered",
#     "SwitchHypoToHyper" = "cornflowerblue"
#)

# Create an alluvial plot
#colnames(freq_alluvial_island)[1]<-"Methylation_dynamics"
#freq_alluvial_island$category<-"Islands"


#Add prefixes to distinguish axes
#freq_alluvial_island$primary <- paste0("P_", freq_alluvial_island$primary)
#freq_alluvial_island$relapse <- paste0("R_", freq_alluvial_island$relapse)

#CGI<-ggplot(freq_alluvial_island, aes(axis1 = primary, axis2 = radiation, y = freq.Freq)) +
#  geom_alluvium(aes(fill = Methylation_dynamics), width = 4/12) +
#  geom_stratum(width = 5/12, fill = "grey") +
#  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#  scale_x_discrete(limits = c("Primary", "Irradiated")) +
#  scale_fill_manual(values = colors) +
#  labs(title = "Alluvial Plot of Status Changes only CGIs", y = "Frequency") +
#  theme_minimal()  +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"),
#        axis.text=element_text(size=17),axis.title=element_text(size=19))


#pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_all_regions_paired_Pvalue1/primary_irradiated/alluvial_all_and_CGIs_primary_irradiated.pdf", height = 8, width= 19)
#combined_plot <- all_cps + CGI + plot_layout(guides = 'collect') & theme(legend.position = "right")
#print(combined_plot)
#dev.off()


######################################################################
######################################################################
# Generate alluvial plots using lapply
alluvial_plots <- lapply(seq_along(annot_types), function(i) {
  freq_alluvial <- prepare_data(annot_types[i], dm_annotated_df, both)
  create_alluvial_plot(freq_alluvial, titles[i], colors)
})

# Combine and display plots
combined_plot <- wrap_plots(alluvial_plots) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)


###############################
###############################

combinations <- list(
  c("N_T1", "N_T2"),
  c("N_T1", "N_T3"),
  c("T1_T2", "T1_T3")
)

for (cc in 1:length(combinations)){  ### loop through each combination

  patterns<-unlist(combinations[cc])
  matching_files <- files[grepl(paste(patterns, collapse = "|"), files)]
  
  processed_files <- lapply(matching_files, function(file) {

    data <- read.table(file, header = TRUE)

    # Apply the mutations
    data <- data %>%
    mutate(id = paste(chr, pos, sep = "_")) %>%
    mutate(methyl_change = ifelse(fdr <= 0.059 & diff < 0, "hyper",
                           ifelse(fdr <= 0.059 & diff > 0, "hypo",
                           ifelse(fdr > 0.059, "no change", "NA")))) %>% 
    filter(chr %in% ChrNames)
  
    # Return the processed data
    return(data)
    }) ### 

   P1 <- processed_files[[1]]  # Processed data for the first file
   P2 <- processed_files[[2]]  # Processed data for the second file
   
   shared_values <- intersect(P1[,"id"], P2[, "id"])   


   P1_shared<- P1 %>% 
    filter(id %in%shared_values) %>% 
    dplyr::select(chr,pos,mu1,mu2,diff,fdr,id, methyl_change)

  P2_shared<- P2 %>% 
    filter(id %in%shared_values) %>% 
    dplyr::select(chr,pos,mu1,mu2,diff,fdr,id, methyl_change)

  both <- P1_shared %>% 
    full_join(P2_shared, by = "id") %>% 
    mutate(pair_change = paste(methyl_change.x,methyl_change.y,sep = "_")) %>% 
    mutate( alluvial= ifelse(pair_change ==  "hypo_hypo", "ConsistHypo",
                           ifelse(pair_change == "hypo_no change" , "LoseHypo",
                           ifelse(pair_change == "hypo_hyper" , "Switch",
                           ifelse(pair_change == "hyper_hyper" , "ConsistHyper",
                           ifelse(pair_change ==  "hyper_no change", "LoseHyper",
                           ifelse(pair_change == "hyper_hypo" , "Switch",
                           ifelse(pair_change ==  "no change_hypo", "GainHypo",
                           ifelse(pair_change == "no change_hyper" , "GainHyper",
                           ifelse(pair_change ==  "no change_no change", "NoChange","NA")
                           )))))))))



freq_alluvial<-data.frame("freq"=table(both$alluvial))
freq_alluvial<-freq_alluvial %>% 
      mutate(status= ifelse(freq.Var1 ==  "ConsistHypo","hypo_hypo",
                           ifelse(freq.Var1 ==  "LoseHypo","hypo_no change",
                           ifelse(freq.Var1 =="Switch", "hypo_hyper" , 
                           ifelse(freq.Var1 == "ConsistHyper","hyper_hyper" , 
                           ifelse(freq.Var1 == "LoseHyper", "hyper_no change", 
                           ifelse(freq.Var1 ==  "Switch","hyper_hypo" ,
                           ifelse(freq.Var1 ==  "GainHypo","no change_hypo", 
                           ifelse(freq.Var1 == "GainHyper","no change_hyper" , 
                           ifelse(freq.Var1 ==  "NoChange","no change_no change","NA")
                           ))))))))) %>% 
      mutate(primary = gsub("_.*$","",status)) %>% 
      mutate(relapse = gsub(".*_","",status))                  
     

  

}  ### cc loop