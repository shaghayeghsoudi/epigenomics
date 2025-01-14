#library(DSS)
library(bsseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2) # For data reshaping
library(pheatmap) # For heatmap visualization
#library(gdata)
#library(kableExtra)
#library(parallel)
#library(foreach)
#library(doParallel)



#make correlation matrix of top 1000 DMLs between R1,R2,...,R10
#Here is the R code to calculate and visualize the correlation matrix of the top 1000 differentially methylated loci (DMLs) between R1, R2, ..., R10 as shown in your example image:
# Code: Correlation Matrix for Top 1000 DMLs

# Simulated data for 1000 DMLs across 10 samples (R1 to R10)
set.seed(123) # For reproducibility
data <- as.data.frame(matrix(runif(1000 * 10, 0, 1), nrow = 1000)) # 1000 rows, 10 columns
colnames(data) <- paste0("R", 1:10) # Sample names: R1, R2, ..., R10

# Compute correlation matrix
cor_matrix <- cor(data, method = "pearson") # Pearson correlation

# Plot heatmap of correlation matrix
pheatmap(cor_matrix, 
         color = colorRampPalette(c("red", "orange", "green"))(50), # Color scale
         display_numbers = TRUE, # Show correlation values
         number_format = "%.2f", # Round to 2 decimal places
         clustering_distance_rows = "euclidean", # Hierarchical clustering
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", # Clustering method
         fontsize_number = 8, # Font size for numbers
         fontsize = 10, # Font size for labels
         main = "DNA Methylation Correlation")





####################################
######## heatmap of top CpGs ######
###################################
### optimized :)

# Define combinations
combinations <- list(
  c("N_T1_T3"),
  c("N_T1_T2"),
  c("T1_T2_T3"),
  c("T1_T2_T3_T4")
)

## cov files
folder_path <- "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/coverage_files_top_CpG_specimens/test2"    #### Path to the folder containing .cov files

dss_path<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_timepoint_effect_TopCpG_regions_per_subject/"
selected_sites_paths <- list(
    paste(dss_path,"Differential_methylation_CpG_timepoint_effect_TopCpG_regions_test3_clean_N_T1_T3.txt"),
    paste(dss_path,"Differential_methylation_CpG_timepoint_effect_TopCpG_regions_test3_clean_N_T1_T2.txt"),
    paste(dss_path,"Differential_methylation_CpG_timepoint_effect_TopCpG_regions_test3_clean_T1_T2_T3.txt"),
    paste(dss_path,"Differential_methylation_CpG_timepoint_effect_TopCpG_regions_test3_clean_N_T1_T2_T3.txt"))
  
### for test    ###
combinations <- list(
  c("N_T1_T2"),
  c("N_T1_T3"))

selected_sites_paths <- list(
    paste(dss_path,"Differential_methylation_CpG_timepoint_effect_TopCpG_regions_test3_clean_N_T1_T2.txt", sep = ""),
    paste(dss_path,"Differential_methylation_CpG_timepoint_effect_TopCpG_regions_test3_clean_N_T1_T3.txt", sep = ""))
#### 

process_combination <- function(patterns, folder_path, selected_sites_path) {
  # List all files in the specified folder
  list_coverage <- list.files(folder_path, full.names = TRUE)
  
  # Read and filter the selected sites based on the given path
  selected_sites <- read.table(file = selected_sites_path, header = TRUE)
  selected_sites <- selected_sites[selected_sites$fdrs <= 0.05, ]
  selected_sites$pos_ID <- paste(selected_sites$chr, selected_sites$pos, sep = "_")
  selected_sites$chr <- paste("chr", selected_sites$chr, sep = "")
  
  # Convert patterns for regex matching
  patty <- gsub("_", "|", patterns)
  
  # Find matching files
  matching_files <- list_coverage[grepl(patty, list_coverage)]
  
  # Read and extract columns from matching files
  attackStats_COV <- lapply(matching_files, function(x) {
    read.table(x, header = FALSE, sep = "\t")[, c(1, 2, 4)]
  })
  
  # Add file paths to the data
  for (i in seq_along(attackStats_COV)) {
    attackStats_COV[[i]] <- cbind(attackStats_COV[[i]], matching_files[i])
  }
  
  # Combine all files into a single data frame
  aa_montecarlo <- do.call("rbind", attackStats_COV)
  names(aa_montecarlo)[4] <- "path"
  
  # Process data
  coverage_good <- aa_montecarlo %>%
    mutate(sample_id = basename(path)) %>%
    mutate(sample_id = sub("(_merged\\.cov|\\.cov)$", "", sample_id)) %>%
    mutate(uniq_ID = paste(V1, V2, sep = "_")) %>%
    filter(uniq_ID %in% selected_sites$pos_ID)
  
  # Extract unique samples and positions
  samples <- unique(coverage_good$sample_id)
  positions <- selected_sites$pos_ID
  
  # Create output matrix
  out_res <- do.call("rbind", lapply(positions, function(pos) {
  # Extract methylation data for the current position
  focal_pos <- t(coverage_good[coverage_good$uniq_ID == pos, "V4"])
  
  # Ensure the number of samples matches the length of focal_pos
  if (length(focal_pos) == length(samples)) {
    rownames(focal_pos) <- pos
    colnames(focal_pos) <- samples
    return(focal_pos)
  } else {
    return(NULL)
  }
}))



# Print the resulting matrix
print(out_res)
  
  # Design data frame for annotation
  design_df <- data.frame("sample" = samples)
  design_df$case <- sapply(design_df$sample, function(x) {
    if (gregexpr("_", x)[[1]][1] == -1) {
      return(NA)
    } else if (length(gregexpr("_", x)[[1]]) == 1) {
      return(sub(".*_", "", x))
    } else {
      return(sub(".*_([^_]+)_.*", "\\1", x))
    }
  })
  
  rownames(design_df) <- design_df$sample
  design_df$case <- sub("^N.*", "N", design_df$case)
  design_df_g <- data.frame("case" = design_df$case)
  rownames(design_df_g) <- rownames(design_df)
  
  # Save the matrix for heatmap input
  write.table(data.frame(out_res), 
              file = paste0("Pheatmap_input_table_", patterns, "_DSS_DMLs.table"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  return(list(matrix = out_res, design = design_df_g))
}


results <- mapply(
  FUN = process_combination,
  patterns = combinations,
  selected_sites_path = selected_sites_paths,
  MoreArgs = list(folder_path = folder_path),
  SIMPLIFY = FALSE
)


##################################
##### simple worked example #####
# Assume `diff_cpgs` is a data frame containing:
# - CpG IDs in rows
# - Beta values for each sample in columns

# Example: Simulated data for demonstration
set.seed(42)
diff_cpgs <- data.frame(
  CpG1 = runif(10, 0, 1),
  CpG2 = runif(10, 0, 1),
  CpG3 = runif(10, 0, 1),
  CpG4 = runif(10, 0, 1),
  CpG5 = runif(10, 0, 1)
)
rownames(diff_cpgs) <- paste0("Sample", 1:10)

> diff_cpgs
              CpG1      CpG2       CpG3        CpG4       CpG5
Sample1  0.9148060 0.4577418 0.90403139 0.737595618 0.37955924
Sample2  0.9370754 0.7191123 0.13871017 0.811055141 0.43577158
Sample3  0.2861395 0.9346722 0.98889173 0.388108283 0.03743103

# Normalize the data (optional, for better contrast)
diff_cpgs_scaled <- scale(diff_cpgs)

# Create the heatmap
pheatmap(
  diff_cpgs_scaled,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,  # Cluster samples
  cluster_cols = TRUE,  # Cluster CpGs
  main = "Heatmap of Top Differentiated CpGs",
  fontsize_row = 10,  # Adjust text size for better readability
  fontsize_col = 10
)


library(ComplexHeatmap)

# Create a heatmap with annotations
Heatmap(
  diff_cpgs_scaled,
  name = "Methylation",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  row_names_gp = gpar(fontsize = 10),  # Adjust font size
  column_names_gp = gpar(fontsize = 10)
)
Customize Annotations
To add annotations for samples or CpGs (e.g., tumor vs. normal), use the annotation_col or annotation_row options in pheatmap or HeatmapAnnotation in ComplexHeatmap.

# Add annotations
annotation_col <- data.frame(
  SampleType = factor(c("Tumor", "Normal", "Tumor", "Normal", "Tumor", "Normal", "Tumor", "Normal", "Tumor", "Normal"))
)
rownames(annotation_col) <- rownames(diff_cpgs)

pheatmap(
  diff_cpgs_scaled,
  annotation_col = annotation_col,
  main = "Heatmap with Annotations"
)



#### very old version
#selected_sites<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_timepoint_effect_TopCpG_regions_per_subject/Differential_methylation_CpG_timepoint_effect_TopCpG_regions_test3_clean_N_T1_T3.txt", header = TRUE)
#selected_sites<-selected_sites[selected_sites$fdrs <= 0.05,]
#selected_sites$pos_ID<-paste(selected_sites$chr , selected_sites$pos, sep = "_")
#selected_sites$chr<-paste("chr",selected_sites$chr,sep = "")


### annotate selectd sites 
## convert DMLs file into Granges
dml_granges <- GRanges(
    seqnames = selected_sites$chr,
    ranges = IRanges(start = selected_sites$pos, end = selected_sites$pos),
    strand = "*",  # Assume no strand information for DMRs
    pos_id<-selected_sites$pos_ID
    #length = dmrs$length,  # Methylation difference
    #nCG = dmrs$nCG ,
    #meanMethy1 = dmrs$meanMethy1,
    #meanMethy2 = dmrs$meanMethy2,
    #diff.Methy = dmrs$diff.Methy,
    #areaStat = dmrs$areaStat

    # p-value
)

#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

### annotate DMRs with annotar ###
annotations_all <- build_annotations(genome = 'hg19', 
                             annotations = c('hg19_basicgenes', # Gene annotations
                                             'hg19_genes_intergenic', # Intergenic regions
                                             'hg19_genes_intronexonboundaries', # Exon-intron boundaries
                                             'hg19_enhancers_fantom', # FANTOM5 Enhancers
                                             'hg19_genes_promoters', # Promoters
                                             'hg19_genes_5UTRs', # 5' UTR regions
                                             'hg19_genes_3UTRs'  # 3' UTR regions
                                             ))  



# Intersect the regions we read in with the annotations
dm_annotated = annotate_regions(
    regions = dml_granges,
    annotations =   annotations_all,
    ignore.strand = TRUE,
    quiet = FALSE)
    # A GRanges object is returned
print(dm_annotated)    
df_dm_annotated = data.frame(dm_annotated)

genesID_focal_annotr <- na.omit(unique(df_dm_annotated$annot.gene_id))    



#############################
### read coverage files #####
list_coverage<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/coverage_files_top_CpG_specimens/test", pattern="*cov", full.names = TRUE)


combinations <- list(
  c("N_T1_T3"),
  c("N_T1_T2"),
  c("T1_T2_T3"),
  c("T1_T2_T3_T4"))


patterns<-unlist(combinations[cc])
  
  patty<-gsub("_","|",patterns)
  matching_files <- list_coverage[grepl(patty, list_coverage)]
  
  attackStats_COV <- lapply(matching_files,function(x) {
     read.table(x, header=FALSE, sep = "\t")[,c(1,2,4)]
     })


for (i in 1:length(attackStats_COV )){
    attackStats_COV[[i]]<-cbind(attackStats_COV[[i]],matching_files[i])
    }
aa_montecarlo<- do.call("rbind", attackStats_COV)      


names(aa_montecarlo)[4]<-"path"
coverage_good<-aa_montecarlo %>% 
    mutate(sample_id=basename(path))%>% 
    mutate(sample_id = sub("(_merged\\.cov|\\.cov)$", "", sample_id)) %>% 
    mutate(uniq_ID=paste(V1,V2, sep = "_"))%>%  
    filter(uniq_ID %in% selected_sites$pos_ID)
    
samples<-unique(coverage_good$sample_id)
positions<-selected_sites$pos_ID


out_res<-NULL
for(jj in 1:length(positions)){

    #for(jj in 1:100){

    focal_pos<-t(coverage_good[coverage_good$uniq_ID == positions[jj],"V4"])   ### V4 is methylation percentage
    if (length(focal_pos) == length(samples)) {

          rownames(focal_pos) <- positions[jj]
          colnames(focal_pos) <- unique(coverage_good$sample_id)

    }   else {

        next
    }
     
    out_res<-rbind(focal_pos,out_res)

}


design_df<-data.frame("sample"=unique(coverage_good$sample_id))
design_df$case <- sapply(design_df$sample, function(x) {
    # Check how many underscores are present
    if (gregexpr("_", x)[[1]][1] == -1) {
    return(NA) # Return NA if no underscores found
     } else if (length(gregexpr("_", x)[[1]]) == 1) {
     return(sub(".*_", "", x))  # Keep characters after the single underscore
     } else {
     return(sub(".*_([^_]+)_.*", "\\1", x))  # Keep characters between two underscores
     }
     })
design_df$case <- sub("^N.*", "N", design_df$case)
 #     sample case
#1 SRC345_T1_D   T1
#2 SRC483_T3_A   T3
#3    SRC486_N    N   

rownames(design_df)<-design_df$sample
#design_df<-dataframe(design_df[,-1])

design_df_g<-data.frame("case"=design_df[,2])
rownames(design_df_g)<-design_df[,1]

ann_colors<-list(
  case = c("N" = "darkgreen",
              "T2" = "blueviolet",
              "T3"="red"))


rownames(out_res)<-NULL   ### format is [1] "matrix" "array" 

### save input for heatmap
write.table(data.frame(out_res), file = "Pheatmap_input_table_three_groups_DSS_DMLs.table", row.names = FALSE, col.names = TRUE, sep = "/t", quote = FALSE)

pdf("Pheatmap_three_groups_DSS_DMLs.pdf", width = 12, height = 14)              

myheat<-pheatmap(out_res,
    main = "DSS three groups comparison",
    cluster_rows = T,
    clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'complete',
         border_color = "black",
         #col = brewer.pal(10, 'RdYlGn'), ## color of heatmap
         cutree_rows = 2, cutree_cols = 4,
         fontsize_col = 7,          # column label font size 
         angle_col = 45,
         annotation_col = design_df_g,
         annotation_colors = ann_colors
         
        
    )
print(myheat)
dev.off()