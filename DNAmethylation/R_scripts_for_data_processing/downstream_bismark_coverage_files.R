#### process coverage files for DNA methylation percentage
# Load the libraries
library(annotatr)
library(ggplot2)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(viridis)
library(patchwork)
#library(rtracklayer)
#library(ChIPseeker)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # TxDb for the human genome hg19
#library(org.Hs.eg.db)  # Annotation database for human genes
#library(clusterProfiler)
#library(enrichplot)
#library(dbplyr)


#out_fol<- list.files(path= "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/coverage_files_top_CpG_specimens", pattern = "*.cov", full.names = TRUE)    #### Path to the folder containing .cov files
out_fol<- list.files(path= "/home/shsoudi/methylation/full_cohort_analysis/coverage_files_top_CpG_specimens", pattern = "*.cov", full.names = TRUE)    #### Path to the folder containing .cov files


### create output dirtory if it does not exist
output_dir <- "/home/shsoudi/methylation/full_cohort_analysis/out_downsteram_coverage_files_top_CpG_specimens/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


dss<-lapply(out_fol,function(x){

     #read.table(x, header = TRUE) [,c("chr","pos","mu1","mu2" "diff" "fdr")]
      read.table(x, header = FALSE)[,1:6]   
})

for (i in 1:length(dss)){
    dss[[i]]<-cbind(dss[[i]],out_fol[i])
    }

type_data <- do.call("rbind", dss) 

colnames(type_data) <- c("Chromosome", "Position_start", "Position_end","Methylation_Percentage", "Methylated_Coverage", "Unmethylated_Coverage") 
colnames(type_data)[ncol(type_data)]<-"path"
ChrNames <- c(1:22,"X","Y")


type_data_good<-type_data%>% 
    filter(Chromosome%in%ChrNames) %>% 
    filter(Methylated_Coverage + Unmethylated_Coverage >= 5) %>% 
    mutate(Chromosome=paste("chr",Chromosome,sep = "")) %>% 
    mutate(sample_info=basename(path)) %>%
    mutate(sample_info=gsub("(_merged\\.cov|\\.cov)$", "", sample_info)) %>% 
    dplyr::select(-path)
    

# Create a new column in the dataframe with extracted information
type_data_good <- type_data_good %>%
  mutate(Sample_Type = case_when(
    grepl("T1", sample_info) ~ "T1",
    grepl("T2", sample_info) ~ "T2",
    grepl("T3", sample_info) ~ "T3",
    grepl("N", sample_info) ~ "N",
    TRUE ~ NA_character_  # Default for unmatched cases
  ))

samples<-unique(type_data_good$sample_info)

out_res<-NULL
for (ii in 1:length(samples)){

   focal_sample<-type_data_good[type_data_good$sample_info==samples[ii],]
   
   new_data<-data.frame("mean_percentage"= (mean(focal_sample$Methylation_Percentage)))
   new_data$sample_info<-samples[ii]
   new_data$Sample_Type<-unique(focal_sample$Sample_Type)
   out_res<-rbind(new_data,out_res)

}

write.table(out_res, file = paste(output_dir, "mean_methylation_percentage_per_top_CpG_region.txt", sep = ""), col.names = TRUE, row.names = FALSE,sep = "\t", quote = FALSE)

# Plot
box_all_plot<-out_res %>%
  ggplot( aes(x=Sample_Type, y=mean_percentage, fill=Sample_Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +  #### dots
  labs(title = "box plot of mean methylation percentage in top CpG regions",  y = "Mean methylation percentage") +
  theme_minimal()  +
        theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=17),
        axis.title=element_text(size=19))
print(box_all_plot)


### violin plot
violin<-out_res %>%
  ggplot(aes(x = Sample_Type, y = mean_percentage, fill = Sample_Type)) +
  geom_violin(trim = FALSE) +  # Creates a violin plot; set trim=FALSE to show full distribution
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add jittered points for individual data
  labs(title = "Violin Plot of Mean Methylation Percentage by Sample Type", 
       y = "Mean methylation percentage",
       x = "Sample Type") +  # Add labels for X and Y axes
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 19)
  ) 
print(violin)  

pdf(file =  paste(output_dir, "patchwork_Plot_mean_methylation_percentage_top_CpG_region.pdf", sep = ""), height = 7, width = 12)
both<-box_all_plot + violin
print(both)
dev.off()

##############################
### per patient plots ####
##############################


##############################
### per patient plots ####
##############################

### optimized ### => need to be tested
sample_types <- unique(type_data_good$Sample_Type)  ### subset to each sample type (N, T1, T2, T3)

# Loop through each sample type and generate the plots
for (ss in 1:length(sample_types)) {
  
  # Filter data for the current sample type
  sample_data <- type_data_good[type_data_good$Sample_Type== sample_types[ss],]
  pdf(file=paste(output_dir,"methylation_percentage_per_",sample_types[ss], "_sample_top_CpG_regions.pdf", sep = ""),height = 7, width = 15)

  # Create the box plot
  box_plot <- sample_data %>%
    ggplot(aes(x = sample_info, y = Methylation_Percentage, fill = sample_info)) +
    geom_boxplot(outlier.shape = NA) +
    labs(
      title = paste("Box plot of methylation percentage per sample in top CpG region,", sample_types[ss]),
      y = "Methylation Percentage"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 17),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
      axis.title = element_text(size = 19)
    ) +
    scale_fill_viridis_d(option = "turbo")
  
  # Save the plot to a PDF
  print(box_plot)
  dev.off()
}

message("Plots generated and saved to ", output_dir)

# Filter data for each Sample_Type and create a list of plots

#plots <- type_data_test %>%
#  split(.$Sample_Type) %>%  # Split the data by Sample_Type
#  lapply(function(df) {
#    ggplot(df, aes(x = sample_info, y = Methylation_Percentage, fill = sample_info)) +
#      geom_boxplot(outlier.shape = NA) +
#      labs(
#        title = paste("Box Plot for Sample Type:", unique(df$Sample_Type)),
#        y = "Methylation Percentage"
#      ) +
#      theme_minimal() +
#      theme(
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.background = element_blank(),
#        axis.line = element_line(colour = "black"),
#        axis.text = element_text(size = 17),
#        axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
#        axis.title = element_text(size = 19)
#      ) +
#      scale_fill_viridis_d(option = "turbo")
#  })

# Combine the plots into one column
#combined_plot <- wrap_plots(plots, ncol = 1)

# Print the combined plot
#print(combined_plot)


################################################
##### annotate each cpg site for each sample ####


samples<-unique(type_data_good$sample_info)

## annotate to find CpG Islands
  annotations <- build_annotations(genome = "hg19", annotations = c(
     "hg19_cpg_islands",
     "hg19_cpg_shores",   # CpG Shores (up to 2kb from islands)
     "hg19_cpg_shelves"    # CpG Shelves (2–4kb from islands)
     # "hg19_intergenic" ## NOTE: annotatr does not support intergenic
  ))


out_res_annot<-NULL
for (ii in 1:length(samples)){

   focal_sample<-type_data_good[type_data_good$sample_info==samples[ii],c("Chromosome" ,"Position_start", "Position_end" ,"Methylation_Percentage", "sample_info", "Sample_Type")]

    ## convert DML file into Granges
    ## convert *** DMLs *** file into Granges
    granges <- GRanges(
      seqnames = focal_sample$Chromosome,
      ranges = IRanges(start = focal_sample$Position_start,end = focal_sample$Position_end),
      strand = "*",  # Assume no strand information for DMRs
      #length = type_data_good$length,  # Methylation difference
      #nCG = type_data_good$nCG ,
      Methylation_Percentage = focal_sample$Methylation_Percentage,
      info = focal_sample$sample_info,
      Type = focal_sample$Sample_Type
    # p-value
     )

  dm_annotated <- data.frame(annotate_regions(
       regions = granges,
       annotations = annotations,
       ignore.strand = TRUE
    ))

  dm_annotated_df <- dm_annotated[,c("seqnames", "start", "end","Methylation_Percentage", "info" ,"Type","annot.id" ,"annot.type")]
  
  mean_methylation_df <- data.frame(dm_annotated_df %>%
  group_by(annot.type) %>%
  summarise(mean_methylation = mean(Methylation_Percentage, na.rm = TRUE)) %>% 
  mutate(info = samples[ii]) %>% 
  mutate(Type = unique(dm_annotated_df$Type))

  )

  out_res_annot<-rbind(mean_methylation_df,out_res_annot)

}

write.table(out_res_annot, file = paste(output_dir, "mean_methylation_per_CpG_annottaion_type_top_CpG_regions.txt", sep = ""), col.names = TRUE, row.names = FALSE,sep = "\t", quote = FALSE)


#> out_res_annot
#         annot.type mean_methylation        info Type
#1  hg19_cpg_islands         4.385336 SRC491_T3_A   T3
#2  hg19_cpg_shelves        36.438876 SRC491_T3_A   T3
#3   hg19_cpg_shores        11.836438 SRC491_T3_A   T3
#4  hg19_cpg_islands        10.858894 SRC490_T2_D   T2
#5  hg19_cpg_shelves        72.366201 SRC490_T2_D   T2

# Create separate boxplots for each annotation type ###
pdf(file =  paste(output_dir, "BoxPlot_mean_methylation_percentage_per_annotation_type_top_CpG_regions.pdf", sep = ""), height = 8, width = 9)
anot_plot<-ggplot(out_res_annot, aes(x = Type, y = mean_methylation, fill = Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, color = "black") +
  facet_wrap(~annot.type, scales = "free_y") + # Separate boxplots for each annotation type
  theme_minimal() +
  labs(
    title = "Methylation Percentage by Annotation Type",
    x = "Sample Type",
    y = "Mean Methylation Percentage"
  ) +
theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 17),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
      axis.title = element_text(size = 19)
    ) 
print(anot_plot)
dev.off()

### plot with notch (optional)
pdf(file =  paste(output_dir, "BoxPlot_notch_mean_methylation_percentage_per_annotation_type_top_CpG_regions.pdf", sep = ""), height = 8, width = 9)
anot_plot2<-ggplot(out_res_annot, aes(x = Type, y = mean_methylation, fill = Type)) +
  geom_boxplot(notch= TRUE) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, color = "black") +
  facet_wrap(~annot.type, scales = "free_y") + # Separate boxplots for each annotation type
  theme_minimal() +
  labs(
    title = "Methylation Percentage by Annotation Type",
    x = "Sample Type",
    y = "Mean Methylation Percentage"
  ) +
theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 17),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
      axis.title = element_text(size = 19)
    ) 
print(anot_plot2)
dev.off()


### violin plot with oxplot overlaid in it
pdf(file =  paste(output_dir, "Violin_BoxPlot_mean_methylation_percentage_per_annotation_type_top_CpG_regions.pdf", sep = ""), height = 8, width = 9)
violin_box<-ggplot(out_res_annot, aes(x = Type, y = mean_methylation, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.7) + # Create the violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9) + # Overlay the boxplot
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, color = "black") + # Add jittered dots
  facet_wrap(~annot.type, scales = "free_y",drop = FALSE) + # Separate plots for each annotation type
  theme_minimal() +
  labs(
    title = "Methylation Percentage by Annotation Type",
    x = "Sample Type",
    y = "Mean Methylation Percentage"
  ) +
  theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 17),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
      axis.title = element_text(size = 19)
    )  +
  scale_fill_brewer(palette = "Set2") # Optional: Set a color palette

print(violin_box)
dev.off()

#################################

out_res_annot2<-NULL
for (ii in 1:length(samples)){

   focal_sample<-type_data_good[type_data_good$sample_info==samples[ii],c("Chromosome" ,"Position_start", "Position_end" ,"Methylation_Percentage", "sample_info", "Sample_Type")]

    ## convert DML file into Granges
    ## convert *** DMLs *** file into Granges
    granges_s <- GRanges(
      seqnames = focal_sample$Chromosome,
      ranges = IRanges(start = focal_sample$Position_start,end = focal_sample$Position_end),
      strand = "*",  # Assume no strand information for DMRs
      #length = type_data_good$length,  # Methylation difference
      #nCG = type_data_good$nCG ,
      Methylation_Percentage = focal_sample$Methylation_Percentage,
      info = focal_sample$sample_info,
      Type = focal_sample$Sample_Type
    # p-value
     )


     dm_annotated_all <- data.frame(annotate_regions(
       regions = granges_s,
       annotations = annotations,
       ignore.strand = TRUE
    ))
  
   dm_annotated_dfall <- dm_annotated_all[,c("seqnames", "start", "end","Methylation_Percentage", "info" ,"Type","annot.id" ,"annot.type")]
   out_res_annot2<-rbind(dm_annotated_dfall,out_res_annot2)

}   
write.table(out_res_annot2, file = paste(output_dir, "annotated_all_top_CpG_regions.txt", sep = ""), col.names = TRUE, row.names = FALSE,sep = "\t", quote = FALSE)


# Loop through each Type
unique_types <- unique(out_res_annot2$Type)
unique_annots <- unique(out_res_annot2$annot.type)

for (tt in 1:length(unique_types)) { ### loop through each type 

      focal1<-out_res_annot2[out_res_annot2$Type == unique_types[tt] ,]

      for (aa in 1:length(unique_annots)) {
     
       #subset_data <- subset(out_res_annot2, Type == type & annot.type == annot)
       focal2<-focal1[focal1$annot.type == unique_annots[aa],]
       
       # Skip if there are no rows for the combination
       #if (nrow(subset_data) == 0) {
       #next
       #}
    



    # Create the plot
    pdf(file=paste(output_dir,"boxplot_methylation_percentage_per_",unique_types[tt], "_sample_top_CpG_regions_",unique_annots[aa],".pdf", sep = ""),height = 7, width = 17)
    p_annot <- ggplot(focal2, aes(x = info, y = Methylation_Percentage, fill = info)) +
      geom_boxplot(outlier.shape = NA) +
      labs(
        title = paste("Box Plot for Type:", unique_types[tt], "| Annot Type:", unique_annots[aa]),
        x = "Info",
        y = "Methylation Percentage"
      ) +
      theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 15),
      axis.title = element_text(size = 21)
    ) +
    scale_fill_viridis_d(option = "turbo")
    
    print(p_annot)
    dev.off()
    # Save the plot
    #ggsave(
    #  filename = paste0("boxplot_", type, "_", annot, ".png"),
    #  plot = p,
    #  width = 8,
    #  height = 6,
    #  dpi = 300
    #)
  }
}





#### optimized ####

#create_annotation_plots <- function(data, output_dir = NULL) {
#  
#  # Ensure the output directory exists, if provided
#  if (!is.null(output_dir)) {
#    if (!dir.exists(output_dir)) {
#      dir.create(output_dir)
#    }
#  }
#  
#  for (type in unique(data$Type)) {
#
#    type_data <- subset(data, Type == type)
#    
#    for (annot in unique(data$annot.type)) {
#      
#      annot_data <- subset(type_data, annot.type == annot)
#      
#
#      p <- ggplot(annot_data, aes(x = info, y = Methylation_Percentage, fill = info)) +
#        geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplots for each sample
#        geom_jitter(width = 0.2, size = 1.5, alpha = 0.8, color = "black") + # Add jittered points
#        theme_minimal() +
#        labs(
#          title = paste("Methylation Percentage for", type, "-", annot),
#          x = "Sample (info)",
#          y = "Methylation Percentage"
#        ) +
#        theme(
#          axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
#          strip.text = element_text(size = 10, face = "bold")
#        ) +
#        scale_fill_brewer(palette = "Set2") # Optional: Set color palette
#      
#    
#      print(p)
#      
#      if (!is.null(output_dir)) {
#        ggsave(
#          filename = paste0(output_dir, "/", type, "_", gsub("hg19_", "", annot), ".png"),
#          plot = p,
#          width = 10,
#          height = 6,
#          dpi = 300
#        )
#      }
#    }
#  }
#}

# Call the function with your dataset
# Replace `out_res_annot2` with your actual data frame
#create_annotation_plots(out_res_annot2, output_dir = "annotation_plots")

########################
######## optimized #####
##########################
# Function to calculate mean CpG percentage for a file
#calculate_mean_cpg_percentage <- function(file) {
#
#  data <- read.table(file, header = FALSE, sep = "\t", col.names = c(
#    "Chromosome", "Position_start",,"Position_end" "Methylated_Percentage", 
#    "Methylated_Coverage", "Unethylated_Coverage"
#  ))
#  
#  # Ensure the 'Methylation_Percentage' column exists and is numeric
#  if (!"Methylation_Percentage" %in% colnames(data) || !is.numeric(data$Methylation_Percentage)) {
#    stop("The file does not have a 'Methylation_Percentage' column or it is not numeric.")
#  }
#  
#  
#  mean_cpg_percentage <- mean(data$Methylation_Percentage, na.rm = TRUE)
#  
#  
#  return(data.frame(File = basename(file), Mean_CpG_Percentage = mean_cpg_percentage))
#}
#mean_cpg_percentages <- lapply(cov_files, calculate_mean_cpg_percentage)

# Combine results into a data frame
#mean_cpg_results <- do.call(rbind, mean_cpg_percentages)

# Print the results
#print(mean_cpg_results)

######################
######## END #########
######################

