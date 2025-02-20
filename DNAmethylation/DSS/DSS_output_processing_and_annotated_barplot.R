### Process DSS output files. This script process DSS putputs files, annotate them based on their location and make bar plots to comapre hypo-heyper methylation in paiwise comparisons
#rm(list = ls())

# Load the libraries
library(annotatr)
library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # TxDb for the human genome hg19
library(org.Hs.eg.db)  # Annotation database for human genes
library(clusterProfiler)
library(enrichplot)
library(dbplyr)
library(stringr)

### load output files 
out_fol<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/DSS_outputs_pair_comparison_merged_regions",pattern 
= "^Differential_methylation_CpG_", full.names= TRUE) 

dss<-lapply(out_fol,function(x){

     #read.table(x, header = TRUE) [,c("chr","pos","mu1","mu2" "diff" "fdr")]
      read.table(x, header = TRUE)   
})

for (i in 1:length(dss)){
    dss[[i]]<-cbind(dss[[i]],out_fol[i])
    }

type_data <- do.call("rbind", dss) 
colnames(type_data)[ncol(type_data)]<-"path"

ChrNames <- c(1:22,"X","Y")

type_data_good<-type_data%>% 
    filter(chr%in%ChrNames) %>% 
    filter(fdr <= 0.059) %>% 
    mutate(chr=paste("chr",chr,sep = "")) %>% 
    mutate(sample_info=basename(path)) %>%
    mutate(sample_info=gsub("Differential_methylation_CpG_|\\.txt", "", sample_info)) %>%
    #dplyr::select(c(chr,start,end,length,nCG,meanMethy1,meanMethy2,diff.Methy,areaStat,sample_info)) ### this is for DMRs
    dplyr::select(c(chr,pos,mu1,mu2,diff,fdr,sample_info)) ### this is for DMRs


## convert *** DMR *** file into Granges
granges <- GRanges(
    seqnames = type_data_good$chr,
    ranges = IRanges(start = type_data_good$start, end = type_data_good$end),
    strand = "*",  # Assume no strand information for DMRs
    length = type_data_good$length,  # Methylation difference
    nCG = type_data_good$nCG ,
    meanMethy1 = type_data_good$meanMethy1,
    meanMethy2 = type_data_good$meanMethy2,
    diff.Methy = type_data_good$diff.Methy,
    areaStat = type_data_good$areaStat,
    info =type_data_good$sample_info

    # p-value
)


## convert *** DMLs *** file into Granges
granges <- GRanges(
    seqnames = type_data_good$chr,
    ranges = IRanges(start = type_data_good$pos,end = type_data_good$pos),
    strand = "*",  # Assume no strand information for DMRs
    #length = type_data_good$length,  # Methylation difference
    #nCG = type_data_good$nCG ,
    meanMethy1 = type_data_good$mu1,
    meanMethy2 = type_data_good$mu2,
    diff.Methy = type_data_good$diff,
    fdr = type_data_good$fdr,
    info =type_data_good$sample_info
    # p-value
)

combinations<-unique(type_data_good$sample_info)

### plot all combinations in one assemble page

plot_list <- list()

for (i in 1:length(combinations)){

    # Extract focal combination
    focal_combination <- subset(type_data_good, sample_info == combinations[i]) ### in dataframe

    # Create plot
    p <- ggplot(focal_combination, aes(x = diff)) +  #### diff needs to be adjusted for DML or DMR
      geom_histogram(binwidth = 0.03, fill = "blue", alpha = 0.7, colour= "black") +
      theme_minimal() +
      labs(
         x = paste("Methylation Difference_", combinations[i], sep = ""),
         y = "Frequency"
      ) +
      ggtitle(combinations[i]) +
      geom_vline(xintercept = c(-0.2, 0.2), col = "red", linetype = "dashed") +
      theme_bw() +
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            text=element_text(size=20),
            plot.margin = unit(c(1,1,1,1), "cm")) 

    # Append the plot to the list using numeric index
    plot_list[[i]] <- p
}

# Combine all plots using patchwork
library(patchwork)
combined_plot <- wrap_plots(plot_list) + plot_layout(ncol = 2)

# Display the combined plot
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/DSS_outputs_pair_comparison_merged_regions/histogram_DMLs_pairwise_merged_regions_CpG_DSS.pdf", height = 18, width=12)
print(combined_plot)
dev.off()

#########################################################
############### annotation with annotatr #################

#view supported annotations by running:
#supported_annotations()
#[1] "hg19_cpg_islands"     
#[2] "hg19_cpg_shores"      
#[3] "hg19_cpg_shelves"     
#[4] "hg19_genes_promoters" 
#[5] "hg19_genes_5UTRs"     
#[6] "hg19_genes_3UTRs"     
#[7] "hg19_genes_introns"   
#[8] "hg19_genes_exons"  

#Load Annotations for CGIs, Shores, Shelves, and open seas
# Build annotations for hg19 genome
annotations <- build_annotations(genome = "hg19", annotations = c(
  "hg19_cpg_islands",
  "hg19_cpg_shores",   # CpG Shores (up to 2kb from islands)
  "hg19_cpg_shelves"    # CpG Shelves (2–4kb from islands)
  # "hg19_intergenic" ## NOTE: annotatr does not support intergenic
))


#annotations_genic <- build_annotations(genome = 'hg19', annotations = c(
    #'hg19_basicgenes'  # Gene annotations
    #'hg19_cpgs'       # CpG islands, shores, shelves
    #'hg19_genes_intersecting', # Genes overlapping with the DMRs
    #'hg19_ensGene'      # Ensembl genes
#))


for (i in 1:length(combinations)){

    gr_focal <- granges[grep(combinations[i], mcols(granges)$info)]
    focal_data<-type_data_good[type_data_good$sample_info ==combinations[i],]

    dm_annotated <- data.frame(annotate_regions(
       regions = gr_focal,
       annotations = annotations,
       ignore.strand = TRUE
    ))

    dm_annotated_df <- as.data.frame(dm_annotated)

    # Merge annotations with original DMRs based on start position
    annotated_dm_full <- merge(focal_data, dm_annotated_df[, c("seqnames", "start", "annot.type")],
       #by.x = c("chr", "start"),  #### if it is DMR
       by.x = c("chr", "pos"),     ### for DML
       by.y = c("seqnames", "start"),
       all.x = TRUE  # Keeps all original rows, adds NA where no overlap
    )

    # Replace NA annotations with 'Open Seas'
   annotated_dm_full <- annotated_dm_full %>%
       mutate(annot.type = ifelse(is.na(annot.type), "Open Seas", annot.type)) %>% 
       mutate(annot.type = gsub("hg19_","",annot.type))


    #write.table((data.frame(dmr_annotated_df)), file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/annotated_CpG_DSS_DMLs_N_T2.table",col.names= TRUE, row.names = FALSE, sep = "\t", quote = FALSE) 

    combination_p1<-gsub("_.*$","",combinations[i])  ### extrcat the first part of the pair
    combination_p2<-sub(".*_", "", combinations[i])  ### extrcat the first part of the pair
    

    annotated_dm_full<- annotated_dm_full %>% 
        mutate(type_p1=ifelse(annotated_dm_full$diff > 0 , "hyper",
                                ifelse(annotated_dm_full$diff < 0 , "hypo",
                                  "NA"))) %>% 
         mutate(type_p1= paste(type_p1,annot.type,sep = "_")) %>% 
         mutate(type_p2=ifelse(annotated_dm_full$diff < 0 , "hyper",
                                ifelse(annotated_dm_full$diff > 0 , "hypo",
                                  "NA"))) %>% 
         mutate(type_p2= paste(type_p2,annot.type,sep = "_"))                         


    colnames(annotated_dm_full)[colnames(annotated_dm_full) == "type_p1"] <- combination_p1[1]
    colnames(annotated_dm_full)[colnames(annotated_dm_full) == "type_p2"] <- combination_p2[1]




key_col<-colnames(annotated_dm_full)[ncol(annotated_dm_full)]  ### find the key column (second part of the combination)
counts<-data.frame("combination"=table(annotated_dm_full[,key_col]))  #### select which part of the pair combination should be used


counts<- counts %>%
mutate(methylation =str_extract(combination.Var1, ".*(?=_[^_]*$)")) %>% 
mutate(combination.Var1 =gsub("_cpg","",combination.Var1)) %>% 
mutate(location=gsub(".*_","",combination.Var1)) %>% 
mutate(methylation =gsub("_cpg","", methylation)) %>% 
mutate(time_point = rep(key_col)) 
colnames(counts)<-c("combination.Var1", "number_of_DMLs" ,   "methylation" , "location", "time_point")



pdf(file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/DSS_outputs_pairwise_comparisons/DSS_outputs_pairwise_comparisons_Pvalue1/DSS_outputs_pair_comparison_merged_regions/barplot_hyper_hypo_",combinations[i],"_DMLs_merged_regions_DSS.pdf",sep = ""))
plotA<-ggplot(counts, aes(fill=methylation, y=number_of_DMLs, x=location)) + 
    geom_bar(stat="identity", colour= "black") +
    scale_fill_brewer(palette = "Dark2") +
    labs(y = "Number of DMLs") +
    ggtitle(paste("#of hypo-hyper DMLs_",key_col,"_in",combinations[i], sep = "")) +
    theme_bw() +
    theme(panel.border = element_blank(), 
    axis.text.y = element_text(angle = 90),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    text=element_text(size=16),
    plot.margin = unit(c(1,1,1,1), "cm"))
print(plotA)
dev.off()    


} ## end of combination loop



### fix integenic issue 
annotations <- build_annotations(genome = "hg19", annotations = c(
  "hg19_cpg_islands",
  "hg19_cpg_shores",
  "hg19_cpg_shelves",
  "hg19_genes_promoters",
  "hg19_genes_introns",
  "hg19_genes_exons"
))


##################
####### END ######
