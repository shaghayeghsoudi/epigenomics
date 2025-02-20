
#rm(list = ls())
library(edgeR)
library(dplyr)
library(tidyr)
library(stringr)


#### Rscript to create volcanoplot with edgeR outputs #### 


edge_out<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/edgeR_outputs/pair_merged_regions/GroupN-GroupT1.topTags_5k_merged_regions_lrt.table", header = TRUE)
cosmic_genes<-read.csv(file = "~/Desktop/Census_all.csv")

edge_out_sig<-edge_out[edge_out$FDR <= 0.055,]
edge_cosmic<-edge_out_sig[edge_out_sig$Symbol %in% cosmic_genes$Gene.Symbol ,]
unique_edge_cos<-unique(edge_cosmic$Symbol)
print(unique_edge_cos)


# Calculate -log10 adjusted P-values
edge_out$negLog10FDR <- -log10(edge_out$FDR)

# Define FDR threshold
fdr_threshold <- 0.05
log10_threshold <- -log10(fdr_threshold)  
print(log10_threshold)

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/edgeR_outputs/pair_merged_regions/Volcano_GroupN-GroupT2.topTags_5k_merged_regions.pdf", height = 9, width = 9)
volcano<-EnhancedVolcano(edge_out,
    lab =  ifelse(duplicated(edge_out$Symbol), "", edge_out$Symbol),                              ###edge_out$Symbol, # Label genes
    x = 'logFC', # X-axis: Log2 Fold Change
    y = 'PValue', # Y-axis: Adjusted p-value
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    title = "Volcano Plot of Differentially Expressed Genes",
    pCutoff = 1.3, # Significance threshold Y -axis
    FCcutoff = 2.0, # Fold change threshold
    pointSize = 1.0, # Size of points
    labSize = 5.0, # Label size
    colAlpha = 0.6, # Transparency of points
    col=c('black', 'black', 'black', 'red3'),
    selectLab = unique_edge_cos,
    drawConnectors = TRUE,  # Enable arrows for labels
    widthConnectors = 0.7,  # Adjust arrow thickness
    colConnectors = "black",  # Color of arrows
    arrowheads = TRUE,  # Display arrowheads at the end of connectors
    max.overlaps = 500,  ### this is jut balc and white
    ylim = c(1, 14),
    hline = NULL
)

print(volcano)
dev.off()