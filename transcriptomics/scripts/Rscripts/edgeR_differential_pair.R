### RNAseq analysis with EdgeR for DMLS

rm(list = ls())
# Load required packages
library(edgeR)  # Alternative: DESeq2
library(dplyr)
library(readr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(DOSE)
library(enrichplot)
library(annotatr)
library(ggplot2)
library(GenomicRanges)
library(stringr)
library(RColorBrewer)
library(viridis)
library(patchwork)


### load tximport count files
tximport_dir<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/rnaSeq/full_cohort/salmon/salmon_counts/"
nova54<-read.table(file = paste(tximport_dir,"salmon.counts.NovaSeq54.tsv", sep = ""), header = TRUE, sep = "\t")
nova54<-nova54[,!(colnames(nova54) %in% c("SRC346.T1.H_quant","SRC486.T3.B_quant" ,"SRC486.T3.C_quant" ,"SRC486.T3.D_quant"))]

nova55<-read.table(file = paste(tximport_dir,"salmon.counts.NovaSeq55.tsv", sep = ""), header = TRUE, sep = "\t")
nova55<-nova55[,!(colnames(nova55) %in% c("SRC125.T2.F_quant","SRC346.T3.C_quant"))]

nova57<-read.table(file = paste(tximport_dir,"salmon.counts.NovaSeq57.tsv", sep = ""), header = TRUE, sep = "\t")
#nova57<-nova57[,!(colnames(nova57) %in% c("SRC486.T3.B_quant", "SRC486.T3.C_quant" ,"SRC486.T3.D_quant"))]
colnames(nova57)[colnames(nova57) == "SRC488.T1.D_quant"] <- "SRC346.T1.H_quant"
colnames(nova57)[colnames(nova57) == "SRC488.T3.C_quant"] <- "SRC346.T3.C_quant"



nova61<-read.table(file = paste(tximport_dir,"salmon.counts.NovaSeq61.tsv", sep = ""), header = TRUE, sep = "\t")
nova62<-read.table(file = paste(tximport_dir,"salmon.counts.NovaSeq62.tsv", sep = ""), header = TRUE, sep = "\t")
nova62<-nova62[,!(colnames(nova62) %in% "SRC343.T1.A_quant")]


df_list <- list(nova54,nova55,nova61,nova62)
nova_all <- Reduce(cbind, df_list)

colnames(nova_all)<- gsub("_quant", "", gsub("NE", "", colnames(nova_all)))
sample_names <- colnames(nova_all)



#nova_all$Symbol<-rownames(nova_all)
#rownames(nova_all)<-NULL
#############################
### prepare to run EdgeR ####
#############################
group<- factor(sub(".*\\.(.*?)\\..*", "\\1", colnames(nova_all)))  ### extract timepoint infomration
y <- DGEList(counts = nova_all, group = group)
num_genes_before <- nrow(y$counts)  ## Count genes before filtering
cat("Number of genes before filtering:", num_genes_before, "\n")  56298 ###

## filter by expression 
#keep <- filterByExpr(y)   # Filter out low-expressed genes , 1760
keep <- filterByExpr(y, min.count = 5, min.total.count = 10)  # Less strict

y <- y[keep, , keep.lib.sizes=FALSE]
cat("Genes retained after filtering:", sum(keep)) #3754


y <- normLibSizes(y) # Normalize counts using TMM (Trimmed Mean of M-values)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design) ###perform quasi-likelihood F-tests

# Perform DE analysis for each comparison
DE_results <- list()

comparisons_char<-c("T1-T2","T1-T3","T2-T3")
comparisons <- list(
  "T1-T2" = c("T1", "T2"),
  "T1-T3" = c("T1", "T3"),
  "T2-T3" = c("T2", "T3")
)

out_fol_edge<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/rnaSeq/full_cohort/edgeR/all_regions/"

for (cc in names(comparisons)) {
  contrast <- makeContrasts(contrasts = paste(comparisons[[cc]], collapse = "-"), levels = design)
  test <- glmQLFTest(fit, contrast = contrast)
  DE_results[[cc]] <- topTags(test, n = Inf)$table  # Store results
  write.table(DE_results[[cc]], file= paste(out_fol_edge, "EdgeR_topTags_DEGs_paired_",comparisons_char[cc],".tsv",sep = ""), sep="\t", quote=FALSE, row.names=TRUE, col.names = TRUE)
  
  DE_results[[cc]]<-DE_results[[cc]][DE_results[[cc]]$FDR<0.055,]
  write.table(DE_results[[cc]], file= paste(out_fol_edge, "EdgeR_topTags_DEGs_paired_",comparisons_char[cc],"_FDR_0.05.tsv",sep = ""), sep="\t", quote=FALSE, row.names=TRUE, col.names = TRUE)

}

#################################################
#### some visualization for RNA-seq analysis#####
cancer_gen<-read.csv(file = "~/Desktop/Census_all.csv", header = TRUE)


rns_fol<-("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/rnaSeq/full_cohort/edgeR/all_regions/")
rns_genes<-list(
  "T1-T2" = read.table(paste(rns_fol,"EdgeR_topTags_DEGs_paired_T1-T2.tsv", sep = ""), header = TRUE),
  "T1-T3" = read.table(paste(rns_fol,"EdgeR_topTags_DEGs_paired_T1-T3.tsv", sep = ""), header = TRUE),
  "T2-T3" = read.table(paste(rns_fol,"EdgeR_topTags_DEGs_paired_T2-T3.tsv", sep = ""), header = TRUE)
)

combinations<-c("T1-T2" ,"T1-T3", "T2-T3")


for(cc in 1:length(combinations)){

    res<-data.frame(rns_genes[grep(combinations[cc],names(rns_genes))])
    colnames(res)<-sub(".*?\\..*?\\.(.*)", "\\1", colnames(res))
    res$combination<-names(rns_genes[cc])

    res$Significant <- res$FDR < 0.05 # Adjust threshold if needed
    #res$Highlight <- ifelse(res$Significant == TRUE, "Significant", "Not Significant")
    sig_genes<-rownames(res[res$FDR < 0.05 & rownames(res) %in% cancer_gen$Gene.Symbol,])
    label_data <- res[rownames(res) %in% sig_genes, ] 
    
    pdf(file = paste(out_fol_edge,"volcano_DEG_all_regions",combinations[cc],".pdf", sep = ""), height = 8, width = 8)
    gvolcano<-ggplot(res, aes(x=logFC, y=-log10(FDR), color=Significant)) +
    geom_point(alpha=0.5, size = 2) +  # Plot points
    scale_color_manual(values=c("black", "red")) +  # Gray for non-significant, Red for significant
    theme_minimal() +
    labs(title=paste("Volcano Plot of DEGs_",combinations[cc]), x="Log2 Fold Change", y="-Log10 FDR") +
    theme(legend.position="top") +
    #theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 19),
          axis.title = element_text(size = 21)) +
          # Label only the genes in sig_genes
          geom_text_repel(data=label_data, aes(label=rownames(label_data)), 
          size=4, max.overlaps=50, color="black")  # Labels in black for contrast
          print(gvolcano)
          dev.off()
}


### TO DO: 
#Heatmap (Top DEGs)

#library(pheatmap)
# Select top 50 DEGs based on FDR
#top_genes <- rownames(res[res$FDR < 0.05 & rownames(res) %in% cancer_gen$Gene.Symbol,])
# Extract normalized expression values (e.g., logCPM)
#expr_matrix <- log2(y$counts[top_genes, ] + 1)

#pheatmap(expr_matrix, scale="row", cluster_rows=TRUE, cluster_cols=TRUE,
#         main="Heatmap of Top 50 DEGs")


#gene_of_interest <- "TP53"  # Replace with a gene from your DEGs
#ggplot(res, aes(x=group, y=y$counts[gene_of_interest, ] + 1)) +
#  geom_boxplot() +
#  scale_y_log10() +
#  labs(title=paste("Expression of", gene_of_interest), y="Normalized Counts", x="Group")


