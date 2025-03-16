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


##############################################
##### Identify DEGs Associated with DMLs #####
##############################################
dml_fol<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/edgeR_outputs/out_edgeR_merged_regions/topTags_5k_merged_regions/"
DML_genes <- list(
  "T1-T2" = read.table(paste(dml_fol,"GroupT1-GroupT2.topTags_5k_merged_regions_lrt.table", sep = ""), header = TRUE)[c("Chr", "Locus", "EntrezID", "Symbol",  "logFC", "FDR")],
  "T1-T3" = read.table(paste(dml_fol,"GroupT1-GroupT3.topTags_5k_merged_regions_lrt.table", sep = ""), header = TRUE)[c("Chr", "Locus", "EntrezID", "Symbol",  "logFC", "FDR")],
  "T2-T3" = read.table(paste(dml_fol,"GroupT2-GroupT3.topTags_5k_merged_regions_lrt.table", sep = ""), header = TRUE)[c("Chr", "Locus", "EntrezID", "Symbol",  "logFC", "FDR")]
)

rns_fol<-("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/rnaSeq/full_cohort/edgeR/all_regions/")
rns_genes<-list(
  "T1-T2" = read.table(paste(rns_fol,"EdgeR_topTags_DEGs_paired_T1-T2_FDR_0.05.tsv", sep = ""), header = TRUE)[c("logFC","FDR")],
  "T1-T3" = read.table(paste(rns_fol,"EdgeR_topTags_DEGs_paired_T1-T3_FDR_0.05.tsv", sep = ""), header = TRUE)[c("logFC","FDR")],
  "T2-T3" = read.table(paste(rns_fol,"EdgeR_topTags_DEGs_paired_T2-T3_FDR_0.05.tsv", sep = ""), header = TRUE)[c("logFC","FDR")]
)


combinations<-c("T1-T2" ,"T1-T3", "T2-T3")
out_enrichment_dir<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/rnaSeq/full_cohort/edgeR/all_regions/enrichment_analysis/DMLs_EdgeR_pairwise/"

for(cc in 1:length(combinations)){

    focal_dml<-data.frame(DML_genes[grep(combinations[cc],names(DML_genes))])
    colnames(focal_dml)<-sub(".*?\\..*?\\.(.*)", "\\1", colnames(focal_dml))
    focal_dml$combination<-names(DML_genes[cc])
    focal_dml<-focal_dml[focal_dml$FDR <= 0.056,]

    focal_rns<-data.frame(rns_genes[grep(combinations[cc],names(rns_genes))])
    colnames(focal_rns)<-sub(".*?\\..*?\\.(.*)", "\\1", colnames(focal_rns))
    focal_rns$combination<-names(rns_genes[cc])

    #focal_rns$Symbol<-rownames(focal_rns)  #### TO Complete
    #DML_DEG_sigtab<-merge(focal_dml,focal_rns, by.x = "Symbol", by.y = "Symbol")

    if(unique(focal_dml$combination) == unique(focal_rns$combination)){

        DML_DEG_overlap<-focal_dml[focal_dml$Symbol %in% rownames(focal_rns),]
        dml_seconpair_hyper<-DML_DEG_overlap[DML_DEG_overlap$logFC <0,]
        dml_seconpair_hyper$status<-"dml_seconpair_hyper"
        
        
        dml_seconpair_hypo<-DML_DEG_overlap[DML_DEG_overlap$logFC >0,]
        dml_seconpair_hypo$status<-"dml_seconpair_hypo"


        both<-rbind(dml_seconpair_hyper,dml_seconpair_hypo)
        status<-unique(both$status)

        for(jj in 1:length(status)){

            both_focal<-both[both$status==status[jj],]
            gene_list <- unique(both_focal$Symbol)

            ##  Enrichment analysis
            gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


            ### Perform GO Enrichment with enrichGO ###
            go_res <- enrichGO(
                gene = gene_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05
                )

            saveRDS(go_res, file = paste(out_enrichment_dir,"enrichGO_",combinations[cc],"_",status[jj],".rds", sep = ""))
            
            pdf(file = paste(out_enrichment_dir,"barplot_count_enrichGO_",combinations[cc],"_",status[jj],".pdf", sep = ""), width = 9, height = 9)
            plota<-barplot(go_res) 
            print(plota)
            dev.off()

            pdf(file = paste(out_enrichment_dir,"barplot_qscore_enrichGO_",combinations[cc],"_",status[jj],".pdf", sep = ""), width = 9, height = 9)
            plotb<-mutate(go_res, qscore = -log(p.adjust, base=10)) %>% 
            barplot(x="qscore")
            print(plotb)
            dev.off()

            
            ## extrcat table 
            go_result<-go_res@result

            significant_go <- go_result %>%
            filter(p.adjust < 0.05) %>%
            arrange(GeneRatio)

            write.table(significant_go, file = paste(out_enrichment_dir,"enrichGO_",combinations[cc],"_",status[jj],"p.adjust0.05.table", sep = ""))

            # Sort data by Fold Enrichment for better visualization
            #go_result_sig<-go_result[go_result$qvalue<0.056,]
            #significant_go  <- go_result_sig %>% arrange(FoldEnrichment)


            ###################
            ### KEGG ##########
            kegg_res <- enrichKEGG(
                gene = gene_entrez$ENTREZID,
                organism = "hsa",  # Human
                pvalueCutoff = 0.05
            )
            saveRDS(kegg_res, file = paste(out_enrichment_dir,"KEGG_",combinations[cc],"_",status[jj],".rds", sep = ""))


            pdf(file = paste(out_enrichment_dir,"barplot_count_KEGG_",combinations[cc],"_",status[jj],".pdf", sep = ""), width = 9, height = 9)
            plotc<-barplot(kegg_res) 
            print(plotc)
            dev.off()

            # View results
            kegg<-kegg_res@result

            significant_kegg <- kegg %>%
            filter(p.adjust < 0.05) %>%
            arrange(GeneRatio)

            write.table(significant_kegg, file = paste(out_enrichment_dir,"KEGG_",combinations[cc],"_",status[jj],"p.adjust0.05.table", sep = ""))

        }
        
 

    }  ## if statement




} ## cc loop



DML_DEG_overlap <- list()

for (comp in names(DML_genes)) {
  DEG_genes <- rownames(DE_results[[comp]])  # Extract DEG gene names
  overlap_genes <- intersect(DML_genes[[comp]], DEG_genes)  # Find overlapping genes
  DML_DEG_overlap[[comp]] <- overlap_genes
  
  # Save overlap results
  write.table(overlap_genes, paste0(comp, "_DML_DEG_overlap.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}


#############
#### END ####
methy<-read.table(file="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/full_cohort/edgeR_outputs/pair_merged_regions/GroupT1-GroupT2.topTags_5k_merged_regions_lrt.table", header = TRUE)
methy_sig<-methy[methy$FDR <= 0.055,]
methy_sig_hyper<-methy_sig[methy_sig$logFC <0,]
methy_sig_hypo<-methy_sig[methy_sig$logFC >0,]

rna<-read.table(file="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/rnaSeq/full_cohort/edgeR/merged_regions/T1-T2_DEGs.tsv", header = TRUE)
sig_DEGs <- rna[rna$FDR <= 0.055 , ]
