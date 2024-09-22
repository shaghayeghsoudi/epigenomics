### Annotate differentially methylated regions (DMRs) and perform gene enrichment anakysis

## Step 1: Annotating Differentially Methylated Regions (DMRs)
rm(list = ls())

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
library(annotatr)
library(dbplyr)


dmls<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/DSS-downstreamDifferential_methylation_CpG_T1_T2.txt", header = TRUE)
dmls$chr<-paste("chr",dmls$chr, sep = "")

dmrs<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/DSS-downstream/NT1/DSS-downstreamDifferential_methylation_CpG_N_T1.txt", header = TRUE)


dmrs$chr <-paste("chr",dmrs$chr, sep = "")
dmrs<-dmrs[dmrs$fdr<= 0.05,]

## convert DMR file into Granges
dmr_granges <- GRanges(
    seqnames = dmrs$chr,
    ranges = IRanges(start = dmrs$start, end = dmrs$end),
    strand = "*",  # Assume no strand information for DMRs
    length = dmrs$length,  # Methylation difference
    nCG = dmrs$nCG ,
    meanMethy1 = dmrs$meanMethy1,
    meanMethy2 = dmrs$meanMethy2,
    diff.Methy = dmrs$diff.Methy,
    areaStat = dmrs$areaStat

    # p-value
)


#### DMLs
dml_granges_DMLs <- GRanges(
    seqnames = dmls$chr,
    ranges = IRanges(start = dmls$pos, end = dmls$pos),
    strand = "*",  # Assume no strand information for DMRs
    #length = dmrs$length,  # Methylation difference
    #nCG = dmrs$nCG ,
    meanMethy1 = dmls$mu1,
    meanMethy2 = dmls$mu2,
    diff.Methy = dmls$diff

    # p-value
)

################# plot the results #################
## Create a simple density plot or frequency plot ##
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/simple_histogram_diff_methylation_T1T2_DMRs_DSS.pdf")
volcano_plot_no_pvalue <- ggplot(dmls, aes(x = diff)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7, colour= "black") +  # Histogram of methylation differences
  theme_minimal() +  # Clean theme
  labs(
    title = "Distribution of Methylation Differences (DMLs)",   ####### change these
    x = "Methylation Difference (T1-T2)",
    y = "Frequency"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), col = "red", linetype = "dashed") +
  theme_bw() +
  theme(panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  axis.line = element_line(colour = "black"),
  text=element_text(size=20),
  plot.margin = unit(c(1,1,1,1), "cm")) # Optional threshold for effect size
print(volcano_plot_no_pvalue)
dev.off()



# Create a density plot (alternative to histogram)
pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/densityplot_diff_methylation_T1T2_DMLs_DSS.pdf")
volcano_plot_density <- ggplot(dmls, aes(x = diff)) +
  geom_density(fill = "blue", alpha = 0.5) +  # Density plot of methylation differences
  theme_minimal() +
  labs(
    #title = "Distribution of Methylation Differences (DMLs)",
    x = "Methylation Difference (T1-T2)",
    y = "Density"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), 
  col = "red", 
  linetype = "dashed")+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  axis.line = element_line(colour = "black"),
  text=element_text(size=20),
  plot.margin = unit(c(1,1,1,1), "cm"))

print(volcano_plot_density)
dev.off()



#Annotate the GRanges Object and plotting
# Build annotations for human (e.g., hg19 genome)
## CpG islands only

annotations_genic <- build_annotations(genome = 'hg19', annotations = c(
    'hg19_basicgenes'  # Gene annotations
    #'hg19_cpgs'       # CpG islands, shores, shelves
    #'hg19_genes_intersecting', # Genes overlapping with the DMRs
    #'hg19_ensGene'      # Ensembl genes
))



annotations_cpg <- build_annotations(genome = 'hg19', annotations = c(
    'hg19_cpgs' ))     # CpG islands, shores, shelves



# Annotate the GRanges DMR object
#dmr_annotated<- annotate_regions(regions =dml_granges_DMLs, annotations = annotations_cpg  , ignore.strand = TRUE)
#dmr_annotated_df<-as.data.frame(dmr_annotated)
#dmr_annotated_df$annot.type[dmr_annotated_df$annot.type=="hg19_cpg_inter"]<-"Open Sea"


dmr_annotated<- annotate_regions(regions =dml_granges_DMLs, annotations =annotations_genic   , ignore.strand = TRUE)
dmr_annotated_df<-as.data.frame(dmr_annotated)


write.table((data.frame(dmr_annotated_df)), file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/annotated_genic_DSS_DMLs_T1T2.table",col.names= TRUE, row.names = FALSE, sep = "\t", quote = FALSE) 

#island<-c("hg19_cpg_islands","hg19_cpg_shores","hg19_cpg_shelves","hg19_cpg_inter")
#dmr_annotated_df_islands<-dmr_annotated_df[dmr_annotated_df$annot.type %in%island, ]

dmr_annotated_df$type <-ifelse(dmr_annotated_df$diff.Methy > 0 , "hypo",
                                ifelse(dmr_annotated_df$diff.Methy < 0 , "hyper",
                                  "NA"))

dmr_annotated_df$type_annote<-paste(dmr_annotated_df$annot.type, dmr_annotated_df$type, sep ="_")                                  
#dmr_annotated_df_islands$id<-paste()

counts<-data.frame("combination"=table(dmr_annotated_df$type_annote))
counts$location <- str_extract(counts$combination.Var1, ".*(?=_[^_]*$)")   
counts$location<-gsub("hg19_","", counts$location)
counts$type <- sub(".*_", "", counts$combination.Var1)
colnames(counts)<-c("combination.Var1", "number_of_DMRs" ,   "genomic_region" , "type")

counts$genomic_region<-gsub("genes_","",counts$genomic_region)


pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/barplot_genic_hyper_hypo_annotated_T1T2_DMLs_DSS.pdf")
plotA<-ggplot(counts, aes(fill=type, y=number_of_DMRs, x=genomic_region)) + 
    geom_bar(stat="identity", colour= "black") +
    scale_fill_brewer(palette = "Dark2") +
    labs(y = "Number of DMLs") +
    ggtitle("#of hypo- hyper methylated DMRs T1-T2") +
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


####### boxplot of 

#types<-unique(dmr_annotated_df_islands$annot.type)

#outres<-NULL
#for (ii in 1:length(types)){
#
#  dmr_type<-dmr_annotated_df_islands[dmr_annotated_df_islands$annot.type == types[ii],]
#  group1<-dmr_type[,c("seqnames","start","meanMethy1","annot.symbol", "annot.type" , "type" , "type_annote")]
#  group1$groupID<-rep("group1")
#  colnames(group1)[3]<-"meanMethy"
#
#  group2<-dmr_type[,c("seqnames","start","meanMethy2","annot.symbol"    ,   "annot.type" , "type"  , "type_annote")]
#  group2$groupID<-rep("group2")
#  colnames(group2)[3]<-"meanMethy"
#
#  both<-rbind(group1, group2)
#
##
#outres<-rbind(both,outres)
#}


#ggplot(outres, aes(x = annot.type, y = meanMethy, fill = type)) +  
#geom_boxplot() +
#scale_fill_brewer(palette = "Dark2") 
#labs(y = "Number of DMRs") +
#    ggtitle("#of hypo- hyper methylated DMRs T1 vs T2") +
#    theme_bw() +
#    theme(panel.border = element_blank(), 
#    axis.text.y = element_text(angle = 90),
#    legend.title = element_blank(),
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(), 
#    axis.line = element_line(colour = "black"),
#    text=element_text(size=20),
#    plot.margin = unit(c(1,1,1,1), "cm"))


################
################

### format annotation file
#annotation<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/hg19_annotation.table", header = TRUE, sep = "\t")
#annot_chr1<-annotation[annotation$seqnames=="chr1",]

#annot_chr1_Granges <- GRanges(
#    seqnames = annot_chr1$seqnames,
#    ranges = IRanges(start = annot_chr1$start, end = annot_chr1$end),
#    strand = "*",
#    id = annot_chr1$id ,
#    tx_id= annot_chr1$tx_id,
#    gene_id= annot_chr1$gene_id,
#    symbol= annot_chr1$symbol,
#    type= annot_chr1$type
#
#    )


# format DMR file
#dmrs_chr1<-dmrs[dmrs$chr=="1",]
#dmrs_chr1$chr<-paste("chr",dmrs_chr1$chr,sep = "")


#dmrs_chr1_Granges <- GRanges(
#    seqnames = dmrs_chr1$chr,
#    ranges = IRanges(start = dmrs_chr1$start, end = dmrs_chr1$end),
#    strand = "*",  # Assume no strand information for DMRs
#    length = dmrs_chr1$length,  # Methylation difference
#    nCG = dmrs_chr1$nCG ,
#    meanMethy1 = dmrs_chr1$meanMethy1,
#    meanMethy2 = dmrs_chr1$meanMethy2,
#    diff.Methy = dmrs_chr1$diff.Methy,
#    areaStat = dmrs_chr1$areaStat
#
    # p-value
#)


#overlaps <- findOverlaps(annot_chr1_Granges, dmrs_chr1_Granges,minoverlap = 10)

#overlapping_dmrs_chr1 <- dmrs_chr1_Granges[subjectHits(overlaps)]
#overlapping_annot_chr1 <- annot_chr1_Granges[queryHits(overlaps)]

#mcols(overlapping_dmrs_chr1)
#mcols(overlapping_annot_chr1)



#result_10bp <- data.frame(
#  gr1_seqnames = seqnames(overlapping_dmrs_chr1),
#  gr1_start = start(overlapping_dmrs_chr1),
#  gr1_end = end(overlapping_dmrs_chr1),
#  gr1_length= mcols(overlapping_dmrs_chr1)$length ,
#  gr1_nCG = mcols(overlapping_dmrs_chr1)$nCG,
#  gr1_meanMethy1 = mcols(overlapping_dmrs_chr1)$meanMethy1,
#  gr1_meanMethy2 = mcols(overlapping_dmrs_chr1)$meanMethy2,
#  gr1_diff.Methy = mcols(overlapping_dmrs_chr1)$diff.Methy,
#  gr1_areaStat = mcols(overlapping_dmrs_chr1)$areaStat,
#  gr2_seqnames = seqnames(overlapping_annot_chr1),
#  gr2_start = start(overlapping_annot_chr1),
#  gr2_end = end(overlapping_annot_chr1),
#  gr2_id = mcols(overlapping_annot_chr1)$id,
#  gr2_gene_id = mcols(overlapping_annot_chr1)$gene_id,
#  gr2_symbol = mcols(overlapping_annot_chr1)$symbol,
#  gr2_type = mcols(overlapping_annot_chr1)$type
# 
#
#)


###########################################################
####### annotating DMRs and enrichment analysis ###########
###########################################################

# Assume `dmr_gr` is your GRanges object containing DMRs
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


## load dmr table
dmrs<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/DSS-downstreamDifferential_methylation_Regions_T1_T2.txt", header = TRUE)
dmrs$chr<-paste("chr",dmrs$chr, sep = "")

dmrs_hyper_T1_T1T2 <-dmrs[dmrs$diff.Methy > 0,]
dmrs_hyper_T2_T1T2 <-dmrs[dmrs$diff.Methy < 0,]

list_names<-c("dmrs_hyper_T1_T1T2","dmrs_hyper_T2_T1T2")  ### adjust this
methyl_types<-list(dmrs_hyper_T1_T1T2,dmrs_hyper_T2_T1T2)

for (cc in 1:length(methyl_types)){

    methyl_types_foc<-data.frame(methyl_types[cc])
    names<-list_names[cc]

    dmrs_focal_Granges <- GRanges(     ### convert DMR table into Granges
    seqnames = methyl_types_foc$chr,
    ranges = IRanges(start = methyl_types_foc$start, end = methyl_types_foc$end),
    strand = "*",                      # Assume no strand information for DMRs
    length = methyl_types_foc$length,  
    nCG = methyl_types_foc$nCG ,
    meanMethy1 = methyl_types_foc$meanMethy1,
    meanMethy2 = methyl_types_foc$meanMethy2,
    diff.Methy = methyl_types_foc$diff.Methy,
    areaStat = methyl_types_foc$areaStat
    # p-value
    )

    ### annotate DMRs with annotar ###
    #annotations_all <- build_annotations(genome = 'hg19', 
    #                             annotations = c('hg19_basicgenes', # Gene annotations
    #                                             'hg19_genes_intergenic', # Intergenic regions
    #                                             'hg19_genes_intronexonboundaries', # Exon-intron boundaries
    #                                             'hg19_enhancers_fantom', # FANTOM5 Enhancers
    #                                             'hg19_genes_promoters', # Promoters
    #                                             'hg19_genes_5UTRs', # 5' UTR regions
    #                                             'hg19_genes_3UTRs'  # 3' UTR regions
    #                                             ))  

    dmrs_annotations <- build_annotations(genome = 'hg19', 
                                 annotations = 'hg19_genes_promoters') # Gene annotations                                               


    # Intersect the regions we read in with the annotations
    dm_annotated = annotate_regions(
    regions = dmrs_focal_Granges,
    annotations =   dmrs_annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
    # A GRanges object is returned
    print(dm_annotated)    

    # Coerce to a data.frame
    df_dm_annotated = data.frame(dm_annotated)
    write.table(df_dm_annotated, paste("Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/gene_enrichment/annotated_DMRs_Promoter_annotar_",list_names[cc],".txt", sep =""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

    genesID_focal_annotr <- na.omit(unique(df_dm_annotated$annot.gene_id))                          
    

    # Perform KEGG pathway enrichment analysis
    kegg_results_annoar <- enrichKEGG(gene         = genesID_focal_annotr ,
                           organism     = "hsa",  # hsa is for Homo sapiens
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)


    kegg_results_annoar_df<-as.data.frame(kegg_results_annoar)

    write.table(kegg_results_annoar_df, paste("Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/gene_enrichment/gene_enrichment_KEGG_annotatted_by_annotar_PromoterRegions_",list_names[cc],".txt", sep =""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    ### plot the results 
    pdf(file = paste("Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/gene_enrichment/dot_plot_gene_enrichment_KEGG_annotatted_by_annotar_PromoterRegions_",list_names[cc],".pdf", sep =""))
    dot_kegg_annotar<-dotplot(kegg_results_annoar, showCategory=30,font.size = 8)
    print(dot_kegg_annotar)
    dev.off()
                       
    # Custom ridge plot with ggridges
    #ggplot(kegg_results_annoar_df, aes(x = GeneRatio, y = Description, fill = p.adjust)) + 
    #geom_density_ridges(scale = 3, rel_min_height = 0.01) + 
    #scale_fill_viridis_c() + 
    #labs(title = "Custom Ridge Plot for Gene Enrichment") + 
    #theme_minimal()
    

    ##########################################################################################
    ####  Alternative Option for annotating DMRs: Annotate DMRs to genes with annotatePeak ####
    dmr_annotated_focal_annotpeak <- annotatePeak(dmrs_focal_Granges , TxDb=txdb, annoDb="org.Hs.eg.db", addFlankGeneInfo = TRUE, tssRegion=c(-3000, 3000))
    #head(head(as.data.frame(dmr_annotated_focal@anno)))


    # Extract gene symbols or Entrez IDs
    genesID_focal <- as.data.frame(dmr_annotated_focal_annotpeak)$geneId   # Use geneId if Entrez IDs are preferred
    #genessymbol_focal <- unique(as.data.frame(dmr_annotated)$SYMBOL)  #
    genes_foc <- unique(genesID_focal)  # Remove duplicates


    # Perform GO enrichment analysis
    go_results <- enrichGO(gene         = genes_foc,         # List of gene IDs (Entrez or SYMBOL)
                       OrgDb        = org.Hs.eg.db,  # The organism's annotation database
                       keyType      = "ENTREZID",    # Set to "SYMBOL" if using gene symbols
                       ont          = "ALL",         # Specify GO category: BP (biological process), CC (cellular component), MF (molecular function), or "ALL"
                       pAdjustMethod = "BH",         # Adjust p-values using the Benjamini-Hochberg method
                       pvalueCutoff = 0.05,          # P-value cutoff
                       qvalueCutoff = 0.05)          # Q-value cutoff (adjusted p-value)

    go_results_df<-as.data.frame(go_results)
    write.table(go_results_df, paste("Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/gene_enrichment/gene_enrichment_enrichGO_",list_names[cc],".txt", sep =""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    


    # Perform KEGG pathway enrichment analysis
    kegg_results <- enrichKEGG(gene         = genes_foc,
                           organism     = "hsa",  # hsa is for Homo sapiens
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
    
    kegg_results_df<-data.frame( kegg_results)
    
    
    pdf(file = paste("Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/gene_enrichment/dot_plot_gene_enrichment_KEGG_annotatted_by_annotPeak_",list_names[cc],".pdf", sep =""))
    dot_kegg_annotar<-dotplot(kegg_results, showCategory=30,font.size = 8)
    print(dot_kegg_annotar)
    dev.off()
    write.table(kegg_results_df, paste("Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/gene_enrichment/gene_enrichment_KEGG_annotated_by_annotatePeak",list_names[cc],".txt", sep =""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
                      

     
                      

} ### for loop



dmrs_all_Granges <- GRanges(
    seqnames = dmrs$chr,
    ranges = IRanges(start = dmrs$start, end = dmrs$end),
    strand = "*",  # Assume no strand information for DMRs
    length = dmrs$length,  # Methylation difference
    nCG = dmrs$nCG ,
    meanMethy1 = dmrs$meanMethy1,
    meanMethy2 = dmrs$meanMethy2,
    diff.Methy = dmrs$diff.Methy,
    areaStat = dmrs$areaStat

    # p-value
)

# Annotate DMRs to genes
dmr_annotated <- annotatePeak(dmrs_all_Granges, TxDb=txdb, annoDb="org.Hs.eg.db")

# View the annotated DMRs




# Extract gene symbols or Entrez IDs
genes <- as.data.frame(dmr_annotated)$geneId  # Use geneId if Entrez IDs are preferred
genes <- unique(genes)  # Remove duplicates

#genes <- unique(as.data.frame(dmr_annotated)$SYMBOL)  # Use SYMBOL if you prefer gene symbols


#Perform Gene Enrichment Analysis


# Perform GO enrichment analysis
go_results <- enrichGO(gene         = genes,         # List of gene IDs (Entrez or SYMBOL)
                       OrgDb        = org.Hs.eg.db,  # The organism's annotation database
                       keyType      = "ENTREZID",    # Set to "SYMBOL" if using gene symbols
                       ont          = "ALL",         # Specify GO category: BP (biological process), CC (cellular component), MF (molecular function), or "ALL"
                       pAdjustMethod = "BH",         # Adjust p-values using the Benjamini-Hochberg method
                       pvalueCutoff = 0.05,          # P-value cutoff
                       qvalueCutoff = 0.05)          # Q-value cutoff (adjusted p-value)
                       


gse <- gseGO(geneList=genes, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             pvalueCutoff = 0.05, 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")



# View the results

head(go_results)

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


# Plot the results using a dotplot
dotplot(go_results, showCategory=20, split=".sign") + facet_grid(.~.sign)

### KEGG Pathway Enrichment
# Perform KEGG pathway enrichment analysis
kegg_results <- enrichKEGG(gene         = genes,
                           organism     = "hsa",  # hsa is for Homo sapiens
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)

# View the results
head(kegg_results)

# Plot the results using a barplot
barplot(kegg_results, showCategory=20)

emapplot(kegg_results, showCategory = 10)
#Visualize the Enrichment Results

#Dotplot of GO Results


dotplot(kegg_results, showCategory=20)



emapplot(kegg_results)

ridgeplot(kegg_results) + labs(x = "enrichment distribution")


########################
###### example data

drosphila_example_de

df = read.csv("~/Downloads/drosphila_example_de.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)