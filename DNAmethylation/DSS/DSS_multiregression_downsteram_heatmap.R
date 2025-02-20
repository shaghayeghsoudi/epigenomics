
#############################
############################

#!/usr/bin/env Rscript

## load libraries
library(DSS)
library(bsseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
#library(gdata)
#library(kableExtra)
#library(knitr)
#library(parallel)
#library(foreach)
#library(doParallel)




## samples to exclude ##
root_dir<-("/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/DSS-downstream/three_groups/DSS_input_shared/")
files<-list.files(root_dir)


## load metadata
meta<-read.csv(file = "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/metadata_grant_updated_f.csv", header = TRUE)
meta_filt <- meta %>% 
        #filter(!(Final_ID%in%target )) %>% 
        #mutate(sample_id=gsub("-","_",sample_id)) %>% 
        mutate(sampleid_status=paste(Final_ID,Tissue_Category, sep = "_"))

       
### select meta rows matching with the sample files
meta_matching_rows <- meta_filt %>%
  filter(sapply(Final_ID, function(id) any(str_detect(files, id))))



samples<-meta_matching_rows$Final_ID
dat.list <- vector(mode = "list", length = length(samples))



for (i in 1:length(samples)){
        
        #for (i in 11:12){ for quick testing
        foc<-files[grep(samples[i],list.files(root_dir))]
        dat.list[[i]] <- read.table(paste(root_dir,foc, sep = ""), header=F, col.names=c("chr", "pos", "N", "X"))

        if(nrow(dat.list[[i]]) > 100){
        
        names(dat.list[[i]])<-c("chr","pos","N","X")

       ### filter cites with less than 5 reads
       dat.list[[i]]<-dat.list[[i]][dat.list[[i]]$N >=5, ]

       #fiter out scafolds
       # remove unplaced contigs from reference

       ChrNames <- c(1:19,"X","Y")
       dat.list[[i]]<-dat.list[[i]][dat.list[[i]]$chr %in% ChrNames,]
       names(dat.list)[[i]] <- samples[i]

        } else {

          print(paste("nrow_",samples[i], "_is not greater than threshold", sep = ""))
        }

}

dat.list <- dat.list[!is.na(names(dat.list))] 
BSobj <- makeBSseqData(dat.list, names(dat.list)) 


design_df<-data.frame("design"=names(dat.list))

design_df$case <- sapply(design_df$design, function(x) {
  # Check how many underscores are present
  if (gregexpr("_", x)[[1]][1] == -1) {
    return(NA) # Return NA if no underscores found
  } else if (length(gregexpr("_", x)[[1]]) == 1) {
    return(sub(".*_", "", x))  # Keep characters after the single underscore
  } else {
    return(sub(".*_([^_]+)_.*", "\\1", x))  # Keep characters between two underscores
  }
})

     
# Fit a linear model using "DMLfit.multiFactor"    
DMLfit = DMLfit.multiFactor(BSobj, design=design_df, formula=~case)   ### adjust based on the multiple factors or interctions
#colnames(DMLfit$X)

# Use DMLtest.multiFactor function to test the cell effect
test3 = DMLtest.multiFactor(DMLfit, term="case")
test3_clean<-na.omit(test3)

ChrNames <- c(1:22,"X","Y") 
test3_clean<-test3_clean[test3_clean$chr %in%ChrNames, ]


ix=sort(test3_clean[,"pvals"], index.return=TRUE)$ix
sorted_test<-test3_clean[ix,]
write.table(sorted_test, file = "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/DSS-downstream/three_groups/DSS-outputs/DMLs_three_groups_tissue_type.table", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

DMLs<-callDMR(test3_clean, p.threshold=0.05)
write.table(sorted_test, file = "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/DSS-downstream/three_groups/DSS-outputs/DMRs_three_groups_tissue_type.table", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

     
#######################
####### heatmap #######
library("pheatmap")
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
library(RColorBrewer)

selected_sites<-read.table(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/three_groups/DMLs_three_groups_tissue_type.table", header = TRUE)
#selected_sites<-read.table(file = "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/DSS-downstream/three_groups/DSS-outputs/DMLs_three_groups_tissue_type.table", header = TRUE)
selected_sites<-selected_sites[selected_sites$fdrs <= 0.05,]
selected_sites$pos_ID<-paste(selected_sites$chr , selected_sites$pos, sep = "_")
selected_sites$chr<-paste("chr",selected_sites$chr,sep = "")

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

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene



### annotate DMRs with annotar ###
annotations_all <- build_annotations(genome = 'hg19', 
                             annotations = c('hg19_basicgenes', # Gene annotations
                                             #'hg19_genes_intergenic', # Intergenic regions
                                             #'hg19_genes_intronexonboundaries', # Exon-intron boundaries
                                             #'hg19_enhancers_fantom', # FANTOM5 Enhancers
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
    

# Perform KEGG pathway enrichment analysis
kegg_results_annoar <- enrichKEGG(gene         = genesID_focal_annotr ,
                           organism     = "hsa",  # hsa is for Homo sapiens
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)


kegg_results_annoar_df<-as.data.frame(kegg_results_annoar)


dmr_annotated_focal_annotpeak <- annotatePeak(dml_granges , TxDb=txdb, annoDb="org.Hs.eg.db", addFlankGeneInfo = TRUE, tssRegion=c(-3000, 3000))

# Extract gene symbols or Entrez IDs
genesID_focal <- as.data.frame(dmr_annotated_focal_annotpeak)$geneId   # Use geneId if Entrez IDs are preferred
#genessymbol_focal <- unique(as.data.frame(dmr_annotated)$SYMBOL)  #
genes_foc <- unique(genesID_focal)  # Remove duplicates


#############################
### read coverage files #####
list_coverage<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/analysis_for_grant_Sep24/DSS-downstream/three_groups/bismark_coverage", pattern="*cov", full.names = TRUE)
attackStats_COV <- lapply(list_coverage,function(x) {
     read.table(x, header=FALSE, sep = "\t")[,c(1,2,4)]
     })


for (i in 1:length(attackStats_COV )){
    attackStats_COV[[i]]<-cbind(attackStats_COV[[i]],list_coverage[i])
    }
aa_montecarlo<- do.call("rbind", attackStats_COV)      

#bb <- aa_montecarlo %>%   #### if shared variants fiels are loaded
#     separate(V1, into = c("chr", "pos", "total_reads","methylated"), sep = " ") %>% 
#     mutate(methylation_percentage=(methylated/total_reads)*100 )

names(aa_montecarlo)[4]<-"path"
coverage_good<-aa_montecarlo %>% 
    mutate(sample_id=basename(path))%>% 
    mutate(sample_id=gsub("^(([^_]*_){2}[^_]*)_.*", "\\1", sample_id)) %>% 
    mutate(sample_id=gsub("_R1","",sample_id), uniq_ID = paste(V1, V2, sep = "_")) %>% 
    filter(uniq_ID %in% selected_sites$pos_ID)
    
samples<-unique(coverage_good$sample_id)
positions<-selected_sites$pos_ID

out_res<-NULL
for(jj in 1:length(positions)){

    #for(jj in 1:1500){

    focal_pos<-t(coverage_good[coverage_good$uniq_ID == positions[jj],"V4"])
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

   
rownames(design_df)<-design_df$sample
#design_df<-dataframe(design_df[,-1])

design_df_g<-data.frame("case"=design_df[,2])
rownames(design_df_g)<-design_df[,1]

ann_colors<-list(
  case = c("N" = "darkgreen",
              "T1" = "blueviolet",
              "T2"="red"))
  
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

