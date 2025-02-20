#!/usr/bin/env Rscript

### edgeR updated (no need for target file)
### load required libraries
#rm(list = ls())
library(edgeR)
library(dplyr)
library(tidyr)
library(stringr)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(limma)
#library(parallel)
#library(foreach)

cov_fol<-list.files("/home/shsoudi/methylation/full_cohort_analysis/coverage_files_merged_regions/zipped",pattern = "*cov.gz", full.names= TRUE) 
cov_fol_df <- data.frame(file_path = cov_fol)


output_dir <- "/home/shsoudi/methylation/full_cohort_analysis/out_edgeR_merged_regions/"

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Directory created:", output_dir, "\n")
} else {
  cat("Directory already exists:", output_dir, "\n")
}

## create target file from coverage file
targets<-cov_fol_df %>% 
     mutate(File = basename(file_path)) %>% 
     mutate(Sample =sub("(_merged\\.cov.gz|\\.cov.gz|\\_merged_DSS_input.txt_bismark.cov.gz)$", "", File)) %>% 
     mutate(Group = str_extract(Sample, "T[1-3]|N")) %>% 
     mutate(GEO = seq(1:length((cov_fol)))) %>% 
     mutate(GEO =paste("ST",GEO, sep = "")) %>% 
     dplyr::select(GEO,Sample, Group , File)
     

yall <- readBismark2DGE(cov_fol, sample.names=targets$Sample)
table(yall$genes$Chr)
#yall$genes$Chr <- factor(yall$genes$Chr, levels=ChrNames)
ChrNames <- as.character(c(1:22, "X", "Y"))  # Ensure it's a character vector
yall <- yall[yall$genes$Chr %in% ChrNames, ]
o <- order(yall$genes$Chr, yall$genes$Locus)
yall <- yall[o,]


### annotating the CpG with the identity of the nearest gene
TSS <- nearestTSS(yall$genes$Chr, yall$genes$Locus, species="Hs")
yall$genes$EntrezID <- TSS$gene_id
yall$genes$Symbol <- TSS$symbol
yall$genes$Strand <- TSS$strand
yall$genes$Distance <- TSS$distance
yall$genes$Width <- TSS$width


### Filtering and normalization 
### suming up the counts of methylated and unmethylated reads to get the total read coverage at each CpG site for each sample
Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))
Me <- yall$counts[, Methylation=="Me"]
Un <- yall$counts[, Methylation=="Un"]
Coverage <- Me + Un
head(Coverage)


#HasCoverage <- rowSums(Coverage >= 5 ) == length(targets$Sample)  ### 20 total number of samples, seems very conservative
HasCoverage <- rowSums(Coverage >= 5 ) == round(length(targets$Sample)/2)  ## at least half of the samples
table(HasCoverage)
true_count <- sum(HasCoverage == TRUE)


# filter out CpGs that are never methylated or always methylated, they provide no information about differential methylation:
HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0
table(HasCoverage, HasBoth)

#The DGEList object is subsetted to retain only the non-filtered loci:
#y <- yall[HasCoverage,, keep.lib.sizes=FALSE]
y <- yall[HasCoverage & HasBoth,, keep.lib.sizes=FALSE]


TotalLibSize <- y$samples$lib.size[Methylation=="Me"] +
            + y$samples$lib.size[Methylation=="Un"]

y$samples$lib.size <- rep(TotalLibSize, each=2)
y$samples


## Data exploration
Me <- y$counts[, Methylation=="Me"]
Un <- y$counts[, Methylation=="Un"]
M <- log2(Me + 2) - log2(Un + 2)
colnames(M) <- targets$Sample

#group <- targets$group
#y <- DGEList(counts = M, group = group)
#y$samples

## sample and RT based
#pdf("M-plot-sample-RTbased-colored.pdf", height = 10, width = 10)
#plotMDS(M, col=c(rep("black",1), rep("firebrick",2),rep("blue",1),rep("green",2),rep("orange",2),rep("chartreuse4",1),rep("chocolate4",2),rep("chartreuse4",1),rep("royalblue4",3),rep("chartreuse4",1),rep("chartreuse4",1),rep("chartreuse4",1),rep("chartreuse4",1),rep("steelblue4",1)), main="M-values",cex=0.8)
#dev.off()

## RT based ##
#pdf("M-plot-RTbased-colored_adjusted.pdf", height = 10, width = 10)
#par(mar = c(7, 7, 7, 7)) # Set the margin on all sides to 6
#plotMDS(M, pch =19,
#     col=c(rep("forestgreen",3), 
#     rep("steelblue4",3),
#     rep("forestgreen",2),
#    rep("tomato",1),
#     rep("forestgreen",2),
#     rep("tomato",1),
#     rep("forestgreen",3),
#     rep("tomato",4),
#     rep("forestgreen",1), alpha = 0.5), 
#     main="M-values",cex=2.5,cex.axis=2,cex.lab=2.8)
#dev.off()

#####################
## Design matrix and estimate disperssion ####
designSL <- model.matrix(~0+Group, data=targets)
design <- modelMatrixMeth(designSL)

#The first XX columns=>the sample coverage effects.
#The last three columns => the methylation levels (in logit units) in the three groups (BeforeRT, AfterRT and normal).

#y1 <- y[y$genes$Chr==1, ]    ### for testing only chrom 1
#y1 <- estimateDisp(y1, design=design, trend="none")
y1 <- estimateDisp(y, design=design, trend="none")

y1$common.dispersion
summary(y1$prior.df)

### save y1 file for easier load later:
saveRDS(y1, file=paste(output_dir,"edgeR_y1_merged_regions.rds", sep = "")) #### y1 file saved
#y1<-readRDS(file = "/home/shsoudi/methylation/full_cohort_analysis/out_edgeR_merged_regions/edgeR_y1_merged_regions.rds")

###################################################
## Differential methylation analysis at CpG loci ##
fit <- glmFit(y1, design)  ### fitting NB GLM for all CpGs

### assing groups for comparisons
#groups<-c("GroupN-GroupT1","GroupN-GroupT2","GroupN-GroupT3","GroupT1-GroupT2","GroupT1-GroupT3","GroupT2-GroupT3")
combinations<-c("GroupN-GroupT1","GroupN-GroupT2","GroupN-GroupT3","GroupT1-GroupT2","GroupT1-GroupT3","GroupT2-GroupT3")

#### Load Cosmic cancer genes
cancer_genes<-read.csv(file = "/home/shsoudi/methylation/full_cohort_analysis//Census_all.csv", header = TRUE, sep = ",")[,c("Gene.Symbol","Role.in.Cancer")]

### activate this part for paralelization
#n.cores <- parallel::detectCores() - 80
#
#my.cluster <- parallel::makeCluster(
#  n.cores, 
#  type = "PSOCK"
#)

#doParallel::registerDoParallel(cl = my.cluster)


## loop through each comparison type
for(cc in 1:length(combinations)){
#foreach(cc = 1:length(combinations))%dopar%{    

        #library(edgeR)
        #library(dplyr)
        #library(tidyr)
        #library(stringr)
        #library(org.Hs.eg.db)
        #library(EnhancedVolcano)
        #library(limma)  ### added

        focal_comb<-combinations[cc]
        print(focal_comb)
        
        #contr <- makeContrasts(GroupNvsB = GroupNormal-GroupBeforeRT, levels=design)
        contr <- makeContrasts(contrasts = focal_comb, levels = design)
        print(contr)

        lrt <- glmLRT(fit, contrast=contr)
        saveRDS(lrt, file=paste(output_dir,combinations[cc],"_merged_regions_lrt.rds",sep = ""))
        
        decide.tab<-data.frame(summary(decideTests(lrt)))
        write.table(decide.tab,file=paste(output_dir,combinations[cc],"_decide.tab_merged_regions_lrt.table",sep = ""),  quote = FALSE)   ### fixed

        topTags_5k<-topTags(lrt,500000)
        write.table(topTags_5k,file=paste(output_dir,combinations[cc],".topTags_5k_merged_regions_lrt.table",sep = ""))

        
        #The total number of DMRs in each direction at a FDR of 5% can be examined with decide Tests
        #summary(decideTests(lrt))
        #plotMD(lrt,bg.pch=0.5, bg.cex=0.5)
    
        #making a volcano plot #
        #pdf(file = paste(output_dir,"enhancedVolcao_edgeR",focal_comb,".pdf",sep = ""), width=14, height = 14)
         
        #aa<-lrt$table
        #EnhancedVolcano(aa,
        #         lab = rownames(aa),
        #         x = 'logFC',
        #         y = 'PValue',
        #         pCutoff = 0.001,
        #         FCcutoff = 5,
        #         pointSize = 1.0,
        #         labSize = 2.0)
        #dev.off()   
        
        #cleaned_comb <- gsub("Group", "", unlist(strsplit(focal_comb, "-")))
        #filtered_targets <- targets %>% filter(Group %in% cleaned_comb)

        InPromoter <- yall$genes$Distance >= -1000 & yall$genes$Distance <= 2000  ## counts in promoter regions
        yIP <- yall[InPromoter,,keep.lib.sizes=FALSE]
        ypr <- rowsum(yIP, yIP$genes$EntrezID, reorder=FALSE)  #count the total counts for each gene promoter
        ypr$genes$EntrezID <- NULL

        zero_genes <- rowSums(ypr$counts) == 0  #Check if you have all-zero rows and remove them
        summary(zero_genes)
        ypr <- ypr[!zero_genes, , keep.lib.sizes=FALSE]  

        summary(ypr$counts)
        any(is.na(ypr$counts))  # TRUE if NAs are present
        any(!is.finite(ypr$counts))  # TRUE if Inf values are present

        ypr$counts[is.na(ypr$counts)] <- 1    ## replace NA or Inf with a small positive value:
        ypr$counts[!is.finite(ypr$counts)] <- 1
        ypr <- calcNormFactors(ypr)
        

        Mepr <- ypr$counts[,Methylation=="Me"]
        Unpr <- ypr$counts[,Methylation=="Un"]
        Coveragepr <- Mepr + Unpr


        HasCoveragepr <- rowSums(Coveragepr >= 4) == round(length(targets$Sample)/2) 
        table(HasCoveragepr)
        #HasBothpr <- rowSums(Mepr) > 0 & rowSums(Unpr) > 0
        #table(HasCoveragepr, HasBothpr)

        ypr <- ypr[HasCoveragepr,,keep.lib.sizes=FALSE]
        

        TotalLibSizepr <- 0.5*ypr$samples$lib.size[Methylation=="Me"] +
        + 0.5*ypr$samples$lib.size[Methylation=="Un"]

        ypr$samples$lib.size <- rep(TotalLibSizepr, each=2)
        ypr$samples
 
       # check lib sizes are set sorrectly
        ypr$samples$lib.size[ypr$samples$lib.size == 0] <- 1  
        ypr <- calcNormFactors(ypr)
        
        ## Differential methylation in gene promoters 
        ypr <- estimateDisp(ypr, design, trend="none")  ### error comes from here
        ypr$common.dispersion

        fitpr <- glmFit(ypr, design)
        lrtpr <- glmLRT(fitpr, contrast=contr)
        saveRDS(lrtpr, file=paste(output_dir,combinations[cc],".promoter_merged_regions_lrt.rds",sep = ""))

        #topTags(lrtpr, n=20)
        aa_promoter<-lrtpr$table
        gg<-lrtpr$genes
        aa_promoter_genes<-cbind(aa_promoter,gg)
        rownames(aa_promoter_genes)<-aa_promoter_genes$Symbol
        write.table(aa_promoter_genes,file=paste(output_dir,combinations[cc],"_aa_promoter_genes_merged_regions_lrtpr.table",sep = ""),  quote = FALSE)   ### fixed


        cancer_promote<-aa_promoter_genes[aa_promoter_genes$Symbol%in%cancer_genes$Gene.Symbol,]
        cancer_promote$Symbol


        top<-data.frame(topTags(lrtpr, n=1000))
        write.table(top,file=paste(output_dir,combinations[cc],"_Promoter_topTags_1k_merged_regions_lrtpr.table",sep = ""),  quote = FALSE)   ### fixed

        top_cancer<-top[top$Symbol%in%cancer_genes$Gene.Symbol,]


    #pdf(file = paste(output_dir,"enhancedVolcao_promoter_annotatted_edgeR",focal_comb,".pdf",sep = ""), width=8, height = 9)
    #EnhancedVolcano(aa_promoter_genes,
    #             lab = rownames(aa_promoter_genes),
    #             x = 'logFC',
    #             y = 'PValue',
    #             pCutoff = 0.001,
    #             FCcutoff = 2,
    #             pointSize = 3.0,
    #             labSize = 5.0, 
    #             selectLab =top_cancer$Symbol,
    #             #selectLab =co_genes,
    #             #boxedLabels = TRUE,
    #             parseLabels = TRUE,
    #             drawConnectors = TRUE,
    #             widthConnectors = 0.75,
    #             colConnectors = 'black',
    #             max.overlaps=30)
    #dev.off()   




} ## end of group loop

#parallel::stopCluster(cl = my.cluster)
#co_genes<-c("CCNB1IP1","NFKB2","ALK","CNOT3","SRC","ZNF479","KRAS","TPM4","MYCN","KLF4","MYO5A","POLE","CHD2","ACVR1B","MUC16","CHD2","FUS","ARHGEF10L","NTRK2","CRTC3","SS18L1","TFPT","XPC","PML")

################
##### END ######
################

