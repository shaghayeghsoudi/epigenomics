############################################
##### process RRBS outputs with EdgeR ######
############################################
rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggforce)
library(reshape2)
library(ggpubr)
library(edgeR)
library(EnhancedVolcano)
#https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/RRBS/downstream/edgeR/pilot1")
targets <- read.delim(paste("targets.csv", sep = ""), stringsAsFactors=FALSE, sep = ",")

files <- paste("coverage_files/",targets$File, sep = "")
yall <- readBismark2DGE(files, sample.names=targets$Sample)

dim(yall)

#We remove the mitochondrial genes as they are usually of less interest, and also unassembled contigs
table(yall$genes$Chr)

#yall <- yall[yall$genes$Chr!="MT", ]


ChrNames <- c(1:19,"X","Y")

yall<-yall[yall$genes$Chr%in%ChrNames, ]

#yall$genes$Chr <- factor(yall$genes$Chr, levels=ChrNames)
o <- order(yall$genes$Chr, yall$genes$Locus)
yall <- yall[o,]


#We now annotate the CpG loci with the identity of the nearest gene
TSS <- nearestTSS(yall$genes$Chr, yall$genes$Locus, species="Hs")
yall$genes$EntrezID <- TSS$gene_id
yall$genes$Symbol <- TSS$symbol
yall$genes$Strand <- TSS$strand
yall$genes$Distance <- TSS$distance
yall$genes$Width <- TSS$width
head(yall$genes)


#################################
## Filtering and normalization ##

### Note1: CpG loci that have low coverage are removed prior to downstream analysis as they provide little information for assessing methylation levels
#summing up the counts of methylated and unmethylated reads to get the total read coverage at each CpG site for each sample
Methylation <- gl(2,1,ncol(yall), labels=c("Me","Un"))
Me <- yall$counts[, Methylation=="Me"]
Un <- yall$counts[, Methylation=="Un"]
Coverage <- Me + Un
head(Coverage)


## As a conservative rule of thumb, we require a CpG site to have a total count (both methylated and unmethylated) of at least 8 in every sample before it is considered in the study.
#nsamples<-length(colnames(Coverage))
HasCoverage <- rowSums(Coverage >= 5) == 20   ### 20 total number of samples, seems very conservative
#table(HasCoverage)

#We also filter out CpGs that are never methylated or always methylated as they provide no information about differential methylation:
#HasBoth <- rowSums(Me) > 0 & rowSums(Un) > 0
#table(HasCoverage, HasBoth)

#The DGEList object is subsetted to retain only the non-filtered loci:
# y_coverage_both <- yall[HasCoverage & HasBoth,, keep.lib.sizes=FALSE] 
y <- yall[HasCoverage,, keep.lib.sizes=FALSE]


TotalLibSize <- y$samples$lib.size[Methylation=="Me"] +
+ y$samples$lib.size[Methylation=="Un"]

y$samples$lib.size <- rep(TotalLibSize, each=2)
y$samples


### NAs??
#i <- is.na(y$genes$EntrezID)
#y <- y[!i, ]

## Data exploration
Me <- y$counts[, Methylation=="Me"]
Un <- y$counts[, Methylation=="Un"]
M <- log2(Me + 2) - log2(Un + 2)
colnames(M) <- targets$Sample

## sample and RT based
#pdf("M-plot-sample-RTbased-colored.pdf", height = 10, width = 10)
#plotMDS(M, col=c(rep("black",1), rep("firebrick",2),rep("blue",1),rep("green",2),rep("orange",2),rep("chartreuse4",1),rep("chocolate4",2),rep("chartreuse4",1),rep("royalblue4",3),rep("chartreuse4",1),rep("chartreuse4",1),rep("chartreuse4",1),rep("chartreuse4",1),rep("steelblue4",1)), main="M-values",cex=0.8)
#dev.off()

## RT based ##
pdf("M-plot-RTbased-colored.pdf", height = 10, width = 10)
plotMDS(M, col=c(rep("springgreen4",3), rep("steelblue4",3),rep("springgreen4",2),rep("tomato4",1),rep("springgreen4",2),rep("tomato4",1),rep("springgreen4",3),rep("tomato4",4),rep("springgreen4",1)), main="M-values",cex=0.8)
dev.off()


####################
##Design matrix ####
#One aim of this study is to identify differentially methylated (DM) loci between the different cell populations.
designSL <- model.matrix(~0+Group, data=targets)
design <- modelMatrixMeth(designSL)

#The first XX columns=>the sample coverage effects.
#The last three columns => the methylation levels (in logit units) in the three groups (BeforeRT, AfterRT and normal).


###Dispersion estimation
#y1 <- y[y$genes$Chr==1, ]    ### for testing only chrom 1
#y1 <- estimateDisp(y1, design=design, trend="none")
y1 <- estimateDisp(y, design=design, trend="none")

y1$common.dispersion
summary(y1$prior.df)

### save y1 file for easier load later:
saveRDS(y1, file="edgeR.y1.rds")  #### y1 file saved


###################################################
## Differential methylation analysis at CpG loci ##
fit <- glmFit(y1, design)

### assing groups for comparisons
groups<-c("GroupBeforeRT-GroupNormal","GroupafterRT-GroupNormal","GroupBeforeRT-GroupafterRT")


## loop through each comparison type
for(jj in 1:length(group)){

        contr <- makeContrasts(GroupNvsB = groups[jj], levels=design)
        #contr <- makeContrasts(GroupNvsB = GroupNormal-GroupBeforeRT, levels=design)

        lrt <- glmLRT(fit, contrast=contr)
        saveRDS(lrt, file=paste(groups[jj],".lrt.rds",sep = ""))
        topTags(lrt)

        topTags_50k<-topTags(lrt,50000)
        write.table(topTags_50k,file=paste(groups[jj],".topTags_50k.table",sep = ""))

        #The total number of DMRs in each direction at a FDR of 5% can be examined with decide Tests
        #summary(decideTests(lrt))
        #plotMD(lrt,bg.pch=0.5, bg.cex=0.5)
    

        #### custom visulaization #####
        #making a volcano plot #
        pdf(file = paste("enhancedVolcao_edgeR",groups[jj],".pdf",sep = ""), width=14, height = 14)
        aa<-lrt$table
        EnhancedVolcano(aa,
                 lab = rownames(aa),
                 x = 'logFC',
                 y = 'PValue',
                 pCutoff = 0.001,
                 FCcutoff = 5,
                 pointSize = 1.0,
                 labSize = 2.0)
        dev.off()   

    
        ## Summarizing counts in promoter regions ##
        InPromoter <- yall$genes$Distance >= -1000 & yall$genes$Distance <= 2000
        yIP <- yall[InPromoter,,keep.lib.sizes=FALSE]

        #compute the total counts for each gene promoter:
        ypr <- rowsum(yIP, yIP$genes$EntrezID, reorder=FALSE)
        ypr$genes$EntrezID <- NULL


        Mepr <- ypr$counts[,Methylation=="Me"]
        Unpr <- ypr$counts[,Methylation=="Un"]
        Coveragepr <- Mepr + Unpr


        HasCoveragepr <- rowSums(Coveragepr >= 5) == 20
        HasBothpr <- rowSums(Mepr) > 0 & rowSums(Unpr) > 0
        table(HasCoveragepr, HasBothpr)
        ypr <- ypr[HasCoveragepr & HasBothpr,,keep.lib.sizes=FALSE]

        TotalLibSizepr <- 0.5*ypr$samples$lib.size[Methylation=="Me"] +
        + 0.5*ypr$samples$lib.size[Methylation=="Un"]
        ypr$samples$lib.size <- rep(TotalLibSizepr, each=2)
        ypr$samples
 
        ## Differential methylation in gene promoters ##
        ypr <- estimateDisp(ypr, design, trend="none")
        ypr$common.dispersion

        fitpr <- glmFit(ypr, design)
        lrtpr <- glmLRT(fitpr, contrast=contr)
        saveRDS(lrtpr, file=paste(groups[jj],".promoter.lrt.rds",sep = ""))

        topTags(lrtpr, n=20)

        aa_promoter<-lrtpr$table

    pdf(file = paste("enhancedVolcao_promoter_edgeR",groups[jj],".pdf",sep = ""), width=14, height = 14)
    EnhancedVolcano(aa_promoter,
                 lab = rownames(aa_promoter),
                 x = 'logFC',
                 y = 'PValue',
                 pCutoff = 0.001,
                 FCcutoff = 2,
                 pointSize = 2.0,
                 labSize = 2.0)
    dev.off()   




} ## end of group loop



################################
####### downstream edgeR #######
################################


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library("colorRamp2")

### load Bismark coverage files
list_coverage<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/RRBS/bismark-calling/trimmed_nondirectional/clipped/coverage", pattern="*.cov", full.names = TRUE)
attackStats_COV <- lapply(list_coverage,function(x) {
     read.delim(x, header=FALSE, sep = "\t")
     })


for (i in 1:length(attackStats_COV )){
    attackStats_COV[[i]]<-cbind(attackStats_COV[[i]],list_coverage[i])
    }
aa_montecarlo<- do.call("rbind", attackStats_COV)      


names(aa_montecarlo)[7]<-"path"
coverage_good<-aa_montecarlo %>% 
    mutate(sample_id=sub('.*/\\s*','', gsub("_R1_val_1_bismark_bt2_pe.bismark.cov","",path)),)%>% 
    select(V1,V2,V3,V4,V5,V6,sample_id)

names(coverage_good)<-c("chr", "start", "end", "methylation_percent", "count_methylated", "count_unmethylated","sample_id")
coverage_good$pos_id<-paste(coverage_good$chr, coverage_good$start , sep = "-")

### coverage file within target files
coverage_target<-coverage_good[coverage_good$sample_id %in% targets$Sample,]

### save it as rds for easier re-loading
saveRDS(coverage_target, file="Bismark_coverage_edgeR_target.rds")

coverage_target<-readRDS(file ="Bismark_coverage_edgeR_target.rds")
###############################################
## laod metadata for the methylation analysis
meta<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/RRBS/metadata/RRBS_methylation_pilot_experiment_samples.csv", header = TRUE, sep = ",")[,c("sample_id","tissue","RT_status")]
names(coverage_good)<-c("chr", "start", "end", "methylation_ratio", "count_methylated", "count_unmethylated","sample_id")
coverage_good_meta<-merge(coverage_good,meta, by.x = "sample_id", by.y = "sample_id")
coverage_good_meta$RT_status<-gsub("nopreopRT","beforeRT",coverage_good_meta$RT_status)       




pair_types<-unique(dss_all_pairgood$pair_type)



######################
##### for TESTING ####

#samples<-c("SRC159-T1","SRC159-T2","SRC161-T2","SRC161-T3","SRC161-T4","SRC163-N","SRC168-N","SRC172-N")
samples_beforeafter<-c("SRC125-T1","SRC127-T1","SRC127-T3","SRC159-T1","SRC159-T2","SRC160-T2","SRC160-T5","SRC161-T2","SRC161-T3","SRC161-T4","SRC168-T4","SRC127-T7","SRC158-T5","SRC158-T8")


cov_samples<-coverage_target[coverage_target$sample_id%in%samples_beforeafter,]


ltr<-read.delim(file = "GroupBeforeRT-GroupafterRT.topTags_50k.table", header = TRUE, sep = "")
ltr_top<-ltr[ltr$FDR<0.05,]

cov_ltr_positions<-cov_samples[cov_samples$pos_id %in%rownames(ltr_top),]


positions<-unique(cov_ltr_positions$pos_id)


out_res<-NULL
for (kk in 1:length(positions)){

    #for (kk in 1:2000){

    cov_ltr_positions_focal<-data.frame(t(cov_ltr_positions[cov_ltr_positions$pos_id==positions[kk],"methylation_percent"]))
    samples_true<-cov_ltr_positions[cov_ltr_positions$pos_id==positions[kk],"sample_id"]

    ### assign hyper vs. hypo methylation
    ltr_focal<-ltr[rownames(ltr)==positions[kk],]
    ltr_focal$direction<-ifelse(ltr_focal$logFC>0 , "hyper", "hypo")


    if(ncol(cov_ltr_positions_focal)==14){    #### length should be adjusted

        rownames(cov_ltr_positions_focal)<-positions[kk]
        cov_ltr_positions_focal$direction<-ltr_focal$direction
        cov_ltr_positions_focal$logFC<-ltr_focal$logFC 
        out_res<-rbind(out_res,cov_ltr_positions_focal)

    } ### if statement loop 


}


out_res_direction<-out_res$direction
out_res_FC<-out_res$logFC
out_res_sample<-out_res[1:(ncol(out_res)-2)] 
#colnames(out_res_sample) <-c("SRC159-T1","SRC159-T2","SRC161-T2","SRC161-T3","SRC161-T4","SRC163-N","SRC168-N","SRC172-N")
colnames(out_res_sample)<-samples_true
out_res_mat<-as.matrix(out_res_sample)


column_tree = hclust(dist(t(out_res_mat)))

length_cpg<-2000
out_res_mat_1000<-out_res_mat[1:length_cpg,]
rownames(out_res_mat_1000)<-NULL


out_res_direction_1000<-out_res_direction[1:length_cpg]
out_res_FC_1000<-out_res_FC[1:length_cpg]


#type<-c("Tumour","Tumour","Tumour","Tumour","Tumour","Normal","Normal","Normal")
type<-c("BeforeRT","BeforeRT","BeforeRT","afterRT","afterRT","afterRT","BeforeRT","BeforeRT","BeforeRT","BeforeRT","BeforeRT","BeforeRT","BeforeRT","BeforeRT")


#ha = HeatmapAnnotation(df = data.frame(type = type), 
#    col = list(type = c("Tumour" = "mediumseagreen", "Normal" = "mediumslateblue")))
#ha2 = HeatmapAnnotation(df = data.frame(type = type), 
#    col = list(type = c("Tumour" = "mediumseagreen", "Normal" = "mediumslateblue")), 
#    show_legend = FALSE)


ha = HeatmapAnnotation(df = data.frame(type = type), 
    col = list(type = c("BeforeRT" = "mediumseagreen", "afterRT" = "mediumslateblue")))
ha2 = HeatmapAnnotation(df = data.frame(type = type), 
    col = list(type = c("BeforeRT" = "mediumseagreen", "afterRT" = "mediumslateblue")), 
    show_legend = FALSE)

#Heatmap(out_res_mat_1000)


Heatmap(out_res_mat_1000, name = "methylation",
                          cluster_columns = column_tree, 
                          column_dend_reorder = FALSE,
                          top_annotation = ha,
                          km = 5, 
                          column_title = "Methylation", 
                          column_title_gp = gpar(fontsize = 2),
                          row_title_gp = gpar(fontsize = 10)) +
                          Heatmap(out_res_direction_1000, name = "direction", col = c("hyper" = "orange2", "hypo" = "#4803bf")) +
                          Heatmap(out_res_FC_1000, name = "log(FC)", col = colorRamp2(c(-12, 0, 12), c("seagreen4", "white", "tomato3")))


###########
### END ###
###########

