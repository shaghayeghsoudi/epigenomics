#!/usr/bin/env Rscript

## load libraries
library(DSS)
library(bsseq)
library(dplyr)
library(ggplot2)
library(tidyr)
#library(gdata)
#library(kableExtra)
#library(knitr)
#library(parallel)
library(foreach)
library(doParallel)


## load metadata
meta<-read.csv(file = "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/pilot1/metadata/RRBS_methylation_pilot_experiment_samples.csv", header = TRUE)
meta$RT_status<-gsub("nopreopRT","beforeRT",meta$RT_status)


## samples to exclude ##
target<-c("SRC125-N8","SRC160-T1","SRC163-T1","SRC127-N7")
meta_filt <- meta %>% 
        filter(!(sample_id%in%target )) %>% 
        #mutate(sample_id=gsub("-","_",sample_id)) %>% 
        mutate(sampleid_status=paste(sample_id,RT_status, sep = "_")) 
        #dplyr::select(sample_id)

samples<-meta_filt$sample_id
#samples<-c("SRC168-N","SRC168-T4","SRC158-T8")  ## for testing only three samples 

##### Preprocessing
# read in methylation data in tab delimited format: <chromosome> <position> <count total> <count methylated>

root_dir<-("/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/pilot1/methylation_calls_pbat_from_nondirectional/CpG.counts/")
dat.list <- vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
        dat.list[[i]] <- read.table(paste(root_dir,samples[i], "_R1_val_1_bismark_bt2_pe.bismark.cov.gz.CpG.counts.txt", sep = ""), header=F, col.names=c("chr", "pos", "N", "X"))
        
        names(dat.list[[i]])<-c("chr","pos","N","X")

       ### filter cites with less than 5 reads
       dat.list[[i]]<-dat.list[[i]][dat.list[[i]]$N >5, ]

       #fiter out scafolds
       # remove unplaced contigs from reference

       ChrNames <- c(1:19,"X","Y")
       dat.list[[i]]<-dat.list[[i]][dat.list[[i]]$chr %in% ChrNames,]

}

### Set up BSobj ###
samples_equ<-meta_filt$sampleid_status
#samples_equ<-c("SRC168-N_Normal","SRC168-T4_beforeRT","SRC158-T8_afterRT")  ## for testing only three samples
samples_equ<-gsub("-","_",samples_equ)
names(dat.list)<-samples_equ

### create list of tables for different comparisons
combinations<-c("Normal|afterRT","Normal|beforeRT","beforeRT|afterRT")
   
#n.cores <- 1

### **create the cluster
#my.cluster <- parallel::makeCluster(
#  n.cores, 
#  type = "PSOCK"
#)

#register it to be used by %dopar%
#doParallel::registerDoParallel(cl = my.cluster)


### assign relevant directory and create a new directory for the final file 
mainDir <- "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/pilot1/downstream"
subDir <- "DSS-downstream"
if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
} else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    
}


for(ll in 1:length(combinations)){

#foreach(ll = 1:length(combinations))%dopar%{     

     #mainDir <- "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/pilot1/downstream"
     #subDir <- "test1"
     #if (file.exists(subDir)){
     #setwd(file.path(mainDir, subDir))
     #} else {
     #dir.create(file.path(mainDir, subDir))
     #setwd(file.path(mainDir, subDir))
     #}   

    #library(DSS)
    #library(bsseq)

      mm.list_focal<-dat.list[grep(combinations[ll], names(dat.list))]
      #BSobj <- makeBSseqData(mm.list_focal, names(mm.list_focal)) [1:100,]
      BSobj <- makeBSseqData(mm.list_focal, names(mm.list_focal)) 
     
      comp<-sub('.*\\_', '', names(mm.list_focal))

      comp<-unique(comp)

      groupA<-names(mm.list_focal)[grep(unique(comp)[1],names(mm.list_focal))]
      groupB<-names(mm.list_focal)[grep(unique(comp)[2],names(mm.list_focal))]
     
      #save(dmlTest, file="DML.RData")
      dmlTest <- DMLtest(BSobj, group1=groupA, group2=groupB, smoothing=TRUE)


      #save(dmlTest, file="DML.RData")
      #load("DML.RData")

      # call differential methylated CpG (DML)
      dmls <- callDML(dmlTest, delta=0.1, p.threshold=0.001)
      # write dmls to file results
      write.table(dmls, file=paste(subDir,"Differential_methylation_CpG_",comp[1],"_",comp[2],".txt", sep=""), col.names=T, row.names=F, quote=F)

      # call differentially methylated regions (DMR)
      dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=0.001)
      write.table(dmrs, file=paste(subDir,"Differential_methylation_Regions_",comp[1],"_",comp[2],".txt", sep=""), col.names=T, row.names=F, quote=F)
      #write.table(dmls, file=paste("Differential_methylation_Regions_",comp[1],"_",comp[2],".txt", sep=""), col.names=T, row.names=F, quote=F)

      #dmls_ordered<-dmrs[order(dmrs$fdr),]

     ## plot top DMRs
     for(jj in 1:200){
     pdf(paste(subDir,jj,"_",comp[1],"_",comp[2],".pdf", sep = ""), width = 20, height = 18)
     showOneDMR(dmrs[jj,], BSobj)
     dev.off()
     }


} ### foreach loop


#parallel::stopCluster(cl = my.cluster)



