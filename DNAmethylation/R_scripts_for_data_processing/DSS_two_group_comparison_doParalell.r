#!/usr/bin/env Rscript

## load libraries
library(DSS)
library(bsseq)
library(dplyr)
#library(ggplot2)
#library(tidyr)
#library(gdata)
#library(kableExtra)
#library(knitr)
library(parallel)
library(foreach)
#library(doParallel)


root_dir<-("/home/shsoudi/methylation/full_cohort_analysis/DSS_input_filt_5reads_merged_per_subject/processed_file_merged_by_addingupp/")
#root_dir<-("/home/shsoudi/methylation/full_cohort_analysis/test/coverage_for_DSS_testing/dss_input_files/")

file_names<-list.files(root_dir)
#dat.list <- vector(mode = "list", length = length(samples))
# Extract sample names based on conditions

#sample_names <- data.frame("sample_name"=sapply(file_names, function(file) {
#  if (grepl("_merged\\.cov\\.DSS_input$", file)) {
#    # Remove '_merged.cov.DSS_input' if it exists
#    gsub("_merged\\.cov\\.DSS_input$", "", file)
#  } else {
#    # Remove '.cov.DSS_input'
#    gsub("\\.cov\\.DSS_input$", "", file)
#  }
#}))

#samples<-sample_names$sample_name
samples<-gsub("_merged_DSS_input.txt","",file_names)
#dat.list <- vector(mode = "list", length = length(samples))
combinations<-c("T1|T3","N|T3","N|T1","N|T2","T1|T2","T2|T3")

### create output dirtory if it does not exist
output_dir <- "/home/shsoudi/methylation/full_cohort_analysis/test_out/"

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Directory created:", output_dir, "\n")
} else {
  cat("Directory already exists:", output_dir, "\n")
}



n.cores <- parallel::detectCores() - 80

### **create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#for (cc in 1:length(combinations)){  ### loop through each combination

foreach(cc = 1:length(combinations))%dopar%{    

    #if (!dir.exists(output_dir)) {
    #  dir.create(output_dir, recursive = TRUE)
    #  cat("Directory created:", output_dir, "\n")
    #} else {
    #  cat("Directory already exists:", output_dir, "\n")
    #}

    library(DSS)
    library(bsseq)

    focal_comb<-combinations[cc]
    focal_files <- grep(focal_comb, file_names, value = TRUE)
    dat.list <- vector(mode = "list", length = length(focal_files))

  for (i in 1:length(focal_files)){  ### lopp trough each file in the focal combination

        dat.list[[i]] <- read.table(paste(root_dir,focal_files[i], sep = ""), header=F, col.names=c("chr", "pos", "N", "X","id"))
        #dat.list[[i]] <- read.table(paste(root_dir,focal_files[i], sep = ""), header=F, col.names=c("chr", "pos", "N", "X"))

        dat.list[[i]]<-dat.list[[i]][,c("chr", "pos", "N", "X")]  ##### NOTE : check and activate this if the column 5 is present and if not deactivate this line

        if(nrow(dat.list[[i]]) > 500000){     #### adjust threshold
        
             names(dat.list[[i]])<-c("chr","pos","N","X")
             ### filter cites with less than 5 reads
             dat.list[[i]]<-dat.list[[i]][dat.list[[i]]$N >=5, ]

             #fiter out scafolds
             #remove unplaced contigs from reference

             ChrNames <- c(1:19,"X","Y")
             dat.list[[i]]<-dat.list[[i]][dat.list[[i]]$chr %in% ChrNames,]
             names(dat.list)[[i]] <- sub("\\.cov_DSS_input$", "", focal_files[i])

         } else {

             print(paste("nrow_",samples[i], "_is not greater than threshold", sep = ""))
        }

   } ## "i" loop" (each file in the focal combination)

    dat.list <- dat.list[!is.na(names(dat.list))] 


    BSobj <- makeBSseqData(dat.list, names(dat.list)) 
    comp<-strsplit(combinations[cc], "\\|")

    
      groupA<-names(dat.list)[grep(comp[[1]][1],names(dat.list))]
      groupB<-names(dat.list)[grep(comp[[1]][2],names(dat.list))]
     
      #save(dmlTest, file="DML.RData")
      dmlTest <- DMLtest(BSobj, group1=groupA, group2=groupB, smoothing=TRUE,smoothing.span = 500)

      #save(dmlTest, file="DML.RData")
      #load("DML.RData")

      # call differential methylated CpG (DML)
      dmls <- callDML(dmlTest, delta=0.1, p.threshold=0.05)
      # write dmls to file results
      write.table(dmls, file=paste(output_dir,"Differential_methylation_CpG_",comp[[1]][1],"_",comp[[1]][2],".txt", sep=""), col.names=T, row.names=F, quote=F)


      # call differentially methylated regions (DMR)
      dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=0.05)
      write.table(dmrs, file=paste(output_dir,"Differential_methylation_Regions_",comp[[1]][1],"_",comp[[1]][2],".txt", sep=""), col.names=T, row.names=F, quote=F)
      #write.table(dmls, file=paste("Differential_methylation_Regions_",comp[1],"_",comp[2],".txt", sep=""), col.names=T, row.names=F, quote=F)

      #dmls_ordered<-dmrs[order(dmrs$fdr),]

     ## plot top DMRs
     #for(jj in 1:50){
     #pdf(paste(subDir,jj,"_",comp[[1]][1],"_",comp[[1]][2],".pdf", sep = ""), width = 20, height = 18)
     #showOneDMR(dmrs[jj,], BSobj)
     #dev.off()
     #}


} ### end of combination loop




#parallel::stopCluster(cl = my.cluster)


