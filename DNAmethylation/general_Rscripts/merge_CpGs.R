### merge CpGs per subject_tissue for downsteram analysis ###
rm(list = ls())
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
#library(foreach)
#library(doParallel)


## load metadata
meta<-read.csv(file = "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/coverage_files/meta_testing_updated.csv", header = TRUE)

## samples to exclude ##
meta_filt <- meta %>% 
        mutate(sampleid_status=paste(Final_ID,Tissue_Category, sep = "_")) 


samples<-meta_filt$Final_ID
subject_tissue<-unique(meta_filt$Subject_Tissue)


root_dir<-("/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/coverage_files/for_test_coveragefiles/top_1000")
files<-list.files(root_dir)
#dat.list <- vector(mode = "list", length = length(samples))


mainDir <- "/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/"
subDir <- "coverage_merged"
if (file.exists(subDir)){
   setwd(file.path(mainDir, subDir))
} else {
   dir.create(file.path(mainDir, subDir))
   setwd(file.path(mainDir, subDir))
   
}





for(ss in 1:length(subject_tissue)){  ### loop through each subject tissue

    focal_files<-files[grep(subject_tissue[ss], files)]
    dat.list <- vector(mode = "list", length = length(focal_files))

    
     
    for (kk in 1:length(focal_files)){   ## loop through each files belonging to the subject tissue to merge
    dat.list[[kk]] <- read.table(paste("/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/coverage_files/for_test_coveragefiles/top_1000/",focal_files[kk], sep = ""), header=F, col.names=c("chr", "start", "enc", "Methylation_Percentage","Count_Methylated","Count_Unmethylated"))

    if(nrow(dat.list[[kk]]) > 100){

        dat.list[[kk]]<-dat.list[[kk]][(dat.list[[kk]]$Count_Methylated  + dat.list[[kk]]$Count_Unmethylated >=5), ]
        dat.list[[kk]]$uniq_ID<-paste(dat.list[[kk]]$chr,dat.list[[kk]]$start , sep = "_")

       #fiter out scafolds
       # remove unplaced contigs from reference

       ChrNames <- c(1:19,"X","Y")
       dat.list[[kk]]<-dat.list[[kk]][dat.list[[kk]]$chr %in% ChrNames,]
       names(dat.list)[[kk]] <- focal_files[kk]

             } else {

             print(paste("nrow_",focal_files[kk], "_is not greater than threshold", sep = ""))
            }
    }    ## kk loop (focal files)
    
         dat.list <- dat.list[names(dat.list)!="NULL"] 

        if (length(dat.list)==1 ){


            aa<-as.data.frame(dat.list[1])
            colnames(aa)<-c("chr","start","end","Methylation_Percentage","Count_Methylated","Count_Unmethylated","ID")
            write.table(aa, file = paste(mainDir,subject_tissue[ss],"_merged_deduped.bismark.cov", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

        } else {
        merged <- Reduce(function(x, y) merge(x, y, by = "uniq_ID", all = TRUE), dat.list)
    

    #} ## kk loop (focal files)

    
    result <- data.frame(chr = NA, start = NA, end = NA , Methylation_Percentage = NA, Count_Methylated = NA, Count_Unmethylated = NA , ID = merged$uniq_ID)


    for (i in 1:nrow(merged)) {
         # Check if the 'ID' exists in all data frames (not NA)
         
         merged_focal<-merged[i,]
         filtNA_merged_focal<-Filter(function(x) !all(is.na(x)), merged_focal)
         result$chr[i]<-filtNA_merged_focal[,grep("chr",colnames(filtNA_merged_focal))][1]
         result$start[i]<-filtNA_merged_focal[,grep("start",colnames(filtNA_merged_focal))][1]
         result$end[i]<-filtNA_merged_focal[,grep("enc",colnames(filtNA_merged_focal))][1]


         foc_column_Methylation_Percentage<-merged_focal[,grep("Methylation_Percentage",colnames(merged_focal))]
         result$Methylation_Percentage[i] <- round(rowSums(foc_column_Methylation_Percentage, na.rm = TRUE)/rowSums(!is.na(foc_column_Methylation_Percentage)))

         foc_column_Count_Methylated<-merged_focal[,grep("Count_Methylated",colnames(merged_focal))]
         result$Count_Methylated[i] <- round(rowSums(foc_column_Count_Methylated, na.rm = TRUE)/rowSums(!is.na(foc_column_Count_Methylated)))

         foc_column_Count_Unmethylated<-merged_focal[,grep("Count_Unmethylated",colnames(merged_focal))]
         result$Count_Unmethylated[i] <- round(rowSums(foc_column_Count_Unmethylated, na.rm = TRUE)/rowSums(!is.na(foc_column_Count_Unmethylated)))
    
    } ## kk loop (focal files)

    df <- apply(result,2,as.character)
    write.table(df, file = paste("/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/analysis_for_grant/coverage_merged/",subject_tissue[ss],"_merged_deduped.bismark.cov", sep = ""),col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
 }
           

}



######################
######## END #########
######################



