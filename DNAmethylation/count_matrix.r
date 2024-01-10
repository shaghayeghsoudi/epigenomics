
library(dplyr)

## load metadata table
setwd("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/RRBS/downstream/count_matrices")


meta<-read.csv(file = "test/metadata.table", header = TRUE, sep = "\t") 
meta_good <- meta %>% 
  mutate(RT_status = str_replace(RT_status, "nopreopRT", "beforeRT"))


####### load bismark coverage files
list_coverage<-list.files("test", pattern="*.cov", full.names = TRUE)

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

#########################################
### make a count matrix ###

### exclude FFPES and low coverage samples
target<-c("SRC125-N8","SRC160-T1","SRC163-T1","SRC127-N7")   ### low coverage samples to be dropped

#coverage_good %>% inner_join(coverage_good, meta_good, by = "sample_id")
coverage_good_table <-coverage_good %>% 
     mutate(pos_id=paste(chr,start, sep ="_")) %>% 
     filter(!(sample_id%in%target))%>% 
     mutate(coverage=count_methylated+count_unmethylated)


## make categories for comparisons 
cats<-c("beforeRT-afterRT","Normal-beforeRT")

for (ll in 1:length(cats)){

  cat_foc<-cats[ll]
  cat1<-as.character(sapply(strsplit(cat_foc, "-"), "[", c(1,2,3)))
  meta_cat1<-meta_good[meta_good$RT_status%in%cat1,"sample_id"]  ### sample we need for each category

  coverage_cat1<-coverage_good_table[coverage_good_table$sample_id%in% meta_cat1, ]
  ChrNames <- c(1:22,"X","Y")

}


#######################
###### up to here #####
## keep only pre-and postRT samples
#coverage_good_tableRT<-coverage_good_table[grep("T", coverage_good_table$sample_id), ]  ## 14 samples only
coverage_good_tableRT<-coverage_good_table


ChrNames <- c(1:22,"X","Y")
#ChrNames <- c(20:22) 
coverage_good_tableRT<-coverage_good_tableRT[coverage_good_tableRT$chr%in%ChrNames, ]


### keep positions shared between all tumour samples
count_14<-table(coverage_good_tableRT$pos_id) 
position14<-subset(names(count_14), count_14 == 17)

coverage_good_14<-coverage_good_tableRT[coverage_good_tableRT$pos_id%in%position14,]


out_res<-NULL
for (i in 1:length(ChrNames)){

  coverage_good_14_chrom<-coverage_good_14[coverage_good_14$chr==ChrNames[i],]
  positions<-unique(coverage_good_14_chrom$pos_id)

  out_res_chrom<-NULL
  for (kk in 1:length(positions)){

    coverage_good_14_final<-data.frame(t(coverage_good_14_chrom[coverage_good_14_chrom$pos_id==positions[kk],c("coverage","sample_id")]))
    colnames(coverage_good_14_final)<-coverage_good_14_final[2,]
    focal_pos_good<-coverage_good_14_final[1,]
    rownames(focal_pos_good)<-positions[kk]
    out_res_chrom<-rbind(focal_pos_good,out_res_chrom)


    #coverage<-data.frame(t(coverage_good_14_chrom[coverage_good_14_chrom$pos_id==positions[kk],c("count_methylated","sample_id")]))
    #colnames(coverage)<-coverage[2,]
    #focal_pos_good_cov<-coverage[1,]
    #rownames(focal_pos_good_cov)<-positions[kk]
    #out_res_chrom<-rbind(focal_pos_good,out_res_chrom)

  }

write.table(out_res_chrom, file = paste("methylation_coverage_normal_beforeRT_chrom",ChrNames[i],".table", sep = ""),quote = FALSE, sep = "\t")
}


