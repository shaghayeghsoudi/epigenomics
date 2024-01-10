## this script process Bismark coverage output files

# Attributions
# original writer: Shaghayegh Soudi
# contributor: NA

## load required packages 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggforce)
library(reshape2)
library(ggpubr)


####### psot-processing quality check ####
####### load bismark coverage files
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


#ggplot(ww, aes(methylation_percent)) +
#  geom_histogram(aes(y = after_stat(density)), color = "#000000", fill = "#0099F8", bins = 100) +
#  geom_density(color = "#000000", fill = "#F85700", alpha = 0.6)

pdf("Histogram_density_of_CpG_percentage_methylation_pilot1.pdf", width = 10, height = 18)
#ggsave("",plot=)
 plot_hist<- ggplot(coverage_good, aes(methylation_percent)) +
  geom_histogram(aes(y = after_stat(density)), color = "#000000", fill = "#0099F8", bins = 20) +
  geom_density(color = "#000000", fill = "#F85700", alpha = 0.6) +facet_wrap( ~sample_id, ncol = 4)
print(plot_hist)
dev.off()


pdf("Histogram_of_CpG_percentage_methylation_pilot1.pdf", width = 10, height = 18)
 plot_hist2<-ggplot(coverage_good, aes(methylation_percent)) +
  geom_histogram(aes(y = after_stat(density)), color = "#000000", fill = "#0099F8", bins = 20) +
  facet_wrap( ~sample_id, ncol = 4)
print(plot_hist2)
dev.off()
 

#########################################
### make a count matrix ###


coverage_good_table <-coverage_good %>% 
     mutate(pos_id=paste())


#########################################
### plot methylation ratio per sample ###
# create methylation stats file
meta<-read.delim(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/RRBS/metadata/RRBS_methylation_pilot_experiment_samples.csv", header = TRUE, sep = ",")[,c("sample_id","tissue","RT_status")]
names(coverage_good)<-c("chr", "start", "end", "methylation_ratio", "count_methylated", "count_unmethylated","sample_id")
overage_good_meta<-merge(coverage_good,meta, by.x = "sample_id", by.y = "sample_id")
        


overage_good_meta_d1 <- overage_good_meta %>% 
         mutate(depth=count_methylated+count_unmethylated)%>% 
         filter(depth  >=1) 
        
mean_d1<-overage_good_meta_d1 %>% 
         group_by(sample_id) %>% 
         summarise_at(vars(methylation_ratio), list(name = mean))


#### keep CpGs with at least depth of 5
overage_good_meta_d5 <- overage_good_meta %>% 
         mutate(depth=count_methylated+count_unmethylated)%>% 
         filter(depth  >=5)

mean_d5<-overage_good_meta_d5 %>% 
     group_by(sample_id) %>% 
     summarise(mean(methylation_ratio)
 )


#########################################
###### load methylation stats file ######

methyl_stats<-read.csv(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/RRBS/downstream/methylation_stats_pilot1.csv", header = TRUE, sep = ",") 
methyl_stats_meta<-merge(methyl_stats,meta, by.x = "Sample_ID", by.y = "sample_id")

methyl_stats_meta$RT_status<-gsub("nopreopRT","beforeRT",methyl_stats_meta$RT_status)



methyl_stats_meta$tissue<-gsub("FFPE_Tumour","FFPE",
                          gsub("FFPE_Normal","FFPE",
                          gsub("To_Confirm","Fresh_Frozen_Tumor",methyl_stats_meta$tissue)))



# preliminary count bar chart of CpG count and
d1<-ggplot(methyl_stats_meta, aes(x=Sample_ID, y=Mean_MethRatio_d1,fill = tissue)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~RT_status, scales="free", space = "free") +
  ggtitle("percentage of methylation ratio minimum depth 1")
 
d5<-ggplot(methyl_stats_meta, aes(x=Sample_ID, y=Mean_MethRatio_d5,fill = tissue)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~RT_status, scales="free", space = "free") +
  ggtitle("percentage of methylation ratio minimum depth 5")

pdf(file = "~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/methylation/RRBS/downstream/barchart_CpG_methylation_ratio_d1_d5.pdf",width = 7, height=10 )
plot_both<-ggarrange(d1,d5, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
print(plot_both)   
dev.off()       

##############################
### make stacked bar plot ####
#d1_data<- methyl_stats_meta %>% 
#    select(Sample_ID,CpG.num_d1,Mean_MethRatio_d1,tissue, RT_status) %>% 
#    mutate(coverage=rep("d1")) 
#names(d1_data)[2:3]<-c("CpG.num", "Mean_MethRatio")   
#
#
#d5_data<- methyl_stats_meta %>% 
#    select(Sample_ID,CpG.num_d5,Mean_MethRatio_d5,tissue, RT_status) %>% 
#   mutate(coverage=rep("d5")) 
#names(d5_data)[2:3]<-c("CpG.num", "Mean_MethRatio")   
#
#both<-rbind(d1_data,d5_data)


