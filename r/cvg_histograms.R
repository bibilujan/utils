#!/usr/bin/env Rscript

library(ggpubr)
library(ggplot2)
library(tidytext)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file)/nUsage:/nRscript cvg_histograms.R sample_bedtools_cvg.txt probe.bed", call.=FALSE)
} 

# Read and add column names to probe bed file
#bed<-read.delim(file="/Users/blujantoro/probeCoverageDev/GBS-2987/bed/LLDM.EXOME.TS.hg19.bed", sep="\t",header=F)
bed<-read.delim(file=args[1], sep="\t",header=F)
colnames(bed)<-c("chrom","start","stop","pool")
bed<-bed[order(bed$pool,bed$chrom,bed$start),]
rownames(bed)<-paste(bed$chrom,":",bed$start,"-",bed$stop,sep="")

# Read and add column names to bedtools histogram
#hist<-read.delim("/Users/blujantoro/probeCoverageDev/GBS-2987/cvg/LLDM_0019_Ln_P_PE_358_TS_210727_A00469_0191_AH7MVVDSX2_3_GCCTATCA-AATGGTCG_R1.fastq.gz.cvghist.txt",sep="\t",as.is=T,header=F)
hist<-read.delim(id,sep="\t",as.is=T,header=F)
colnames(hist)<-c("chrom","start","stop","pool","depth","bases","size","proportion")
hist<-hist[hist$chrom != "all",]
hist$interval<-paste(hist$chrom,":",hist$start,"-",hist$stop,sep="")

# Rename id 
#id <- "/Users/blujantoro/probeCoverageDev/GBS-2987/cvg/LLDM_0019_Ln_P_PE_358_TS_210727_A00469_0191_AH7MVVDSX2_3_GCCTATCA-AATGGTCG_R1.fastq.gz.cvghist.txt"
#id<-sub("^.*cvg/","",id)
id<-sub("^.*cvg/","",args[1])
id<-sub("_TS.*$","",id)
cat(id)
#hist.files=dir("cvg",full.names=T)
#id<-"cvg/LLDM_0019_Ln_P_PE_358_TS_210727_A00469_0191_AH7MVVDSX2_4_GCCTATCA-AATGGTCG_R1.fastq.gz.bam.cvg"
#id<-sub("cvg/","",args[1])
#id<-gsub("TS_210727_A00469_0191_AH7MVVDSX2_","",id,perl=T)
#id<-gsub("(_[[:alpha:]]+-[[:alpha:]]+_R1.*)","",id,perl=T)
#id<-gsub("_R1.*","",id,perl=T) #
#hist<-read.delim(file=args[1],sep="\t",as.is=T,header=F)
#id <- "LLDM_0095_Ln_P_PE_394_TS"
df.all<-NULL

# CALCULATE MEAN COVERAGE
mean_coverage<-by(hist,hist$interval,function(x){sum(x$depth*x$bases)/as.numeric(x$size[1])})
intervals<-as.vector(names(mean_coverage))
mean_coverage<-as.vector(mean_coverage)
df<-data.frame(id=id,interval=intervals,metric="cvg_mean",value=mean_coverage)
df.all<-rbind(df.all,df)

# CALCULATE NO COVERAGE INTERVALS
intervals<-unique(hist$interval)
no_coverage<-rep(0,length(intervals))
names(no_coverage)<-intervals
hist.nocoverage<-hist[hist$depth==0,]
no_coverage[hist.nocoverage$interval]<-hist.nocoverage$proportion*100
df<-data.frame(id=id,interval=intervals,no_coverage=no_coverage)
df<-data.frame(id=id,interval=intervals,metric="cvg0",value=no_coverage)
df.all<-rbind(df.all,df)

# CALCULATE PERCENT COVERED
hist$covered<-ifelse(hist$depth>0,1,0)
hist$proportion2<-hist$proportion*hist$covered
percent_covered <- by(hist,hist$interval,function(x){sum(x$proportion2)})
intervals<-as.vector(names(percent_covered))
percent_covered<-as.vector(percent_covered)
#df<-data.frame(id=id,interval=intervals,percent_covered=percent_covered)
df<-data.frame(id=id,interval=intervals,metric="pct_cvd",value=percent_covered)
df.all<-rbind(df.all,df)

### set interval as factor to order by intervals in the bed file
df.all$interval<-factor(df.all$interval,levels=rownames(bed))
df.all$pool<-bed[df.all$interval,]$pool

#save.image(file="image.RData")

############# Pool and Subsampled Exome PLOT
set1<-rownames(bed[bed$pool!="EXOME",])
set2<-sample(rownames(bed[bed$pool=="EXOME",]),450)
set<-c(set1,set2)
df.set<-df.all[df.all$interval %in% set,]

############# Percent coverage
cvg0<-df.all[df.all$metric=="cvg0",]
cvg_mean<-df.all[df.all$metric=="cvg_mean",]
percent_intervals_with_coverage<-aggregate(value~id + pool, data=cvg_mean,function(x){length(x[x>0])/length(x)*100})
#percent_intervals_withno_coverage<-aggregate(value~id + pool, data=cvg0,function(x){length(x[x>0])/length(x)*100})

g1<-ggplot(percent_intervals_with_coverage,aes(y=value,x=pool,col=pool)) + 
  geom_bar(stat="identity") + 
  #facet_wrap(~id,ncol=3) +
  theme(axis.text.x = element_blank(), legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.4, 'cm')) +
  labs(title=paste("                                                                                                  ", id,"\nPercent of Intervals with coverage", sep = "")) +
  xlab("pool") + 
  ylab("percent")

g2<-ggplot(df.set[df.set$metric=="pct_cvd",], aes(x=as.factor(pool), y=value, col=pool)) + 
  #geom_bar(stat="identity") + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  labs(title=paste ("Proportion of Interval covered", sep = "")) + xlab("pool") + ylab("percent")+
  #ylab("Proportion of interval covered") + 
  theme(axis.text.x = element_blank(), legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.4, 'cm')) #+
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+ 
#geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)
#geom_jitter(shape=16, position=position_jitter(0.2))


pct_cvd<-df.all[df.all$metric=="pct_cvd",]
percent_intervals_with_90_coverage<-aggregate(value~id + pool, data=pct_cvd,function(x){length(x[x>0.9])/length(x)*100})
g3<-ggplot(percent_intervals_with_90_coverage,aes(y=value,x=pool,col=pool)) + 
  geom_bar(stat="identity") + 
  #facet_wrap(~id,ncol=3) +
  theme(axis.text.x = element_blank(), legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.4, 'cm')) +
  labs(title=paste ("Percent of Intervals >90% covered", sep = "")) + 
  xlab("pool") + ylab("percent")
#ggsave(g3,file=paste("percent_90covered.png", sep=""),dev="png",height=30,width=15)

percent_intervals_with_75_coverage<-aggregate(value~id + pool, data=pct_cvd,function(x){length(x[x>0.75])/length(x)*100})
g4<-ggplot(percent_intervals_with_75_coverage,aes(y=value,x=pool,col=pool)) + 
  geom_bar(stat="identity") + 
  #facet_wrap(~id,ncol=3) +
  theme(axis.text.x = element_blank(), legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.4, 'cm')) +
  labs(title=paste ("Percent of Intervals >75% covered", sep = "")) + 
  xlab("pool") + ylab("percent")

gall <-ggarrange(g1, g2, g3, g4,
                 #labels = c("A", "B", "C"),
                 ncol = 1, nrow =4)
#gall <-annotate_figure(gall,top = text_grob(id, size = 10))
ggsave(gall,file=paste(id, "_pct_covered.png", sep = ""),dev="png",height=10,width=15)


############# Sorted coverage
df2<-df.all
#df2<-df2[df2$metric == "cvg_mean",]
df2$set<-ifelse(df2$pool=="EXOME","EXOME","LLDM_TARGETS")

subset1<-rownames(bed[bed$pool!="EXOME",])
subset2<-sample(rownames(bed[bed$pool=="EXOME",]),4000)
subset<-c(subset1,subset2)

df2<-df2[df2$interval %in% subset,]

g6<-ggplot(df2[df2$metric == "cvg_mean",],aes(x=reorder_within(interval,value,list(id,set)),y=value)) + geom_point() + 
  facet_wrap(~set,ncol=2,scales="free_x") +
  scale_y_log10()+
  theme(axis.text.x = element_blank()) +
  guides(x = "none") +
  labs(title=paste ("Sorted Mean Interval coverage-log scale", sep = "")) + xlab("interval") + ylab("mean coverage")
ggsave(g6,file=paste(id, "_sorted_mean_intv_cov.png", sep = ""),dev="png",height=10,width=15)


