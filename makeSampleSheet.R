library(dplyr)
library(tidyr)

serverPath="/Volumes/meister.data"
workDir=paste0(serverPath,"/jsemple/20260728_ma_10x_multi_PMW941_1/ahringer_pipeline_claude")
setwd(workDir)

df<-read.delim(paste0(workDir,"/config/samples.tsv"))
fq<-read.delim(paste0(workDir,"/fastqfiles.txt"),header=F)

gex<-data.frame(fastq=fq[grepl("_GEX_",fq$V1),])
gex$gex_fastq_dir<-dirname(gex$fastq)
gex$gex_fastq_prefix<-"TEV0_GEX"
gex$batch<-"batch1"
gex$stage_group<-"l3"
gex$idr_rep<-1
gex$utr_ext<-"yes"


gex$fastq<-NULL
gex<-unique(gex)

atac<-data.frame(fastq=fq[grepl("_ATAC_",fq$V1),])
atac$atac_fastq_dir<-dirname(atac$fastq)
atac$atac_fastq_prefix<-gsub("_S.*_L.*fastq.gz","",basename(atac$fastq))
atac$batch<-"batch1"
atac$stage_group<-"l3"
atac$idr_rep<-1
atac$utr_ext<-"yes"

atac$fastq<-NULL
atac<-unique(atac)
atac$atac_fastq_prefix<-paste0(atac$atac_fastq_prefix,collapse=",")
atac<-unique(atac)

ss<-inner_join(gex, atac)
ss$sample_id<-"PMW941_1"
ss$barhop<-NA
ss$ambient<-NA

ss<-ss[,colnames(df)]

ss

write.table(ss,paste0(workDir,"/config/samples.tsv"),sep="\t",row.names=F,quote=F)
