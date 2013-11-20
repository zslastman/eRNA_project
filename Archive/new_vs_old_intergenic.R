


#New vs old intergenic numbers

sum(distanceToNearest(crmgrs,tss.gr)$distance>intergenic_dist)
tss.old<-'/g/furlong/project/11_Histone_ChIP_seq/BiTS-ChIP/analysis/used-TSSs.txt'
tss.old<-read.table(tss.old,header=T)
tss.old<-GRanges(paste0('chr',tss.old[,1]),IRanges(tss.old[,2],tss.old[,3]),strand=tss.old[,4],seqinfo=si)

transcripts.old<-'/g/furlong/project/11_Histone_ChIP_seq/BiTS-ChIP/analysis/TSS_description.txt'
transcripts.old<-read.table(transcripts.old,header=T)
transcripts.old<-GRanges(paste0('chr',transcripts.old$CHR),IRanges(transcripts.old$GeneStart,transcripts.old$GeneStop),
                         strand=transcripts.old$STRAND,seqinfo=si)

nontssnum<-sum(distanceToNearest(crmgrs,tss.gr)$distance>intergenic_dist)
oldnontssnum<-sum(distanceToNearest(crmgrs,tss.old)$distance>intergenic_dist)


sum(distanceToNearest(crmgrs,transcripts)$distance>intergenic_dist)
sum(distanceToNearest(crmgrs,transcripts.old)$distance>intergenic_dist)

transcripts.old<-'/g/furlong/project/11_Histone_ChIP_seq/BiTS-ChIP/analysis/'