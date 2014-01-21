#script to load RNAseq data, produce binned coverage wigs etc.
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')

# #Specify the bam files
# rnabams<-list('data/RNAseq/U6_PE2_thout/accepted_hits.nodup.sort.bam',
# 'data/RNAseq/mesoRNA_jul2013/U19_nodup.bam',
# 'data/RNAseq/S3_PE_thout/accepted_hits.nodup.sort.bam',
# 'data/RNAseq/S6_PE_thout/accepted_hits.nodup.sort.bam',
# 'data/RNAseq/mesoRNA_jul2013/S19_nodup.bam',
# '/g/furlong/Harnett/RNASEQ_download_20aug/S12_g1/accepted_hits_sorted.bam',
# '/g/furlong/Harnett/RNASEQ_download_20aug/S2/accepted_hits_sorted.bam')


#function to add up list of lists of rles.
sum.list.of.rle.lists <- function (sll) {
  #this takes in a list of lists of simplerlelists
  #(e.g. one of the structure - acc,strand,chr)
  summed<-sll[[1]]
  expect_SRL(summed[[1]])
  for(n in 2:length(sll)){
    for(strand in names(summed)){
     summed[[strand]]<-sll[[n]][[strand]]+summed[[strand]]
    }
  }
  return(summed)
}

wig2rle<-function(wigfile,strandsplit=F){
  out=import(wigfile)
  out=as(out,'GRanges')
  out=keepSeqlevels(out,chrs.keep)
  seqinfo(out)=si
  if(strandsplit){
    out=list(
      pos=coverage(out[out$score>=0],weight='score'),
      neg=coverage(out[out$score<=0],weight='score'))
    sapply(out,abs)
  }else{
    out$score=abs(out$score)
    return(coverage(out,weight='score'))
  }
}



################################################################################
#1 Get our RNAseq data
rnabams<-list.files(recursive=T,full.names=T,'data/RNAseq_links/',pattern='.bam$')
is.meso<-grepl('meso',rnabams)
#get coverage rles
rna.seq<-mclapply(mc.cores=10,rnabams,function(bam)bam2coverage(bam,doshift=F,doresize=F,stranded=T,is.ProperPair=NA,is.FirstMateRead=F,use.size=F,maxsize=NA))
#merge the embryo and mesoderm rles
rna.seq = c(meso = sum.list.of.rle.lists(rna.seq[is.meso]) ,embryo = sum.list.of.rle.lists(rna.seq[!is.meso]) )
#now produce the coverage wigs
export(rna.seq[[embryo]]$pos,'data/RNAseq_links/secondmate/all_rnaseq_emb_6.8.pos.bw')
export(rna.seq[[embryo]]$neg,'data/RNAseq_links/secondmate/all_rnaseq_emb_6.8.neg.bw')
export(rna.seq[[meso]]$pos,'data/RNAseq_links/secondmate/all_rnaseq_meso_6.8.pos.bw')
export(rna.seq[[meso]]$neg,'data/RNAseq_links/secondmate/all_rnaseq_meso_6.8.neg.bw')
save(rna.seq,file='data/objects/rna.seq.object.R')


################################################################################
#2 Now get the coverage wigs from teh Celniker total RNAseq experiments
c.rnawigs=list.files('data/graveley_totalrnaseq/',full.names=T,pattern='wig',recursive=T)
c.rna.seq=sapply(c.rnawigs,wig2rle)

i=seq(from=2,by=2,to=length(c.rna.seq))
c.rna.seq=lapply(i,function(i){list(pos=c.rna.seq[[i]],neg=c.rna.seq[[i-1]])})
names(c.rna.seq) = paste0('tp',gsub(x=c.rnawigs[i],pattern='.*(\\d\\d?)-(\\d\\d?h).*',rep='\\1\\2'))


#Now let's add the Zeitlinger RNAseq as well
zrnabams = list.files('/g/furlong/Harnett/Pioneer_TFs/data/zeitlinger_etal_2012/',pattern='.bam$',full.names=T)
zrnabams =   zrnabams[grepl('06to08h|10to12h|02to04h',zrnabams)]#only 
rna.seq<-mclapply(mc.cores=10,zrnabams,function(bam)bam2coverage(bam,doshift=F,doresize=F,stranded=F,is.ProperPair=NA,is.FirstMateRead=NA,use.size=F,maxsize=500))


#check our RNA

# bam='data/RNAseq_links/samp.bam'
# cov=bam2coverage(bam,doshift=F,doresize=F,stranded=T,is.ProperPair=NA,is.FirstMateRead=NA,use.size=F,maxsize=NA)
# export(cov$pos+cov$neg,'data/RNAseq_links/justreads.bw')

# cov=bam2coverage(bam,doshift=F,doresize=F,stranded=T,is.ProperPair=T,is.FirstMateRead=NA,use.size=F,maxsize=NA)
# export(cov$pos+cov$neg,'data/RNAseq_links/justpairedreads.bw')


# cov=bam2coverage(bam,doshift=F,doresize=F,stranded=T,is.ProperPair=T,is.FirstMateRead=T,use.size=F,maxsize=NA)
# export(cov$pos+cov$neg,'data/RNAseq_links/firstmate.bw')

# cov=bam2coverage(bam,doshift=F,doresize=F,stranded=T,is.ProperPair=T,is.FirstMateRead=F,use.size=F,maxsize=NA)
# export(cov$pos+cov$neg,'data/RNAseq_links/secmate.bw')


# coverageplot(peaks.pos$chr10[wpeaks[1]], peaks.neg$chr10[wpeaks[1]])


# #load the celniker data from tibors bams
# c.rna.bams=list.files('/g/tier2/furlong/pakozdi/alignment/celniker/rna/merged/',pattern='bam$',full.names=T)[3]
# #try different methods of leading it first of all
# tmp1 = mclapply(mc.cores=10,c.rna.bams,function(bam)bam2coverage(bam,doshift=F,doresize=F,stranded=F,is.ProperPair=NA,is.FirstMateRead=NA,use.size=F,maxsize=NA))
# tmp2 = mclapply(mc.cores=10,c.rna.bams,function(bam)bam2coverage(bam,doshift=F,doresize=F,stranded=F,is.ProperPair=NA,is.FirstMateRead=NA,use.size=T,maxsize=600))
# tmp3 = mclapply(mc.cores=10,c.rna.bams,function(bam)bam2coverage(bam,doshift=F,doresize=F,stranded=F,is.ProperPair=NA,is.FirstMateRead=T,use.size=F,maxsize=NA))
# # locatio of twist chr2R:18,931,631-18,937,849
# pdf('coverageplots.celniker.pdf')
# coverageplot(Views(c.rna.seq[[8]][['chr2R']],18931631,18937849),Views(c.rna.seq[[7]][['chr2R']],18931631,18937849))
# coverageplot(Views(tmp1[['chr2R']],18931631,18937849))
# coverageplot(Views(tmp2[['chr2R']],18931631,18937849))
# coverageplot(Views(tmp3[[1]][['chr2R']],18931631,18937849))
# coverageplot(Views(tmp4[[1]][['chr2R']],18931631,18937849))
# dev.off()



