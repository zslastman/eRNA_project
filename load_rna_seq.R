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
rna.seq = list(
  meso = sum.list.of.rle.lists(rna.seq[is.meso]) ,
  embryo = sum.list.of.rle.lists(rna.seq[!is.meso]) 
  )
#now produce the coverage wigs
export(rna.seq[[embryo]]$pos,'data/RNAseq_links/secondmate/all_rnaseq_emb_6.8.pos.bw')
export(rna.seq[[embryo]]$neg,'data/RNAseq_links/secondmate/all_rnaseq_emb_6.8.neg.bw')
export(rna.seq[[meso]]$pos,'data/RNAseq_links/secondmate/all_rnaseq_meso_6.8.pos.bw')
export(rna.seq[[meso]]$neg,'data/RNAseq_links/secondmate/all_rnaseq_meso_6.8.neg.bw')
save(rna.seq,file='data/objects/rna.seq.object.R')
#now produce the coverage wigs
export(embsum$pos,'data/RNAseq_links/secondmate/all_rnaseq_emb_6.8.pos.bw')
export(embsum$neg,'data/RNAseq_links/secondmate/all_rnaseq_emb_6.8.neg.bw')
export(mesosum$pos,'data/RNAseq_links/secondmate/all_rnaseq_meso_6.8.pos.bw')
export(mesosum$neg,'data/RNAseq_links/secondmate/all_rnaseq_meso_6.8.neg.bw')



################################################################################
#2 Now get the coverage wigs from teh Celniker total RNAseq experiments
c.rnawigs=list.files('data/graveley_totalrnaseq/',full.names=T,pattern='wig',recursive=T)
# coverageplot(peaks.pos$chr10[wpeaks[1]], peaks.neg$chr10[wpeaks[1]])
c.rna.seq.wig=sapply(c.rnawigs,wig2rle)
#now sum the strands
i=seq(from=2,by=2,to=length(c.rna.seq.wig))
c.rna.seq.wig=lapply(i,function(i){list(pos=c.rna.seq.wig[[i]],neg=c.rna.seq.wig[[i-1]])})
#and fix the names
names(c.rna.seq.wig) = paste0('tp',gsub(x=c.rnawigs[i],pattern='.*(\\d\\d?)-(\\d\\d?h).*',rep='\\1\\2'))
######################################
#load the celniker data from tibors bams
c.rna.bams=list.files('/g/tier2/furlong/pakozdi/alignment/celniker/rna/merged/',pattern='bam$',full.names=T)
c.rna.seq = mclapply(mc.cores=10,c.rna.bams,function(bam)bam2coverage(bam,doshift=F,doresize=F,stranded=F,is.ProperPair=NA,is.FirstMateRead=T,use.size=T,maxsize=600))
names(c.rna.seq) = paste0('tp',gsub(x=c.rna.bams,pattern='.*(\\d\\dh).*',rep='\\1'))
c.rna.seq = sapply(c.rna.seq,function(x){list(both=x)}


#export bigwigs for these
for(set in names(c.rna.seq)){
    export( c.rna.seq[[set]],paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/celniker.rnaseq.',set,'.bw'))
}


################################################################################
#3 Now let's add the Zeitlinger RNAseq as well
zrnabams = list.files('/g/furlong/Harnett/Pioneer_TFs/data/zeitlinger_etal_2012/',pattern='.bam$',full.names=T)
zrnabams =   zrnabams[grepl('06to08h|10to12h|02to04h',zrnabams)]#only 
z.rna.seq<-mclapply(mc.cores=10,zrnabams,function(bam)bam2coverage(bam,doshift=F,doresize=F,stranded=F,is.ProperPair=NA,is.FirstMateRead=NA,use.size=F,maxsize=500))
#save it all
save(rna.seq,c.rna.seq.wig,c.rna.seq,file='data/objects/rna.seq.object.R')
z.rna.rles


