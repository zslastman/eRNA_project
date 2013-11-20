setwd('~/Harnett/TSS_CAGE_myfolder')
source('src/tss_cage_functions.R')
genomesize <- sum(seqlengths(Dmelanogaster))#genome length
file.chrom.rles<-'data/objects/chrom.rles.object.R'
file.chrom_noctl.rles<- 'data/objects/chrom_noctl.rles.object.R'
file.input.rles<-'data/objects/input.rles.object.R'
file.chrom.rles.rpgc<-'data/objects/chrom.rles.rpgc.R'
file.input.rles.rpgc<-'data/objects/input.rles.rpgc.R'
file.chrom.rles.rpgc.sub.merge<-'data/objects/chrom.rles.rpgc.sub.merge.object.R'

fraglength=190



#list.files('/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/',pattern='Pol')

#list.files(recursive=T,'data/chromatin/bam')

chrom_bam_files =  list(
    PolII_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-PolII-Rpb3-6-8h.bam',
    PolII_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B51_PolII-Rpb3_6-8h.bam',
    H3K27ac_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B23-K27ac.bam',
    H3K27ac_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-A-K27ac.bam',
    H3K4me3_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B24-K4me3.bam',
    H3K4me3_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B53_K4me3_6-8h.bam',
    H3K4me1_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-E-K4me1.bam',
    H3K4me1_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-F-K4me1.bam',
    H3K79me3_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B38_K79me3.bam',
    H3K79me3_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B52_K79me3_6-8h.bam',
    H3K36me3_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B02-K36me3.bam',
    H3K36me3_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-B-K36me3.bam',
    H3K27me3_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-K27me3-6-8h.bam',
    H3K27me3_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B54_K27me3_6-8h.bam'
  #  H3K36me3_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-K27me3-6-8h.bam'
    
    )
input_bam_files =  list(
    Input_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B03-input.bam',
    Input_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-D-input.bam'
)



if(file.exists(file.chrom.rles.rpgc.sub.merge)){
  load(file.chrom.rles.rpgc.sub.merge)
  message('Chromatin rle files already present, loaded')
}else{
  
  message(paste0('\nloading the following files  -  \n',chrom_bam_files, '\n', input_bam_files))
  #let's make this nice and neat with rapply later
  
  chrom.rles<-mclapply(mc.cores=10,chrom_bam_files,function(x)bam2coverage(x,doshift=T,fragment.length=180,use.size=F))
  input.rles<-mclapply(mc.cores=10,input_bam_files,function(x)bam2coverage(x,doshift=F,doresize=F,use.size=F))
  
  stopifnot(class(chrom.rles[[1]])=='SimpleRleList')
  
  names(chrom.rles)<-names(chrom_bam_files)#name them
  names(input.rles)<-names(input_bam_files)#name them
  
  save(chrom.rles,file=file.chrom.rles)
  save(input.rles,file=file.input.rles)
  
  load(file.input.rles)
  load(file.chrom.rles)
  
  message('rpgc normalization')
  chrom.gcov <- lapply(chrom.rles,function(j) sum(sapply(j, function(i) sum(as.numeric(i)))) / genome.size)
  input.gcov <- lapply(input.rles,function(j) sum(sapply(j, function(i) sum(as.numeric(i)))) / genome.size)
  
  chrom.rles.rpgc<-mapply(chrom.rles,chrom.gcov,FUN=function(x,y){x/y})
  input.rles.rpgc<-mapply(input.rles,input.gcov,FUN=function(x,y){x/y})
  
  message('saving rpgc coverage')
  
#    load(file.chrom.rles.rpgc)
#    load(file.input.rles.rpgc)
  save(chrom.rles.rpgc,file=file.chrom.rles.rpgc)
  save(input.rles.rpgc,file=file.input.rles.rpgc)
  
#      load(file.chrom_noctl.rles.)
  
  message('Input subtraction')
  input.rlelist.rpgc.merge<-(input.rles.rpgc[[1]]+input.rles.rpgc[[2]])/2#combine the input replicates
  chrom.rles.rpgc.sub<-sapply(chrom.rles.rpgc,function(srlelist){srlelist - input.rlelist.rpgc.merge})#subtract
  
  #we'll just do this by indexes, assumes pairs of replicates in order
  message('combining replicates')
  chrom.rles.rpgc.sub.merge<-lapply(seq(from=1,to=length(chrom.rles),by=2),function(i){chrom.rles.rpgc.sub[[i]]+chrom.rles.rpgc.sub[[i+1]]/2})
  #we may want the unsubtracted ones too
  chrom.rles.rpgc.merge<-lapply(seq(from=1,to=length(chrom.rles),by=2),function(i){chrom.rles.rpgc[[i]]+chrom.rles.rpgc[[i+1]]/2})
  
  #name
  names(chrom.rles.rpgc.sub.merge)<-(names(chrom_bam_files))[seq(from=1,to=length(chrom.rles),by=2)]
  #get rid of R numbers
  names(chrom.rles.rpgc.sub.merge)<-lapply(names(chrom.rles.rpgc.sub.merge),function(name){gsub(x=name,pattern='_R\\d',replacement='')})
  names(chrom.rles.rpgc.merge)<-lapply(names(chrom.rles.rpgc.sub.merge),function(name){gsub(x=name,pattern='_R\\d',replacement='')})
  
  save(chrom.rles.rpgc.sub.merge,file=file.chrom.rles.rpgc.sub.merge)
  
#   
#   #lets just look at these in igv
#   #bothnoctl<-(chrom_nctl.rlelists[[1]]+chrom_nctl.rlelists[[2]])/2
#   export(chrom.rles$PolII_6.8_R1,'data/unprocessedpolII68.bigWig')
#   export(chrom.rles.rpgc.sub.merge,'data/processedpolII68.Bigwig')
#   
#   
#   

}


# Load the MACs Peaks for our Chromatin data ------------------------------
chrompeakfiles<-list.files('data/chromatin/',pattern='.bed|.txt',full.names=T)

chrompeaks<-sapply(chrompeakfiles,function(filen){
  cat(filen)
  suppressWarnings({
      #import the bed files
    if(grepl('.bed',filen)){
      x=import(filen,asRangedData=F)
      if(!grepl('chr',seqnames(x)[1])){
        seqlevels(x)=paste0('chr',seqlevels(x))
      }
      x=keepSeqlevels(x,chrs.keep)
      seqlevels(x)=seqlevels(si)
      seqinfo(x)=si
      x
    }else{
    x=read.delim(header=T,filen,comment.char='#')  
    x[,1]=paste0('chr',x[,1])
    x=x[ x$chr %in% chrs.keep, ]
    x.gr<-GRanges(x[,1],IRanges(x[,2],x[,3]),seqinfo=si)
    mcols(x.gr)<-x[,-c(1,2,3)]
    start(x.gr)<-pmax(start(x.gr),1)
    end(x.gr)<-pmin(end(x.gr),seqlengths(x.gr)[as.vector(seqnames(x.gr))])
    x.gr
    }
  })
})

#format the names of the peakfiles
names(chrompeaks)<-gsub('.*//(\\w+?)_.*(\\d\\-\\dh)_.*','\\1_\\2',names(chrompeaks))




# Get the peaks from the modencode data as well ---------------------------
#For now I'll use the peaks on our server which Nico made using genomic input DNA
chrompeakfiles<-list.files('data/modencode/peaks/',pattern='.txt|.bed|.tab',full.names=T)
chrompeaks.modencode<-sapply(chrompeakfiles,function(x){
  suppressWarnings({
      #import the bed files
    if(grepl('.bed',x)){
      x=import(x,asRangedData=F)
      if(!grepl('chr',seqnames(x)[1])){
        seqlevels(x)=paste0('chr',seqlevels(x))
      }
      x=keepSeqlevels(x,chrs.keep)
      seqlevels(x)=seqlevels(si)
      seqinfo(x)=si
      x
    }else{
    x=read.delim(header=T,x)  
    x[,1]=paste0('chr',x[,1])
    x=x[ x$chr %in% chrs.keep, ]
    x.gr<-GRanges(x[,1],IRanges(x[,2],x[,3]),seqinfo=si)
    mcols(x.gr)<-x[,-c(1,2,3)]
    start(x.gr)<-pmax(start(x.gr),1)
    end(x.gr)<-pmin(end(x.gr),seqlengths(x.gr)[as.vector(seqnames(x.gr))])
    x.gr
    }
  })
})
#format the names of the peakfiles
names(chrompeaks.modencode)<-  gsub('.*//(\\w+?)_.*(\\d\\-\\dh)_.*','\\1_\\2',names(chrompeaks.modencode))


# Modencode continous signal ----------------------------------------------

#modencode data
chrom.rles.modencode<-list(  
  K27ac='data/modencode/coverage_wigs/4-8hr/H3K27ac-Developmental-Stage=Embryos-4-8-hr#Strain=Y-cn-bw-sp-ChIP-seq-Rep-1-ChIP-Dmel_r5.32-modENCODE_835.wig',
  PolII='data/modencode/coverage_wigs/4-8hr/PolII-E4-8__density.wig'
               )

chrom.rles.modencode<-sapply(chrom.rles.modencode,function(wigfile){
  mark.modencode<-import(format='wig',wigfile)
  mark.modencode<-as(mark.modencode[chrs.keep],'GRanges')
  seqlevels(mark.modencode)<-seqlevels(si)
  seqinfo(mark.modencode)<-si
  coverage(mark.modencode,weight='score')
})


################################################################################
#Dnase data
###############################################################################
message("loading dnase from bam files")
if(!file.exists(file.dnase.rles)){
  dnase_filelist<-list.files(full.names=T,'data/dnase_thomas_etal/', pattern='R\\d.combined.bam$')
  dnase_filelist<-sort(dnase_filelist)
  #actually we should just pull this out
  dnase.rles<-mclapply(mc.cores=10,dnase_filelist,function(x)bam2coverage(x,doshift=F,doresize=F))
  names(dnase.rles)<-gsub('.*(STG\\d\\d?)_.*','\\1',dnase_filelist)
  #add together the dnase rle lists with the same name(stage)
  dnase.rles<-sapply(names(dnase.rles)[!duplicated(names(dnase.rles))],function(stg){do.call('+',unname(dnase.rles[names(dnase.rles)==stg]))})
  message('saving this very large chromdata file')
  save(dnase.rles,file=file.dnase.rles)
}else{load(file.dnase.rles)}
