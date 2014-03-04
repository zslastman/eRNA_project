setwd('~/Harnett/TSS_CAGE_myfolder')
source('src/tss_cage_functions.R')
genomesize <- sum(seqlengths(Dmelanogaster))#genome length

kchromatinFolder = '/g/tier2/furlong/project/11_Histone_ChIP_seq/BiTS-ChIP/alignment/'
fraglength=190
chrom_bam_files =  list(
    PolII_6.8_R1= paste0( kchromatinFolder, '6-8h/bam/B-PolII-Rpb3-6-8h.bam'),
    PolII_6.8_R2= paste0( kchromatinFolder, '6-8h/bam/B51_PolII-Rpb3_6-8h.bam'),
    H3K27ac_6.8_R1= paste0( kchromatinFolder, '6-8h/bam/B23-K27ac.bam'),
    H3K27ac_6.8_R2= paste0( kchromatinFolder, '6-8h/bam/B-A-K27ac.bam'),
    H3K4me3_6.8_R1= paste0( kchromatinFolder, '6-8h/bam/B24-K4me3.bam'),
    H3K4me3_6.8_R2= paste0( kchromatinFolder, '6-8h/bam/B53_K4me3_6-8h.bam'),
    H3K4me1_6.8_R1= paste0( kchromatinFolder, '6-8h/bam/B-E-K4me1.bam'),
    H3K4me1_6.8_R2= paste0( kchromatinFolder, '6-8h/bam/B-F-K4me1.bam'),
    H3K79me3_6.8_R1= paste0( kchromatinFolder, '6-8h/bam/B38_K79me3.bam'),
    H3K79me3_6.8_R2= paste0( kchromatinFolder, '6-8h/bam/B52_K79me3_6-8h.bam'),
    H3K36me3_6.8_R1= paste0( kchromatinFolder, '6-8h/bam/B02-K36me3.bam'),
    H3K36me3_6.8_R2= paste0( kchromatinFolder, '6-8h/bam/B-B-K36me3.bam'),
    H3K27me3_6.8_R1= paste0( kchromatinFolder, '6-8h/bam/B-K27me3-6-8h.bam'),
    H3K27me3_6.8_R2= paste0( kchromatinFolder, '6-8h/bam/B54_K27me3_6-8h.bam'))

expect_true(all(sapply(chrom_bam_files,file.exists)))
  #  H3K36me3_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-K27me3-6-8h.bam'

  #  H3K36me3_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-K27me3-6-8h.bam'

input_bam_files =  list(
    Input_6.8_R1='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B03-input.bam',
    Input_6.8_R2='/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/bam/B-D-input.bam'
)

if(file.exists(file.chrom.rles.rpgc.sub.merge)){
  load(file.chrom.rles.rpgc.sub.merge)
  message('Chromatin rle files already present, loaded')
}else{
  message(paste0('\nloading the following files  -  \n',chrom_bam_files, '\n', input_bam_files))
  #read chromatin files
  chrom.rles<-mclapply(mc.cores=10,chrom_bam_files,function(x)bam2coverage(x,doshift=T,fragment.length=fraglength,use.size=F))
  input.rles<-mclapply(mc.cores=10,input_bam_files,function(x)bam2coverage(x,doshift=F,doresize=F,use.size=F))
  sapply(chrom.rles,expect_SRL)
  sapply(input.rles,expect_SRL)

  
  names(chrom.rles)<-names(chrom_bam_files)#name them
  names(input.rles)<-names(input_bam_files)#name them
  
  save(chrom.rles,file=file.chrom.rles)
  save(input.rles,file=file.input.rles)
  
  # load(file.input.rles)
  # load(file.chrom.rles)
  
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
  #get rid of R numbers since they're merged
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
  message(filen)
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
  K27ac_0.4h= 'data/modencode/coverage_wigs/0-4hr/H3K27ac-Developmental-Stage=Embryos-0-4-hr#Strain=Y-cn-bw-sp-ChIP-seq-Rep-1-ChIP-Dmel_r5.32-modENCODE_834.wig',
  K27ac_4.8h='data/modencode/coverage_wigs/4-8hr/H3K27ac-Developmental-Stage=Embryos-4-8-hr#Strain=Y-cn-bw-sp-ChIP-seq-Rep-1-ChIP-Dmel_r5.32-modENCODE_835.wig',
  K27ac_8.12h='data/modencode/coverage_wigs/8-12hr/H3K27ac-Developmental-Stage=Embryos-8-12-hr#Strain=Y-cn-bw-sp-ChIP-seq-Rep-1-ChIP-Dmel_r5.32-modENCODE_836.wig',
  PolII_4.8h="/g/furlong/Harnett/TSS_CAGE_myfolder/data/modencode/coverage_wigs/4-8hr/PolII-E4-8__density.wig"
)

expect_true( all( sapply(chrom.rles.modencode,file.exists)))

chrom.rles.modencode<-sapply(chrom.rles.modencode,function(wigfile){
  if( grepl('.gz$',wigfile) ){
    system(paste0('gunzip ', wigfile))
    wigfile = gsub(wigfile,pat='(.*)\\.gz$',rep='\\1')
  }
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



################################################################################################
# Also load dnase peaks ---------------------------------------------------
#first just load the peaks for stages 10,11
if(file.exists(file.dnase.peaks)){
  load(file.dnase.peaks)
}else{
  dnase.files<-list.files(full.names=T,'/g/furlong/Harnett/TSS_CAGE_myfolder/data/dnase_thomas_etal/',pattern='.bed$')
  dnase.peaks<-sapply(dnase.files,function(f){
    tmp<-import(f,asRangedData=F)  
    tmp$stage = gsub('.*bdtnpDnaseAccS(\\d\\d?)\\.bed','stg\\1',f)
    tmp
  })
  stopifnot(all(grepl('bdtnpDnaseAccS\\d\\d?.bed$',dnase.files)))
  names(dnase.peaks) =   gsub('.*bdtnpDnaseAccS(\\d\\d?)\\.bed','stg\\1',names(dnase.peaks))
 # dnase.peaks = GRangesList(dnase.peaks)
  dnase.peaks<-reduce(do.call('c',unname(dnase.peaks)))
  # dnase.peaks$name<-paste0('dnase_10.11',1:length(dnase.peaks))
  dnase.peaks<-keepSeqlevels(dnase.peaks,chrs.keep[chrs.keep %in% seqlevels(dnase.peaks)])
  seqinfo(dnase.peaks) <- si
  save(dnase.peaks,file=file.dnase.peaks)
}
