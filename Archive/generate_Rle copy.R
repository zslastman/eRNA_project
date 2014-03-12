###This is a script to generate RLE objects out of the bam files for our CAGE data

# @ author Dermot Harnett, EMBL Heidelberg
# @date 16/5/2013
# @title Export coverage into wig format
########################################
setwd(dir='/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
library(VGAM,lib.loc='~/Harnett/R')
library(igraph,lib.loc='~/Harnett/R')
library(ggplot2)
#variables
refsize=100000000


################################################################################
#1 read in the bam files as rles
cage_bam_files.68 = list.files('/g/furlong/Harnett/24_TSSCAGE/data/solexa/bam/',full.names=T,pattern='_68h.*trimmed.*quality.bwa.sort.bam$')
cage_bam_files.24 = list.files('/g/furlong/Harnett/24_TSSCAGE/data/solexa/bam/',full.names=T,pattern='24h.*trimmed.*quality.bwa.sort.bam$')
cage_bam_files.1012 = list.files('/g/furlong/Harnett/24_TSSCAGE/data/solexa/bam/',full.names=T,pattern='1012h.*trimmed.*quality.bwa.sort.bam$')
cage_bam_files.meso= list.files('/g/furlong/Harnett/24_TSSCAGE/data/solexa/bam/',full.names=T,pattern='splitmeso.*sort.*.bam$')



#2 perpare metadata dataframe from the filenames
accession.df<-data.frame(
  accession=gsub(pattern='.*(split[meso]?\\w\\w\\w?_\\w+_\\d\\d?\\d\\d?h.*?).trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  sample=gsub(pattern='.*split(\\d\\d?\\w)_\\w+_\\d\\d?\\d\\d?h.*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  run=gsub(pattern='.*split(\\d\\d?)\\w_\\w+_\\d\\d?\\d\\d?h.*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  lane=gsub(pattern='.*split\\d\\d?(\\w)_\\w+_\\d\\d\\d?\\d?h.*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  reseq= grepl(pattern='reseq',x=cage_bam_files,perl=T),
  line=gsub(pattern='.*split\\w\\d?\\w_(\\d+).*?\\d\\d\\d?\\d?h.*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  timepoint=gsub(pattern='.*split\\w\\d?\\w_\\w+_(\\d\\d\\d?\\d?h).*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  tissue='embryo',
  stringsAsFactors=F
  )  
#add data for mesodermal cage libraries
accession.df=rbind(accession.df,data.frame(
  accession=gsub(pattern='.*splitmeso_(\\w\\d\\d?_r\\d).*trimmed.*',replacement='\\1',x=cage_bam_files.meso,perl=T),
  sample = NA,
  run=NA,
  lane=NA,
  reseq=F,
  line='BiTS_line',
  timepoint=gsub(pattern='.*meso(\\d\\dh).*',replacement='\\1',x=cage_bam_files.meso,perl=T),
  tissue='meso',
  stringsAsFactors=F
  ))
#convert all our line information to the RAL #
accession.df$RAL=!grepl('25\\d\\d\\d',x=accession.df$line)
linetable = read.table('data/line_name_table.txt',header=T)
bloom2ral = linetable[,2]
names(bloom2ral) = linetable[,1]
accession.df$line[ !accession.df$RAL ] = bloom2ral[accession.df$line[ !accession.df$RAL]]
#now insert the correct line numbers into our accessions
inds=!is.na(accession.df$sample)
accession.df$acc[inds] = paste0('split',accession.df$sample[inds],'_',accession.df$line[inds],'_',accession.df$timepoint[inds],sep='')
# #now filter out duplicates
# dup=duplicated(accession.df)
accession.df[dup,]$acc = paste0(accession.df[dup,]$acc,'_reseq')
accession.df[dup,]$reseq = T
# cage_bam_files=cage_bam_files[!dup]



#5 Read in data on replicates
message('Reading in data on replicates')
#tp - line - col - prep - seq
accession.df$collection = 1:nrow(accession.df)
accession.df$prep = 1:nrow(accession.df)
accession.df$seq = 1:nrow(accession.df)
#get info on bad samples and replicates
replicatetable = read.table('data/replicate_pairs.txt',header=T,sep='\t')
badsamples=as.character(unlist(replicatetable[as.logical(replicatetable$IsBad),1:2]))
accession.df$isbadsample=accession.df$sample %in% badsamples
#go through each accession, and if it's a 'reseq', change it's collection and prep to the match
for(i in 1:nrow(accession.df)){
  if(accession.df[i,]$reseq){
      reseq<-accession.df$accession[i]
      orig<-gsub(reseq,pattern='_reseq',replacement='')
      o=which(accession.df$accession==orig)
      stopifnot(length(o)>0 & accession.df$line[i]==accession.df$line[o])
      accession.df$collection[i]=accession.df$collection[o]
      accession.df$prep[i]=accession.df$collection[o]
      # accession.df$
  }
}
#go through each technical replicate and change the collections to be the same
for(i in 1:nrow(replicatetable)){
  if(replicatetable[i,]$ReplicateType=='Technical'){
      r1<-replicatetable$Sample.ID.r1[i]
      r2<-replicatetable$Sample.ID.r1[i]
      i=which(accession.df$sample==r1)
      o=which(accession.df$sample==r2)
      stopifnot(length(o)>0 & accession.df$line[i]==accession.df$line[o])
      accession.df$collection[i]=accession.df$prep[o]
  }
}
save(accession.df,file='data/objects/accession.df.full.object.R')





# Read in the bam files as RLE lists --------------------------------------
cage_bam_files=c(cage_bam_files.68,cage_bam_files.24,cage_bam_files.1012)
message('reading bam files....')
cage.tag.rles<-mclapply(mc.cores=5,mc.cleanup=T,c(cage_bam_files),function(x)bam2coverage(x,doshift=F,doresize=T,stranded=T,fragment.length=1))
save(cage.tag.rles,file='data/objects/all.cage.unprocessed.object.R')
#now for the mesodermal
cage.tag.rles<-c(cage.tag.rles,mclapply(mc.cores=5,mc.cleanup=T,c(cage_bam_files.meso),function(x)bam2coverage(x,doshift=F,doresize=T,stranded=T,fragment.length=1)))
save(cage.tag.rles,file='data/objects/all.cage.unprocessed.object.R')

#and name our Rle object appropriately
names(cage.tag.rles)<-accession.df$accession  
#create a third strand for each library with the other two summed
for(acc in names(cage.tag.rles)){cage.tag.rles[[acc]][['both']]<-cage.tag.rles[[acc]]$pos+cage.tag.rles[[acc]]$neg}
#calculate size of  library
accession.df$library.size<-sapply(accession.df$accession,function(acc){sum(as.numeric(sum(cage.tag.rles[[acc]]$both)))})
accession.df$genome.coverage<-(accession.df$library.size/sum(seqlengths(si)))
#and save
save(cage.tag.rles,file='data/objects/all.cage.unprocessed.object.R')
save(accession.df,file='data/objects/accession.df.full.object.R')
#load('data/objects/all.cage.unprocessed.object.R')
#load('data/objects/accession.df.full.object.R')








#6
#create a summed alltags object for each timepoint
#system.time({
message('summing cage tags')
#split our accs into a list of grouped libraries
splitaccs=with(accession.df,split(accession,paste0(timepoint,tissue)))
#now use this list to sum up the libraries
allcage=mclapply(mc.cores=4,splitaccs,function(accs){
  accs=accs
  cage.tag.rles=cage.tag.rles[accs]
  alltags=list(
      pos=Reduce('+',sapply(cage.tag.rles,'[[','pos')),
      neg=Reduce('+',sapply(cage.tag.rles,'[[','neg'))
    )
  alltags$both=alltags$pos+alltags$neg#and add a both rle for each one.
  alltags
})
#export bigwigs for these
for(set in names(allcage)){
    export(allcage[[set]]$pos,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/allcage.',set,'.pl.pos.bw'))
    export(allcage[[set]]$neg,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/allcage.',set,'.pl.neg.bw'))
}
#Also one which just has the sum of ALL sites

save(allcage,file='data/objects/allcage.object.R')


## 4  do the power law normalization 
fits<-sapply(accession.df$accession,function(acc){
  sitecounts=sort(c(unlist(unname(cage.tag.rles[[acc]][['pos']])),unlist(unname(cage.tag.rles[[acc]][['neg']]))))
  sitecounts=sitecounts[sitecounts!=0]
  fit=power.law.fit(as.vector(sitecounts),xmin=200)#not that power.law.fit excludes the lower count numbers
  total.tags=sum(sitecounts)
  o=getOffset(fit$alpha,total.tags)#calculate the offset (second parameter)
  c(fit,offset=o,tagcount=total.tags)#now tack it onto the fit list and return
})
accession.df$alpha<-fits['alpha',]
accession.df$offset<-fits['offset',]
save(accession.df,file='data/objects/accession.df.full.object.R')
save(cg.68.pl,file='data/objects/cg.68.pl')
#4a now finally normalize all of our cage libraries to a common power law reference
cg.pl<-mapply(SIMPLIFY=F,names(cage.tag.rles),accession.df$alpha,accession.df$offset,FUN=function(acc,alpha,offset){
   lapply(cage.tag.rles[[acc]],function(srle){ 
    pl.norm(srle,x.alpha=alpha,x.offset=offset[[1]])
  })
})
save(cg.pl,file='data/objects/cg.all.pl')
#4b normalize the middle time points ot eachother 
emb.68<-accession.df$timepoint=='68h' & accession.df$tissue=='embryo'
mean.alpha=mean(unlist(accession.df$alpha[emb.68]))
r.offset<-getOffset(mean.alpha,refsize)
cg.68.pl<-mapply(SIMPLIFY=F,accession.df$accession[emb.68],accession.df$alpha[emb.68],accession.df$offset[emb.68],FUN=function(acc,alpha,offset){
 lapply(cage.tag.rles[[acc]],function(srle){ 
  pl.norm(srle,x.alpha=alpha,x.offset=offset[[1]])
})
})






################7
#now create rles with our replciates and reseqs merged
#now create lists, but 
cg.merge=cage.tag.rles
message('CAGE Rle objects calculated, now summing reseqs')
reseqs<-accession.df$acc[grepl(accession.df$acc,pattern='reseq')]
for(reseq in reseqs){
  orig<-gsub(reseq,pattern='_reseq',replacement='')
  cg.merge[[orig]]$pos<-cg.merge[[orig]]$pos+cg.merge[[reseq]]$pos
  cg.merge[[orig]]$neg<-cg.merge[[orig]]$pos+cg.merge[[reseq]]$neg
}
#now get rid of the reseqs  
cg.merge<-cg.merge[!names(cg.merge) %in% reseqs]
#now the replicates
message('... now summing replicates')
replicates<-accession.df$accession[grepl(accession.df$accession,pattern='(+\\d)_?r2')]  
for(replicate in replicates){
  orig<-gsub(reseq,pattern='(+\\d)_?r2',replacement='\\1')
  cg.merge[[orig]]$pos<-cg.merge[[orig]]$pos+cg.merge[[replicate]]$pos
  cg.merge[[orig]]$neg<-cg.merge[[orig]]$pos+cg.merge[[replicate]]$neg
}
#now get rid of the replicates  
cg.merge<-cg.merge[!names(cg.merge) %in% replicates]
save(cg.merge,file='data/objects/cg.merge.object.R')


#7b
###and for the normalized data

message('CAGE Rle objects calculated, now summing reseqs')

reseqs<-names(cg.pl)[grepl(names(cg.pl),pattern='reseq')]

for(reseq in reseqs){
  orig<-gsub(reseq,pattern='_reseq',replacement='')
  cg.pl[[orig]]$pos<-cg.pl[[orig]]$pos+cg.pl[[reseq]]$pos
  cg.pl[[orig]]$neg<-cg.pl[[orig]]$pos+cg.pl[[reseq]]$neg
}

  #now get rid of the reseqs  
cg.pl<-cg.pl[!names(cg.pl) %in% reseqs]

message('CAGE Rle objects calculated, now summing replicates')

replicates<-names(cg.pl)[grepl(names(cg.pl),pattern='(+\\d)_?r2')]  
for(replicate in replicates){
  orig<-gsub(reseq,pattern='(+\\d)_?r2',replacement='\\1')
  cg.pl[[orig]]$pos<-cg.pl[[orig]]$pos+cg.pl[[replicate]]$pos
  cg.pl[[orig]]$neg<-cg.pl[[orig]]$pos+cg.pl[[replicate]]$neg
}

  #now get rid of the replicates  
cg.pl<-cg.pl[!names(cg.pl) %in% replicates]

save(cg.pl,file='data/objects/cg.pl.merge.object.R')





save(allcage,file='data/objects/allcage.object.R')
save(accession.df,file=file.accession.df)






