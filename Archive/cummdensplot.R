#quick s
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/generate_Rle.R')
#install.packages('VGAM')

library(ggplot2)


cg=cage.tag.rles
accs=names(cage.tag.rles)

strand='pos'
acc=names(cage.tag.rles)[1]


bin.size=500

#create genomic bins
bin_ranges <- IRangesList(#initially just create disjoint regions 
  sapply(seqlevels(si)[!grepl('gl|hap', seqlevels(si))], function(tmp_chr) {
    IRanges(start = seq(1, seqlengths(si)[[tmp_chr]] - bin.size, bin.size), width = bin.size)
  })
)
bin_ranges<-as(bin_ranges,'GRanges')
seqinfo(bin_ranges)<-si


srle2vector<-function(srlelist){as.vector(unlist(unname(srlelist)))}
log10nz<-function(x){x=log10(x);x[x==-Inf]<-0;stopifnot(all(x>=0));x}

#cumsum.df<-do.call(rbind,lapply(accs[seq(1,length(accs),by=10)],function(acc){

#for each accesion produce the summed tags in each bin
cumsum.df<-do.call(rbind,lapply(accs,function(acc){
	#get the summed reads in each genomic bin for this library
	# sitecounts<-unlist(viewSums(Views.gr(cage.tag.rles[[acc]][[3]],bin_ranges)))
	sitecounts<-c(srle2vector(cg[[acc]][[1]]),srle2vector(cg[[acc]][[2]]))	
	freq<-cumsum(table(sitecounts))
	freq<-max(freq)-freq
	level<-names(freq)
	data.frame(freq,level=as.numeric(level),acc)

	
}))
qplot(data=cumsum.df,color=acc,y=freq,x=as.numeric(level),log='xy',geom='line')


p.df<-sapply(c(),function(p){
	s=rpareto(100000,location=p,shape=1.25)


})

qplot(data=p.df,y=freq,x=as.numeric(level),log='xy',geom='line')

#let's look at the reproducibility scatter plot
load('data/objects/unprocessed.cage.tag.object.R')
u.cg<-cage.tag.rles
load(file.cage.tag.rles)



  reseqs<-names(u.cg)[grepl(names(u.cg),pattern='reseq')]
  
  sapply(reseqs,function(reseq){
    orig<-gsub(reseq,pattern='_reseq',replacement='')
    eitherpos<-u.cg[[reseq]]$pos+u.cg[[orig]]$pos > 0
    eitherneg<-u.cg[[reseq]]$neg+u.cg[[orig]]$neg > 0

    reseqvals<-c(srle2vector(u.cg[[reseq]][[1]][eitherpos]),srle2vector(u.cg[[reseq]][[2]][eitherneg]))
    origvals<-c(srle2vector(u.cg[[orig]][[1]][eitherpos]),srle2vector(u.cg[[orig]][[2]][eitherneg]))
    
    reseqvals<-log10nz(reseqvals)
    origvals<-log10nz(origvals)

    qplot(reseqvals,origvals,geom='point',log='',alpha=0.2)

})



names(u.cg)


