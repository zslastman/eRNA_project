#We want to address some quick questions re mapability, that might bear on this weird peak of correlation I'm
#Seeing at the 5' vs 3' correlation

setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R')

######1 load up rles and metadata on them
#annotation and metadata
load(file.transcripts)
load('data/objects/accession.df.full.object.R')
load('data/objects/tagseq.df.object.R')
load('/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/cg.object.R')
#cage and tagseq data
cg=cg.pl.map
load('data/objects/all.tagseq.unprocessed.R')
ts=ts
rm(ts)
iterations=1000#number of iterations when doing empirical pvalues
#names of timepoints
tps=c('tp24h','tp68h','tp1012h')
load('data/objects/allmap.object.R')



#1 autocorrelation of mappability - how far does it extend?
acfdists=10^(seq(1,5,length.out=200))#distances to calculate acf at
acf.list=sapply(simplify=F,names(allmap),function(chr){acf(lag=acfdists,as.vector(as.numeric(allmap[[chr]])),plot=F)})
dir.create('analysis/mappability_5prime_3prime')
pdf('analysis/mappability_5prime_3prime/autocor.pdf',h=14,w=14)
for(chr in names(acf.list)){
	plot(acf.list[[chr]])
}
dev.off()

#2 correlation in mappability of 5' end and 3' end
	#does it exist?
	#is it related to the correlation 

######3 Load up the peaks for the 3' data as well.
tagseq.peaks=list(
	tp24=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.2h.default.gff.features.gff',header=F,comment.char='#'),
	tp68=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.6h.default.gff.features.gff',header=F,comment.char='#'),
	tp1012=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.10h.default.gff.features.gff',header=F,comment.char='#')
)
tagseq.peaks=sapply(tagseq.peaks,function(f){
	GRanges(paste0('chr',f$V1),IRanges(f$V4,f$V5),strand=f$V7,seqinfo=si)
})
#create 500bp windows around the end of trnascripts for finding our tagseq peaks
filt.transcripts.gr = transcripts.gr[countOverlaps(transcripts.gr,transcripts.gr)==1]
ends.gr=resize(filt.transcripts.gr,width=500,fix='end')
ends.gr=shift(ends.gr,250)
#just 1 tp for now
tp=tps[2]
#select tss with only 1 peak at each tp
fov = findOverlaps(ends.gr,tagseq.peaks[[tp]])
fov = fov[!duplicated(fov@queryHits) & !duplicated(fov@subjectHits),]
#get the tss
utss=resize(filt.transcripts.gr[ fov@queryHits ],width=1,fix='start')
utss=sort(resize(utss,width=500,fix='center'))
#and en of transcripts from the peak data
uends = sort(tagseq.peaks[[tp]][ fov@subjectHits ])
stopifnot(length(uends)==length(utss))

#Now get mappability data as means
utss$map=unlist(viewMeans(GRViews(allmap,utss))
uends$map=unlist(viewMeans(GRViews(allmap,uends))


