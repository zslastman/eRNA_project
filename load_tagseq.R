###This is a script to generate RLE objects out of the bam files for our CAGE data

# @ author Dermot Harnett, EMBL Heidelberg
# @date 18/10/2013
# @title Load the 3' tag data
########################################
setwd(dir ='/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')

#1 Read in files
tagseq.files = list.files('/g/furlong/mdavis/projects/enrico/data/local_bams',full.names=T,pattern='Tagseq.*.bam$')
linetable    = read.table('data/line_name_table.txt',header=T)

	
tagseq.df<-data.frame(
	accession        = gsub('.*/(.*?)\\..*','\\1',tagseq.files),
	replicate        = gsub(pattern='.*\\dh_((RAL)?\\d+)(_(\\d+))?\\..*',replacement='\\4',x=tagseq.files,perl=T),  
	line             = gsub(pattern='.*\\dh_(RAL)?(\\d+)(_\\d+)?\\..*',replacement='\\2',x=tagseq.files,perl=T),
	timepoint        = gsub(pattern='.*_(\\d+_\\d+)h_(RAL)?\\d+(_\\d)?\\..*',replacement='\\1',x=tagseq.files,perl=T),
	tissue           = 'embryo',
	RAL              = grepl('h_RAL\\d+',x=tagseq.files),#vector describing if the RAL were used (so 0 for bloomington lin number)
	stringsAsFactors = F
)


#convert all our line information to the RAL #
bloom2ral                      = linetable[,2]
names(bloom2ral)               = linetable[,1]
tagseq.df$line[!tagseq.df$RAL] = bloom2ral[tagseq.df$line[ !tagseq.df$RAL]]
#make timepoint info comparable between tagseq and cage tables
tagseq.df$timepoint            =gsub('(\\d+)_(\\d+)','\\1\\2h',tagseq.df$timepoint)
#save
save(tagseq.df,file            ='data/objects/tagseq.df.object.R')  

#3 Process the BAM files, name them with the dataframe, creating and naming our Rle list of stucture name->strand->chr
ts<-mclapply(mc.cores=20,mc.cleanup=T,tagseq.files,function(x)bam2coverage(x,doshift=F,doresize=F,stranded=T))
names(ts)<-tagseq.df$accession#name them
 #create a third strand for each library with the other two summed
for(acc in names(ts)){ts[[acc]]$both<-ts[[acc]]$pos+ts[[acc]]$neg}
#calculate size of  library and put this in the dataframe
tagseq.df$library.size<-sapply(as.character(tagseq.df$accession),function(acc){
   sum(as.numeric(sum(ts[[as.character(acc)]]$both)))
})
tagseq.df$genome.coverage<-(tagseq.df$library.size/sum(seqlengths(si)))
#save
save(tagseq.df,file= 'data/objects/tagseq.df.object.R' )  
save(ts,file='data/objects/ts.object.R')  
#load('data/objects/tagseq.df.object.R')
#load('data/objects/all.tagseq.unprocessed.R')



## 4  do the power law normalization 
fits<-sapply(tagseq.df$accession,function(acc){
  sitecounts=sort(c(unlist(unname(ts[[acc]][['pos']])),unlist(unname(ts[[acc]][['neg']]))))
  sitecounts=sitecounts[sitecounts!=0]
  fit=power.law.fit(as.vector(sitecounts),xmin=200)#not that power.law.fit excludes the lower count numbers
  total.tags=sum(sitecounts)
  o=getOffset(fit$alpha,total.tags)#calculate the offset (second parameter)
  c(fit,offset=o,tagcount=total.tags)#now tack it onto the fit list and return
})
tagseq.df$alpha<-as.numeric(fits['alpha',])
mean.alpha=mean(tagseq.df$alpha)
r.offset <-getOffset(mean.alpha,refsize)
tagseq.df$offset<-as.numeric(fits['offset',])
save(tagseq.df,file='data/objects/tagseq.df.full.object.R')
#4a now finally normalize all of our cage libraries to a common power law reference
ts.pl<-mapply(SIMPLIFY=F,names(ts),tagseq.df$alpha,tagseq.df$offset,FUN=function(acc,alpha,offset){
   lapply(ts[[acc]],function(srle){ 
    pl.norm(srle,x.alpha=alpha,x.offset=offset[[1]])
  })
})
save(ts.pl,file='data/objects/ts.pl.object.R')


#5 create a summed alltags object for each timepoint
 message('summing tagseq tags')
 alltagseq=sapply(simplify=F,split(ts,as.character(tagseq.df$timepoint)),function(alltags){
 	alltags=alltags
 	list(
 			pos=Reduce('+',sapply(alltags,'[[','pos')),
 			neg=Reduce('+',sapply(alltags,'[[','neg')),
 			both=Reduce('+',sapply(alltags,'[[','both'))
 		)
 })
#export bigwigs for these
 for(set in names(alltagseq)){
 	export(alltagseq[[set]]$pos,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/alltagseq.',set,'.pl.pos.bw'))
 	export(alltagseq[[set]]$neg,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/alltagseq.',set,'.pl.neg.bw'))
 }
 save(alltagseq,file ='data/objects/alltagseq.object.R')

#6 load up the peaks and creat files for these as well.
tagseq.peaks=list(
	tp24=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.2h.default.gff.features.gff',header=F,comment.char='#'),
	tp68=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.6h.default.gff.features.gff',header=F,comment.char='#'),
	tp1012=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.10h.default.gff.features.gff',header=F,comment.char='#')
)
tagseq.peaks=sapply(tagseq.peaks,function(f){
	GRanges(paste0('chr',f$V1),IRanges(f$V4,f$V5),strand=f$V7,seqinfo=si)
})
sapply(names(tagseq.peaks),function(npks){
	export(tagseq.peaks[[npks]],paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/tagseqpeaks_',npks,'.bed'))
})
save(tagseq.peaks,file='data/objects/tagseq.peaks.object.R')

