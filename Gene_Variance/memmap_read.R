setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
rootfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/'
outfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/mapfilterplnorm/'
dir.create(outfolder,showWarnings=F)

#load two cage datasets
load(paste0(rootfolder,'data/objects/cg.pl.object.R'))

library(gplots)
library('mmap', lib.loc='/g/furlong/Degner/R/')
library(gplots)

#append chromosome to our mapability file
lowmaplim=1#minimum number of reads before we call a region as mapped
lineallowance=0#number of lines that are allowed to be unmapped
#collect all the mappability info
#collect list of mappability files
filelist <- scan(pipe("du /g/furlong/project/24_TSSCAGE/analysis/MAPPABILITY/MAP_BIN/DNAmapRAL_* | awk '$1 > 1000000 {print $2}'"), what='char')

#create template SRL with all zeroes
mapsum = cg.pl[[1]][[1]]
mapsum[mapsum>0] <- 0
goodcount=0

for(i in (1:length(filelist))){
	# mmfile=filelist[[1]]
	mlist=sapply(chrs.keep,function(chr){#for each chr
		mmchrfile=paste0(mmfile, paste0('/',chr,'.memmap'))#define the file for each chr
		test <- mmap(mmchrfile, mode=int64())#get the accessor
		Rle(test[])#access it and create an Rle
	})
	mlist = as( mlist , 'SimpleRleList' )
	if(any(sum(mlist)==0)){next}
	goodcount=goodcount+1
	mlist=mlist>=lowmaplim
	mapsum = mapsum + mlist
	cat('.')
}

#now create a logical vector describing which locations have non-mapping lines
allmap = mapsum => goodcount-lineallowance

#now save our allmap object
save(allmap,file='data/objects/allmap.object.R')



load('data/objects/cg.pl.object.R')
load('data/objects/allmap.object.R')
# #now make a version with only sites that are totally mappable
# cg.pl.map=rapply(cg.pl,how='replace',function(SRL){
# 	SRL[!allmap]<-0
# # 	SRL
# # })

cg.pl.map=cg.pl

# load('data/objects/cg.pl.object.R')
# load('data/objects/allmap.object.R')
# #now make a version with only sites tha

nallmap=!allmap
cg.pl.map=mclapply(mc.cores=5,names(cg.pl),function(acc){
	sapply(c('pos','neg'),function(s){
		cg.pl[[acc]][[s]][nallmap]<-0
		cat('.')
		cg.pl[[acc]][[s]]
	})
})
names(cg.pl.map) = names(cg.pl)

save(cg.pl.map,file='data/objects/cg.pl.map.object.R')









#all this was assuming they fit in memory, which they probably don't even as Rles
#for each file create a SimpleRleList of coverage
uplim=min(length(filelist),20)
acc.map.list=mclapply(mc.cores=20,filelist[1:uplim],function(mmfile){
	
})
filelist = filelist[ (uplim+1):length() ]
#eliminate the ones that are all zeros on at least one chr
all.chrs.cov = function(SRL){}
acc.map.list = Filter(all.chrs.cov,acc.map.list)
#convert all of them to logical Rles describing if they are above our limit
acc.map.list.l = mclapply(mc.cores=20,acc.map.list,function(SRL){
})
#and reduce the whole thing down to a single logical SRL,counting the lines that map
acc.map.list.r = Reduce('+',acc.map.list.l)
#
linenumlim=length(acc.map.list.l)




 