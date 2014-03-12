setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
library('mmap', lib.loc='/g/furlong/Degner/R/')
rootfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/'
outfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/mapfilterplnorm/'
dir.create(outfolder,showWarnings=F)

klowmaplim=1#minimum number of reads before we call a region as mapped
klineallowance=1#number of lines that are allowed to be unmapped
#collect all the mappability info
#collect list of mappability files
filelist <- scan(pipe("du /g/furlong/project/24_TSSCAGE/analysis/MAPPABILITY/MAP_BIN/DNAmapRAL_* | awk '$1 > 1000000 {print $2}'"), what='char')
#create template SRL with all zeroes
mapsum <- coverage(GRanges(seqinfo=si))
linecount=0#counter
for(i in (1:length(filelist))){
	mmfile=filelist[[i]]
	mlist=sapply(chrs.keep,function(chr){#for each chr
		mmchrfile=paste0(mmfile, paste0('/',chr,'.memmap'))#define the file for each chr
		test <- mmap(mmchrfile, mode=int64())#get the accessor
		Rle(test[])#access it and create an Rle
	})
	mlist = as( mlist , 'SimpleRleList' )#
	if(any(sum(mlist)==0)){next}#skip empty or nearly files
	linecount=linecount+1#incrememnt counter
	mlist=mlist>=klowmaplim#get logical Rle List
	mapsum = mapsum + mlist#add this to our mapsum object
	cat('.')
}
#now create a logical vector describing which locations have non-mapping lines
allmap = mapsum >= linecount-klineallowance

#check it's worked
test_that("the allmap object is correct",{
	expect_is(allmap,'SimpleRleList')
  expect_equal(names(allmap),chrs.keep)
	expect_true(all(mean(allmap)>0))
	expect_true(mean(allmap)[['chr2L']]>mean(allmap)[['chr4']])
})

#now save our allmap object
save(allmap,file='data/objects/allmap.object.R')
save(mapsum,file='data/objects/mapsum.object.R')
save(allmap,file='/g/furlong/project/24_TSSCAGE/data/allmap.object.R')
save(mapsum,file='/g/furlong/project/24_TSSCAGE/data/mapsum.object.R')
export(allmap*1,'/g/furlong/Harnett/TSS_CAGE_myfolder/data/mappability.bw')

