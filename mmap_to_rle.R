library(GRanges)
library('mmap', lib.loc='/g/furlong/Degner/R/')
mfile <- commandArgs(trailingOnly = TRUE)
#filelist <- scan(pipe("du /g/furlong/project/24_TSSCAGE/analysis/MAPPABILITY/MAP_BIN/DNAmapRAL_* | awk '$1 > 1000000 {print $2}'"), what='char')
mlist=sapply(chrs.keep,function(chr){#for each chr
    mmchrfile=paste0(mmfile, paste0('/',chr,'.memmap'))#define the file for each chr
    test <- mmap(mmchrfile, mode=int64())#get the accessor
    Rle(test[])#access it and create an Rle
})
mlist = as( mlist , 'SimpleRleList' )#
if(any(sum(mlist)==0)){next}#skip empty or nearly files
mlist=mlist>=0#get logical Rle List
export(mlist,file=paste0(getwd(),'/',basename(arg),'lrle.object.R')
