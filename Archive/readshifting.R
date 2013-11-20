
#I need a strategy for spawning parallel processes from Rstudio

setwd('~/Harnett/TSS_CAGE_myfolder/')
load('./Rimage')




save('.Rimage')

inputsubpol<-import('/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/wig-merged/InputSubtracted/PolII_B-PolII-Rpb3-B51_6-8h_input-subtracted.wig.gz')
rpgcpol<-import('/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/wig-merged/RPGC/PolII_B-PolII-Rpb3-B51_6-8h.wig.gz')
rpgcinput<-import('/g/tier2/furlong/delhomme/projects/11_Histone_ChIP_seq/BiTS-ChIP/alignment/6-8h/wig-merged/RPGC/input_D-3_6-8h.wig.gz')




#tibors comverage
load('/g/tier2/furlong/lab_resource/histmod/DataVector_68h_1bp.rda')
#delete unnecessary ones
rm(list=ls(pattern='[^It]_68h'))

#now load the rpgc vectors
load('/g/tier2/furlong/lab_resource/histmod/DataVector_68h_10bp_RPGC.rda')

#now compare to our ones just made from tibors, and from my own data


#now do binning

tmp1<-chrom.rles$PolII_6.8_R1[['chr2L']][10000:10100]
tmp2<-PolII_68h_R1_B_10bp[['chr2L']][1000:1010]#+PolII_68h_R2_B51_10bp[['chr2L']][1000:1010])/2


which(!PolII_68h_R1_B_1bp[['chr3L']][1:100000]==0)[1]
which(!as.vector(chrom.rles$PolII_6.8_R1[['chr3L']][1:100000]==0))[1]

as.vector(chrom_nctl.rlelists$PolII_6.8_R1$chr2L[4998:5098])

#and compare again.









# #trying to compare to tibors
# t.g<-as(PolII_68h_R1_B_1bp,'SimpleRleList')
# tiborchrs<-names(gcov)[names(gcov) %in% names(t.g)]
# gl<-as(g,'RangesList')[tiborchrs]
# 
# t.g<-t.g[names(gcov)[names(gcov) %in% names(t.g)]]
# 
# 
# gcov<-coverage()
# 
# #plotting our coverage
#  cutoff<-
# 
#  peaks <- peakSummary(g.peaks)
#  peaks<-peaks[1:10]
# 
#  peak.depths <- viewMaxs(g.peaks)
#  cov.pos <- coverage(g[strand(g) == "+"])
#  cov.neg <- coverage(g[strand(g) == "-"])
#  peaks.pos <- Views(cov.pos, ranges(g.peaks))
#  peaks.neg <- Views(cov.neg, ranges(g.peaks))
#  wpeaks <- tail(order(peak.depths[['chr2L']]), 4)
#  coverageplot(peaks.pos[['chr2L']][wpeaks[4]], peaks.neg[['chr2L']][wpeaks[4]])




#bin our coverage vector
tmp<-as(viewSums(Views(chrom.rles$PolII_6.8,bins$chr2L)),'SimpleRleList')


# 
# 
# chr<-'chr2L'
# i<-1
# s=10000
# e=20000
# plotrange<-ranges(GRanges(chr,IRanges(s,e)))

g.peaks <- slice(coverage(g), lower = peakCutoff(gcov, fdr = 0.00001) )

rlelist1<-gcov
rlelist2<-gscov

g<-g.noshift
g<-g.shifted
g<-g.resize
#get coverage on strands 
rlelist1<- coverage(g[strand(g)=='+'])
rlelist2<- coverage(g[strand(g)=='-'])
#plot the views
coverageplot(opposite=T,
   Views(rlelist1[[chr]],ranges(GRanges(chr,IRanges(22958793,22962230)))),
   Views(rlelist2[[chr]],ranges(GRanges(chr,IRanges(22958793,22962793)))))


coverageplot(opposite=T,
             Views(rlelist1[[chr]],ranges(crmgrs[30])),
             Views(rlelist2[[chr]],ranges(crmgrs[30])))









?Views
?ranges

# #only those on tibors chromosomes
# g.peaks<-g.peaks[names(gcov)]
# peaks.mine <- Views(gcov, ranges(g.peaks))
# 
# g.peaks<-g.peaks[names(t.g)]
# peaks.tibors <- Views(t.g, ranges(g.peaks))
# 
# 
# 
# 
# coverageplot(peaks.mine[['chr2L']][wpeaks[4]], peaks.tibors[['chr2L']][wpeaks[4]])


coverageplot(peaks.mine[['chr2L']][wpeaks[4]], peaks.tibors[['chr2L']][wpeaks[4]])










# Functions for binning, under work ---------------------------------------


rlelist<-PolII_68h_10bp_Merged_ControlSubtracted

exprlelist<-mclapply(mc.cores=10,names(rlelist),function(chr){
  rle<-rlelist[[chr]]
  l<-length(rle)
  exprle<-t(matrix(as.vector(rle),ncol=10,nrow=l))[1:seqlengths(Dmelanogaster)[[chr]]]
  Rle(exprle)
})
names(exprlelist)<-names(rlelist)
seqlengths(crmgrs)
sapply(exprlelist,length)
sapply(PolII_68h_R1_B_1bp,length)
sapply(PolII_68h_10bp_Merged_ControlSubtracted,length)
seqlengths(Dmelanogaster)


length(PolII_68h_R1_B_1bp$chr2L)
length(PolII_68h_10bp_Merged_ControlSubtracted$chr2L)
length(exprle)


coverageplot(Views(rlelist=PolII_68h_R1_B_1bp,gr=crmgrs[2])$chr2R)
coverageplot(Views(as(PolII_68h_R1_B_1bp,'SimpleRleList'),as(crmgrs[2],'RangesList')))

seqlengths(crmgrs)