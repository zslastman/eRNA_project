# choosing regions to look at ---------------------------------------------

#get the CAGE tag info
#need to reorder the chromsomes I think...

gr<-resize(crms.all,width=500,fix='center')
gr<-sort(gr)
strandviews<-sapply(alltags,function(strand){Views(strand,as(gr,'RangesList'))[unique(seqnames(gr))]})
gr$strandsums<-sapply(strandviews,function(strand){unlist(viewSums(strand))})
gr$strandmaxs<-sapply(strandviews,function(strand){unlist(viewMaxs(strand))})
gr$cagesum<-rowSums(sapply(strandviews,function(strand){unlist(viewSums(strand))}))
gr<-gr[rev(order(gr$cagesum))]
gr$K27ac<-countOverlaps(gr,chrompeaks[['K27Ac_6-8h']])>0
gr$PolII<-countOverlaps(gr,chrompeaks[['PolII_6-8h']])>0
gr$K4me3<-
  0< (countOverlaps(gr,modK4me3)+
        countOverlaps(gr,chrompeaks[['K4me3_6-8h']]))
gr$bimod<-rowMax(gr$strandsums)/rowSums(gr$strandsums)
gr$bimod[gr$bimod==NaN]<-NA
#for enhancers
gr<-gr[gr$intergenic]
gr<-gr[gr$K27ac&gr$PolII]


#now looking at the more bimodal ones
gr<-gr[ (gr$bimod<0.75) %in%T  ]
gr$name

#for genes
#get the GENIC ones
gr<-gr[!gr$intergenic]
#ones not DNase
gr<-gr[!grepl(pattern='dnase',gr$name)]
#ones with many TRs
gr<-gr[countOverlaps(gr,tss.gr)>3]


r=gr[1]