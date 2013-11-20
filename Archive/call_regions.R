source('src/generate_Rle.R')
source('src/generate_chromdata_cage.R')
source('src/load_annotations.R')

bin.size=50
bin.size=500
# Call our enriched regions -----------------------------------------------
#filter windows near genes etc.
act<-list(
  transcripts.gr,
  cad3.gr,
#  lincRNA.gr,
  crmgrs
)

act<-sapply(act,  function(gr){ mcols(gr) <-NULL;gr})
act<-do.call('c',act)
act<-resize(act,width(act)+(2*intergenic_dist),fix='center')
act<-reduce(act)
act<-keepSeqlevels(act,chrs.keep)

#now we select windows of size w

bin_ranges <- IRangesList(#initially just create disjoint regions 
  sapply(seqlevels(si)[!grepl('gl|hap', seqlevels(si))], function(tmp_chr) {
    IRanges(start = seq(1, seqlengths(si)[[tmp_chr]] - bin.size, bin.size), width = bin.size)
  })
)
bin_ranges<-resize(bin_ranges,width=500,fix='start')#expand them to size so they overlap
bin_ranges<-as(bin_ranges,'GRanges')#to GR
bin_ranges<-keepSeqlevels(bin_ranges,chrs.keep)#get rid of weird chrs
bin_ranges<-bin_ranges[countOverlaps(bin_ranges,act)==0]#get rid of those near our predefined active regions

#we cycle through our windows, say 10000 at a time
#we can start by eliminating all of the ranges which don't have a certain max read number, say 3
load(file.alltags.rpgc)

cut=428

v<-Views(alltags.rpgc$pos,as(bin_ranges,'RangesList'))
vn<-Views(alltags.rpgc$neg,as(bin_ranges,'RangesList'))
bin_ranges$maxpos<-unlist(viewSums(v))
bin_ranges$maxneg<-unlist(viewSums(vn))
#give each one a max score pick best strand
bin_ranges<-bin_ranges[bin_ranges$maxpos>cut |bin_ranges$maxneg>cut]
# 
# #just checking on the reduce functionality
# bin_ranges<-bin_ranges[countOverlaps(bin_ranges,GRanges('chr2L',IRanges(706000,760000)))>0]
# bin_ranges_reduced<-bin_ranges_reduced[countOverlaps(bin_ranges_reduced,GRanges('chr2L',IRanges(706000,760000)))>0]



#we can now view our intergenic bins in IGV
bin_ranges_reduced<-reduce(bin_ranges)
bin_ranges_reduced$max<-pmax(unlist(viewSums(Views(alltags$pos,as(bin_ranges_reduced,'RangesList')))),unlist(viewSums(Views(alltags$neg,as(bin_ranges_reduced,'RangesList')))))
bin_ranges_reduced$name<-paste0('intbin_',1:length(bin_ranges_reduced))
tail(bin_ranges_reduced[order(bin_ranges_reduced$max)],n=50)
export(bin_ranges_reduced,'analysis/call_regions/intergenic_bins.bed')




mean(countOverlaps(lincRNA.gr,bin_ranges_reduced)>0)


windchunk=1:10000
current_ranges<-bin_ranges[windchunk] 

#idea 1 - simply count rpgc tags, cutoff with a simple cutoff

#idea 2 - simply set a single cutoff, then count the number of lines with/without tags

#idea 3 - individual cutoffs.

#parameters of our call
fname=paste0('called_regions_',methodstring,'_',cutoffstring,'.bed')

#export the regions, finaly
export(called.windows,con=paste0('analysis/call_regions/enriched_windows/',fname))

#now count how many overlap DHS, and the lincRNAs