
#checking out RNAP at crms
rlecovplot(rlelist1=cage.tag.rles[[2]][[2]]*100,rlelist2=chrom.rles.rpgc.sub.merge[[1]]*100,g=resize(crmgrs[1],20000,fix='center'))

rlecovplot(rlelist1=cage.tag.rles[[2]][[2]]*100,rlelist2=chrom.rles.rpgc.sub.merge[[1]]*100,g=resize(crmgrs[1],20000,fix='center'))



#okay, so charles suggests plotting the distribution of reproducibility for peaks with peak height 1,2 etc.
#do I do this only once for the lines, for all lines, and with or without rpgc?
#I can do it over 1 line to begin with,later maybe over the averaged rpgcs
#may can later do this for high and low RNA-seq zones, illustrating the need for background modelling

#parameters here are
#1)The distance we expect clusters to move between our datasets
w<-10#for now

#so we can perhaps check this

#get the top 100 highest windows from our cage data
chr<-'chr2R'
cutoff<-sort(unique(as.vector(alltags$pos[[chr]])),decreasing=T)[5]
toppeaks<-slice(alltags$pos[[chr]],lower=cutoff)
toppeaks<-GRanges(seqnames=chr,IRanges(start(toppeaks),end(toppeaks)),name=paste0('toppeak',1:length(toppeaks)),seqinfo=si)

logpnratio<-abs(log2(pnratio))

!(logpnratio) %in% c(NaN,Inf) && logpnratio
length(allsums)
#can just pick one



