#load stuff
########################
setwd(dir='/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
library(VGAM,lib.loc='~/Harnett/R')
library(igraph,lib.loc='~/Harnett/R')

#load our mesodermal cage libraries
#this gives the four of them

#How much bigger are these guys? already have that info

#Do we sum them all - yes to begin with

#scatter plot of our meso data.vs our means int he other libraries
#to do this we can use our count cage list
mesocountmats=lapply(cagecountmatlist[[,'68']][1:4],function(regmat){rowSums(regmat[[1]])})

meso_crm_counts=rowSums(cagecountmatlist['crm8008','68'][[1]])
emb_crm_counts=rowSums(cagecountmatlist['crm8008','68h'][[1]])

qplot(meso_crm_counts,emb_crm_counts,log='xy')


meso_tss_counts=rowSums(cagecountmatlist['tss','68'][[1]])
emb_tss_counts=rowSums(cagecountmatlist['tss','68h'][[1]])

qplot(meso_tss_counts,emb_tss_counts,log='xy')


crm24=rowSums(cagecountmatlist['crm8008','24h'][[1]])
crm68=rowSums(cagecountmatlist['crm8008','1012h'][[1]])

qplot(crm24,crm68,log='xy')

load('data/objects/all.cage.unprocessed.object.R')
load(file.tss)
load('/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/allcage.object.R')
load('data/objects/accession.df.full.object.R')
tss.gr=sort(tss.gr)


mesosums=unlist(viewSums(GRViews(allcage[['meso68h']]$both,tss.gr)))
embsums=unlist(viewSums(GRViews(allcage[['embryo68h']]$both,tss.gr)))


# mesosums=mesosums[mesosums>quantile(mesosums,0.1)]
# embryosums=mesosums[mesosums>quantile(mesosums,0.1)]


mesoranks=rank(mesosums)
embranks=rank(embsums)

library(LSD)

heatscatter(mesoranks,embranks,main='Quantile normalized Cage 6-8hours, Embryo vs. Mesoderm',xlab='Quantile Normalized Mesodermal Cage',ylab='Quantile Normalized Whole Embryo Cage')
#comparison with a large emb library.
replicate(3,{
	
	s=sample(1:sum(accession.df$timepoint=='68h'),1)
	subset = cage.tag.rles[ accession.df$timepoint=='68h'][[s[1]]]$both
	for(i in s[-1]){	subset=subset+cage.tag.rles[ accession.df$timepoint=='68h'][[i]]$both}
	subembsum=unlist(viewSums(GRViews(subset,tss.gr)))
	subembranks=rank(subembsum)
	heatscatter(subembranks,embranks,main='Quantile normalized Cage 6-8hours, Embryo vs. Embryo subset',xlab='Quantile Normalized Subset of Mesodermal Cage',ylab='Quantile Normalized Whole Embryo Cage')
})
#What percentage of the eRNAs above our cutoff are in this biggest library?

#question - do we see more genuine eRNAs 

#What if we include only our crms with eRNAs in the 


#can we get a hold of some genes we know should 'NOT be in the mesoderm
bdgp<-read.csv2('/g/furlong/Harnett/data_general/bdgp_insitu.csv')