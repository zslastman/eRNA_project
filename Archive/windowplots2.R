setwd('~/Harnett/TSS_CAGE_myfolder/')
rm(list=ls())

source('src/tss_cage_functions.R')
library(reshape2)
library(ggplot2)
library(ROCR,lib.loc='~/Harnett/R')

rpgc.step.num<-100
pval.step.num<-20

#load the window table as 'm'
load('windowtable.robject.R')
load(file.accession.df)
load(file.cage.tag.rles)
cg<-cage.tag.rles
regs=list(
  pos=import(con='analysis/positive_8008.bed',asRangedData=F,seqinfo=si),
  neg=import(con='analysis/negative_8008.bed',asRangedData=F,seqinfo=si),
  int=import(con='analysis/random.intergenic.bed',asRangedData=F,seqinfo=si),
  cad3pos=import(con='analysis/cad3.pos.bed',asRangedData=F,seqinfo=si),
  cad3neg=import(con='analysis/cad3.neg.bed',asRangedData=F,seqinfo=si)
)




# Plot of distance between best windows of techreps -----------------------
# 
# dists.df<-sapply(simplify=F,unique(z$windowsize),function(w){
#     z<-z[ z$windowsize ==w ,]
#     #for each replicated line
#     sapply(simplify=F,1:length(multiaccs[[1]]),function(n){
#       #get our distance
#       dists<-
#         (z[ z$acc==multiaccs[[1]][[n]] ,]$pos)-
#         (z[ z$acc==multiaccs[[2]][[n]] ,]$pos)
#       dists<-dists[!is.na(dists)]
#     })
#   })
# #melt for ggplot2
# d<-melt(dists.df)
# h(d)
# colnames(d)<-c('value','acc','windowsize')
# 
# #strip the ws of our windowsizes
# d$value<-abs(d$value)
# #now we plot the distance
# ggplot(d,aes(x=factor(windowsize),y=as.numeric(value)))+stat_summary(fun.data='mean_cl_boot')+
#   scale_y_continuous(name='distance in bp between center of best windows')+
#   scale_x_discrete(name='Size of best Window')+
#   ggtitle('mean distance between the best windows in Replicates')
# 


# names(dists.df[[1]})

# h(d)  
# dim(dists.df)
# d$value


# reproducibility plots for the technically replicated lines --------------

#put the position entires in a new column
for(l in 1:length(multiaccs[[1]])){
  
  accs=sapply(multiaccs,'[[',l)
  tmp<- z[z$acc %in% accs[1],]
  tmp$val2<-z[z$acc %in% accs[2],]$val.libnorm
  #excludeing outliers
  tmp<-tmp[ tmp$val.libnorm<1000 & tmp$val2<1000,]
  title=paste0('Correspondance between technical replicates for line ',strsplit(multiaccs[[1]][l],'_')[[1]][2])
  ggsave(h=25,w=25,filename=paste0(accs[1],' tmp.jpeg'),
         qplot(data=tmp,x=val.libnorm,y=val2,facets=windowsize~region,main=title,geom='jitter')+geom_abline(intercept=0)
  )
  mean(tmp[ tmp$value<3,]$val2==0)
  
}
