
load('data/objects/accession.df.object.R')
library(reshape2)
gn=200
gsigs=log(10^(seq(1,3,length.out=gn)))
#vector of library sizes
reptable=accession.df
reptable=reptable[reptable$tissue=='embryo',]
libsizes=reptable$library.size
samplegroups=c('timepoint','line','collection','prep','seq')
reptable=reptable[,c('acc','seq','prep','collection','line','timepoint')]

#vector of noise magnitudes for each different factor (seq, techrep, collection(biorep),line)
noisemags=c('seq'=0.1,'prep'=0.1,'collection'=0.2,'line'=0.5,'timepoint'=1)
#create matrix of counts
simcountmat=matrix(gsigs,nrow=length(gsigs),ncol=length(libsizes))
colnames(simcountmat)=reptable$acc
#for(reptype in names(noisemags) ){

for(reptype in samplegroups ){
  SD=sqrt(noisemags[[reptype]])
  repvect=reptable[,reptype]
  for(repn in unique(repvect)){
    simcountmat[,repvect==repn] = rnorm(n=gn,
                                        mean=simcountmat[,repvect==repn],
                                        sd=SD
    )
  }
}
colnames(simcountmat)=reptable$acc
head(sim.melt)
sim.melt = melt(simcountmat)
colnames(sim.melt) <- c('Gene','lib','value')
sim.melt$Gene = as.factor(sim.melt$Gene)
var= 'line'

for(var in samplegroups){
  sim.melt[[var]] = reptable[[var]][ match(sim.melt$lib,reptable$acc) ]
}

#let's fit a series of progressively more complicated models....
library(nlme)
attach(sim.melt)

lmfit = lm(formula = value ~ as.factor(Gene) * timepoint,
           data=sim.melt)

lmefit = lme(fixed = value ~ Gene * timepoint,
  random = ~ timepoint | Gene,
  data=sim.melt)

lmefit = lme(fixed = )
