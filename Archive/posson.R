
# A more probabalistic approach... ----------------------------------------


# First we test if the two distributions are different in a line ----------

# But this test would be underpowered in one line -------------------------


#optimizing a poisson

simdata<-rpois(n=300000,0.25)

lf<-function(x,sdata=simdata){sum(dpois(sdata,x,log=T))}
lf(0.1)
lf(0.25)
lf(0.3)
#optim is behaving badly

optim(c(0.1,0.5,1),lf,lower=0.0000001,upper=1,method='CG')
optim(c(0.1,0.5,1),lf,lower=0.0000001,upper=1,method='L-BFGS-B')


printlocation<-function(gr){paste0(seqnames(gr),':',start(gr),'-',end(gr))}

intergenic<-import('analysis/intergenic.inactive.regions.bed')
intvect<-c(unlist(unname(p)),unlist(unname(n)))
ilen<-length(intvect)
tab<-table(intvect)
vals<-as.numeric(names(tab))
zero_one_r<-(tab[2])/(tab[1]+tab[2])

r<-0.25
d<-rpois(100000,r)
tab<-table(d)
mean(d)

zero_one_r<-(tab[2])/(tab[1]+tab[2])
zero_one_r
a<-dpois(1,x)
b<-dpois(0,lambda=x)

lfp<-function(x,zone_ratio=zero_one_r){ abs(  (dpois(1,x)/dpois(0,x))-zone_ratio   ) }
lfp(0.23)
lfp(0.19)

optimise(lf,interval=c(0,3),maximum=F)

curve(dpois(x,lambda=3),0,3)



p[mappability.rle==1]
length(mappability.rle[[1]])
length(cage.tag.rles[[1]][[1]][[1]])
sapply(intergenic,length)
sapply(mappability.rle,length)
sapply(cage.tag.rles[[1]][[1]],length)

length(intergenic[[1]])

length(p[[1]])

       length(mappability.rle[[1]])

n=cage.tag.rles[[1]][[2]][intergenic]

