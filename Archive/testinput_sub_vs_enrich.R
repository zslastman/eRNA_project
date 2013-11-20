#testing input subtraction and binning

input.sig<-1+rnorm(400)
input.sig[100:300]<-input.sig[100:300]+3
sig=rnorm(400)
sig[100:300]<-sig[100:300]+9

points(input.sig)
plot(sig)

library(zoo)
enrich<-sig/input.sig
enrich.bin<-
plot(enrich)

dev.off()


binvector<-function(v,binsize){
  l=length(v)
  binnum<-floor(l/binsize)
  binstarts<-(seq(1:binnum)-1)*binsize+1
  sapply(binstarts,function(binstart)v[binstart:(binstart+binsize-1)])
}

n=10000

input<-c(rnorm(n),0+rnorm(n))
sig<-c(rnorm(n),0+rnorm(n))
input<-input+abs(min(c(input,sig)))
sig<-sig+abs(min(c(input,sig)))
plot(input,ylim=c(-3,20))
plot(sig,ylim=c(-3,20))
sub<-sig-input
enr<-sig/input
b.input<-colMeans(binvector(input,10))
b.sig<-colMeans(binvector(sig,10))
b.sub<-(b.sig-b.input)
b.en<-(b.sig/b.input)


plot(sub)
plot(enr)
plot(b.sub)
plot(b.en)

sd(input)
sd(sig)
sd(sub)
sd(enr)
sd(b.input)
sd(b.sig)
sd(b.sub)
sd(b.en)