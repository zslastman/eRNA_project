


#we COULD also calculate the library size like they do in DESeq - 
#DESeq uses the median of the ratios, between the samples counts at the (tss say)
#and the geometric mean across samples

cm<-matrix(rnbinom(size=15,prob=0.6,n=100000),ncol=100)
plot(density(cm))
libsizes=runif(100,1,10)
plot(libsizes)
cm<-t(apply(cm,1,function(v){v*libsizes}))

regnum<-nrow(cm)
getBestWindowMat(reg=gr[10000:11000,],w=100,chrs.keep,bigchrs,cage=cage.tag.rles)
cm<-gr$cagemat[10000:11000,]
#now we get the geometric mean of each row
gmean<-function(x){exp(mean(log(x)))}
#row vector of geometric means across sampels for different regions
gmeans<-apply(cm,1,gmean)
#now get each genes ratio to the median across all sapmles
ratios<-t(sapply(1:regnum,function(n){cm[n,]/gmeans[n]}))
#Now for each library(column) get the median ratio
scales<-apply(ratios,2,median)

#Now we can work out the variance for each of our eRNAs (I think)
#Now we can derived negative binomdial distributions fo reach one
#Which allows us to derive a null distribution over the correlation?
#We 

m<-cagecountmatlist