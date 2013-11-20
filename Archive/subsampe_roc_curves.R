load(file.alltags)

#function that takes an RleList and returns a new one with the counts resampled
#according to binomial probabilities 
sampleRleList<-function(srle,p){
  as(
    sapply(names(srle),function(chr){      
      chr<-srle[[chr]]
      sites<-chr>0
      subsample<-rbinom(n=length(as.vector(chr[sites])),prob=p,size=as.vector(chr[sites]))
      chr[sites]<-subsample
      chr
    }),'SimpleRleList')
}

sample.alltags<-sapply(simplify=F,names(alltags),function(strand){
  sampleRleList(alltags[[strand]],0.2)
})
