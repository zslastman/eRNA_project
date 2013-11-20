
#  concatenate RleLists for easy access ------------------------
load(file.cage.tag.rles)
lines.cat.pos<-unlist(sapply(cage.tag.rles, "[[", 1))
lines.cat.neg<-unlist(sapply(cage.tag.rles, "[[", 2))
rm(cage.tag.rles)

for(n in seq_along(lines.cat.pos)){
  names(lines.cat.pos[[n]])<-NULL
  lines.cat.pos[[n]]<-do.call('c',as.list(lines.cat.pos[[n]]))
  #lines.cat.pos[[n]]<-Matrix(lines.cat.pos[[n]],sparse=T)
  
  
  names(lines.cat.neg[[n]])<-NULL
  lines.cat.neg[[n]]<-do.call('c',as.list(lines.cat.neg[[n]]))
  #lines.cat.neg[[n]]<-Matrix(lines.cat.neg[[n]],sparse=T)
  
}
multilines <- names(table(accession.df$line))[table(accession.df$line)>2]


lines.cat<-mapply(c,lines.cat.pos,lines.cat.neg)
n=2
macc1<-accession.df$accession[grepl( multilines[n],accession.df$accession)][1]
macc2<-accession.df$accession[grepl( multilines[n],accession.df$accession)][2]

l1<-lines.cat[[macc1]]
l2<-lines.cat[[macc2]]
s1<-accession.df$library.size[ accession.df$accession==macc1 ]
s2<-accession.df$library.size[ accession.df$accession==macc2 ]

#I have 1 tag. How often do I expect non reproduction, given a certain libray
#size?
actual.rep.rate<-mean(l2[l1==1]>0)
#now what's the expect rate for that?
#we can calculate the odds of a given true rate outputting a zero, but then we need to know the true rate
#for that we could in principle get a posterior for the true rate from our other rate
#for this I think we can just use a uniform prior


#this simulates the results 
score=2
res=1000
rate.post<-rgamma(res,score+1,1)
rate.post<-rate.post*0.2
sim<-sapply(rate.post,function(rate){rpois(res,rate)})
hist(rate.post,breaks=100)
#and here's our simulated rep. rate
mean(sim==0)




#build the predicted poisson
n=1
i=10000
r1=n/s1
r2=(r1*s1)/s2
x1=log2(n)
x2<-log2(rpois(i,r2))
m=x1-x2
a=0.5*(x1+x2)

plot(a,m)

names(cag)[grepl('7')]
