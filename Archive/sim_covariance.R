#simulate eRNA QTLs

x<-runif(1000,0,100)
y<-x+rnorm(1000,0,100)
plot(x,y)
cor(x,y)

corrnorm<-function(n,rho,mean1,mean2,sd1,sd2){
  # desired correlation = cos(angle)
  theta <- acos(rho)             # corresponding angle
  x1    <- rnorm(n, mean1, sd1)        # fixed given data
  x2    <- rnorm(n, mean2, sd2)      # new random data
  X     <- cbind(x1, x2)         # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
  
  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
  
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  return(list(x1, x))  
  
}





#define techrepped lines
multiaccs<-sort(accession.df$accession)[1:40]
multiaccs<-as.character(multiaccs[!grepl('reseq',multiaccs)])
multiaccs<-list(multiaccs[1:10],multiaccs[11:20])

#our null to begin with - a poisson

floor(10^(rnorm(1000)))
sample(1:1000,size=10000,prob=)

#gamma is our 
curve({dgamma(x,20,30)},from=0,to=20)
rgamma(1000,)



#paramaters
#dispersion of distirubiton 0=poisson
#effect size distribution
#read numbers
#number of candidate genes/enhancers

#generate bunch of signals for our crms

#matrix of crms
crmmat<-matrix(0,ncol=linenum,nrow=crmnum)
#matrix of TSS
tssmat<-matrix(0,ncol=linenum,nrow=crmnum)
#to encode crm tss linkages we have a list of integers - the index each crm is linked to
#for now just have 1 linkage
crm2tsslist<-list()
for(crm in 1:nrow(crmmat)){crm2tsslist[[crm]]<-crm}
#we populate our list of crms with data. We do this by sampling from 

#first we generate a distribution of eRNA signals. We can later do this realistically
#for now we just assume a normal distribution of eRNA and TSS signal

crmn<-100
tssn<-100
crmmean<-60
tssmean<-6000
#our null hypothesis is that the vector of TSS and CRM
libsizes=accession.df$library.size
#values to adjust rates by
libadj=libsizes/max(libsizes)
#underlying isgnals of eRNAs
crmsigs=as.integer(rpois(crmn,crmmean))
#underlying signals of TSS
tssigs=as.integer(rpois(tssn,tssmean))

#now add correlated fracion

####__
plot(rpois(100,1000),ylim=c(0,1200))
plot(rnorm(100,1000,sd=sqrt(1000)),ylim=c(0,1200))

#generate data - rows as librarys, columns as libraries
crmmat<-sapply(crmsigs,function(sig){sapply(libadj,function(lam){rpois(1,sig*lam)})})
tssmat<-sapply(tssigs,function(sig){sapply(libadj,function(lam){rpois(1,sig*lam)})})
dim(crmmat)
#now normalize
crmmat.n<-t(vapply(1:nrow(crmmat),function(linenum){crmmat[linenum,]/libsizes[linenum]},FUN.VALUE=crmmat[1,]))
tssmat.n<-t(vapply(1:nrow(tssmat),function(linenum){tssmat[linenum,]/libsizes[linenum]},FUN.VALUE=tssmat[1,]))

pvals<-sapply(1:nrow(crmmat),function(crm){
  sapply(crm2tsslist[[crm]],function(tss){
    tmp<-lm(tssmat.n[tss,]~crmmat.n[crm,])#perform linear regression on eRNAs and their TSS
    x<-summary(tmp,)#get the summary
    pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)#get pvalue off that
    })
})
#so we're getting what look a lot like uniform pvalues.
plot(density(pvals),xlim=0:1)



# Now repeat the above but with signal ------------------------------------

prop.cor<-0.5
cor.dist
#our null hypothesis is that the vector of TSS and CRM
libsizes=accession.df$library.size
#values to adjust rates by
libadj=libsizes/max(libsizes)
#underlying isgnals of eRNAs
crmsigs=as.integer(rnorm(100,60,sd=10))
#underlying signals of TSS
tssigs=as.integer(rnorm(100,6000,sd=1000))

#generate data - rows as librarys, columns as libraries
crmmat<-sapply(crmsigs,function(sig){sapply(libadj,function(lam){rpois(1,sig*lam)})})
tssmat<-sapply(tssigs,function(sig){sapply(libadj,function(lam){rpois(1,sig*lam)})})
dim(crmmat)
#now normalize

crmmat.n<-t(vapply(1:nrow(crmmat),function(linenum){crmmat[linenum,]/libsizes[linenum]},FUN.VALUE=crmmat[1,]))
tssmat.n<-t(vapply(1:nrow(tssmat),function(linenum){tssmat[linenum,]/libsizes[linenum]},FUN.VALUE=tssmat[1,]))

#find d2nearest for each crm
#assign the closest TSS for each crm
#now remove those from the 
tss
tss<-coverage(resize(tss.gr,width=500,fix='center'))>0
tss<-tss[chrs.keep]


#for every element of crm2tsslist, do a regression from that crm to those tss.

pvals<-sapply(1:nrow(crmmat),function(crm){
  sapply(crm2tsslist[[crm]],function(tss){
    tmp<-lm(tssmat.n[tss,]~crmmat.n[crm,])#perform linear regression on eRNAs and their TSS
    x<-summary(tmp,)#get the summary
    pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)#get pvalue off that
  })
})

#so we're getting what look a lot like uniform pvalues.
plot(density(pvals),xlim=0:1)


  #now we do our linear regressions
  
  reps=1000
  maxsig=200
  m<-a<-matrix(NA,nrow=maxsig,ncol=reps)
  sigs=c(1:maxsig)
  libadj=c(1,0.8)
  
  for(i in 1:reps){
  
    mat<-sapply(sigs,function(sig){sapply(libadj,function(lam){rpois(1,sig*lam)})})
    mat<-log10(mat)
    m[,i]<-mat[1,]-mat[2,]
    a[,i]<-0.5*(mat[1,]+mat[2,])
  
  }

qplot(as.vector(a),as.vector(m))


















#define techrepped lines
multiaccs<-sort(accession.df$accession)[1:40]
multiaccs<-as.character(multiaccs[!grepl('reseq',multiaccs)])
multiaccs<-list(multiaccs[1:10],multiaccs[11:20])




#playing with glm.nb
sigs<-runif(100,1,20)
sigs2<-runif()
out<-rnbinom(100,mu=sigs,size=sigs^1.2)

