# @ author Dermot Harnett, EMBL Heidelberg
# @date 17/5/2013
# @title normalize_by_RNASEQ.R
########################################
###Script to subtract out the background from our cage data
#So in Hoskins et al they normalized their cage data using the RNASEQ data from the encode project
#We model our CAGE data as coming from a mixture of RNASEQ data and the signal
#We estimate the  mixture proportion empirically first, then assign explicity probabilites to
#of our bins, estimate the mixture again, iterations shouldn't really help


load(cage.rle.file)



yl<-c(0,max(cage.sig))
yllog<-c(0,log10(yl[2]))
#first simuate:
#our cage data looks like this
cage.sig<-rep(0,100)#zeroes
cage.sig[c(1,25,50,75)]<-c(3000,6000,8000,6000)#with some very high peaks
plot(cage.sig,ylim=yl)
cage.sig<-Rle(cage.sig)
noise.p<-sapply( exp(rnorm(100)),function(x){rpois(1,x)})
true.mix<-0.2# we want to recover this quantity
plot(noise.p*true.mix,ylim=yl)
cage.sim<-(cage.sig)+noise.p*true.mix

true.lambda=sum(true.mix*noise.p)/sum(cage.sig+true.mix*noise.p)

#here's our simulated dataset,we also know the true noise
plot(log10(as.vector(cage.sim)))

N=sum(cage.sim)
#now, state the function to minimize
noiseprobs<-true.noise/sum(true.noise)
est.mix<-function(m,data=cage.sim,noise=noise.p){sum(abs(m*noise-data))}
#now optimize(works)
est.m<-optimize(f=est.mix,interval=c(0,1))
est.m<-est.m$minimum
#get our estimated lambda
est.lambda<-sum(est.m*noise.p)/sum(cage.sim)

plot(as.vector(noise.p*est.mix))

#now go through our vector, calculate likelihood of the total at each bin
#for each bin, we get the beta distribution of bin frequencies





P-dat1

mix.exp <- function(c(x)){1/sum(abs((rep(x,3)*c(3,2,4))-c(9,6,120)))}
mix.exp<-function(x){3*x*}

plot(dat1)






X <- cbind(1, runif(1000))
theta.true <- c(2,4,6) # error variance = 2, intercept = 4, slope = 6.
y <- X %*% theta.true[-1] + sqrt(theta.true[1]) * rnorm(1000)






ols.lf1 <- function(theta, y, X) {
  beta <- theta[-1]
  sigma2 <- theta[1]
  if (sigma2 <= 0) return(NA)
  n <- nrow(X)
  e <- y - X%*%beta                                  # t() = matrix transpose
  logl <- ((-n/2)*log(2*pi)) - ((n/2)*log(sigma2)) - ((t(e)%*%e)/(2*sigma2))
  return(-logl) # since optim() does minimisation by default.
}
optim(c(0.5,6,8), method="L-BFGS-B", fn=ols.lf1,
      lower=c(1e-6,-Inf,-Inf), upper=rep(Inf,3), y=y, X=X)

optimize(c(-4,4),f= function(x){abs(x^2-4)})
?optimize

?optim
D(trig.exp,'x')

D({3*x^5},x)