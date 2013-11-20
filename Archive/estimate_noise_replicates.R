#using the equtation from Balwierz et al and the normalized counts, let's try and calculate
#the noise component as they did in their paper.


#the full expression, with sigma 2 expanded
L<-function(sig,m,n,n.n,m.n){
	log(  ( 2*sig^2 + 1/n + 1/m )  ) + (		  (log(n.n/m.n)   ^   2)   /    (2*(2*  ( 2*sig^2 + 1/n + 1/m )  )   ) )
}
#test this expr with dummy data

f=0.001 #fraction in input pool
x=log(f)




#function to generate our tag distribution given a variance, a tag magnitude and a library size
###This works fine FOR SMALL SIGMA (less than ~ 0.5)

tag.dist<-function(n,sigma,x,N){

		exp( -(( log(n/N) - x )^2)     /  (2*(sigma^2+1/n)^2) ) / 
		(	n * sqrt(2 * pi) * (sigma^2+ 1/n )  )


}
#so with this we can generate dummy data by scanning through values of X
#generating data points with a frequency determined by the power law distribution
#then generating the actual tag numbers for those.
signals = 10^seq(1,10,by=0.01)
sig=0.085
libsizes=c(1,0.5)

#our model posits that signals are subject ot a multiplicative law and then a poisson law
signals.mulnoise=exp(log(signals)+rnorm(length(signals),0,sig))
signals.mulshot=rpois(length(signals.mulnoise),signals.mulnoise*libsizes[1])

signals.mulnoise.b=exp(log(signals)+rnorm(length(signals),0,sig))
signals.mulshot.b=rpois(length(signals.mulnoise),signals.mulnoise*libsizes[2])

N=sum(signals.mulshot)
M=sum(signals.mulshot.b)

sapply(seq(0.001,0.2,by=0.01),function(sigma){
	sum(mapply(sigma,signals.mulshot,signals.mulshot.b,signals.mulshot/N,signals.mulshot.b/M,FUN=L))
})

#this is outputting numbers that are way too low, peaking at around 4
tag.dist(4,4,log(0.001),1000000)
curve(tag.dist(x,4,log(0.001),1000000),from=0,to=2000)



#maximize correlations with enrico's data, with timepoints