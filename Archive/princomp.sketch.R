#PCA


#function extract.pcas
subtract.pcas<-function(m,num){
 stopifnot(is.matrix(m)) 
 stopifnot(is.numeric(m))
 
 #we want to tget pcas to reduce the variability between columns
 p=prcomp(t(m))
 
 for(i in 1:ncol(m)){
   for(j in 1:num){   
      m[,i] = p$x[i,j]*p$rotation[,j]
      
   }
 }
 m
}

#test feature extraction
#matrix to test on
m=matrix(ncol=10,rnorm(300,10))
#columns account for equal amount of each row
qplot(data=melt(m),x=factor(Var2),y=value,geom='boxplot')
#Add systematic variation
for(i in 1:3){m[,i]= m[,i]+(30*10)}
#columns 1:3 now account for too muhc of the variationi in each row
qplot(data=melt(m),x=factor(Var2),y=value,geom='boxplot')

#subtract it out
m=subtract.pcas(m,2)
qplot(data=melt(m),x=factor(Var2),y=value,geom='boxplot')


dim(m)

dim(m)
p=prcomp(t(m))
p$rotation #matrix where column j describes PC j
p$x #matrix where column j describes amount of PCj in column i of the matrix

plot(p$rotation[,1])

#each principle component is available as a column of p$loadings

#p$scores[1,1] tells us how muhc of the first component is in the first row


p$scores[3,1]*p$loadings[,1]



melt(m)
