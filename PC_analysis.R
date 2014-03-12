##A Load libraries,functions,data
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
#load the annotation data with scores
load('data/objects/scored_regions.object.R')
load('data/objects/gene.annotation.object.R')
load('data/objects/accession.df.object.R')
#load('data/objects/cg.mapfilt.pl.object.R')
outfolder = 'analysis/PC_analysis'
dir.create(outfolder,showWarnings=F)





#see this page
#http://stats.stackexchange.com/questions/57467/how-to-perform-dimension-reduction-after-doing-pca-in-r
#simulate some data to make sure this is working properly
#we'll have four columns, with the second two being dependant on the first two

c1 = rbinom(n=1000,1,0.5)
c2 = rbinom(n=1000,1,0.5)
c1c = rbinom(n=1000,1,0.95*c1)
c2c = rbinom(n=1000,1,0.95*c2)
#somewhat laboriously print the correlations for all pairwise combinations of the above
#hopefully this code comes in useful someday...
clist=list(c1,c2,c1c,c2c)
names(clist)=c('c1','c2','c1c','c2c')
combnames=apply(combn(names(clist),2),2,function(x){do.call(paste0,as.list(x))})
setNames(combnames,apply(combn(clist,2),2,function(x){cor(x[[1]],x[[2]])}))

mat = as.matrix(do.call(cbind,clist))
#using some copied code from stack exchange...
res <- prcomp(mat, center = F, scale = FALSE)
res$sdev
names(res)
length(res$sdev)
cumsum(res$sdev^2/sum(res$sdev^2))
res$sdev^2/sum(res$sdev^2)
res$rotation
dim(res$rotation)
res$x
dim(res$x)
pcs = split (pcob$loadings , col(pcob$scores)) #from matrix to list

# Save plot as pdf
pdf(file=paste0(outfolder,'/',"pctest.pdf"))
plot(cumsum(res$sdev^2/sum(res$sdev^2))) #cumulative explained variance
dev.off()





