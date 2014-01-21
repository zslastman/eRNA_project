#system for creating mmap files

setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
load('data/objects/cg.pl.map.object.R')





#we want a function that creates a big 

#we want to first creat a csv file for some rle objects (later maybe also GRanges?)

#needs to be in a system of fodlers for each chr



#maybe use bigmemory?



library(bigmemory)
load('data/objects/cg.pl.map')
cg.pl.map-cg.pl.map[1:2]
lines.cat.pos<-unlist(sapply(cg.pl.map, "[[", 1))
lines.cat.neg<-unlist(sapply(cg.pl.map, "[[", 2))
rm(cg.pl.map)
for(n in seq_along(lines.cat.pos)){
  names(lines.cat.pos[[n]])<-NULL
  lines.cat.pos[[n]]<-do.call('c',as.list(lines.cat.pos[[n]]))
  lines.cat.pos[[n]]<-Matrix(lines.cat.pos[[n]],sparse=T)
  
  
  names(lines.cat.neg[[n]])<-NULL
  lines.cat.neg[[n]]<-do.call('c',as.list(lines.cat.neg[[n]]))
  lines.cat.neg[[n]]<-Matrix(lines.cat.neg[[n]],sparse=T)
  
}

#first create a suitably huge filebacked bigmatrix

###THis won't work - the matrices are too big to fit all of the cage data on them.
tracks=cg.pl.map
tracks=unlist(tracks)
stopifnot(sapply(tracks,function(x)'SimpleRleList' %in% is(x)))
stopifnot(sapply(tracks,function(x)names(x) == chrs.keep))
stopifnot(!grepl('.both$',names(tracks)))
tracknum=length(tracks)
cagemats=sapply(chrs.keep,function(chr){
	matfile=paste0(chr,'.genome.matrix.bin')
	descfile=paste0(chr,'.genome.matrix.desc')
	chrlength=seqlengths(si)[chr]
	cmat=filebacked.big.matrix(type='double',nrow=chrlength,ncol=tracknum,dimnames=c(NULL,names(tracks)),
	backingfile=matfile,
	descriptorfile=descfile)
	for(i in seq_along(tracks)){
		cmat[,i]=tracks[[i]][[chr]]
	}
	scmat=Matrix(sparse=T,cmat[,],ncol=tracknum)

	file.remove(matfile)
	scmat
})

tracks=cg.pl.map
tracks=unlist(tracks)
stopifnot(sapply(tracks,function(x)'SimpleRleList' %in% is(x)))
stopifnot(sapply(tracks,function(x)names(x) == chrs.keep))
stopifnot(!grepl('.both$',names(tracks)))
tracknum=length(tracks)
cagemats=sapply(chrs.keep,function(chr){
	chrlength=seqlengths(si)[chr]
	cmat=Matrix(sparse=T,0,nrow=chrlength,ncol=tracknum)
	for(i in seq_along(tracks)){
		cmat[,i]=tracks[[i]][[chr]]
	}
	file.remove(matfile)
	scmat
})

#Let's try a gigantic sparse matrix instead


# #check our matrix works
# fromtrack=tracks[[1]][[chr]][(start(loc)-10):(end(loc)+10)]
# frommat=cmat[(start(loc)-10):(end(loc)+10),1]
# fromsmat=scmat[(start(loc)-10):(end(loc)+10),1]
# find a location with some info in one of our tracks
# loc=slice(lower=100,tracks[[1]])[[1]][1]
# as.vector(fromtrack)
# as.vector(fromsmat)
# frommat
