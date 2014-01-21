#set of methods for working with lots of tracks at once


# concatenate our tag data, make sparse matrices --------------------------


load(file.cage.tag.rles)
cage.tag.rles<-cage.tag.rles[eachline]
lines.cat.pos<-unlist(sapply(cage.tag.rles, "[[", 1))
lines.cat.neg<-unlist(sapply(cage.tag.rles, "[[", 2))
rm(cage.tag.rles)

for(n in seq_along(lines.cat.pos)){
  names(lines.cat.pos[[n]])<-NULL
  lines.cat.pos[[n]]<-do.call('c',as.list(lines.cat.pos[[n]]))
  lines.cat.pos[[n]]<-Matrix(lines.cat.pos[[n]],sparse=T)
  
  
  names(lines.cat.neg[[n]])<-NULL
  lines.cat.neg[[n]]<-do.call('c',as.list(lines.cat.neg[[n]]))
  lines.cat.neg[[n]]<-Matrix(lines.cat.neg[[n]],sparse=T)
  
}

#make sparse matrices with rows as samples
#sparse matrix approach won't work
lines.cat.pos<-do.call(rBind,lines.cat.pos)
lines.cat.neg<-do.call(rBind,lines.cat.neg)






# try with big.matrix instead ---------------------------------------------


install.packages('bigmemory',lib='/g/furlong/Harnett/R/')
install.packages('bigmemory.sri',lib='/g/furlong/Harnett/R/')
install.packages('BH',lib='/g/furlong/Harnett/R/')
library(BH)
library(bigmemory)
library('bigmemory',lib.loc='/g/furlong/Harnett/R/')


big.matrix((tmp),nrow=1,ncol=length(tmp))



load(file.cage.tag.rles)
cage.tag.rles<-cage.tag.rles[eachline]
lines.cat.pos<-unlist(sapply(cage.tag.rles, "[[", 1))
lines.cat.neg<-unlist(sapply(cage.tag.rles, "[[", 2))
rm(cage.tag.rles)

for(n in seq_along(lines.cat.pos)){
  names(lines.cat.pos[[n]])<-NULL
  lines.cat.pos[[n]]<-do.call('c',as.list(lines.cat.pos[[n]]))
  lines.cat.pos[[n]]<-Matrix(lines.cat.pos[[n]],sparse=T)
  
  
  names(lines.cat.neg[[n]])<-NULL
  lines.cat.neg[[n]]<-do.call('c',as.list(lines.cat.neg[[n]]))
  lines.cat.neg[[n]]<-Matrix(lines.cat.neg[[n]],sparse=T)
  
}

#make sparse matrices with rows as samples
#sparse matrix approach won't work
lines.cat.pos<-do.call(rBind,lines.cat.pos)
lines.cat.neg<-do.call(rBind,lines.cat.neg)






c.sum<-c(0,cumsum(seqlengths(si)[chrs.keep])[])
c.sum<-c.sum[-length(c.sum)]
names(c.sum)<-chrs.keep
#we now have a named vector we can use to convert to and from the coordinates of
#the concatenated chromosomes







# Concatenate our tag rles ------------------------------------------------






#lets check the positioning is okay
crm8008.gr[7009]



RleViewsList(rleList = cage.tag.rles[[1]], rangesList = crm8008.gr[1:3])


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

lines.cat.pos<-do.call(rBind,lines.cat.pos)
lines.cat.neg<-do.call(rBind,lines.cat.neg)

lines.cat<-mapply(c,lines.cat.pos,lines.cat.neg)

lines.cat[[1]]




# Function to convert GRange objects to the concatenated genome -----------


# function to carry out with out concatenated objects ---------------------



#we have to kill the names in order to concatenate rle lists
function()

do.call('c',lapply(cage.tag.rles$split1A_732_68h$pos,function(x)x))

do.call('c',))

(
  (runsum(cage.tag.rles[[1]][[1]]==17,k=5))
  + (cage.tag.rles[[1]][[1]]==17 ))


  
  
  

# concatenation runs afoul of memory problems -----------------------------



# sparse matrix.. ---------------------------------------------------------
library(Matrix)

Matrix(lines.cat.pos[[n]],sparse=T)

tmp<-lapply(names(cage.tag.rles),function(x){
  Matrix(unlist(cage.tag.rles[[x]][['pos']],use.names=F),ncol=1,sparse=T)})
#now create a sparse matrix for all of that, with our lines as columns
posmat<-do.call(cBind,tmp)


write.table(x=lines.cat.neg[[1]],'tmp.txt')




install.packages('gdsfmt',lib='/g/furlong/Harnett/R/')
library('gdsfmt',lib.loc='/g/furlong/Harnett/R/')
  )





install.packages('GWASTools',lib='/g/furlong/Harnett/R/')
library('gdsfmt',lib.loc='/g/furlong/Harnett/R/')
)


source("http://bioconductor.org/biocLite.R")
biocLite(pkgs='GWASTools',lib.loc='~/Harnett/R')
library('GWASTools',lib.loc='~/Harnett/R/')
