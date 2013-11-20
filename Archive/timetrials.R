
cage.tag.rles->cage.tag.rles.bak
cage.tag.rles<-cage.tag.rles.bak











#whats' the quickest way to pull out all the data for a specific location?
#How does views compare to pulling out with an rle list?, to just using locations?
#with the rle, does it make a difference how many chromosomes the tags are on?
#How many there are?

#create small object to pull out 
tmp<-cage.tag.rles[[1]][[1]]
tmp<-tmp==17



system.time(
{
  sapply(cage.tag.rles,function(srl){
    srl[[1]][[1]][80000:9000]
  })
})

system.time({
  for(srl in cage.tag.rles){
    srl[[1]][[1]][80000:9000]
    
  }
})