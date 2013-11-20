

#  concatenate RleLists for easy access ------------------------
load(file.cage.tag.rles)
#cage.tag.rles<-cage.tag.rles[eachline]
lines.cat.pos<-unlist(sapply(cage.tag.rles, "[[", 1))
lines.cat.neg<-unlist(sapply(cage.tag.rles, "[[", 2))
#concatenate them
for(n in seq_along(lines.cat.pos)){
  names(lines.cat.pos[[n]])<-NULL
  lines.cat.pos[[n]]<-do.call('c',as.list(lines.cat.pos[[n]]))
  #lines.cat.pos[[n]]<-Matrix(lines.cat.pos[[n]],sparse=T)
  
  
  names(lines.cat.neg[[n]])<-NULL
  lines.cat.neg[[n]]<-do.call('c',as.list(lines.cat.neg[[n]]))
  #lines.cat.neg[[n]]<-Matrix(lines.cat.neg[[n]],sparse=T)
  
}
cg.cat<-mapply('c',lines.cat.pos,lines.cat.neg)

#now acf on a selection of lines

accs<-accession.df$accession[c(1,2,78,79)]
acc<-accs[1]
acf.tmp<-acf(as.vector(as.numeric(cg.cat[[acc]]>0)))
acc<-accs[4]
acf.tmp2<-acf(as.vector(cg.cat[[acc]]))
plot(acf.tmp,xlim=c(0,200),main='cross correlation of Cage signal for Library 3E')


#also produce Jack's cumulative density plots
cat<-cg.cat[[1]]
m<-sapply(cg.cat,function(cat){
  cat<-rev(sort(cat[cat!=0]))
  rle<-cumsum(cat)/sum(cat)
  s<-seq(from=1,to=length(rle),length.out=1000)
   as.vector(rle[s])                       
})

matplot(m,type='l')
