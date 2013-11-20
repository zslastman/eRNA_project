#Normalization-
# If it's really just a library size difference between our different data sets we should expect to see a very good
# correlation between the heights of the tallest peaks and the library size
# if on the other hand there's a signal/noise difference we'll see another factor, relating ot how much of the
# library is bullshit. We can Use the number of lone tags as a proxy for this....


source('src/tss_cage_functions.R')
source('src/generate_Rle.R')
source('src/load_annotations.R')

# calculate size of 'intergenic library' ----------------------------------
#get rle with the intergenic genome
transcripts.gr<-keepSeqlevels(transcripts.gr,chrs.keep)
int<-coverage(resize(transcripts.gr,width=width(transcripts.gr)+1000,fix='center'))
int<-int==0

accession.df$intergenic.lib.size <-sapply(accession.df$accession,function(acc){
sum(sum(cage.tag.rles[[acc]]$neg[int]+cage.tag.rles[[acc]]$pos[int]))
})
tmp<-sapply(accession.df$accession,function(acc){
  sum(sum(cage.tag.rles[[acc]]$neg[!int]+cage.tag.rles[[acc]]$pos[!int]))
})

#calculate size of library at all TSS
tss<-coverage(resize(tss.gr,width=500,fix='center'))>0
tss<-tss[chrs.keep]
accession.df$tss.lib.size <-sapply(accession.df$accession,function(acc){
  sum(sum(cage.tag.rles[[acc]]$neg[tss]+cage.tag.rles[[acc]]$pos[tss]))
})
save(accession.df,file=file.accession.df)
exit()

attach(accession.df)



ggplot(accession.df,aes(x=accession.df$library.size,y=tss.lib.size,color=line))+
  geom_point()+
  scale_x_continuous(limits=c(0,20000000))+
  scale_y_continuous(limits=c(0,20000000))+
  guides(col=guide_legend(ncol=4))+
  ggtitle('Library Size, vs. Total Tags in TSS ')
  
ggplot(accession.df,aes(x=library.size,y=intergenic.lib.size,color=line))+
  geom_point()+
  scale_x_continuous(limits=c(0,20000000))+
  scale_y_continuous(limits=c(0,100000))+
  guides(col=guide_legend(ncol=4))+
  ggtitle('Library Size, vs. Total Tags in intergenic space')

# ggplot(accession.df,aes(x=library.size,y=intergenic.lib.size,color=line))+
#   geom_point()+
#   scale_x_continuous(limits=c(0,20000000))+
#   scale_y_continuous(limits=c(0,20000000))+
#   guides(col=guide_legend(ncol=4))+
#   ggtitle('Total Tags in Intergenic Space vs. Total Tags in TSS')



multilines <- names(table(accession.df$line))[table(accession.df$line)>2]
cg<-cage.tag.rles[]

#get a matrix with the all
r<-tss.gr
r<-resize(r,width=500,fix='center')
r<-as(r,'RangesList')[ chrs.keep]
#now go through all our libraries and get a matrix of the totals at those regions
#cage.tag.rles[[2]]->rlel
all.p<-sapply(cg,function(rlel){
  unlist(viewSums(Views(rlel$pos,r)))  
})
all.n<-sapply(cg,function(rlel){
  unlist(viewSums(Views(rlel$neg,r)))  
})
#now we have a matrix with rows as TSS, columns as lines, and cells as tag counts
all<-all.p+all.n

line=multilines[1]
accs<-line.accessions[grepl(line,line.accessions)]
accs<-accs[!grepl('reseq',accs)]
M<-log2(all[,accs[1]])-log2(all[,accs[2]])
A<-0.5*(log2(all[,accs[1]])+log2(all[,accs[2]]))
plot(A,M)

#cumulative desnity of reads vs 





library(ggplot2)
a.df<-accession.df[ accession.df$line%in%multilines & !accession.df$reseq,]
#now go through from line to line and ask if the variability is larger 
for(line in multilines){
  accs<-a.df[ a.df$line==line, ]$accession 
  v<-sd(3,2)/2.5
}





