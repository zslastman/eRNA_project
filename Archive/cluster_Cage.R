#hierarchically cluster cage peaks
source('src/tss_cage_functions.R')


clusts.list<-list()

#for the best library
clusts.big.5<-Hcluster.cage.rle(cutoff=5,cage.tag.rles$split3E_375_68h$pos+cage.tag.rles$split3E_375_68h$neg)
clusts.small.5<-Hcluster.cage.rle(cutoff=5,cage.tag.rles$split7J_491_68h$pos+cage.tag.rles$split7J_491_68h$neg)

clusts.big.10<-Hcluster.cage.rle(cutoff=10,cage.tag.rles$split3E_375_68h$pos+cage.tag.rles$split3E_375_68h$neg)
clusts.small.10<-Hcluster.cage.rle(cutoff=10,cage.tag.rles$split7J_491_68h$pos+cage.tag.rles$split7J_491_68h$neg)

clustlist<-list()
for(multiline in multilines){
  #get the cage tag rles
  set.df<-accession.df[grepl(x=accession.df$accession,multiline),]
  #exclude reseqs
  set.df<-set.df[ set.df$reseq==F,]
  accs<-as.character(set.df$accession)#data sets to use 
  
  clustlist[[accs[[1]]]]<-Hcluster.cage.rle(cutoff=0,cage.tag.rles[[accs[[1]]]]$pos+cage.tag.rles[[accs[[1]]]]$neg)
  clustlist[[accs[[2]]]]<-Hcluster.cage.rle(cutoff=0,cage.tag.rles[[accs[[2]]]]$pos+cage.tag.rles[[accs[[2]]]]$neg)
  
}



# cluster on a single line ------------------------------------------------


#for viewing
export(clustlist[[1]],'hclusters_test.bed')
export(cage.tag.rles[[1]]$pos,'hclusts_cagesig_pos.bedGraph')
export(cage.tag.rles[[1]]$neg,'hclusts_cagesig_neg.bedGraph')

#Now filter
tmp<-clustlist[[1]]
tmp<-tmp[width(tmp)!=1]
tmp<-tmp[ tmp$score>3 ]
export(tmp,'hclusters_test_filt.bed')



# Track describing Biological Reproducibility -----------------------------

#####Creat a track describing biological reproducibility
w=1
#tracks describing presence or absence of 2 or more cage tags
cage.tag.rles.windowed.log<-rapply(cage.tag.rles,how='list',f=function(rl)runsum(rl,w)>1)
cage.tag.rles.windowed.log<-cage.tag.rles.windowed.log[eachline]
#sum them
rep.rle<-cage.tag.rles.windowed.log[[1]]
lapply(cage.tag.rles.windowed.log[-1],function(x){
  rep.rle[['pos']]<<-x[['pos']]+rep.rle[['pos']]
  rep.rle[['neg']]<<-x[['neg']]+rep.rle[['neg']]
  NULL
})
#and export
export(rep.rle[['pos']],'LineCount.pos.bedGraph')
export(rep.rle[['neg']],'LineCount.neg.bedGraph')




# cluster on all ----------------------------------------------------------
allclust.pos=Hcluster.cage.rle(cutoff=0,alltags$pos)
allclust.neg=Hcluster.cage.rle(cutoff=0,alltags$neg)
allclust.pos
allclust.neg
export(allclust.pos,'hclusts_alllines_pos.bed')
export(allclust.neg,'hclusts_alllines_neg.bed')

#filter using reproducibility
Views(rep.rle$pos,as(allclust.pos,'RangesList'))

