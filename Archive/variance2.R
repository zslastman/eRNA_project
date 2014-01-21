#Variance analysis on cage data
source('src/generate_chromdata_cage.R')
source('src/load_crm_data.R')
#library(Matrix)
#library(zoo)
library(reshape2)
library(ggplot2)
load(file.cage.tag.rles)
load(file.alltags)
#load(file.cage.tag.rles.rpgc)
load(file.accession.df)


w<-1
replim=0.95


# Picking a simple dataset for each line ----------------------------------
acc.tn.rep<-list()

#for each line with 3 replicates
multilines<-names(table(accession.df$line))[table(accession.df$line)==4]




RleViewsList(rleList=cage.tag.rles[[1]],rangesList=crm8008.gr)


repmats<-mclapply(multilines,mc.cores=10,function(line){
  #get the cage tag rles
  set.df<-accession.df[grepl(x=accession.df$accession,line),]
  #exclude reseqs
  set.df<-set.df[ set.df$reseq==F,]
  accs<-as.character(set.df$accession)#data sets to use
  #now transform to a list of logical vectors so we can sum.
  cage.tag.rles.windowed.log<-rapply(cage.tag.rles[accs],how='list',f=function(rl)runsum(rl,w)!=0)
  
  
  #thisi s a little tricky- do.call doesn't work properly on a list of SimpleRleLists
  all.windowed.log<-cage.tag.rles.windowed.log[[1]]
  #so we initilialize from the firs 
  all.windowed.log[['both']]<-cage.tag.rles.windowed.log[[1]][['pos']]|cage.tag.rles.windowed.log[[1]][['neg']]
  #then loop over the rest
  all.windowed.log[['pos']]<-cage.tag.rles.windowed.log[[1]][['pos']]+cage.tag.rles.windowed.log[[2]][['pos']]
  all.windowed.log[['neg']]<-cage.tag.rles.windowed.log[[1]][['neg']]+cage.tag.rles.windowed.log[[2]][['neg']]
  all.windowed.log[['both']]<-all.windowed.log[['pos']]+  all.windowed.log[['neg']]  
#   
#   lapply(cage.tag.rles.windowed.log[-1],function(x){
#     x=cage.tag.rles.windowed.log[[2]]
#     #incrementing by 1 whenever a line has tags present
#     all.windowed.log[['pos']]<<-x[['pos']]+all.windowed.log[['pos']]
#     all.windowed.log[['neg']]<<-x[['neg']]+all.windowed.log[['neg']]
#     #not that we don't want to count lines twice if they have plus and negative strand activity
#     all.windowed.log[['both']]<<-(x[['pos']] | x[['neg']])+all.windowed.log[['both']]
#     NULL
#   })
  
  
  #now find the frequency of reproducibility for each tag height
  blanktable<-rep(0,times=)
  names(blanktable)<-as.character(1:length(accs))
  
  
  sapply(simplify=F,accs,function(acc){#for each line
    m<-sapply(1:25,function(tn){#return a matrix whose columns correspond to tag hieghts
      #and whose rows correspond to the reproduciblity ('1','2', etc.)
      t=colSums(table(
        all.windowed.log[['pos']][cage.tag.rles[[acc]][['pos']]==tn]
        ))[names(blanktable)]
      names(t)=names(blanktable)
      t[t%in%NA]=0
      t
    })
    
   
  })
  
})
repmats<-unlist(repmats,recursive=F)

#we want to find the column of m for which some percentage are fully reproducible
  # fullrep=nrow(set.df)#will depend on if we include the bad run
  #rep.fraction<-m[fullrep,]/colSums(m)#portion reproducible in each column
  #which(rep.fraction>replim)[1]#minimum tag height meeting reproducibility standards
cutoffs<- sapply(repmats,function(repmat){
   fullrep=nrow(repmat)#will depend on if we include the bad run
   rep.fraction<-repmat[fullrep,]/colSums(repmat)#portion reproducible in each column
   which(rep.fraction>replim)[1]#minimum tag height meeting reproducibility standards
})

nonrep<-sapply(repmats,function(m)m[1,]/colSums(m))
matplot(nonrep,type='l',xlab='Number of Tags',ylab='Fraction not reproducible')
title('fraction of nonreproducible tagsites by tagnumber')



#calculate cutoffs seperately
accession.df$cutoff<-unlist(cutoffs)[ accession.df$accession ]

accession.df.tmp<-accession.df[!accession.df$cutoff%in%NA,]
plot(accession.df.tmp$cutoff,accession.df.tmp$library.size,
     main='library size vs. 95% reproducibility limit',ylab='library size',xlab='Min tag height with 95% repr')

ggplot(accession.df.tmp,aes(x=cutoff,y=library.size,color=line))+geom_point(size=4)







# Trying to set window for 'reproducibility' -------------------------------


dists<-mclapply(multilines,mc.cores=10,function(line){
  #get the cage tag rles
  set.df<-accession.df[grepl(x=accession.df$accession,line),]
  #exclude reseqs
  set.df<-set.df[ set.df$reseq==F,]
  accs<-as.character(set.df$accession)#data sets to use
  
  pos1<-as(cage.tag.rles[[accs[[1]]]]$pos,'GRanges')
  pos2<-as(cage.tag.rles[[accs[[2]]]]$pos,'GRanges')
  neg1<-as(cage.tag.rles[[accs[[1]]]]$neg,'GRanges')
  neg2<-as(cage.tag.rles[[accs[[2]]]]$neg,'GRanges')
  
  pos1<-pos1[pos1$score!=0]
  pos2<-pos2[pos2$score!=0]
  neg1<-neg1[neg1$score!=0]
  neg2<-neg2[neg2$score!=0]
  
  pos1$dist<-distanceToNearest(pos1,pos2)$distance
  pos2$dist<-distanceToNearest(pos2,pos1)$distance
  neg1$dist<-distanceToNearest(neg1,neg2)$distance
  neg2$dist<-distanceToNearest(neg2,neg1)$distance
  
  d=pos1[pos1$score==1]$dist
  plot(
    hist(d[d<200],breaks=200)
      )
  
  
  (pos1[pos1$score==1]$dist,xlim=c(0,250))
  
  #thisi s a little tricky- do.call doesn't work properly on a list of SimpleRleLists
  all.windowed.log<-cage.tag.rles.windowed.log[[1]]
  #so we initilialize from the firs 
  all.windowed.log[['both']]<-cage.tag.rles.windowed.log[[1]][['pos']]|cage.tag.rles.windowed.log[[1]][['neg']]
  #then loop over the rest
  
  
  all.windowed.log[['pos']]<-cage.tag.rles.windowed.log[[1]][['pos']]+cage.tag.rles.windowed.log[[2]][['pos']]
  all.windowed.log[['neg']]<-cage.tag.rles.windowed.log[[1]][['neg']]+cage.tag.rles.windowed.log[[2]][['neg']]
  all.windowed.log[['both']]<-all.windowed.log[['pos']]+  all.windowed.log[['neg']]  
  #   
  #   lapply(cage.tag.rles.windowed.log[-1],function(x){
  #     x=cage.tag.rles.windowed.log[[2]]
  #     #incrementing by 1 whenever a line has tags present
  #     all.windowed.log[['pos']]<<-x[['pos']]+all.windowed.log[['pos']]
  #     all.windowed.log[['neg']]<<-x[['neg']]+all.windowed.log[['neg']]
  #     #not that we don't want to count lines twice if they have plus and negative strand activity
  #     all.windowed.log[['both']]<<-(x[['pos']] | x[['neg']])+all.windowed.log[['both']]
  #     NULL
  #   })
  
  
  #now find the frequency of reproducibility for each tag height
  blanktable<-rep(0,times=length(accs))
  names(blanktable)<-as.character(1:length(accs))
  
  
  sapply(simplify=F,accs,function(acc){#for each line
    m<-sapply(1:25,function(tn){#return a matrix whose columns correspond to tag hieghts
      #and whose rows correspond to the reproduciblity ('1','2', etc.)
      t=colSums(table(all.windowed.log[['pos']][cage.tag.rles[[acc]][['pos']]==tn]))[names(blanktable)]
      names(t)=names(blanktable)
      t[t%in%NA]=0
      t
    })
    
    
  })
  
})





