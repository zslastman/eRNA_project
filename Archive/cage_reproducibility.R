setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/generate_Rle.R')
source('src/generate_chromdata_cage.R')
source('src/load_annotations.R')

# look at reproducibility in different classes of sequence ----------------

##define active tss for the first class
tss.gr$active68<-countOverlaps(tss.gr,active68.gr)>0

##define a group of intergenic CRMs with low H3K4me3
#have to remove metadata to concatenate them
crms.filtered =  crms.all
#put the names back in
#crms.filtered$name<-c(crm8008.gr$id,cad3.gr$name)
crms.filtered$intergenic<-distanceToNearest(crms.filtered,transcripts.gr)$distance>intergenic_dist

#filter out those near genes
crms.filtered<-crms.filtered[crms.filtered$intergenic]

#now filter for H3K4me3 using just peaks, for now
modK4me3<-chrompeaks.modencode[['K4me3_4-8h']]
modK4me3<- modK4me3[modK4me3$FDR...  < 0.50]#accept even low confidence peaks

crms.filtered$H3K4me3_peak<-
  0< (countOverlaps(crms.filtered,modK4me3)+
  countOverlaps(crms.filtered,chrompeaks[['K4me3_6-8h']])+
  countOverlaps(crms.filtered,chrompeaks[['K4me3_4-6h']]))

#we can also do it with the actual signal and views...



##define our 'intergenic' windows - regions outside of genes, crms, K4me3 peaks,

intergenic.regions <-  gaps(do.call(c,lapply(
  list(
  crm8008.gr,
  cad3.gr,
  transcripts.gr,
  chrompeaks.modencode[['K4me3_4-8h']],
  chrompeaks[['K4me3_6-8h']],
  chrompeaks[['K4me3_4-6h']]),
  function(gr){mcols(gr)<-NULL;gr})))

#and our 'inactive' windows, which are the same, only without even K4me1
intergenic.inactive.regions <- gaps(do.call(c,lapply(
  list(
    crm8008.gr,
    cad3.gr,
    transcripts.gr,
    chrompeaks.modencode[['K4me3_4-8h']],
    chrompeaks[['K4me3_6-8h']],
    chrompeaks[['K4me3_4-6h']],
    chrompeaks[['K4me1_4-6h']],
    chrompeaks[['K4me1_6-8h']],
    chrompeaks.modencode[['K4me1_4-8h']]
    ),
  function(gr){mcols(gr)<-NULL;gr})))

genome.classlist<-list(
  activetss=resize(tss.gr[tss.gr$active68],width=500,fix='center'),
 # crms=crms.filtered[!crms.filtered$H3K4me3_peak],
   crms=crms.filtered[],
  
  intergenic=intergenic.regions,
  intergenic.inactive=intergenic.inactive.regions
)

grles<-sapply(genome.classlist,function(x)(coverage(x)>0)[chrs.keep])





# Assessing variability within regions ------------------------------------


# Technical Reproducibility -----------------------------------------------
#choose the lines with 2 replicates

multilines <- names(table(accession.df$line))[table(accession.df$line)>2]

#gclass=grles[[1]]
bigchrs<-chrs.keep[1:6]
heightrange=1:25
dists=seq(1,21,2)

repmats<-mclapply(multilines,mc.cores=10,function(mline){
    dists=seq(1,21,2)
    #get the cage tag rles
    set.df<-accession.df[grepl(x=accession.df$accession,mline),]
    #exclude reseqs
    set.df<-set.df[ set.df$reseq==F,]
    accs<-as.character(set.df$accession)#data sets to use
    
    #we want to know when a location is 'reproducible' 
    #what we do is cycle through different window sizes, creating a smoothed rle list each time,
    #we annotate a 'reproducibility' rle with the first window size it becomes reproducible on
    #those which are never reproducible will take a value of -1 when we subtract 1
    
    #make srlelist of zeros
    repdist.pos<-cage.tag.rles[[1]][[1]]==-1
    repdist.neg<-cage.tag.rles[[1]][[1]]==-1
    
    for(w in dists){
      wind.log<-rapply(cage.tag.rles[accs],how='list',f=function(rl)runsum(rl,w,e='c')!=0)
      
      rep.pos<-wind.log[[1]][['pos']]&wind.log[[2]][['pos']]
      rep.neg<-wind.log[[1]][['neg']]&wind.log[[2]][['neg']]
      
      repdist.pos[repdist.pos==0 & rep.pos] <- w
      repdist.neg[repdist.neg==0 & rep.neg] <- w
    }
    #change to show the distance to the nearest, rather than the window size
    repdist.neg<-floor((repdist.neg-1)/2)
    repdist.pos<-floor((repdist.pos-1)/2)
    dists<-(dists-1)/2
    
    
    #now we have a simple rle list with a number over each tag position
#     #telling us how close the nearest tag in the replicate is
#       sapply(simplify=F,accs,function(acc){
#         pos.p<-cage.tag.rles[[acc]]$pos<=20  #locations of our tags
#         pos.n<-cage.tag.rles[[acc]]$neg<=20
#         suppressWarnings({vals<-unlist((cage.tag.rles[[acc]]$pos[pos.p])[ bigchrs ])})
#         suppressWarnings({reps<-unlist((cage.tag.rles[[acc]]$pos[pos.p])[ bigchrs ])})
#         list(
#           tagnum=suppressWarnings({vals<-c(vals,unlist((cage.tag.rles[[acc]]$neg[pos.p])[ bigchrs ]))}),
#         repdist=suppressWarnings({reps<-c(reps,unlist((cage.tag.rles[[acc]]$neg[pos.p])[ bigchrs ]))})
#         )
#       })
#     
   # this produces tables, but maybe better off with the rles for now
    #and whose rows correspond to the reproduciblity ('1','2', etc.)
    sapply(simplify=F,grles,function(gclass){
      sapply(simplify=F,accs,function(acc){#for each line
        m<-sapply(heightrange,function(tn){#return a matrix whose columns correspond to tag hieghts
          tab=table(repdist.pos[cage.tag.rles[[acc]]$pos==tn & gclass])#only those in our region
          if(length(tab)==0){tab=rep(0,length(dists))
          }else{
            tmp<-colnames(tab)
            tab=as.matrix(tab[bigchrs,])
            colnames(tab)<-tmp
            tab=colSums(tab)[as.character(dists-1)]#problems here, if table has only one columns
          }
          names(tab)=as.character(dists-1)
          tab[tab%in%NA]=0
          tab
      })
    })
  })
})
    
names(repmats)<-multilines
# Now produce reproducibility curves for each class and for each  ---------
#we go through each distance



cutoffs<-list()
nonrepnum<-list()

dists=(dists-1)/2

for(gclass in names(genome.classlist)){
  cutoffs[[gclass]]<-list()
  nonrepnum[[gclass]]<-list()
  for(d in 0:9){
  m<-sapply(repmats,function(repmatl){
    repmatl<-repmatl[[gclass]][[2]]
    rep<-rownames(repmatl)[ rownames(repmatl) %in% 0:d  ]
    rep<-repmatl[rep,]
    if(is.matrix(rep)){rep<-colSums(rep)}
    1-(rep/colSums(repmatl))
  })
  colnames(m)<-multilines
  main.tit=paste('Technical Reproducibility at distance',d,'for',gclass)
  
  
  jpeg(paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/plots/reproducibility/',main.tit,'.jpeg'),height=1200,width=1920)
  
  matplot(m,type='l',main=main.tit,cex.main=3,col=rainbow(ncol(m)),ylim=c(0,1))
  abline(h=0.05)
  cutoffs[[gclass]][[as.character(d)]]<-apply(m,2,function(v){which(v<0.05)[1]})
  legend(x='topright',legend=paste(colnames(m),'-',cutoffs[[gclass]][[as.character(d)]]),fill=rainbow(ncol(m)),cex=3)
  dev.off()
  
}
}

# get total number of nonreproducible tags --------------------------------
for(gclass in names(genome.classlist)){
  nonrepnum[[gclass]]<-  sapply(  0:9, function(d){
    m<-sapply(repmats,function(repmatl){
      repmatl<-repmatl[[gclass]][[2]]
      nonrep<-!rownames(repmatl) %in% 0:d  
      nonrep<-repmatl[nonrep,]
     
    })
   colSums(m)
  })

}




# look at cutoffs ---------------------------------------------------------

# each lib not line -------------------------------------------------------
nonrepnum<-list()

mt<-lapply(names(genome.classlist),function(gclass){
  sapply(0:9,function(d){
    m<-lapply(repmats,function(repmat){
      #each library, not each line
      repmatl<-repmat[[gclass]][[2]]
      rep<-rownames(repmatl)[ rownames(repmatl) %in% 0:d  ]
      rep<-repmatl[rep,]
      if(is.matrix(rep)){rep<-colSums(rep)}
      a=1-(rep/colSums(repmatl))
      repmatl1<-repmat[[gclass]][[1]]
      rep<-rownames(repmatl1)[ rownames(repmatl1) %in% 0:d  ]
      rep<-repmatl1[rep,]
      if(is.matrix(rep)){rep<-colSums(rep)}
      b=1-(rep/colSums(repmatl1))
      m<-cbind(b,a)
      colnames(m)<-names(repmat[[gclass]])
      
      v<-list(sum(repmatl1[!rownames(repmatl1) %in% 0:d ,]),sum(repmatl[!rownames(repmatl) %in% 0:d ,]))
      names(v)<-names(repmat[[gclass]])
      v
      #return(m)
    })
    
    unlist(unname(m))
  })
})

Non_Reproducible_Sites<-rowSums(sapply(mt,function(m)m[,1]))











