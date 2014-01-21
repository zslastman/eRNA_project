setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/generate_chromdata_cage.R')
source('src/generate_Rle.R')
#source('src/normalization.R')
library(reshape2)
library(ggplot2)

#name cage data for convenience
cg<-cage.tag.rles
cgr<-cage.tag.rles.rpgc

# CAD3 SET ----------------------------------------------------------------
cad3.gr$name<-gsub('[#\\(\\)\\s/\\-]','',cad3.gr$name)
cad3.gr$name<-sub(' ','',cad3.gr$name)
cad3.pos<-cad3.gr[ cad3.df$active68 %in% T & cad3.df$intergenic ]
cad3.neg<-cad3.gr[ cad3.df$inactive68 %in% T & cad3.df$intergenic ]
cad3.pos<-keepSeqlevels(cad3.pos,chrs.keep)
cad3.neg<-keepSeqlevels(cad3.neg,chrs.keep)
cad3.neg<-cad3.neg[countOverlaps(cad3.neg,chrompeaks[['PolII_6-8h']])==0]

#maybe filter out GJ guys?
cad3.neg<-cad3.neg[ !grepl('GJ\\d+',cad3.neg$name) ]
cad3.pos<-cad3.pos[ !grepl('GJ\\d+',cad3.pos$name) ]
cad3.neg$name<-paste0('cadneg_',cad3.neg$name)
cad3.pos$name<-paste0('cadpos_',cad3.pos$name)


export(cad3.pos,con='analysis/make_regions_bedfiles/cad3.pos.bed')
export(cad3.neg,con='analysis/make_regions_bedfiles/cad3.neg.bed')
# import('analysis/make_regions_bedfiles/cad3.neg.bed')



# Now 8008 CRMs -----------------------------------------------------------

#crm8008.gr<-crm8008.gr[ countOverlaps(crm8008.gr,cad3.gr)==0 ]#get rid of ones overlapping hte CAD3 enhancers
#get rid of the intergenic ones
crm8008.gr$intergenic<-distanceToNearest(crm8008.gr,transcripts.gr)$distance>intergenic_dist
#And add the K4me3 info from our chromatin and the modencode projects'
modK4me3<-chrompeaks.modencode[['K4me3_4-8h']]
modK4me3<- modK4me3#accept even low confidence peaks
#mark those with any H3K4me3 overlap
crm8008.gr$H3K4me3_peak<- 0< (countOverlaps(crm8008.gr,modK4me3)+ countOverlaps(crm8008.gr,chrompeaks[['K4me3_6-8h']]))
#And add the K4me1 info from our chromatin and the modencode projects'
modK4me1<-chrompeaks.modencode[['K4me1_4-8h']]
# modK4me1<- modK4me1 #accept even low confidence peaks
modK4me1<- keepSeqlevels(modK4me1,chrs.keep)#
#mark those with H3K4me1 overlaps
crm8008.gr$H3K4me1_peak<- 0< (countOverlaps(crm8008.gr,modK4me1)+  countOverlaps(crm8008.gr,chrompeaks[['K4me1_6-8h']]))

#And add the K27ac info from our chromatin and the modencode projects'
modK27ac<-chrompeaks.modencode[['K27ac_4-8h']]
modK27ac<- modK27ac#accept even low confidence peaks
#mark those with H3K4me1 overlaps
crm8008.gr$H3K27ac<- 0< (countOverlaps(crm8008.gr,modK27ac)+  countOverlaps(crm8008.gr,chrompeaks[['K27Ac_6-8h']]))
#and our PolII
crm8008.gr$polII<-   0 < countOverlaps(crm8008.gr,chrompeaks[['PolII_6-8h']])

crm8008.gr$chrom.neg.set<-crm8008.gr$H3K4me1_peak & !crm8008.gr$H3K27ac & !crm8008.gr$H3K4me3_peak& crm8008.gr$intergenic & ! crm8008.gr$polII


#now construct an even more negative set by taking the bottom of the H3K27ac

sum.K27ac<-unlist(viewSums(Views(chrom.rles.rpgc.sub.merge[['H3K27ac_6.8']],as(crm8008.gr,'RangesList'))))
sum.K27.modencode<-unlist(viewSums(Views(chrom.rles.modencode[['H3K27ac']],as(crm8008.gr,'RangesList'))))

#all 8008
qplot(sum.K27ac,sum.K27.modencode,color=crm8008.gr$chrom.neg.set,shape=crm8008.gr$chrom.neg.set,xlab='Mesodermal K27ac',ylab='whole embryo ac',
      main='Comparison or our mesodermal K27ac data with\nModencode data - 6.8hrs, all 8008 crms')

neg.crms<-crm8008.gr[crm8008.gr$chrom.neg.set,]
neg.crms$islowk27<-sum.K27.modencode[crm8008.gr$chrom.neg.set]<quantile(sum.K27ac[crm8008.gr$chrom.neg.set],0.2)

#just negatives
qplot(sum.K27ac[crm8008.gr$chrom.neg.set],sum.K27.modencode[ crm8008.gr$chrom.neg.set],color=neg.crms$islowk27,xlab='Mesodermal K27ac',ylab='whole embryo ac',
      main='Comparison or our mesodermal K27ac data with\nModencode data - 6.8hrs, all negative 8008 crms')

#put the lowk27 status in as the score so we can export to a bed
neg.crms$name<-paste0('neg',1:length(neg.crms))
neg.crms$score<-as.numeric(neg.crms$islowk27)


# Now do the positives with more stringent FDR ----------------------------


#And add the K4me3 info from our chromatin and the modencode projects'
modK4me3<-chrompeaks.modencode[['K4me3_4-8h']]
modK4me3<- modK4me3#accept even low confidence peaks
#mark those with any H3K4me3 overlap
crm8008.gr$H3K4me3_peak<- 0< (countOverlaps(crm8008.gr,modK4me3)+ countOverlaps(crm8008.gr,chrompeaks[['K4me3_6-8h']]))

#And add the K4me1 info from our chromatin and the modencode projects'
modK4me1<-chrompeaks.modencode[['K4me1_4-8h']]
modK4me1<- modK4me3#accept even low confidence peaks
#mark those with H3K4me1 overlaps
crm8008.gr$H3K4me1_peak<- 0< (countOverlaps(crm8008.gr,modK4me1)+  countOverlaps(crm8008.gr,chrompeaks[['K4me1_6-8h']]))

#And add the K27ac info from our chromatin and the modencode projects'
modK27ac<-chrompeaks.modencode[['K27ac_4-8h']]
modK27ac<- modK27ac#accept even low confidence peaks

#mark those with H3K4me1 overlaps
crm8008.gr$H3K27ac<- 0< (countOverlaps(crm8008.gr,modK27ac)+  countOverlaps(crm8008.gr,chrompeaks[['K27Ac_6-8h']]))


#load sets
pos.crms<-crm8008.gr[crm8008.gr$H3K4me1_peak & crm8008.gr$H3K27ac & !crm8008.gr$H3K4me3_peak & crm8008.gr$intergenic,]

pos.crms<-keepSeqlevels(pos.crms,chrs.keep)
neg.crms<-keepSeqlevels(neg.crms,chrs.keep)



#resize them, for now, to 500bp
pos.crms<-resize(pos.crms,500,fix='center')
neg.crms<-resize(neg.crms,500,fix='center')
#possibly filter out the noninergenic ones again
#we'll also get rid on an outlier - neg 376
neg.crms<-neg.crms[-376]


#and our 'inactive' windows, which are the same, only without even K4me1
act<-list(
  crm8008.gr,
  cad3.gr,
  transcripts.gr,
  chrompeaks[['K4me3_6-8h']],
  chrompeaks[['K4me3_4-6h']],
  chrompeaks[['K4me1_4-6h']],
  chrompeaks[['K4me1_6-8h']],
  chrompeaks.modencode[['K4me3_4-8h']],
  chrompeaks.modencode[['K4me1_4-8h']],
  chrompeaks.modencode[['K27ac_4-8h']],
  chrompeaks[['PolII_6-8h']],
  chrompeaks[['PolII_4-6h']],
  lincRNA.gr
)
act<-sapply(act,  function(gr){ mcols(gr) <-NULL;gr})
act<-do.call('c',act)
act<-resize(act,width(act)+1000,fix='center')
act<-reduce(act)
act<-keepSeqlevels(act,chrs.keep)
strand(act)<-'*'
intergenic.inactive.regions<-gaps(act)
intergenic.inactive.regions<-intergenic.inactive.regions[strand(intergenic.inactive.regions)=='*']
#export tracks to make sure this all looks okay
export(transcripts.gr,con='analysis/make_regions_bedfiles/transcripts.gr.bed')
#export(gene.model.data,'analysis/make_regions_bedfiles/gene_models.bed')
export(act,con='analysis/make_regions_bedfiles/act.bed')
export(intergenic.inactive.regions,con='analysis/make_regions_bedfiles/intergenic.inactive.regions.bed')

#function that imports binned UCSC wig files and exports as an unbinned rle 
import.UCSC.binned<-function(trackfile){
  # Now use mappability and GC to filter out regions ------------------------
  trackfile<-import('~/Harnett/data_general/map_dm3_Bin10.wig.gz')#import as UCSC track
  trackfile<-trackfile[chrs.keep]#only good chrs
  seqinfo(trackfile)<-si
  #unbin
  trackfile<-coverage(x=trackfile,weight='score',width=sapply(cage.tag.rles[[1]][[1]],length))
  
}

#now we need to start picking random windows in our regions
#function to pick windows of size x in a set 
windowGRange<-function(gr,x){
  gr<-gr[width(gr)>=x]
  gr<-resize(gr,width(gr)-(x-1),fix='start')
  c<-coverage(gr)
  c<-which(c>0)
  c<-unlist(c)
  return(
    (GRanges(names(c),IRanges(c,width=x)))
  )
}

# Load mappability and gc, or make if they don't exist --------------------

file.mappability.rle<-'data/objects/mappability.rle.object.R'
file.gc.rle<-'data/objects/gc.rle.object.R'
if(!file.exists(file.mappability.rle)){
  mappability.rle<-import.UCSC.binned('~/Harnett/data_general/map_dm3_Bin10.wig.gz')
  save(mappability.rle,file=file.mappability.rle)
}else{load(file.mappability.rle)}


if(!file.exists(file.gc.rle)){
  gc.rle<-import.UCSC.binned('~/Harnett/data_general/gc_dm3_Bin10.wig.gz')
  save(gc.rle,file=file.gc.rle)
}else{load(file.gc.rle)}

int<-windowGRange(intergenic.inactive.regions,500)
int<-sample(int,size=20000)

#only mappable regions
int<-int[unlist(viewMeans(  Views(mappability.rle[chrs.keep],as(int,'RangesList')[chrs.keep])  ))==1]
#gc content
#int<-int[abs(unlist(viewMeans(  (Views(gc.rle[chrs.keep],as(int,'RangesList')[chrs.keep])  )))-0.5)<0.2 ]

#exclude polII
int<-keepSeqlevels(int,chrs.keep)
#get rid of the top 1%, likely to be anomalously high
#tmp<-unlist(viewSums(  (Views(alltags.rpgc$pos[chrs.keep],as(int,'RangesList')[chrs.keep]))))
#tmp.n<-unlist(viewSums(  (Views(alltags.rpgc$neg[chrs.keep],as(int,'RangesList')[chrs.keep]))))
#tmp<-pmax(tmp,tmp.n)  
#order(tmp)[1:(length(tmp)*0.99)]
int<-sample(int,size=10000)    
int<-sort(int)
seqinfo(int)<-si



# export as beds ----------------------------------------------------------


pos.crms$name<-paste0('pos',1:length(pos.crms))
export(pos.crms,con='analysis/make_regions_bedfiles/positive_8008.bed')

export(neg.crms,con='analysis/make_regions_bedfiles/negative_8008.bed')

export(int,con='analysis/make_regions_bedfiles/random.intergenic.bed')

intergenic<-import('analysis/make_regions_bedfiles/intergenic.inactive.regions.bed')
intergenic<-coverage(intergenic,width=seqlengths(si))>0






# Now let's try defining negative sets using our TF data ------------------


#get the rpgc values for our 

#get the tf data 
tffiles<-list.files('data/TFs/peaks/',full.names=T)
tffiles.tfnames<-gsub('peaks_(.*?)_.*','\\1',list.files('data/TFs/peaks/',full.names=F))
tffiles.tps<-gsub('.*?_(\\d+\\-\\d+).bed','\\1',tffiles)
tffile.df<-data.frame(tfname=tffiles.tfnames,tp=tffiles.tps,file=tffiles)
tfgrlist<-sapply(simplify=F,as.character(unique(tffile.df$tp)),function(tp){
  sapply(simplify=F,as.character(unique(tffile.df$tfname)),function(tfname){
    if(!any( tffile.df$tp==tp & tffile.df$tfname==tfname)){return(NULL)}
    
    gr<-import(as.character(tffile.df$file[tffile.df$tp==tp & tffile.df$tfname==tfname]),asRangedData=F)
    
  })
})
#now let's go through our 6-8 hour set and add the overlap info for the tfs at various timepoints

#for our tf,tp, do countoverlaps four our crm8008.gr
tmp<-sapply(simplify=F,as.character(unique(tffile.df$tp)),function(tp){
  sapply(simplify=F,as.character(unique(tffile.df$tfname)),function(tf){
    if(is.null(tfgrlist[[tp]][[tf]])){return(NULL)}   
    mcols(crm8008.gr)[[paste(tf,tp,sep='_')]]<<-countOverlaps(crm8008.gr,tfgrlist[[tp]][[tf]])>0 
    NULL
  })
}) 

#now lets classify our crms according to four groups
#All heart bound, all 5 mesobound,2 heart and 2 8008
#we now have columns in our gr matching the ts and tps
#get logical matrix specifying tfs for heart
heartmat<-as.matrix(mcols(crm8008.gr)[,c('tin_6.8','doc2_6.8','dTCF_6.8','mef2_6.8','pnr_6.8')])
mesomat<-as.matrix(mcols(crm8008.gr)[,c('tin_6.8','twi_6.8','bap_6.8','bin_6.8','mef2_6.8')])

mcols(crm8008.gr)$Activity24   <- (apply(as.matrix(mcols(crm8008.gr)[,colnames(mcols(crm8008.gr))[grepl("2.4", colnames(mcols(crm8008.gr)), fixed=T)]]) , 1, sum))
mcols(crm8008.gr)$Activity46   <- (apply(as.matrix(mcols(crm8008.gr)[,colnames(mcols(crm8008.gr))[grepl("4.6", colnames(mcols(crm8008.gr)), fixed=T)]]) , 1, sum))
mcols(crm8008.gr)$Activity68   <- (apply(as.matrix(mcols(crm8008.gr)[,colnames(mcols(crm8008.gr))[grepl("6.8", colnames(mcols(crm8008.gr)), fixed=T)]]) , 1, sum))
mcols(crm8008.gr)$Activity810   <- (apply(as.matrix(mcols(crm8008.gr)[,colnames(mcols(crm8008.gr))[grepl("8.10", colnames(mcols(crm8008.gr)), fixed=T)]]) , 1, sum))
mcols(crm8008.gr)$Activity1012   <- (apply(as.matrix(mcols(crm8008.gr)[,colnames(mcols(crm8008.gr))[grepl("10.12", colnames(mcols(crm8008.gr)), fixed=T)]]) , 1, sum))

crm8008.gr$Meso5<-apply(mesomat,1,all)
crm8008.gr$Heart5<-apply(heartmat,1,all)
crm8008.gr$Meso2<-apply(mesomat,1,function(x)sum(x)==2)
crm8008.gr$Heart2<-apply(heartmat,1,function(x)sum(x)==2)

#Now get the rpgc values for them all
alltags.rpgc$both<-alltags.rpgc$pos+alltags.rpgc$neg
crm8008.gr$allsum.rpgc<-unlist(viewSums(Views(alltags.rpgc$both,as(crm8008.gr,'RangesList'))))

#classify our crms using the tf data
#get the mean normalized rpgc figures for each crm in each class
tmp<-list(Heart_5TF=crm8008.gr$allsum.rpgc[ crm8008.gr$Heart5 & crm8008.gr$intergenic],
          Meso_5TF=crm8008.gr$allsum.rpgc[ crm8008.gr$Meso5 & crm8008.gr$intergenic],
          Heart_2TF=crm8008.gr$allsum.rpgc[ crm8008.gr$Meso2 & crm8008.gr$intergenic],
          Meso_2TF=crm8008.gr$allsum.rpgc[ crm8008.gr$Heart2 & crm8008.gr$intergenic])
tmp<-stack(tmp)
#do boxplot
jpeg('analysis/make_regions_bedfiles/Heart_Meso_5bound_2bound_TF_boxplot.jpeg')
qplot(y=tmp$values,log='y',x=tmp$ind,color=tmp$ind,geom='boxplot',main='RPGC Cage signal for 4 classes of CRMs',ylab='normalized CAGE signal',xlab='')
dev.off()
#show sizes of sets
table(tmp$ind)

crm8008.gr$in.tf.pos<- !crm8008.gr$H3K4me3_peak & crm8008.gr$intergenic & crm8008.gr$Activity68
crm8008.gr$in.tf.neg<- !crm8008.gr$H3K4me3_peak & crm8008.gr$intergenic & !crm8008.gr$Activity68 & !crm8008.gr$Activity46

pos.tf.crms<-crm8008.gr[ crm8008.gr$in.tf.pos]
neg.tf.crms<-crm8008.gr[ crm8008.gr$in.tf.neg]

export(pos.tf.crms,con ='analysis/make_regions_bedfiles/pos.tf.8008crms.bed')
export(pos.tf.crms,con='analysis/make_regions_bedfiles/neg.tf.8008crms.bed')

pos.tf.ac.crms<-crm8008.gr[ crm8008.gr$in.tf.pos & crm8008.gr$H3K27ac]
neg.tf.ac.crms<-crm8008.gr[  crm8008.gr$in.tf.neg & !crm8008.gr$H3K27ac]

export(pos.tf.ac.crms,con ='analysis/make_regions_bedfiles/pos.tf.ac.8008crms.bed')
export(pos.tf.ac.crms,con='analysis/make_regions_bedfiles/neg.tf.ac.8008crms.bed')


# Now use DNase as well. --------------------------------------------------
#get the summed dnase reads over 6-8hrs for our crms
crm.dnase.68<-Views(dnase.rles$STG10+dnase.rles$STG11,as(crm8008.gr,'RangesList'))
crm.dnase.68<-unlist(viewSums(crm.dnase.68))/width(crm8008.gr)

#does dnase sensitivity correlate with cage signal?
#scatter plot, coloring positives, negatives
#scatterplot showing dnase vs cage for intergenic and non
qplot(x=crm8008.gr$allsum.rpgc,y=crm.dnase.68,log='xy',xlab='RPGC cage signal',ylab='Mean DNase Sensitivity',
      color= crm8008.gr$intergenic ,main='Dnase Sensitivity vs Cage signal for 8008 crms' )+
  scale_color_discrete(name='Intergenic')

setvector=rep('Neither',length(crm8008.gr))
setvector[ crm8008.gr$in.tf.pos ]<-'Positive'
setvector[ crm8008.gr$in.tf.neg ]<-'Negative'
int<-crm8008.gr$intergenic

#scatterplot showing this for our positive negative tf set
qplot(x=crm8008.gr$allsum.rpgc[int],y=crm.dnase.68[int],log='xy',xlab='RPGC cage signal',ylab='Mean DNase Sensitivity',
      color= setvector[int] ,main='Dnase Sensitivity vs Cage signal for Intergenic 8008 crms' )+
  scale_color_discrete(name='Set - TF binding')


qplot(x=crm8008.gr$allsum.rpgc[int],color=setvector[int],geom='density',log='x',xlab='Total Normalized Cage Signal',
      main='Distribution of Cage Signal for Sets')+    scale_color_discrete(name='Set - TF binding')

setvector=rep('Neither',length(crm8008.gr))
setvector[ crm8008.gr$in.tf.pos & crm8008.gr$H3K27ac ]<-'Positive'
setvector[ crm8008.gr$in.tf.neg & ! crm8008.gr$H3K27ac ]<-'Negative'

#scatterplot showing this for our positive negative tf set
qplot(x=crm8008.gr$allsum.rpgc[int],y=crm.dnase.68[int],log='xy',xlab='RPGC cage signal',ylab='Mean DNase Sensitivity',
      color= setvector[int] ,main='Dnase Sensitivity vs Cage signal for Intergenic 8008 crms' )+
  scale_color_discrete(name='Set - TF binding, K27ac')

qplot(x=crm8008.gr$allsum.rpgc[int],color=setvector[int],geom='density',log='x',xlab='Total Normalized Cage Signal',
      main='Distribution of Cage Signal for Sets')+scale_color_discrete(name='Set - TF binding + K27ac')





#now define the positives and negatives with dnase as well.
crm8008.gr$low.dnase<-crm.dnase.68<quantile(crm.dnase.68,0.2)
crm8008.gr$high.dnase<-crm.dnase.68>quantile(crm.dnase.68,0.8)

setvector=rep('Neither',length(crm8008.gr))
setvector[ crm8008.gr$in.tf.pos & crm8008.gr$high.dnase ]<-'Positive'
setvector[ crm8008.gr$in.tf.neg & crm8008.gr$low.dnase ]<-'Negative'
int<-crm8008.gr$intergenic

#scatterplot showing this for our positive negative tf set, with dnase
qplot(x=crm8008.gr$allsum.rpgc[int],y=crm.dnase.68[int],log='xy',xlab='RPGC cage signal',ylab='Mean DNase Sensitivity',
      color= setvector[int] ,main='Dnase Sensitivity vs Cage signal for Intergenic 8008 crms' )+
  scale_color_discrete(name='Set - TF binding + Dnase')

qplot(x=crm8008.gr$allsum.rpgc[int],color=setvector[int],geom='density',log='x',xlab='Total Normalized Cage Signal',
      main='Distribution of Cage Signal for Sets')+  scale_color_discrete(name='Set - TF binding + Dnase')


##And finally do it with K27ac as well
setvector=rep('Neither',length(crm8008.gr))
setvector[ crm8008.gr$in.tf.pos & crm8008.gr$high.dnase & crm8008.gr$H3K27ac ]<-'Positive'
setvector[ crm8008.gr$in.tf.neg & crm8008.gr$low.dnase & !crm8008.gr$H3K27ac ]<-'Negative'
int<-crm8008.gr$intergenic


#scatterplot showing this for our positive negative tf set
qplot(x=crm8008.gr$allsum.rpgc[int],y=crm.dnase.68[int],log='xy',xlab='RPGC cage signal',ylab='Mean DNase Sensitivity',
      color= setvector[int] ,main='Dnase Sensitivity vs Cage signal for Intergenic 8008 crms' )+
  scale_color_discrete(name='Set - TF binding + Dnase + K27ac')

qplot(x=crm8008.gr$allsum.rpgc[int],color=setvector[int],geom='density',log='x',xlab='Total Normalized Cage Signal',
      main='Distribution of Cage Signal for Sets')+scale_color_discrete(name='Set - TF binding + Dnase + K27ac')



#now finally produce the most negative set by putting a cutoff on dnase - say the bottom 20%
#We can do an additional filter by filtering these again using the previous stages dnase


pos.tf.dnase.crms<-crm8008.gr[ crm8008.gr$in.tf.pos & crm8008.gr$high.dnase]
neg.tf.dnase.crms<-crm8008.gr[ crm8008.gr$in.tf.neg & crm8008.gr$low.dnase]

export(pos.tf.dnase.crms,con ='analysis/make_regions_bedfiles/pos.tf.8008crms.bed')
export(pos.tf.dnase.crms,con='analysis/make_regions_bedfiles/neg.tf.8008crms.bed')

pos.tf.ac.dnase.crms<-crm8008.gr[ crm8008.gr$in.tf.pos & crm8008.gr$H3K27ac & crm8008.gr$high.dnase ]
neg.tf.ac.dnase.crms<-crm8008.gr[  crm8008.gr$in.tf.neg & !crm8008.gr$H3K27ac & crm8008.gr$low.dnase ]

export(pos.tf.ac.crms,con ='analysis/make_regions_bedfiles/pos.tf.ac.8008crms.bed')
export(pos.tf.ac.crms,con='analysis/make_regions_bedfiles/neg.tf.ac.8008crms.bed')


save(crm8008.gr,file=file.crm8008.gr)


# Some summary graphs of the positive and negative sets -------------------



accs<-as.character(accession.df$accession[ !accession.df$reseq])
cagemat<-get.cage.matrix(crm8008.gr,sapply(cage.tag.rles[accs][],'[[','both'))
rownames(cagemat)<-crm8008.gr$id
tmp<-melt(cagemat)
tmp$Var2 <-factor(tmp$Var2, levels = as.character(unique(tmp$Var2)))
tmp$Var1 <-factor(tmp$Var1, levels = as.character(unique(tmp$Var1)))

#plot showing mean value for all crms across libraries
plot1<-ggplot(tmp,aes(x=as.factor(tmp$Var2),y=tmp$value))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",width=0.1,colour = "red")+
  scale_x_discrete(name='Library',labels='')+ggtitle('Mean Value for all CRMs, for each CAGE library')+scale_y_continuous(name='Mean value in unnormalized Tags')
plot1

#for normalized libraries
cagemat<-get.cage.matrix(crm8008.gr,sapply(cage.tag.rles.rpgc[accs][],'[[','both'))
rownames(cagemat)<-crm8008.gr$id
tmp<-melt(cagemat)
tmp$Var2 <-factor(tmp$Var2, levels = as.character(unique(tmp$Var2)))
tmp$Var1 <-factor(tmp$Var1, levels = as.character(unique(tmp$Var1)))

plot1<-ggplot(tmp,aes(x=as.factor(tmp$Var2),y=tmp$value))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",width=0.1,colour = "red")+
  scale_x_discrete(name='Library',labels='')+ggtitle('Mean Value for all CRMs, for each CAGE library')+scale_y_continuous(name='Mean value in Normalized Tags')
plot1


#now normalize cagemat by row so we can see if our libraries account for a similiar amount of each one
cagemat.rownorm<-t(apply(cagemat,1,function(x){x/median(x)}))
cagemat.rownorm[is.nan(cagemat.rownorm)]<-0
tmp<-melt(cagemat.rownorm)
tmp$Var2 <-factor(tmp$Var2, levels = as.character(unique(tmp$Var2)))
tmp$Var1 <-factor(tmp$Var1, levels = as.character(unique(tmp$Var1)))
plot1<-ggplot(tmp,aes(x=as.factor(tmp$Var2),y=tmp$value))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",width=0.1,colour = "red")+
  scale_x_discrete(name='crm')+stat_summary(fun.data="mean",geom="point",width=0.1,colour = "blue")
plot1

tmp<-melt(cagemat.rownorm)

#results of normalizing the 
plot1<-ggplot(tmp,aes(x=as.character(tmp$Var2),y=tmp$value))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",width=0.1,colour = "red")+
  scale_x_discrete(name='crm',breaks=as.character(unique(tmp$Var2)))





#Now for crms instead of libraries
setvector=rep('Neither',length(crm8008.gr))
setvector[ crm8008.gr$in.tf.pos & crm8008.gr$high.dnase & crm8008.gr$H3K27ac ]<-'Positive'
setvector[ crm8008.gr$in.tf.neg & crm8008.gr$low.dnase & !crm8008.gr$H3K27ac ]<-'Negative'

cagemat<-get.cage.matrix(crm8008.gr,sapply(cage.tag.rles.rpgc[accs][],'[[','both'))
rownames(cagemat)<-1:length(crm8008.gr)
tmp<-melt(cagemat)
tmp$set<-setvector[as.numeric(as.character(tmp$Var1))]

table(tmp$set)
table(setvector)


tmp$logval<-log10(tmp$value)
tmp$logval[tmp$logval==-Inf]<-0
tmp.s<-tmp[ tmp$set!='Neither',]
by(data=
tmp.s<-tmp.s[order(tmp.s$value) ,]
tmp.s<-tmp.s[order(tmp.s$set) ,]
#NOW reorder the factor levels
tmp.s$Var1 <-factor(tmp.s$Var1, levels = as.character(unique(tmp.s$Var1)))
tmp.s$Var2 <-factor(tmp.s$Var2, levels = as.character(unique(tmp.s$Var2)))

plot1<-ggplot(tmp.s,aes(x=factor(tmp.s$Var1),y=tmp.s$logval,color=tmp.s$set))+stat_summary(fun.data="mean_cl_boot",geom="errorbar")+
  scale_x_discrete(labels='',breaks=reorder(unique(tmp.s$Var2),setvector[unique(tmp.s$Var2)]))

jpeg('tmp.jpeg',h=2000,w=4000);print(plot1);dev.off()




jpeg('tmp.jpeg',h=2000,w=4000)
plot1<-ggplot(tmp,aes(x=factor(tmp$Var1),y=tmp$logval,color=tmp$set))+stat_summary(fun.data="mean_cl_boot",geom="errorbar")+scale_x_discrete()
dev.off()

sampcrms<-sample(unique(tmp.s$Var1),20)
tmp.s<-tmp.s[ tmp.s$Var1 %in% sampcrms ,]




















r=full.neg.gr[7]
b=start(r)
e=end(r)
width(r)

   c=as.character(seqnames(r))[1]
s='both'
#does alltag object match igv? Yes
plot(alltags.rpgc[[s]][[c]][b:e])
ch#do the sums match the alltags object?


#checking alltags object
sum(sapply(1:80,function(n){sum(cage.tag.rles[[n]]$pos[[as.character(seqnames(r))]][start(r):end(r)])}))
       

# calculate size of 'intergenic library' ----------------------------------
#get rle with the intergenic genome
transcripts.gr<-keepSeqlevels(transcripts.gr,chrs.keep)
int<-coverage(resize(transcripts.gr,width=width(transcripts.gr)+1000,fix='center'))
int<-int==0


#calculate size of intergenic library
accession.df$intergenic.lib.size <-sapply(accession.df$accession,function(acc){
  sum(sum(cage.tag.rles[[acc]]$neg[int]+cage.tag.rles[[acc]]$pos[int]))
})
tmp<-sapply(accession.df$accession,function(acc){
  sum(sum(cage.tag.rles[[acc]]$neg[!int]+cage.tag.rles[[acc]]$pos[!int]))
})
qplot(accession.df$library.size,accession.df$intergenic.lib.size,xlim=c(0,25000000),ylim=c(0,50000))
#calculate size of library at all TSS
tss<-coverage(resize(tss.gr,width=500,fix='center'))>0
tss<-tss[chrs.keep]
accession.df$tss.lib.size <-sapply(accession.df$accession,function(acc){
  sum(sum(cage.tag.rles[[acc]]$neg[tss]+cage.tag.rles[[acc]]$pos[tss]))
})
save(accession.df,file=file.accession.df)








