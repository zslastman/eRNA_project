
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/generate_Rle.R')
source('src/generate_chromdata_cage.R')
source('src/load_annotations.R')

# save.image('tmp.image.R')
# load('tmp.image.R')
#source('src/normalization.R')
library(reshape2)
library(ggplot2)
library(ROCR,lib.loc='~/Harnett/R/')
#define window sizes
window.sizes<-c('w25'=25,'w50'=50,'w100'=100,'w200'=200,'w300'=300,'w500'=500)
#exclude the reseq libraries.
cg<-cage.tag.rles[!grepl('reseq',names(cage.tag.rles))]
#windowsize
w=100
#precission cutoff when setting cutoff
prec.cutoff<-0.95













#####First define our inactive regions ---------------------------------------

#define 'inactive space'
act<-list(
  crmgrs,
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
#concatenate grs
act<-sapply(act,  function(gr){ mcols(gr) <-NULL;gr})
act<-do.call('c',act)
#resize them so inactive space is 1k from active
act<-resize(act,width(act)+1000,fix='center')
act<-reduce(act)#reduce to same strand nonoverlapping regions
act<-keepSeqlevels(act,chrs.keep)#only good chrs
strand(act)<-'*'
#define inactive space 
intergenic.inactive.regions<-gaps(act)
intergenic.inactive.regions<-intergenic.inactive.regions[strand(intergenic.inactive.regions)=='*']
#export tracks to make sure this all looks okay
export(transcripts.gr,con='analysis/make_regions_bedfiles/transcripts.gr.bed')
#export(gene.model.data,'analysis/make_regions_bedfiles/gene_models.bed')
export(act,con='analysis/make_regions_bedfiles/act.bed')
export(intergenic.inactive.regions,con='analysis/make_regions_bedfiles/intergenic.inactive.regions.bed')


#Now we filter out the unmappable or High GC content regions

#function that imports binned UCSC wig files and exports as an unbinned rle 
import.UCSC.binned<-function(trackfile){
  # Now use mappability and GC to filter out regions ------------------------
  trackfile<-import('~/Harnett/data_general/map_dm3_Bin10.wig.gz')#import as UCSC track
  trackfile<-trackfile[chrs.keep]#only good chrs
  seqinfo(trackfile)<-si
  #unbin
  trackfile<-coverage(x=trackfile,weight='score',width=sapply(cage.tag.rles[[1]][[1]],length))
  
}
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


int<-windowGRange(intergenic.inactive.regions,500)
int<-sample(int,size=20000)


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

#only mappable regions
int<-int[unlist(viewMeans(  Views(mappability.rle[chrs.keep],as(int,'RangesList')[chrs.keep])  ))==1]
#gc content
int<-int[abs(unlist(viewMeans(  (Views(gc.rle[chrs.keep],as(int,'RangesList')[chrs.keep])  )))-0.5)<0.2 ]
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
export(int,con='analysis/make_regions_bedfiles/random.intergenic.bed')














# Now CAD regions ---------------------------------------------------------


# CAD3 SET ----------------------------------------------------------------
cad3.gr<-sort(cad3.gr)
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













#3) Now 8008 CRMs -----------------------------------------------------------

# Chromatin data ----------------------------------------------------------


##mark intergenic CRMS
#also including lincRNAs
trancripts.lincs.gr<-do.call('c',sapply(list(lincRNA.gr,transcripts.gr),  function(gr){ mcols(gr) <-NULL;gr}))
#concatenate our transcripts and lincRNAs.
crmgrs$intergenic<-distanceToNearest(crmgrs,trancripts.lincs.gr)$distance>intergenic_dist

##And add the binary information from peaks
#mark those with any H3K4me3 overlap
crmgrs$H3K4me3_peak<- 0 < 
  (countOverlaps(crmgrs,chrompeaks.modencode[['K4me3_4-8h']])+ 
     countOverlaps(crmgrs,chrompeaks[['K4me3_6-8h']]))+countOverlaps(crmgrs,chrompeaks[['K4me3_4-6h']])
#mark those with H3K4me1 overlaps
crmgrs$H3K4me1_peak<- 0 < 
  (countOverlaps(crmgrs,chrompeaks.modencode[['K4me1_4-8h']])+  
     countOverlaps(crmgrs,chrompeaks[['K4me1_6-8h']]))+countOverlaps(crmgrs,chrompeaks[['K4me1_4-6h']])
#And add the K27ac info from our chromatin and the modencode projects'
crmgrs$H3K27ac<- 0 < 
  (countOverlaps(crmgrs,chrompeaks.modencode[['K27ac_4-8h']])+  
     countOverlaps(crmgrs,chrompeaks[['K27Ac_6-8h']]))+countOverlaps(crmgrs,chrompeaks[['K27Ac_4-6h']])
#And add the K36me3 info from our chromatin and the modencode projects'
crmgrs$H3K79me3<- 0 < 
#  (countOverlaps(crmgrs,chrompeaks.modencode[['K27ac_4-8h']])+  
     countOverlaps(crmgrs,chrompeaks[['K79me3_6-8h']])+countOverlaps(crmgrs,chrompeaks[['K79me3_4-6h']])
#And add the K79ac info from our chromatin and the modencode projects'
crmgrs$H3K36me3<- 0 < 
     countOverlaps(crmgrs,chrompeaks[['K36me3_6-8h']])+countOverlaps(crmgrs,chrompeaks[['K36me3_4-6h']])
#and our PolII
crmgrs$polII <-   0 < countOverlaps(crmgrs,chrompeaks[['PolII_6-8h']]) +countOverlaps(crmgrs,chrompeaks[['PolII_4-6h']])


##And continous information 
chrom.mean.mats<-sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crmgrs,'RangesList'))))})
chrom.mean.mats.mod<-sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crmgrs,'RangesList'))))})


# chrom.mean.mats2<-get.best.chr.window.mat(resize(crmgrs,width=750,fix='center'),w,cage=chrom.rles.rpgc.sub.merge)
# chrom.mean.mats2<-chrom.mean.mats2/w
# poschrom<- chrom.mean.mats2[crmgrs$abovecut & crmgrs$intergenic,]
# negchrom<-chrom.mean.mats2[!crmgrs$abovecut & crmgrs$intergenic,]
# print(Chrom.boxplots(poschrom , negchrom))

# comparing our datasets to the modencode datasets on the crmgrs ----------


cor(chrom.mean.mats[,'PolII_6.8'],chrom.mean.mats.mod[,'PolII'])
cor(chrom.mean.mats[,'H3K27ac_6.8'],chrom.mean.mats.mod[,'K27ac'])

qplot(x=chrom.mean.mats[,'PolII_6.8'],y=chrom.mean.mats.mod[,'PolII'])
qplot(x=chrom.mean.mats[,'H3K27ac_6.8'],y=chrom.mean.mats.mod[,'K27ac'])









# Now TF data -------------------------------------------------------------


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
#for our tf,tp, do countoverlaps four our crmgrs
tmp<-sapply(simplify=F,as.character(unique(tffile.df$tp)),function(tp){
  sapply(simplify=F,as.character(unique(tffile.df$tfname)),function(tf){
    if(is.null(tfgrlist[[tp]][[tf]])){return(NULL)}   
    mcols(crmgrs)[[paste(tf,tp,sep='_')]]<<-countOverlaps(crmgrs,tfgrlist[[tp]][[tf]])>0 
    NULL
  })
}) 

#now lets classify our crms according to four groups
#All heart bound, all 5 mesobound,2 heart and 2 8008
#we now have columns in our gr matching the ts and tps
#get logical matrix specifying tfs for heart
heartmat<-as.matrix(mcols(crmgrs)[,c('tin_6.8','doc2_6.8','dTCF_6.8','mef2_6.8','pnr_6.8')])
mesomat<-as.matrix(mcols(crmgrs)[,c('tin_6.8','twi_6.8','bap_6.8','bin_6.8','mef2_6.8')])

mcols(crmgrs)$Activity24   <- (apply(as.matrix(mcols(crmgrs)[,colnames(mcols(crmgrs))[grepl("2.4", colnames(mcols(crmgrs)), fixed=T)]]) , 1, sum))
mcols(crmgrs)$Activity46   <- (apply(as.matrix(mcols(crmgrs)[,colnames(mcols(crmgrs))[grepl("4.6", colnames(mcols(crmgrs)), fixed=T)]]) , 1, sum))
mcols(crmgrs)$Activity68   <- (apply(as.matrix(mcols(crmgrs)[,colnames(mcols(crmgrs))[grepl("6.8", colnames(mcols(crmgrs)), fixed=T)]]) , 1, sum))
mcols(crmgrs)$Activity810   <- (apply(as.matrix(mcols(crmgrs)[,colnames(mcols(crmgrs))[grepl("8.10", colnames(mcols(crmgrs)), fixed=T)]]) , 1, sum))
mcols(crmgrs)$Activity1012   <- (apply(as.matrix(mcols(crmgrs)[,colnames(mcols(crmgrs))[grepl("10.12", colnames(mcols(crmgrs)), fixed=T)]]) , 1, sum))

crmgrs$Meso5<-apply(mesomat,1,all)
crmgrs$Heart5<-apply(heartmat,1,all)
crmgrs$Meso2<-apply(mesomat,1,function(x)sum(x)==2)
crmgrs$Heart2<-apply(heartmat,1,function(x)sum(x)==2)

#Now get the rpgc values for them all
alltags.rpgc$both<-alltags.rpgc$pos+alltags.rpgc$neg
crmgrs$allsum.rpgc<-unlist(viewSums(Views(alltags.rpgc$both,as(crmgrs,'RangesList'))))

crmgrs$in.tf.pos <- !crmgrs$H3K4me3_peak & crmgrs$intergenic & crmgrs$Activity68
crmgrs$in.tf.neg<- !crmgrs$H3K4me3_peak & crmgrs$intergenic & !crmgrs$Activity68 & !crmgrs$Activity46

export(crmgrs[ crmgrs$in.tf.pos],con ='analysis/make_regions_bedfiles/pos.tf.8008crms.bed')
export(crmgrs[ crmgrs$in.tf.neg],con='analysis/make_regions_bedfiles/neg.tf.8008crms.bed')



# Plots comparing CAGE signal to PolII signal -----------------------------



qplot(crmgrs$allsum.rpgc,y={a=chrom.mean.mats.mod[,'PolII'];a[a<1]=1;a},log='xy',color=crmgrs$intergenic,
      xlab='Normalized CAGE tags',ylab='Modencode Whole Embryo PolII 4-8h',main='CAGE signal Vs. Modencode PolII, 8008 CRMs')
cor(crmgrs$allsum.rpgc,chrom.mean.mats.mod[,'PolII'],method='s')



qplot(crmgrs$allsum.rpgc,y={a=chrom.mean.mats.mod[,'PolII'];a[a<1]=1;a},log='xy',color=crmgrs$intergenic,
      xlab='Normalized CAGE tags',ylab='Meso PolII 6_8h ',main='CAGE signal Vs. Modencode PolII, 8008 CRMs')
cor(crmgrs$allsum.rpgc,chrom.mean.mats[,'PolII_6.8'],method='s')


qplot(crmgrs$allsum.rpgc[ crmgrs$intergenic],y={a=chrom.mean.mats.mod[,'PolII'][ crmgrs$intergenic];a[a<1]=1;a},log='xy',
      xlab='Normalized CAGE tags',ylab='Modencode Whole Embryo PolII 4-8h',main='CAGE signal Vs. Modencode PolII, Intergenic 8008 CRMs')
cor(crmgrs$allsum.rpgc[ crmgrs$intergenic],chrom.mean.mats.mod[,'PolII'][ crmgrs$intergenic],method='s')


qplot(crmgrs$allsum.rpgc[ crmgrs$intergenic],y={a=chrom.mean.mats.mod[,'PolII'][ crmgrs$intergenic];a[a<1]=1;a},log='xy',
      xlab='Normalized CAGE tags',ylab='Meso PolII 6_8h ',main='CAGE signal Vs. Modencode PolII, Integergenic 8008 CRMs')
cor(crmgrs$allsum.rpgc[ crmgrs$intergenic],chrom.mean.mats[,'PolII_6.8'][ crmgrs$intergenic],method='s')






# Now  DNase as well. --------------------------------------------------
#get the summed dnase reads over 6-8hrs for our crms
crm.dnase.68<-Views(dnase.rles$STG10+dnase.rles$STG11,as(crmgrs,'RangesList'))
crmgrs$dnase.density.68<-unlist(viewSums(crm.dnase.68))/width(crmgrs)
crmgrs$low.dnase<-crmgrs$dnase.density.68<quantile(crmgrs$dnase.density.68,0.2)

save(crmgrs,file=file.crmgrs)











# Now get cage information ------------------------------------------------
cadpos.nc<-sort(cadpos.nc)#sort the cads
cadpos.nc<-cad3.pos
mcols(cadpos.nc)<-NULL#and make a version with no mcols

cagecountmatlist<-list(#get the best window of w in each one.
  crm8008=get.best.window.mat(crmgrs,w,cage=cg),
  cadpos=get.best.window.mat(cadpos.nc,w,cage=cg)
)



####create the matrix of normalized counts.
v<-accession.df$genome.coverage
names(v)<-accession.df$accession

cagecountmatlist.r<-sapply(simplify=F,names(cagecountmatlist),function(reg){
    sapply(1:ncol(cagecountmatlist[[reg]]),function(n){
      cagecountmatlist[[reg]][,n]<-cagecountmatlist[[reg]][,n]/v[names(cg)[n]]
    })
})









#define our datasets, a list of matrices
setlist=list(
  tfset=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos,],
  tf.27ac.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K27ac,],
  tf.pol.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$polII,],
  tf.K79.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K79me3,],
 # tf.K36.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K36me3,],
  tf.K4me1.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K4me1_peak,],
  cad3pos=cagecountmatlist.r[['cadpos']]
)
#and the negative dataset we'll compare them to
crmgrs$in.full.neg <- crmgrs$in.tf.neg & !crmgrs$H3K27ac & crmgrs$low.dnase 
negmat=cagecountmatlist.r[['crm8008']][ crmgrs$in.full.neg ,]
nrow(negmat)

# full.neg.gr<-crmgrs[crmgrs$in.full.neg]
# full.neg.gr$name<-paste0('full.neg.',1:length(full.neg.gr))
# export(full.neg.gr,'analysis/crm_chrom_analysis/full.neg.bed')





# Define functions for plotting each set ----------------------------------



ROCfunc <- function (pos,neg,main,sub) {
  scores<-stack(list(pos=pos,neg=neg))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T,lwd=6)
  auc=slot(performance(sumpred,'auc'),'y.values')[[1]]
  aucstring=paste0('AUC = ',round(auc,4))
  text(x=0.90,y=0,labels=aucstring)
  
  posstring=paste0('Positive Set - ', length(pos))
  negstring=paste0('Negative Set - ', length(neg))
  text(x=0.85,y=0.2,labels=posstring)
  text(x=0.85,y=0.1,labels=negstring)
  
  title(main,sub)
  abline(coef=c(0,1),lty=3)
  grid(lwd=2)
}

PreRecfunc <- function (pos,neg,main,sub) {
  scores<-stack(list(pos=pos,neg=neg))
  sumpred<-prediction(scores$values,scores$ind)
  perf=performance(sumpred,'prec','rec')
  plot(perf,colorize=T,ylim=c(0.5,1))
  
  posstring=paste0('Positive Set - ', length(pos))
  negstring=paste0('Negative Set - ', length(neg))
  text(x=0.2,y=0.4,labels=posstring)
  text(x=0.2,y=0.35,labels=negstring)
  
  title(main,sub)
  cutoff<-max(perf@x.values[[1]][perf@y.values[[1]]>prec.cutoff],na.rm=T)
  
  grid(lwd=2)
  abline(h=prec.cutoff,lty=3)
  abline(v=cutoff,lty=3)
  
  perf<-performance(sumpred,'prec')
  cutoff<-min(perf@x.values[[1]][perf@y.values[[1]]>prec.cutoff],na.rm=T)
  cutstring=paste0('Cutoff at precision ',prec.cutoff,' - ',round(cutoff,1),' RPGC tags')
  text(x=0.4,y=0.55,labels=cutstring)
  return(cutoff)
}

Barplot.from.matrices<-function(posmat,negmat, yname='Normalized Cage Signal 95% cl', xname='CRM', barsize=2,cut=1,  bar=T  ){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Positive',nrow(posmelt)),rep('Negative',nrow(negmelt)))
  
  
  #NOW reorder the factor levels
  bothmelt$Var1 <-factor(bothmelt$Var1, levels = as.character(unique(bothmelt$Var1)))
  bothmelt$Var2 <-factor(bothmelt$Var2, levels = as.character(unique(bothmelt$Var2)))
  bothmelt$set <-factor(bothmelt$set, levels = as.character(unique(bothmelt$set)))
  
  means=by(bothmelt,bothmelt$Var1,FUN=function(df){mean(df$value)})
  names(means)<-unique(bothmelt$Var1)
  
  bothmelt<-bothmelt[order(means[as.character(bothmelt$Var1)]) ,]
  bothmelt<-bothmelt[order(bothmelt$set) ,]
  #NOW reorder the factor levels
  bothmelt$Var1 <-factor(bothmelt$Var1, levels = as.character(unique(bothmelt$Var1)))
  bothmelt$Var2 <-factor(bothmelt$Var2, levels = as.character(unique(bothmelt$Var2)))
  bothmelt$set <-factor(bothmelt$set, levels = as.character(unique(bothmelt$set)))
  bothmelt$logval<-bothmelt$value
  bothmelt$logval[ bothmelt$logval==0]<-1#note that we take the logarithm of the logvals with ggplot later
  
  
  if(bar==T){
    return(
  #line graph of our se
  ggplot(data=bothmelt,aes(x=Var1,y=logval,color=set))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=barsize)+
    scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000))+
    #scale_y_continuous(name=yname)+
    scale_x_discrete(name=xname,labels='')+
    ggtitle('95% cl for mean of RPGC normalized Signal')
    )
  }
  if(bar==F){
    return(
    ggplot(data=bothmelt,aes(x=Var1,y=value,color=set))+stat_summary(fun.data="mean_cl_boot",geom="pointrange",size=1)+
      scale_y_log10(name=yname)+
      #scale_y_continuous(name=yname)+
      scale_x_discrete(name=xname,labels='')
    )
  }  
}
#Barplot.from.matrices(posmat,negmat,bar=T)
Chrom.boxplots<-function(posmat,negmat, lowlim=F,yname='Mean Inp-Sub Chromatin Signals for crms above/below cutoff ', xname='Set'){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Above',nrow(posmelt)),rep('Below',nrow(negmelt)))
  #we need to do something iwth our zero/negative values to do a log scale....
  if(lowlim==F){
    lowlim=min(bothmelt$value[bothmelt$value>0])
  }
  bothmelt$value[bothmelt$value<=0] <- lowlim
  

    #line graph of our se
    ggplot(data=bothmelt,aes(x=set,facet=Var2,y=value,color=Var2))+geom_boxplot()+
      scale_y_log10(name=yname)+    
      facet_wrap(~Var2)+
      #scale_y_continuous(name=yname)+
      scale_x_discrete(name=xname,labels=c('Above','Below'))     
}
# chrom.mean.mats[crmgrs$abovecut,] ->posmat
# chrom.mean.mats[!crmgrs$abovecut,]->negmat
# Chrom.boxplots( chrom.mean.mats[crmgrs$abovecut,] ,chrom.mean.mats[!crmgrs$abovecut,],lowlim=0.000000001)
# Chrom.boxplots( chrom.mean.mats[crmgrs$abovecut,] ,chrom.mean.mats[!crmgrs$abovecut,],)

#for now we'll just score the matrices by adding up the columns (libraries)
scorefunc<-rowSums


setn=names(setlist)[6]


# Now go through each of our sets and produce our plots -------------------

for(setn in names(setlist)){
  #dir.create('analysis/crm_chrom_analysis/')
  
  posmat=setlist[[setn]]#get our data
  negmat=negmat
  pos=scorefunc(posmat)#for now just sum across lines
  neg=scorefunc(negmat)
 
  pdf(paste0('analysis/crm_chrom_analysis/crmbarplot',setn,'.pdf'))
  print(Barplot.from.matrices(posmat,negmat,bar=T))
  dev.off()
  
  pdf(paste0('analysis/crm_chrom_analysis/crmrocplot',setn,'.pdf'))
  ROCfunc(pos,neg,
          main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,setn),
          sub=paste0('windowsize ',w))
  dev.off()
  
  pdf(paste0('analysis/crm_chrom_analysis/crm_prec_rec_plot',setn,'.pdf'))
  cutoff<-PreRecfunc(pos,neg,
          main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,setn),
          sub=paste0('windowsize ',w))
  dev.off()
  
  #now define our above and below cutoff sets.
  crmgrs$abovecut<-scorefunc(cagecountmatlist.r[['crm8008']])>cutoff
  
  #and produce the boxplots for the other chromatin marks
  pdf(paste0('analysis/crm_chrom_analysis/crm_chrom_boxplot',setn,'.pdf'))
  
  poschrom<- chrom.mean.mats[crmgrs$abovecut & crmgrs$intergenic,]
  negchrom<-chrom.mean.mats[!crmgrs$abovecut & crmgrs$intergenic,]
  print(Chrom.boxplots(poschrom , negchrom))
  
  dev.off()
  # poschrom<- chrom.mean.mats2[crmgrs$abovecut & crmgrs$intergenic,]
  # negchrom<-chrom.mean.mats2[!crmgrs$abovecut & crmgrs$intergenic,]
  # print(Chrom.boxplots(poschrom , negchrom))  
  
  cat(setn)
  cat(' ... ')
  
}  
  
# exporting tracks for IGBY viewing ---------------------------------------

#gr with our negatives
full.neg.gr<-crmgrs[crmgrs$in.full.neg]
#order them by their cage signal
neg=scorefunc(negmat)
full.neg.gr<-full.neg.gr[order(neg,decreasing=T)]
#appropriate names
full.neg.gr$name<-paste0('full.neg.',1:length(full.neg.gr))
mcols(full.neg.gr)<-mcols(full.neg.gr)[,c('name')]
full.neg.gr$name<-full.neg.gr$value
export(full.neg.gr,'analysis/crm_chrom_analysis/full.neg.bed')


#gr with our negatives
above.cutoff.gr<-crmgrs[ crmgrs$abovecut & crmgrs$intergenic]
#order them by their cage signal
above.cutoff.gr<-above.cutoff.gr[order(above.cutoff.gr$allsum.rpgc,decreasing=T)]
#appropriate names
above.cutoff.gr$name<-paste0('above.cutoff.',1:length(above.cutoff.gr))
mcols(above.cutoff.gr)<-mcols(above.cutoff.gr)[,c('name')]
above.cutoff.gr$name<-above.cutoff.gr$value
export(above.cutoff.gr,'analysis/crm_chrom_analysis/above.cutoff.bed')


make.IGV.snapshot.batchfile <- function(snapranges,snapshotdir,igvfile="tmp.igvbatch.txt",windowsize=20000, wait.milliseconds=5000){
  dir.create(snapshotdir, showWarnings = F)
  write(igvfile,x="")
  write(igvfile,x=paste0("snapshotDirectory ",snapshotdir),append=T)
  for(i in 1:length(snapranges)){
    chr=snapranges@seqnames[i]
    start=start(snapranges)[i]-(windowsize/2)
    end=end(snapranges)[i]+(windowsize/2)
    write(igvfile,x=paste0('goto ',chr,':',start,'-',end),append=T)
    write(igvfile,x=paste0('snapshot'),append=T)
    write(igvfile,x=paste0("setSleepInterval ",wait.milliseconds),append=T)        
  }
}



make.IGV.snapshot.batchfile(snapranges=full.neg.gr[1:20],snapshotdir='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/browser_screenshots/',igvfile='negbatch.tmp.igvbatch.txt',windowsize=10000)
make.IGV.snapshot.batchfile(snapranges=above.cutoff.gr[1:20],snapshotdir='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/browser_screenshots/',igvfile='abovecut.tmp.igvbatch.txt',windowsize=10000)

Views( alltags above.cutoff.gr[10])









