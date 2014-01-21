
# source('src/generate_Rle.R')
# source('src/load_annotations.R')

 setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')

source('src/tss_cage_functions.R')
#load(file.cage.tag.rles)
load(file.cg.pl)
load(file.crm8008.gr)
source('src/generate_chromdata_cage.R')
load(file.transcripts)
load(file.cad3)
load(file.tss)
load(file.alltags)
# save.image('tmp.image.R')
# load('tmp.image.R')
#source('src/normalization.R')
library(reshape2)
library(ggplot2)
library(ROCR,lib.loc='~/Harnett/R/')
library(LSD)

#cg<-cage.tag.rles
cg<-cg.pl
#windowsize
w=100
#precission cutoff when setting cutoff
prec.cutoff<-0.95


# Now CAD regions ---------------------------------------------------------


# CAD3 SET ----------------------------------------------------------------
cad3.gr$polII = countOverlaps(cad3.gr,chrompeaks[['PolII_6-8h']])>0
cad3.gr$pos <- cad3.gr$active68 %in% T & cad3.gr$intergenic 
cad3.gr$neg <- cad3.gr$inactive68 %in% T & cad3.gr$intergenic & !cad3.gr$GJ & !cad3.gr$polII

save(cad3.gr,file=file.cad3)


center.gr<-function(rle,gr){
  centred.gr<-nomcols(gr)
  cageviews<-Views(alltags$both,as(centred.gr,'RangesList'))[unique(seqnames(centred.gr))] 
  maxpos<- unlist( viewWhichMaxs(cageviews,na.rm=T) )
  start(centred.gr)<-maxpos-250
  centred.gr<-resize(centred.gr,width=500,fix='start')
  centred.gr<-sort(centred.gr)
}

cad.centered.gr<-center.gr(alltags$both,cad3.gr)
    #we need to make versions of our cad3 range centred on the cage signal
  # cad.centered.gr<-nomcols(cad3.gr)
  # cageviews<-Views(alltags$both,as(cad.centered.gr,'RangesList'))[unique(seqnames(cad.centered.gr))] 
  # maxpos<- unlist( viewWhichMaxs(cageviews,na.rm=T) )
  # start(cad.centered.gr)<-maxpos-250
  # cad.centered.gr<-resize(cad.centered.gr,width=500,fix='start')
  # cad.centered.gr<-sort(cad.centered.gr)

#export oour positives and negatives for viewing
export(keepSeqlevels(cad3.gr[ cad3.gr$pos ],chrs.keep),'analysis/cad3.pos.bed')
export(keepSeqlevels(cad3.gr[ cad3.gr$neg ],chrs.keep),'analysis/cad3.neg.bed')








#3) Now 8008 CRMs -----------------------------------------------------------

# Chromatin data ---------------------------------------------------------

overlapsAny<-function(gr,peaklist,maxgap=50){
  0 < rowSums(sapply(peaklist,function(peaks){
    countOverlaps(gr,peaks,maxgap)
  }))
}

#columns describind overlap with peak
crm8008.gr$H3K4me3_peak   <-overlapsAny(crm8008.gr,list(chrompeaks.modencode[['K4me3_4-8h']], chrompeaks[['K4me3_4-6h']] ,chrompeaks[['K4me3_6-8h']]),maxgap=50) 
crm8008.gr$H3K4me1_peak   <-overlapsAny(crm8008.gr,list(chrompeaks.modencode[['K4me1_4-8h']], chrompeaks[['K4me1_4-6h']] ,chrompeaks[['K4me1_6-8h']]),maxgap=50)
crm8008.gr$H3K27ac_peak   <-overlapsAny(crm8008.gr,list(chrompeaks.modencode[['K27ac_4-8h']], chrompeaks[['K27Ac_4-6h']] ,chrompeaks[['K27Ac_6-8h']]),maxgap=50)
crm8008.gr$H379me3_peak   <-overlapsAny(crm8008.gr,list(chrompeaks[['K79me3_4-6h']],chrompeaks[['K79me3_6-8h']]),maxgap=50)
crm8008.gr$H3K36me3_peak  <-overlapsAny(crm8008.gr,list(chrompeaks[['K36me3_4-6h']],chrompeaks[['K36me3_6-8h']]),maxgap=50)
crm8008.gr$polII          <-overlapsAny(crm8008.gr,list(chrompeaks[['PolII_4-6h']],chrompeaks[['PolII_6-8h']]),maxgap=50)

##And continous information for graphing later


# DNase as well. --------------------------------------------------
#get the summed dnase reads over 6-8hrs for our crms
crm.dnase.68<-Views(dnase.rles$STG10+dnase.rles$STG11,as(crm8008.gr,'RangesList'))
crm8008.gr$dnase.density.68<-unlist(viewSums(crm.dnase.68))/width(crm8008.gr)
crm8008.gr$low.dnase<-crm8008.gr$dnase.density.68<quantile(crm8008.gr$dnase.density.68,0.2)


chrom.rles.rpgc.sub.merge=chrom.rles.rpgc.sub.merge[!grepl('K36me3|K4me3',names(chrom.rles.rpgc.sub.merge))]
chrom.rles.rpgc.sub.merge=c(chrom.rles.rpgc.sub.merge,'dnase_6.8'=(dnase.rles$STG10+dnase.rles$STG11))


chrom.mean.mats<-list(
  crm8008=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crm8008.gr,'RangesList'))))}),
  cad=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
)
save(chrom.mean.mats,file='chrom.mean.mats.object.R')
load('chrom.mean.mats.object.R')
chrom.mean.mats.modencode<-list(
  crm8008=sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crm8008.gr,'RangesList'))))}),
  cad=sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
)


# chrom.mean.mats.mod<-sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crm8008.gr,'RangesList'))))})


# chrom.mean.mats2<-get.best.chr.window.mat(resize(crm8008.gr,width=750,fix='center'),w,cage=chrom.rles.rpgc.sub.merge)
# chrom.mean.mats2<-chrom.mean.mats2/w
# poschrom<- chrom.mean.mats2[crm8008.gr$abovecut & crm8008.gr$intergenic,]
# negchrom<-chrom.mean.mats2[!crm8008.gr$abovecut & crm8008.gr$intergenic,]
# print(Chrom.boxplots(poschrom , negchrom))

# comparing our datasets to the modencode datasets on the crm8008.gr ----------


# cor(chrom.mean.mats[,'PolII_6.8'],chrom.mean.mats.mod[,'PolII'])
# cor(chrom.mean.mats[,'H3K27ac_6.8'],chrom.mean.mats.mod[,'K27ac'])

# qplot(x=chrom.mean.mats[,'PolII_6.8'],y=chrom.mean.mats.mod[,'PolII'])
# qplot(x=chrom.mean.mats[,'H3K27ac_6.8'],y=chrom.mean.mats.mod[,'K27ac'])









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

#finally construct our negative and positve sets from TF data
crm8008.gr$in.tf.pos <- !crm8008.gr$H3K4me3_peak & crm8008.gr$intergenic & crm8008.gr$Activity68
crm8008.gr$in.tf.neg<- !crm8008.gr$H3K4me3_peak & crm8008.gr$intergenic & !crm8008.gr$Activity68 & !crm8008.gr$Activity46

export(crm8008.gr[ crm8008.gr$in.tf.pos],con ='analysis/make_regions_bedfiles/pos.tf.8008crms.bed')
export(crm8008.gr[ crm8008.gr$in.tf.neg],con='analysis/make_regions_bedfiles/neg.tf.8008crms.bed')


save(crm8008.gr,file=file.crm8008.gr)

cagecountmatlist<-list(#get the best window of w in each one.
  crm8008=getBestWindowMat(crm8008.gr,w,cage=cg),
  cad=getBestWindowMat(  nomcols(cad3.gr)  ,w,cage=cg)
#  tss=getBestWindowMat(  nomcols(tss.gr)  ,w,cage=cg)
  # cadpos=getBestWindowMat(  nomcols(cad3.gr[cad3.gr$pos])  ,w,cage=cg),
  # cadneg=getBestWindowMat(  nomcols(cad3.gr[ cad3.gr$neg ])  ,w,cage=cg)
)

load(file.tss)
tss.gr[!duplicated(tss.gr)]
tss.gr<-tss.gr[seqnames(tss.gr)%in% bigchrs]
# Now get cage information ------------------------------------------------

tss.gr=sort(tss.gr)
tssmat<-getBestWindowMat(  nomcols(tss.gr)  ,w,cage=cg)
save(tssmat,file='tssmat.object.R')
save(cagecountmatlist,file='cagecountmatlist.object.R')

load('tssmat.object.R')
load('cagecountmatlist.object.R')


cagecountmatlist[[3]]<-tssmat
names(cagecountmatlist)[[3]]<-'tss'


#######################################use all this data to define our datasets
#define our datasets, a list of matrices
possetlist=list(
  tfset=cagecountmatlist[['crm8008']][crm8008.gr$in.tf.pos,],
  tf.27ac.set=cagecountmatlist[['crm8008']][crm8008.gr$in.tf.pos & crm8008.gr$H3K27ac,],
  tf.pol.set=cagecountmatlist[['crm8008']][crm8008.gr$in.tf.pos & crm8008.gr$polII,],
  tf.K79.set=cagecountmatlist[['crm8008']][crm8008.gr$in.tf.pos & crm8008.gr$H3K79me3,],
  # tf.K36.set=cagecountmatlist[['crm8008']][crm8008.gr$in.tf.pos & crm8008.gr$H3K36me3,],
  tf.K4me1.set=cagecountmatlist[['crm8008']][crm8008.gr$in.tf.pos & crm8008.gr$H3K4me1_peak,],
  cad3pos=cagecountmatlist[['cad']][cad3.gr$pos,]
)

negsetlist=list(
  full.neg= cagecountmatlist[['crm8008']][crm8008.gr$in.tf.neg & !crm8008.gr$H3K27ac & crm8008.gr$low.dnase & ! crm8008.gr$polII,],
  cad.neg=  cagecountmatlist[['cad']][cad3.gr$neg,]
)




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
  cutstring=paste0('Cutoff at precision ',prec.cutoff,' - ',round(cutoff,4),' Normalized Tags')
  text(x=0.4,y=0.55,labels=cutstring)
  return(cutoff)
}
# list(pos=c(rep(T,30),rep(F,9)),neg=c(rep(F,100),rep(T,80),rep(F,7)))->above.cut.vects

Barplot.from.matrices<-function(posmat,negmat, logy=T,yname='Normalized Cage Signal 95% cl', xname='CRM',above.cut.vects=NA, barsize=2,cut=1,  bar=T,reorder=F  ){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Positive',nrow(posmelt)),rep('Negative',nrow(negmelt)))
  if(length(above.cut.vects)!=1){

    stopifnot(is(above.cut.vects,'list'))
    stopifnot(is(above.cut.vects,'vector'))

    bothmelt$abovecut<-c(above.cut.vects[[1]],above.cut.vects[[2]])

  }



  if(reorder==T){
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
  }
  
  bothmelt$logval <-bothmelt$value
  if(logy ){
    # m<-min(bothmelt$logval[ bothmelt$logval != 0 ])
    # bothmelt$logval <- bothmelt$logval / m 
    bothmelt$logval[ bothmelt$logval==0]<-1

  }

  #note that we take the logarithm of the logvals with ggplot later

  yscale=ifelse(logy,scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000)),scale_y_continuous(name=yname))

  if(length(above.cut.vects)!=1){
    return(
      #line graph of our se
      ggplot(data=bothmelt,aes(x=Var1,y=logval,color=set))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=barsize)
      +
        yscale+
        #scale_y_continuous(name=yname)+
        scale_x_discrete(name=xname,labels='')+
        ggtitle('95% cl for mean of  normalized Signal')
    )
  }else{
    return(
       #line graph of our se
      ggplot(data=bothmelt,aes(x=Var1,y=logval,color=set))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=barsize)+
        scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000))+
        #scale_y_continuous(name=yname)+
        scale_x_discrete(name=xname,labels='')+
        ggtitle('95% cl for mean of normalized Signal')
    )
  }  
}

# Barplot.from.matrices(posmat,negmat,logy=T)

# Barplot.from.matrices(posmat[1:10,],negmat[1:10,],bar=T)
Chrom.boxplots<-function(posmat,negmat, lowlim=F,yname='Mean Inp-Sub Chromatin Signals for crms above/below cutoff ', xname='Set',tit=''){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Above',nrow(posmelt)),rep('Below',nrow(negmelt)))
  bothmelt$nset<-c(rep(nrow(posmat),nrow(posmelt)),rep(nrow(negmat),nrow(negmelt)))

  #we need to do something iwth our zero/negative values to do a log scale....
  if(lowlim==F){
    lowlim=min(bothmelt$value[bothmelt$value>0])
  }
  bothmelt$value[bothmelt$value<=0] <- lowlim
  
  num.df <-  data.frame(set=c('Above','Below'),number=c(nrow(posmat),nrow(negmat)))
  num.df <- do.call(rbind, lapply(unique(bothmelt$Var2),function(v) cbind(num.df,v) ))
  num.df$Var2<-num.df$v
  num.df$yval<-min(bothmelt$value)


  #line graph of our se
  ggplot(data=bothmelt,aes(x=set,facet=Var2,y=value,color=Var2))+geom_boxplot()+
    scale_y_log10(name=yname)+    
    facet_wrap(~Var2)+
    #scale_y_continuous(name=yname)+
    scale_x_discrete(name=xname,labels=c('Above','Below'))+     
    geom_text(data=num.df ,aes(x=set,y=yval,label=number))+
    ggtitle(tit)+
    theme(strip.text=element_text(size=10))
}



dotplot.from.matrices<-function(posmat,negmat,yname='Normalized Cage Signal', xname='CRM'){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)  
  bothmelt$CRM<-bothmelt$Var1
  bothmelt$set<-c(rep('Positive',nrow(posmelt)),rep('Negative',nrow(negmelt)))
  
      #line graph  our se
      ggplot(data=bothmelt,aes(x=CRM,y=value,color=set))+geom_point()+
      scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000))
       
}
dotplot.from.matrices(posmat,negmat)



simple.barplot.matrices<-function(posmat,negmat,yname='Normalized Cage Signal', xname='CRM',main.tit='tit'){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Positive',nrow(posmelt)),rep('Negative',nrow(negmelt)))
  bothmelt$CRM<-bothmelt$Var1
  bothmelt[['Normalized_Signal']]<-bothmelt$value

  
      #line graph  our se
      ggplot(data=bothmelt,aes(x=CRM,y=Normalized_Signal,color=set))+
      stat_summary(fun.data="mean_cl_boot",geom="errorbar")+
      ggtitle(main.tit)

      # scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000))
       
}
simple.barplot.matrices(posmat,negmat,main.tit='Mean Cage Signal 95% bootstrap cl - CAD enhancers')


number.not.zero<-function(mat){apply(mat,1,function(x)sum(x!=0))}

#for now we'll just score the matrices by adding up the columns (libraries)
scorefunc<-function(x){rowSums(x)}
#setn=names(possetlist)[6]







#what if we have a list of lists, the bottom lists are going to contain
#the type of scoring used, the cutoff, vectors describing pos and neg for cad and crms
#the descriminances and accuracy.
negn=names(negsetlist)[2]
setn=names(possetlist)[6]
posmat=possetlist[[setn]]#get our data
negmat=negsetlist[[negn]]


dotplot.from.matrices(posmat,negmat)
simple.barplot.matrices(posmat,negmat)





cutoffs<-sapply(negsetlist,function(x)sapply(possetlist,function(z)list()))
ncutoffs<-sapply(negsetlist,function(x)sapply(possetlist,function(z)list()))

# Now go through each of our sets and produce our plots -------------------
for(negn in names(negsetlist)){
  for(setn in names(possetlist)){
    #dir.create('analysis/crm_chrom_analysis2/')
    
    posmat=possetlist[[setn]]#get our data
    negmat=negsetlist[[negn]]
    pos=scorefunc(posmat)#for now just sum across lines
    neg=scorefunc(negmat)
    
    pdf(paste0('analysis/crm_chrom_analysis2/crmbarplot',setn,negn,'.pdf'))
    print(Barplot.from.matrices(posmat,negmat,bar=T))
    dev.off()
    
    pdf(paste0('analysis/crm_chrom_analysis2/crmrocplot',setn,negn,'.pdf'))
    ROCfunc(pos,neg,
            main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,setn),
            sub=paste0('windowsize ',w))
    dev.off()
    
    pdf(paste0('analysis/crm_chrom_analysis2/crm_prec_rec_plot',setn,negn,'.pdf'))
    cutoff<-PreRecfunc(pos,neg,
                       main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,setn),
                       sub=paste0('windowsize ',w))
    cutoff
    dev.off()
    
    cutoffs[[negn]][[setn]][['crm8008']]  <-scorefunc(cagecountmatlist[[1]])>cutoff
    cutoffs[[negn]][[setn]][['cad3']]     <-scorefunc(cagecountmatlist[[2]])>cutoff

    #and produce the boxplots for the other chromatin marks
    pdf(paste0('analysis/crm_chrom_analysis2/crm_chrom_boxplot',setn,negn,'.pdf'))
    poschrom<- chrom.mean.mats[[1]][cutoffs[[negn]][[setn]][['crm8008']] & crm8008.gr$intergenic,]
    negchrom<-chrom.mean.mats[[1]][!cutoffs[[negn]][[setn]][['crm8008']] & crm8008.gr$intergenic,]
    print(Chrom.boxplots(poschrom , negchrom,tit=paste0('All Intergenic CRMs, Split by cutoff ',cutoff)))
    dev.off()
    
    cat(setn)
    cat(negn)
    cat(' ... ')
  }  
}





#Now for cads and crms
#for cads we use the 500 bp around cage center
#rank by cage signal
#and then produce heatmaps for various marks
# chrom.rles.rpgc.sub.merge->chrom.rles.rpgc.sub.merge.bak
# crm8008.gr->crm8008.gr.bak


profile.matrix<-function(rle,gr){
  stopifnot(length(unique(width(gr)))==1)
  views<-GRViews(rle,gr)#get the views
  matlist<-lapply(views,function(x){#turn them into a matrix
      suppressWarnings(as.matrix(x))
  })
  do.call(rbind,matlist)
}




rankedheatmap<-function(x,col1='white',col2='black',brks='quantile',filt=NA,rankvect=NA){
    if(!is.na(rankvect[[1]])){x=x[rankvect,]}
   
    if(!is.na(filt)){z=x[x>filt]
    }else{z<-x}

    if(brks=='quantile'){brks= quantile(z, probs = seq(0, 1, 0.01))}
    else if(brks=='uniform'){brks= seq(min(z),max(z),length.out=100)
    }else(warn('brks should be quantile, uniform'))

    image(t(x[nrow(x):1, ]), 
          breaks =brks,
          col = colorRampPalette(c(col1, col2))(length(brks ) - 1L),
          xaxt = 'none',
          yaxt = 'none')

}





load(file.alltags.pl)

chrom.profile.mats<-list(
  crm8008=lapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){
    profile.matrix(chrom.rle,center.gr(alltags.pl$both,resize(crm8008.gr,wi=500,fi='center')))
  }),
  cad=lapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){
    profile.matrix(chrom.rle,center.gr(alltags.pl$both,cad.centered.gr))
  })
 )











#go through the different marks and make heatmaps with our enhancers ranked by cage signal
for(setname in names(chrom.profile.mats)){
  jpeg(w=1920,h=1200,paste0('analysis/crm_chrom_analysis2/ranked_heatmap',setname,'.jpeg'))
  par(mfrow=c(1,6))

  for(markname in names(chrom.profile.mats[[setname]]) ){


    chrdata=chrom.profile.mats[[setname]][[markname]]
    
    cagerank=order(scorefunc(cagecountmatlist[[setname]]),decreasing=T)

    chrcol='red'

    rankedheatmap(chrdata,col2=chrcol,brks='quantile',filt=0,rankvect=cagerank)
    
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ',setname,'\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)
   

  }
 dev.off()
}







#using our final annotations of above and below cutoff, how these divide up the positive sets

#boxplots for the tf binding set should be our abovecut guys, vs. guys with lower cage signal
cutoff=113
rs=rowSums(cagecountmatlist[[1]])
crm8008.gr$abovecut<-rs>cutoff
cad3.gr$abovecut<-rowSums(cagecountmatlist[[2]])>cutoff
crm8008.gr$wellbelowcut<-rs<quantile(rs[rs<cutoff],0.25)#those rows which are in the bottom half of those below the cutoff
rsc<-rowSums(cagecountmatlist[[2]])
cad3.gr$wellbelowcut<-rsc<quantile(rsc[rsc<cutoff],0.25)





poschrom<-chrom.mean.mats[[1]][crm8008.gr$in.tf.pos & crm8008.gr$abovecut,]
negchrom<-chrom.mean.mats[[1]][crm8008.gr$in.tf.pos & !crm8008.gr$abovecut,] 
vnegchrom<-chrom.mean.mats[[1]][crm8008.gr$in.tf.pos & !crm8008.gr$abovecut & crm8008.gr$wellbelowcut,] 
posrank=order(scorefunc(cagecountmatlist[[1]][crm8008.gr$in.tf.pos & crm8008.gr$abovecut,]),decreasing=T)
negrank=order(scorefunc(cagecountmatlist[[1]][crm8008.gr$in.tf.pos & !crm8008.gr$abovecut,]),decreasing=T)

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_tf_boxplots','.jpeg'))

print(Chrom.boxplots(poschrom , negchrom))

dev.off()

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_tf_boxplots_vneg','.jpeg'))

print(Chrom.boxplots(poschrom , vnegchrom))

dev.off()

for(markname in names(chrom.profile.mats[[1]]) ){
    poschrom<-chrom.profile.mats[[1]][[markname]][crm8008.gr$in.tf.pos & crm8008.gr$abovecut,]
    negchrom<-chrom.profile.mats[[1]][[markname]][crm8008.gr$in.tf.pos & !crm8008.gr$abovecut,] 

    jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_tf_heatmaps',markname,'.jpeg'))
    par(mfrow=c(1,2))

 
    rankedheatmap(poschrom,col2=chrcol,brks='quantile',filt=0,rankvect=posrank)
    
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ','tf_bound above cutoff','\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)
   

    rankedheatmap(negchrom,col2=chrcol,brks='quantile',filt=0,rankvect=negrank)
 
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ','tf_bound below cutoff','\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)


    dev.off()

}



#boxplots for the cad binding set




poschrom<-chrom.mean.mats[[2]][ cad3.gr$pos & cad3.gr$abovecut,]
negchrom<-chrom.mean.mats[[2]][cad3.gr$pos & !cad3.gr$abovecut,] 
vnegchrom<-chrom.mean.mats[[2]][cad3.gr$pos & !cad3.gr$abovecut & cad3.gr$wellbelowcut,] 

posrank=order(scorefunc(cagecountmatlist[[2]][cad3.gr$pos & cad3.gr$abovecut,]),decreasing=T)
negrank=order(scorefunc(cagecountmatlist[[2]][cad3.gr$neg & !cad3.gr$abovecut,]),decreasing=T)

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_cad_boxplots','.jpeg'))

print(Chrom.boxplots(poschrom , negchrom))
dev.off()

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_cad_boxplots_vneg','.jpeg'))

print(Chrom.boxplots(poschrom , vnegchrom))
dev.off()

for(markname in names(chrom.profile.mats[[2]][[setname]]) ){
      poschrom<-chrom.profile.mats[[2]][[markname]][ cad3.gr$pos & cad3.gr$abovecut,]
      negchrom<-chrom.profile.mats[[2]][[markname]][ cad3.gr$pos & !cad3.gr$abovecut,] 

    par(mfrow=c(1,2))
    jpeg(h=1200,w=1920,filename=paste0('analysis/crm_chrom_analysis2/possplit_cad_heatmaps',markname,'.jpeg'))

 
    rankedheatmap(poschrom,col2=chrcol,brks='quantile',filt=0,rankvect=posrank)
    
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ','tf_bound above cutoff','\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)
   

    rankedheatmap(negchrom,col2=chrcol,brks='quantile',filt=0,rankvect=negrank)
 
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ','tf_bound below cutoff','\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)


    dev.off()

}



#boxplots for crms with PollII 




poschrom<-chrom.mean.mats[[1]][crm8008.gr$polII & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$abovecut,]
negchrom<-chrom.mean.mats[[1]][crm8008.gr$polII & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & !crm8008.gr$abovecut,] 
vnegchrom<-chrom.mean.mats[[1]][crm8008.gr$polII & !crm8008.gr$abovecut & crm8008.gr$wellbelowcut,] 
posrank=order(scorefunc(cagecountmatlist[[1]][crm8008.gr$polII & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$abovecut,]),decreasing=T)
negrank=order(scorefunc(cagecountmatlist[[1]][crm8008.gr$polII & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & !crm8008.gr$abovecut,]),decreasing=T)

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_polII_boxplots','.jpeg'))

print(Chrom.boxplots(poschrom , negchrom))

dev.off()

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_polII_boxplots_vneg','.jpeg'))

print(Chrom.boxplots(poschrom , vnegchrom))

dev.off()






for(markname in names(chrom.profile.mats[[1]]) ){
    poschrom<-chrom.profile.mats[[1]][[markname]][crm8008.gr$polII & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$abovecut,]
    negchrom<-chrom.profile.mats[[1]][[markname]][crm8008.gr$polII & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & !crm8008.gr$abovecut,] 

    jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_polII_heatmaps',markname,'.jpeg'))
    par(mfrow=c(1,2))

 
    rankedheatmap(poschrom,col2=chrcol,brks='quantile',filt=0,rankvect=posrank)
    
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ','tf_bound above cutoff','\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)
   

    rankedheatmap(negchrom,col2=chrcol,brks='quantile',filt=0,rankvect=negrank)
 
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ','tf_bound below cutoff','\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)


    dev.off()

}


#or  K27ac


poschrom<-chrom.mean.mats[[1]][crm8008.gr$H3K27ac & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$abovecut,]
negchrom<-chrom.mean.mats[[1]][crm8008.gr$H3K27ac & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & !crm8008.gr$abovecut,] 
vnegchrom<-chrom.mean.mats[[1]][crm8008.gr$H3K27ac & !crm8008.gr$abovecut & crm8008.gr$wellbelowcut,] 
posrank=order(scorefunc(cagecountmatlist[[1]][crm8008.gr$H3K27ac & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$abovecut,]),decreasing=T)
negrank=order(scorefunc(cagecountmatlist[[1]][crm8008.gr$H3K27ac & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & !crm8008.gr$abovecut,]),decreasing=T)

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_H3K27acboxplots','.jpeg'))

print(Chrom.boxplots(poschrom , negchrom))

dev.off()

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_H3K27acboxplots_vneg','.jpeg'))

print(Chrom.boxplots(poschrom , vnegchrom))

dev.off()

for(markname in names(chrom.mean.mats[[setname]]) ){
    poschrom<-chrom.profile.mats[[1]][[markname]][crm8008.gr$H3K27ac & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$abovecut,]
    negchrom<-chrom.profile.mats[[1]][[markname]][crm8008.gr$H3K27ac & crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & !crm8008.gr$abovecut,] 

    jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_H3K27acheatmaps',markname,'.jpeg'))
    par(mfrow=c(1,2))

 
    rankedheatmap(poschrom,col2=chrcol,brks='quantile',filt=0,rankvect=posrank)
    
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ','tf_bound above cutoff','\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)
   

    rankedheatmap(negchrom,col2=chrcol,brks='quantile',filt=0,rankvect=negrank)
 
    abline(v=0.5)
    title(paste0('Input Subtracted ',markname,' for ','tf_bound below cutoff','\n','ranked by cage signal'))
    axis(side=1,at=c(0,0.25,0.5,0.75,1),labels=c('-250bp','-125bp','0','+125bp','+250bp'),las=0,cex.axis=0.8)


    dev.off()

}


#Take most extreme comparison possible - 

intset <- crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak 
vals<-rs[intset]
uq<-quantile(vals,0.95)
bq<-quantile(vals,0.05)
poschrom<-chrom.mean.mats[[1]][intset & rs>=uq,]
negchrom<-chrom.mean.mats[[1]][intset & rs<=bq,] 
print(Chrom.boxplots(poschrom , negchrom))



#Show signal at tss compared to signal at Enhancers
tss.rs<-rowSums(cagecountmatlist[[3]])
posmat<-cagecountmatlist[[1]][intset & rs>=uq,]
negamt<-chrom.mean.mats[[2]][intset & rs<=bq,] 
print(Chrom.boxplots(poschrom , negchrom))



crmmat=cagecountmatlist[[1]][intset,]
tss.s.mat=cagecountmatlist[[3]][ tss.gr$active68 ,]
crmmat=crmmat[rank(rowSums(crmmat)),]
tss.s.mat=tss.s.mat[rank(rowSums(tss.s.mat)),]
crmmat<-crmmat[seq(1,nrow(crmmat),length.out=25),]
tss.s.mat<-tss.s.mat[seq(1,nrow(tssmat),length.out=25),]

#Dotplot comparing 25 CRMs to 25 TSS


