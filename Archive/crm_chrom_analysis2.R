
# source('src/generate_Rle.R')
# source('src/load_annotations.R')

 setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')

source('src/tss_cage_functions.R')
load(file.cage.tag.rles)
load(file.crmgrs)
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
#define window sizes
window.sizes<-c('w25'=25,'w50'=50,'w100'=100,'w200'=200,'w300'=300,'w500'=500)
#exclude the reseq libraries.
cg<-cage.tag.rles[!grepl('reseq',names(cage.tag.rles))]
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

overlapsAny<-function(gr,peaklist){
  0 < rowSums(sapply(peaklist,function(peaks){
    countOverlaps(gr,peaks)
  }))
}

crmgrs$H3K4me3_peak   <-overlapsAny(crmgrs,list(chrompeaks.modencode[['K4me3_4-8h']],chrompeaks[['K4me3_4-6h']],chrompeaks[['K4me3_6-8h']])) 
crmgrs$H3K4me1_peak   <-overlapsAny(crmgrs,list(chrompeaks.modencode[['K4me1_4-8h']],chrompeaks[['K4me1_4-6h']],chrompeaks[['K4me1_6-8h']]))
crmgrs$H3K27ac_peak   <-overlapsAny(crmgrs,list(chrompeaks.modencode[['K27ac_4-8h']],chrompeaks[['K27ac_4-6h']],chrompeaks[['K27ac_6-8h']]))

crmgrs$H379me3_peak   <-overlapsAny(crmgrs,list(chrompeaks[['K79me3_4-6h']],chrompeaks[['K79me3_6-8h']]))
crmgrs$H3K36me3_peak  <-overlapsAny(crmgrs,list(chrompeaks[['K36me3_4-6h']],chrompeaks[['K36me3_6-8h']]))
crmgrs$polII          <-overlapsAny(crmgrs,list(chrompeaks[['PolII_4-6h']],chrompeaks[['PolII_6-8h']]))





##And continous information for graphing later


# DNase as well. --------------------------------------------------
#get the summed dnase reads over 6-8hrs for our crms
crm.dnase.68<-Views(dnase.rles$STG10+dnase.rles$STG11,as(crmgrs,'RangesList'))
crmgrs$dnase.density.68<-unlist(viewSums(crm.dnase.68))/width(crmgrs)
crmgrs$low.dnase<-crmgrs$dnase.density.68<quantile(crmgrs$dnase.density.68,0.2)


chrom.rles.rpgc.sub.merge=chrom.rles.rpgc.sub.merge[!grepl('K36me3|K4me3',names(chrom.rles.rpgc.sub.merge))]
chrom.rles.rpgc.sub.merge=c(chrom.rles.rpgc.sub.merge,'dnase_6.8'=(dnase.rles$STG10+dnase.rles$STG11))


chrom.mean.mats<-list(
  crm8008=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crmgrs,'RangesList'))))}),
  cad=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
)

# chrom.mean.mats.mod<-sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crmgrs,'RangesList'))))})


# chrom.mean.mats2<-get.best.chr.window.mat(resize(crmgrs,width=750,fix='center'),w,cage=chrom.rles.rpgc.sub.merge)
# chrom.mean.mats2<-chrom.mean.mats2/w
# poschrom<- chrom.mean.mats2[crmgrs$abovecut & crmgrs$intergenic,]
# negchrom<-chrom.mean.mats2[!crmgrs$abovecut & crmgrs$intergenic,]
# print(Chrom.boxplots(poschrom , negchrom))

# comparing our datasets to the modencode datasets on the crmgrs ----------


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

#finally construct our negative and positve sets from TF data
crmgrs$in.tf.pos <- !crmgrs$H3K4me3_peak & crmgrs$intergenic & crmgrs$Activity68
crmgrs$in.tf.neg<- !crmgrs$H3K4me3_peak & crmgrs$intergenic & !crmgrs$Activity68 & !crmgrs$Activity46

export(crmgrs[ crmgrs$in.tf.pos],con ='analysis/make_regions_bedfiles/pos.tf.8008crms.bed')
export(crmgrs[ crmgrs$in.tf.neg],con='analysis/make_regions_bedfiles/neg.tf.8008crms.bed')


save(crmgrs,file=file.crmgrs)





# Plots comparing CAGE signal to PolII signal -----------------------------



# qplot(crmgrs$allsum.rpgc,y={a=chrom.mean.mats.mod[,'PolII'];a[a<1]=1;a},log='xy',color=crmgrs$intergenic,
#       xlab='Normalized CAGE tags',ylab='Modencode Whole Embryo PolII 4-8h',main='CAGE signal Vs. Modencode PolII, 8008 CRMs')
# cor(crmgrs$allsum.rpgc,chrom.mean.mats.mod[,'PolII'],method='s')



# qplot(crmgrs$allsum.rpgc,y={a=chrom.mean.mats.mod[,'PolII'];a[a<1]=1;a},log='xy',color=crmgrs$intergenic,
#       xlab='Normalized CAGE tags',ylab='Meso PolII 6_8h ',main='CAGE signal Vs. Modencode PolII, 8008 CRMs')
# cor(crmgrs$allsum.rpgc,chrom.mean.mats[,'PolII_6.8'],method='s')


# qplot(crmgrs$allsum.rpgc[ crmgrs$intergenic],y={a=chrom.mean.mats.mod[,'PolII'][ crmgrs$intergenic];a[a<1]=1;a},log='xy',
#       xlab='Normalized CAGE tags',ylab='Modencode Whole Embryo PolII 4-8h',main='CAGE signal Vs. Modencode PolII, Intergenic 8008 CRMs')
# cor(crmgrs$allsum.rpgc[ crmgrs$intergenic],chrom.mean.mats.mod[,'PolII'][ crmgrs$intergenic],method='s')


# qplot(crmgrs$allsum.rpgc[ crmgrs$intergenic],y={a=chrom.mean.mats.mod[,'PolII'][ crmgrs$intergenic];a[a<1]=1;a},log='xy',
#       xlab='Normalized CAGE tags',ylab='Meso PolII 6_8h ',main='CAGE signal Vs. Modencode PolII, Integergenic 8008 CRMs')
# cor(crmgrs$allsum.rpgc[ crmgrs$intergenic],chrom.mean.mats[,'PolII_6.8'][ crmgrs$intergenic],method='s')














# Now get cage information ------------------------------------------------

#cage data
# crmgrs$allsum.rpgc<-unlist(viewSums(Views(alltags.rpgc$both,as(crmgrs,'RangesList'))))



cagecountmatlist<-list(#get the best window of w in each one.
  crm8008=get.best.window.mat(crmgrs,w,cage=cg),
  cad=get.best.window.mat(  nomcols(cad3.gr)  ,w,cage=cg)
  # cadpos=get.best.window.mat(  nomcols(cad3.gr[cad3.gr$pos])  ,w,cage=cg),
  # cadneg=get.best.window.mat(  nomcols(cad3.gr[ cad3.gr$neg ])  ,w,cage=cg)
)
save(cagecountmatlist,file='cagecountmatlist.object.R')


####create the matrix of normalized counts.
load(file.accession.df)
v<-accession.df$genome.coverage
names(v)<-accession.df$accession

cagecountmatlist.r<-sapply(simplify=F,names(cagecountmatlist),function(reg){
  sapply(1:ncol(cagecountmatlist[[reg]]),function(n){
    cagecountmatlist[[reg]][,n]<-cagecountmatlist[[reg]][,n]/v[names(cg)[n]]
  })
})


#take in a matrix and perform DEseqs normalization on it. The columns should be big,
#if there's a list of matrices require they're all the same size, concatenate them,
#then spit them again at the end
load(file.cage.tag.res.de)
deseqNormalize<-function(m){
  islist=F
  if(is(m,'list')){
    islist=T
    nrows<-sapply(m,nrow)
    splitfacts<-unlist(mapply(names(m),nrows,FUN=function(name,num){rep(name,num)}))
    m<-do.call(rbind,m)
  }



  require(DESeq)
  oldrownames<-rownames(m)
  rownames(m)<-1:nrow(m)#no duplicate names

  conditions=colnames(m)#
  cds=newCountDataSet(m,conditions)#create the DESeq object
  
  cds=estimateSizeFactors(cds)#returns object with size factors
  m=counts(cds,normalized=T)#returns the normalized matrix

  rownames(m) =oldrownames#put the old rownames back in

  #join the list again if necessary
  if(islist){
    nrows=cumsum(nrows)
    starts=c(0,nrows)[1:length(nrows)]+1
    m=sapply(seq_along(nrows),function(n){  m[starts[n]:nrows[n],]})
    names(m)=names(nrows)
  }
  #return normalized matrix or list thereof
  return(m)
}

cagecountmatlist.de<-deseqNormalize(cagecountmatlist.de)









#######################################use all this data to define our datasets

#define our datasets, a list of matrices
possetlist=list(
  tfset=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos,],
  tf.27ac.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K27ac,],
  tf.pol.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$polII,],
  tf.K79.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K79me3,],
  # tf.K36.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K36me3,],
  tf.K4me1.set=cagecountmatlist.r[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K4me1_peak,],
  cad3pos=cagecountmatlist.r[['cad']][cad3.gr$pos,]

)

negsetlist=list(
  full.neg= cagecountmatlist.r[['crm8008']][crmgrs$in.tf.neg & !crmgrs$H3K27ac & crmgrs$low.dnase & ! crmgrs$polII,],
  cad.neg=  cagecountmatlist.r[['cad']][cad3.gr$neg,]
)

#and the negative dataset we'll compare them to
nrow(negsetlist[[1]])

# full.neg.gr<-crmgrs[crmgrs$in.full.neg]
# full.neg.gr$name<-paste0('full.neg.',1:length(full.neg.gr))
# export(full.neg.gr,'analysis/crm_chrom_analysis2/full.neg.bed')





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
# list(pos=c(rep(T,30),rep(F,9)),neg=c(rep(F,100),rep(T,80),rep(F,7)))->above.cut.vects
Barplot.from.matrices<-function(posmat,negmat,logy=T, yname='Normalized Cage Signal 95% cl', xname='CRM',above.cut.vects=NA, barsize=2,cut=1,  bar=T,reorder=F  ){
  
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
  
  bothmelt$logval<-bothmelt$value
  bothmelt$logval[ bothmelt$logval==0]<-1#note that we take the logarithm of the logvals with ggplot later

  yscale=ifelse(logy, scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000)) ,scale_y_continuous(name=yname))   


  if(length(above.cut.vects)!=1){
    return(
      #line graph of our se
      ggplot(data=bothmelt,aes(x=Var1,y=logval,color=set,alpha=abovecut))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=barsize)+
        yscale+
        #scale_y_continuous(name=yname)+
        scale_x_discrete(name=xname,labels='')+
        ggtitle('95% cl for mean of RPGC normalized Signal')
    )
  }else{
    return(
       #line graph of our se
      ggplot(data=bothmelt,aes(x=Var1,y=logval,color=set))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=barsize)+
        yscale+
        #scale_y_continuous(name=yname)+
        scale_x_discrete(name=xname,labels='')+
        ggtitle('95% cl for mean of RPGC normalized Signal')
    )
  }  
}
#limits=c(0,100000)

# Barplot.from.matrices(posmat[1:10,],negmat[1:10,],bar=T)
Chrom.boxplots<-function(posmat,negmat,ylog=T ,lowlim=F,yname='Mean Inp-Sub Chromatin Signals for crms above/below cutoff ', xname='Set',tit=''){
  
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

  yscale=ifelse(logy,scale_y_log10(name=yname),scale_y_continuous(name=yname))   

  #line graph of our se
  ggplot(data=bothmelt,aes(x=set,facet=Var2,y=value,color=Var2))+geom_boxplot()+
    yscale+
    facet_wrap(~Var2)+
    #scale_y_continuous(name=yname)+
    scale_x_discrete(name=xname,labels=c('Above','Below'))+     
    geom_text(data=num.df ,aes(x=set,y=yval,label=number))+
    ggtitle(tit)
}


number.not.zero<-function(mat){apply(mat,1,function(x)sum(x!=0))}

#for now we'll just score the matrices by adding up the columns (libraries)
scorefunc<-rowSums
#setn=names(possetlist)[6]

#what if we have a list of lists, the bottom lists are going to contain
#the type of scoring used, the cutoff, vectors describing pos and neg for cad and crms
#the descriminances and accuracy.
negn=names(negsetlist)[2]
setn=names(possetlist)[6]
cutoffs<-sapply(negsetlist,function(x)sapply(possetlist,function(z)list()))
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
    
    #now define our above and below cutoff sets.
    cutoffs[[negn]][[setn]]=list(
      crm8008=scorefunc(cagecountmatlist.r[['crm8008']])>cutoff,
      cad=scorefunc(cagecountmatlist.r[['cad']])>cutoff
    )


    #create vectors showing which crms are above and below the threshold
    pos.abovecut<-pos>cutoff
    neg.abovecut<-neg>cutoff
    

    #and produce the boxplots for the other chromatin marks
    pdf(paste0('analysis/crm_chrom_analysis2/crm_chrom_boxplot',setn,negn,'.pdf'))
    
    poschrom<- chrom.mean.mats[[1]][cutoffs[[negn]][[setn]][['crm8008']] & crmgrs$intergenic,]
    negchrom<-chrom.mean.mats[[1]][!cutoffs[[negn]][[setn]][['crm8008']] & crmgrs$intergenic,]

    print(Chrom.boxplots(poschrom , negchrom,tit=paste0('All Intergenic CRMs, Split by cutoff ',cutoff)))
    
    dev.off()
    # poschrom<- chrom.mean.mats2[crmgrs$abovecut & crmgrs$intergenic,]
    # negchrom<-chrom.mean.mats2[!crmgrs$abovecut & crmgrs$intergenic,]
    # print(Chrom.boxplots(poschrom , negchrom))  
    
    cat(setn)
    cat(negn)
    cat(' ... ')
    
  }  

}



#Does it make a difference to use the total tags vs the number of lines with something?
scorefunc<-number.not.zero
poschrom<- chrom.mean.mats[[1]][scorefunc(cagecountmatlist.r[[1]])>40 & crmgrs$intergenic,]
negchrom<-chrom.mean.mats[[1]][!scorefunc(cagecountmatlist.r[[1]])>2 & crmgrs$intergenic,]
print(Chrom.boxplots(poschrom , negchrom,tit='tmp'))



rs=rowSums(cagecountmatlist.r[[1]])
rs.nn=rowSums(cagecountmatlist[[1]])

nnz=number.not.zero(cagecountmatlist.r[[1]])

qplot(rs,nnz,log='x')
qplot(rs.nn,nnz,log='x')

(rank(rs),rank(nnz),log='',xlab='rank of normalized tag total',ylab='rank of linecount (lines with 1+ tags)',main='Rank in total tags vs. Rank in number of lines with tags')
qplot(rank(rs.nn),rank(nnz),log='')

#for each line, how many of our crms have anything?



crmgrs$abovecut<-scorefunc(cagecountmatlist.r[[1]])>cutoff
cad3.gr$abovecut<-scorefunc(cagecountmatlist.r[[2]])>cutoff
#are we getting problems with single libraries boosting negatives into significance?
#find lines where 80% or more of our signal comes from a single library....
rs<-rowSums(cagecountmatlist.r[[1]])
rm<-apply(cagecountmatlist.r[[1]],1,max)
single.lib.dominant=  rm/rs>0.8
single.lib.dominant=single.lib.dominant%in%T
meanaboveone<-apply(cagecountmatlist[[1]],1,mean)>1
twolibswithtwo<-apply(cagecountmatlist[[1]],1,function(x) sum(x>0)>1 )
libnum<-apply(cagecountmatlist[[1]],1,function(x) sum(x>0) )
abovecut<-rs>350 &

#no cases in which a single library is 80 of the signal 
sum(crmgrs$abovecut & single.lib.dominant)
#no cases above our cutoff without at least two libraries above 2 tags
sum(crmgrs$abovecut & !twolibswithtwo)
#checking what the maximum is for these two groups
max(rs[single.lib.dominant])
max(rs[!twolibswithtwo])



#Now for cads and crms
#for cads we use the 500 bp around cage center
#rank by cage signal
#and then produce heatmaps for various marks
# chrom.rles.rpgc.sub.merge->chrom.rles.rpgc.sub.merge.bak
# crmgrs->crmgrs.bak


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





chrom.profile.mats<-list(
  crm8008=lapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){
    profile.matrix(chrom.rle,center.gr(alltags$both,resize(crmgrs,wi=500,fi='center')))
  }),
  cad=lapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){
    profile.matrix(chrom.rle,center.gr(alltags$both,cad.centered.gr))
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

crmgrs$wellbelowcut<-rs<quantile(rs[rs<cutoff],0.25)#those rows which are in the bottom half of those below the cutoff
rsc<-rowSums(cagecountmatlist.r[[2]])
cad3.gr$wellbelowcut<-rsc<quantile(rsc[rsc<cutoff],0.25)





poschrom<-chrom.mean.mats[[1]][crmgrs$in.tf.pos & crmgrs$abovecut,]
negchrom<-chrom.mean.mats[[1]][crmgrs$in.tf.pos & !crmgrs$abovecut,] 
vnegchrom<-chrom.mean.mats[[1]][crmgrs$in.tf.pos & !crmgrs$abovecut & crmgrs$wellbelowcut,] 
posrank=order(scorefunc(cagecountmatlist[[1]][crmgrs$in.tf.pos & crmgrs$abovecut,]),decreasing=T)
negrank=order(scorefunc(cagecountmatlist[[1]][crmgrs$in.tf.pos & !crmgrs$abovecut,]),decreasing=T)

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_tf_boxplots','.jpeg'))

print(Chrom.boxplots(poschrom , negchrom))

dev.off()

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_tf_boxplots_vneg','.jpeg'))

print(Chrom.boxplots(poschrom , vnegchrom))

dev.off()

for(markname in names(chrom.profile.mats[[1]]) ){
    poschrom<-chrom.profile.mats[[1]][[markname]][crmgrs$in.tf.pos & crmgrs$abovecut,]
    negchrom<-chrom.profile.mats[[1]][[markname]][crmgrs$in.tf.pos & !crmgrs$abovecut,] 

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




poschrom<-chrom.mean.mats[[1]][crmgrs$polII & crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$abovecut,]
negchrom<-chrom.mean.mats[[1]][crmgrs$polII & crmgrs$intergenic & !crmgrs$H3K4me3_peak & !crmgrs$abovecut,] 
vnegchrom<-chrom.mean.mats[[1]][crmgrs$polII & !crmgrs$abovecut & crmgrs$wellbelowcut,] 
posrank=order(scorefunc(cagecountmatlist[[1]][crmgrs$polII & crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$abovecut,]),decreasing=T)
negrank=order(scorefunc(cagecountmatlist[[1]][crmgrs$polII & crmgrs$intergenic & !crmgrs$H3K4me3_peak & !crmgrs$abovecut,]),decreasing=T)

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_polII_boxplots','.jpeg'))

print(Chrom.boxplots(poschrom , negchrom))

dev.off()

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_polII_boxplots_vneg','.jpeg'))

print(Chrom.boxplots(poschrom , vnegchrom))

dev.off()

for(markname in names(chrom.profile.mats[[1]]) ){
    poschrom<-chrom.profile.mats[[1]][[markname]][crmgrs$polII & crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$abovecut,]
    negchrom<-chrom.profile.mats[[1]][[markname]][crmgrs$polII & crmgrs$intergenic & !crmgrs$H3K4me3_peak & !crmgrs$abovecut,] 

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
poschrom<-chrom.mean.mats[[1]][crmgrs$H3K27ac & crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$abovecut,]
negchrom<-chrom.mean.mats[[1]][crmgrs$H3K27ac & crmgrs$intergenic & !crmgrs$H3K4me3_peak & !crmgrs$abovecut,] 
vnegchrom<-chrom.mean.mats[[1]][crmgrs$H3K27ac & !crmgrs$abovecut & crmgrs$wellbelowcut,] 
posrank=order(scorefunc(cagecountmatlist[[1]][crmgrs$H3K27ac & crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$abovecut,]),decreasing=T)
negrank=order(scorefunc(cagecountmatlist[[1]][crmgrs$H3K27ac & crmgrs$intergenic & !crmgrs$H3K4me3_peak & !crmgrs$abovecut,]),decreasing=T)

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_H3K27acboxplots','.jpeg'))

print(Chrom.boxplots(poschrom , negchrom))

dev.off()

jpeg(h=1200,w=1920,paste0('analysis/crm_chrom_analysis2/possplit_H3K27acboxplots_vneg','.jpeg'))

print(Chrom.boxplots(poschrom , vnegchrom))

dev.off()

for(markname in names(chrom.mean.mats[[setname]]) ){
    poschrom<-chrom.profile.mats[[1]][[markname]][crmgrs$H3K27ac & crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$abovecut,]
    negchrom<-chrom.profile.mats[[1]][[markname]][crmgrs$H3K27ac & crmgrs$intergenic & !crmgrs$H3K4me3_peak & !crmgrs$abovecut,] 

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












# Views( alltags above.cutoff.gr[10])



# v1<-rnorm(100)
# v2<-rnorm(100)*0.8 +v1
# plot(1:3)
# vectlist=list(v1,v2)
# rank.order.plot<-function(vectlist){
#   par(mfrow=c(1,2))
#   v1=sort(v1,d=T)
#   v2=sort(v2,d=F)
#   r1=rank(v1)
#   r2=rank(v2)
#   image(matrix(r1,ncol=length(r1),nrow=1))
#   image(matrix(r2[r1],ncol=length(r2),nrow=1))

# }








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
export(full.neg.gr,'analysis/crm_chrom_analysis2/full.neg.bed')


#gr with our negatives
above.cutoff.gr<-crmgrs[ crmgrs$abovecut & crmgrs$intergenic]
#order them by their cage signal
above.cutoff.gr<-above.cutoff.gr[order(above.cutoff.gr$allsum.rpgc,decreasing=T)]
#appropriate names
above.cutoff.gr$name<-paste0('above.cutoff.',1:length(above.cutoff.gr))
mcols(above.cutoff.gr)<-mcols(above.cutoff.gr)[,c('name')]
above.cutoff.gr$name<-above.cutoff.gr$value
export(above.cutoff.gr,'analysis/crm_chrom_analysis2/above.cutoff.bed')


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







