##Todo on this script
#Incorporate modencode stuff
#Modencode RNAseq
#1 load functions,datasets
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
source('src/generate_chromdata_cage.R')
load(file.cg.pl)
load(file.
library(ROCR,lib.loc='~/Harnett/R/')
cg<-cg.pl.map
rm(cg.pl)
w=300
#precission cutoff when setting cutoff
prec.cutoff<-0.95

##2   get cage data summed over windwos for our regions, at various timestages
tss.gr<-tss.gr[!duplicated(tss.gr)]
tss.gr<-tss.gr[seqnames(tss.gr)%in% bigchrs]
tss.gr=sort(tss.gr)
tps=unique(accession.df$timepoint)[1:3]
cagecountmatlist<-sapply(simplify=F,tps,function(tp){
    list(
        crm8008=getBestWindowMat(crm8008.gr[],w,cage=cg[accession.df$timepoint==tp][]),
        cad=getBestWindowMat(  nomcols(cad3.gr[])  ,w,cage=cg[accession.df$timepoint==tp][]),
    ) 
})
for(tp in tps){
  cagecountmatlist[[tp]]$tss<-tss.gr$tsscage[,accession.df$time==tp]
}

save(cagecountmatlist,file='cagecountmatlist.object.R')
# load('cagecountmatlist.object.R')
#calculate means for our matrices
meanlist<-lapply(cagecountmatlist,function(mat){row})
message('cage info collected')







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

crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak

##And continous information for graphing later




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








# DNase as well. --------------------------------------------------
#get the summed dnase reads over 6-8hrs for our crms
crm.dnase.68<-Views(dnase.rles$STG10+dnase.rles$STG11,as(crm8008.gr,'RangesList'))
crm8008.gr$dnase.density.68<-unlist(viewSums(crm.dnase.68))/width(crm8008.gr)
crm8008.gr$low.dnase<-crm8008.gr$dnase.density.68<quantile(crm8008.gr$dnase.density.68,0.2)
chrom.rles.rpgc.sub.merge=chrom.rles.rpgc.sub.merge[!grepl('K36me3',names(chrom.rles.rpgc.sub.merge))]
chrom.rles.rpgc.sub.merge=c(chrom.rles.rpgc.sub.merge,'dnase_6.8'=(dnase.rles$STG10+dnase.rles$STG11))
chrom.mean.mats<-list(
  crm8008=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crm8008.gr,'RangesList'))))}),
  cad=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
)
chrom.mean.mats$crm8008 <- cbind(chrom.mean.mats$crm8008,data.frame(NearestGeneExpr=GetNearestExpression(crm8008.gr)))
chrom.mean.mats$cad <- cbind(chrom.mean.mats$cad,data.frame(NearestGeneExpr=GetNearestExpression(cad3.gr)))

save(chrom.mean.mats,file='chrom.mean.mats.object.R')

load('chrom.mean.mats.object.R')
chrom.mean.mats.modencode<-list(
  crm8008=sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crm8008.gr,'RangesList'))))}),
  cad=sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
)

save(chrom.mean.mats.modencode,file='chrom.mean.mats.modencode.object.R')
load('chrom.mean.mats.modencode.object.R')








# Barplot.from.matrices(posmat[1:10,],negmat[1:10,],bar=T)
Chrom.boxplots <-function(posmat,negmat, lowlim=F,yname='Mean Inp-Sub Chromatin Signals for crms above/below cutoff ', xname='Set',tit=''){
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Above',nrow(posmelt)),rep('Below',nrow(negmelt)))
  bothmelt$nset<-c(rep(nrow(posmat),nrow(posmelt)),rep(nrow(negmat),nrow(negmelt)))
  bothmelt$zscore=bothmelt$value
  #we need to do something iwth our zero/negative values to do a log scale....
  if(lowlim==F){
    lowlim=min(bothmelt$zscore)
  }
  # bothmelt$value[bothmelt$value<=0] <- lowlim
  
  num.df <-  data.frame(set=c('Above','Below'),number=c(nrow(posmat),nrow(negmat)))
  num.df <- do.call(rbind, lapply(unique(bothmelt$Var2),function(v) cbind(num.df,v) ))
  num.df$Var2<-num.df$v
  num.df$yval<-min(bothmelt$zscore)

  #for each pair, do a wilcox test and append an asterix if significant
  pvals=sapply(1:ncol(posmat),function(i){wilcox.test(posmat[,i],negmat[,i])$p.value})
  pvals=pvals<0.05
  pvals=ifelse(pvals,'*','')
  names(pvals)=colnames(posmat)
  num.df$number=paste0(num.df$number,pvals[as.character(num.df$Var2)])
  #line graph of our se
  ggplot(data=bothmelt,aes(x=set,facet=Var2,y=zscore,color=Var2))+geom_boxplot()+
    # scale_y_log10(name=yname)+    
    #scale_y_continuous(name=yname)+
    scale_x_discrete(name=xname,labels=c('Above','Below'))+     
    geom_text(data=num.df ,aes(x=set,y=yval,label=number))+
    ggtitle(tit)+
    theme(strip.text=element_text(size=10))+
    facet_wrap(~Var2, scale = 'free_y')
}

chrom.mean.mats = sapply(chrom.mean.mats,function(x){apply(x,2,FUN=qnormvect)})
chrom.mean.mats.modencode = sapply(chrom.mean.mats.modencode,function(x){apply(x,2,FUN=qnormvect)})



chromscatters<-function(chrommat,eRNAvect,qnorm.erna=T,qnorm.chrom=F,tit){
    if(qnorm.erna){eRNAvect<-qnormvect(eRNAvect)}
    if(qnorm.chrom){chrommat=apply(chrommat,2,qnormvect)}
    chrommelt <- melt(chrommat)
    chrommelt$eRNA <- eRNAvect
    #
    colnames(chrommelt)<-c('CRM','Mark','Level','eRNA')
    #add correlation figure to plot
    cors <- sapply(unique(chrommelt$Mark),function(mark){cor(chrommelt$Level[chrommelt$Mark==mark],eRNAvect,meth='s')})
    cors = sapply(cors,round,4)
    names(cors)=unique(chrommelt$Mark)
    chrommelt$correlation <- cors[as.character(chrommelt$Mark)]

    stopifnot('numeric' %in%is(chrommelt$eRNA))
    ggplot(data=chrommelt,aes(x=eRNA,facet=Mark,y=Level,color=Mark,alpha=0.1))+geom_point()+
    facet_wrap(~Mark, scale = 'free_y')+
    scale_y_continuous(name='Zscore')+
    geom_text(data=chrommelt ,aes(x=0,y=-2,label=paste('rho =',correlation)))+
    ggtitle(tit)
   }
# m=matrix(1:24,ncol=4)
# colnames(m)=letters[1:4]
# eRNAvect=seq(0,1,length.out=6)
# pdf('tmp.pdf');chromscatters(m,eRNAvect,tit='test');dev.off()
# # #


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





#######################################use all this data to define our datasets
#define our datasets, a list of matrices
possetlist=list(
  tfset=cagecountmatlist[['crm8008','68h']][crm8008.gr$in.tf.pos,],
  tf.27ac.set=cagecountmatlist[['crm8008','68h']][crm8008.gr$in.tf.pos & crm8008.gr$H3K27ac,],
  tf.pol.set=cagecountmatlist[['crm8008','68h']][crm8008.gr$in.tf.pos & crm8008.gr$polII,],
  tf.K79.set=cagecountmatlist[['crm8008','68h']][crm8008.gr$in.tf.pos & crm8008.gr$H3K79me3,],
  # tf.K36.set=cagecountmatlist[['crm8008']][crm8008.gr$in.tf.pos & crm8008.gr$H3K36me3,],
  tf.K4me1.set=cagecountmatlist[['crm8008','68h']][crm8008.gr$in.tf.pos & crm8008.gr$H3K4me1_peak,],
  cad3pos=cagecountmatlist[['cad','68h']][cad3.gr$pos,]
)
negsetlist=list(
  full.neg= cagecountmatlist[['crm8008','68h']][crm8008.gr$in.tf.neg & !crm8008.gr$H3K27ac & crm8008.gr$low.dnase & ! crm8008.gr$polII,],
  cad.neg=  cagecountmatlist[['cad','68h']][cad3.gr$neg,]
)





