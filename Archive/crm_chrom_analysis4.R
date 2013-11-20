##Todo on this script
#Incorporate modencode stuff
#Modencode RNAseq
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')

source('src/generate_chromdata_cage.R')

load(file.cg.pl)
load(file.accession.df)
load(file.crmgrs)
load(file.cad3)
load(file.tss)
library(ROCR,lib.loc='~/Harnett/R/')

#cg<-cage.tag.rles
cg<-cg.pl
#windowsize
w=100
#precission cutoff when setting cutoff
prec.cutoff<-0.95





#get cage data summed over windwos for our regions, at various timestages
load(file.tss)
tss.gr<-tss.gr[!duplicated(tss.gr)]
tss.gr<-tss.gr[seqnames(tss.gr)%in% bigchrs]
tss.gr=sort(tss.gr)
tps=accession.df$timepoint

cagecountmatlist<-list(#get the best window of w in each one.
  sapply(unique(tps),function(tp){
    list(crm8008=get.best.window.mat(crmgrs[],w,cage=cg[accession.df$timepoint==tp][]),
        cad=get.best.window.mat(  nomcols(cad3.gr[])  ,w,cage=cg[accession.df$timepoint==tp][]),
        tss=get.best.window.mat(  nomcols(tss.gr[])  ,w,cage=cg[accession.df$timepoint==tp][])
        )
  # cadpos=get.best.window.mat(  nomcols(cad3.gr[cad3.gr$pos])  ,w,cage=cg),
  # cadneg=get.best.window.mat(  nomcols(cad3.gr[ cad3.gr$neg ])  ,w,cage=cg)
  })
)
save(cagecountmatlist,file='cagecountmatlist.object.R')

load('cagecountmatlist.object.R')
cagecountmatlist=cagecountmatlist[[1]]#can fix this later

meanlist<-lapply(cagecountmatlist,function(mat){row})

#load('cagecountmatlist.object.R')


#first do a dotplot comparing levels at CRMs vs TSS at 68 hours
#

tp=unique(accession.df$timepoint)[[3]]








#CRMs vs TSS.

for(tp in unique(accession.df$timepoint)){

  mag=7#orders of magnitude to plot

  crmorder=order(rowMeans(cagecountmatlist['crm8008',tp][[1]][crmgrs$intergenic,]),decreasing=F)
  tssorder=order(rowMeans(cagecountmatlist['tss',tp][[1]]),decreasing=F)

  crmrows=crmorder[seq(from=1,to=length(crmorder),length.out=50)]
  tssrows=tssorder[seq(from=1,to=length(tssorder),length.out=50)]

  m = rbind( cagecountmatlist['crm8008',tp][[1]][crmgrs$intergenic,][crmrows,],cagecountmatlist['tss',tp][[1]][tssrows,])
  rownames(m)<-1:nrow(m)
  m=melt(m)
  m$value[m$value==0]<-1
  #now color are CRMs vs TSS
  colvect=rep(c('Intergenic CRM','TSS'),each=50)
  # setwd('analysis/crm_chrom_analysis4')
  # pdf(paste0('crms_v_tss_levels',tp,'.pdf'))


  ggplot(data=m,aes(x=Var1,y=value,color=colvect[Var1]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='Log10 normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
  scale_color_discrete(name='set',labels=unique(colvect))+scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
  ggtitle(paste0('Cage Signal at 50 quantiles for CRMs vs TSS at ',tp ))
  # dev.off()
  # setwd('../..')

}


#Intergenic vs. Genic TSS


tp=unique(accession.df$timepoint)[[1]]

  mag=7#orders of magnitude to plot

  intcrmorder=order(rowMeans(cagecountmatlist['crm8008',tp][[1]][crmgrs$intergenic,]),decreasing=F)
  crmorder=order(rowMeans(cagecountmatlist['crm8008',tp][[1]][!crmgrs$intergenic,]),decreasing=F)

  intcrmrows=intcrmorder[seq(from=1,to=length(intcrmorder),length.out=50)]
  crmrows=crmorder[seq(from=1,to=length(crmorder),length.out=50)]

  m = rbind( cagecountmatlist['crm8008',tp][[1]][crmgrs$intergenic,][intcrmrows,],
    cagecountmatlist['crm8008',tp][[1]][!crmgrs$intergenic,][crmrows,])

  rownames(m)<-1:nrow(m)

  m=melt(m)
  m$value[m$value==0]<-1
  #now color are CRMs vs TSS
  colvect=rep(c('Intergenic CRM','IntraGenic CRM'),each=50)
  # setwd('analysis/crm_chrom_analysis4')
  # pdf(paste0('crms_v_tss_levels',tp,'.pdf'))


  ggplot(data=m,aes(x=Var1,y=value,color=colvect[Var1]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='Log10 normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
  scale_color_discrete(name='set',labels=unique(colvect))+scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
  ggtitle(paste0('Cage Signal at 50 quantiles for CRMs vs TSS at ',tp ))
  # dev.off()
  # setwd('../..')




#mixing them


tp=unique(accession.df$timepoint)[[1]]

mag=7#orders of magnitude to plot
mat=cagecountmatlist['crm8008',tp][[1]]
o=order(rowMeans(mat),decreasing=F)
crmrows=o[seq(from=1,to=length(o),length.out=50)]
mat=mat[o,]#sort the matrix
rows=7907:8007#take rows from the matrix
mat=mat[rows,]
rownames(mat)<-1:nrow(mat)
m=melt(mat)
m$value[m$value==0]<-1
#now color are CRMs vs TSS
colvect=crmgrs$intergenic[o][rows]
# setwd('analysis/crm_chrom_analysis4')
# pdf(paste0('crms_v_tss_levels',tp,'.pdf'))
ggplot(data=m,aes(x=Var1,y=value,color=colvect[Var1]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='Log10 normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
scale_color_discrete(name='set',labels=unique(colvect))+scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
ggtitle(paste0('Cage Signal CRMs',tp ))





#Cad positives vs Cad negatives
tp=unique(accession.df$timepoint)[[1]]

mag=3#orders of magnitude to plot
mat=cagecountmatlist['cad',tp][[1]]
mat=rbind(mat[cad3.gr$pos,],mat[cad3.gr$neg,])
rownames(mat)<-1:nrow(mat)
  m=melt(mat)
  m$value[m$value==0]<-1
#now color are CRMs vs TSS
colvect=c(rep('cadpos',sum(cad3.gr$pos)),rep('cadneg',sum(cad3.gr$neg)))
# setwd('analysis/crm_chrom_analysis4')
dev.off()
pdf(paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/cadpos_vs_cadneg','.pdf'))
ggplot(data=m,aes(x=Var1,y=value,color=colvect[Var1]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='Log10 normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
scale_color_discrete(name='set',labels=unique(colvect))+scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
ggtitle(paste0('Cage Signal at CAD',tp ))
dev.off()










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
crmgrs$H3K4me3_peak   <-overlapsAny(crmgrs,list(chrompeaks.modencode[['K4me3_4-8h']], chrompeaks[['K4me3_4-6h']] ,chrompeaks[['K4me3_6-8h']]),maxgap=50) 
crmgrs$H3K4me1_peak   <-overlapsAny(crmgrs,list(chrompeaks.modencode[['K4me1_4-8h']], chrompeaks[['K4me1_4-6h']] ,chrompeaks[['K4me1_6-8h']]),maxgap=50)
crmgrs$H3K27ac_peak   <-overlapsAny(crmgrs,list(chrompeaks.modencode[['K27ac_4-8h']], chrompeaks[['K27Ac_4-6h']] ,chrompeaks[['K27Ac_6-8h']]),maxgap=50)
crmgrs$H379me3_peak   <-overlapsAny(crmgrs,list(chrompeaks[['K79me3_4-6h']],chrompeaks[['K79me3_6-8h']]),maxgap=50)
crmgrs$H3K36me3_peak  <-overlapsAny(crmgrs,list(chrompeaks[['K36me3_4-6h']],chrompeaks[['K36me3_6-8h']]),maxgap=50)
crmgrs$polII          <-overlapsAny(crmgrs,list(chrompeaks[['PolII_4-6h']],chrompeaks[['PolII_6-8h']]),maxgap=50)

##And continous information for graphing later






# DNase as well. --------------------------------------------------
#get the summed dnase reads over 6-8hrs for our crms
crm.dnase.68<-Views(dnase.rles$STG10+dnase.rles$STG11,as(crmgrs,'RangesList'))
crmgrs$dnase.density.68<-unlist(viewSums(crm.dnase.68))/width(crmgrs)
crmgrs$low.dnase<-crmgrs$dnase.density.68<quantile(crmgrs$dnase.density.68,0.2)


chrom.rles.rpgc.sub.merge=chrom.rles.rpgc.sub.merge[!grepl('K36me3',names(chrom.rles.rpgc.sub.merge))]
chrom.rles.rpgc.sub.merge=c(chrom.rles.rpgc.sub.merge,'dnase_6.8'=(dnase.rles$STG10+dnase.rles$STG11))


chrom.mean.mats<-list(
  crm8008=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crmgrs,'RangesList'))))}),
  cad=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
)
save(chrom.mean.mats,file='chrom.mean.mats.object.R')

load('chrom.mean.mats.object.R')
chrom.mean.mats.modencode<-list(
  crm8008=sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crmgrs,'RangesList'))))}),
  cad=sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
)
save(chrom.mean.mats.modencode,file='chrom.mean.mats.modencode.object.R')
load('chrom.mean.mats.modencode.object.R')









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


#we can now do the boxplots for our quantiles based on the normalized cage data


#for 8008 crms

intset=crmgrs$intergenic

rowranks=rank(rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,]))

#define two logical vectors for the postivie and negative
top=rowranks>=quantile(rowranks,0.95)
bottom=rowranks<=quantile(rowranks,0.05)

poschrom<-chrom.mean.mats[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats[['crm8008']][intset,][bottom,]

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/crm_5quant_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()



#using the modencode data



#for CADs

rowranks=rank(rowMeans(cagecountmatlist['cad','68h'][[1]][cad3.gr$intergenic,]))

#define two logical vectors for the postivie and negative

top=rowranks>quantile(0.95,rowranks)
bottom=rowranks<quantile(0.05,rowranks)













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





#######################################use all this data to define our datasets
#define our datasets, a list of matrices
possetlist=list(
  tfset=cagecountmatlist[['crm8008','68h']][crmgrs$in.tf.pos,],
  tf.27ac.set=cagecountmatlist[['crm8008','68h']][crmgrs$in.tf.pos & crmgrs$H3K27ac,],
  tf.pol.set=cagecountmatlist[['crm8008','68h']][crmgrs$in.tf.pos & crmgrs$polII,],
  tf.K79.set=cagecountmatlist[['crm8008','68h']][crmgrs$in.tf.pos & crmgrs$H3K79me3,],
  # tf.K36.set=cagecountmatlist[['crm8008']][crmgrs$in.tf.pos & crmgrs$H3K36me3,],
  tf.K4me1.set=cagecountmatlist[['crm8008','68h']][crmgrs$in.tf.pos & crmgrs$H3K4me1_peak,],
  cad3pos=cagecountmatlist[['cad','68h']][cad3.gr$pos,]
)

negsetlist=list(
  full.neg= cagecountmatlist[['crm8008','68h']][crmgrs$in.tf.neg & !crmgrs$H3K27ac & crmgrs$low.dnase & ! crmgrs$polII,],
  cad.neg=  cagecountmatlist[['cad','68h']][cad3.gr$neg,]
)







#what if we have a list of lists, the bottom lists are going to contain
#the type of scoring used, the cutoff, vectors describing pos and neg for cad and crms
#the descriminances and accuracy.
negn=names(negsetlist)[1]
setn=names(possetlist)[6]
posmat=possetlist[[setn]]#get our data
negmat=negsetlist[[negn]]


dotplot.from.matrices(posmat,negmat)
simple.barplot.matrices(posmat,negmat)






#first let's just produce 





cutoffs<-sapply(negsetlist,function(x)sapply(possetlist,function(z)list()))
ncutoffs<-sapply(negsetlist,function(x)sapply(possetlist,function(z)list()))

# Now go through each of our sets and produce our plots -------------------
for(negn in names(negsetlist)){
  for(setn in names(possetlist)){
    #dir.create('analysis/crm_chrom_analysis2/')
    mag=3
    
    posmat=possetlist[[setn]]#get our data
    negmat=negsetlist[[negn]]
    pos=scorefunc(posmat)#for now just sum across lines
    neg=scorefunc(negmat)
    

    m=rbind(posmat,negmat)
    rownames(m)<-1:nrow(m)
    colvect=c(rep(setn,nrow(posmat)),rep(negn,nrow(negmat)))
    m=melt(m)
    m$value[m$value==0]<-1
  # setwd('analysis/crm_chrom_analysis4')
   pdf(paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/posneg_barplot_',setn,negn,'.pdf'))
   par(asp=2)
    ggplot(data=m,aes(x=Var1,y=value,color=colvect[Var1]))+coord_trans(y = "log10")+scale_y_continuous(name='Log10 normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
    scale_color_discrete(name='set',labels=unique(colvect))+scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
    ggtitle(paste0('Cage Signal for ',setn,' vs ',negn ))
  dev.off()

    
    pdf(paste0('analysis/crm_chrom_analysis4/crmrocplot',setn,negn,'.pdf'))
    ROCfunc(pos,neg,
            main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,setn),
            sub=paste0('windowsize ',w))
    dev.off()
    
    pdf(paste0('analysis/crm_chrom_analysis4/crm_prec_rec_plot',setn,negn,'.pdf'))
    cutoff<-PreRecfunc(pos,neg,
                       main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,setn),
                       sub=paste0('windowsize ',w))
    cutoff
    dev.off()
    
    cutoffs[[negn]][[setn]][['crm8008']]  <-scorefunc(cagecountmatlist[[1]])>cutoff
    cutoffs[[negn]][[setn]][['cad3']]     <-scorefunc(cagecountmatlist[[2]])>cutoff

    #and produce the boxplots for the other chromatin marks
    pdf(paste0('analysis/crm_chrom_analysis2/crm_chrom_boxplot',setn,negn,'.pdf'))
    poschrom<- chrom.mean.mats[[1]][cutoffs[[negn]][[setn]][['crm8008']] & crmgrs$intergenic,]
    negchrom<-chrom.mean.mats[[1]][!cutoffs[[negn]][[setn]][['crm8008']] & crmgrs$intergenic,]
    print(Chrom.boxplots(poschrom , negchrom,tit=paste0('All Intergenic CRMs, Split by cutoff ',cutoff)))
    dev.off()
    
    cat(setn)
    cat(negn)
    cat(' ... ')
  }  
}

#2.95 - cutoff for now

mean(rowMeans(cagecountmatlist[['crm8008','68h']][crmgrs$intergenic & ! crmgrs$H3K4me3_peak,])>cutoff)

sum(rowMeans(cagecountmatlist[['crm8008','68h']][crmgrs$intergenic & ! crmgrs$H3K4me3_peak,])>cutoff)



mean(rowMeans(cagecountmatlist[['cad','68h']][cad3.gr$neg,])>cutoff)

sum(rowMeans(cagecountmatlist[['cad','68h']][cad3.gr$pos,])>cutoff)


intset=crmgrs$intergenic & !crmgrs$H3K4me3_peak

#hunt out some for browser viewing
crmgrs$name<-gsub('^(\\d+).*','\\1',x=crmgrs$name)
save(crmgrs,file=file.crmgrs)

crmgrs$name[crmgrs$countmean=rowMeans(cagecountmatlist[['crm8008','68h']])]


cad3.gr$countmeans=rowMeans(cagecountmatlist[['cad','68h']])






intset=crmgrs$intergenic & !crmgrs$H3K4me3_peak
rs=rowMeans(cagecountmatlist['crm8008',tp][[1]][,])
rstofind=rs
rstofind[!intset]=0
topdudes=order(decreasing=T,rstofind)
topcrms=crmgrs[topdudes]
topcrms$name<-paste0('topcrm',1:length(topcrms))
export(topcrms,'/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/topcrms.bed')



  #show the difference between above cutoff and below cutoff guys in a dotplot
  tp='68h'
  mag=4#orders of magnitude to plot
  abovecut=rowMeans(cagecountmatlist['crm8008',tp][[1]][,])>cutoff

  abovemat=cagecountmatlist['crm8008',tp][[1]][intset & abovecut,]
  belowmat=cagecountmatlist['crm8008',tp][[1]][intset & !abovecut,]

  abovemat=abovemat[order(rowMeans(abovemat)),]
  belowmat=belowmat[order(rowMeans(belowmat)),]

  abovemat=abovemat[seq(1,nrow(abovemat),length.out=50),]
  belowmat=belowmat[seq(1,nrow(belowmat),length.out=50),]

  m = rbind( abovemat,belowmat)

  rownames(m)<-1:nrow(m)

  m=melt(m)
  m$value[m$value==0]<-1
  #now color are CRMs vs TSS
  colvect=rep(c('Above Cutoff','Below Cutoff'),each=50)
  # setwd('analysis/crm_chrom_analysis4')
   pdf(paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/above_v_belowcut_dotplot.pdf'))

  ggplot(data=m,aes(x=Var1,y=value,color=colvect[Var1]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='Log10 normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
  scale_color_discrete(name='set',labels=unique(colvect))+scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
  ggtitle(paste0('Samples of Cage signal for above and below cutoff CRMs'))
  dev.off()







#show the seperation in terms of the chromatin 
intset=crmgrs$intergenic & !crmgrs$H3K4me3_peak

top=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])>cutoff
bottom=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])<cutoff


poschrom<-chrom.mean.mats[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats[['crm8008']][intset,][bottom,]

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/crm_cutoff_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()



#now take the bottom 10%
intset=crmgrs$intergenic & !crmgrs$H3K4me3_peak

top=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])>cutoff
bottom=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])<quantile(rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,]),0.10)


poschrom<-chrom.mean.mats[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats[['crm8008']][intset,][bottom,]

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/crm_cutoff_lowqunat__boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()





#Now conditioning on various marks


#show the seperation in terms of the chromatin 

#pol
intset=crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$polII

top=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])>cutoff
bottom=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])<cutoff


poschrom<-chrom.mean.mats[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats[['crm8008']][intset,][bottom,]

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/crm_cutoff_polpeak_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()



allseqlist<- getSeq(Dmelanogaster,crmgrs@seqnames,start(crmgrs),end(crmgrs),as.character=T)
names(allseqlist)<-crmgrs$name


#Now let's output the sequences of those guys
seqlist<-allseqlist[intset & top]
outfile<-paste("pol_crms_high_eRNA.fasta",sep="")
message(outfile)
write.table(paste(">", names(seqlist), "\n", as.character(seqlist), sep=""), file=outfile, row.names=F, col.names=F, sep = "\n", quote=F)

#Now let's output the sequences of those guys
seqlist<-allseqlist[intset & bottom]
outfile<-paste("pol_crms_low_eRNA.fasta",sep="")
message(outfile)
write.table(paste(">", names(seqlist), "\n", as.character(seqlist), sep=""), file=outfile, row.names=F, col.names=F, sep = "\n", quote=F)







#ac
intset=crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$H3K27ac

top=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])>cutoff
bottom=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])<cutoff


poschrom<-chrom.mean.mats[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats[['crm8008']][intset,][bottom,]

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/crm_cutoff_K27Acpeak_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()




allseqlist<- getSeq(Dmelanogaster,crmgrs@seqnames,start(crmgrs),end(crmgrs),as.character=T)
names(allseqlist)<-crmgrs$name


#Now let's output the sequences of those guys
seqlist<-allseqlist[intset & top]
outfile<-paste("ac_crms_high_eRNA.fasta",sep="")
message(outfile)
write.table(paste(">", names(seqlist), "\n", as.character(seqlist), sep=""), file=outfile, row.names=F, col.names=F, sep = "\n", quote=F)

#Now let's output the sequences of those guys
seqlist<-allseqlist[intset & bottom]
outfile<-paste("ac_crms_low_eRNA.fasta",sep="")
message(outfile)
write.table(paste(">", names(seqlist), "\n", as.character(seqlist), sep=""), file=outfile, row.names=F, col.names=F, sep = "\n", quote=F)





#K4me1
intset=crmgrs$intergenic & !crmgrs$H3K4me3_peak & crmgrs$polII

top=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])>cutoff
bottom=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])<cutoff


poschrom<-chrom.mean.mats[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats[['crm8008']][intset,][bottom,]

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis4/crm_cutoff_polpeak_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()






