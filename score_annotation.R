#This script integrates the various datasets on genes, crms, chroamtin, rnaseq, cage, tagseq and dnase
#it outputs a list of GRange objects with the data stored as vectors and matrices in the meta data columns
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
source('src/load_chromdata.R')
load('data/objects/gene.annotation.object.R')
load('data/objects/crm.annotation.object.R')
load('data/objects/cg.object.R')
load('data/objects/cg.mapfilt.object.R')
load('data/objects/cg.pl.object.R')
load('data/objects/cg.mapfilt.pl.object.R')
load('data/objects/ts.object.R')
load('data/objects/ts.pl.object.R')
load('data/objects/rna.seq.object.R')
library(ROCR,lib.loc='~/Harnett/R/')
w=500#window size for determining CRM scores


# getStrandedExprMat(gr=transcripts.gr[c(2,9,10000,14000)],cg[1:10])


####TODO - turn the above into loops - as below
gene.annotations = c('transcripts.gr','tss.gr','genes.gr')
crm.annotations = c('crm8008.gr','cad3.gr')
#expr.datasets = c('cg','cg.pl','cg.mapfilt','cg.mapfilt.pl','ts','ts.pl','rna.seq','c.rna.seq','z.rna.seq')
expr.datasets = c('cg','cg.pl','cg.mapfilt','cg.mapfilt.pl','ts','ts.pl','rna.seq','c.rna.seq')


file.remove('analysis/expr.mats.hdf5')
h5createFile('analysis/expr.mats.hdf5')
#now iterate over our features and datasets to get the scores
expr.mats=sapply(simplify=F,gene.annotations,function(region){
  h5createGroup(file='analysis/expr.mats.hdf5',group=region)
  # browser()
  region=get(region)[]
  region$id=1:length(region)
	sapply(simplify=F,expr.datasets,function(data){
		suppressWarnings({wideregion = resize(region,width=width(region)+250,fix='center')})
		wideregion=sort(wideregion)
     #  browser()
		m=getStrandedExprMat(wideregion,get(data))
		rownames(m)=wideregion$id
 #   h5write(file='analysis/expr.mats.hdf5',name=paste0(region,'/',data),obj= m[as.character(region$id),])
  #  cat('.')
    m[as.character(region$id),]
  })
})

str(expr.mats)
save(expr.mats,file = 'data/objects/expr.mats.object.R')




#shortened data works okay
expr.crm.mats=sapply(simplify=F,crm.annotations,function(region){
  cat(region)
  region=get(region)
  region$id=1:length(region)
  sapply(simplify=F,expr.datasets,function(data){
    cat(data)
    suppressWarnings({wideregion = resize(region,width=width(region)+0,fix='center')})
    wideregion= sort (wideregion)
    expect_is(wideregion,'GRanges')
    expect_true(all(wideregion==sort(wideregion)))
    cat('.')
      # browser()
    # tryCatch({
      m=getBestWindowMat(reg=wideregion,windowsize=500,chrs.keep=seqnames(si),cage=get(data))
      expect_true(is.numeric(m))
      expect_is(m,'matrix')
      rownames(m)=wideregion$id
      m[as.character(region$id),]
    # },warning = function(w)browser(),error = function(e)browser)
  })
})


###
reg=wideregion
windowsize=500
chrs.keep=seqnames(si)
cage=get(data)[1:2]
##

str(expr.crm.mats)


save(expr.crm.mats,file = 'data/objects/expr.crm.mats.object.R')


########################################################################################
#Now collect 

genes.gr$geneDensity = getGeneDensity(genes.gr,w=200000)


save(transcripts.gr,genes.gr,tss.gr,file='data/objects/scored.annotation.object.R')
load('data/objects/scored.annotation.object.R')



#3) Now 8008 CRMs -----------------------------------------------------------

# Chromatin data ---------------------------------------------------------

overlapsAny<-function(gr,peaklist,maxgap=50){
  0 < rowSums(sapply(peaklist,function(peaks){
    countOverlaps(gr,peaks,maxgap)
  }))
}

center.gr<-function(rle,gr){
  centred.gr<-nomcols(gr)
  cageviews<-Views(alltags$both,as(centred.gr,'RangesList'))[unique(seqnames(centred.gr))] 
  maxpos<- unlist( viewWhichMaxs(cageviews,na.rm=T) )
  start(centred.gr)<-maxpos-250
  centred.gr<-resize(centred.gr,width=500,fix='start')
  centred.gr<-sort(centred.gr)
}

cad.centered.gr<-center.gr(alltags$both,cad3.gr)

#columns describind overlap with peak
crm8008.gr$H3K4me3  <-overlapsAny(crm8008.gr,list(chrompeaks.modencode[['K4me3_4-8h']], chrompeaks[['K4me3_4-6h']] ,chrompeaks[['K4me3_6-8h']]),maxgap=50) 
crm8008.gr$H3K4me1   <-overlapsAny(crm8008.gr,list(chrompeaks.modencode[['K4me1_4-8h']], chrompeaks[['K4me1_4-6h']] ,chrompeaks[['K4me1_6-8h']]),maxgap=50)
crm8008.gr$H3K27ac   <-overlapsAny(crm8008.gr,list(chrompeaks.modencode[['K27ac_4-8h']], chrompeaks[['K27Ac_4-6h']] ,chrompeaks[['K27Ac_6-8h']]),maxgap=50)
crm8008.gr$H3K79me3   <-overlapsAny(crm8008.gr,list(chrompeaks[['K79me3_4-6h']],chrompeaks[['K79me3_6-8h']]),maxgap=50)
crm8008.gr$H3K36me3  <-overlapsAny(crm8008.gr,list(chrompeaks[['K36me3_4-6h']],chrompeaks[['K36me3_6-8h']]),maxgap=50)
crm8008.gr$polII          <-overlapsAny(crm8008.gr,list(chrompeaks[['PolII_4-6h']],chrompeaks[['PolII_6-8h']]),maxgap=50)
crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3

#columns describind overlap with peak
cad3.gr$H3K4me3  <-overlapsAny(cad3.gr,list(chrompeaks.modencode[['K4me3_4-8h']], chrompeaks[['K4me3_4-6h']] ,chrompeaks[['K4me3_6-8h']]),maxgap=50) 
cad3.gr$H3K4me1   <-overlapsAny(cad3.gr,list(chrompeaks.modencode[['K4me1_4-8h']], chrompeaks[['K4me1_4-6h']] ,chrompeaks[['K4me1_6-8h']]),maxgap=50)
cad3.gr$H3K27ac   <-overlapsAny(cad3.gr,list(chrompeaks.modencode[['K27ac_4-8h']], chrompeaks[['K27Ac_4-6h']] ,chrompeaks[['K27Ac_6-8h']]),maxgap=50)
cad3.gr$H3K79me3   <-overlapsAny(cad3.gr,list(chrompeaks[['K79me3_4-6h']],chrompeaks[['K79me3_6-8h']]),maxgap=50)
cad3.gr$H3K36me3  <-overlapsAny(cad3.gr,list(chrompeaks[['K36me3_4-6h']],chrompeaks[['K36me3_6-8h']]),maxgap=50)
cad3.gr$polII          <-overlapsAny(cad3.gr,list(chrompeaks[['PolII_4-6h']],chrompeaks[['PolII_6-8h']]),maxgap=50)
cad3.gr$intset=cad3.gr$intergenic & !cad3.gr$H3K4me3
cad3.gr$pos <- cad3.gr$active68 %in% T & cad3.gr$intergenic 
cad3.gr$neg <- cad3.gr$inactive68 %in% T & cad3.gr$intergenic & !cad3.gr$GJ & !cad3.gr$polII


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
save(tfgrlist,file='data/objects/tfgrlist.object.R')
#now let's go through our 6-8 hour set and add the overlap info for the tfs at various timepoints
#for our tf,tp, do countoverlaps four our crm8008.gr
tmp<-sapply(simplify=F,as.character(unique(tffile.df$tp)),function(tp){
  sapply(simplify=F,as.character(unique(tffile.df$tfname)),function(tf){
    if(is.null(tfgrlist[[tp]][[tf]])){return(NULL)}   
    mcols(crm8008.gr)[[paste(tf,tp,sep='_')]]<<-countOverlaps(crm8008.gr,tfgrlist[[tp]][[tf]])>0 
    NULL
  })
}) 
rm(tmp)
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
#and export bedfiles
export(crm8008.gr[ crm8008.gr$in.tf.pos],con ='analysis/make_regions_bedfiles/pos.tf.8008crms.bed')
export(crm8008.gr[ crm8008.gr$in.tf.neg],con='analysis/make_regions_bedfiles/neg.tf.8008crms.bed')


save(cad3.gr,file=file.cad3)


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

#Get matrices of chromatin overlap data
if(file.exists('data/objects/chrom.mean.mats.object.R')){
  load('data/objects/chrom.mean.mats.object.R')
  message('data/objects/chrom.mean.mats.object.R already present, loaded')
}else{
  chrom.rles.rpgc.sub.merge=chrom.rles.rpgc.sub.merge[!grepl('K36me3',names(chrom.rles.rpgc.sub.merge))]
  chrom.rles.rpgc.sub.merge=c(chrom.rles.rpgc.sub.merge,'dnase_6.8'=(dnase.rles$STG10+dnase.rles$STG11))
  chrom.mean.mats<-list(
    crm8008=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crm8008.gr,'RangesList'))))}),
    cad=sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
  )
  chrom.mean.mats$crm8008 <- cbind(chrom.mean.mats$crm8008,data.frame(NearestGeneExpr=GetNearestExpression(crm8008.gr)))
  chrom.mean.mats$cad <- cbind(chrom.mean.mats$cad,data.frame(NearestGeneExpr=GetNearestExpression(cad3.gr)))
  save(chrom.mean.mats,file='data/objects/chrom.mean.mats.object.R')
}  #read chromati
#Now for the modencode data
if(file.exists('data/objects/chrom.mean.mats.modencode.object.R')){
  load('data/objects/chrom.mean.mats.modencode.object.R')
  message('data/objects/chrom.mean.mats.modencode.object.R already present, loaded')
}else{
  chrom.mean.mats.modencode<-list(
    crm8008=sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crm8008.gr,'RangesList'))))}),
    cad=sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(cad.centered.gr,'RangesList'))))})
  )
  chrom.rles.rpgc.sub.merge=chrom.rles.rpgc.sub.merge[!grepl('K36me3',names(chrom.rles.rpgc.sub.merge))]
  chrom.rles.rpgc.sub.merge=c(chrom.rles.rpgc.sub.merge,'dnase_6.8'=(dnase.rles$STG10+dnase.rles$STG11))
  chrom.mean.mats.modencode$crm8008 <- cbind(chrom.mean.mats.modencode$crm8008,data.frame(NearestGeneExpr=GetNearestExpression(crm8008.gr)))
  chrom.mean.mats.modencode$cad <- cbind(chrom.mean.mats.modencode$cad,data.frame(NearestGeneExpr=GetNearestExpression(cad3.gr)))
  save(chrom.mean.mats.modencode,file='data/objects/chrom.mean.mats.modencode.object.R')
}  #read chromati



load('cagecountmatlist.object.R')


simplify2array(expr.crm.mats$cg.pl)

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

negn='full.neg'
posn='cad3pos'
prec.cutoff=0.95
cutoffs=list()
# Now go through each of our sets and produce our plots -------------------
for(negn   in names(negsetlist)){
  for(posn in names(possetlist)){
    #dir.create('analysis/crm_chrom_analysis2/')
    mag=3
    
    posmat=possetlist[[posn]]#get our data
    negmat=negsetlist[[negn]]
    pos=scorefunc(posmat)#for now just sum across lines
    neg=scorefunc(negmat)
    

    m=rbind(posmat,negmat)
    rownames(m)<-1:nrow(m)
    colvect=c(rep(posn,nrow(posmat)),rep(negn,nrow(negmat)))
    m=melt(m)
    m$value[m$value==0]<-1
  # setwd('analysis/crm_chrom_analysis5')
   pdf(paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/posneg_barplot_',posn,negn,'.pdf'),w=14)
    qplot(data=m,x=Var1,y=value,color=colvect[Var1])+
    coord_trans(y = "log10")+
    scale_y_continuous(name=' normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
    scale_color_discrete(name='set',labels=unique(colvect))+
    scale_x_discrete(name='CRM',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
    ggtitle(paste0('Cage Signal for ',posn,' vs ',negn ))
  dev.off()

    
    pdf(paste0('analysis/crm_chrom_analysis/crmrocplot',posn,negn,'.pdf'))
    ROCfunc(pos,neg,
            main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,posn),
            sub=paste0('windowsize ',w))
    dev.off()
    
    pdf(paste0('analysis/crm_chrom_analysis/crm_prec_rec_plot',posn,negn,'.pdf'))
    cutoff<-PreRecfunc(pos,neg,
                       main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,posn),
                       sub=paste0('windowsize ',w))
    cutoff
    dev.off()
    
    cutoffs[[negn]][[posn]][['crm8008']]  <-scorefunc(cagecountmatlist[[1]])>cutoff
    cutoffs[[negn]][[posn]][['cad3']]     <-scorefunc(cagecountmatlist[[2]])>cutoff

    #and produce the boxplots for the other chromatin marks
    pdf(paste0('analysis/crm_chrom_analysis/crm_chrom_boxplot',posn,negn,'.pdf'))
    poschrom<- chrom.mean.mats[[1]][cutoffs[[negn]][[posn]][['crm8008']] & crm8008.gr$intergenic,]
    negchrom<-chrom.mean.mats[[1]][!cutoffs[[negn]][[posn]][['crm8008']] & crm8008.gr$intergenic,]
    print(Chrom.boxplots(poschrom , negchrom,tit=paste0('All Intergenic CRMs, Split by cutoff ',cutoff)))
    dev.off()
    
    cat(posn)
    cat(negn)
    cat(' ... ')
  }  
}

#save cutoff info on crm8008.gr
crm8008.gr$abovecadfullcut <- cutoffs$full.neg$cad3pos[['crm8008']]
cad3.gr$abovecadfullcut <- cutoffs$full.neg$cad3pos[['cad3']]
crm8008.gr$abovecadfullcut

#So now we can integrate all our information to pick the ones which are above the cutoff, have low
crm8008.gr$cagesum.68 = rowSums(cagecountmatlist[['crm8008','68h']])
crm8008.gr$startmotif.overlap = motif.overlap$crm8008.gr
crm8008.gr$any.promotor.motifs = apply(crm8008.gr$startmotif.overlap,1,any)
for_testing.gr = crm8008.gr[ crm8008.gr$abovecadfullcut & crm8008.gr$in.tf.pos ]
for_testing.gr = for_testing.gr[ order(decreasing=T,for_testing.gr$cagesum.68)]
for_testing.gr$startmotif.overlap
for_testing.gr$name = paste0('ex',1:length(for_testing.gr),'_eRNA_crm')
export(for_testing.gr,'olga_eRNA_test_crms.bed')
write.csv(as.data.frame(mcols(for_testing.gr)),row.names=F,quote=F,file='olga_eRNA_test_crms.annotation.csv')



#So now we can integrate all our information to pick the ones which are above the cutoff, have low
cad3.gr$cagesum.68 = rowSums(cagecountmatlist[['cad','68h']])
cad3.gr$startmotif.overlap = motif.overlap$cad3.gr
cad3.gr$any.promotor.motifs = apply(cad3.gr$startmotif.overlap,1,any)
cadfor_testing.gr = cad3.gr[ cad3.gr$abovecadfullcut & cad3.gr$intergenic & !cad3.gr$H3K4me3 ] 
cadfor_testing.gr = cadfor_testing.gr[ order(decreasing=T,cadfor_testing.gr$cagesum.68)]
cadfor_testing.gr$startmotif.overlap
cadfor_testing.gr$name = paste0('ex',1:length(cadfor_testing.gr),'_cad_eRNA_crm')
export(cadfor_testing.gr,'olga_eRNA_test_cad_crms.bed')
tab=mcols(cadfor_testing.gr)
tab[,2] = gsub(' ','_',x=tab[,2])
write.table(as.data.frame(tab),sep='\t',row.names=F,col.names=T,quote=F,file='olga_eRNA_test_cad_crms.annotation.txt')
read.table('olga_eRNA_test_cad_crms.annotation.txt',sep='\t')
write.(as.data.frame(tab),sep='\t',row.names=F,col.names=F,quote=F,file='olga_eRNA_test_cad_crms.annotation.csv')
read.table('olga_eRNA_test_cad_crms.annotation.csv',sep='\t')


#Now let's pick some negatives for the 8008
crm8008.gr$startmotif.overlap = motif.overlap$crm8008.gr
crm8008.gr$any.promotor.motifs = apply(crm8008.gr$startmotif.overlap,1,any)
neg_examples_8008.gr = crm8008.gr[ !crm8008.gr$abovecadfullcut & crm8008.gr$in.tf.pos & crm8008.gr$H3K27ac  ]
neg_examples_8008.gr = neg_examples_8008.gr[ order(decreasing=F,neg_examples_8008.gr$cagesum.68)]
neg_examples_8008.gr$startmotif.overlap
neg_examples_8008.gr$name = paste0('anex',1:length(neg_examples_8008.gr),'_eRNA_crm')
export(neg_examples_8008.gr,'olga_neg_eRNA_test_crms.bed')
write.csv(as.data.frame(mcols(neg_examples_8008.gr)),row.names=F,quote=F,file='olga_neg_eRNA_test_crms.annotation.csv')

#and negatives for the CAD
cad3.gr$startmotif.overlap = motif.overlap$cad3.gr
cad3.gr$any.promotor.motifs = apply(cad3.gr$startmotif.overlap,1,any)
neg_examples_cad3.gr = cad3.gr[ !cad3.gr$abovecadfullcut & cad3.gr$active68 & cad3.gr$active24 & cad3.gr$intergenic ]
neg_examples_cad3.gr = cad3.gr[ !cad3.gr$abovecadfullcut & cad3.gr$active68 & cad3.gr$intergenic ]

neg_examples_cad3.gr = neg_examples_cad3.gr[ order(decreasing=F,neg_examples_cad3.gr$cagesum.68)]
neg_examples_cad3.gr$startmotif.overlap
neg_examples_cad3.gr$name = paste0('znex',1:length(neg_examples_cad3.gr),'_eRNA_cd')
export(neg_examples_cad3.gr,'olga_neg_eRNA_test_cad.bed')
write.csv(as.data.frame(mcols(neg_examples_cad3.gr)),row.names=F,quote=F,file='olga_neg_eRNA_test_cad.annotation.csv')

