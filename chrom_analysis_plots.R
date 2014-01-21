setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
chrom.quantile.cutoff = 0.75 #cutoff for conditioning on chromatin marks


crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak
sum(crm8008.gr$intset)

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_chrom_scatter.pdf',h=14,w=14)
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 2-4h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','68h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 6-8h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','1012h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 10-12h')
dev.off()

dim(chrom.mean.mats$crm8008[crm8008.gr$intset,])
length(rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,]))
############
#Start Plotting
###########
#first do a dotplot comparing levels at CRMs vs TSS at 68 hours
tp=unique(accession.df$timepoint)[[3]]
#CRMs vs TSS.
for(tp in unique(accession.df$timepoint)){
  mag=7#orders of magnitude to plot
  crmorder=order(rowMeans(cagecountmatlist['crm8008',tp][[1]][crm8008.gr$intset,]),decreasing=F)
  tssorder=order(rowMeans(cagecountmatlist['tss',tp][[1]]),decreasing=F)
  crmrows=crmorder[seq(from=1,to=length(crmorder),length.out=50)]
  tssrows=tssorder[seq(from=1,to=length(tssorder),length.out=50)]

  m = rbind( cagecountmatlist['crm8008',tp][[1]][crm8008.gr$intset,][crmrows,],cagecountmatlist['tss',tp][[1]][tssrows,])
  rownames(m)<-1:nrow(m)
  m=melt(m)
  m$value[m$value==0]<-1
  #now color are CRMs vs TSS
  colvect=rep(c('Intergenic CRM','TSS'),each=50)
  pdf(paste0('analysis/crm_chrom_analysis/','crms_v_tss_levels',tp,'.pdf'),w=14)
  print(
    ggplot(data=m,aes(x=Var1,y=value,color=colvect[Var1]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='Normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
    scale_color_discrete(name='set',labels=unique(colvect))+scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
    ggtitle(paste0('Cage Signal at 50 quantiles for CRMs vs TSS at ',tp ))
  )
  dev.off()
}


#Cad positives vs Cad negatives
tp=unique(accession.df$timepoint)[[1]]
mag=3#orders of magnitude to plot
mat=cagecountmatlist['cad',tp][[1]]
pmat = mat[cad3.gr$pos,]
pmat = pmat[order(rowMeans(pmat)),]
nmat = mat[cad3.gr$neg,]
nmat = nmat[order(rowMeans(nmat)),]
mat=rbind(pmat,nmat)
rownames(mat)<-1:nrow(mat)
m=melt(mat)
m$value[m$value==0]<-1
#now color are CRMs vs TSS
colvect=c(rep('cadpos',sum(cad3.gr$pos)),rep('cadneg',sum(cad3.gr$neg)))
# setwd('analysis/crm_chrom_analysis5')
pdf(paste0('analysis/crm_chrom_analysis/cadposvneg.pdf'),w=14)
  ggplot(data=m,aes(x=Var1,y=value,color=colvect[as.numeric(as.character(Var1))]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
  scale_color_discrete(name='set')+
  scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
  ggtitle(paste0('Cage Signal at CAD',tp ))
dev.off()









#first let's just produce plots for various possible positive/negative sets

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
  # setwd('analysis/crm_chrom_analysis5')
   pdf(paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/posneg_barplot_',setn,negn,'.pdf'),w=14)
    qplot(data=m,x=Var1,y=value,color=colvect[Var1])+
    coord_trans(y = "log10")+
    scale_y_continuous(name=' normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
    scale_color_discrete(name='set',labels=unique(colvect))+
    scale_x_discrete(name='CRM',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
    ggtitle(paste0('Cage Signal for ',setn,' vs ',negn ))
  dev.off()

    
    pdf(paste0('analysis/crm_chrom_analysis/crmrocplot',setn,negn,'.pdf'))
    ROCfunc(pos,neg,
            main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,setn),
            sub=paste0('windowsize ',w))
    dev.off()
    
    pdf(paste0('analysis/crm_chrom_analysis/crm_prec_rec_plot',setn,negn,'.pdf'))
    cutoff<-PreRecfunc(pos,neg,
                       main=paste0('ROC curve for Summed Normalized Tags, Set: ' ,setn),
                       sub=paste0('windowsize ',w))
    cutoff
    dev.off()
    
    cutoffs[[negn]][[setn]][['crm8008']]  <-scorefunc(cagecountmatlist[[1]])>cutoff
    cutoffs[[negn]][[setn]][['cad3']]     <-scorefunc(cagecountmatlist[[2]])>cutoff

    #and produce the boxplots for the other chromatin marks
    pdf(paste0('analysis/crm_chrom_analysis/crm_chrom_boxplot',setn,negn,'.pdf'))
    poschrom<- chrom.mean.mats[[1]][cutoffs[[negn]][[setn]][['crm8008']] & crm8008.gr$intergenic,]
    negchrom<-chrom.mean.mats[[1]][!cutoffs[[negn]][[setn]][['crm8008']] & crm8008.gr$intergenic,]
    print(Chrom.boxplots(poschrom , negchrom,tit=paste0('All Intergenic CRMs, Split by cutoff ',cutoff)))
    dev.off()
    
    cat(setn)
    cat(negn)
    cat(' ... ')
  }  
}




#what if we have a list of lists, the bottom lists are going to contain
#the type of scoring used, the cutoff, vectors describing pos and neg for cad and crms
#the descriminances and accuracy.
negn=names(negsetlist)[1]
setn=names(possetlist)[6]
posmat=possetlist[[setn]]#get our data
negmat=negsetlist[[negn]]

#Cad positives vs Cad negatives
tp=unique(accession.df$timepoint)[[1]]
mag=3#orders of magnitude to plot
pmat=possetlist[[setn]]#get our data
nmat=negsetlist[[negn]]
pmat = pmat[order(rowMeans(pmat),dec=T),]
nmat = nmat[order(rowMeans(nmat),dec=T),]
mat=rbind(pmat,nmat)
rownames(mat)<-1:nrow(mat)
m=melt(mat)
m$value[m$value==0]<-1
#now color are CRMs vs TSS
colvect=c(rep('Active Cad3 Crms',sum(cad3.gr$pos)),rep('Negative 8008 Crms',sum(cad3.gr$neg)))
#
pdf(paste0('analysis/crm_chrom_analysis/cutoff_dotplot.pdf'),w=14)
  ggplot(data=m,aes(x=Var1,y=value,color=colvect[as.numeric(as.character(Var1))]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
  scale_color_discrete(name='set')+
  scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
  ggtitle(paste0('Cage Signal at CAD',tp ))

dev.off()





#save cutoff info on crm8008.gr
crm8008.gr$abovecadfullcut <- cutoffs$full.neg$cad3pos[['crm8008']]
cad3.gr$abovecadfullcut <- cutoffs$full.neg$cad3pos[['cad3']]

#hunt out some for browser viewing
crm8008.gr$name<-gsub('^(\\d+).*','\\1',x=crm8008.gr$name)
save(crm8008.gr,file=file.crm8008.gr)

#export for viewing in genome browser
intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak
rs=rowMeans(cagecountmatlist['crm8008',tp][[1]][,])
rstofind=rs
rstofind[!intset]=0
topdudes=order(decreasing=T,rstofind)
topcrms=crm8008.gr[topdudes]
topcrms$name<-paste0('topcrm',1:length(topcrms))
export(topcrms,'/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/topcrms.bed')

#print some relevant numbers
sum(crm8008.gr$abovecadfullcut & crm8008.gr$intset)
sum(crm8008.gr$abovecadfullcut & crm8008.gr$intset & crm8008.gr$polII)
sum(crm8008.gr$abovecadfullcut & crm8008.gr$intset & crm8008.gr$H3K27ac)
sum(crm8008.gr$abovecadfullcut & crm8008.gr$intset& crm8008.gr$H3K27ac & crm8008.gr$polII)


#
cad3.gr$H3K4me3_peak <- overlapsAny(cad3.gr,list(chrompeaks.modencode[['K4me3_4-8h']], chrompeaks[['K4me3_4-6h']] ,chrompeaks[['K4me3_6-8h']]),maxgap=50) 
sum(cad3.gr$abovecadfullcut & cad3.gr$intergenic & !cad3.gr$H3K4me3_peak )




#we can now do the boxplots for our quantiles based on the normalized cage data
#for 8008 crms
intset=crm8008.gr$intergenic& !crm8008.gr$H3K4me3_peak
rowranks=rank(rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,]))
#define two logical vectors for the postivie and negative
top=rowranks>=quantile(rowranks,0.95)
bottom=rowranks<=quantile(rowranks,0.05)
poschrom<-chrom.mean.mats[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats[['crm8008']][intset,][bottom,]

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_5quant_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()

poschrom<-chrom.mean.mats.modencode[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats.modencode[['crm8008']][intset,][bottom,]

#with modencode data
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_5quant_modchrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()








  #show the difference between above cutoff and below cutoff guys in a dotplot
  tp='68h'
  mag=4#orders of magnitude to plot
  abovecut=crm8008.gr$abovecadfullcut

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
  # setwd('analysis/crm_chrom_analysis5')
   pdf(paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/above_v_belowcut_dotplot.pdf'),w=14)

  ggplot(data=m,aes(x=Var1,y=value,color=colvect[Var1]))+geom_point()+coord_trans(y = "log10")+scale_y_continuous(name='Log10 normalized Cage Tags',limits=c(1,10^mag),breaks=c(10^(0:mag)))+
  scale_color_discrete(name='set',labels=unique(colvect))+scale_x_discrete(name='Region',labels='')+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=1)+
  ggtitle(paste0('Samples of Cage signal for above and below cutoff CRMs'))
  dev.off()











#show the seperation in terms of the chromatin 
intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak
top=crm8008.gr$abovecadfullcut
bottom=!crm8008.gr$abovecadfullcut
poschrom<-chrom.mean.mats[['crm8008']][intset&top,]
negchrom<-chrom.mean.mats[['crm8008']][intset&bottom,]
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_cutoff_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()



#now take the bottom 10%
intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak
top=crm8008.gr$abovecadfullcut
bottom=rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,])<quantile(rowMeans(cagecountmatlist['crm8008','68h'][[1]][intset,]),0.10)
poschrom<-chrom.mean.mats[['crm8008']][intset,][top,]
negchrom<-chrom.mean.mats[['crm8008']][intset,][bottom,]
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_cutoff_lowqunat__boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()








#Now conditioning on various marks
#show the seperation in terms of the chromatin 


#pol
crm8008.gr$intset = crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$polII
sum(crm8008.gr$intset)

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_chrom_scatter_polpeak.pdf',w=14,h=14)
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 2-4h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','68h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 6-8h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','1012h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 10-12h')
dev.off()

top=crm8008.gr$abovecadfullcut
bottom=!crm8008.gr$abovecadfullcut
poschrom <-chrom.mean.mats[['crm8008']][crm8008.gr$intset&top,]
negchrom <-chrom.mean.mats[['crm8008']][crm8008.gr$intset&bottom,]


pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_cutoff_polpeak_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()

#now with continous values instead of MACs peak
crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak
 crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & chrom.mean.mats$crm8008[,"PolII_6.8"]>quantile(chrom.mean.mats$crm8008[crm8008.gr$intset,"PolII_6.8"],chrom.quantile.cutoff)

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_chrom_scatter_polcut.pdf',w=14,h=14)
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 2-4h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','68h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 6-8h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','1012h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 10-12h')
dev.off()

top=crm8008.gr$abovecadfullcut
bottom=!crm8008.gr$abovecadfullcut
poschrom<-chrom.mean.mats[['crm8008']][crm8008.gr$intset&top,]
negchrom<-chrom.mean.mats[['crm8008']][crm8008.gr$intset&bottom,]
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_cutoff_polcut_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()


allseqlist<- getSeq(Dmelanogaster,crm8008.gr@seqnames,start(crm8008.gr),end(crm8008.gr),as.character=T)
names(allseqlist)<-crm8008.gr$name
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





dev.off()

#ac
crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$H3K27ac
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_chrom_scatter_acpeak.pdf',w=14,h=14)
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 2-4h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','68h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 6-8h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','1012h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 10-12h')
dev.off()

top=crm8008.gr$abovecadfullcut
bottom=!crm8008.gr$abovecadfullcut
poschrom<-chrom.mean.mats[['crm8008']][crm8008.gr$intset&top,]
negchrom<-chrom.mean.mats[['crm8008']][crm8008.gr$intset&bottom,]
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_cutoff_K27Acpeak_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()


#now with continous values instead of MACs peak
crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak
crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & chrom.mean.mats$crm8008[,"H3K27ac_6.8"]>quantile(chrom.mean.mats$crm8008[crm8008.gr$intset,"H3K27ac_6.8"],chrom.quantile.cutoff)
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_chrom_scatter_ac_cutoff.pdf',w=14,h=14)
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 2-4h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','68h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 6-8h')
chromscatters(chrom.mean.mats$crm8008[crm8008.gr$intset,],rowMeans(cagecountmatlist['crm8008','1012h'][[1]][crm8008.gr$intset,]),tit='chromatin marks vs eRNA at 10-12h')
dev.off()
top=crm8008.gr$abovecadfullcut
bottom=!crm8008.gr$abovecadfullcut
poschrom<-chrom.mean.mats[['crm8008']][crm8008.gr$intset&top,]
negchrom<-chrom.mean.mats[['crm8008']][crm8008.gr$intset&bottom,]
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_cutoff_K27Ac_cutoff_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()


allseqlist<- getSeq(Dmelanogaster,crm8008.gr@seqnames,start(crm8008.gr),end(crm8008.gr),as.character=T)
names(allseqlist)<-crm8008.gr$name
#Now let's output the sequences of those guys
seqlist<-allseqlist[crm8008.gr$intset & top]
outfile<-paste("ac_crms_high_eRNA.fasta",sep="")
message(outfile)
write.table(paste(">", names(seqlist), "\n", as.character(seqlist), sep=""), file=outfile, row.names=F, col.names=F, sep = "\n", quote=F)
#Now let's output the sequences of those guys
seqlist<-allseqlist[crm8008.gr$intset & bottom]
outfile<-paste("ac_crms_low_eRNA.fasta",sep="")
message(outfile)
write.table(paste(">", names(seqlist), "\n", as.character(seqlist), sep=""), file=outfile, row.names=F, col.names=F, sep = "\n", quote=F)




#K4me1
crm8008.gr$intset=crm8008.gr$intergenic & !crm8008.gr$H3K4me3_peak & crm8008.gr$polII
top=crm8008.gr$abovecadfullcut
bottom=!crm8008.gr$abovecadfullcut
poschrom<-chrom.mean.mats[['crm8008']][crm8008.gr$intset&top,]
negchrom<-chrom.mean.mats[['crm8008']][crm8008.gr$intset&bottom,]

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/crm_chrom_analysis/crm_cutoff_polpeak_mesochrom_boxplots.pdf')
print(Chrom.boxplots(poschrom , negchrom))
dev.off()




chrommat=chrom.mean.mats$crm8008[crm8008.gr$intset,]
eRNA24=rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,])
eRNA68=rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,])
eRNA1012=rowMeans(cagecountmatlist['crm8008','24h'][[1]][crm8008.gr$intset,])

lm.ac.pol<-lm(chrommat[,'H3K27ac_6.8']~chrommat[,'PolII_6.8'])
lm.pol.ac<-lm(chrommat[,'PolII_6.8' ]~chrommat[,'H3K27ac_6.8'])
cor(eRNA68,lm.ac.pol$residuals)
cor(eRNA68,lm.pol.ac$residuals)


