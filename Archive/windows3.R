crm(list=ls())

source('src/tss_cage_functions.R')
library(reshape2)
library(ggplot2)
library(ROCR,lib.loc='~/Harnett/R')

rpgc.step.num<-100
pval.step.num<-20
prec.cutoffs<-c(0.8,0.85,0.9,0.95)

load(file.cage.tag.rles)
load(file.accession.df)
regs=list(
  pos=import(con='analysis/positive_8008.bed',asRangedData=F,seqinfo=si),
  neg=import(con='analysis/negative_8008.bed',asRangedData=F,seqinfo=si),
  int=import(con='analysis/random.intergenic.bed',asRangedData=F,seqinfo=si),
  cad3pos=import(con='analysis/cad3.pos.bed',asRangedData=F,seqinfo=si),
  cad3neg=import(con='analysis/cad3.neg.bed',asRangedData=F,seqinfo=si)
)
#define window sizes
window.sizes<-c('w25'=25,'w50'=50,'w100'=100,'w200'=200,'w300'=300,'w500'=500)
window.sizes<-window.sizes[6]
#exclude the reseq libraries.
cg<-cage.tag.rles[!grepl('reseq',names(cage.tag.rles))]

#rm(cage.tag.rles)

cagecountmatlist<-sapply(simplify=F,window.sizes,function(w){
  sapply(simplify=F,regs,function(reg){   
    # sapply(simplify=F,c(10,100,500),function(w){
    #make sure all our regions are at least the window size
    reg<-resize(reg,width=pmax(width(reg),w),fix='center')
    #make sure th
    #convert them to a Rangelist for Views
    winds<-as(reg,'RangesList')[chrs.keep]
    winds<-winds[bigchrs]
    starts<-unlist(start(winds))
    chrs<-names(starts) 
    #for each line  
    maxwinds<-mclapply(mc.cores=10,cg,function(acc){
      #now get the windowed views for each window
      v<-unlist(viewApply(Views(acc[['pos']][bigchrs],winds[bigchrs]),FUN=function(x){ (runsum(x,w)) }))
      v.n<-unlist(viewApply(Views(acc[['neg']][bigchrs],winds[bigchrs]),FUN=function(x){ (runsum(x,w)) }))
      #get the bigger of the two strands
      ispos<-(sapply(v,max)>sapply(v.n,max))
      v<-ifelse( ispos,v,v.n)
      sapply(v,max)
    })
    simplify2array(maxwinds)
  })
})

#
names(cagecountmatlist)<-window.sizes
#add a column of normalized values
v<-accession.df$intergenic.lib.size
names(v)<-accession.df$accession
#windowsize should be a numeric variable

save(cagecountmatlist,file='data/objects/cagecountmatlist.object.R')

cagecountmatlist.r<-sapply(simplify=F,as.character(window.sizes),function(w){
  sapply(simplify=F,names(regs),function(reg){
    sapply(1:ncol(cagecountmatlist[[w]][[reg]]),function(n){
      cagecountmatlist[[w]][[reg]][,n]<-cagecountmatlist[[w]][[reg]][,n]/v[names(cg)[n]]
      })
  })
})
 



#create matrices for recording parameter thresholds
rpgc.wind.matrix<-matrix(NA,nrow=rpgc.step.num,ncol=length(window.sizes))
pval.wind.matrix<-matrix(NA,nrow=pval.step.num,ncol=length(window.sizes))
wind.pval.prec<-array(dim=c(length(window.sizes),pval.step.num,length(prec.cutoffs)))
#window sizes as columns
colnames(rpgc.wind.matrix)<-as.character( window.sizes)
colnames(pval.wind.matrix)<-as.character( window.sizes)

#duplicate these matrices for accuracies instead of aucs
rpgc.wind.matrix.acc<-rpgc.wind.matrix
pval.wind.matrix.acc<-pval.wind.matrix





message('Done with data matrices')

# first kind of fdr - total reads -----------------------------------------
for (w in as.character(window.sizes)){
  
  #select relevant matrix
  cagecountmats<-cagecountmatlist[[w]]
  cagecountmats.r<-cagecountmatlist.r[[w]]
  
  stopifnot(!NA%in%cagecountmats.r )
  #calculate the summed 
  sums<-sapply(simplify=F,cagecountmats,function(reg){rowSums(reg)})#total in all lines
  sums.r<-sapply(simplify=F,cagecountmats.r,function(reg){rowSums(reg)})#normalized
  
  pdf(paste0('analysis/summ.roc.',w,'.pdf'))
  
  #8008 positive plus negative set
  #par(mfrow=c(2,2))
  scores<-stack(list(pos=sums[[1]],neg=sums[[2]]))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T)
  aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
  text(x=0.90,y=0,labels=aucstring)
  title('ROC curve for Summed Unnormalized Tags,8008 set',sub=paste0('windowsize ',w))
  # plot(performance(sumpred,'prec','rec'),colorize=T)
  # plot(performance(sumpred,'cal',window.size=5))
  # performance(sumpred,'auc')@y.values[[1]]
  abline(coef=c(0,1),lty=3)
  
  #8008 positive vs the intergenic set
  #par(mfrow=c(2,2))
  scores<-stack(list(pos=sums[[1]],neg=sums[[3]]))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T)
  aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
  text(x=0.90,y=0,labels=aucstring)
  title('ROC curve for Summed Unnormalized Tags,crms vs random inactive',sub=paste0('windowsize ',w))
  # plot(performance(sumpred,'prec','rec'),colorize=T)
  # plot(performance(sumpred,'cal',window.size=5))
  # performance(sumpred,'auc')@y.values[[1]]
  abline(coef=c(0,1),lty=3)
  
  #crm set with low H3k27ac
  scores<-stack(list(pos=sums[[1]],neg=sums[[2]][ regs[[2]]$score==1 ]))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T)
  aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
  text(x=0.90,y=0,labels=aucstring)
  title('ROC curve for Summed Tags,8008 low k27ac negative set',sub=paste0('windowsize ',w))
  abline(coef=c(0,1),lty=3)
  
  #cad3
  scores<-stack(list(pos=sums[[4]],neg=sums[[5]]))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T)
  aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
  text(x=0.90,y=0,labels=aucstring)
  title('ROC curve for Summed Tags, CAD3 set',sub=paste0('windowsize ',w))
  abline(coef=c(0,1),lty=3)
  
  #crms with normalized reads
  scores<-stack(list(pos=sums.r[[1]],neg=sums.r[[2]]))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T)
  aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
  text(x=0.90,y=0,labels=aucstring)
  title('ROC curve for Summed RPGC Tags,8008 set',sub=paste0('windowsize ',w))
  abline(coef=c(0,1),lty=3)
  
  #crms vs intergenic 
  scores<-stack(list(pos=sums.r[[1]],neg=sums.r[[3]]))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T)
  aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
  text(x=0.90,y=0,labels=aucstring)
  title('ROC curve for Summed RPGC Tags,crms vs random inactive',sub=paste0('windowsize ',w))
  abline(coef=c(0,1),lty=3)
  
  #cad
  scores<-stack(list(pos=sums.r[[4]],neg=sums.r[[5]]))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T)
  aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
  text(x=0.90,y=0,labels=aucstring)
  title('ROC curve for Summed Normalized Tags,CAD3 set',sub=paste0('windowsize ',w))
  abline(coef=c(0,1),lty=3)
  
  dev.off()
  
  #let's also try using our 'very negative' set
  
  
  
  
}

message('Done with Summed tag plots')


# Second type of FDR - lines above flat cutoff ----------------------------

message('Done with flat cutoff plots')


for (w in as.character(window.sizes)){
  cagecountmats<-cagecountmatlist[[w]]
  cagecountmats.r<-cagecountmatlist.r[[w]]
  
  
  #define a range to cut off our peaks at
  
  
  highnum<-20#arbitrarily set a hight limit on our rpgc cutoff
  all<-do.call(rbind,cagecountmats)#get the scores in all regions
  all.r<-do.call(rbind,cagecountmats.r)#get the scores in all regions
  highnum<-all.r[which(all[,1]>20)[1]]
  lowvals<-sort((unique(as.vector(all.r[,1]))))#for only the biggest library
  lowvals<-lowvals[lowvals < highnum]
  lowvals<-lowvals[seq(from=1,to=length(lowvals),length.out=rpgc.step.num)]
  lnrange<-1:ncol(cagecountmats.r[[1]])
  
  
  #first go through our cutoffs and give each crm a linescore
  cutoff.linenum<-sapply(simplify=F,as.character((lowvals)),function(val){#now pick an RPGC value
    #now go through our normalized libraries 
    lnmats<-sapply(simplify=F,cagecountmats.r,function(reg){
      tmp<-apply(reg,1,function(r){#for each locations
        sum(r>as.numeric(val))#count the lines with tagnum over our cutoff
      })
    })  
  })
  
  
  #we now have a matrix with four rows and columns for each cutoff - this can give us an fdr with confidence
  
  #this function takes in a matrix with four columns true positive, false positive, total true, total false.
  #it outputs three columns - an fdr, and the 95% jeffreys credible intervals for these.
  
  aucs<-sapply(names(cutoff.linenum),function(cutoff){
    scores<-stack(list(pos=cutoff.linenum[[cutoff]][[1]],neg=cutoff.linenum[[cutoff]][[2]]))
    sumpred<-prediction(scores$values,scores$ind)
    performance(sumpred,'auc')@y.values[[1]]
  })
  
  lowcuts<-names(cutoff.linenum)
  rownames(rpgc.wind.matrix)<-lowcuts
  rpgc.wind.matrix[,as.character(w)]<-aucs
  
  
  #now do t
  
  #logical vector describing if the cutoff gives at least a nonzero for each region.
  cutlog<-sapply(lowcuts,function(bestcut){ all( sapply(cutoff.linenum[[bestcut]][c(1,2,4,5)],function(reg){any(reg>0)}))})
  #only show cutoffs with some non zeros   
  cuts<-lowcuts[cutlog]
  #we need to exlclude 
  cuts<-cuts[seq(length(cuts)/3,length(cuts),length.out=3)]
  best<-lowcuts[which.max(aucs)]
  
  #let's do four different cutoffs
  pdf(paste0('analysis/rpgc.roc.',w,'.pdf'))
  for (bestcut in c(best,cuts)){
    
    #crm set
    cutlist<-    cutoff.linenum[[bestcut]]
    scores<-stack(list(pos=cutlist[[1]],neg=cutlist[[2]]))
    lcpred<-prediction(scores$values,scores$ind)
    perf<-performance(lcpred,'tpr','fpr')
    plot(perf,colorize=T)
    title(paste0('ROC curve for Number of lines above cutoff ',round(as.numeric(bestcut),4),'  Tags,crm set'),sub=paste0('windowsize ',w))
    abline(coef=c(0,1),lty=3)
    aucstring=paste0('AUC = ',round(performance(lcpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    
    plot(performance(lcpred,'prec','rec'),colorize=T)
    title(paste0('Precision/Recall curve for Number of lines above cutoff\n',round(as.numeric(bestcut),4),'  Tags,crm set'),sub=paste0('windowsize ',w))
    
    
    #crm set with low H3k27ac
    scores<-stack(list(pos=cutlist[[1]],neg=cutlist[[2]][ regs[[2]]$score==0 ]))
    sumpred<-prediction(scores$values,scores$ind)
    perf<-performance(sumpred,'tpr','fpr')
    plot(perf,colorize=T)
    aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    title('ROC curve for Number of lines above cutoff,8008 low k27ac negative set',sub=paste0('windowsize ',w))
    abline(coef=c(0,1),lty=3)
    
    
    #crms vs the randoms
    cutlist<-    cutoff.linenum[[bestcut]]
    scores<-stack(list(pos=cutlist[[1]],neg=cutlist[[3]]))
    lcpred<-prediction(scores$values,scores$ind)
    perf<-performance(lcpred,'tpr','fpr')
    plot(perf,colorize=T)
    title(paste0('ROC curve for Number of lines above cutoff ',round(as.numeric(bestcut),4),'  Tags,crms vs random'),sub=paste0('windowsize ',w))
    abline(coef=c(0,1),lty=3)
    aucstring=paste0('AUC = ',round(performance(lcpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    
    plot(performance(lcpred,'prec','rec'),colorize=T)
    title(paste0('Precision/Recall curve for Number of lines above cutoff\n',round(as.numeric(bestcut),4),'  Tags,crms vs random '),sub=paste0('windowsize ',w))
    
    
    
    
    #cad 
    scores<-stack(list(pos=cutlist[[4]],neg=cutlist[[5]]))
    lcpred<-prediction(scores$values,scores$ind)
    perf<-performance(lcpred,'tpr','fpr')
    plot(perf,colorize=T)
    title(paste0('ROC curve for Number of lines above cutoff ',round(as.numeric(bestcut),4),'  Tags,CAD set'),sub=paste0('windowsize ',w))
    abline(coef=c(0,1),lty=3)
    aucstring=paste0('AUC = ',round(performance(lcpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    
    plot(performance(lcpred,'prec','rec'),colorize=T)
    title(paste0('Precision/Recall curve for Number of lines above cutoff\n',bestcut,'  Tags,CAD set'),sub=paste0('windowsize ',w))
  }
  dev.off()
}

message('Done with flat cutoff plots')





#first calculate the 

# lambdas for third kind of fdr -------------------------------------------


inter=import('analysis/intergenic.inactive.regions.bed')
inter=coverage(inter)>0
#first we'll define lambda for each line
#we can do this easily enough by just defining a 
if(file.exists(file.lambdas)){load(file.lambdas)
}else{
  lambdas<-sapply(names(cg),function(acc){
    zeros<-sum(as.numeric(sum(cg[[acc]][[1]][inter]==0))+  as.numeric(sum(cg[[acc]][[1]][inter]==0)))
    ones<-sum(as.numeric(sum(cg[[acc]][[1]][inter]==1))+  as.numeric(sum(cg[[acc]][[2]][inter]==1)))
    ones/(zeros+ones)#The poisson's density will be mostly concentrated on zero and 1
  })
  save(lambdas,file=file.lambdas)
}

# v=accession.df$library.size
# names(v)<-accession.df$accession
# qplot(v[names(lambdas)],lambdas,xlab='Library Size',main='Calculated Poisson Rate Vs. Library Size')
# v=accession.df$intergenic.lib.size
# names(v)<-accession.df$accession
# qplot(v[names(lambdas)],lambdas,xlab='Intergenic Library Size',main='Calculated Poisson Rate Vs. Intergenic Library Size')




# Have accidentally deleted a section here - the part for calculat --------






# save stuff --------------------------------------------------------------



write.table(rpgc.wind.matrix,'analysis/rpgc.window.mat.txt')

write.table(pval.wind.matrix,'analysis/pval.window.mat.txt')



read.table( 'analysis/')

save.image('tmp.image.R')
jpeg('analysis/AUC_over_windows_pvalues.jpeg')
# x=pval.wind.matrix
# image(t(x[nrow(x):1, ]),xlab='WindowSize', ylab='-ve Log10 Pvalue',main='AUC  for 8008 positive vs. negative set over \ndifferent parameter values',
#       breaks = quantile(x, probs = seq(0, 1, 0.01)),
#       col = colorRampPalette(c('white', 'orange'))(length(quantile(x, probs = seq(0, 1, 0.01)) ) - 1L),
#       xaxt = 'none',
#       yaxt = 'none')
# axis(1,at=seq(0,1,length.out=ncol(x)),labels=colnames(x))
# axis(2,at=seq(0,1,length.out=nrow(x)),labels=nrow(x):1)


jpeg('analysis/AUC_over_windows_pvalues.jpeg')

x=pval.wind.matrix
levelplot(t(x[nrow(x):1, ]),col.regions=colorRampPalette(c('white', 'orange'))(length(quantile(x, probs = seq(0, 1, 0.01)) ) - 1L),xlab='WindowSize', ylab='-ve Log10 Pvalue',main='AUC  for 8008 positive vs. negative set over \ndifferent parameter values',
      breaks = quantile(x, probs = seq(0, 1, 0.01)))
          
dev.off()
          
dev.off()
     

jpeg('analysis/AUC_over_windows_rpgccutoffs.jpeg')

x=rpgc.wind.matrix
rownames(x)<-NULL
levelplot(aspect=3,t(x[nrow(x):1, ]),col.regions=colorRampPalette(c('white', 'orange'))(length(quantile(x, probs = seq(0, 1, 0.01)) ) - 1L),xlab='WindowSize', ylab='-ve Log10 Pvalue',main='AUC  for 8008 positive vs. negative set over \ndifferent parameter values',
          breaks = quantile(x, probs = seq(0, 1, 0.01)))
dev.off()

dev.off()





# recall at given precision --------------------------------------------

rec.cutoffs<-sapply(as.character(window.sizes),function(w){
  cagecountmats<-cagecountmatlist[[w]]
  cagecountmats.r<-cagecountmatlist.r[[w]]
  
  pvalrange<-1/(10^(1:pval.step.num))
  #Now we can go through each line and call the windows in a binary way
  #first go through our cutoffs and give each crm a linescore
  cutoff.linenum<-sapply(simplify=F,as.character(pvalrange),function(pval){#for each pval cutoff
    #derive the cutoffs for each library
    cutoffs<-sapply(lambdas,function(lam){ qpois(1-as.numeric(pval),lambda=lam*as.numeric(w)) })
    cutoffs<-cutoffs[colnames(cagecountmats[[1]])]
    #now go through our normalized libraries 
    lnmats<-sapply(simplify=F,cagecountmats,function(reg){
      tmp<-apply(reg,1,function(r){#for each locations
        sum(r>cutoffs)#count the lines with tagnum over our cutoff
      })
    })  
  })
  sapply(names(cutoff.linenum),function(cutoff){
        scores<-stack(list(pos=cutoff.linenum[[cutoff]][[1]],neg=cutoff.linenum[[cutoff]][[2]]))
        pred<-prediction(scores$values,scores$ind)
        perf<-performance(pred,'prec','rec')
        if(!any(perf@y.values[[1]]>prec.cutoff)%in%T){return(NA)}
        max(perf@x.values[[1]][perf@y.values[[1]]>prec.cutoff],na.rm=T)
  })
    
})



jpeg('analysis/recall_at_90pc_prec.jpeg')

x=rec.cutoffs
rownames(x)<-20:1
levelplot(t(x[nrow(x):1, ]),aspect=3,main='recall at 90pc precision for pvalue cutoff and windowsize combinations',xlab='windowsize',,ylab='Pval cutoff -ve Log10')
          
          
          
          dev.off()