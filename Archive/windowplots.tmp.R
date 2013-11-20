setwd('~/Harnett/TSS_CAGE_myfolder/')
rm(list=ls())

source('src/tss_cage_functions.R')
library(ggplot2)
#Load ROCR library for roc plots
library(ROCR,lib.loc='~/Harnett/R')
#define the number of rpgc and pvalue thresholds we'll try
rpgc.step.num<-100
pval.step.num<-20

#load the window table as 'z'
load('windowtable.robject.R')
load(file.accession.df)#load our data on the cage lines
load(file.cage.tag.rles)#load the actual cage data
cg<-cage.tag.rles

#load the regions we've defined in the makeregions script
regs=list(
  pos=import(con='analysis/positive_8008.bed',asRangedData=F,seqinfo=si),
  neg=import(con='analysis/negative_8008.bed',asRangedData=F,seqinfo=si),
  int=import(con='analysis/random.intergenic.bed',asRangedData=F,seqinfo=si),
  cad3pos=import(con='analysis/cad3.pos.bed',asRangedData=F,seqinfo=si),
  cad3neg=import(con='analysis/cad3.neg.bed',asRangedData=F,seqinfo=si)
)

#backup z for when we subset
z.bak<-z

#create matrices for recording parameter thresholds
rpgc.wind.matrix<-matrix(NA,nrow=rpgc.step.num,ncol=length(unique(z.bak$windowsize)))
pval.wind.matrix<-matrix(NA,nrow=pval.step.num,ncol=length(unique(z.bak$windowsize)))

#window sizes as columns
colnames(rpgc.wind.matrix)<-as.character( unique(z.bak$windowsize))
colnames(pval.wind.matrix)<-as.character( unique(z.bak$windowsize))

#duplicate these matrices for accuracies instead of aucs
rpgc.wind.matrix.acc<-rpgc.wind.matrix
pval.wind.matrix.acc<-pval.wind.matrix

#lists for when we record out matrix for results
cagecountmatlist<-cagecountmatlist.r<-list()




stopifnot(is.data.frame(z))
stopifnot(is.numeric(unique(z.bak$windowsize)))
stopifnot((unique(z.bak$acc))%in%names(cage.tag.rles))


# data in matrix format ---------------------------------------------------


for (w in unique(z.bak$windowsize)){
  
  cagecountmatlist[[w]]<-sapply(unique(z.bak$region),function(reg){
    sapply(unique(z.bak$acc),function(acc){
      z.bak[ z.bak$windowsize==w & z.bak$region ,]
    })
  })
  
  cagecountmatlist.r[[w]]<-sapply(unique(z.bak$region),function(reg){
    sapply(unique(z.bak$acc),function(acc){
      z.bak[ z.bak$windowsize==w & z.bak$region ,]$val.libnorm
    })
  })
  
}



message('Done with data matrices')
# first kind of fdr - total reads -----------------------------------------
for (w in unique(z.bak$windowsize)){
  
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
  scores<-stack(list(pos=sums[[1]],neg=sums[[2]][ regs[[2]]$score==0 ]))
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


for (w in unique(z.bak$windowsize)){
  cagecountmats<-cagecountmatlist[[w]]
  cagecountmats.r<-cagecountmatlist.r[[w]]
  
# Second kind of fdr - number of lines above fixed cutoff -----------------

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
  #or a roc
  
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




# #third kind of fdr - number of lines above fixed pvalue cutoff ----------


#first calculate the 

inter=import('analysis/intergenic.inactive.regions.bed')
inter=coverage(inter)>0
#first we'll define lambda for each line
#we can do this easily enough by just defining a 
lambdas<-sapply(names(cg),function(acc){
  zeros<-sum(as.numeric(sum(cg[[acc]][[1]][inter]==0))+  as.numeric(sum(cg[[acc]][[1]][inter]==0)))
  ones<-sum(as.numeric(sum(cg[[acc]][[1]][inter]==1))+  as.numeric(sum(cg[[acc]][[2]][inter]==1)))
  ones/(zeros+ones)#The poisson's density will be mostly concentrated on zero and 1
})

v=accession.df$library.size
names(v)<-accession.df$accession
qplot(v[names(lambdas)],lambdas,xlab='Library Size',main='Calculated Poisson Rate Vs. Library Size')
v=accession.df$intergenic.lib.size
names(v)<-accession.df$accession
qplot(v[names(lambdas)],lambdas,xlab='Intergenic Library Size',main='Calculated Poisson Rate Vs. Intergenic Library Size')


for (w in unique(z.bak$windowsize)){
  cagecountmats<-cagecountmatlist[[w]]
  cagecountmats.r<-cagecountmatlist.r[[w]]
    
  pvalrange<-1/(10^(1:pval.step.num))
  
  #Now we can go through each line and call the windows in a binary way
  #first go through our cutoffs and give each crm a linescore
  cutoff.linenum<-sapply(simplify=F,as.character(pvalrange),function(pval){#for each pval cutoff
    #derive the cutoffs for each library
    cutoffs<-sapply(lambdas,function(lam){ qpois(1-as.numeric(pval),lambda=lam*w) })
    cutoffs<-cutoffs[colnames(cagecountmats[[1]])]
    #now go through our normalized libraries 
    lnmats<-sapply(simplify=F,cagecountmats,function(reg){
      tmp<-apply(reg,1,function(r){#for each locations
        sum(r>cutoffs)#count the lines with tagnum over our cutoff
      })
    })  
  })

  ######Or using ROCR
  aucs<-sapply(names(cutoff.linenum),function(cutoff){
    scores<-stack(list(pos=cutoff.linenum[[cutoff]][[1]],neg=cutoff.linenum[[cutoff]][[2]]))
    pred<-prediction(scores$values,scores$ind)
    performance(pred,'auc')@y.values[[1]]
  })
  accs<-sapply(names(cutoff.linenum),function(cutoff){
    scores<-stack(list(pos=cutoff.linenum[[cutoff]][[1]],neg=cutoff.linenum[[cutoff]][[2]]))
    pred<-prediction(scores$values,scores$ind)
    performance(pred,'acc')@y.values[[1]]
  })
  
  auc.df<-as.data.frame(aucs)
  bestcut<-rownames(auc.df)[which.max(auc.df$aucs)]
  bestcut
  
  
  lowcuts<-names(cutoff.linenum)
  rownames(pval.wind.matrix)<-lowcuts
  pval.wind.matrix[,as.character(w)]<-aucs
  #pval.wind.matrix.acc[,as.character(w)]<-accs
  
  #logical vector describing if the cutoff gives at least a nonzero for each region.
  cutlog<-sapply(lowcuts,function(bestcut){ all( sapply(cutoff.linenum[[bestcut]][c(1,2,4,5)],function(reg){any(reg>0)}))})
  #only show cutoffs with some non zeros   
  cuts<-lowcuts[cutlog]
  #we need to exlclude 
  cuts<-cuts[seq(length(cuts)/3,length(cuts),length.out=3)]
  best<-lowcuts[which.max(auc.df$aucs)]
  
  #let's do four different cutoffs
  pdf(paste0('analysis/pois_pval.roc.',w,'.pdf'))
  
  
  #let's do four different cutoffs
  for (bestcut in c(best,cuts)){
      
    #  par(mfrow=c(2,2))
      
      
      
    #crms
    scores<-stack(list(pos=cutoff.linenum[[bestcut]][[1]],neg=cutoff.linenum[[bestcut]][[2]]))
    lcpred<-prediction(scores$values,scores$ind)
    perf<-performance(lcpred,'tpr','fpr')
    plot(perf,colorize=T)
    title(paste0('ROC curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),4),'  Tags,crm set'),sub=paste0('windowsize ',w))
    aucstring=paste0('AUC = ',round(performance(lcpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    abline(coef=c(0,1),lty=3)
    
    plot(performance(lcpred,'prec','rec'),colorize=T)
    title(paste0('Precision/Recall curve for Number of lines above Pval cutoff\n',round(as.numeric(bestcut),4),'  Tags,crm set'),sub=paste0('windowsize ',w))
    
    
    #crms
    scores<-stack(list(pos=cutoff.linenum[[bestcut]][[1]],neg=cutoff.linenum[[bestcut]][[2]][regs[[2]]$score==1]))
    lcpred<-prediction(scores$values,scores$ind)
    perf<-performance(lcpred,'tpr','fpr')
    plot(perf,colorize=T)
    title('ROC curve for Number of lines above Pval cutoff,8008 low k27ac negative set',sub=paste0('windowsize ',w))
    aucstring=paste0('AUC = ',round(performance(lcpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    abline(coef=c(0,1),lty=3)
    
    plot(performance(lcpred,'prec','rec'),colorize=T)
    title(paste0('Precision/Recall curve for Number of lines above Pval cutoff\n',round(as.numeric(bestcut),4),'  Tags,crm set'),sub=paste0('windowsize ',w))
    
    
    #crms vs random
    scores<-stack(list(pos=cutoff.linenum[[bestcut]][[1]],neg=cutoff.linenum[[bestcut]][[3]]))
    lcpred<-prediction(scores$values,scores$ind)
    perf<-performance(lcpred,'tpr','fpr')
    plot(perf,colorize=T)
    title(paste0('ROC curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),4),'  Tags,crm vs random'),sub=paste0('windowsize ',w))
    aucstring=paste0('AUC = ',round(performance(lcpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    abline(coef=c(0,1),lty=3)
    
    plot(performance(lcpred,'prec','rec'),colorize=T)
    title(paste0('Precision/Recall curve for Number of lines above Pval cutoff\n',round(as.numeric(bestcut),4),'  Tags,crm vs random'),sub=paste0('windowsize ',w))
    
    
    #cad
    scores<-stack(list(pos=cutoff.linenum[[bestcut]][[4]],neg=cutoff.linenum[[bestcut]][[5]]))
    lcpred<-prediction(scores$values,scores$ind)
    perf<-performance(lcpred,'tpr','fpr')
    plot(perf,colorize=T)
    title(paste0('ROC curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),4),'  Tags,CAD set'))
    aucstring=paste0('AUC = ',round(performance(lcpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    abline(coef=c(0,1),lty=3)
    
    plot(performance(lcpred,'prec','rec'),colorize=T)
    title(paste0('Precision/Recall curve for Number of lines above Pval cutoff\n',bestcut,'  Tags,CAD set'),sub=paste0('windowsize ',w))
  }
  dev.off()
  
}

write.table(rpgc.wind.matrix,'analysis/rpgc.window.mat.txt')
write.table(pval.wind.matrix,'analysis/pval.window.mat.txt')

dev.off()
save.image('tmp.image.R')

# experimenting with negative binomail distributions ----------------------




