rm(list=ls())
setwd('~/Harnett/TSS_CAGE_myfolder/')
# Another windows script. Now just do all 8008 so we can label lat --------
source('src/tss_cage_functions.R')
source('src/generate_Rle.R')
library(reshape2)
library(ggplot2)
library(ROCR,lib.loc='~/Harnett/R')


load(file.crmgrs)
load(file.cage.tag.rles)
load(file.accession.df)



cagemat<-get.best.window.mat(crmgrs,100,cage=cg)







# Get the score in the best window for each region ------------------------
sapply(simplify=F,regs,function(reg){ 
    



# create rpgc data and save --------------------------------------------------------

names(cagecountmatlist)<-window.sizes
#add a column of normalized values

####create the matrix of normalized counts.

v<-accession.df$genome.coverage
names(v)<-accession.df$accession

cagecountmatlist.r<-sapply(simplify=F,names(cagecountmatlist),function(w){
  sapply(simplify=F,names(cagecountmatlist[[w]]),function(reg){
    sapply(1:ncol(cagecountmatlist[[w]][[reg]]),function(n){
      cagecountmatlist[[w]][[reg]][,n]<-cagecountmatlist[[w]][[reg]][,n]/v[names(cg)[n]]
    })
  })
})

setlist=list(
  tfset=list(crmgrs$in.tf.pos,crmgrs$in.tf.neg),
  tf.27ac.set=list(crmgrs$in.tf.pos & crmgrs$H3K27ac,crmgrs$in.tf.neg & crmgrs$low.dnase & !crmgrs$H3K27ac),
  tf.pol.set=list(crmgrs$in.tf.pos & crmgrs$,crmgrs$in.tf.neg & crmgrs$low.dnase & !crmgrs$H3K27ac),
  tf.K79.set=list(crmgrs$in.tf.pos & crmgrs$high.dnase,crmgrs$in.tf.neg & crmgrs$low.dnase & !crmgrs$H3K27ac),
  tf.K36.set=list(crmgrs$in.tf.pos & crmgrs$high.dnase,crmgrs$in.tf.neg & crmgrs$low.dnase & !crmgrs$H3K27ac),
)


ROCfunc <- function (pos,neg,main,sub) {
  scores<-stack(list(pos=pos,neg=neg))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T,lwd=6)
  length(posset)
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

ROCfunc(pos=rnorm(100,2),
        neg=rnorm(100,1),
        main=paste0('ROC curve Testing'),
        sub=paste0('windowsize '))


PreRecfunc <- function (pos,neg,main,sub) {
  scores<-stack(list(pos=pos,neg=neg))
  sumpred<-prediction(scores$values,scores$ind)
  plot(performance(sumpred,'prec','rec'),colorize=T)
    
  posstring=paste0('Positive Set - ', length(pos))
  negstring=paste0('Negative Set - ', length(neg))
  text(x=0.85,y=0.2,labels=posstring)
  text(x=0.85,y=0.1,labels=negstring)
  
  title(main,sub)
  grid(lwd=2)
}


PreRecfunc(pos=rnorm(100,2),
        neg=rnorm(100,1),
        main=paste0('ROC curve Testing'),
        sub=paste0('windowsize '))





sum.cutoff<-list()

# first kind of fdr - total reads -----------------------------------------
for (w in as.character(window.sizes)){
  sum.cutoff[[w]]=list()
               
  for(setn in names(setlist)){
    
    posset=setlist[[setn]][[1]]
    negset=setlist[[setn]][[2]]
    #select relevant matrix
    cagecountmats<-cagecountmatlist[[w]]
    cagecountmats.r<-cagecountmatlist.r[[w]]
    
    #we need the annotation from the 'makeregions' script which should be in crmgrs
    if( is.null(crmgrs$high.dnase) | is.null(crmgrs$in.tf.pos)){stop('Missing annotation of crm object - run makeregions.')}
    stopifnot(!NA%in%cagecountmats.r[[1]][,1] )#check cagecountmats.r looks okay
    
    #calculate the summed 
    sums<-sapply(simplify=F,cagecountmats,function(reg){rowSums(reg)})#total in all lines
    sums.r<-sapply(simplify=F,cagecountmats.r,function(reg){rowSums(reg)})#normalized
    
    pdf(paste0('analysis/rocplots_pos_neg/summ.roc.',setn,'.',w,'.pdf'))
    
#     
#     ROCfunc(pos=crmgrs$allsum.rpgc[ crmgrs$in.tf.pos],neg=crmgrs$allsum.rpgc[ crmgrs$in.tf.neg],
#             main='ROC curve for Summed Unnormalized Tags,Pos crms vs Neg Crms',
#             sub=paste0('windowsize ',w))
#       ROCfunc(pos=cagecountmatlist[[w]][[1]][posset],neg=cagecountmatlist[[w]][[1]][negset],
#                   main='ROC curve for Summed Unnormalized Tags,Pos crms vs Neg Crms',
#                   sub=paste0('windowsize ',w))
#     
#     
    
    ROCfunc(pos=sums[[1]][posset],neg=sums[[1]][negset],
            main='ROC curve for Summed Unnormalized Tags,Pos crms vs Neg Crms',
            sub=paste0('windowsize ',w))
    
    ROCfunc(pos=sums.r[[1]][posset],neg=sums.r[[1]][negset],
            main='ROC curve for Summed Normalized Tags,Pos crms vs Neg Crms',
            sub=paste0('windowsize ',w))
    
    ROCfunc(pos=sums[[1]][posset],neg=sums[[2]],
            main='ROC curve for Summed Unnormalized Tags,Pos crms vs random inactive',
            sub=paste0('windowsize ',w))
    
    ROCfunc(pos=sums.r[[1]][posset],neg=sums.r[[2]],
            main='ROC curve for Summed Normalized Tags,Pos crms vs random inactive',
            sub=paste0('windowsize ',w))
    
    #8008 negative vs. the intergenic set
    ROCfunc(pos=sums[[1]][negset],neg=sums[[2]],
    main='ROC curve for Summed Unnormalized Tags,Neg crms vs random inactive',
            sub=paste0('windowsize ',w))
    
    #8008 negative vs. the intergenic set
    ROCfunc(pos=sums.r[[1]][negset],neg=sums.r[[2]],
            main='ROC curve for Summed Normalized Tags,Neg crms vs random inactive',
            sub=paste0('windowsize ',w))
    
    PreRecfunc(pos=sums[[1]][posset],neg=sums[[1]][negset],
            main='Precision/Recall curve for Summed Unnormalized Tags,Pos crms vs Neg Crms',
            sub=paste0('windowsize ',w))
    
    PreRecfunc(pos=sums.r[[1]][posset],neg=sums.r[[1]][negset],
            main='Precision/Recall for Summed Normalized Tags,Pos crms vs Neg Crms',
            sub=paste0('windowsize ',w))
    
    PreRecfunc(pos=sums[[1]][posset],neg=sums[[2]],
            main='Precision/Recall for Summed Unnormalized Tags,Pos crms vs random inactive',
            sub=paste0('windowsize ',w))
    
    PreRecfunc(pos=sums.r[[1]][posset],neg=sums.r[[2]],
            main='Precision/Recall for Summed Normalized Tags,Pos crms vs random inactive',
            sub=paste0('windowsize ',w))
    
    #8008 negative vs. the intergenic set
    PreRecfunc(pos=sums[[1]][negset],neg=sums[[2]],
            main='Precision/Recall for Summed Unnormalized Tags,Neg crms vs random inactive',
            sub=paste0('windowsize ',w))
    
    #8008 negative vs. the intergenic set
    PreRecfunc(pos=sums.r[[1]][negset],neg=sums.r[[2]],
            main='Precision/Recall for Summed Normalized Tags,Neg crms vs random inactive',
            sub=paste0('windowsize ',w))
    
    dev.off()
    
    
    
  } 
  
  scores<-stack(list(pos=sums.r[[1]][posset],neg=sums.r[[1]][negset]))
  pred<-prediction(scores$values,scores$ind)
  perf.pred<-performance(pred,"prec") #find list of accuracies
  sum.cutoff[[w]][[setn]]=   min(perf.pred@x.values[[1]][ perf.pred@y.values[[1]]>prec.cutoff ],na.rm=T)    
  
}






# Cutoffs using poisson distribution --------------------------------------



# lambdas for third kind of fdr -------------------------------------------

#get the intergenic regions
inter=import('analysis/intergenic.inactive.regions.bed')
inter=coverage(inter)>0
#first we'll define lambda for each line
if(file.exists(file.lambdas)){load(file.lambdas)
}else{
  lambdas<-sapply(names(cg),function(acc){
    zeros<-sum(as.numeric(sum(cg[[acc]][[1]][inter]==0))+  as.numeric(sum(cg[[acc]][[1]][inter]==0)))
    ones<-sum(as.numeric(sum(cg[[acc]][[1]][inter]==1))+  as.numeric(sum(cg[[acc]][[2]][inter]==1)))
    ones/(zeros+ones)#The poisson's density will be mostly concentrated on zero and 1
  })
  save(lambdas,file=file.lambdas)
}





# #third kind of fdr - number of lines above fixed pvalue cutoff ----------

pval.wind.matrix<-list()
recall=list()
cutoffs=list()

for(setn in names(setlist)){
  pval.wind.matrix[[setn]]<-list()
  recall[[setn]]<-list()
  cutoffs[[setn]]<-list()
  
  for (w in as.character(window.sizes)){
    pval.wind.matrix[[setn]][[w]]<-list()
    recall[[setn]][[w]]<-list()
    cutoffs[[setn]][[w]]<-list()
    
    posset=setlist[[setn]][[1]]
    negset=setlist[[setn]][[2]]
      
    #get the relevant data
    cagecountmats<-cagecountmatlist[[w]]
    cagecountmats.r<-cagecountmatlist.r[[w]]
    #define our range of pvalues
    pvalrange<-1/(10^(1:pval.step.num))
    
    #Now we can go through each line and call the windows in a binary way
    #first go through our cutoffs and give each crm a linescore
    cutoff.linenum<-sapply(simplify=F,as.character(pvalrange),function(pval){#for each pval cutoff
    #derive the cutoffs for each library from our lambdas
      cutoffs<-sapply(lambdas,function(lam){ qpois(1-as.numeric(pval),lambda=lam*as.numeric(w)) })
      #in the right order, using only the libraries in our matrix
      cutoffs<-cutoffs[colnames(cagecountmats[[1]])]
      #now go through our NONnormalized libraries 
      lnmats<-sapply(simplify=F,cagecountmats,function(reg){
        tmp<-apply(reg,1,function(r){#for each locations
          sum(r>cutoffs)#count the lines with tagnum over our cutoff
        })
      })  
   })
    #we now have an object which, for each pval cutoff, gives us the number of libraries above this threshold 
    
    
    #Now we'll go through each cutoff and get the auc for it
    aucs<-sapply(names(cutoff.linenum),function(cutoff){
      scores<-stack(list(pos=cutoff.linenum[[cutoff]][[1]][posset],neg=cutoff.linenum[[cutoff]][[1]][negset]))
      pred<-prediction(scores$values,scores$ind)
      performance(pred,'auc')@y.values[[1]]
    })
    #And the accuracy at our chosen threshold
    accs<-sapply(names(cutoff.linenum),function(cutoff){
      scores<-stack(list(pos=cutoff.linenum[[cutoff]][[1]][posset],neg=cutoff.linenum[[cutoff]][[1]][negset]))
      pred<-prediction(scores$values,scores$ind)
      performance(pred,'acc')@y.values[[1]]
    })
    
    #and the recall at 95% precision
    recall[[setn]][[w]]<-
      sapply(names(cutoff.linenum),function(cutoff){
      scores<-stack(list(pos=cutoff.linenum[[cutoff]][[1]][posset],neg=cutoff.linenum[[cutoff]][[1]][negset]))
      pred<-prediction(scores$values,scores$ind)
      perf<-performance(pred,'prec','rec')
      plot(perf)
      if(!any(perf@y.values[[1]]>prec.cutoff)%in%T){return(NA)}
      max(perf@x.values[[1]][perf@y.values[[1]]>prec.cutoff],na.rm=T)
    })
    
      
    
    #and the 
    cutoffs[[setn]][[w]]<-
      sapply(names(cutoff.linenum),function(cutoff){
        scores<-stack(list(pos=cutoff.linenum[[cutoff]][[1]][posset],neg=cutoff.linenum[[cutoff]][[1]][negset]))
        pred<-prediction(scores$values,scores$ind)
        perf.pred<-performance(pred,"prec") #find list of accuracies
        min(perf.pred@x.values[[1]][ perf.pred@y.values[[1]]>prec.cutoff ],na.rm=T)    
    })
    
    #get the cutoff with the best auc
    pval.wind.matrix[[setn]][[w]] <-as.data.frame(aucs)
    rownames( pval.wind.matrix[[setn]][[w]] )<-names(cutoff.linenum)
    bestcut<-rownames( pval.wind.matrix[[setn]][[w]])[which.max( pval.wind.matrix[[setn]][[w]]$aucs)]
      
    lowcuts<-names(cutoff.linenum)
    #logical vector describing if the cutoff gives at least a nonzero for each region.
    cutlog<-sapply(lowcuts,function(bestcut){ all( sapply(cutoff.linenum[[bestcut]][1],function(reg){any(reg>0)}))})
    #only show cutoffs with some non zeros   
    cuts<-lowcuts[cutlog]
    #we need to exlclude 
    cuts<-cuts[seq(length(cuts)/3,length(cuts),length.out=3)]
    best<-lowcuts[which.max( pval.wind.matrix[[setn]][[w]]$aucs)]
    
    
    #let's do four different cutoffs
    pdf(paste0('analysis/rocplots_pos_neg/pois_pval.roc.',setn,'.',w,'.pdf'))
    
    #let's do four different cutoffs
    for (bestcut in c(best,cuts)){
      
      #  par(mfrow=c(2,2))
            
      ROCfunc(pos=cutoff.linenum[[bestcut]][[1]][posset],
              neg=cutoff.linenum[[bestcut]][[1]][negset],
              main=paste0('ROC curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),4),'  Tags,crm set'),
              sub=paste0('windowsize ',w))
      
    
      ROCfunc(pos=cutoff.linenum[[bestcut]][[1]][posset],
              neg=cutoff.linenum[[bestcut]][[2]],
              main=paste0('ROC curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),4),'  Tags, Pos crm vs. random'),
              sub=paste0('windowsize ',w))
      
      ROCfunc(pos=cutoff.linenum[[bestcut]][[1]][negset],
              neg=cutoff.linenum[[bestcut]][[2]],
              main=paste0('ROC curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),6),'  Tags, Neg crm vs. random'),
              sub=paste0('windowsize ',w))
      
      PreRecfunc(pos=cutoff.linenum[[bestcut]][[1]][posset],
              neg=cutoff.linenum[[bestcut]][[1]][negset],
              main=paste0('Prec/Recall curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),4),'  Tags,crm set'),
              sub=paste0('windowsize ',w))
      
      
      PreRecfunc(pos=cutoff.linenum[[bestcut]][[1]][posset],
              neg=cutoff.linenum[[bestcut]][[2]],
              main=paste0('Prec/Recall curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),4),'  Tags, Pos crm vs. random'),
              sub=paste0('windowsize ',w))
      
      PreRecfunc(pos=cutoff.linenum[[bestcut]][[1]][negset],
              neg=cutoff.linenum[[bestcut]][[2]],
              main=paste0('Prec/Recall curve for Number of lines above Pval cutoff ',round(as.numeric(bestcut),4),'  Tags, Neg crm vs. random'),
              sub=paste0('windowsize ',w))
      
      
    
    }
  dev.off()
  
  }

}





# Plots of parameter scanning ---------------------------------------------







setn=names(setlist)[1]

for(setn in names(setlist)){
  
  #now create our heatmap for different pvalues
  pdf(paste0('analysis/rocplots_pos_neg/AUC_over_windows_pvalues.',setn,'.pdf'))
  
    x<-as.matrix(do.call(cbind,pval.wind.matrix[[setn]]))
    colnames(x)=names(pval.wind.matrix[[setn]])
    levelplot(t(x[nrow(x):1, ]),col.regions=colorRampPalette(c('white', 'orange'))(length(quantile(x, probs = seq(0, 1, 0.01)) ) - 1L),xlab='WindowSize', ylab='-ve Log10 Pvalue',
              main=paste0('AUC  for 8008 positive vs. negative set over \ndifferent parameter values - ',setn),
            breaks = quantile(x, probs = seq(0, 1, 0.01)))
  
  dev.off()
  
}

for(setn in names(setlist)){
  
  #now create our heatmap for different pvalues
  pdf(paste0('analysis/rocplots_pos_neg/Recall_at_95prec',setn,'.pdf'))
  
  x<-as.matrix(do.call(cbind,recall[[setn]]))
  colnames(x)=names(recall[[setn]])
  levelplot(t(x[nrow(x):1, ]),col.regions=colorRampPalette(c('white', 'lightgreen'))(length(quantile(x, probs = seq(0, 1, 0.01)) ) - 1L),xlab='WindowSize', ylab='-ve Log10 Pvalue',
            main=paste0('Recall at 95pc Prec  for 8008 positive vs. negative set over \ndifferent parameter values - ',setn),
            breaks = quantile(x, probs = seq(0, 1, 0.01)))
  
  dev.off()
  
}


