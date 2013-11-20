setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
library(reshape2)
library(ggplot2)

load(file.cage.tag.rles)
load(file.accession.df)
#eliminate reseqs since they're already added to their originals
regs=list(
  pos=import(con='analysis/positive_8008.bed',asRangedData=F,seqinfo=si),
  neg=import(con='analysis/negative_8008.bed',asRangedData=F,seqinfo=si),
  int=import(con='analysis/random.intergenic.bed',asRangedData=F,seqinfo=si),
  cad3pos=import(con='analysis/cad3.pos.bed',asRangedData=F,seqinfo=si),
  cad3neg=import(con='analysis/cad3.neg.bed',asRangedData=F,seqinfo=si)
)
# 
# #for speed on rstudio
# regs[[3]]<-regs[[3]][1:500]
# cg<-cage.tag.rles[!grepl('reseq',names(cage.tag.rles))]
# cg<-cg[1:5]
# accession.df<-accession.df[1:5,]



#get matrix with windows etc. 
cagecountmats<-sapply(regs,function(reg){
  reg.r<-as(reg,'RangesList')
  sapply(names(cg),function(acc){
    unlist(viewSums(Views(cg[[acc]]$pos[chrs.keep],reg.r[chrs.keep])))
  })+
    sapply(names(cg),function(acc){
      unlist(viewSums(Views(cg[[acc]]$pos[chrs.keep],reg.r[chrs.keep])))
    })   
})

  ####__####stopped here
  
# density plot to illustrate cutoff ---------------------------------------
# The more complex method - we try a cutoff, and a nmber of lines ---------
#we go through a selection of RPGC cutoffs

highnum<-1000#arbitrarily set a hight limit on our rpgc cutoff
all<-do.call(rbind,cagecountmats.r)
lowvals<-sort((unique(as.vector(all))))
lowvals<-lowvals[lowvals <  highnum/rev(accession.df$genome.coverage)[1]]
lowvals<-lowvals[seq(from=1,to=length(lowvals),length.out=100)]

#now we can cycle through each cutoff, and for each of our sets, determine
#if it's above the criteria
#for(lfrac in seq(0,1,by=0.1)){
#  lnum=floor(length(cg)*lfrac)#get actual number of libraries we need


cutoff.linenum<-
sapply(simplify=F,as.character(unique(lowvals)),function(val){#now pick an RPGC value
  #now go through our normalized libraries 
  sapply(simplify=F,cagecountmats.r,function(reg){
    tmp<-apply(reg,1,function(r){#for each locations
      sum(r>as.numeric(val))#count the lines with tagnum over our cutoff
    })
  }) 
})

c.df<-melt(cutoff.linenum)
colnames(c.df)<-c('library.count','region','cutoff')
#tail(c.df)
c.df$line.fraction<-c.df$library.count/length(cg)

#now, we can plot the distribution of reproducibility score above certain cutoffs
#first find the combination which gives use the most positives at 95% spec
bestpos<-0
bestlf<-NULL
bestcut<-NULL
bestpos<-NULL
#we don't want to calculate for every possible fraction line number cutoff
#just a selection of them
fract.cutoffs<-unique(c.df$line.fraction)
fract.cutoffs<-fract.cutoffs[seq(1,length(fract.cutoffs),length.out=5)]
rpgccutoffs<-unique(c.df$cutoff)[seq(1,length(unique(c.df$cutoff)),length.out=10)]#actually just integers for now

senslist<-
sapply(simplify=F,rpgccutoffs,function(lowval){
  #sapply(fract.cutoffs,function(linefraction){
    cat('.')
    #we want to find the 95% quantile of the linecount for the negatives
    p=c.df[c.df$cutoff ==lowval & c.df$region == 'pos',]
    n=c.df[c.df$cutoff ==lowval & c.df$region == 'neg',]
    cut=quantile(n$library.count,0.95)
    posnum=sum(p$library.count>cut)
    posnum
 # })
})

maxlibcount<-max(c.df$library.count)
sapply(simplify=F,rpgccutoffs,function(lowval){
  #sapply(fract.cutoffs,function(linefraction){
  cat('.')
  #we want to find the 95% quantile of the linecount for the negatives
  p=c.df[c.df$cutoff ==lowval & c.df$region == 'pos',]
  n=c.df[c.df$cutoff ==lowval & c.df$region == 'neg',]
  cut=quantile(n$library.count,0.95)
  posnum=sum(p$library.count>cut)
  posnum
  # })
})

#list of sensitivies at different cutoffs for rpgc and then linenumber
senslist


lowval=lowvals[99]
n=c.df[c.df$cutoff ==lowval & c.df$region == 'neg',]
n$library.count



#plots of the reproducibility distribution
lowval=rpgccutoffs[4]
jpeg('tmp.jpeg',width=1000,height=1000)
libcountdata<-unlist(list(pos=c.df[c.df$cutoff ==lowval & c.df$region == 'pos',]$library.count,
            neg=c.df[c.df$cutoff ==lowval & c.df$region == 'neg',]$library.count,
            int=c.df[c.df$cutoff ==lowval & c.df$region == 'int',]$library.count
          #  cad3pos=c.df[c.df$cutoff ==lowval & c.df$region == 'cad3pos',]$library.count,
           # cad3neg=c.df[c.df$cutoff ==lowval & c.df$region == 'cad3neg',]$library.count)
))
libcountdata<-data.frame(lib.count=libcountdata,region= gsub('(.*?)\\d+','\\1',names(libcountdata)))
ggplot(libcountdata,aes(x=lib.count,color=factor(region)))+
  geom_density()+
  scale_x_discrete(name='proportion of libraries with tag number above cutff')+
  ggtitle(paste0('Density Plot of "Line Number" using RPGC cutoff of ',lowval))  
dev.off()

#now produce some density plots showing the distribution of lscore for pos neg and int, when
#considering various cutoffs.

# plots -------------------------------------------------------------------


plot(density(log10(sums$pos)),ylim=c(0,0.60),col='red',xlab='log 10 tag count',
     main='log10 density of summed cage librarys')
points(density(log10(sums$neg)),col='blue',type='l')
points(density(log10(sums$int)),col='green',type='l')
legend(x='topright',fill=c('red','blue','green'),legend=c('positive','negative','intergenic'))


fill_density_plot(density(log10(sums$pos)),cutoff.all.neg,'pink')


for(regn in seq_along(regs)){
  regs[[regn]]$called<-regs[[regn]]$cutoff.all.neg  
}

# sensitivity and specificity ---------------------------------------------
#print specificity according to the 8008
#prop of negatives correctly classified
spec.neg.8008<-mean(!regs[['neg']]$called.n.8)
sens.neg.8008<-mean(!regs[['pos']]$called.n.8)

spec.neg.8008<-mean(!regs[['neg']]$called.n.8)
sens.neg.8008<-mean(!regs[['pos']]$called.n.8)

sens.neg.8008<-mean(!regs[['neg']]$called)
sens.int.8008<-mean(!regs[['int']]$called)
sens.neg.cad3<-mean(!cad3.neg[['neg']]$called)



cat('Specificity:')


#and to CAD



#and now the rpgc
#normalized plots now
sums<-sapply(simplify=F,allmats.r,function(reg){
  rowSums(reg)
})
plot(density(log10(sums$pos)),col='red',ylim=c(0,0.8),xlab='log10 normalized tag count',
     main='log10 density of summed cage librarys')
points(density(log10(sums$neg)),col='blue',type='l')
points(density(log10(sums$int)),col='green',type='l')
legend(x='topright',fill=c('red','blue','green'),legend=c('positive','negative','intergenic'))








# density plots for individual lines --------------------------------------

m<-(melt(list(pos=allmats$pos,neg=allmats$neg,intergenic=allmats$int)))
ggplot(m,aes(x=value,color=L1))+geom_density()+scale_x_log10()
ggplot(m[m$L1=='pos',],aes(x=value,color=Var2))+geom_density()+scale_x_log10()



m<-(melt(list(pos=allmats.r$pos,neg=allmats.r$neg,intergenic=allmats.r$int)))
ggplot(m,aes(x=value,color=L1))+geom_density()+scale_x_log10()
ggplot(m[m$L1=='pos',],aes(x=value,color=Var2))+geom_density()+scale_x_log10()




#do individual plots for the lines
ggplot(m[m$L1=='pos'&m$Var2%in%accession.df$accession[c(1:3,70:73)],],aes(x=value,color=Var2))+geom_line(stat='density')+scale_x_log10()+
  ggtitle('Reads in Positive CRM set')

ggplot(m[m$L1=='neg'&m$Var2%in%accession.df$accession[c(1:3,70:73)],],aes(x=value,color=Var2))+geom_density()+scale_x_log10()
ggplot(m[m$L1=='int'&m$Var2%in%accession.df$accession[c(1:3,70:73)],],aes(x=value,color=Var2))+geom_density()+scale_x_log10()
ggplot(m[m$Var2%in%accession.df$accession[c(1:3,70:73)],],aes(x=value,color=L1))+geom_density()+scale_x_log10()

              
#see if these normalize any better for rpgc versions
plot(density(rowSums(allmats.r$pos)),col='red')
points(density(rowSums(allmats.r$neg)),col='blue',type='l')
points(density(rowSums(allmats.r$int)),col='green',type='l')






#now melt into a dataframe
w.df<-melt(window.sum.distributions)
colnames(w.df)<-c('value','acc','windowsize','region')
h(w.df)

v=accession.df$genome.coverage
names(v)<-as.character(accession.df$accession)
w.df$val.libnorm<- w.df$value/ v[w.df$acc]
# 
# v=accession.df$intergenic.lib.size
# names(v)<-as.character(accession.df$accession)
# w.df$val.intlibnorm<- w.df$value/ v[w.df$acc]
# 
# v=accession.df$tss.lib.size
# names(v)<-as.character(accession.df$accession)
# w.df$val.tsslibnorm<- w.df$value/ v[w.df$acc]


w.df$nonzero<-as.numeric(w.df$value>0)
#accs<-names(cage.tag.rles)[c(1,5,10,15,20,25,30,35,40,45,50)]
w.df.bak<-w.df
write.table(quote=F,row.names=F,col.names=T,w.df,file='analysis/windowsize.crm.scores.table.txt')
#w.df<-read.table('analysis/windowsize.crm.scores.table.txt')
#w.df<-w.df[w.df$acc%in%accs,]



ggplot(w.df,aes(x=acc,fill=windowsize,y=nonzero,facet=region))+
  stat_summary(fun.y='mean',geom='bar',position='dodge')+
  facet_wrap(~region)+
  ggtitle('proportion of regions with a maximum score of more than zero, for\n intergenic, negative and positive regions')


ggplot(w.df[ w.df$windowsize=='w500', ],aes
    (y=value,x=acc,fill=acc,facet=region))+geom_boxplot()+
  facet_wrap(~region)+scale_y_log10()+ggtitle('')
  

ggplot(w.df[ w.df$windowsize=='w500' , ],aes(y=val.libnorm,x=acc,fill=acc,facet=region))+geom_boxplot()+
  facet_wrap(~region)+scale_y_log10()


ggplot(w.df[ w.df$windowsize=='w500', ],aes(y=val.libnorm,fill=region,x=acc))+geom_violin()+
  scale_y_log10()



regs$int$alltagsum<-sums$int
regs$pos$alltagsum<-sums$pos
regs$neg$alltagsum<-sums$neg

regs$neg[order(regs$neg$alltagsum,decreasing=T),]
regs$int[order(regs$int$alltagsum,decreasing=T),]$alltagsum[1:100]










#now with a cutoff and window size decided on, we can ask, 
#we go through each line
#95% cutoff for different sets

#sapply(1:20,function(cut){mean(w.df[  w.df$region=='int'   ,]$value>cut)})

tmp<-sapply(accession,function(acc){
  w.df<-w.df[ w.df$acc==acc &w.df$windowsize=='w500',]
  sapply(1:30,function(cut){sum(w.df[  w.df$region=='neg'   ,]$value >=cut) /( sum(w.df[  w.df$region=='neg'   ,]$value >=cut) +sum(w.df[  w.df$region=='pos'   ,]$value>=cut))})
  
  #sapply(1:30,function(cut){mean(w.df[  w.df$region=='neg'   ,]$value >=cut) / mean(w.df[  w.df$region=='pos'   ,]$value>=cut)})
})
matplot(tmp[,1:20],type='l')
#fdrs at various cutoffs
tmp<-sapply(accession,function(acc){
  w.df<-w.df[ w.df$acc==acc &w.df$windowsize=='w500',]
  sapply(1:30,function(cut){sum(w.df[  w.df$region=='neg'   ,]$val.libnorm >=cut) /( sum(w.df[  w.df$region=='neg'   ,]$value >=cut) +sum(w.df[  w.df$region=='pos'   ,]$value>=cut))})
  
  #sapply(1:30,function(cut){mean(w.df[  w.df$region=='neg'   ,]$value >=cut) / mean(w.df[  w.df$region=='pos'   ,]$value>=cut)})
})
w.df$
matplot(tmp[,1:20],type='l')

icut<-quantile(w.df[  w.df$region=='int'   ,]$value,0.95)
ncut<-quantile(w.df[ w.df$region=='neg'   ,]$value,0.95)




