setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/generate_Rle.R')
source('src/generate_chromdata_cage.R')
source('src/load_annotations.R')

cg<-cage.tag.rles
cgr<-cage.tag.rles.rpgc
# now create an object with all potential crms ----------------------------

#concatenate the crm8008.gr and cad3.gr objects
cad3.pos<-cad3.gr[ cad3.df$active68 %in% T & cad3.df$intergenic ]
cad3.neg<-cad3.gr[ cad3.df$inactive68 %in% T & cad3.df$intergenic ]
pos.crms<-keepSeqlevels(cad3.pos,chrs.keep)
neg.crms<-keepSeqlevels(cad3.neg,chrs.keep)



#resize them, for now, to 500bp
pos.crms<-resize(pos.crms,500,fix='center')
neg.crms<-resize(neg.crms,500,fix='center')
#possibly filter out the noninergenic ones again


#and our 'inactive' windows, which are the same, only without even K4me1
act<-list(
  crm8008.gr,
  cad3.gr,
  transcripts.gr,
  chrompeaks[['K4me3_6-8h']],
  chrompeaks[['K4me3_4-6h']],
  chrompeaks[['K4me1_4-6h']],
  chrompeaks[['K4me1_6-8h']],
  chrompeaks.modencode[['K4me3_4-8h']],
  chrompeaks.modencode[['K4me1_4-8h']],
  chrompeaks.modencode[['K27ac_4-8h']],
  chrompeaks[['PolII_6-8h']]
)

act<-sapply(act,  function(gr){ mcols(gr) <-NULL;gr})
act<-do.call('c',act)
act<-resize(act,width(act)+1000,fix='center')
act<-reduce(act)
act<-keepSeqlevels(act,chrs.keep)
intergenic.inactive.regions<-gaps(act)


#now we need to start picking random windows in our regions
#function to pick windows of size x in a set 
windowGRange<-function(gr,x){
  gr<-gr[width(gr)>=x]
  gr<-resize(gr,width(gr)-(x-1),fix='start')
  c<-coverage(gr)
  c<-which(c>0)
  c<-unlist(c)
  return(
    (GRanges(names(c),IRanges(c,width=x)))
  )
}


int<-windowGRange(intergenic.inactive.regions,500)
int<-sample(int,size=10000)
regs<-list(pos=pos.crms,neg=neg.crms,int=int)

# Basic state for max and min ---------------------------------------------
#make sizes the same, make intergenic windows
#measure 

#get matrix with windows etc. 
allmats<-sapply(regs,function(reg){
  reg.r<-as(reg,'RangesList')
  sapply(names(cg),function(acc){
    unlist(viewSums(Views(cg[[acc]]$pos,reg.r)[seqlevels(reg.r)]))
  })+
    sapply(names(cg),function(acc){
      unlist(viewSums(Views(cg[[acc]]$pos,reg.r)[seqlevels(reg.r)]))
    })   
})
allmats.r<-sapply(allmats,function(reg){
  t(apply(reg,1,function(x)x/accession.df$library.size))
})

sums<-sapply(allmats,function(reg){
  rowSums(reg)
})
plot(density(log10(sums$pos)),ylim=c(0,0.60),col='red',xlab='log 10 tag count',
     main='log10 density of summed cage librarys - cad3 set')
points(density(log10(sums$neg)),col='blue',type='l')
points(density(log10(sums$int)),col='green',type='l')
legend(x='topright',fill=c('red','blue','green'),legend=c('positive','negative','intergenic'))


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


m<-melt(sums)
h(m)
ggplot(m,aes(color=factor(L1),x=as.numeric(value)))+geom_density()+scale_y_log10()
# allmats.r<-sapply(regs,function(reg){
#   reg.r<-as(reg,'RangesList')
#   sapply(names(cg),function(acc){
#     unlist(viewSums(Views(cgr[[acc]]$pos,reg.r)[seqlevels(reg.r)]))
#   })+
#     sapply(names(cg),function(acc){
#       unlist(viewSums(Views(cgr[[acc]]$pos,reg.r)[seqlevels(reg.r)]))
#     })   
# })

#now plot distributions for them all

m<-(melt(list(pos=allmats$pos,neg=allmats$neg,intergenic=allmats$int)))
ggplot(m,aes(x=value,color=L1))+geom_density()+scale_x_log10()
ggplot(m[m$L1=='pos',],aes(x=value,color=Var2))+geom_density()+scale_x_log10()



m<-(melt(list(pos=allmats.r$pos,neg=allmats.r$neg,intergenic=allmats.r$int)))
ggplot(m,aes(x=value,color=L1))+geom_density()+scale_x_log10()
ggplot(m[m$L1=='pos',],aes(x=value,color=Var2))+geom_density()+scale_x_log10()



# 
# plot(density(rowSums(allmats$pos)),col='red',ylim=c(0,0.00003))
# points(density(rowSums(allmats$neg)),col='blue',type='l')
# points(density(rowSums(allmats$int)),col='green',type='l')


m<-melt(allmats)

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







# Function to pick the best window in a set of ranges ---------------------
window.sizes<-c('w5'=5,'w10'=10,'w25'=25,'w50'=50,'w100'=100,'w150'=150,'w250'=250,'w500'=500)
cg<-cg[1:10]
window.sizes=window.sizes[c(3,5,8)]

window.sum.distributions<-sapply(simplify=F,regs,function(reg){
  sapply(simplify=F,window.sizes,function(w){
    # sapply(simplify=F,c(10,100,500),function(w){
    cat('.')
    #make sure all our regions are at least the window size
    reg<-resize(reg,width=pmax(width(reg),w),fix='center')
    winds<-windowGRange(reg,w)
    #winds<-sample(winds,size=min(10000,length(winds)))#sample of our windows
    windregs=findOverlaps(winds,reg,select='first')
    f=findOverlaps(winds,reg)
    regnums<-1:length(reg)
    stopifnot(length(winds)>100)
    winds<-as(winds,'RangesList')
    #for each line  
    maxwinds<-mclapply(mc.cores=10,names(cg),function(acc){
      #now get the total in each window
      v<-unlist(viewSums(Views(cg[[acc]]$pos,winds)[seqlevels(winds)]))
      z<-!v==0
      v1<-sapply(split(x=v[z],windregs[z]),max)
      v1<-v1[as.character(regnums)]
      names(v1)<-as.character(regnums)
      v1[is.na(v1)]<-0
      #now for negative strand
      v<- unlist(viewSums(Views(cg[[acc]]$neg,winds)[seqlevels(winds)]))#sum strands
      z<-!v==0
      v2<-sapply(split(x=v,windregs),max)
      v2<-v2[as.character(regnums)]
      names(v2)<-as.character(regnums)
      v2[is.na(v2)]<-0
      #and get the max window for each
      pmax(v1,v2) 
    })
    names(maxwinds)<-names(cg)
    maxwinds
  })
})

#now melt into a dataframe
w.df<-melt(window.sum.distributions)
h(w.df)
colnames(w.df)<-c('value','acc','windowsize','region')

v=accession.df$library.size
names(v)<-as.character(accession.df$accession)
w.df$val.libnorm<- w.df$value/ v[w.df$acc]

v=accession.df$intergenic.lib.size
names(v)<-as.character(accession.df$accession)
w.df$val.intlibnorm<- w.df$value/ v[w.df$acc]

v=accession.df$tss.lib.size
names(v)<-as.character(accession.df$accession)
w.df$val.tsslibnorm<- w.df$value/ v[w.df$acc]


w.df$nonzero<-as.numeric(w.df$value>0)
accs<-names(cage.tag.rles)[c(1,5,10,15,20,25,30,35,40,45,50)]
w.df.bak<-w.df

w.df<-w.df[w.df$acc%in%accs,]



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



# now with techreps -------------------------------------------------------
mline<-multilines[3]
accs<-accession[grepl(mline,accession)]
accs<-accs[!grepl('reseq',accs)]
w.df<-w.df.bak

ggplot(w.df,aes(x=w.df[w.df$acc==accs[1],]$value,y=w.df[w.df$acc==accs[2],]$value))+geom_point()
#which window size maximizes agreement between replicates?
lapply(unique(w.df$windowsize)function(w){
  
  
  
})


# now do density plots ----------------------------------------------------
df<-melt(rapply(window.sum.distributions,how='replace',function(x){v<-(as.vector(x));v[v<quantile(v,0.999)]}))
h(df)
df<-df[df$L2==w & df$L3 %in% accession.df$accession[1:5],]

#now produce density plots for the lines, showing distribution of counts for the inactive, negative and postiive regions
ggplot(df,aes(x=value,color=L3,facet=L1))+geom_density()+facet_grid(.~L1)+scale_x_log10()#+scale_y_log10()

df<-melt(rapply(window.sum.distribution ,how='replace',function(x){v<-(as.vector(x));v[v<quantile(v,0.999)]}))
df<-df[df$L2==w & df$L3 %in% accession.df$accession[1:10],]
#now produce density plots for the lines, showing distribution of counts for the inactive, negative and postiive regions
ggplot(df,aes(x=value,color=L3,facet=L1))+geom_density()+facet_grid(.~L1)+scale_x_log10()#+scale_y_log10()


#we want to see if dividing by library size normalizes the distributions. (I suspect not)


# Now do discriminance analysis -------------------------------------------
#supposing we have a our cutoff set for each 



