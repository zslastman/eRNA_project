
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/generate_Rle.R')
source('src/generate_chromdata_cage.R')
source('src/load_annotations.R')
#We need our tss information, 
#our crm information, 
#our pollII peak information with summits.
#later we can just use view maxs to find the max of poll II where we don't have peaks




#we want a heatmap with the cage signal, centered on our pollII (and teh reverse?)
#and we want a heatmap with the poll signal, centered on our cage


#get pol peaks with summit location 
# pol.peaks<-read.delim('data/chromatin/with_summits/PolII_B-PolII-Rpb3-B51_6-8h_peaks.txt',comment.char='#')

# Polymerase summits from the MACS summits --------------------------------

pol_summits<- chrompeaks[['PolII_6-8h']]#start with the peak
start(pol_summits)<-pol_summits$summit+start(pol_summits)#shift to the summit
end(pol_summits)<-start(pol_summits)



# polymerase summits from the local maximum -------------------------------
pol_maxs<-combinegrs(list(tss=tss.gr,crm=crm8008.gr))

find_polmax<-function(gr,polrle,size=500){
  gr<-resize(gr,size,fix='center')
  gr<-gr[width(gr)==size]
  v<-Views(polrle,as(gr,'RangesList'))[unique(seqnames(gr))]
  m<-do.call(rbind,(lapply(v,function(x){suppressWarnings(as.matrix(x))})))
  m<-t(apply(m,1,function(x)smooth.spline(x,spar=0.4)$y))
  #  plot(m[sample(1:8007,1),],type='l')#spar4 looks about right
  #   maxs<-apply(m,1,function(x){#function returning NA if the max isn't 2x as high as the mean
  #     w=which.max(x)
  #     if(all(x==0)){return(NA)}
  #     if(  (x[w]/mean(x,na.rm=T)) < 2 ){return(NA)
  #     }else{return(w)}
  #     
  #     })
  maxs=apply(m,1,which.max)
  start(gr)=start(gr)+maxs-1
  end(gr)=start(gr)
  return(gr)
}

pol_maxs<-find_polmax(pol_maxs,chrom.rles.rpgc.sub.merge[['PolII_6.8']])


for(pol.peak.gr in list(pol_summits,pol_maxs)){

  pol.peak.gr<-resize(pol.peak.gr,width=500,fix='center')#now expand to a 500bp window around the summit
  pol.peak.gr=pol.peak.gr[width(pol.peak.gr)==500]#eliminate those near the edges of chrs so they are a uniform length
  pol.peak.gr=sort(pol.peak.gr)
  
#get polymerase signal
pol.peak.gr$polsum<-unlist(viewSums(GRViews(chrom.rles.rpgc.sub.merge[['PolII_6.8']],pol.peak.gr)))




#annotate overlap with crm or tss
pol.peak.gr$tss <- 0 < countOverlaps(pol.peak.gr,tss.gr,maxgap=500)
pol.peak.gr$crm <- 0 < countOverlaps(pol.peak.gr,crm8008.gr,maxgap=500)

pol.peak.gr$name<-paste0('pol_summit',1:length(pol.peak.gr))
export(pol.peak.gr,'polsummits.tmp.bed')






#this takes in a grange object and a list iwth positive/negative/both srls. It puts the matrix of 
gr.add.cagemats<-function(gr,rles){
	#get the CAGE tag info
	stopifnot(length(unique(width(gr)))==1)#all should be the same width
	strandviews<-sapply(rles,function(strand){Views(strand,as(gr,'RangesList'))[unique(seqnames(gr))]})
	#summaries over our regions
	strandsums<-sapply(strandviews,function(strand){unlist(viewSums(strand))})
	gr$bothsum<-strandsums[,'both']
	#strandmaxs<-sapply(strandviews,function(strand){unlist(viewMaxs(strand))})
	#get matrices over our windows
	gr$plusmat<-do.call(rbind,(lapply(strandviews[[1]],function(x){suppressWarnings(as.matrix(x))})))
	gr$negmat<- do.call(rbind,(lapply(strandviews[[2]],function(x){suppressWarnings(as.matrix(x))})))
	#defining dominant strand - have a boolean and a score.
	gr$plusisgreater<-apply(strandsums,1,function(x)x[1]>=x[2])
	gr$strandscore<-apply(strandsums,1,function(x)(x[1]-x[2])/(x[1]+x[2]))
	return(gr)
}
  #include only those tss which overlap a pol peak(redudant for pol_summits)
pol.peak.gr<-pol.peak.gr[ pol.peak.gr$crm | pol.peak.gr$active68]
#add the cage info to the pol summits object
pol.peak.gr<- gr.add.cagemats(pol.peak.gr,alltags.rpgc)
dim(pol.peak.gr$negmat)
#get matrices with the dominant strand
pol.peak.gr$dommat<-pol.peak.gr$negmat
pol.peak.gr$dommat[pol.peak.gr$plusisgreater,]<-pol.peak.gr$plusmat[pol.peak.gr$plusisgreater,]
pol.peak.gr$dommat[!pol.peak.gr$plusisgreater,]<-t(apply(pol.peak.gr[!pol.peak.gr$plusisgreater]$dommat,1,rev))

  
  
  
# 
# #do clustering seperately from the heatmap
# # cluster rows
# hc.rows <- hclust(dist(m))
# plot(hc.rows)
# # transpose the matrix and cluster columns
# hc.cols <- hclust(dist(t(mtscaled)))
# # draw heatmap for first cluster
# heatmap(mtscaled[cutree(hc.rows,k=2)==1,], Colv=as.dendrogram(hc.cols), scale='none')
# 
# #this will give you a tall, thin plot - just use mfrow, easier than screwing with margins
# par(mfrow=c(1,2),mar=c(2,1,3,1))
# 
# #plot all of our pol peaks
# m<-pol.peak.gr$dommat[1:50,]
# cage.heatmap(m)
# plot(1:3)#since we're using mfrow



#now plot the crms and the tss side by side
m.tss<-pol.peak.gr$dommat[pol.peak.gr$tss,]
m.crm<-pol.peak.gr$dommat[pol.peak.gr$crm & ! pol.peak.gr$tss,]

par(mfrow=c(1,2),mar=c(4,1,3,1))

#produce heatmaps of two sets
cage.heatmap<-function(m,tit,maxpoint=F){
  
  if(maxpoint){
    m<-t(apply(m,1,function(x){
    z=rep(0,length(x))
    z[which.max(x)]<-1
    z
    }))
  }
  
  heatmap.simple(m[order(rowSums(m),decreasing=T),])

  title(main=tit,xlab='bp from polymerase summit')
  axis(side=1,at=seq(0,1,l=11),labels=seq(-250,250,by=50),cex=0.5,)
  abline(v=0.5,lty=3,lwd=3)

}

cage.heatmap(m.tss[,],tit='CAGE tag distribution all TSS Summits \n Mesoderm 6-8 hrs')
cage.heatmap(m.crm[,],tit='CAGE tag distribution all CRM Summits \n Mesoderm 6-8 hrs')

cage.heatmap(m.tss[,],tit='Maximum CAGE tag site - all TSS Summits \n Mesoderm 6-8 hrs',maxpoint=T)
cage.heatmap(m.crm[,],tit='Maximum CAGE tag site - all CRM Summits \n Mesoderm 6-8 hrs',maxpoint=T)


#produce density plot (averagogram)
ggplot(data=melt(list(m.tss[,],m.crm[,])),aes(x=Var2,y=value,color=factor(L1)))+
  scale_color_manual(name='Type',values=c('red','blue'), labels=c('tss','crm'))+
  stat_smooth(stat_summary='mean_cl_boot')+
#  stat_smooth(stat_summary='mean_cl_boot',data=melt(m.crm[1:100,]))+
  scale_y_continuous(name='Mean RPGC Tag count')+
  scale_x_discrete(name='Distance from Pol Summit',labels=seq(-250,250,by=50),breaks=floor(seq(1,500,l=11)))+
  ggtitle('density of CAGE reads centered on Polymerase Summit')
  



}
#
heatmap(m.tss,Colv=NA,Rowv=NULL,labCol=F,labRow=F,col=rgb((1:100)/100,0,0))
heatmap(m.crm,Colv=NA,Rowv=NULL,labCol=F,labRow=F,col=rgb((1:100)/100,0,0))



#exlude low-scoring poll peaks
#exclude those which aren't clearly stranded


heatmap(m.max,Colv=NA,Rowv=NULL,labCol=F,labRow=F,col=rgb((1:100)/100,0,0))


#now make a matrix of zeroes with ones at these points


#now we can do heatmaps of the cage signal 

#for those which are above our cutoff.

#for just genes

#for just crm8008.gr


gr <- granges()

#Plot a gene if the gene region has any coverage
gene_plot <- if(length(width(gr)) >= 1){
#exons contains the start and stop regions for each exon in the gene
exons <- data.frame(xmin = exon_starts, xmax = exon_ends, ymin = 0, ymax = Inf)

cov_plot <- autoplot (
gr, 
stat = "coverage",
geom = "area",
colour = "black",
fill = "green",
xlab = "Position",
ylab = "Coverage",
main = "Coverage across gene"

#Marks the exon regions out so that you can see how much of the sequencing is on target	
) + geom_rect(
data = exons, 
aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
fill = "yellow", 
colour = "black", alpha = 0.5

)
}

#plot the graph	
gene_plot



#Calculate the point of maximum cage singal
#e.g. by defining the point where one basepair accounts for at least 50% of the signal.
