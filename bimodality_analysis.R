# bimodality_analysis.R
# Script to look at bimodality in our cage data
#  author:Dermot Harnett      
#  Requires load_cage,load_tagseq,load_annotation,load_chromatin

setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )     
#A run scripts, functions etc.
source ( '/g/furlong/Harnett/TSS_CAGE_myfolder/src/tss_cage_functions.R' )
source('src/load_chromdata.R')
#Load libraries
##load data
load('/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/allcage.object.R')

load('data/objects/allcage.object.R')
load('data/objects/gene.annotation.object.R')

outfolder = 'analysis/bimodality_analysis'
dir.create(outfolder,showWarnings=F)


#############################################################################
#TSS cage score by strand
#############################################################################
tss.gr = tss.gr[!duplicated(tss.gr)]
# posscores = GRViewsums(resize(tss.gr,width=500,fix='center'),rle=allcage$tp68hembryo$pos)
# negscores = GRViewsums(resize(tss.gr,width=500,fix='center'),rle=allcage$tp68hembryo$neg)
posscores = GRViewsums(resize(tss.gr,width=500,fix='center'),rle=allcage$tp68hembryo$pos)
negscores = GRViewsums(resize(tss.gr,width=500,fix='center'),rle=allcage$tp68hembryo$neg)
posgreater = posscores > negscores
domscores = ifelse(posgreater, posscores, negscores)
minorscores = ifelse(!posgreater, posscores, negscores)
logcagetotal=log10(domscores+minorscores+1)
#logcagetotal=logcagetotal[is.finite(logcagetotal)]
for(i in 1:10){ try({cagecutoff=getBimodalCutoff(logcagetotal,2,lower = -1,upper=10)} ) }

posscores.meso = GRViewsums(resize(tss.gr,width=500,fix='center'),rle=allcage$tp68hmeso$pos)
negscores.meso = GRViewsums(resize(tss.gr,width=500,fix='center'),rle=allcage$tp68hmeso$neg)
posgreater.meso = posscores.meso > negscores.meso
domscores.meso = ifelse(posgreater.meso, posscores.meso, negscores.meso)
minorscores.meso = ifelse(!posgreater.meso, posscores.meso, negscores.meso)
logcagetotal.meso=log10(domscores.meso+minorscores.meso+1)

# Density plot
plot(density(logcagetotal,na.rm=TRUE, data=dataName, legend=T, xlab="xLabel", ylab="yLabel", main="Total Cage Tag density on TSS - Embryo 68h"))
plot(density(logcagetotal.meso,na.rm=TRUE, data=dataName, legend=T, xlab="xLabel", ylab="yLabel", main="Total Cage Tag density on TSS - Meso 68h"))
plot(density(domscores.meso,na.rm=TRUE, data=dataName, legend=T, xlab="xLabel", ylab="yLabel", main="Total Cage Tag density on TSS - Dominant strand 68h meso"))
plot(density(negscores.meso,na.rm=TRUE, data=dataName, legend=T, xlab="xLabel", ylab="yLabel", main="Total Cage Tag density on TSS - weak strand 68h meso"))


abline(v=cagecutoff)
#Now use a bimodality cutoff on the total expression to pick expressed and nonexpressed genes.
activegenes = domscores+minorscores > cagecutoff



#############################################################################
#load our chromatin information
#############################################################################
#columns describind overlap with peak
tss.gr$H3K4me3  <-overlapsAny(tss.gr,list(chrompeaks.modencode[['K4me3_4-8h']], chrompeaks[['K4me3_4-6h']] ,chrompeaks[['K4me3_6-8h']]),maxgap=50) 
tss.gr$H3K4me1  <-overlapsAny(tss.gr,list(chrompeaks.modencode[['K4me1_4-8h']], chrompeaks[['K4me1_4-6h']] ,chrompeaks[['K4me1_6-8h']]),maxgap=50)
tss.gr$H3K27ac  <-overlapsAny(tss.gr,list(chrompeaks.modencode[['K27ac_4-8h']], chrompeaks[['K27Ac_4-6h']] ,chrompeaks[['K27Ac_6-8h']]),maxgap=50)
tss.gr$H3K79me3 <-overlapsAny(tss.gr,list(chrompeaks[['K79me3_4-6h']],chrompeaks[['K79me3_6-8h']]),maxgap=50)
tss.gr$H3K36me3 <-overlapsAny(tss.gr,list(chrompeaks[['K36me3_4-6h']],chrompeaks[['K36me3_6-8h']]),maxgap=50)
tss.gr$polII    <-overlapsAny(tss.gr,list(chrompeaks[['PolII_4-6h']],chrompeaks[['PolII_6-8h']]),maxgap=50)
#matrix of chromatin data
tss.chrom.mat<-sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(tss.gr,'RangesList'))))})
tss.mod.chrom.mat<-sapply(chrom.rles.modencode,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(tss.gr,'RangesList'))))})





#############################################################################
#Graphs showing Data on H3K4me3 peak present
#############################################################################
dir.create('analysis/bimodality_analysis')
#Scatterplot showing both strands of (bimodality)
pdf(paste0(outfolder,'/','Dom_vs_Min.heatscatter.pdf'))
heatscatter(cor=F,x= log10( domscores ) ,y= log10( minorscores ) ,log='',main='Heatmap of Main vs Minor Strand for all TSS at 68 Hours\n Whole Embryo Cage',xlab = 'Log10 Dominant Strand Cage' , ylab = 'Log10 Minor Strand Cage')
dev.off()



pdf(paste0(outfolder,'/','Dom_vs_Min.heatscatter.pdf'))
heatscatter(cor=F,x= log10( domscores.meso ) ,y= log10( minorscores.meso ) ,log='',main='Heatmap of Main vs Minor Strand for all TSS at 68 Hours\n Whole Mesodermal Cage',xlab = 'Log10 Dominant Strand Cage' , ylab = 'Log10 Minor Strand Cage')
dev.off()

mean( log10( domscores.meso ) > 2 & log10( minorscores.meso ) >2 ) * 100
mean( log10( domscores ) > 3 & log10( minorscores ) >3 ) * 100



#Density plot asking if we have seperate population of bimodal genes and splitting if so
# Density plot
# Save plot as pdf
pdf(file="bimodality_density_plot.pdf")
qplot(geom=bimodality,na.rm=TRUE,ggdf,legend=T, xlab="Ratio of Dominant to Negative Strand", ylab="Density", main="Ratio of Dominant to Negative Strand\n For Expressed Genes")
dev.off()

#barplot showing numbers with/without peak
# Save plot as pdf
pdf(file="Chromatin_Barplots.pdf")
plotHere
dev.off()


#Scatterplot showing expression vs H3K4me3


ggdf = data.frame(x=rnorm(100),a=sample(letters[1:2],rep=T,100))

qplot(stat='bin',data=ggdf,x=a,fill=a)+stat_bin(geom='text',aes(label=..count..),vjust=5)+theme_bw()

qplot(stat='bin',data=ggdf,x=a,fill=a)+stat_bin(geom='text',aes(label=..count..),vjust=5)+theme_minimal()

qplot(stat='bin',data=ggdf,x=a,fill=a)+stat_bin(geom='text',aes(label=..count..),vjust=5)
