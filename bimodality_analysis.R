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
posscores = unlist(viewSums(GRViews(gr=tss.gr,allcage$tp68embryo$pos)))
negscores = unlist(viewSums(GRViews(gr=tss.gr,allcage$tp68embryo$neg)))
posgreater = posscores > negscores
domscores = ifelse(posgreater, posscores, negscores)
minorscores = ifelse(!posgreater, posscores, negscores)
getBim

#Now use a bimodality cutoff on the total expression to pick expressed and nonexpressed genes.


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
tss.chrom.mat<-sapply(chrom.rles.rpgc.sub.merge,function(chrom.rle){unlist(viewMeans(Views(chrom.rle,as(crm8008.gr,'RangesList'))))})

#############################################################################
#Graphs showing Data on H3K4me3 peak present
#############################################################################
dir.create('analysis/bimodality_analysis')
#Scatterplot showing both strands of (bimodality)
pdf(paste0(outfolder,'/','Dom_vs_Min.heatscatter.pdf'))
heatscatter(x= log10( domscores ) ,y= log10( minorscores ) ,log='',main='Heatmap of Main vs Minor Strand for all TSS at 68 Hours\n Whole Embryo Cage',xlab = 'Log10 Dominant Strand Cage' , ylab = 'Log10 Negative Strand Cage')
dev.off()

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
