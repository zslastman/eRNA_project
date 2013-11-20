

hg19_refGene_txdb <- makeTranscriptDbFromUCSC(genome="hg19", tablename="refGene")
ttrack.h<-GeneRegionTrack(hg19_refGene_txdb)
gtrack.h<-GenomeAxisTrack()
itrack.h<-IdeogramTrack(genome='hg19')

#74509725

totRNA.kimetal.h1<-import('data/kim_etal_2010_eRNA_mm/GSE21161_RAW/GSM530211_hr1c+.bw',asRangedData=F)
max(totRNA.kimetal.h1$score)

totRNAp_h0<-DataTrack(import('data/kim_etal_2010_eRNA_mm/GSE21161_RAW/GSM530210_hr0c+.bw',asRangedData=F), genome = "hg19", type = "l",
           name = "Coverage", window = -1, chromosome = "chr15",
           importFunction = function(file){import(con=file,asRangedData=F)}, stream = TRUE)
totRNAp_h1<-DataTrack(import('data/kim_etal_2010_eRNA_mm/GSE21161_RAW/GSM530211_hr1c+.bw',asRangedData=F), genome = "hg19", type = "l",
                      name = "Coverage", window = -1, chromosome = "chr15",
                      importFunction = function(file){import(con=file,asRangedData=F)}, stream = TRUE)

save(list=c('totRNAp_h1','totRNAp_h0','ttrack.h','gtrack.h','itrack.h'),file='data/objects/human.gviz.objects.R')


windowsize=10000
ch='chr15'
n=74509725

jpeg(paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/plots/Gviz/',main.tit,'.h.jpeg'),height=1200,width=1920)

plotTracks(trackList=list(totRNAp_h1,totRNAp_h0,ttrack.h,gtrack.h),type='p',chr=ch,from=n-windowsize,to=n+windowsize)

dev.off()