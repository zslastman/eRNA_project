#script to load motif data, integrate it with the Dnase sensitivity data, and then the gene data
#A Load libraries,functions,data
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
load( 'data/objects/synonyms.object.R' )
windowsize=20000
#get the relevant files


#First get the names of all the motif files
motif.files=list.files('/g/furlong/Harnett/TSS_CAGE_myfolder/data/scanPWMoutput/p002t3/',full.name=T,pattern='output\\.clean\\.bed$')
#only those files which are reasonably big
motif.files = Filter(function(x)file.info(x)$size > 200 ,motif.files)
motif.names=list.files('/g/furlong/Harnett/TSS_CAGE_myfolder/data/scanPWMoutput/p002t3/',full.name=F,pattern='output\\.clean\\.bed$')
motif.names=gsub(x=motif.names,pattern='.pfm.*',rep='\\1')
#get the FBgn numbers and the names
FBgns=gsub(x=motif.names,pattern='.*(FBgn\\d+).*',rep='\\1')
nofbinds=!grepl(x=FBgns,'FBgn')
FBgns[nofbinds]<-NA
motif.names=gsub(x=motif.names,pattern='(.*)_FBgn\\d+.*',rep='\\1')
MAnames = motif.names[grepl(x=motif.names,'MA\\d+')]
#MAnames = gsub(x=MAnames,pat='(MA\\d+)\\..*',rep='\\1')
###MA00001 means that the pwm is from Jaspar, have to translate these
jasparnames = read.delim(header=F,'data/jaspar_names.txt')
colnames(jasparnames) = c('MAnum','symbol')
MAnames = as.character(jasparnames$symbol[match( MAnames, jasparnames$MAnum )])
motif.names[grepl(x=motif.names,'MA\\d+')] <- MAnames


#now read our files, using bedtools and the dnase peaks to subset each one
export(dnase.peaks,con='dnase.peaks.bed')
#read the motif files and get only those motifs within dnase peaks
motifs.gr=sapply( motif.files[],function(f){
  cat('.')
  #use bedtools to subset bed  
  system( paste0('intersectBed -a ',f,' -b dnase.peaks.bed',' >tmp2.bed') ) 
  gr=import('tmp2.bed',format='bed',asRangedData=F) #import as GRanges
  strand(gr) <- gr$name
  gr$name=NULL
  seqlevels(gr)=seqlevels(si)
  seqinfo(gr)=si
  gr
#  import('tmp.bed',format='bed','asRangesList'=F)#import as GRanges
})
names(motifs.gr)<-motif.names[]

#Rles with the position score at the start of the motif (shouldn't make much difference
#Over the 50kb window anyway)
motif.rles=sapply( motifs.gr,function(gr){
  grthin=resize(gr,width=1,fix='start')
  grthin=grthin[order(gr$score,decreasing=T)]
  grthin=grthin[!duplicated(gr)]
  # gr=gr[strand(gr)=='+']
  score=round(grthin$score,digits=6)
  cov=coverage(grthin,weight=score)#this is doing something weird with the score values....
})

#Now threshold the motif rles
####PUT SOMETHING HERE

#Now 
cutoff=0.001#the coverage method seems to create rounding errors, so a cutoff is needed
motif.rles.bin = sapply( motif.rles,function(srle){srle[srle>cutoff]<-1;srle[srle<cutoff]<-0;srle})


save(motif.rles.bin,'data/objects/motif.rles.bin.object.R')
save(motif.rles,'data/objects/motif.rles.bin.object.R')


#now get a matrix, with columns for each PWM and rows for each tss
sapply(grlist,function(gr){
  sapply(motif.rles.bin,function(motif){
    grwindow = resize(grwindow,width=windowsize,fix='center')



  })
})


reglist = list(genes.gr,tss.gr)

#shortened data works okay
motif.windowcount=sapply(simplify=F,reglist,function(region){
  # cat(region)
  # region=get(region)
  region$id=1:length(region)
  suppressWarnings({wideregion = resize(region,width=windowsize,fix='center')})
  wideregion= sort (wideregion)
  m=getStrandedExprMat(wideregion,motif.rles.bin)  
  expect_is(m,'matrix')
  rownames(m)=wideregion$id
  m[as.character(region$id),]
  })
})


#Let's look at our mef2 peaks and see how they relate to motif strength
#get mef2peaks








#I can now score the genes with respect to the number of motifs near them above some threshold do one of two things


#bug: some of the values in motif.rles aren't in the score column of the Granges object - I think this is just numerical









#now create a list of motif Rle objects
motifs.gr[[1]]->gr


#why these low values, what's wrong with these coverage vectors?
cov=motif.rles[[1]][[1]]
covals=runValue(cov)
unique(covals) %in% gr$score
unique(covals[covals>0.001]) %in% gr$score
countOverlaps(resize(gr,width=1,fix='start'),resize(gr,width=1,fix='start')
grthin=resize(gr,width=1,fix='center')
length(grthin)
length(grthin[!duplicated(grthin)])
#how does the duplication method work on GRanges???
sum( duplicated(start(grthin)) && duplicated(seqnames(grthin)) && duplicated(strand(grthin))       )
sum( duplicated(paste0(start(grthin),seqnames(grthin))) )
sum(duplicated(grthin))
gr[seqnames(gr)==chr]

min(runValue(cov[[1]]))

min(gr$score)

save(motifs.gr,file='data/objects/motifs.gr.object.R')

str(motifs.gr)
length(motifs.gr)
names(motifs.gr)

str(motif.rles)
length(motif.rles)
names(motif.rles)
motif.rles[[1]]






#now for each of our 

#We can get the motif content around a gene by expanding it and counting overlaps with each motif set
#Or we could make a load of Rles and use views...



#Plots at the end of this

#Distribution of motif number near genes - genes in different expr groups
#Fall-off of motif number as function of distance
#Degree to which motif strength + Dnase predicts e.g. Mef2 occupancy



#nearby mef2 occupancy vs 
#Does the number of nearby motifs effect trans variance more than cis?


export(GRanges('a',IRanges(c(1,10,20),c(5,17,30))),format='bed',con='tmp1.bed')
export(GRanges('a',IRanges(c(9),c(22))),format='bed',con='tmp2.bed')



GRanges('a,',IRanges(c(9),c(22)))



MAnames[which(is.na(match( MAnames, jasparnames$MAnum)))]


 jasparnames$name 

head(jasparnames)
#Granges list
motif.files[1]->mfile




tab = read.delim(header=F,mfile,stringsAsFactors=F)
tab$V2<-as.numeric(tab$V2)
tab$V3<-as.numeric(tab$V3)
tab=tab[tab$V4 %in% c('+','-'),]
tab=tab[tab$V1 %in% chrs.keep,]
expect_true(mean( tab$V2 <tab$V3,na.rm=T ) > 0.99)
GRanges(tab$V1,IRanges(as.numeric(tab$V2),as.numeric(tab$V3)),strand=tab$V4,score=tab$V5,seqinfo=si)


readLines(mfile)[which(as.numeric(tab$V2)-as.numeric(tab$V3) != -6)]

as.numeric(tab$V3)-as.numeric(tab$V2)
import(mfile)

#Now combine
#Now filter for score
#Filter out those not in open chromatin
#Now attribute a score to each gene using the motifs, let's say the number within 


####Also load the genetic interaction data from Thomas Horn et al

