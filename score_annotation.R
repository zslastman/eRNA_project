#This script integrates the various datasets on genes, crms, chroamtin, rnaseq, cage, tagseq and dnase
#it outputs a list of GRange objects with the data stored as vectors and matrices in the meta data columns
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
source('src/load_chromdata.R')
load('data/objects/gene.annotation.object.R')
load('data/objects/cg.object.R')
load('data/objects/cg.mapfilt.object.R')
load('data/objects/cg.pl.object.R')
load('data/objects/cg.mapfilt.pl.object.R')
load('data/objects/ts.object.R')
load('data/objects/ts.pl.object.R')
w=500#window size for determining CRM scores

################################################################################
#1 Get Expression Scores for the RNAseq, Cage and Tagseq datasets
getStrandedExprMat <- function(gr,SRLs){
	#function that takes in some SimpleRleLists (e.g. RNAseq data)
	#which have a 'positive' and 'negative' sublist structure
	#and then gets uses the GRViews function to get the strand specific
	#score for each of the input regions in each of the SRLs 
	# browser()
	  require(testthat)
	  # browser()
	  simplify2array(mclapply(mc.cores=10,SRLs,function(SRL){
	  	  if(identical(names(SRL)[1:2],c('pos','neg'))) {
		  v=rep(0,length(gr))
		  ispos  = as.vector(strand(gr)=='+')
		  isneg  = as.vector(strand(gr)=='-')
		  isstar = as.vector(strand(gr)=='*')
		  v[ isneg ] = unlist(viewSums(GRViews(SRL$neg,gr[ isneg ])))
		  v[ ispos ] = unlist(viewSums(GRViews(SRL$pos,gr[ ispos ] )))
		  v[ isstar ] = unlist(viewSums(GRViews(SRL$pos,gr[ isstar ] ))) +
		     unlist(viewSums(GRViews(SRL$neg,gr[ isstar ] ))) 
		  v
		  }else{
			unlist(viewSums(GRViews(SRL,gr)))
		  }
	  }))
}

#The below is so obviously repetitive that it should clearly be looped.
#But fuck it.

#first, the transcripts - we face a problem when widening our grange object
#doing this upsets the sort order, and thus the output of GRViews Becomes wrong
#
transcripts.gr=resize(transcripts.gr,width=width(transcripts.gr)+250,fix='center')
transcripts.gr=sort(transcripts.gr)
transcripts.gr$cg = getStrandedExprMat(transcripts.gr,cg)
transcripts.gr$cg.pl = 	getStrandedExprMat(transcripts.gr,cg.pl)
transcripts.gr$cg.mapfilt = getStrandedExprMat(transcripts.gr,cg.mapfilt)
transcripts.gr$cg.mapfilt.pl = getStrandedExprMat(transcripts.gr,cg.mapfilt.pl)
transcripts.gr$ts = getStrandedExprMat(transcripts.gr,ts)
transcripts.gr$ts.pl = getStrandedExprMat(transcripts.gr,cts.pl)
transcripts.gr$rna.seq = getStrandedExprMat(transcripts.gr,rna.seq)
transcripts.gr$c.rna.seq = getStrandedExprMat(transcripts.gr,c.rna.seq)
transcripts.gr$z.rna.seq = getStrandedExprMat(transcripts.gr,z.rna.seq)

maxrange=which.max(transcripts.gr$cg[,1])
transcripts.gr$cg[maxrange,1]
maxg=transcripts.gr[maxrange,1]
sum(cg[[1]][[1]][[as.character(seqnames(maxg))]][start(maxg) :end(maxg) ])
max(transcripts.gr$cg[,1])





#now individual TSS
tss.gr.wide=resize(tss.gr,width=width(tss.gr)+250,fix='center')
tss.gr$cg = getStrandedExprMat(tss.gr.wide,cg)
tss.gr$cg.pl = getStrandedExprMat(tss.gr.wide,cg.pl)
tss.gr$cg.mapfilt = getStrandedExprMat(tss.gr.wide,cg.mapfilt)
tss.gr$cg.mapfilt.pl = getStrandedExprMat(tss.gr.wide,cg.mapfilt.pl)
tss.gr$ts = getStrandedExprMat(tss.gr.wide,ts)
tss.gr$ts.pl = getStrandedExprMat(tss.gr.wide,cts.pl)
tss.gr$rna.seq = getStrandedExprMat(tss.gr.wide,rna.seq)
tss.gr$c.rna.seq = getStrandedExprMat(tss.gr.wide,c.rna.seq)
tss.gr$z.rna.seq = getStrandedExprMat(tss.gr.wide,z.rna.seq)

#####Now genes
genes.gr=resize(genes.gr,width=width(genes.gr)+250,fix='center')
genes.gr=sort(genes.gr)
genes.gr$cg = getStrandedExprMat(genes.gr,cg)
genes.gr$cg.pl = 	getStrandedExprMat(genes.gr,cg.pl)
genes.gr$cg.mapfilt = getStrandedExprMat(genes.gr,cg.mapfilt)
genes.gr$cg.mapfilt.pl = getStrandedExprMat(genes.gr,cg.mapfilt.pl)
genes.gr$ts = getStrandedExprMat(genes.gr,ts)
tmp = getStrandedExprMat(genes.gr,ts.pl[1:2])
genes.gr$ts.pl = getStrandedExprMat(genes.gr,ts.pl)

genes.gr$rna.seq = getStrandedExprMat(genes.gr,rna.seq)
genes.gr$c.rna.seq = getStrandedExprMat(genes.gr,c.rna.seq)
genes.gr$z.rna.seq = getStrandedExprMat(genes.gr,z.rna.seq)

#annotate those transcripts which are the only transcript of for their gene
tsnumtable = table(transcripts.gr$Gene)
single.tr.genes = as.integer(names(tsnumtable)[tsnumtable==1])
transcripts.gr$only.tr = transcripts.gr$Gene %in% single.tr.genes

#### Annotate the crms with expression data, but using the 'getBestWindow' function
crm8008.gr.wide=resize(crm8008.gr,width=width(crm8008.gr)+250,fix='center')
crm8008.gr$cg = getBestWindowMat(crm8008.gr.wide,cg)
crm8008.gr$cg.pl = getBestWindowMat(crm8008.gr.wide,cg.pl)
crm8008.gr$ts = getBestWindowMat(crm8008.gr.wide,ts)
crm8008.gr$ts.pl = getBestWindowMat(crm8008.gr.wide,cts.pl)
crm8008.gr$rna.seq = getBestWindowMat(crm8008.gr.wide,rna.seq)
crm8008.gr$c.rna.seq = getBestWindowMat(crm8008.gr.wide,c.rna.seq)
crm8008.gr$z.rna.seq = getBestWindowMat(crm8008.gr.wide,z.rna.seq)

#and the CAD enhancer
cad3.gr.wide=resize(cad3.gr,width=width(cad3.gr)+250,fix='center')
cad3.gr$cg = getBestWindowMat(cad3.gr.wide,cg)
cad3.gr$cg.pl = getBestWindowMat(cad3.gr.wide,cg.pl)
cad3.gr$ts = getBestWindowMat(cad3.gr.wide,ts)
cad3.gr$ts.pl = getBestWindowMat(cad3.gr.wide,cts.pl)
cad3.gr$rna.seq = getBestWindowMat(cad3.gr.wide,rna.seq)
cad3.gr$c.rna.seq = getBestWindowMat(cad3.gr.wide,c.rna.seq)
cad3.gr$z.rna.seq = getBestWindowMat(cad3.gr.wide,z.rna.seq)



####TODO - turn the above into loops - as below
gene.annotations = c('transcripts.gr','tss.gr','genes.gr')[1:2]
crm.annotations = c('crm8008.gr','cad3.gr')
expr.datasets = c('cg','cg.pl','cg.mapfilt','cg.mapfilt.pl','ts','ts.pl','rna.seq','c.rna.seq','z.rna.seq')[1:2]


expr.mats=sapply(simplify=F,gene.annotations,function(region){
  # browser()
  region=get(region)[1:10]
	sapply(simplify=F,expr.datasets,function(data){
		suppressWarnings({wideregion = resize(transcripts.gr,width=width(transcripts.gr)+250,fix='center')})
		wideregion=sort(wideregion)
		m=getStrandedExprMat(wideregion,get(data)[1:2])
		rownames(m)=wideregion$id
    m[as.character(region$id),]
  })
})

Q
sapply(simplify=F,expr.datasets,function(data){
	suppressWarnings({wideregion = resize(transcripts.gr,width=width(transcripts.gr)+250,fix='center')})
	mcols(transcripts.gr)[[data]]<<-getStrandedExprMat(wideregion,get(data))
    NULL
})


tmp<-sapply(simplify=F,gene.annotations,function(region){
  regionget(region)
  sapply(simplify=F,expr.datasets,function(data){
  	suppressWarnings({wideregion = resize(region,width=width(region)+250,fix='center')})
    mcols(crm8008.gr)[[data]]<<-getStrandedExprMat(wideregion,get(data))
    NULL
  })
})
tmp<-sapply(simplify=F,crm.annotations,function(region){
  sapply(simplify=F,expr.datasets,function(data){
    suppressWarnings({wideregion = resize(region,width=width(region)+250,fix='center')})
    mcols(crm8008.gr)[[data]]<<-getBestWindowMat(wideregion,data)
    NULL
  })
})

genes.gr$geneDensity = getGeneDensity(genes.gr,w=200000)


save(transcripts.gr,genes.gr,tss.gr,file='data/objects/scored.annotation.object.R')
load('data/objects/scored.annotation.object.R')