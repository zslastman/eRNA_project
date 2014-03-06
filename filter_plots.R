#This script will analyze the effects of the various filters we're applying to our genes for the 
#analysis of expression variability
#This script integrates the various datasets on genes, crms, chroamtin, rnaseq, cage, tagseq and dnase
#it outputs a list of GRange objects with the data stored as vectors and matrices in the meta data columns
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
source('src/load_chromdata.R')
load('data/objects/scored.annotation.object.R')
load('data/objects/cg.object.R')
load('data/objects/cg.mapfilt.object.R')
load('data/objects/cg.pl.object.R')
load('data/objects/cg.mapfilt.pl.object.R')
load('data/objects/ts.object.R')
load('data/objects/ts.pl.object.R')
load('data/objects/rna.seq.object.R')
load('data/objects/synonyms.object.R')
#Apply mappability filter - FInd out which genes change a lot under mappability filter
#So genes whose total changes by more than 5% 
nofilttot = rowSums(genes.gr$cg.pl)
filttot   = rowSums(genes.gr$cg.mapfilt.pl)
genes.gr$filter.changed = (nofilttot - filttot) / nofilttot 
#This ends up being very high
mean(genes.gr$filter.changed > 0.1, na.rm = T)
#try number of reasonably high tagged sites 


insitu.stages = c('stage1-3','stage4-6','stage7-8','stage9-10','stage11-12','stage13-16')
##Apply BDGP ubiquitous/non ubiquitous filter
#read bdgp data
insitu.df = read.delim(header=F,sep='\t','data/BDGP/insitu_annot.csv')
insitu.df2 = read.delim(header=T,sep='\t','data/BDGP/bdgp.txt')


colnames(insitu.df) = c('symbol','name','FBgn','stage','tissue')
insitu.df$stage = insitu.stages[insitu.df$stage]
dim(insitu.df)
head(insitu.df)

#first, for each of our genes, get all synonyms, then find these in the in situe tabl
load('data/objects/synonyms.object.R')
symbolmatch = insitu.df$symbol %in% synonyms$syn
insitu.df$gene_id[symbolmatch] = synonyms$gene_id[ match(insitu.df$symbol[symbolmatch],synonyms$syn) ]
#shows how we hardly gain anything by also matching FBgn or name...
mean(insitu.df$FBgn[!symbolmatch] %in% synonyms$syn)
mean(insitu.df$name[!symbolmatch] %in% synonyms$syn)
genes.gr$insitu.annotated = genes.gr$id %in% insitu.df$gene_id
sum(genes.gr$insitu.annotated)
#get rid of the annotation with no gene_id
insitu.df = insitu.df [ symbolmatch, ]
#now to each gene attach a table of expression data
genes.gr$bdgpanno=lapply(genes.gr$id,function(x){ insitu.df[insitu.df$gene_id == x , c('stage','tissue') ]})
genes.gr$nostain = sapply(genes.gr$bdgpanno,function(x){all(grepl(x$tissue,pat='no_staining',fixed=T))})#this returns true if the df is empty as well
genes.gr$insitu.annotated = genes.gr$insitu.annotated &  ! genes.gr$nostain
sum(genes.gr$insitu.annotated )
#how many genes have only ubiquitous as their tissue annotation?
genes.gr$just.ub=F
genes.gr$just.ub[genes.gr$insitu.annotated] = sapply(genes.gr$bdgpanno[genes.gr$insitu.annotated],function(x){
  all(grepl(x$tissue,pat = 'ubiquitous|maternal'))
})
genes.gr$any.ub=F
genes.gr$any.ub[genes.gr$insitu.annotated] = sapply(genes.gr$bdgpanno[genes.gr$insitu.annotated],function(x){
  any(grepl(x$tissue,pat = 'ubiquitous',fixed=T))
})
sum(genes.gr$just.ub)
sum(is.na(genes.gr$just.ub))
sum(genes.gr$any.ub)
sum(is.na(genes.gr$any.ub))

#get Davids germlayer annotation
germterms=readLines('data/BDGP/germlayers_noOverlappingTerms.txt')
germterms= strsplit(germterms,split=';|,')
names(germterms) =sapply(germterms,'[[',1) 
germterms = sapply(germterms,function(x){x[-1]})
germterms = melt(germterms)
germterms$gene_id = synonyms$gene_id[match(germterms$value,synonyms$syn)]
#and 
ov.germterms=readLines('data/BDGP/germlayers_wOverlappingTerms.txt')
ov.germterms= strsplit(ov.germterms,split=';|,')
names(ov.germterms) =sapply(ov.germterms,'[[',1) 
ov.germterms = sapply(ov.germterms,function(x){x[-1]})
ov.germterms = melt(ov.germterms)
ov.germterms$gene_id = synonyms$gene_id[match(ov.germterms$value,synonyms$syn)]

genes.gr$endoderm = with(germterms,genes.gr$id %in% gene_id[L1=='endoderm'])
genes.gr$ectoderm = with(germterms,genes.gr$id %in% gene_id[L1=='ectoderm'])
genes.gr$mesoderm = with(germterms,genes.gr$id %in% gene_id[L1=='mesoderm'])
genes.gr$endoderm.ov = with(ov.germterms,genes.gr$id %in% gene_id[L1=='endoderm'])
genes.gr$ectoderm.ov = with(ov.germterms,genes.gr$id %in% gene_id[L1=='ectoderm'])
genes.gr$mesoderm.ov = with(ov.germterms,genes.gr$id %in% gene_id[L1=='mesoderm'])
with(germterms,genes.gr$insitu.annotated[genes.gr$id %in% gene_id] <- T)

sum(genes.gr$bdgp.germlayer.nov == 'mesoderm',na.rm=T)
sum(genes.gr$bdgp.germlayer.nov == 'ectoderm',na.rm=T)
sum(genes.gr$bdgp.germlayer.nov == 'endoderm',na.rm=T)
sum(genes.gr$bdgp.germlayer.ov == 'mesoderm',na.rm=T)
sum(genes.gr$bdgp.germlayer.ov == 'ectoderm',na.rm=T)
sum(genes.gr$bdgp.germlayer.ov == 'endoderm',na.rm=T)


sum(genes.gr$endoderm , na.rm=T)
sum(genes.gr$ectoderm , na.rm=T)
sum(genes.gr$mesoderm , na.rm=T)
sum(genes.gr$endoderm.ov , na.rm=T)
sum(genes.gr$ectoderm.ov , na.rm=T)
sum(genes.gr$mesoderm.ov , na.rm=T)

library(VennDiagram)
#annotated,ub,meso,ecto,endo
setlist = list(
  Annotated=genes.gr$insitu.annotated,
  Ubiquitous=genes.gr$just.ub,
  Mesodermal=genes.gr$mesoderm.ov,
  Ectoodermal=genes.gr$ectoderm.ov,
  Endodermal=genes.gr$endoderm.ov
 )
setlist = sapply(setlist,FUN=which)
#Venn diagrams showing size of expression catagories
venn.diagram(setlist,'filter_Venn_Diagrams.tiff',fill=rainbow(length(setlist)))
venn.diagram(setlist[-1],'filter_Venn_Diagrams_4.tiff',fill=rainbow(length(setlist[-1])))


save(genes.gr,tss.gr,transcripts.gr,file='data/objects/scored.annotation.R')
