source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt",lib='~/Harnett/R')

listMarts()
bm=useMart('ensembl')
listDatasets(bm)

dm3_refGene_txdb <- makeTranscriptDbFromBiomart(biomart="bio",
                                                dataset="dmelanogaster_gene_ensembl",
                                                id_prefix="ensembl_",
                                                host="www.biomart.org",
                                                port=80)
dm3_refGene_txdb.bak = dm3_refGene_txdb


cols(dm3_refGene_txdb)
select(dm3_refGene_txdb, keys = keys, cols=cols(dm3_refGene_txdb), keytype="GENEID")


symbol(ttrack)

#put in the symbol names to our GeneRegionTrack
library(biomaRt)
mart = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
newsymb=getGene(id = symbol(ttrack), type = "refseq_mrna", mart = mart)
newsymb.v=newsymb$external_gene_id
names(newsymb.v)=newsymb$refseq_mrna
symbol(ttrack) = newsymb.v[symbol(ttrack)]
displayPars(ttrack)$showId=T




# 
# ######Attempting to build a transcript track....
# library(GenomicFeatures)
# dm3_refGene_txdb <- makeTranscriptDbFromBiomart(biomart="ensembl",
#                         dataset="dmelanogaster_gene_ensembl",
#                         id_prefix="ensembl_",
#                         host="www.biomart.org",
#                         port=80)
# 
# 
# transcripts(dm3_refGene_txdb)<-transcripts(dm3_refGene_txdb)
# dm3_refGene_txdb.bak=dm3_refGene_txdb
# 
# dm3_refGene_txdb<-dm3_refGene_txdb.bak
# 
# keepSeqlevels(dm3_refGene_txdb,bigchrs)
# dm3_refGene_txdb$.chrom<-paste0('chr',dm3_refGene_txdb$.chrom)
# makeTranscriptDb(transcripts=data.frame(tx_id=1:length(transcripts.gr),
#                               tx_name=transcripts.gr$TrID,
#                               tx_chrom=as.vector(seqnames(transcripts.gr)),
#                               tx_strand=as.vector(strand(transcripts.gr)),
#                               tx_start=start(transcripts.gr),
#                               tx_end=end(transcripts.gr)),
#                  genes=data.frame(gene_id=transcripts.gr$symbol)
#                  )
# transcripts.gr$symbol <- fbtr2symb.v[as.character(transcripts.gr$TrID)]
# with(mcols(transcripts.gr),{TrID})
#      
# seqinfo(dm3_refGene_txdb)
# chrs.keep(dm3_refGene_txdb,bigchrs)
# ttrack<-GeneRegionTrack(dm3_refGene_txdb,showId=T,geneSymbols=T)