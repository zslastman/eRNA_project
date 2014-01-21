###This script loads the information on our crms for the CAGE project

message('loading Gene, transcript, CRM data ')
setwd(dir="/g/furlong/Harnett/TSS_CAGE_myfolder/")

intergenic_dist=500
# constants ---------------------------------------------------------------

sl<-seqlengths(Dmelanogaster)

crm_list.binary.file<-"data/TF8008.txt"
file.transcripts<-'data/objects/transcripts.object.R'
file.tssgr<-'data/objects/tssgr.R'
file.gene.model.data<-'data/objects/gene.model.data.object.R'


# load transcript data from sql database ----------------------------------
if(   all(sapply(c(file.transcripts,file.tssgr,file.gene.model.data),FUN=function(f){file.exists(f)}))  ){
  sapply(c(file.transcripts,file.tssgr,file.gene.model.data),function(f){load(f,envir=.GlobalEnv)})
  message('Transcript files already present, loaded')
}else{
  
  library(RMySQL,lib.loc='~/Harnett/R')
  library(RPostgreSQL, lib.loc='~/Harnett/R')
  
  dbase = dbConnect(MySQL(), user='furlong', password='gnolruf', host='base3', dbname='genomedb2')
  
  p.transcripts  = dbGetQuery(dbase, "select c.name, t.start, t.stop + 1, t.acc, t.name, t.strand from transcript_6 t, chromosome_6 c where t.strand = '+' AND t.chromosome_id = c.id")
  n.transcripts  = dbGetQuery(dbase, "select c.name, t.start, t.stop + 1, t.acc, t.name, t.strand from transcript_6 t, chromosome_6 c where t.strand = '-' AND t.chromosome_id = c.id")
  
  colnames(p.transcripts) <-c('Chr', 'S1', 'S2', 'TrID', 'Gene', 'S')
  colnames(n.transcripts) <- c('Chr', 'S1', 'S2', 'TrID', 'Gene', 'S')
  transcripts  = rbind(p.transcripts, n.transcripts)
  transcripts<-GRanges(paste0('chr',transcripts[,1]),
                       IRanges(transcripts[,2]+1,transcripts[,3]+1),#note the indexes are off by 1 in the db
                       strand=Rle(transcripts[,'S']),
                       TrID=transcripts[,'TrID'],
                       Gene=transcripts[,'Gene'],
                       seqinfo=seqinfo(Dmelanogaster)
  )
  #dbListTables(dbase)
  dbListFields(dbase,'gene_model_feature_6')
  #this gets us introns, UTRs etc, but currentl doesn't seem to have the chromosome name...
  gene.model.data  = dbGetQuery(dbase, "select g.id, g.gene_id, g.ftype, g.chromosome_id, g.start, g.stop, g.strand from gene_model_feature_6 g")
  colnames(gene.model.data)<-c('id','gene_id','ftype','chromosome_id','start','stop','strand')
  gene.model.data <- GRanges(seqnames=gene.model.data$chromosome_id,
                             IRanges(gene.model.data$start,gene.model.data[,'stop']),
                             strand=gene.model.data$strand,
                             id=gene.model.data$id,
                             ftype=gene.model.data$ftype,
                             gene_id=gene.model.data$gene_id)
  
  #create tss gr
  tss.gr<-transcripts
  end(tss.gr)[as.logical(transcripts@strand=='+')]<-start(tss.gr)[as.logical(transcripts@strand=='+')]
  start(tss.gr)[as.logical(transcripts@strand=='-')]<-end(tss.gr)[as.logical(transcripts@strand=='-')]
  #Checked these in Gbrowse, they're correct
  #tss.gr[300:320]
  
  save(transcripts,file=file.transcripts)
  save(tss.gr,file=file.tssgr)
  save(gene.model.data,file=file.gene.model.data)
  

}

# load crm binary data ----------------------------------------------------

message('load the crm_binary data')
crm_list <- read.delim(comment="#",file.path(crm_list.binary.file))
tfcols<-colnames(crm_list)[grepl("_", colnames(crm_list), fixed=T)]
tfs<-unique(gsub(tfcols,pattern='^(.*?)_.*',replacement='\\1'))

#fix the chr names
crm_list$CHR<-paste0('chr',as.character(crm_list$CHR))

#and the indexing (one of the CRMS starts at 0)
crm_list$START<-crm_list$START+1
crm_list$STOP<-crm_list$STOP+1
crm_list$middle <- crm_list$START + round((crm_list$STOP-crm_list$START)/2)

#sumarrized data for the time points
crm_list$Activity24   <- (apply(crm_list[,colnames(crm_list)[grepl("2.4", colnames(crm_list), fixed=T)]], 1, sum))
crm_list$Activity46   <-(apply(crm_list[,colnames(crm_list)[grepl("4.6", colnames(crm_list), fixed=T)]], 1, sum))
crm_list$Activity68   <- (apply(crm_list[,colnames(crm_list)[grepl("6.8", colnames(crm_list), fixed=T)]], 1, sum))
crm_list$Activity810  <- (apply(crm_list[,colnames(crm_list)[grepl("8.10", colnames(crm_list), fixed=T)]], 1, sum))
crm_list$Activity1012 <- (apply(crm_list[,colnames(crm_list)[grepl("10.12", colnames(crm_list), fixed=T)]], 1, sum))

#summarized data for the tfs
#lapply(tfs,function(tf){
#  apply(crm_list[,colnames(crm_list)[grepl(tf, colnames(crm_list))]],1,function(x){any(as.logical(x))})
#  })

crm_list$bin <- as.logical(apply(crm_list[,colnames(crm_list)[grepl("bin", colnames(crm_list), fixed=T)]], 1, function(x) {if (sum(x)>0) {return(1)}; return(0)}))
crm_list$tin <- as.logical(apply(crm_list[,colnames(crm_list)[grepl("tin", colnames(crm_list), fixed=T)]], 1, function(x) {if (sum(x)>0) {return(1)}; return(0)}))
crm_list$mef2 <- as.logical(apply(crm_list[,colnames(crm_list)[grepl("mef2", colnames(crm_list), fixed=T)]], 1, function(x) {if (sum(x)>0) {return(1)}; return(0)}))
crm_list$twi <- as.logical(apply(crm_list[,colnames(crm_list)[grepl("twi", colnames(crm_list), fixed=T)]], 1, function(x) {if (sum(x)>0) {return(1)}; return(0)}))
crm_list$bap <- as.logical(crm_list$bap_6.8)
crm_list$Ntf <- as.logical(apply(crm_list[,c('bin','tin','twi','mef2','bap')],1,sum))

crm_list$activitycode<-paste0(crm_list$mef2_2.4,crm_list$tin_2.4,crm_list$twi_2.4,'-',
       crm_list$mef2_4.6,crm_list$tin_4.6,crm_list$twi_4.6,'-',
       crm_list$bap_6.8,crm_list$bin_6.8,crm_list$mef2_6.8,crm_list$tin_6.8,crm_list$twi_6.8,'-',
       crm_list$bin_8.10,crm_list$mef2_8.10,'-',
       crm_list$bin_10.12,crm_list$mef2_10.12)

#convert the binary values to logical values, prevents accidents with indexing
for(i in 1:ncol(crm_list)){
  if(all(crm_list[,i] %in% c(0,1))){
    crm_list[,i]<-as.logical(crm_list[,i])  
  }
}

crm_list<-crm_list[order(crm_list$CRMID),]
crm_list<-crm_list[
  crm_list$START>250 & 
    crm_list$STOP < seqlengths(Dmelanogaster)[as.character(crm_list$CHR)] ,]#or end

#make granges object
crm8008.gr <- GRanges(seqnames=Rle(crm_list$CHR),
                  IRanges(crm_list$START,crm_list$STOP),
                  id=crm_list$CRMID,
                  seqinfo=seqinfo(Dmelanogaster))

# Define Intergenic CRMs
crm_list$intergenic <- distanceToNearest(crm8008.gr,transcripts)$distance>intergenic_dist

#check our crmgr object is the right length
stopifnot(nrow(crm_list)==length(crm8008.gr))

mark_rgbs=list(
  "Dnase"=col2rgb('orange'),
  "Mnase"=col2rgb('brown'),
  "H3"=col2rgb('gray'),
  "H3K4me3"=col2rgb('darkorange'),
  "H3K4me1"=col2rgb('orange'),
  "H3K9me2"=col2rgb('palegreen'),
  "H3K9me3"=col2rgb('darkgreen'),
  "H3K27ac"=col2rgb('darkred'),
  "H3K27me3"=col2rgb('blue'),
  "PolII"=col2rgb('yellow'),
  "SBP"=col2rgb('gray'),
  "H3K79me3"=col2rgb('purple'),
  "H3K36me3"=col2rgb('palegreen'),
  "Pho"=col2rgb('darkblue'))
tf_rgbs=list(
  "bap"=col2rgb('purple'),
  "bin"=col2rgb('darkorange'),
  "mef2"=col2rgb('blue'),
  "tin"=col2rgb('red'),
  "twi"=col2rgb('green'))
  
crm_list$itemRgb<-rgb(0.5,0.5,0)

#give our crms a color value for each combination of tfs, by adding the rgb values
#for each crm
for(tf in tfs){
  cols<-crm_list$itemRgb[crm_list[,tf]]
  cols<-lapply(cols,col2rgb)
  cols<-lapply(cols,function(col){(col+tf_rgbs[[tf]])})
  cols<-lapply(cols,function(col){(col)/max(col)})
  crm_list$itemRgb[crm_list[,tf]] <- unlist( lapply(cols,function(color){do.call(rgb,as.list(color))})) 
}
    
crm8008.gr$name<-paste0(crm8008.gr$id,'...',crm_list$activitycode)
crm8008.gr$thickStart<- crm_list$middle-24
crm8008.gr$thickEnd<- crm_list$middle+25
crm8008.gr$score<-crm_list$Activity68
crm8008.gr$itemRgb<-crm_list$itemRgb

f<-'data/crm8008.gr.color.bed'
write(x='track name="ItemRGBDemo" description="Item RGB demonstration" itemRgb="On"',file=f)

export(crm8008.gr,con='data/crm8008.gr.color.bed',format='bed',append=T)
readLines('data/crm8008.gr.color.bed')[1:2]


seqlevels(crm8008.gr)<-seqlevels(si)

message('crm data loaded')


