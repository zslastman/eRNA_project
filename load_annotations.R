###This script loads the information on our crms for the CAGE project
source('src/tss_cage_functions.R')
message('loading Gene, transcript, CRM data ')
setwd(dir="/g/furlong/Harnett/TSS_CAGE_myfolder/")

intergenic_dist=1000
# constants ---------------------------------------------------------------

sl<-seqlengths(Dmelanogaster)


################################################################################################
# load transcript data from sql database ----------------------------------
if(   all(sapply(c(file.transcripts,file.tssgr,file.gene.model.data),FUN=function(f){file.exists(f)}))  ){
  sapply(c(file.transcripts,file.tssgr,file.gene.model.data,file.synonyms),function(f){load(f,envir=.GlobalEnv)})
  message('Transcript files already present, loaded')
}else{
  #Load from the labs databse
  library(RMySQL)
  # library(RPostgreSQL, lib.loc='/g/furlong/Degner/R')

  #first get info on accesions/ids/symbols (synonyms)
  dbase = dbConnect(MySQL(), user='furlong', password='gnolruf', host='gbcs', dbname='genomedb2')  
  synonyms  = dbGetQuery(dbase, "select id, gene_id, syn from genesynonyms_6 ")
  #we can leave out all but the FBgn ids
  #synonyms<-synonyms[grepl(x=synonyms$syn,'FBgn'),]
  save(synonyms,file=file.synonyms)

  #Now load transcripts from the databse
  p.transcripts  = dbGetQuery(dbase, "select c.name, t.start, t.stop + 1, t.acc, gene_id, t.strand id from transcript_6 t, chromosome_6 c where t.strand = '+' AND t.chromosome_id = c.id")
  n.transcripts  = dbGetQuery(dbase, "select c.name, t.start, t.stop + 1, t.acc, gene_id, t.strand id from transcript_6 t, chromosome_6 c where t.strand = '-' AND t.chromosome_id = c.id")
  colnames(p.transcripts) <-c('Chr', 'S1', 'S2', 'TrID', 'Gene', 'S')
  colnames(n.transcripts) <- c('Chr', 'S1', 'S2', 'TrID', 'Gene', 'S')
 #actually better not to put FBGN   over id, as the FBgn accession aren't very reliable over time.
#   p.transcripts$FBgn<- as.vector(gid2FBgn[ as.numeric(p.transcripts$Gene) ])
#   n.transcripts$FBgn<- as.vector(gid2FBgn[ as.numeric(n.transcripts$Gene) ])
  #now into a GRange object 
  transcripts  = rbind(p.transcripts, n.transcripts)
  transcripts.gr<-GRanges(paste0('chr',transcripts[,1]),
                       IRanges(transcripts[,2]+1,transcripts[,3]+1),#note the indexes are off by 1 in the db
                       strand=Rle(transcripts[,'S']),
                       TrID=transcripts[,'TrID'],
                       Gene=transcripts[,'Gene'],
                       seqinfo=seqinfo(Dmelanogaster)
  )
  #remove duplicates, unassigned chromosomes etc.
  transcripts.gr<-keepSeqlevels(transcripts.gr,chrs.keep)
  transcripts.gr<-transcripts.gr[!duplicated(transcripts.gr)]
  # transcripts.gr<-transcripts.gr[seqnames(transcripts.gr)%in% bigchrs] #not sure why this is here...
  transcripts.gr<- sort(transcripts.gr)
  # #Load gene Activity Calls -----------------------------------------------
  active68.gr<- as(import('data/ActiveGenesPrediction_68h.bed'),'GRanges')
  active46.gr<- as(import('data/ActiveGenesPrediction_46h.bed'),'GRanges')
  transcripts.gr$active68<-countOverlaps(transcripts.gr,active68.gr)>0
  transcripts.gr$active46<-countOverlaps(transcripts.gr,active46.gr)>0

  #create tss gr
  tss.gr<-transcripts.gr
  end(tss.gr)[as.logical(transcripts.gr@strand=='+')]<-start(tss.gr)[as.logical(transcripts.gr@strand=='+')]
  start(tss.gr)[as.logical(transcripts.gr@strand=='-')]<-end(tss.gr)[as.logical(transcripts.gr@strand=='-')]
  tss.gr<-tss.gr[!duplicated(tss.gr)]


  #dbListTables(dbase)
  #dbListFields(dbase,'gene_model_feature_6')
  #Get the names of the chrs from a seperate table
  chrnamelist = dbGetQuery(dbase, "select g.name ,g.id from chromosome_6 g") 
  chrnamelist = as.vector( chrnamelist[,1][order(chrnamelist$id)])
  chrnamelist = paste0('chr', chrnamelist)
  #now the genes
  gene.model.data  = dbGetQuery(dbase, "select g.id, g.gene_id, g.ftype, g.chromosome_id, g.start, g.stop, g.strand from gene_model_feature_6 g")
  colnames(gene.model.data)<-c('id','gene_id','ftype','chromosome_id','start','stop','strand')
  #edit chrnames
  gene.model.data$chromosome_id = chrnamelist[as.numeric(as.character(gene.model.data$chromosome_id))]
  gene.model.data = gene.model.data[gene.model.data$chromosome_id %in% seqnames(si),]
  #   To check data
#   'FBgn0003900'
#   'RAFBtr0071953'
#   transcripts[ transcripts$TrID=='FBtr0071953' ,]
#   gene.model.data[which(gene.model.data$gene_id=='1866'),]
#   gene.model.data[which(gene.model.data$gene_id=='4359'),]
  
  gene.model.data <- GRanges(seqnames = gene.model.data$chromosome_id,
                             IRanges(gene.model.data$start,gene.model.data[,'stop']),
                             strand=gene.model.data$strand ,
                             id=gene.model.data$id ,
                             ftype=gene.model.data$ftype,
                             gene_id=gene.model.data$gene_id,
                             seqinfo=si)
  
  
  #Now let's just get the genes as a gr
  genes.gr = dbGetQuery(dbase, "select g.id, g.acc, g.name,g.symbol,g.chromosome_id,g.start,g.stop,g.strand from gene_6 g")
  genes.gr$chromosome_id = chrnamelist[as.numeric(as.character(genes.gr$chromosome_id))]
  genes.gr = genes.gr[genes.gr$chromosome_id %in% seqnames(si),]
  genes.gr <- GRanges(seqnames=genes.gr$chromosome_id,
                             IRanges(genes.gr$start,genes.gr[,'stop']),
                             strand=genes.gr$strand,
                             id=genes.gr$id,
                             acc=genes.gr$acc,
                             name=genes.gr$name,
                            symbol=genes.gr$symbol,
                             seqinfo=si)

  #add the info from the gene to the transcripts

  mcols(transcripts.gr) = cbind(mcols(transcripts.gr),genes.gr[match(transcripts.gr$Gene,genes.gr$id),c('name','symbol','acc')])
  
  save(transcripts.gr,tss.gr,gene.model.data,genes.gr,file=file.gene.annotation)
    

}

# Also Load the lincRNA data from Young et al -----------------------------
lincRNA.df<-read.table('data/young_etal_2012/Young_etal_2012.txt',header=T,sep='\t')
lincRNA.df<-lincRNA.df[ grepl(':',lincRNA.df$Position) ,]
gregexpr(pattern='(\\w+):()()',as.character(lincRNA.df[,2]))
tmp<-strsplit(as.character(lincRNA.df[,2]),split=':')
chrs<-sapply(tmp,function(s){paste0('chr',s[[1]])})
pos<-sapply(tmp,function(s){ as.numeric(unlist(strsplit(s[[2]],'-')) ) })
lincRNA.gr<-GRanges(chrs,IRanges(pos[1,],pos[2,]))
#mcols(lincRNA.gr)<-lincRNA.df
lincRNA.gr$name<-paste0('lincRNA_',1:length(lincRNA.gr))
export(lincRNA.gr,'data/lincRNAs.youngetal.bed')

#concatenate our transcripts and lincRNAs.
trancripts.lincs.gr<-do.call('c',sapply(list(lincRNA.gr,transcripts.gr),  function(gr){ mcols(gr) <-NULL;gr}))



################################################################################################
# load crm binary data ----------------------------------------------------
message('load the crm_binary data')
crm_list <- read.delim(comment="#",file.path(file.crm.list.binary))
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
                  seqinfo=si)
crm8008.gr<-sort(crm8008.gr)

# Define Intergenic CRMs
crm8008.gr$intergenic <- distanceToNearest(crm8008.gr,trancripts.lincs.gr)$distance>intergenic_dist

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

message('8008 data loaded')




####################################################################################
# Load CAD data -----------------------------------------------------------
#load the file I've formatted
cad3.df<-read.delim('data/CAD3.txt',comment='',header=T,sep='\t',quote='',dec='.',)
for(col in colnames(cad3.df)){
  cad3.df[,col][cad3.df[,col]=='ni']<-NA
}

#make a granges object
cad3.gr<-GRanges(paste0('chr',cad3.df$Chr),IRanges(cad3.df$Start,cad3.df$Stop),seqinfo=si,name=cad3.df$Name)
                      
#add in the dataframe to the Granges, and delete it
for(c in colnames(cad3.df)[!colnames(cad3.df)%in%c('Start','Stop','Name','Chr')] ){mcols(cad3.gr)[[c]]<-cad3.df[[c]]}
rm(cad3.df)
cad3.gr=sort(cad3.gr)

#fix the names, the have a lot of weird characters in them
cad3.gr$name<-gsub('[#\\(\\)\\s/\\-]','',cad3.gr$name)


#do gene assingment - first get indices in synonym list, looking up the FBgns
#format the FBGn column
cad3.gr$FBGn<-sapply(cad3.gr$FBGn,function(s){
  s=gsub(x=s,'\\s','')#remove whitespace
  if(grepl(',',s)){
    strsplit(s,',')[[1]]#split on commas
  }else{s}
})
inds<-sapply(cad3.gr$FBGn,function(FBgn){which(synonyms$syn == FBgn)[1]})
#matching our FBgns to gene ids.
cad3.gr$g.id<-synonyms$gene_id[inds]
#match gene ids to transcripts
cad3.gr$transcripts<-sapply(cad3.gr$g.id,function(g.id){
    if(is.na(g.id)){return(NA)
    }else{
    transcripts.gr$TrID[ transcripts.gr$Gene==g.id]}
 })
#Now let's mark those cad3 enhancers active in mesoderm at 10 or 11
cad3.gr$meso68<-!cad3.gr$M10 %in% c(0,NA) | !cad3.gr$M11 %in% c(0,NA)
cad3.gr$other68<-!cad3.gr$O10 %in% c(0,NA) | !cad3.gr$O11 %in% c(0,NA)
#fix some colnames
colnames(mcols(cad3.gr))<-gsub(x=colnames(mcols(cad3.gr)),pattern='(.*)\\.(.*)',rep='\\1\\2')
#define crms inactive at tmepoints
cad3.gr$inactive24<- (cad3.gr$O5==0 & cad3.gr$M5==0 & cad3.gr$O6==0 & cad3.gr$M6==0 & cad3.gr$O7==0 & cad3.gr$M7==0) %in% T
cad3.gr$active24<-cad3.gr$O5%in%c(1,'S,V') | cad3.gr$M5%in%c(1,'S,V')  | cad3.gr$O6%in%c(1,'S,V') | cad3.gr$M6%in%c(1,'S,V')  | cad3.gr$O7%in%c(1,'S,V') | cad3.gr$M7%in%c(1,'S,V') 
#
cad3.gr$inactive68<- (cad3.gr$O10==0 & cad3.gr$O11==0 & cad3.gr$M10==0 & cad3.gr$M11==0) %in% T
cad3.gr$active68<-cad3.gr$O10%in%c(1,'S,V') | cad3.gr$O11%in%c(1,'S,V') | cad3.gr$M10%in%c(1,'S,V')  | cad3.gr$M11%in%c(1,'S,V') 
#
cad3.gr$inactive1012<- cad3.gr$O13%in%0 & cad3.gr$M13%in%0 & cad3.gr$O14%in%0 & cad3.gr$M14%in%0 & cad3.gr$O15 %in%0 & cad3.gr$M15%in%0
cad3.gr$active1012<-cad3.gr$O13%in%c(1,'S,V') | cad3.gr$M13%in%c(1,'S,V')  | cad3.gr$O14%in%c(1,'S,V') | cad3.gr$M14%in%c(1,'S,V')  | cad3.gr$O15%in%c(1,'S,V') | cad3.gr$M15%in%c(1,'S,V') 
# Define Intergenic CRMs
cad3.gr$intergenic <-  distanceToNearest(cad3.gr,trancripts.lincs.gr)$distance>intergenic_dist
#catagory stating if te enhancer is one of guillames
cad3.gr$GJ<-grepl('GJ\\d+',cad3.gr$name)
#now export a bed file with our cad3 enhancers and names
export(cad3.gr,'analysis/cad3.bed')


#put a unique gene number on each tss and cad object
#now load the table linking fbgns and fbtrs
fbgn_fbtr_fbpp.df <- read.delim(header=F, comment.char='#',stringsAsFactors=F,sep='\t','/g/furlong/Harnett/TSS_CAGE_myfolder/data/fbgn_fbtr_fbpp_fb_2013_06.tsv.gz' )
colnames(fbgn_fbtr_fbpp.df) <- c('FBgn','FBtr','FBpp')
trs = tss.gr$TrID


save(cad3.gr,crm8008.gr,file=file.crm.annotations)



################################################################################
#load phastcons
read.table('/g/furlong/Harnett/TSS_CAGE_myfolder/data/phastcons.sga.gz' )
phastconstable<-phastconstable[,-2]
