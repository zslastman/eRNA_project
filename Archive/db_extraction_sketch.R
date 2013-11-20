# 
dbase = dbConnect(MySQL(), user='furlong', password='gnolruf', host='base3', dbname='genomedb2')
# dbListTables(dbase)
# dbColumnInfo(dbase,'transcript_6')
# dbColumnInfo(dbase,'genesynonyms_6')
# dbColumnInfo(dbase,'accessionset')
# dbColumnInfo(dbase,'geneaccession_6')
# dbColumnInfo(dbase,'transcript_model_feature_6')
# dbColumnInfo(dbase,'gene_6')

#   
#   
  this gives numeric ids and symbols
   syn  = dbGetQuery(dbase, "select id, gene_id, syn from genesynonyms_6 ")
   h(syn,n=100)

#go to transcript model features, get the mapping from FBgn to transcript
gene  = dbGetQuery(dbase, "select id, acc from gene_6 ")
colnames(gene)<- c('id','acc')
gene<-gene[order(gene$id),]
gid2FBgn<-gene$acc
#   
#   #this gives FBgns and CGs
#   accset  = dbGetQuery(dbase, "select id, name, description from accessionset ")
#   h(accset,n=100)
#   
#   #this gives gene_ids and accessions - which are eitehr CGs or FBGns
#   gacc  = dbGetQuery(dbase, "select id, gene_id, acc from geneaccession_6 ")
#   h(gacc,n=100)
#   


#go to transcript model features, get the mapping from FBgn to transcript
tscpt_mod  = dbGetQuery(dbase, "select transcript_id, acc from transcript_model_feature_6 ")
colnames(tscpt_mod)<- c('transcript_id','acc')
tscpt_mod<-  tscpt_mod[grepl('FBgn',tscpt_mod$acc,fixed=T),]#only fbgn rows
tscpt_mod$acc<- gsub('.*(FBgn\\d+).*','\\1',  tscpt_mod$acc)#edit them down to just the FBgns
tscpt_mod<-tscpt_mod[! duplicated(tscpt_mod$transcript_id),]#one for each combo
tscpt_mod<-tscpt_mod[!duplicated(tscpt_mod,MARGIN=1),]
sapply(unique(tscpt_mod$transcript_id),function(tid){
  tid2FBgns<-by(simplify=F,tscpt_mod[1:50,],tscpt_mod$transcript_id[1:50],function(x)x$acc)
  #this is a one to one map bar a single transcript.
  tid2FBgns<-tid2FBgns[length(tid2FBgns)==1]#just get rid of the single multimap
  