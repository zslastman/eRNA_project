


###############################################################################
#Plots comparing individual features to our variances
###############################################################################
load('data/objects/scored.annotation.R')
load('data/objects/tfgrlist.object.R')
load('data/objects/motif.overlap.object.R')#get motif overlap
load('data/objects/motifs.gr.object.R')
load('data/objects/dnase.peaks.object.R')



# functions for boxplot labelling
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
mean.n <- function(x){
  return(c(y = median(x)*0.97, label = round(mean(x),2))) 
  # experiment with the multiplier to find the perfect position
}

genetic.vars = list(
  'Proportional' = (var.all$Cis[,tp] + var.all$Trans[,tp]) / var.all$Noise[,tp],
  'Cis' = var.all$Cis[,tp] ,
  'Trans' = var.all$Trans[,tp] ,
  'Genetic' = var.all$Cis[,tp] + var.all$Trans[,tp],
  'Non developmental' = var.all$Cis[,tp] + var.all$Trans[,tp] + var.all$Noise[,tp],
  'Limma' = apply(fit$coeff,1,sd)#limma only uses replicated ones anyway...
  # 'Limma_replicated_only' = apply(fit$coeff[,repped_coeffs],1,sd)
)
windowsizes = c('50kb'=50000,'25kb'=25000,'12.5kb'=12500)

gn = 'Limma'
w=windowsizes[[2]]

for(gn in names(genetic.vars)){
genetic.var = genetic.vars[[gn]]


pdf(paste0('analysis/variance_rankings/indiv_features_vs_Var_',gn,' w=',w,'_.pdf'))
# pdf('analysis/variance_rankings/density_vs_variance_heatscatter.pdf')
h5sumgenes.gr$genedensity50 = getGeneDensity(h5sumgenes.gr,genes=genes.gr,windowsize=w)
heatscatter(x=h5sumgenes.gr$genedensity50,y=log10(genetic.var),main=paste0('Gene Density vs ',gn,' windowsize ',w),xlab='gene density',ylab=paste0('Log10',gn))
# dev.off()
# Nearby TF peaks
h5sumgenes.gr$tfpeakdensity = rowSums(sapply( tfgrlist[['6-8']] , function(tfpeaks.gr){
  getGeneDensity(h5sumgenes.gr,genes=tfpeaks.gr,windowsize=w)
}))
# pdf('analysis/variance_rankings/tfpeakdensity_vs_variance_heatscatter.pdf')
heatscatter(x=h5sumgenes.gr$tfpeakdensity,y=log10(genetic.var),main=paste0('Local big5 TF peak density vs ',gn,' windowsize ',w),xlab='tf peak density',method='s',ylab=gn)
# dev.off()
# Nearby Dnase Peaks
# pdf('analysis/variance_rankings/tfpeakdensity_vs_variance_heatscatter.pdf')
# load('data/objects/dnase.peaks.object.R')
h5sumgenes.gr$dnasepeakdensity = getGeneDensity(h5sumgenes.gr,genes=dnase.peaks,windowsize=w)
heatscatter(x=h5sumgenes.gr$dnasepeakdensity,y=log10(genetic.var),main=paste0('Local Dnase Peak Density vs ',gn,' windowsize ',w),xlab='Dnase Peak density',method='s',ylab=gn)
# dev.off()
dev.off()

allmotifs.gr<-import('/g/furlong/Harnett/TSS_CAGE_myfolder/data/fimo_out_JASPAR_CORE_insects_02/fimo.gff',asRangedData=F)
allmotifs.gr = allmotifs.gr[countOverlaps(allmotifs.gr,dnase.peaks)>0]

# Promotor Type
gene.motif.overlap = matrix(0,ncol=ncol(motif.overlap$tss.gr),nrow=length(genes.gr))
colnames(gene.motif.overlap) = colnames(motif.overlap$tss.gr)
rownames(gene.motif.overlap) = genes.gr$id
# the below code sums all tss for the gene
for(i in length(tss.gr$Gene)){
  tss.ov = motif.overlap$tss.gr[i,]#get the motifs overlapping that gene
  gid = tss.gr$Gene[i]#get the corresponding gene
  gene.motif.overlap[gid,]=gene.motif.overlap[gid,] + tss.ov
}

pdf(paste0('analysis/variance_rankings/motifs_vs_Var_',gn,' w=',w,'_.pdf'))

#or just calculate overlap directly on the peaks....
motifs=unique(motifs.gr$name)
thresh=0.0001
h5sumgenes.gr$motif.overlap=sapply(motifs,function(motif){
      mots=motifs.gr[motifs.gr$name==motif & motifs.gr$pvalue < thresh]
      countOverlaps(resize(h5sumgenes.gr,width=500,fix='center'),mots)
})
for(motif in motifs){
h5sumgenes.gr$mot.overlap = h5sumgenes.gr$motif.overlap[ ,motif ] > 0 
# qplot(geom='errorbar',log='y', y = genetic.var, x = h5sumgenes.gr$TATA.overlap , fill= h5sumgenes.gr$TATA.overlap )
d=data.frame(motif=factor(h5sumgenes.gr$mot.overlap) , 'GeneticVariance' = genetic.var)
print(ggplot(d, aes(factor(motif),GeneticVariance,color=motif)) +
  stat_summary(fun.data='mean_cl_boot',geom='errorbar',size=4) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
  ggtitle(paste0(motif,' presence vs ',gn)))
  scale_y_continuous(limits=c(0,2))
}

dev.off()



#could also skip tss with no expression, or just use the main one for each gene.
#now show boxplots for various sets.
expr.pattern = rep(NA,length(genes.gr))
expr.pattern[genes.gr$endoderm] = 'endoderm'
expr.pattern[genes.gr$mesoderm] = 'mesoderm'
expr.pattern[genes.gr$ectoderm] = 'ectoderm'
expr.pattern[genes.gr$just.ub] = 'ubiquitious'
h5sumgenes.gr$expr.pattern = expr.pattern[match(h5sumgenes.gr$Gene,genes.gr$id)]

pdf(paste0('analysis/variance_rankings/indivexpr_features_vs_Var_',gn,' w=',w,'_.pdf'))

# Expression catagory (from BDGP)
qplot(log = 'y', main = 'expression pattern vs. Genetic variability' ,
 x = h5sumgenes.gr$expr.pattern, fill = h5sumgenes.gr$expr.pattern, y = genetic.var, geom = 'Violin' )

expr.tab=c(table(h5sumgenes.gr$expr.pattern),'NA'=sum(is.na(h5sumgenes.gr$expr.pattern)))
expr.tab=data.frame(pattern=names(expr.tab),freq=expr.tab)

d=data.frame('Pattern'=factor(h5sumgenes.gr$expr.pattern) , 'GeneticVariance' = genetic.var)

print(
ggplot(d, aes(factor(Pattern),GeneticVariance)) +
  stat_summary(data=d,fun.data='mean_cl_boot',geom='errorbar',size=4, colour = c('yellow','red','orange','blue','black')) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) 
)

dev.off()

}
#Finally, take the number of motifs (any )

#Now let's output files for Davids analysis suite.

for(gn in names(genetic.vars)){
  genetic.var = genetic.vars[[gn]]

  hasID = !is.na(h5sumgenes.gr$Gene)
  IDs = genes.gr[match(h5sumgenes.gr$Gene[hasID],genes.gr$id)]$acc
  var.GOtab.df = data.frame(gene_ID = IDs, score = genetic.var[hasID])
  write.table(var.GOtab.df[hasID,],file=paste0('analysis/GO_BDGP_analysis/go_var_table_',gn,'_w',w,'_.txt'),sep='\t',col.names=T,row.names=F,quote=F)
}
