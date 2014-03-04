#this script uses the fimo motif calls (at least for now) to look at motif occurences within tss and crms
load('data/objects/gene.annotation.object.R')
load('data/objects/crm.annotation.object.R')
load('data/objects/motifs.gr.object.R')

motifs.gr$name = 

motnames=unique(motifs.gr$name)

motname='TATA-Box'

#First - what percentage of promotors have various motifs
threshs = 10^( -(seq(2,6,length.out=10)) )
sapply(threshs,function(thresh){
  mots=motifs.gr[motifs.gr$name==motname & motifs.gr$pvalue < thresh]
  mean(countOverlaps(resize(tss.gr,width=500,fix='center'),mots))
})

sapply(motnames,function(motname){
  sapply(threshs,function(thresh){
    mots=motifs.gr[motifs.gr$name==motname & motifs.gr$pvalue < thresh]
    mean(countOverlaps(resize(crm8008.gr,width=500,fix='center'),mots))
  })
})

thresh=0.01
motif.overlap=sapply(simplify=F,c('tss.gr','crm8008.gr','cad3.gr'),function(gr){
  gr=get(gr)
  ov=sapply(motnames,function(motname){
      mots=motifs.gr[motifs.gr$name==motname & motifs.gr$pvalue < thresh]
      countOverlaps(resize(gr,width=500,fix='center'),mots)
  })
})
mean(apply(motif.overlap$tss.gr,1,any))
mean(apply(motif.overlap$crm8008.gr,1,any))

canyprom.motif = 
mcols(crm8008.gr)[,]

#Print Bed files for our motifs
dir.create('analysis/motifbeds')
for(motname in motnames){
  mots=motifs.gr[motifs.gr$name==motname]
  export(mots,con=paste0('analysis/motifbeds/',motname,'.bed'))
}




save(motif.overlap,file='data/objects/motif.overlap.object.R')