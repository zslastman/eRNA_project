#

ohlerpwms = '/g/furlong/Harnett/TSS_CAGE_myfolder/data/Ohler/ohler_2002_PWMs.txt'
pwmlines=readLines(ohlerpwms)
matstarts = grep('^A',pwmlines)
namelines = sapply(1:length(pwmlines),function(l){ grepl('^-+$',pwmlines[l-1]) &  grepl('^-+$',pwmlines[l+1])  })
namelines = which(unlist(namelines))+1



numlines = which(grepl('[-\\d]+.*\\s[-\\d+]',pwmlines) & ! grepl('[a-z]',pwmlines))
pwmlines[numlines]
pwmlines[namelines]

motiflist=list()


for(i in 1:length(namelines)){
  cat('.')
  matstart=matstarts[i] + 1
  if(i == length(matstarts)){matend=length(pwmlines)
  }else{ matend=max(numlines[numlines <matstarts[i+1]] ) }
  matname=pwmlines[max(namelines[namelines < matstart])]
  mat=pwmlines[matstart:matend]
  mat
  mat=sapply(USE.NAMES=F,mat,strsplit,'\\s+')
  mat=sapply(USE.NAMES=F,mat,'[',-1)
  mat=apply(mat,2,as.numeric)
  rownames(mat) = c('A|','C|','G|','T|')
  motiflist[[matname]]=mat
}
names(motiflist) = paste0('ohler_',names(motiflist))

motiflist

length(motiflist)
dir.create('data/Ohler/pfmfiles')

for(motif in names(motiflist)){
  write.table(col.names=F,row.names=F,x=motiflist[[motif]],file=paste0('data/Ohler/pfmfiles/',motif,'.pfm'))
}

system('jaspar2meme -pfm data/Ohler/pfmfiles > data/Ohler/ohler.meme')

system('meme2meme ~/Harnett/data_general/motif_databases/JASPAR_POLII_2008.meme Ohler/ohler.meme >Ohler_and_Jaspar_promotor_motifs.meme')

system('')

fimo --max-stored-scores 100000 motif_databases/*POL* genome/dm3allchrs.fa
fimo --max-stored-scores 1000000 --qv-thresh --thresh 0.01 motif_databases/*POL* genome/dm3allchrs.fa


file.info(motifs.file)$size

object.size(motifs.gr)







runFimo<-function(
genomefile = '~/Harnett/data_general/genome/dm3allchrs.fa',#genome file (unzipped)
motiffile = '~/Harnett/data_general/motif_databases/big5.meme', #pwm file
thresh = 0.005,
outdir = '~/Harnett/data_general/motif_databases/fimo_out_big5'){
  #make erase the outdir if it exists already
  system(paste0('rm -r ', outdir) )
  #system(paste0('fimo --max-stored-scores 1000000 --qv-thresh --thresh ',qvalthresh,' --o ',outdir,' ',motiffile,' ',genomefile))
  system(paste0('fimo --max-stored-scores 10000000 --thresh ',thresh,' --o ',outdir,' ',motiffile,' ',genomefile))
  #now import the results as a GRanges object
  motifs.file = paste0(outdir,'/fimo.gff')
  motifs.gr = import(asRangedData=F,motifs.file)
  #Now we have to Map the names retaine by FIMO to the names in the pwm file, which are readable
  #Read the original file
  pwmlines = readLines(motiffile)
  motifs.gr$name = gsub('^[+-]','',motifs.gr$Name )
  motifs.gr$Name=NULL
  motifs.gr$type=NULL
  motifs.gr$source=NULL
  pwmlines = pwmlines[ grepl('MOTIF' , pwmlines  ) ]#Get the lines with names
  polnames = gsub('MOTIF\\s+(.*?)\\s+','',pwmlines)#Make a named vector for the conversion
  names(polnames) = gsub('MOTIF\\s+(\\w*?)\\s+.*','\\1',pwmlines)
  #convert hose which have readable names
  havename=grepl('\\w+',polnames[motifs.gr$name])
  motifs.gr$name[havename]=polnames[motifs.gr$name[havename]]
  return(motifs.gr)
}



#Now verify a threshold by reproducing the figures in zinzen supp fig. 1 
#Get the peaks 
#get the tf data 
tffiles<-list.files('/g/furlong/girardot/meso_chip_migration_dm3/results/export/peaks/bed/',pattern='bed',full.names=T)
tffiles
tffiles.tfnames<-gsub('peaks_(.*?)_.*','\\1',list.files('/g/furlong/girardot/meso_chip_migration_dm3/results/export/peaks/bed/',pattern='bed',full.names=F))
tffiles.tps<-paste0('tp',gsub('.*?_(\\d+)\\-(\\d+).bed','\\1\\2',tffiles))
tffile.df<-data.frame(tfname=tffiles.tfnames,tp=tffiles.tps,file=tffiles)
tfgrlist<-sapply(simplify=F,as.character(unique(tffile.df$tp)),function(tp){
  sapply(simplify=F,as.character(unique(tffile.df$tfname)),function(tfname){
    if(!any( tffile.df$tp==tp & tffile.df$tfname==tfname)){return(NULL)}    
    gr<-import(as.character(tffile.df$file[tffile.df$tp==tp & tffile.df$tfname==tfname]),asRangedData=F)
  })
})
fimo --max-stored-scores 1000000 --o fimo_out_big5 ~/Harnett/data_general/motif_databases/big5.meme  ~/Harnett/data_general/genome/dm3allchrs.fa
#Get the motif data from our 5 pwms
#big5motifs.gr=import('/g/furlong/Harnett/data_general/motif_databases/fimo_out_big5/fimo.gff',asRangedData=F)
big5motifs.gr = runFimo()#use default settings
#And a selection of pvalue thresholds

#make sure the tf names match in the motif and peak objects
stopifnot(identical(sort(unique(tffiles.tfnames)),sort(unique(motifs.gr$name))))

motif=unique(motifs.gr$name)[[1]]
dists=c('100bp'=100,'200bp'=200,'500bp'=500)
dists=c(50,100,150,200,250,300,350,400,450,500)
names(dists)=paste0(as.character(dists),'bp')
pvals=c(10^-3,10^-4,10^-5,10^-6,10^-7)
pvals=10( -seq(3,7,length.out = 12) )
names(pvals)=as.character(pvals)
#And for the 5 motifs
peak.motifs.data=sapply(simplify=F,pvals,function(pval){
  sapply(simplify=F,unique(motifs.gr$name),function(motif){
    #At various bp distances
    sapply(simplify=F,dists,function(distance){
      #We want to compute the overlap with the correct chip peak type
      #just at 68 hours for now
      tp = names(tfgrlist)[[1]]    
      peaks=resize(tfgrlist[[tp]][[motif]],width=distance,fix='center')
      mots.gr=motifs.gr[motifs.gr$name==motif & motifs.gr$pvalue < pval]
      peaksoverlapping = mean(countOverlaps(peaks,mots.gr) > 0)
      genomic.density =  length(mots.gr)/genome.size
      total.motifs.inpeaks = sum(countOverlaps(mots.gr,peaks)>0)
      enr= (total.motifs.inpeaks/sum(width(peaks))) / genomic.density#for now, genome wide enrichment, later we can do non-repeat genome wide or something
      list(enrichment=enr,peakswithmotifs=peaksoverlapping)
    })
  })
})
peak.motifs.data.m=melt(peak.motifs.data)
head(peak.motifs.data.m)
colnames(peak.motifs.data.m)=c('value','measurement','distance','motif','pvalue')
peak.motifs.data.m$pvalue = as.numeric(peak.motifs.data.m$pvalue)
peak.motifs.data.m$lpvalue = -log10(as.numeric(peak.motifs.data.m$pvalue))

#now with the log -ve pvalues.
pdf('analysis/load_motifs/threshold.plots.pdf')
qplot(data=peak.motifs.data.m,x=lpvalue,y=value,color=numdist,group=numdist,geom='point')+
facet_grid(measurement~motif,scale='free')+
geom_line()
dev.off()




#Now plot the number of peaks with TFs as a function of the pvalue threshold
dir.create('analysis/load_motifs/')
pdf('analysis/load_motifs/threshold.plots.pdf')
dat=peak.motifs.data.m[peak.motifs.data.m$measurement=='peakswithmotifs',]
ggplot(data=dat,aes(x=as.numeric(lpvalue),y=value,group=distance,facet=motif))+
facet_wrap(~motif)+
geom_point()+
geom_line()
dev.off()

pdf('analysis/load_motifs/threshold.plots.pdf')
dat=peak.motifs.data.m[peak.motifs.data.m$measurement=='peakswithmotifs',]
qplot(data=dat,x=lpvalue,y=value,color=distance,group=distance,facet=motif,log='x',geom='point')+
facet_wrap(~motif)+
geom_line()
dev.off()

peak.motifs.data.m$numdist=as.numeric(substr(peak.motifs.data.m$distance,start=1,stop=3))

pdf('analysis/load_motifs/threshold.plots.pdf')
dat=peak.motifs.data.m[peak.motifs.data.m$measurement=='peakswithmotifs',]
qplot(data=dat,x=pvalue,y=value,color=distance,group=distance,facet=motif,log='x',geom='point')+
facet_wrap(~motif)+
geom_line()
dev.off()



pdf('analysis/load_motifs/threshold.plots.pdf')
qplot(data=peak.motifs.data.m,x=lpvalue,y=value,color=numdist,group=numdist,geom='point')+
facet_grid(measurement~motif,scale='free')+
geom_line()
dev.off()
#And enrichment relative to genomic
#Collect all this data into a dataframe

pdf('analysis/load_motifs/threshold.plots.pdf')
qplot(data=peak.motifs.data.m,x=numdist,y=value,color=pvalue,group=pvalue,log='x',geom='point')+
scale_color_gradient(name = "Pvalue", trans = "log")+
facet_grid(measurement~motif,scale='free')+
geom_line()
dev.off()
#Now plot in GGPLOT




#Now, to callibrate the number of motifs we call for each TF, we'll make sure we're getting about the same numbers as Ohler et al
#


#Export our motifs to a bed folder
for(motif in unique(motifs.gr$name)){
  tmp=motifs.gr[motifs.gr$name == motif]
  export(tmp,file='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/motifbeds/')
}



