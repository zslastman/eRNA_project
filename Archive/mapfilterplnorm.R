setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
rootfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/'
outfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/mapfilterplnorm'
dir.create(outfolder,showWarnings=F)

load('data/objects/cg.pl.object.R')
load('data/objects/allmap.object.R')
# #now make a version with only sites that are totally mappable
# cg.pl.map=rapply(cg.pl,how='replace',function(SRL){
# 	SRL[!allmap]<-0
# # 	SRL
# # })

cg.pl.map=cg.pl

# load('data/objects/cg.pl.object.R')
# load('data/objects/allmap.object.R')
# #now make a version with only sites tha

nallmap=!allmap
cg.pl.map=mclapply(mc.cores=5,names(cg.pl),function(acc){
	sapply(c('pos','neg'),function(s){
		cg.pl[[acc]][[s]][nallmap]<-0
		cat('.')
		cg.pl[[acc]][[s]]
	})
})
names(cg.pl.map) = names(cg.pl)

save(cg.pl.map,file='data/objects/cg.pl.map.object.R')