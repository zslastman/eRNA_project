setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
rootfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/'

#we'll perform steps 5 and onwards once for each normalization method
#b)sqrt(lib size)
#c)power law
#c.2)power law with the original equation
#d)qnorm((rank(reads) - 0.5)/length(reads), ties.method='average')
#e)qunatile normalization but counting sites with reads in ANY library


#a)library size
cg.libnorm=sapply(simplify=F,names(cg),function(acc){
	sapply(c('pos','neg'),function(strand){
		cg[[acc]][[strand]]/with(accession.df,genome.coverage[accession==acc])
	})
})
length(cg.libnorm)

#bsqrt(libsize)
cg.sqrt=sapply(simplify=F,names(cg),function(acc){
	sapply(c('pos','neg'),function(strand){
		sqrt(cg.libnorm[[acc]][[strand]])
	})
})
save(cg.libnorm,file='data/objects/cg.libnorm.object.R')
save(cg.sqrt,file='data/objects/cg.sqrt.object.R')
rm(cg.libnorm)
rm(cg.sqrt)


#c power law
cg.pl=fit.and.plnorm(cg)
save(cg.pl,file='data/objects/cg.pl.object.R')
rm(cg.pl)

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


#d)qnorm((rank(reads) - 0.5)/length(reads), ties.method='average')
cg.qn=mclapply(mc.cores=10,unique(unlist(accs)),function(acc){
	sites=sapply(cg[[acc]][1:2],'>',0)
	readvect = lapply(c('pos','neg'),function(st){
		alls = cg[[acc]][[st]]
		sites.s = sites[[st]]
		c(unlist(unname( alls[sites.s] )))
	})
	readvect = as.vector(do.call('c',readvect))
	readvect = qnormvect(readvect)
	#split by strand
	readvect=chunkvect(readvect,sapply(sites,function(s){sum(sum(s))}) )
	#now split by chr and return
	sapply(names(readvect),function(st){
		chrsites=sapply(sites[[st]],function(s){ sum(s) })
		sitesasrleslist=as(chunkvect(readvect[[st]],chrsites,as.rle=T),'SimpleRleList')
		cg[[acc]][[st]][sites[[st]]]<-sitesasrleslist#and enter normalized vals back into our sites
	})
})
names(cg.qn) <- unique(unlist(accs))
save(cg.qn,file='data/objects/cg.qn.object.R')
rm(cg.qn)


#e)qunatile normalization but counting sites with reads in ANY library
load('/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/alltaglist.object.R')
sites=sapply(c('pos','neg'),function(st){
	sites=lapply(alltaglist,'[[',st)
	sites=sapply(sites,'>',0)
	Reduce('|',sites)
})#make pos/neg list 
cg.qn.allsites=cg
cg.qn.allsites=mclapply(mc.cores=10,unique(unlist(accs)),function(acc){
	readvect = lapply(c('pos','neg'),function(st){
		alls = cg[[acc]][[st]]
		sites.s = sites[[st]]
		c(unlist(unname( alls[sites.s] )))
	})
	readvect = as.vector(do.call('c',readvect))
	readvect = qnormvect(readvect)
	#split by strand
	readvect=chunkvect(readvect,sapply(sites,function(s){sum(sum(s))}) )
	#now split by chr and return
	sapply(names(readvect),function(st){
		chrsites=sapply(sites[[st]],function(s){ sum(s) })
		sitesasrleslist=as(chunkvect(readvect[[st]],chrsites,as.rle=T),'SimpleRleList')
		cg[[acc]][[st]][sites[[st]]]<-sitesasrleslist#and enter normalized vals back into our sites
	})
})
names(cg.qn.allsites) <- unique(unlist(accs))
save(cg.qn.allsites,file='data/objects/cg.qn.allsites.object.R')
rm(cg.qn.allsites)

