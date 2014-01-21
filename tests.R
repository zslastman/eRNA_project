#This file contains tests of the various objects used in these scripts, put here for easy access
library(pryr)

expect_SRL<-function(x,s=si){
	# This function tests a SimpleRLeList object to make sure it
	# is correct
	 expect_is(x,'SimpleRleList')
	 expect_equal(names(x),seqnames(s))
	 expect_equal(sapply(x,length),seqlengths(s))
	 expect_false(all(sum(x)==0))
}


#test the cage data
test_that('Cage data read and processed correctly',{
  #The metadata
  expect_false(any(grepl(accession.df,pat='/'))) #no filenames unparsed
  #check data rle
  expect_is(cg,'list')
  expect_true(all(sapply(cg, function(x){names(x)==c('pos','neg','both')})))
  expect_equal(length(cg),nrow(accession.df))
  rapply(cg,f(expect_equal(x[[!allmap]],0)))
  rapply(cg,expect_SRL)
  rapply(cg,f(expect_true( runValue(x)[[1]][[2]]%%1 == 0 ))) #test for integerness
  #check map filtering
  expect_true(!all(sum(cg[[1]][[1]][!allmap])==0))#cg not filtered
  expect_false(!all(sum(cg.mapfilt[[1]][[1]][!allmap])==0))#cg not filtered
  #check power law
  #check they have same structure, names
  identical(rapply(cg,is,how='replace'),rapply(cg.pl,is,how='replace'))
  #alphas should be in the 1.5 to 3 range
  expect_true(all(accession.df$alpha > 1.4))
  expect_true(all(accession.df$alpha < 3))
  #
  rapply(cg.pl,f(expect_true( runValue(x)[[1]][[2]]%%1 != 0 ))) #test for integerness
  expect_true(!all(sum(cg.pl[[1]][[1]][!allmap])==0))#cg.pl not filtered
  expect_false(!all(sum(cg.mapfilt.pl[[1]][[1]][!allmap])==0))#cg.mapfitl filtered
  #check 
}
)

#test the tagseq data
test_that('Tagseq data read and processed correctly',{
  #The metadata
  expect_false(any(grepl(tagseq.df,pat='/'))) #no filenames unparsed
  #check data rle
  expect_is(ts,'list')
  expect_true(all(sapply(ts, function(x){names(x)==c('pos','neg','both')})))
  expect_equal(length(ts),nrow(tagseq.df))
  rapply(ts,expect_SRL)
  rapply(ts,f(expect_true( runValue(x)[[1]][[2]]%%1 == 0 ))) #test for integerness
  #check power law
  #check they have same structure, names
  identical(rapply(ts,is,how='replace'),rapply(ts.pl,is,how='replace'))
  #alphas should be in the 1.5 to 3 range
  expect_true(all(tagseq.df$alpha > 1.4))
  expect_true(all(tagseq.df$alpha < 3))
  #
  rapply(ts.pl,f(expect_true( runValue(x)[[1]][[2]]%%1 != 0 ))) #test for integerness
  #check 
}
)

#test the annotation data
test_that('Annotation data read and processed correctly',{
	#correct number of transcripts
	#few tss, all tss are starts of transcript
  act_id=synonyms$gene_id[synonyms$syn == 'FBgn0000042']
  expect_true(is.integer(act_id))
  expect_true()
  expect_true(length(transcripts.gr)>length(tss.gr))
  expect_true(length(tss.gr)>12000)
  expect_true(all(width(tss.gr))==1 )
  expect_true( all( start(tss.gr) %in% start(transcripts.gr) | 	start(tss.gr) %in% end(transcripts.gr) ) )
  #check twist is there in the correct place

}
)

#test the chromatin data
test_that('Cage data read and processed correctly',{

  expect_is(cg,'list')
  expect_is(cg.pl,'list')
  expect_equal(length(cg),length(cg.pl))
  expect_equal(length(cg),nrow(accession.df))
  rapply(cg,expect_SRL)
}
)


#test the RNA data
test_that('RNA data read and processed correctly',{
	#Our RNA data
  expect_equal(length(rna.seq),2)#length
  rapply(rna.seq,expect_SRL)#correct SRLs
  expect_true( all(names(rna.seq) %in% c('meso','embryo')) )#named as timepoints
  #pos/neg structure
  length(grep("pos", names(rapply(rna.seq,names)))) == length(grep("neg", names(rapply(rna.seq,names))))
  #Celniker data
  expect_equal(length(c.rna.seq),3)#length
  a=sapply(c.rna.seq,expect_SRL)#correct SRLs
  expect_true(all(grepl(x=names(c.rna.seq),pat='tp\\d+h')))#named as timepoints
}
)
save(rna.seq,c.rna.seq.wig,c.rna.seq,z.rna.rles,file='data/objects/rna.seq.object.R')

test_that('test',{expect_SRL(cg[[1]][[1]]}))


#Test the annotations have been scored correctly
test_that('Datasets and CRM/Gene Annotation correctly integrated',{
  #
  expect_true( 'cg' %in% colnames(mcols(transcripts.gr)) )
  expect_true( 'cg.pl' %in% colnames(mcols(transcripts.gr)) )
  expect_true( 'ts' %in% colnames(mcols(transcripts.gr)) )
  expect_true( 'ts.pl' %in% colnames(mcols(transcripts.gr)) )
  expect_true( 'rna.seq' %in% colnames(mcols(transcripts.gr)) )
  expect_true( 'c.rna.seq' %in% colnames(mcols(transcripts.gr)) )
  expect_true( 'z.rna.seq' %in% colnames(mcols(transcripts.gr)) )
  expect_true( 'cg' %in% colnames(mcols(transcripts.gr)) )



})



#This is for visually inspecting SRLs
#genomic location of Actin C - +/- 20kb
SRL = cg[[1]][[1]]


actgr=genes.gr[ which(genes.gr$symbol =='Act5C') ]
actsum = unlist(viewSums(GRViews(SRL,actgr)))
expect_true(actsum > 0)

actgr=genes.gr[ which(genes.gr$symbol =='twi') ]
actsum = unlist(viewSums(GRViews(SRL,actgr)))
expect_true(actsum > 0)

export(actgr,'/g/furlong/Harnett/TSS_CAGE_myfolder/tmp.bed')
export(SRL,'/g/furlong/Harnett/TSS_CAGE_myfolder/tmp.bw')

makeInspectionCovPlot<-function(x){
	require(chipseq)
	expect_SRL(x)
	actchr = 'chr2R'
	actstart = 5794894 - 10000 
	actstop = 5794894 + 10000
	pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/inspection_coverage_plot_tmp.pdf')
	coverageplot(Views(x[[actchr]],actstart,actstop))
	dev.off()


}
makeInspectionCovPlot(c.rna.seq[[1]])

makeGvizCovPlot<-function(x){
  require(chipseq)
  expect_SRL(x)
  actchr = 'chr2R'
  actstart = 5794894 - 10000 
  actstop = 5794894 + 10000
  expect_
  pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/inspection_coverage_plot_tmp.pdf')
  coverageplot(Views(x[[actchr]],actstart,actstop))
  dev.off()
  
  
}

# library(chipseq)
# pdf('coverageplots.celniker.pdf')
# coverageplot(Views(c.rna.seq[[8]][['chr2R']],18931631,18931631),Views(c.rna.seq[[7]][['chr2R']],18931631,18937849),main='wig')
# coverageplot(Views(tmp1[['chr2R']],18931631,18937849))
# coverageplot(Views(tmp2[['chr2R']],18931631,18937849))
# coverageplot(Views(tmp3[[1]][['chr2R']],18931631,18937849))
# coverageplot(Views(tmp4[[1]][['chr2R']],18931631,18937849))
# dev.off()







