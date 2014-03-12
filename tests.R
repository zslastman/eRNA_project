#This file contains tests of the various objects used in these scripts, put here for easy access
library(pryr)

expect_SRL<-function(x,s=si){
	# This function tests a SimpleRLeList object to make sure it
	# is correct
	 expect_is(x,'SimpleRleList')
	 expect_equal(names(x),seqnames(si))
	 expect_equal(sapply(x,length),seqlengths(si))
	 expect_false(all(sum(x)==0))
   NULL
}


expect_GR<-function(x,s=si){
  # This function tests a SimpleRLeList object to make sure it
  # is correct
   expect_is(x,'GRanges')
   expect_identical(seqinfo(x),si)
   expect_false(length(x)==0)
   NULL
}

#check we have a list of numeric matrices
expect_num_mat_list<-function(matlist){
  #recursively check structre, values not zero
  rapply(matlist,function(x){expect_true(identical(is(x),is(matrix(0))))})
  rapply(matlist,function(x){expect_false({all(x==0)})})
  NULL
}


load('data/objects/cg.object.R')
load('data/objects/cg.mapfilt.object.R')
load('data/objects/cg.pl.object.R')
load('data/objects/cg.mapfilt.pl.object.R')
load('data/objects/accession.df.object.R')

expect_is(cg[[1]]$neg,'SimpleRleList')

#test the cage data
test_that('Cage data read and processed correctly',{
  #The metadata
  expect_false(any(grepl(accession.df,pat='/'))) #no filenames unparsed
  expect_false(any(duplicated(accession.df$accession)))
  #check data rle
  expect_is(cg,'list')#is a list
  expect_true(all(sapply(cg, function(x){names(x)==c('pos','neg','both')})))#all three strands
  expect_equal(length(cg),nrow(accession.df))#correct length
  expect_identical(accession.df$accession,names(cg))#correct names
  tmp= rapply(cg,expect_SRL)#right type of objects
  tmp = rapply(cg,function(x){expect_true( runValue(x)[[1]][[2]]%%1 == 0 )}) #test for integerness
  expect_false(any(duplicated(sapply(cg,f(sum(x[[1]][[1]]))))))#not duplicated
  #mapfiltered object
  expect_identical(names(cg.mapfilt),names(cg))#names okay
  expect_true(all(sapply(cg.mapfilt, function(x){names(x)==c('pos','neg','both')})))#all strands
  tmp = rapply(cg.mapfilt,function(x){ expect_true( all(all(x[!allmap]==0))) ;NULL})#mapfilter worked
  #power law norm
  expect_identical(names(cg.pl),names(cg))#names okay
  expect_true(all(sapply(cg.pl, function(x){names(x)==c('pos','neg','both')})))#all strands
  tmp = rapply(cg.pl,function(x){ expect_false( all(all(x[!allmap]==0)));NULL })#mapfilter not done yet
  tmp = rapply(cg.pl[1:2],f(expect_true( any(runValue(x)[[1]][[2]]%%1 != 0)) )) #test for NOTintegerness
  #mapfiltered , power law norm object
  expect_identical(names(cg.mapfilt.pl),names(cg))#names okay
  expect_true(all(sapply(cg.mapfilt.pl, function(x){names(x)==c('pos','neg','both')})))#all strands
  tmp = rapply(cg.mapfilt.pl,function(x){ expect_true( all(all(x[!allmap]==0))) })#mapfilter worked
  tmp = rapply(cg,f(expect_true( any(runValue(x)[[1]][[2]]%%1 != 0) ))) #test for NOTintegerness
  #test the allcage object
  splitaccs=with(accession.df,split(accession,paste0(timepoint,tissue)))
  expect_equal(length(cg),length(splitacces))#correct length
  expect_identical(names(allcage),names(splitaccs))#names okay
  expect_true(all(sapply(allcage, function(x){names(x)==c('pos','neg','both')})))#all strands  
  tmp = rapply(allcage,function(x){expect_true( runValue(x)[[1]][[2]]%%1 == 0 )}) #test for integerness
  #nothing duplicated
  expect_false(any(duplicated(sapply(cg.pl,f(sum(x[[1]][[1]]))))))
  expect_false(any(duplicated(sapply(cg,f(sum(x[[1]][[1]]))))))
  expect_false(any(duplicated(sapply(cg,f(sum(x[[1]][[1]]))))))

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
  # expect_true( 'cg' %in% colnames(mcols(transcripts.gr)) )
  # expect_true( 'cg.pl' %in% colnames(mcols(transcripts.gr)) )
  # expect_true( 'ts' %in% colnames(mcols(transcripts.gr)) )
  # expect_true( 'ts.pl' %in% colnames(mcols(transcripts.gr)) )
  # expect_true( 'rna.seq' %in% colnames(mcols(transcripts.gr)) )
  # expect_true( 'c.rna.seq' %in% colnames(mcols(transcripts.gr)) )
  # expect_true( 'z.rna.seq' %in% colnames(mcols(transcripts.gr)) )
  # expect_true( 'cg' %in% colnames(mcols(transcripts.gr)) )

  # sapply(gene.annotations,function(g){expr.da})
  # identical(rapply(expr.mats,is,how='replace'),rapply(cg.pl,is,how='replace'))
  expect_num_mat_list(expr.mats)
  expect_num_mat_list(expr.crm.mats)



})
  #chose a random set from each gene region

  totest=sample(1:length(region),10)
  chrs= seqnames(region[totest])
  starts=start(region)[totest]
  ends=end(region)[totest]
  #now for each region get the coordinates and chr
  #use this to access our objects
  #check our matrix matches the results

  #get the score in these random sets according to our matrices

  expect_true(expr.mats)




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

makeInspectionCovPlot<-function(x,  actchr = 'chr2R', actstart = 5794894 - 10000 ,  actstop = 5794894 + 10000){
  require(testthat)
	require(chipseq)
	expect_SRL(x)
	pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/inspection_coverage_plot_tmp.pdf')
	coverageplot(Views(x[[actchr]],actstart,actstop))
	dev.off()
}

pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/inspection_coverage_plot_tmp.pdf')
actchr = 'chrX'
  actstart = 18101918 - 1000 
  actstop = 18102417 + 1000
x = allcage.pl[['24hembryo']]$both
coverageplot(Views(x[[actchr]],actstart,actstop))
dev.off()



viewSums(GRViews(gr=cad3.gr[569],rle=allcage.pl[[2]]$both))
mean(expr.crm.mats$cad3.gr$cg.pl[569,accession.df$time=='24h'])


coverageplot(Views(x[[actchr]],actstart,actstop))
makeInspectionCovPlot(allcage.pl[[2]]$both)

makeInspectionCovPlot(c.rna.seq[[1]])

makeGvizCovPlot<-function(x){
  require(testthat)
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

#Testing Snippet
test_that('Test peakcounts object',{
  peakcounts
  #test length of object
  #test names of object
  #type
  expect_is(peakcounts,'matrix')
  expect_true(is.numeric(peakcounts))
  expect_false(is.character(peakcounts))
  expect_false(is.factor(peakcounts))
  expect_false(is.list(peakcounts))
  #test values of object
  expect_true( all(peakcounts != 0))
  #test dimensions of object
  expect_equal( dim(peakcounts) , matrix  )
  expect_identical(colnames(  peakcounts),  accession.df$accession)
  expect_identical(rownames(  peakcounts),  mypeaks$id)
}
)


