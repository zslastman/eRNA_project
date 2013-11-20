# checking out negatives with high cage tags ------------------------------

#just verifying my data is okay in IGV

unlist(viewSums(Views(alltags$neg, as(cad3.neg,'RangesList'))))

cad3.neg<-sort(cad3.neg)
cad3.pos<-sort(cad3.pos)
unlist(viewSums(Views(alltags$pos[bigchrs],as(cad3.neg,'RangesList')[bigchrs])))
unlist(viewSums(Views(alltags$neg[bigchrs],as(cad3.neg,'RangesList')[bigchrs])))

unlist(viewSums(Views(cage.tag.rles[[1]]$pos[bigchrs],as(cad3.neg,'RangesList')[bigchrs])))
unlist(viewSums(Views(cage.tag.rles[[1]]$neg[bigchrs],as(cad3.neg,'RangesList')[bigchrs])))

findOverlaps(cad3.neg,chrompeaks[['PolII_6-8h']])

r<-cad3.neg[8]
sum(alltags$pos[[as.character(seqnames(r))]][start(r):end(r)])
sum(alltags$neg[[as.character(seqnames(r))]][start(r):end(r)])
plot(alltags$pos[[as.character(seqnames(r))]][start(r):end(r)])
plot(alltags$neg[[as.character(seqnames(r))]][start(r):end(r)])


#now (using stuff from the windowplots script) display names of positives so we can view
pos<-sums[[1]]
names(pos)<-regs[[1]]$name
sort(pos)

neg<-sums[[2]]
names(neg)<-regs[[2]]$name
sort(neg)

ne#quick thing to produce a 