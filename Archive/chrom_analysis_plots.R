
#classify our crms using the tf data
#get the mean normalized rpgc figures for each crm in each class
tmp<-list(Heart_5TF=crm8008.gr$allsum.rpgc[ crm8008.gr$Heart5 & crm8008.gr$intergenic],
          Meso_5TF=crm8008.gr$allsum.rpgc[ crm8008.gr$Meso5 & crm8008.gr$intergenic],
          Heart_2TF=crm8008.gr$allsum.rpgc[ crm8008.gr$Meso2 & crm8008.gr$intergenic],
          Meso_2TF=crm8008.gr$allsum.rpgc[ crm8008.gr$Heart2 & crm8008.gr$intergenic])
tmp<-stack(tmp)
#do boxplot
jpeg('analysis/make_regions_bedfiles/Heart_Meso_5bound_2bound_TF_boxplot.jpeg')
qplot(y=tmp$values,log='y',x=tmp$ind,color=tmp$ind,geom='boxplot',main='RPGC Cage signal for 4 classes of CRMs',ylab='normalized CAGE signal',xlab='')
dev.off()
#show sizes of sets
table(tmp$ind)
