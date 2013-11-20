
#classify our crms using the tf data
#get the mean normalized rpgc figures for each crm in each class
tmp<-list(Heart_5TF=crmgrs$allsum.rpgc[ crmgrs$Heart5 & crmgrs$intergenic],
          Meso_5TF=crmgrs$allsum.rpgc[ crmgrs$Meso5 & crmgrs$intergenic],
          Heart_2TF=crmgrs$allsum.rpgc[ crmgrs$Meso2 & crmgrs$intergenic],
          Meso_2TF=crmgrs$allsum.rpgc[ crmgrs$Heart2 & crmgrs$intergenic])
tmp<-stack(tmp)
#do boxplot
jpeg('analysis/make_regions_bedfiles/Heart_Meso_5bound_2bound_TF_boxplot.jpeg')
qplot(y=tmp$values,log='y',x=tmp$ind,color=tmp$ind,geom='boxplot',main='RPGC Cage signal for 4 classes of CRMs',ylab='normalized CAGE signal',xlab='')
dev.off()
#show sizes of sets
table(tmp$ind)
