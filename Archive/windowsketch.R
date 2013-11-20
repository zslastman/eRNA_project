

# and the log10 plots -----------------------------------------------------


#for each techrepped line
for(l in 1:length(multiaccs[[1]])){
  #for each of 5 window sizes
  w.df<-m[ ,]
  accs=sapply(multiaccs,'[[',l)
  
  jpeg(width=2000,height=2000,filename=paste0(accs[1],' tmp.log10.jpeg'))
  par(mfrow=c(length(windows.toplot),3))
  for(w in windows.toplot){
    #for each region type
    # max(w.df<-m[m$windowsize== w&m$region==reg &m$acc %in% accs& m$windowsize%in%windows.toplot & m$variable=='score',]$value)->maxval
    #now have max val for plotting limits
    for(reg in c('pos','neg','int')){
      
      w.df<-        m[m$windowsize== w & m$region==reg &m$acc %in% accs& m$windowsize%in%names(windows.toplot) & m$variable=='score',]
      #now normalize rpgc
      #plots showin correspondance of windows in technical replicates
      maxval.n<-quantile(m[m$region =='pos' & m$windowsize== w &m$acc %in% accs& m$windowsize%in%windows.toplot & m$variable=='score',]$val.libnorm,0.99)
      
      x=w.df[w.df$acc==accs[1],]$val.libnorm
      y=w.df[w.df$acc==accs[2],]$val.libnorm
      xmax<-100/v[accs[1]]
      ymax<-100/v[accs[2]]
      # x=w.df[w.df$acc==accs[1],]$value
      #y=w.df[w.df$acc==accs[2],]$value
      
      title=paste0('Correspondance between technical replicates for line ',strsplit(multiaccs[[1]][l],'_')[[1]][2],
                   '\nRegion ',reg,' ,windowsize ',w)
      plot(x,y,main=title,xlim=c(0,maxval.n),ylim=c(0,maxval.n))
      
      plot(x,y,main=title,xlim=c(0,maxval.n),ylim=c(0,maxval.n))
      plot(log10(x),log10(y),xlim=c(0,log10(maxval.n)),ylim=c(0,log10(maxval.n)),main=title)
      # print( ggplot(w.df,aes(x,y))+geom_point()+scale_y_log10()+scale_x_log10()+ggtitle(title)  )
      
    }
  }
  dev.off()
}



# Now fdrs for individual lines -------------------------------------------
#
#so now we have our window size

#assing a p value to each window based on a poisson distribution for that windowsize and line
m$pval<-NA
lamdas<-list()
for(w in c('w500')){
  lamdas[[w]]<-list()
  mclapply(mc.cores=10,unique(m$acc),function(acc)){
    
    w.df<-m[ m$variable=='score' & m$windowsize==w & m$region=='int'& m$acc == acc,]
    
    trimint<-w.df$value[w.df$value<quantile(w.df$value,0.99)]
    lam<-mean(trimint)
    lam
    cat('.')
    #  m[ m$variable=='score' & m$windowsize==w & m$acc == accession[1],]$pval<-
    #   dpois(m[ m$variable=='score' & m$windowsize==w & m$acc == accession[1],]$value,lam)
  })
}
save.image('tmp.image.R')
w.df$called<- w.df$value>cutoffs[w.df$acc]



#window.sum.distributions.rpgc<-window.sum.distributions.rpgc.int<-window.sum.distributions
# 
# 
# #normalize by library size
# for(reg in names(window.sum.distributions)){
#  for(w in names(window.sum.distributions[[reg]])){
#    for(n in names(window.sum.distributions[[reg]][[w]])){
#     i=which(accession.df$accession==n)
#    window.sum.distributions.rpgc<-window.sum.distributions[[reg]][[w]][[n]]/accession.df$library.size[i]
#  }
# }
# }
# #and by intergenic library size
# #normalize by library size
# for(reg in names(window.sum.distributions)){
#   for(w in names(window.sum.distributions[[reg]])){
#     for(n in names(window.sum.distributions[[reg]][[w]])){
#       i=which(accession.df$accession==n)
#       window.sum.distributions.rpgc.int<-window.sum.distributions[[reg]][[w]][[n]]/accession.df$intergenic.lib.size[i]
#     }
#   }
# }

# 
# #how many windows are zero?
# df<-rapply(window.sum.distributions,how='replace',function(x){mean(x==0)})
# #boxplot showing proportion of zero values for all 
# ggplot(melt(df),aes(x=L3,y=value,fill=factor(L2),group=factor(L2)))+stat_summary(fun.y='median',geom='bar',position='dodge')+
#   facet_wrap(~L1)





# now do density plots ----------------------------------------------------
df<-melt(rapply(window.sum.distributions,how='replace',function(x){v<-(as.vector(x));v[v<quantile(v,0.999)]}))
h(df)
df<-df[df$L2==w & df$L3 %in% accession.df$accession[1:5],]

#now produce density plots for the lines, showing distribution of counts for the inactive, negative and postiive regions
ggplot(df,aes(x=value,color=L3,facet=L1))+geom_density()+facet_grid(.~L1)+scale_x_log10()#+scale_y_log10()

df<-melt(rapply(window.sum.distribution ,how='replace',function(x){v<-(as.vector(x));v[v<quantile(v,0.999)]}))
df<-df[df$L2==w & df$L3 %in% accession.df$accession[1:10],]
#now produce density plots for the lines, showing distribution of counts for the inactive, negative and postiive regions
ggplot(df,aes(x=value,color=L3,facet=L1))+geom_density()+facet_grid(.~L1)+scale_x_log10()#+scale_y_log10()


#we want to see if dividing by library size normalizes the distributions. (I suspect not)


# Now do discriminance analysis -------------------------------------------
#supposing we have a our cutoff set for each 




for(l in 1:length(multiaccs[[1]])){
  #for each of 5 window size
  jpeg(width=2000,height=2000,filename=paste0(accs[1],' tmp.jpeg'))
  par(mfrow=c(length(windows.toplot),3))
  for(w in windows.toplot){
    #for each region type
    w.df<-m[m$windowsize== w &m$acc %in% accs& m$windowsize%in%names(windows.toplot) & m$variable=='score',]
    #now have max val for plotting limits
    maxval.n<-quantile(w.df$val.libnorm,0.99)
    
    for(reg in c('pos','neg','int')){
      
      #now normalize rpgc
      #plots showin correspondance of windows in technical replicates
      x=w.df[w.df$region==reg&w.df$acc==accs[1],]$val.libnorm
      y=w.df[w.df$region==reg&w.df$acc==accs[2],]$val.libnorm
      # x=w.df[w.df$acc==accs[1],]$value
      #y=w.df[w.df$acc==accs[2],]$value
      maxval.n<-500
      
      
      plot(x,y,main=title,xlim=c(0,maxval.n),ylim=c(0,maxval.n))
      
      print(qplot(x,y,main=title,xlim=c(0,maxval.n),ylim=c(0,maxval.n),geom='jitter'))
      #plot(log10(x),log10(y),xlim=c(0,log10(maxval.n)),ylim=c(0,log10(maxval.n)),main=title)
      # print( ggplot(w.df,aes(x,y))+geom_point()+scale_y_log10()+scale_x_log10()+ggtitle(title)  )
      
    }
  }
  dev.off()
}



#plots make no sense, inspecting data
tmp.bak<-tmp
tmp<-tmp[ tmp$windowsize==500,]
tmp<-tmp[ tmp$region=='pos',]


#my code does seem to be returning the correct values. e.g.
h(tmp,n=30)
regs[[1]][30]
cage.tag.rles[[accs[1]]]$pos$chr2L[8846672:8847171]
#or
h(tmp,n=399)
sum(cage.tag.rles[[accs[1]]]$neg[[as.character(seqnames(regs[[1]][399]))]][start(regs[[1]][399]):end(regs[[1]][399])])

tmp[tmp$region=='pos',][ which(tmp[ tmp$region=='pos',]$val.libnorm==max(tmp[tmp$region=='pos',]$val.libnorm)), ]
tmp[tmp$region=='int',][ which(tmp[ tmp$region=='int',]$val.libnorm==max(tmp[tmp$region=='int',]$val.libnorm)), ]
#just 8 values in tmp now
tmp<-tmp[tmp$val.libnorm>500,]
ggplot(data.frame(),aes(x=tmp$val.libnorm,y=tmp$val2))+geom_jitter(h=0.1)
qplot(tmp$val.libnorm,tmp$val2,geom='jitter')+geom_abline(intercept=0)
dev.off()
# #inspect m
# z<-m[ m$variable=='score',]
# zp<-m[ m$variable=='pos',]
# z<-z[ m$region=='int'&z$acc %in% accs& z$windowsize%in%c('w100','w500'),]
# zp<-zp[ m$region=='int'&z$acc %in% accs& z$windowsize%in%c('w100','w500'),]
# 
# z100<-z[ z$windowsize=='w100' ,]
# z500<-z[ z$windowsize=='w500' ,]       
# max( z$value )