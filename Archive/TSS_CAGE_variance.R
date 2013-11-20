#Variance analysis on cage data
source('src/generate_chromdata_cage.R')
source('src/load_crm_data.R')
#library(Matrix)
#library(zoo)
library(reshape2)
library(ggplot2)
load(file.cage.tag.rles)
load(file.alltags)
#load(file.cage.tag.rles.rpgc)
load(file.accession.df)


replim=0.95

# Picking a simple dataset for each line ----------------------------------




#we need to select each line only once
#order each accession by library size
accession.df<-accession.df[order(accession.df$library.size,decreasing=T),]
#now pick a unique dataset for each line, prioritizing the largest


cage.tag.rles<-cage.tag.rles[eachline]
  




for(w in c(1,5)){
  message(paste0('computing windowed rles for windowsize ', w ))
  #get our rle lists, take the rolling sum of each 
  
  
  #go from the tag score to a logical indicating presence or absence, also flatten the list of rles
  cage.tag.rles.windowed.log<-rapply(cage.tag.rles,how='list',f=function(rl)runsum(rl,w)!=0)
  #cage.tag.rles.windowed.log<-rapply(cage.tag.rles,how='list',f=function(rl)rl!=0)
  
  #now we can just do a logical OR on our Rle lists to give us an Rle saying how many of our
  #lines have a score of any magnitude within a distance w each each site, on either strand
  message(paste0('computing windowed rles for windowsize ', w ))
  #thisi s a little tricky- do.call doesn't work properly on a list of SimpleRleLists
  all.windowed.log<-cage.tag.rles.windowed.log[[1]]
  #so we initilialize from the firs 
  all.windowed.log[['both']]<-cage.tag.rles.windowed.log[[1]][['pos']]|cage.tag.rles.windowed.log[[1]][['neg']]
  #then loop over the rest
  message('calculating reproducibility track..')
  lapply(cage.tag.rles.windowed.log[-1],function(x){
    #incrementing by 1 whenever a line has tags present
    all.windowed.log[['pos']]<<-x[['pos']]+all.windowed.log[['pos']]
    all.windowed.log[['neg']]<<-x[['neg']]+all.windowed.log[['neg']]
    #not that we don't want to count lines twice if they have plus and negative strand activity
    all.windowed.log[['both']]<<-(x[['pos']] | x[['neg']])+all.windowed.log[['both']]
    NULL
  })
  
    for(tn in 1:10){
    #now iterate over each line, and whenever we have a single tag, ask what the number of lines
    #is with another tag nearby - we don't need to do all the tags, 10k  from each line will do
      
    message(paste0('sampling locations with tag number', tn) ) 
    
    onelinecount<-mclapply(mc.cores=20,(cage.tag.rles),function(rl){
        sites<-c(unlist(use.names=F,all.windowed.log[['both']][rl$pos==tn]),unlist(use.names=F,all.windowed.log[['both']][rl$neg==tn])) 
        sample(sites,size=min(10000,length(sites)))
        })   
    names(onelinecount)<-names(cage.tag.rles)
    onelinecount.df<-melt(sapply(onelinecount,as.vector))
    if(ncol(head(onelinecount.df))==2){linecol=c(F,T)
    }else{ linecol<-sapply(1:ncol(head(onelinecount.df)),function(n)is.factor(head(onelinecount.df[,n])))}
    #melt seems to be giving erratic behavior in sorting the columns, have to do this manually
    names(onelinecount.df)[ linecol]<-'line'
    names(onelinecount.df)[ names(onelinecount.df)=='value']<-'count'
    onelinecount.df<-onelinecount.df[,c('count','line')]
    onelinecount.df$count<-as.numeric(onelinecount.df$count)/length(eachline)
    onelinecount.df$tn<-tn
    onelinecount.df$w<-w
    
    message('writing table')
    write.table(quote=F,sep=' ',row.names=F,col.names=T,
                onelinecount.df,
                file=paste0('analysis/tagvartables/tagvartable.w',w,'_tn',tn,'.txt'))
    }         
  
  
} 


tagvar.df<-do.call(rbind,lapply(list.files('analysis/tagvartables/',full.names=T),function(f){
  x<-read.table(header=T,f)
}))
head(tagvar.df)

quit()


tmp=sampledfrows(tagvar.df,num=10000)

colnames(tagvar.df)[1]<-'reproducibility'

#plot showing reproducibility distribution over peak heights
ggplot(tagvar.df[tagvar.df$w==1&tagvar.df$tn==1,]
       ,aes(x=reproducibility,color=factor(tn)))+geom_density()
ggplot(tagvar.df[tagvar.df$w==5&tagvar.df$tn==1,]
       ,aes(x=reproducibility,color=factor(tn)))+geom_density()

ggplot(tagvar.df[tagvar.df$w==1&tagvar.df$tn==3,]
       ,aes(x=reproducibility,color=factor(tn)))+geom_density()
ggplot(tagvar.df[tagvar.df$w==5&tagvar.df$tn==3,]
       ,aes(x=reproducibility,color=factor(tn)))+geom_density()


ggplot(tagvar.df[tagvar.df$w==1&tagvar.df$tn==5,]
       ,aes(x=reproducibility,color=factor(tn)))+geom_density()
ggplot(tagvar.df[tagvar.df$w==5&tagvar.df$tn==5,]
       ,aes(x=reproducibility,color=factor(tn)))+geom_density()

ggplot(tagvar.df[tagvar.df$w==1&tagvar.df$tn==10,]
       ,aes(x=reproducibility,color=factor(tn)))+geom_density()
ggplot(tagvar.df[tagvar.df$w==5&tagvar.df$tn==10,]
       ,aes(x=reproducibility,color=factor(tn)))+geom_density()









#plot showing reproducibility distribution over lines
lines<-accession.df$accession[ accession.df$run%in%c(1,2) ]
tagvar.df.sample<-tagvar.df[as.character(tagvar.df$line) %in% as.character(lines) ,]

unique(tagvar.df.sample$line)

ggplot(tagvar.df.sample,aes(x=count,color=factor(line)))+
  geom_density()+ 
  guides(colour = guide_legend(title.hjust = 0.20,ncol=2))+scale_x_log10()


lines<-accession.df$accession[ accession.df$run %in% c(1,2)  ]

tagvar.df.sample<-tagvar.df[tagvar.df$line %in% (lines),]

ggplot(tagvar.df.sample,aes(x=count,color=line))+geom_density()+scale_x_log10()+
  ggtitle('Distribution of reproducibility for heigth 1 tagsites')



#now for different tag numbers...
tagvar.df<-do.call(rbind,lapply(list.files('analysis/tagvartables/',pattern='w1',full.names=T),function(f){
  x<-read.table(header=T,f)
  x<-x[sample(1:nrow(x),size=100000,replace=F),]
  x$tn<-gsub(x=f,pattern='.*tn(\\d+)\\..*',replacement='\\1')
  x
}))







#already sorted the df by library size so...
smalllib=rev(accession.df$accession)[1]
any(tagvar.df$line==biglib)

tmp=tagvar.df[(tagvar.df$line %in% c(biglib,smalllib)),]

ggplot(tagvar.df[tagvar.df$line %in% c(biglib,smalllib),]
       ,aes(x=count,color=line))+geom_density()+scale_x_log10()







# 
# 
# 
# #actually if we're going ot sample we can do this much easier with a simple GRanges/views approach
# onelinecount<-lapply((cage.tag.rles),function(rl){
#   taglocs<-(c(as(rl$neg==1,'GRanges'),as(rl$pos==1,'GRanges')))
#   taglocs<-sample(taglocs[taglocs$score==T],10000)
#   taglocs<-resize(taglocs,width=w,fix='center')#resize them to our window size
#   
#   sample(which(rl$pos==1)$chr2R,10000)
#   
#   #now just 
#   cat('.')
#   
#   mclapply(mc.cores=20,cage.tag.rles,function(acc){
#     v1<-Views(acc$pos,as(taglocs,'RangesList'))
#     v2<-Views(acc$neg,as(taglocs,'RangesList'))
#     unlist(viewSums(v1))!=0 | unlist(viewSums(v2))!=0
#   })
# })
# 
#   
# names(onelinecount)<-names(cage.tag.rles)
# onelinecount<-melt(sapply(onelinecount,as.vector))
# onelinecount<-onelinecount[,-1]
# colnames(onelinecount)<-c('line','count')
# onelinecount$count<-onelinecount$count/length(eachline)
# 
# message('writing table')
# write.table(quote=F,sep=' ',row.names=F,col.names=T,
#             onelinecount,
#             file=paste0('analysis/tagvartables/tagvartable.w',w,'_tn',tn,'.txt'))
# }       
# 



# 
# 
# 
# #could i just feed all of the data into a huge, sparse matrix?
# #in theory this would be fine - 
# 
# alltagsboths<-alltags$pos+alltags$neg
# taggedsites<-alltagsboths!=0
# 
# 

# 
# #get an object containing the counts at each site for each line
# #using a logical Simple Rle List is *much* quicker
# countmats<-list()
# 
# 
# onesites<-mclapply(mc.cores=20,unlist(cage.tag.rles),function(rl){rl==1})
# onesites<-do.call(sum,onesites)
# onesites[[1]] + onesites[[2]]
# 
# #create a list of Rles - data for each of our lines, with chromosomes appended
 tmp<-lapply(eachline,function(x){
   Matrix(unlist(cage.tag.rles[[x]][['pos']],use.names=F),ncol=1,sparse=T)})
 #now create a sparse matrix for all of that, with our lines as columns
 countmats[['pos']]<-do.call(cBind,tmp)

# 
# 
# #create a list of Rles - data for each of our lines, with chromosomes appended
# tmp<-lapply(eachline,function(x){
#   Matrix(unlist(cage.tag.rles[[x]][['neg']],use.names=F),ncol=1,sparse=T)})
# #now create a sparse matrix for all of that, with our lines as columns
# countmats[['neg']]<-do.call(cBind,tmp)
# 
# rm(tmp)
# 
# countmats[['both']]<-countmats[['pos']]+countmats[['neg']]
# 
# save(countmats,file='data/objects/countmat.object.R')
# 
# 

# Compute variance at single tag locations --------------------------------


#Now I can easily roll along the matrix and get the windowed variances
#rows with 1 in them
onerows<-Rle(rowSums(countmats$pos==1)+rowSums(countmats$neg==1))

#create a new version of countmat by doing rolling sums over the columns (a little involved)
countmats.w<-countmats
countmats.w$pos<-do.call(cBind,mclapply(mc.cores=20,1:ncol(countmats$pos),function(i){
  Matrix(runsum(Rle(countmats$pos[,i]),5),sparse=T,ncol=1)}))
countmats.w$neg<-do.call(cBind,mclapply(mc.cores=20,1:ncol(countmats$neg),function(i){
  Matrix(runsum(Rle(countmats$neg[,i]),5),sparse=T,ncol=1)}))

#now we can simply ask - for our one rows, what is the distribution of the other lines 
countmats.w

# 
# 
# l# Transforming Rle to GR object -------------------------------------------
# 
# #rather than working with Rles, we can also just store the range that 'do' have tags in 'any' of our libraries
# alltags.gr<-as(alltags$pos+alltags$neg,'GRanges')
# alltags.gr<-alltags.gr[alltags.gr$score!=0]
# #now we need to turn all our GRange objects that have a width greater than 1 into width 1 ranges
# alltagrs.gr<-width1gr(alltags.gr)
# 
# #now get individual tag counts at each location
# #tagscoremat




# getting individual counts at tag sites ----------------------------------
# 
# alltags.gr.bak<-alltags.gr
# alltags.gr<-alltags.gr.bak
# alltags.resize.gr<-resize(alltags.gr,width=5,fix='center')
# 

#   #Okay so I now have a sparse matrix I can query
#   findOverlaps(alltags.gr.bak,ignoreSelf=T,maxgap=5)
#   
#   system.time({
#   
#   })
#   
#   
#   
#   
#   
#   
#   
# 
# # linecounting ------------------------------------------------------------
# 
#   
#   
#   #turn all values in our Cage Rles into ones
#   cage.tag.ones<-sapply(cage.tag.rles,simplify=F,function(acc){
#     sapply(acc,simplify=F,function(strand){
#       as(sapply(strand,function(rle){
#         runValue(rle)[runValue(rle)!=0]<-1
#         rle
#       }),'SimpleRleList')
#     })
#   })
#   
#   #now add up all our ones to give us a count of the ones with tags at each position
#   linecount<-cage.tag.ones[eachline][[1]]
#   lapply(cage.tag.rles[eachline][-1],function(x){
#     linecount[['pos']]<<-(x[['pos']])+linecount[['pos']]
#     linecount[['neg']]<<-(x[['neg']])+linecount[['neg']]
#   })
#   
#   
# 
# # Transforming linecount Rle to GR object (may not be useful)------------------------------------------
# 
# #rather than working with Rles, we can also just store the range that 'do' have tags in 'any' of our libraries
# linecount.gr.pos<-as(linecount$pos,'GRanges')
# linecount.gr.neg<-as(linecount$neg,'GRanges')
# linecount.gr<-as(linecount$pos+linecount$neg,'GRanges')
# 
# linecount.gr.pos<-linecount.gr.pos[linecount.gr.pos$score!=0]
# linecount.gr.neg<-linecount.gr.neg[linecount.gr.neg$score!=0]
# linecount.gr<-linecount.gr[linecount.gr$score!=0]
# 
# #now we need to turn all our GRange objects that have a width greater than 1 into width 1 ranges
# #for now we'll skip this
# linecount.gr.pos<- width1gr(linecount.gr.pos)
# linecount.gr.neg<- width1gr(linecount.gr.neg)
# linecount.gr<- width1gr(linecount.gr)
# 
# #now we can find how many are present in only 1 line
# oneline<-linecount.gr$score==1
# 
# #now we find the ones of those further than 300bp away from any other tag, in any other line
# d<-distanceToNearest(linecount.gr[oneline])
# badtags<-linecount.gr[which(oneline)[d$distance>50]]
# #turns out to be a farily small fraction - just 15k locations for 300bp, and 212k locations for 50bp
# 
# #how dense are our tags actually?
# d<-distanceToNearest(linecount.gr)
# ggplot(as.data.frame(d),aes(x=queryHits,y=distance))+geom_density()
# z=log10(sample(d$distance,size=1000))
# ggplot(data.frame(z=log10(sample(d$distance,size=1000)),w=0),aes(y=))+geom_density()
# 
# 
# 
# 
# #we could maybe just make a spare matrix out of this whole thing..
# #go through each line, get each tag that's equal to one, see what the percentage 
# #overlap is for each window size
# 
# tagn=1
# linenum<-length(eachline)
# windowsize=1
# acc<-eachline[1]
# strands<-c('pos','neg')
# 
# #We want plots showing the distriubtion of the linecount for tag locations with magnitude 1, 2,3 etc.
# #we want to go through each line
# cage.tag.grs<-mclapply(mc.cores=10,eachline,function(line){
#   sapply(strands,function(strand){
#     cat('.')
#     g<-as(cage.tag.rles[[line]][[strand]],'GRanges')#as a grange object
#     g<-g[g$score==tagn]#get the single tag locations
#     g<-width1gr(g)
#     v<-Views(linecount[[strand]],as(g,'RangesList'))#get views on the linecount
#     unlist(viewSums(v))#out as percentage.
#   })
# })
# 
# 
# tagnoverlap<-mclapply(mc.cores=10,eachline,function(line){
#   sapply(strands,function(strand){
#     cat('.')
#     g<-as(cage.tag.rles[[line]][[strand]],'GRanges')#as a grange object
#     g<-g[g$score==tagn]#get the single tag locations
#     g<-width1gr(g)
#     v<-Views(linecount[[strand]],as(g,'RangesList'))#get views on the linecount
#     unlist(viewSums(v))#out as percentage.
#   })
# })
# names(tagnoverlap)<-eachline
# library(ggplot2)
# library(reshape2)
# tagoverlap.df <-  melt(tagnoverlap)
# colnames(tagoverlap.df)<-c('overlap','strand','line')
# #now we melt that list down into a dataframe
# #...and do a density plot on it.
# ggplot(tagoverlap.df,aes(y=overlap))+geom_density()
# 
# #overlap them with the linecount GRange object
# 
# 
# #then we do views on our linecount object with a given radius and ask how many of them have something,
# #in any line, at a given radius.
# 
# 
# 
# 
# 
# 
# 
# #this should be easy to view with Gviz, and to do statistics on
# #First we can check the number of peaks which are not in any other dataset for one.
# #Then we can do this, again, for each 
# 
# 
# #We want a function that can display a regions's cage tages as a cloud of points over the right bp.
# #The easiest way of doing this would just be a matrix
# #probably we'll get views first.
# 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # #first get some of the peaks
# #   #use slice on the alltagset
# # #we can do this by just picking the 1000 highest value on chr2R, for now
# # cutoff<-sort(unique(as.vector(alltags[[strand]][[chr]])),decreasing=T)[20]
# # toppeakviews<-slice(alltags$pos,lower=cutoff)
# # 
# # #toppeakviews<-toppeakviews
# # #as a granges object
# # t<-as(toppeakviews,'GRanges')
# # seqinfo(t)<-si
# # toppeaks.gr<-GRanges(chr,as(toppeaks,'IRanges'),seqinfo=si)
# # toppeaks.gr<-resize(toppeaks.gr,width=500,fix='center')
# # toppeakviews<-toppeakviews[[chr]]
# # toppeaks.gr<-toppeaks.gr[1:length(toppeakviews)]
# # 
# # #peaks with more than one part above the cutoff
# # broadpeaks<-sapply(toppeakviews,function(x)x[width(x)>1])
# #                
# # #get the views for our peaks in all the accesions
# # all.acc.views<-sapply(accession.df$accession,function(acc){
# #   Views(as(cage.tag.rles.rpgc[[acc]][[strand]],'SimpleRleList'),as(toppeaks.gr,'RangesList') )[[chr]]})
# # 
# # #maximums in each line
# # maxmatrix<-sapply(all.acc.views,function(x){viewMaxs(x)})
# # #positions of max in each line
# # alltagpos<- viewWhichMaxs(toppeakviews)
# # #get matrix with columns as positions of these peaks, rows in each acc
# # posmatrix<-sapply(all.acc.views,function(acc) viewWhichMaxs(acc))
# # #turn each of these rows into a deviation from the max in alltags
# # devmat<-sapply(1:nrow(posmatrix),function(x)posmatrix[x,]-alltagpos[x])
# # #we seem to get large number of deviations even on the very high peaks
# # apply(devmat,2,table)[2]
# # 
# #                
# # #check if deviations tend to be on the same 
# # 
# #                
# #                
# #                
# #                
# # 
# # 
# # coverageplot(toppeaks
# #   Views(as(cage.tag.rles[[acc]][[strand]],'SimpleRleList'),as(toppeaks.gr[1],'RangesList') )[[chr]]
# #              ,
# # matplot(type='line',
# #   unlist(sapply(all.acc.views,function(x){(x[1])}))
# #         )
# # 
# # coverage(toppeaks[[1]])
# # 
# #              #get the locations of these 5 peaks
# # 
# # #plot the max and the others a little
# # coverageplot(toppeaks,Views()i=3)
# # 
# # 
# # 
# # 
# #   #go through each site and ask 
# #     #much it's moved
# # 
# # #having determined that, let's pick our biggest and smallest library
# # #for the bigger, go through each peak, collecting locations and height, and number of peaks within 300bp
# # #now go through all other datasets and ask - how often does our peak 
# # #Do certain libraries have more shifting from place to place?
# # #How does the
# # 
# # 
# # ggplot(as.data.frame(y=maxmatrix[1,]))+geom_density()
# #              
# 
# 
# 
# 
# # 
# # #could maybe have groupings with ggplot, sequencing runs etc. etc.
# # accession.df
# # ggplot(accession.df,aes(y=accession.df$library.size,x=line,fill=run))+
# #     geom_bar(stat='identity')
# # 
# # acc<-names(cage.tag.rles)[[1]]
# # strand<-'pos'
# #lets produce an Rle list which counts the number of lines with tags at a site, rather than the number of tags


