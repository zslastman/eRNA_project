# load the library
library(spp);

# The following section shows how to initialize a cluster of 8 nodes for parallel processing
# see "snow" package manual for details.
library(snow)
cluster <- makeCluster(8);


chip.data <- read.bam.tags(chrom_bam_files[[1]]);
input.data <- read.bam.tags(input_bam_files[[1]]);

binding.characteristics <- get.binding.characteristics(chip.data,srange=c(50,500),bin=5,cluster=cluster);

print(paste("binding peak separation distance =",binding.characteristics$peak$x))

par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation");
abline(v=binding.characteristics$peak$x,lty=2,col=2)

select.informative.tags(reads,binding.characteristics)
smoothed.density <- get.smoothed.tag.density(reads,bandwidth=200,step=100,tag.shift=tag.shift);
writewig(smoothed.density,"example.density.wig","Example smoothed, background-subtracted tag density");

?get.smoothed.tag.density


tag.shift <- round(binding.characteristics$peak$x/2)
