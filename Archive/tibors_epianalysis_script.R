# working directory
setwd('/g/furlong/pakozdi/')

# required packages
library(IRanges)
library(multicore)

# genomic annotation (FBv547)
genedb <- read.delim('database/genomic_annotation/FlyBase_dm3_v5.47.txt')

# select reasonable chromosomes
genedb <- genedb[!grepl('Het|U|M', genedb$chr), ]

# mesoderm specific RNA-seq data (Zeitlinger lab)
exprs <- read.delim('database/gene_expression/Zeitlinger_RNA_FBv547_RPKM.txt')

# summarize 6-8h replicates
exprs <- rowMeans(exprs[,c('Zeitlinger_Muscle_RNA_6.8h_NA_R1', 'Zeitlinger_Muscle_RNA_6.8h_NA_R2')])

# convert to log scale
exprs <- log10(exprs)

# calculate probability density of the rpkm values
pdens <- density(exprs)

# find local maxima and select 2 highest pts.
lmax <- which(diff(sign(diff(pdens$y)))==-2) + 1
lmax <- lmax[which(pdens$y[lmax] %in% tail(sort(pdens$y[lmax]), 2))]

# select local minima between the two maximas
threshold <- pdens$x[which.min(pdens$y[lmax[1]:lmax[2]]) + lmax[1]]

# visualize gene expression bimodality and local minima splitting the active from non-active genes
plot(pdens, xlab = expression(paste(Log[10], ' RPKM')), ylab = 'Relative frequency', main = 'RNA-seq 6-8h')
abline(v = threshold, col = 2)

# define all active genes within reasonable chromosomes
agene <- names(exprs[exprs >= threshold])
agene <- agene[agene %in% genedb$gid]

# BiTS-ChIP H3K4me3 peaks at 6-8h
hpeaks <- read.table('01_Polycomb_Repression/results/HistMod_MACS/6-8h/H3K4me3/K4me3_B24-B53_6-8h_peaks.bed', skip = 1)
hpeaks[,1] <- paste0('chr', hpeaks[,1])

# ranges of H3K4me3 peaks
hranges <- IRangesList(sapply(unique(hpeaks[,1]), simplify = F, function(x) IRanges(start = hpeaks[hpeaks[,1] == x,2], end = hpeaks[hpeaks[,1] == x,3])))

# total isoform ranges for all genes
isoranges <- IRangesList(sapply(unique(genedb[,'chr']), simplify = F, function(x) {

	IRanges(start = genedb[genedb[, 'chr'] == x, 'trstart'], end = genedb[genedb[, 'chr'] == x, 'trstop'], names = genedb[genedb[, 'chr'] == x, 'trid'])
}))

# all-to-all overlap between isoforms and H3K4me3 peaks
total_overlap <- unlist(countOverlaps(isoranges, hranges[names(isoranges)], type = 'any'))
names(total_overlap) <- unlist(lapply(names(total_overlap), function(x) strsplit(x, '[.]')[[1]][2]))

# check the isoform overlap of active genes
goverlap <- unlist(mclapply(agene, mc.cores = 8, function(x) {

	# select isoforms
	sel <- genedb[genedb$gid == x, 'trid']

	any(total_overlap[sel] == 1)
}))

table(goverlap)['TRUE'] / sum(table(goverlap))

overlap_dist <- genedb[genedb$gid %in% agene[goverlap], c('gid', 'glen')]
no_overlap_dist <- genedb[genedb$gid %in% agene[!goverlap], c('gid', 'glen')]

# remove redundancy
overlap_dist <- overlap_dist[!duplicated(overlap_dist), ]
no_overlap_dist <- no_overlap_dist[!duplicated(no_overlap_dist), ]

plot(density(log10(overlap_dist[,2])), xlab = expression(paste(Log[10], ' gene length (bp)')), ylab = 'Relative frequency', main = '')
lines(density(log10(no_overlap_dist[,2])), col = 2)
legend('topright', c('With H3K4me3', 'Without H3K4me3'), col = 1:2, lwd = 3, bty = 'n')

wilcox.test(overlap_dist[,2], no_overlap_dist[,2])

# script end