## This script is used to get plots for preliminary quality assesment on your data. 
## Set the appropriate variables below and run on cluster environment (large memory requirements).
## The number of samples in your count file should be equal to the number of samples in your metafile

### Set Variables
count.file <- 'combined.counts'
meta.file <- 'ginty_neuron.csv'
job.name <- 'nmf'
baseDir="/n/home02/mmistry/ginty_neuron"
dataDir=paste(baseDir, "/htseq-count", sep="")
resultsDir=paste(baseDir, "/results", sep="")
metaDir=paste(baseDir, "/meta", sep="")

### Number of iterations 
n.runs= 50 

### Rank factor range
rank.range <- 4:8

### Variability cutoff
var.cutoff = 0.95 #

########### BEGIN SCRIPT ###########
library(CHBUtils)
library(grid)
library(gridExtra)
library(Biobase)
library(NMF)
library(genefilter)

# Load data
counts <- read.delim(file.path(dataDir, count.file), sep="\t", row.names=1, header=T)
meta <- read.delim(file.path(metaDir, meta.file), sep=",", row.names=1, header=T)
meta <- meta[colnames(counts),]

# Create eset
eset <- new("ExpressionSet", exprs=as.matrix(counts))
pData(eset) <- meta

# Filter by variance
# calculates a measure of variability of each row (row inter-quartile range)
# select a fraction (the var.cutoff argument) with most variability
var.eset <- varFilter(eset, var.cutoff=var.cutoff)


## NMF Preparation
### Estimating the factorization rank
var.eset$replicate <- NULL
nmf.brunet <- nmf(var.eset, rank.range, nrun=n.runs, method="brunet", seed = "random", .options='tp')
png(file.path(resultsDir, paste(job.name, '.quality.png', sep="")), width = 10, height = 10, units = 'in', res = 300)
plot(nmf.brunet)
dev.off()

# shuffle original data to look for overfitting
eset.rand <- randomize(var.eset)

# estimate quality measures from the shuffled data (use default NMF algorithm)
nmf.brunet.rand <- nmf(eset.rand, rank.range, nrun = n.runs, seed = 'random', .options='tp')

# plot measures on same graph
png(file.path(resultsDir, paste(job.name, '.overfitting.png', sep="")), width = 10, height = 10, 
    units = 'in', res = 300)
plot(nmf.brunet, nmf.brunet.rand)
dev.off()


## consensus matrix for each value of the factorization rank to see if the clusters (or consensus blocks) obtained correspond to the known cell types.
png(file.path(resultsDir, paste(job.name, '.consensus.png', sep="")), 
    width = 10, height = 10, units = 'in', res = 300)
consensusmap(nmf.brunet, annCol=var.eset,labCol=NA, labRow=NA)
dev.off()

# Algorithm compare: do this for each rank 

for (r in rank.range){
  nmf.multi.method <- nmf(var.eset, r, list("brunet", "KL", "lee","nsNMF"), nrun=n.runs, seed = "random", 
                          .options = "tp")
  png(file.path(resultsDir, paste('Rank', r, '.nmf.algcompare.png', sep="")), width = 10, 
      height = 10, units = 'in', res = 300)
  plot(nmf.multi.method, main="")
  dev.off()
}




