## This script is used to get plots and feature lists once PARAMETERS HAVE BEEN SET
## Set the appropriate variables below and run on cluster environment (large memory requirements)
## The number of samples in your count file should be equal to the number of samples in your metafile (csv)

### Set Variables
count.file <- 'combined.counts'
meta.file <- 'ginty_neuron.csv'
job.name <- 'nmf'
baseDir="/n/home02/mmistry/ginty_neuron"
dataDir=paste(baseDir, "/htseq-count", sep="")
resultsDir=paste(baseDir, "/results", sep="")
metaDir=paste(baseDir, "/meta", sep="")

### Number of iterations 
n.runs= 1000

### Rank factor
rank <- 6

### Specify algorithm
algorithm <- 'lee'

### Variability cutoff
var.cutoff = 0.90


########### BEGIN SCRIPT ###########
library(NMF)
library(CHBUtils)
library(Biobase)
library(genefilter)
library(RColorBrewer)


# Load data
counts <- read.delim(file.path(dataDir, count.file), sep="\t", row.names=1, header=T)
meta <- read.delim(file.path(metaDir, meta.file), sep=",", row.names=1, header=T)
meta <- meta[colnames(counts),]

# Create eset
eset <- new("ExpressionSet", exprs=as.matrix(counts))
pData(eset) <- meta
fData(eset) <- data.frame(GeneSymbol=ann.counts[,'symbol'], row.names=row.names(ann.counts))

# Filter by variance
# calculates a measure of variability of each row (row inter-quartile range)
# select a fraction (the var.cutoff argument) with most variability
var.eset <- varFilter(eset, var.cutoff=var.cutoff)
var.eset$replicate <- NULL

# final run
nmf.final <- nmf(var.eset, rank, algorithm, nrun=n.runs, seed = "random", .options = "tp")

png('results/consensus_final.png', width = 10, height = 10, units = 'in', res = 300)
consensusmap(nmf.final, annCol=var.eset,labCol=NA, labRow=NA)
dev.off()


# show individual genes in relation to metagenes and samples
png('results/coefmap_final.png', width = 10, height = 10, units = 'in', res = 300)
coefmap(nmf.final, labCol=NA, labRow=NA, annCol=var.eset)
dev.off()

# show metagenes in relation to samples
png('results/basismap_final.png', width = 10, height = 10, units = 'in', res = 300)
basismap(nmf.final, subsetRow=TRUE)
dev.off()

# Get fit
bestfit <- fit(nmf.final) #extract the best fit NMf model
fS <- featureScore(bestfit) # get all the feature scores which measure specificity to which a gene contributes to a metagene
featureList <- extractFeatures(bestfit) # extract the features with the most specifciity for each of the metagenes

## get unique annotations for genes (NMF returns number referenced IDs) in metagenes
unique.metagenesymbols <- lapply(featureList, function(x) {
  genenames <- unique(unlist(fData(var.eset)[x,"GeneSymbol"]))
  return(genenames)
  })

## get number of unique genes in each metagene
numgenes <- unlist(lapply(unique.metagenesymbols, length))

# dataframe to present the metagene features that includes both the metagene annotations and the featurescores
for (n in 1:length(featureList)){
  metagene <- cbind(fData(var.eset)[featureList[[n]],], fS[featureList[[n]]]) 
  names(metagene)[ncol(metagene)] <- "featurescore"
  
  write.table(as.matrix(metagene), 
              file=file.path(resultsDir, paste("Rank", rank, "metagenes.", n, ".txt", sep="")), sep="\t", quote=F)
  
}


