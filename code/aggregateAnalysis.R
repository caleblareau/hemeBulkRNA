library(tidyverse)
library(GenomicRanges)

g10X <- read.table("../annotation/genes.tsv", stringsAsFactors = FALSE)
metaData <- read.table("../annotation/lookup.txt", header = TRUE, stringsAsFactors = FALSE)

# Add a ordered, unique sample ID
v <- metaData$Source
idxn <- lapply(unique(v), function(ele){
  nvec <- which(ele == v)
  names(nvec) <- as.character(1:length(nvec))
  nvec
}) %>% unlist() %>% sort() %>% names()

metaData$UID <- paste0(metaData$Source, "_", idxn)

# Import the HT-Seq counts data
shortNames <- gsub('.{17}$', '', list.files("../htseq", full.names = FALSE))
fullnames <- list.files("../htseq", full.names = TRUE)
rawCounts <- sapply(fullnames, function(file){read.table(file)[,2]})

# Do ENSG gene names; cut out some junk
n <- dim(g10X)[1]
mdf <- merge(read.table(fullnames[1])[1:n,], g10X, by.x = "V1", by.y = "V1", sort = FALSE)
rawCounts <- rawCounts[1:n,]
colnames(rawCounts) <- shortNames
rownames(rawCounts) <- mdf[,1]
df <- rawCounts[g10X[,1],]

# Translate names 
namesvec <- metaData$UID
names(namesvec) <-  metaData$Base
colnames(df) <- namesvec[colnames(df)]
write.table(df, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE,
            file = "../output/individualHemeRNAcounts.txt")

# Collapse counts across replicates
namesvec <- metaData$Source
names(namesvec) <-  metaData$UID
df2 <- df
colnames(df2) <- namesvec[colnames(df)]
cc <- sapply(unique(colnames(df2)), function(x) rowSums(df2[,colnames(df2) == x]))
write.table(cc, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE, 
            file = "../output/groupHemeRNAcounts.txt")

# Import GTF; make data frame of spanning gene coordinates 
gtf <-rtracklayer::import("../annotation/hg19_10X.gtf.gz")
gs <- tapply(start(gtf), list(mcols(gtf)@listData$gene_id), min)
ge <- tapply(start(gtf), list(mcols(gtf)@listData$gene_id), max)
chrs <- data.frame(chr = seqnames(gtf), gene = mcols(gtf)@listData$gene_id) %>% unique()
remove(gtf)
coords <- data.frame(gs, ge); coords$gene <- rownames(coords)
mergedDF1 <- merge(chrs, coords, by = "gene")

# Pairwise Limma-Voom
library(limma)
library(edgeR)

doLimmaGetBed <- function(inGroup, prop = 0.1, pad = 100000){
  
  # Run standard limma - voom code
  groups <- as.numeric(namesvec[colnames(df)] == inGroup)
  z <- DGEList(counts = data.matrix(df2), group = groups)
  z <- calcNormFactors(z)
  design <- model.matrix(~groups)
  v <- voom(z,design,plot=FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)
  
  # Make a final data frame for export as bed file 
  results <- as.data.frame(topTable(fit, number = round(nrow(df2)*prop), sort.by = "t"))
  rdf <- merge(results, mergedDF1, by.x = "row.names", by.y  = "gene")
  rdf$start <- rdf$gs - pad
  rdf$end <- rdf$ge + pad
  rdf <- merge(rdf, g10X, by.x = "Row.names", by.y = "V1")
  rdf$chrchr <- paste0("chr", as.character(rdf$chr))
  write.table(rdf[,c("chrchr", "start", "end", "V2")], paste0("../bedFiles/", inGroup, ".25June2017.bed"), 
                    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  inGroup
}

lapply(unique(namesvec), doLimmaGetBed)

