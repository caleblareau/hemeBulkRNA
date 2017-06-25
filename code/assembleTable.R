library(tidyverse)

g10X <- read.table("../data/genes.tsv", stringsAsFactors = FALSE)
metaData <- read.table("../data/lookup.txt", header = TRUE, stringsAsFactors = FALSE)

# Add a ordered, unique sample ID
v <- metaData$Source
idxn <- lapply(unique(v), function(ele){
  nvec <- which(ele == v)
  names(nvec) <- as.character(1:length(nvec))
  nvec
}) %>% unlist() %>% sort() %>% names()

metaData$UID <- paste0(metaData$Source, "_", idxn)

# Import the HT-Seq counts data
shortNames <- gsub('.{17}$', '', list.files("../data/htseq/unpairedOK", full.names = FALSE))
fullnames <- list.files("../data/htseq/unpairedOK", full.names = TRUE)
rawCounts <- sapply(fullnames, function(file){read.table(file)[,2]})

# Do common gene names; cut out some junk
n <- dim(g10X)[1]
mdf <- merge(read.table(fullnames[1])[1:n,], g10X, by.x = "V1", by.y = "V1", sort = FALSE)
rawCounts <- rawCounts[1:n,]
colnames(rawCounts) <- shortNames
rownames(rawCounts) <- mdf[,3]
df <- rawCounts[g10X[,2],]

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
