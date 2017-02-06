# load packages
library(lattice)
library(ggplot2)
library(reshape2)

# load data to be analyzed
prDat <- read.table("./data/GSE4051_data.tsv.txt") 
str(prDat, max.level = 0)

prDes <- readRDS("./data/GSE4051_design.rds")
str(prDes)

# write function that takes probeset IDs as input and gives output as data.frame
library(plyr)

(luckyGenes <- c("1419655_at","1438815_at"))
prepareData <- subset(prDat, rownames(prDat) %in% luckyGenes)
prepareData <- data.frame(gExp = as.vector(t(as.matrix(prepareData))),
                          gene = factor(rep(rownames(prepareData), each = ncol(prepareData)),
                                        levels = luckyGenes))
jDat <- suppressWarnings(data.frame(prDes, prepareData))

str(jDat)

head(jDat)
tail(jDat)

# generate figure to check if function is working:
stripplot(gExp ~ devStage | gene, jDat,
          group = gType, jitter.data = TRUE,
          auto.key = TRUE, type = c('p', 'a'), grid = TRUE)

stripplot(gExp ~ devStage | gene, jDat, pch = 17, cex = 3,
          group = gType, jitter.data = TRUE,
          auto.key = TRUE, type = c('p', 'a'), grid = TRUE)

# make stripplot using ggplot2:


suppressWarnings(data.frame(prDes, prepareData))