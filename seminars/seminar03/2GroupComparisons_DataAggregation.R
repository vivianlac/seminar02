# load packages
library(knitr)
library(lattice)
library(ggplot2)
library(plyr)

prDat <- read.table("./data/GSE4051_data.tsv.txt")
str(prDat, max.level = 0)

prDes <- readRDS("./data/GSE4051_design.rds")
str(prDes)

###############################
## Two sample tests - one gene ##

# extract data for one gene and put in data.frame with experimental info
set.seed(987)
(theGene <- sample(1:nrow(prDat), 1))

pDat <- data.frame(prDes, gExp = unlist(prDat[theGene, ]))
str(pDat)

# data aggregation
aggregate(gExp ~ gType, pDat, FUN = mean)

# doing same as above but with 'plyr' function
ddply(pDat, ~ gType, summarize, gExp = mean(gExp))

# make stripplot to test our t-test result (using 'lattice')
stripplot(gType ~ gExp, pDat)

# same as above but with 'ggplot2'
ggplot(pDat, aes(x = gExp, y = gType)) + geom_point()


# Let's conduct two-sample t-test comparing WT to Nrl knockouts
t.test(gExp ~ gType, pDat)

# save this t-test result in an object
ttRes <- t.test(gExp ~ gType, pDat)
str(ttRes)

# extract t-statistic from saved t-test
ttRes$statistic

# extract p-value
ttRes$p.value

######

# Trying the above myself!
# extract data for one gene and put in data.frame with experimental info
set.seed(500)
(OneGene <- sample(1:nrow(prDat), 1))

pDat2 <- data.frame(prDes, gExp = unlist(prDat[OneGene, ]))
str(pDat2)

# data aggregation
aggregate(gExp ~ gType, pDat2, FUN = mean)

# make stripplot (lattice)
stripplot(gType ~ gExp, pDat2)

# two-sample t-test
ttRes2 <- t.test(gExp ~ gType, pDat2)
str(ttRes2)

ttRes2$statistic
ttRes2$p.value

##########################
#apply() for computing on rows and columns of matrices

kDat <- readRDS("./data/GSE4051_MINI.rds")

kMat <- as.matrix(kDat[c('crabHammer', 'eggBomb', 'poisonFang')])
str(kMat)

# computing median expression for specific gene by hand + using apply()
median(kMat[ ,1])

median(kMat[ , 'eggBomb'])

apply(kMat, 2, median)

# here's an alternative way to compute gene-specific medians:
apply(kMat, 2, quantile, probs = 0.5)

apply(kMat, 2, quantile, probs = c(0.25, 0.75))

# take min gene expression for each sample across these 3 genes
apply(kMat, 1, min)

colnames(kMat)[apply(kMat, 1, which.min)]

# special functions for computing row- and column-wise sums
rowSums(kMat)

all.equal(rowSums(kMat), apply(kMat, 1, sum))

colMeans(kMat)

all.equal(colMeans(kMat), apply(kMat, 2, mean))

#########################
#aggregate() for computing on groups of observations

# compute average expression of 'eggBomb' for different levels of 'devStage'
aggregate(eggBomb ~ devStage, kDat, FUN = mean)

# split data into grups based on combination of factors
aggregate(eggBomb ~ gType * devStage, kDat, FUN = mean)

# report range
aggregate(eggBomb ~ gType * devStage, kDat, FUN = range)

## DOING SAME TASKS with 'plyr' ##

ddply(kDat, ~ devStage, summarize, avg = mean(eggBomb))

ddply(kDat, ~ gType * devStage, summarize, avg = mean(eggBomb))

ddply(kDat, ~ gType * devStage, summarize,
      min = min(eggBomb), max = max(eggBomb))

###########################
## Two sample tests - a handful of genes
keepGenes <- c("1431708_a_at", "1424336_at", "1454696_at",
               "1416119_at", "1432141_x_at", "1429226_at" )
miniDat <- subset(prDat, rownames(prDat) %in% keepGenes)
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                      gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                    levels = keepGenes))
miniDat <- suppressWarnings(data.frame(prDes, miniDat))
str(miniDat)

stripplot(gType ~ gExp | gene, miniDat,
          scales = list(x = list(relation = "free")),
          group = gType, auto.key = TRUE)

# make same plots using 'ggplot2'
ggplot(miniDat, aes(x = gExp, y = gType, color = gType)) +
  facet_wrap(~ gene, scales="free_x") +
  geom_point(alpha = 0.7) +
  theme(panel.grid.major.x = element_blank())

# conducting t-test...
someDat <- droplevels(subset(miniDat, gene == keepGenes[1]))
t.test(gExp ~ gType, someDat)

# but how to scale this up to ALL SIX genes?? (capability of aggregate() has been outgrown...)
# use 'plyr'!!!!

#load 'plyr' package
library(plyr)
d_ply(miniDat, ~ gene, function(x) t.test(gExp ~ gType, x), .print = TRUE)

# we can use dlply() to retain everything in a new list with one component per probeset
ttRes <- dlply(miniDat, ~ gene, function(x) t.test(gExp ~ gType, x))
names(ttRes)

# if we knew in advance that we only wanted, say, the test stat and p-val, then...
ttRes <- ddply(miniDat, ~ gene, function(z) {
  zz <- t.test(gExp ~ gType, z)
  round(c(tStat = zz$statistic, pVal = zz$p.value), 4)
})
ttRes