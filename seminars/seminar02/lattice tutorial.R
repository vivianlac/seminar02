# load 'lattice'
library(lattice)

# load mini dataset from saved R object
kDat <- readRDS('./data/GSE4051_MINI.rds')
str(kDat)

table(kDat$devStage)

table(kDat$gType)

with(kDat, table(devStage, gType))

#############

## Scatterolots ##

#############

# NOTE: lattice uses general form: ' Y ~ x ' ("y twiddle x")

# lattice function for scatterplots is 'xyploit()'
xyplot(eggBomb ~ crabHammer, kDat)

xyplot(poisonFang ~ crabHammer, kDat)

# how to see two y variables plotted against x var?  Use "extended formula interface"
xyplot(eggBomb + poisonFang ~ crabHammer, kDat,
       auto.key = TRUE)

# above, 'auto.key = TRUE' produces legend automatically at the top!!

# putting each response, y var, in its own scatterplot for side-by-side comparison:
xyplot(eggBomb + poisonFang ~ crabHammer, kDat,
       outer = TRUE, grid = TRUE)

# above, 'outer = TRUE' causes info to be presented in two pannels.  'Grid = TRUE' adds a background grid for the plots

# identifying which points are from WT mice vs Nrl knockouts:
xyplot(eggBomb + poisonFang ~ crabHammer, kDat,
       outer = TRUE, gird = TRUE,
       groups = gType, auto.key = TRUE)

# reshaping the data from the previous figures such that the panels correspond to a level of a factor/unique combo of levels of two or more factors
nDat <-
  with(kDat,
       data.frame(sidChar, sidNum, devStage, gType, crabHammer,
                  probeset = factor(rep(c("eggBomb", "poisonFang"), each = nrow(kDat))),
                  geneExp = c(eggBomb, poisonFang)))
str(nDat)

# remaking of previous plot with the newly reshaped data
xyplot(geneExp ~ crabHammer | probeset, nDat,
       grid = TRUE,
       groups = gType, auto.key = TRUE)

# NOTE: syntaax ' y ~ x | z ' requests a scatterplot of y versus x for each level of factor z

# pracice plot: show developmental stage via color
xyplot(geneExp ~ crabHammer | probeset, nDat,
       grid = TRUE,
       groups = devStage, auto.key = TRUE)

###################

## Stripplot ##

###################

oDat <-
  with(kDat,
       data.frame(sidChar, sidNum, devStage, gType,
                  probeset = factor(rep(c("crabHammer", "eggBomb",
                                          "poisonFang"), each = nrow(kDat))),
                  geneExp = c(crabHammer, eggBomb, poisonFang)))
str(oDat)

# generate stripplot (a univariate scatterplot)
stripplot(~ geneExp, oDat)

# let's split things out for the different genes
stripplot(probeset ~ geneExp, oDat)

# add jitter in horizontal position
stripplot(probeset ~ geneExp, oDat, jitter.data = TRUE)

# put genes into different panels
stripplot(~ geneExp | probeset, oDat,
          layout = c(nlevels(oDat$probeset), 1))

# distplay info about WT vs knockout:
stripplot(~ geneExp | probeset, oDat,
          layout = c(nlevels(oDat$probeset), 1),
          groups = gType, auto.key = TRUE)

# exploring gene expression changes over the course of development
stripplot(geneExp ~ devStage, oDat)

# separate panel per gene
stripplot(geneExp ~ devStage | probeset, oDat,
          layout = c(nlevels(oDat$probeset), 1))

# add back the genotype info
stripplot(geneExp ~ devStage | probeset, oDat,
          layout = c(nlevels(oDat$probeset), 1),
          groups = gType, auto.key = TRUE)

# adding averages
stripplot(geneExp ~ devStage | probeset, oDat,
          layout = c(nlevels(oDat$probeset), 1),
          groups = gType, auto.key = TRUE, grid = TRUE,
          type = c('p', 'a'))

#above, 'p' specifies the data as points on the plot / 'a' referes to getting the average of each category and joining them by a line


##################

## Density Plot ##

#################

densityplot(~ geneExp, oDat)

#separate panels
densityplot(~ geneExp | gType, oDat,
            grid = TRUE)

# grouping by genotype
densityplot(~ geneExp, oDat,
            groups = gType, auto.key = TRUE)

# further customization of the plot with 'bw' and 'n'
jBw <- 0.2
jn <- 400
densityplot(~ geneExp, oDat,
            groups = gType, auto.key = TRUE,
            bw = jBw, n = jn,
            main = paste("bw =", jBw, "n = ", jn))

#################

## Boxplot ##

################

# boxplots in lattice from 'bwplot()' function for "box and whiskers plot"
bwplot(geneExp ~ devStage, oDat)

# separate panels
bwplot(geneExp ~ devStage | gType, oDat)

# violin plot
bwplot(geneExp ~ devStage, oDat,
       panel = panel.violin)

##############

## Heat Maps ##

##############

prDat <- read.table("./data/GSE4051_data.tsv.txt")
str(prDat, max.level = 0)

prDes <- readRDS("./data/GSE4051_design.rds")
str(prDes)

set.seed(1)
(yo <- sample(1:nrow(prDat), size = 50))

hDat <- prDat[yo, ]
str(hDat)

# heatmap expects a matrix, not a data.frame... some must convert 'hDat' and transpose
hDat <- as.matrix(t(hDat))
rownames(hDat) <- with(prDes,
                       paste(devStage, gType, sidChar, sep="_"))
str(hDat)

# heatmap!
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8))

# fix the color scheme
heatmap(hDat, Rowv = NA, Colv = NA, col = cm.colors(256),
        scale="none", margins = c(5, 8))

# load RColorBrewer
library(RColorBrewer)
display.brewer.all()

jGraysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))
jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))

# heatmapping in grey with RColorBrewer
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8),
        col = jGraysFun(256))

# heatmapping with blue-purple palette
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8),
        col = jBuPuFun(256))

# default rendering for heatmao
heatmap(hDat, margins = c(5, 8), col = jBuPuFun(256))

# allowing scaling within column
heatmap(hDat, col = jBuPuFun(256), margins = c(5, 8), scale=c("column"))

# using 'heatmap.2()' function from 'gplots' package
install.packages("gplots")

library(gplots)
heatmap.2(hDat, col = jGraysFun, trace = "none")

# adjust color of heatmap
heatmap.2(hDat, col = jBuPuFun, trace = "none")

##################

Overplotting

#################

set.seed(924)
(yo <- sample(1:ncol(prDat), size = 2))

y <- prDat[[yo[1]]]
z <- prDat[[yo[2]]]
str(y)

str(z)

xyplot(y ~ z, asp = 1)

# use 'smoothScatter()' function
smoothScatter(y ~ z, asp = 1)

# 'xyplot() function' in lattice can produce a similar plot by specifying a smoothScatter-type of panel function
xyplot(y ~ z, asp = 1, panel = panel.smoothScatter, nbin = 150)

# add-on package 'hexbin' implements hexagonal binning
library(hexbin)
hexbinplot(y ~ z)

###############

## Plot Matrix ##

##############

# take slightly larger sample of columns
set.seed(3)
(yo <- sample(1:ncol(prDat), size = 4))

pairDat <- subset(prDat, select = yo)
str(pairDat)

# using base function 'pairs()'
pairs(pairDat)

# combing 'pairs()" with 'smoothScatter()' for better result
pairs(pairDat,
      panel = function(...) smoothScatter(..., add=TRUE))

# using 'splom()' function from "lattice"
splom(pairDat)

# using "splom()' from lattice + 'smoothScatter'-type panel function
splom(pairDat, panel = panel.smoothScatter, raster = TRUE)

# using hexplom()
hexplom(pairDat)

