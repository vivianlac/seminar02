# load ggplot2 package
library(ggplot2)

# Check out the 'geometries' functions avaialble!
apropos("^geom_")

# Check out the 'stat' functions available:
apropos("^stat_")

# Check out the 'scale' functions available:
apropos("^scale_")

# Separating factors:
# one factor:
facet_wrap

#two factor:
facet_grid


# load mini dataset from saved R object
kDat <- readRDS('./data/GSE4051_MINI.rds')
str(kDat)

table(kDat$devStage)
table(kDat$gType)

with(kDat, table(devStage, gType))

#########

## qplot function ##

#########

#qplot is a function in 'ggplot2' for quick plots
qplot(crabHammer, eggBomb, data = kDat)

#########

## Scatterplots ##

#########

# specify layer with ggplot()
p <- ggplot(kDat, aes(x = crabHammer, y = eggBomb))
str(p)

# now let's add geometries layer for simple scatter plot
(p <- p + geom_point())

# adding statistic layer to plot... a smoothing line
(p <- p + stat_smooth())

# the above input generates a plot using the default 'loess' function, a local regression model, and shaded ribbon for standard error

# we can also customize some default settings by adding more layers
(p <- p + theme_bw() + 
    xlab("Expression of crabHammer") +
    ylab("Expression of eggBomb") +
    ggtitle("Scatterplot for expression levels"))

## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.

# now let's plot both eggBomb and poisonFang against crabHammer as we did with 'lattice'... this involves some reshaping of the data
nDat <- with(kDat,
             data.frame(sidChar, sidNum, devStage, gType, crabHammer,
                        probeset = factor(rep(c("eggBomb", "poisonFang"),
                                             each = nrow(kDat))),
                        geneExp = c(eggBomb, poisonFang)))
str(nDat)

# specifying color in geometries (applies to ALL data!)
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) + geom_point())

# NOTE: this statement will work as well (uses 'aes' function)
(p <- ggplot(nDat, aes(crabHammer, geneExp)) + geom_point(aes(color = probeset)))

# Try adding a smoothing lune now
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) +
    geom_point() +
    stat_smooth(se = F))
# (note: se = F will turn off the display of standard error ribbon)

# overruling inheritance of previous defining of the groups by specifying another aesthetics in new layer
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) +
    geom_point() +
    stat_smooth(se = F, aes(group = 1)))

# if we want to plot 'poisonFang ~ crabHammer' and 'eggBomb ~ crabHammer' in separated pamnels, we can use facetting:
(p <- ggplot(nDat, aes(crabHammer, geneExp)) +
    geom_point() +
    facet_wrap(~ probeset))
# Note: for separating panels with two variables, see 'facet_grid()

# use different colors for the two genotypes...
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = gType)) +
    geom_point() +
    facet_wrap( ~ probeset))

# use different colors for the diff development stages:
(p <- ggplot(nDat, aes(crabHammer, geneExp, color = devStage)) +
    geom_point() +
    facet_wrap( ~ probeset))

#######

## Stripplot ##

#######

# in ggplot2, Stripplots still use 'geom_point()' ... just one of its coordinates is mapped to a 'factor' rather than a quantitative variable
oDat <-
  with(kDat,
       data.frame(sidChar, sidNum, devStage, gType,
                  probeset = factor(rep(c("crabHammer", "egbomb", 
                                          "poisonFang"), each = nrow(kDat))),
                  geneExp = c(crabHammer, eggBomb, poisonFang)))
str(oDat)

# plot expression level of each gene
(p <- ggplot(oDat, aes(geneExp, probeset)) +
    geom_point())

# we can also add jitter
(p <- ggplot(oDat, aes(geneExp, probeset)) + 
    geom_point(position = position_jitter(height = 0.1)))

# let's explore gene expression changes over the course of development
(p <- ggplot(oDat, aes(devStage, geneExp)) + 
    geom_point())

# Show different genes in separate panels
(p <- p + facet_wrap(~ probeset))

# Add genotype information 
(p <- p + aes(color = gType))

# Add averages with 'stat_summary"
(p <- p + stat_summary(fun.y = mean, geom = "point", shape = 4, size =4))

# the argument 'fun.y' will use the function you feed in... in this case mean' to summarize y for every x.  Alternatively, you can use fun.data if you want a summary on the entire dataset

#########

## Density Plots ##

#########

# thre are two functions 'ggplot2' uses for density plots: 'geom_density' and 'stat_density'

# geom_density method:
(p <- ggplot(oDat, aes(geneExp)) + geom_density())

# stat_density method:
(p <- ggplot(oDat, aes(geneExp)) +
    stat_density(geom = "line", position = "identity"))

# if you want a more similar presentation to the one we created with "lattice", i.e. adding data piunts at the bottom etc.:
(p <- ggplot(oDat, aes(geneExp)) +
    stat_density(geom = "line", position = "identity") +
    geom_point(aes(y = 0.05), position = position_jitter(height = 0.005)))

# Change bandwidth with 'adjust' argument in 'stat_density'
(p <- ggplot(oDat, aes(geneExp)) +
    stat_density(geom = "line", position = "identity", adjust = 0.5) +
    geom_point(aes(y = 0.05), position = position_jitter(height = 0.005)))

# separate panels for different genotype
(p <- p + facet_wrap(~ gType))

# use different colors for different genotyoe
(p <- ggplot(oDat, aes(geneExp, color = gType)) +
    stat_density(geom = "line", position = "identity") +
    geom_point(aes(y = 0.05), position = position_jitter(height = 0.005)))

#############

## Boxplot ##

#############

# 'geom_boxplot' function is available for boxplots in ggplot2
(p <- ggplot(oDat, aes(devStage, geneExp)) +
   geom_boxplot())

# separate by genotypes
(p <- p + facet_wrap(~ gType))

# violinplot is a hybrid of densityplot and histogram
(p <- ggplot(oDat, aes(devStage, geneExp)) +
    geom_violin())


###################################################

##### OVERPLOTTING AND PLOT MATRIX ######

##################################################

# load full data matrix
prDat <- read.table("./data/GSE4051_data.tsv.txt")
str(prDat, max.level = 0)

# loads an object named 'prDes'
prDes <- readRDS("./data/GSE4051_design.rds")
str(prDes)

# pick two samples at random to plot against each other
set.seed(2)
(yo <- sample(1:ncol(prDat), size = 2))

bDat <- data.frame(y = prDat[[yo[1]]], z = prDat[[yo[2]]])
str(bDat)

(p <- ggplot(bDat, aes(z,y)) +
    geom_point())

# reducing transparency of the data point with 'alpha argument'
(P <- ggplot(bDat, aes(z, y)) +
    geom_point(alpha = 0.1))

(p <- ggplot(bDat, aes(z, y)) +
    stat_density2d())

(p <- ggplot(bDat, aes(z, y)) +
    stat_density2d(geom = "tile", contour = F, aes (fill = ..density..)) +
    scale_fill_gradient(low = "white", high = "blue"))

# 'stat_binhex' function
install.packages("hexbin")
library(hexbin)
(p <- ggplot(bDat, aes(z,y)) +
    stat_binhex())

# for pairwise scatterplots, you can use the new 'ggpairs' function in 'GGally'
set.seed(3)
(yo <- sample(1:ncol(prDat), size = 4))

pairDat <- subset(prDat, select = yo)
str(pairDat)

# install 'GGally"
install.packages("GGally")
library(GGally)
(p <- ggpairs(pairDat))

#######################

## Heatmap ##

#######################

library(RColorBrewer)

# set seed so that we have exactly reproducible results
set.seed(1)

# choose 50 probes out of the 30k to work with
yo <- sample(1:nrow(prDat), size = 50)
hDat <- prDat[yo, ]
colnames(hDat) <- with(prDes, paste(devStage, gType, sidChar, sep = "_"))

# transform the data to tall format
prDatTall <- data.frame(sample = rep(colnames(hDat), each = nrow(hDat)),
                        probe = rownames(hDat),
                        expression = unlist(hDat))

# create a blue -> purple palette
jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
paletteSize <- 256
jBuPuPalette <- jBuPuFun(paletteSize)

# heatmap!
ggplot(prDatTall, aes(x = probe, y = sample, fill = expression)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_tile() +
  scale_fill_gradient2(low = jBuPuPalette[1],
                       mid = jBuPuPalette[paletteSize/2],
                       high = jBuPuPalette[paletteSize],
                       midpoint = (max(prDatTall$expression) +
                                     min(prDatTall$expression)) / 2,
                       name = "Expression")

####################

## HOMEWORK

####################

set.seed(1)

probe20 <- sample(1:nrow(prDat), size = 20)
heatDat <- prDat[probe20, ]
colnames(heatDat) <- with(prDes, paste(devStage, gType, sidChar, sep = "_"))

# transform data to tall format
prDataTall <- data.frame(sample = rep(colnames(heatDat), each = nrow(heatDat)),
                         probe = rownames(heatDat),
                         expression = unlist(heatDat))

# create a blue -> purple palette
jBuPurFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
paletteSizee <- 256
jBuPuPalettee <- jBuPurFun(paletteSizee)

# heatmap!
ggplot(prDataTall, aes(x = probe, y = sample, fill = expression)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_tile() +
  scale_fill_gradient2(low = jBuPuPalettee[1],
                       mid = jBuPuPalettee[paletteSizee/2],
                       high = jBuPuPalettee[paletteSizee],
                       midpoint = (max(prDataTall$expression) +
                                     min(prDataTall$expression)) /2,
                       name = "Expression")