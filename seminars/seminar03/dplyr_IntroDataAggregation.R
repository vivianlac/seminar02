# load dplyr package
library(dplyr)

# load Gapminder data
gd_url <- "http://tiny.cc/gapminder"
gdf <- read.delim(file = gd_url)
str(gdf)