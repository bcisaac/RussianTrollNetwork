#Construct mention network from Russian Troll Tweets Dataset (via fivethirtyeight.com
#Author: Chris Marcum @csmarcum <christopher.steven.marcum@gmail.com>
#Date: 1 August 2018
#Last modified: 2 August 2018

#Script assumes that you have downloaded the data
# from https://github.com/fivethirtyeight/russian-troll-tweets

# install.packages("ggplot2")
library(stringr)
library(sna)
###############################################################################
# Load Trolls Rdata
###############################################################################
load("trolls.Rdata")
filter.

###############################################################################
# Plot
###############################################################################
set.seed(20818)
png("RTRT.png",width=3000,height=2500,res=100,pointsize=24)
par(mar=c(3,1,3,1),xpd=TRUE)
gplot(trolls,vertex.col=as.color(trolls%v%"accountcategory"),displayisolates=FALSE,vertex.cex=abs(log(trolls%v%"maxfollowers"+2,base=10)-1),edge.col=rgb(0,0,0,.5),main="Russian Troll to Russian Troll Twitter Mention Network (n=1245) \n fivethirtyeight.com Data 08/01/18 \ by @csmarcum")
legend("bottom",legend=unique(trolls%v%"accountcategory"),pch=19,col=unique(as.color(trolls%v%"accountcategory")),title="Account Type",horiz=TRUE)
dev.off()

###############################################################################
# ALAAM
###############################################################################

library(networks)

trolls (network)
# create subgraph and filter trolls
name <- get.vertex.attribute(trolls, "vertex.names")
cat <- get.vertex.attribute(trolls, "accountcategory")
df <- data.frame(name = name, cat = cat)
# get a list of vertices of specific categories
# pass that list of names/vertex ids to get induced Subgraph
v <- df[(df$cat%in%c("RightTroll","LeftTroll","FearMonger")), "name"]
trolls1 <- get.inducedSubgraph(trolls, v)
trolls2 <- get.inducedSubgraph(trolls,v=which(trolls %v% "accountcategory"%in%c("RightTroll","LeftTroll","FearMonger","NewsFeed","HashtagGamer")))
# make an adjacency matrix of network
adj <- as.matrix.network.adjacency(trolls2)

gplot(trolls2)
# make a cov table
changeinfoll <- get.vertex.attribute(trolls2, "maxfollowers") - get.vertex.attribute(trolls2, "minfollowers")

colnames(covs) <- c("minfollowers",
                   "maxfollowers", ...)

head(covs)



# install.packages("mvtnorm")
library(mvtnorm)
# install.packages("xtable")
library(xtable)