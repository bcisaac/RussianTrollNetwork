#Construct mention network from Russian Troll Tweets Dataset (via fivethirtyeight.com
#Author: Chris Marcum @csmarcum <christopher.steven.marcum@gmail.com>
#Date: 1 August 2018
#Last modified: 2 August 2018

#Script assumes that you have downloaded the data
# from https://github.com/fivethirtyeight/russian-troll-tweets
###############################################################################
# Install and Load Trolls Libraries
###############################################################################
# install.packages("ggplot2")
library(ggplot2)
# install.packages("network")
library(network)
# install.packages("igraph")
library(igraph)
# install.packages("sna")
# library(sna)
# install.packages("mvtnorm")
library(mvtnorm)
# install.packages("xtable")
library(xtable)

# install.packages("CINNA")
# library(CINNA)
# install.packages("intergraph")
# library(intergraph)


source('MultivarALAAMalt.R')

###############################################################################
# Load Trolls Rdata
###############################################################################
load("trolls.Rdata")

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
# Creating RightTroll ALAAM Model
###############################################################################
# Extract righttrolls from troll data
##### 
trollssample.righttroll <- get.inducedSubgraph(
  trolls,
  v=which(trolls %v% "accountcategory"%in%c("RightTroll"))
  )

##### 
# Create component
##### 
#delete attributes with variable rows
delete.vertex.attribute(trollssample.righttroll,"region")
delete.vertex.attribute(trollssample.righttroll,"language")

# get largest component using igraph
rt <- asIgraph(trollssample.righttroll)
comp <- clusters(rt)
comp
component.righttroll <- rt %>%induced.subgraph(
  .,
  v=which(comp$membership == which.max(comp$csize))
)

# return to network package
trollssample.righttroll <- asNetwork(component.righttroll)

##### 
# Create adjacency matrix component (plot and count trolls)
##### 
# make an adjacency matrix of network
adj.righttroll <- as.matrix.network.adjacency(trollssample.righttroll)
m <- nrow(adj.righttroll) # check how many rows are there 
m # this number is important, it is different from n in upper data
gplot(adj.righttroll)

##### 
# Prepare ALAAM covariates
##### 
# Creating context-specific covariates for covs table that will be used for ALAAM ()
accountcategory <- get.vertex.attribute(trollssample.righttroll,"accountcategory")
activity <- get.vertex.attribute(trollssample.righttroll,"activity")
language <- get.vertex.attribute(trollssample.righttroll,"language")
maxfollowing <- get.vertex.attribute(trollssample.righttroll,"maxfollowing")
minfollowing <- get.vertex.attribute(trollssample.righttroll,"minfollowing")
maxfollowers <- get.vertex.attribute(trollssample.righttroll,"maxfollowers")
minfollowers <- get.vertex.attribute(trollssample.righttroll,"minfollowers")

# Creating monadic covariates for covs.righttroll table that will be used for ALAAM
out.degree <-matrix( rowSums(adj.righttroll), m, 1) # number of ties sent
in.degree <- matrix( colSums(adj.righttroll) , m, 1 ) # number of ties received
rec.ties <-  matrix( rowSums(adj.righttroll * t(adj.righttroll) ), m , 1) # number of ties that are mutual
in.two.star <- matrix( choose(in.degree,2),m,1) #  in-stars refecting dispersion in popularity
out.two.star <- matrix( choose(out.degree,2),m,1) #  out-stars refecting dispersion in activity
mix.two.star <- in.degree*out.degree - rec.ties # correlation between indegree and outdegree  
in.three.star <- matrix( choose(in.degree,3),m,1) # furhter measure of in-degree heterogeneity
out.three.star <- matrix( choose(out.degree,3),m,1) # furhter measure of out-degree heterogeneity
triangles <- rowSums( adj.righttroll* (adj.righttroll %*% t(adj.righttroll) )  ) # embedded in transitive triads

## Forming a new covs
covs.righttroll <- cbind(
  activity,
  maxfollowing, 
  minfollowing,
  maxfollowers,
  minfollowers,
  in.degree,
  out.degree, 
  rec.ties,
  in.two.star,
  out.two.star,
  mix.two.star,
  in.three.star,
  out.three.star,
  triangles)
colnames(covs.righttroll) <- c(
  "activity",
  "maxfollowing", 
  "minfollowing",
  "maxfollowers",
  "minfollowers",
  "indegree",
  "outdegree",
  "reciprocation" ,
  "instar",
  "outstar",
  "twopath",
  "in3star",
  "out3star",
  "transitive")
head(covs.righttroll)

##### 
# Create dependent variable "following growth"
##### 
# creating measure of dv
fg.righttroll <- maxfollowers - minfollowers
summary(fg.righttroll) # the mean is 2874,but the median is 30. 
str(fg.righttroll)
fg.righttroll[fg.righttroll<2874] <- 0
fg.righttroll[fg.righttroll>=2874] <- 1
summary(fg.lefttroll)
str(fg.righttroll)

##### 
# Run and Test ALAAM Models
##### 
res.rt <- BayesALAAM(y = fg.righttroll,           # dependent variable
                    ADJ = adj.righttroll,           # network
                    covariates = covs.righttroll[,c(1,6,7,8,11,14)],   # covariates
                    directed = TRUE,    # directed / undirected network
                    Iterations = 6850,   # number of iterations
                    saveFreq = 500,     # print and save frequency
                    PropSigma = Propsigma )

## you can improve the mixing by using a better proposal covariance
## which can be used as an argument PropSigma to BayesALAAM. 
Propsigma <- cov(res.rt$Thetas)

##### 
# Temporary Model Storage
##### 
# good models w/ isolates and median
res.rt1 <- res.rt #6500
  # ESS - intercept, outdegree, twopath
  # SACF10 - intercept, contagion, outdegree, twopath
  # SACF30 -
res.rt2 <- res.rt #6800
# ESS - 
# SACF10 - outdegree, reciprocation
# SACF30 - intercept, indegree, outdegree, reciprocation
res.rt3 <- res.rt #6900
# ESS - 
# SACF10 - intercept, activity
# SACF30 - intercept, activity, outdegree, twopath
res.rt4 <- res.rt #6850
# ESS - 
# SACF10 - reciprocation
# SACF30 - activity, reciprocation

# good models w/ component and mean
res.rt5 <- res.rt
# ESS - intercept, twopath, transitive
# SACF10 - intercept, indegree, twopath, transitive
# SACF30 - intercept, activity, indegree, twopath, transitive
res.rt6 <- res.rt
# ESS -
# SACF10 - contagion, outdegree, reciprocation, twopath, transitive
# SACF30 - contagion, indegree, outdegree, reciprocation, twopath, transitive
res.rt7 <- res.rt # 6600
# ESS - twopath, transitive
# SACF10 - activity, outdegree, twopath, transitive
# SACF30 - contagion, activity, indegree, outdegree, reciprocation, twopath, transitive
res.rt8 <- res.rt # 6850
# ESS - activity
# SACF10 - 
# SACF30 - contagion, activity, transitive

rtmodel <- res.rt8

#####
saveRDS(rtmodel, file = "./.rtmodel.rds")
###############################################################################
# Creating LeftTroll ALAAM Model
###############################################################################
# Extract lefttrolls from troll data
##### 
trollssample.lefttroll <- get.inducedSubgraph(
  trolls,
  v=which(trolls %v% "accountcategory" %in% c("LeftTroll"))
  )

##### 
# Create component
#####
delete.vertex.attribute(trollssample.lefttroll,"region")
delete.vertex.attribute(trollssample.lefttroll,"language")
lt <- asIgraph(trollssample.lefttroll)

comp <- components(lt)
comp
component.lefttroll <- lt %>%induced.subgraph(
  .,
  v=which(comp$membership == which.max(comp$csize))
  )

#return to a network package
trollssample.lefttroll <- asNetwork(component.lefttroll)

##### 
# Create adjacency matrix component (plot and count trolls)
#####
adj.lefttroll <- as.matrix.network.adjacency(trollssample.lefttroll)
m <- nrow(adj.lefttroll) # check how many rows are there 
m # this number is important, it is different from n in upper data
gplot(adj.lefttroll)

##### 
# Prepare ALAAM covariates
#####
# Creating context-specific covariates for covs table that will be used for ALAAM
library(network)
accountcategory <- get.vertex.attribute(trollssample.lefttroll,"accountcategory")
names <- get.vertex.attribute(trollssample.lefttroll,"vertex.names")
activity <- get.vertex.attribute(trollssample.lefttroll,"activity")
language <- get.vertex.attribute(trollssample.lefttroll,"language")
maxfollowing <- get.vertex.attribute(trollssample.lefttroll,"maxfollowing")
minfollowing <- get.vertex.attribute(trollssample.lefttroll,"minfollowing")
maxfollowers <- get.vertex.attribute(trollssample.lefttroll,"maxfollowers")
minfollowers <- get.vertex.attribute(trollssample.lefttroll,"minfollowers")

# Creating monadic covariates for covs.righttroll table that will be used for ALAAM
out.degree <-matrix( rowSums(adj.lefttroll), m, 1) # number of ties sent
in.degree <- matrix( colSums(adj.lefttroll) , m, 1 ) # number of ties received
rec.ties <-  matrix( rowSums(adj.lefttroll * t(adj.lefttroll) ), m , 1) # number of ties that are mutual
in.two.star <- matrix( choose(in.degree,2),m,1) #  in-stars refecting dispersion in popularity
out.two.star <- matrix( choose(out.degree,2),m,1) #  out-stars refecting dispersion in activity
mix.two.star <- in.degree*out.degree - rec.ties # correlation between indegree and outdegree  
in.three.star <- matrix( choose(in.degree,3),m,1) # furhter measure of in-degree heterogeneity
out.three.star <- matrix( choose(out.degree,3),m,1) # furhter measure of out-degree heterogeneity
triangles <- rowSums( adj.lefttroll* (adj.lefttroll %*% t(adj.lefttroll) )  ) # embedded in transitive triads

## Forming a new covs
covs.lefttroll <- cbind(
  activity,
  maxfollowing, 
  minfollowing,
  maxfollowers,
  minfollowers,
  in.degree,
  out.degree, 
  rec.ties,
  in.two.star,
  out.two.star,
  mix.two.star,
  in.three.star,
  out.three.star,
  triangles)
colnames(covs.lefttroll) <- c(
  "activity",
  "maxfollowing", 
  "minfollowing",
  "maxfollowers",
  "minfollowers",
  "indegree",
  "outdegree",
  "reciprocation" ,
  "instar",
  "outstar",
  "twopath",
  "in3star",
  "out3star",
  "transitive")
head(covs.lefttroll)

##### 
# Create dependent variable "following growth"
#####
fg.lefttroll <- maxfollowers - minfollowers
summary(fg.lefttroll) # the mean is 2330/2750, but the median is 702/736 
str(fg.lefttroll)
fg.lefttroll[fg.lefttroll<2750] <- 0
fg.lefttroll[fg.lefttroll>=2750] <- 1
summary(fg.lefttroll)
str(fg.lefttroll)

##### 
# Run and Test ALAAM Models
#####
res.lt <- BayesALAAM(y = fg.lefttroll, # median-based dependent variable
                     ADJ = adj.lefttroll,           # network
                     covariates = covs.lefttroll[,c(1,6,7,8,11,14)],   # covariates
                     directed = TRUE,     # directed / undirecred network
                     Iterations = 4000,   # number of iterations
                     saveFreq = 500)   # print and save frequency

## you can improve the mixing by using a better proposal covariance
## which can be used as an argument PropSigma to BayesALAAM. 
Propsigma <- cov(res.lt$Thetas)

##### 
# Temporary Model Storage
#####
# good models 
res.lt1 <- res.lt #6700
# ESS - twopath, transitive
# SACF10 - activity, indegree, twopath, transitive
# SACF30 - everything but contagion
res.lt2 <- res.lt #6800
# ESS - 
# SACF10 - outdegree, reciprocation
# SACF30 - intercept, indegree, outdegree, reciprocation
res.lt3 <- res.lt #6900
# ESS - 
# SACF10 - intercept, activity
# SACF30 - intercept, activity, outdegree, twopath
res.lt4 <- res.rt #6850
# ESS - 
# SACF10 - reciprocation
# SACF30 - activity, reciprocation

ltmodel <- res.lt #6000

#####
saveRDS(ltmodel, file = "./.ltmodel.rds")
###############################################################################
# Load and Analyze Models
###############################################################################
# Read models
model <- readRDS("./.rtmodel.rds")
model <- readRDS("./.ltmodel.rds")

# Plot the MCMC output in trace plots 
plot(ts(model$Thetas))

# Obtain the model results
write.res.table(burnin=1, # should be set sufficiently high
                datamat=model$Thetas, # the result from BayesALAAM
                thin=1, # should be set so that SACF is sufficiently low, 
                # important for Confidence Intervals
                tabname=NULL) # the name appended to the table that is saved

# Draw posterior samples from model for GOF tests
sim.rt <- get.gof.distribution(NumIterations=500, # number of vectors to draw
                               res=model, # the ALAAM estimation object that contains model and results
                               burnin=100, # no. iterations discarded from GOF distribution
                               thinning = 1000) #, # no. iterations between sample points
# contagion ='none') # should be the same as for model fitted

# Create GOF table from simulated sample object
gof.table(obs.stats= sim.rt$stats, # observed statistics included not 
          # fitted statistics
          sim.stats= sim.rt$Sav.gof, # simulated goodness-of-fit statistics
          name.vec= sim.rt$gof.stats.names, # names of statistics calculate, 
          # not all will be used if undirected
          tabname='ALAAMGoFRightTroll', # name of file saved
          pvalues=TRUE, # posterior predictive p-values
          # save.tab ='csv', # save a csv file or a LaTex file
          directed=TRUE)
