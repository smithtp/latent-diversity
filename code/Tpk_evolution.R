#################################################
## This script is to estimate evolution of     ##
## traits across a phylogeny (Tpk, size, E)    ##
##
## Reproduces figure 3 and related analyses
##
## Tom Smith 2022
##
#################################################

graphics.off()
rm(list = ls())
## load libraries
library(phytools)
library(ggplot2)
library(geiger)
library(coda)

# read the time calibrated tree
time_tree <- read.tree("data/dated_tree_20220912_combined_new")

plot(time_tree)

# then drop the outgroup tip
tree <- drop.tip(time_tree, "Methanospirillum_hungatei_00_00_00")

plot.phylo(tree, type = "fan", show.tip.label = FALSE)

# now lets rearrange the tree somewhat so the branches are more sensibly placed

tree <- rotateNodes(tree, 61)
tree <- rotateNodes(tree, 62)
tree <- rotateNodes(tree, 63)
tree <- rotateNodes(tree, 70)
tree <- rotateNodes(tree, 89)
tree <- rotateNodes(tree, 92)
tree <- rotateNodes(tree, 93)
tree <- rotateNodes(tree, 95)
tree <- rotateNodes(tree, 100)
tree <- rotateNodes(tree, 101)
tree <- rotateNodes(tree, 108)
tree <- rotateNodes(tree, 110)
tree <- rotateNodes(tree, 113)


# now onto the exciting stuff
# of trying to estimate ancestral states of Tpk

# change the tip labels to be less messy
tiplabels <- tree$tip.label
# use gsub to seperate the species names and strain codes
species <- gsub('.{9}$', '', tiplabels)
species_nice <- gsub('_', ' ', species)
IDcode <- gsub('.*(?=.{8}$)', '', tiplabels, perl = T)
# similarly we can take the isolation and incubation temps from these codes
IncuTemp <- gsub('.{6}$', '', IDcode)
Iso1 <- gsub('.{3}$', '', IDcode)
IsoTemp <- gsub('.*(?=.{2}$)', '', Iso1, perl = T)

#bind into dataframe
metadata <- data.frame(species, IDcode, IncuTemp, IsoTemp)
metadata$species <- as.character(metadata$species)

tree$tip.label <- as.character(metadata$IDcode)

plot(tree)

# Create a data frame for the traits
# TPC data a little complicated as we should add Streptomyces in
summary_data <- read.csv("data/growth_summary.csv")
summary_data <- summary_data[!(summary_data$strain %in% c("30_30_01", "30_30_03", "30_30_06", "30_30_07")),]
strep_fits <- read.csv("data/summary_OD_R.csv")

summary_data$iso_temp <- as.character(summary_data$iso_temp)
strep_fits$iso_temp <- as.character(strep_fits$iso_temp)

TPC_data <- bind_rows(summary_data, strep_fits)

Tpk <- c()
Ea <- c()

for(i in 1:length(IDcode)){
  TPC_subs <- TPC_data[TPC_data$strain == IDcode[i],]
  
  # add a device to catch times when there was no model fit
  if (nrow(TPC_subs) > 0){
    Tpk_subs <- TPC_subs$T_pk_est_sch - 273.15 # fuck it, lets stick it back to Celsius
    Ea_subs <- TPC_subs$E_sch
    
    Tpk <- c(Tpk, Tpk_subs)
    Ea <- c(Ea, Ea_subs)
  }
  else {
    Tpk <- c(Tpk, NA)
    Ea <- c(Ea, NA)
  }
}

trait_frame <- data.frame(Tpk, Ea, IsoTemp, IncuTemp, row.names = IDcode)

# are any missing?

trait_frame[is.na(trait_frame$Tpk),]
trait_frame

# we'll take the Tpk by eye from these strains

trait_frame[20,]$Tpk <- 30 # 30_30_03
trait_frame[60,]$Tpk <- 25 # 50_RT_02
trait_frame[11,]$Tpk <- 25 # 30_RT_04
trait_frame[12,]$Tpk <- 25 # 21_RT_01
trait_frame[15,]$Tpk <- 25 # 30_RT_05

# also need a fix on 50_50_04 because its got a Tpk fit with no points after the curve
trait_frame[24,]$Tpk <- 55

# map Tpk to tree

Tpk_matrix <- data.matrix(trait_frame)[,1]

plot(tree)

# used fastAnc to estimate ancestral states
fit<-fastAnc(tree,Tpk_matrix,vars=TRUE,CI=TRUE)
fit
fit$CI95[1,]

lims<-c(floor(min(Tpk_matrix)*10)/10,ceiling(max(Tpk_matrix)*10)/10)

## projection of the reconstruction onto the edges of the tree
obj<-contMap(tree,Tpk_matrix,plot=FALSE)
plot(obj, type = "fan",legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.7,0.9))

## change colours to blue -> red
n<-length(obj$cols)
obj$cols[1:n]<-colorRampPalette(c("cyan", "orange", "red"), space="Lab")(n) 

# plot out a nice fan version
svg("results/Tpk_evolution.svg", width = 10, height = 10)
plot(obj, type = "fan", ftype = "off", legend = FALSE, xlim = c(-4000, 4000), lwd = 8)
add.color.bar(1500, obj$cols, title = "",
              lims = obj$lims, digits = 0, prompt = FALSE, x = -3800,
              y = -3000, lwd = 8, fsize = 1.5, subtitle = "1.5 BY")
# add clade labels
arc.cladelabels(text = "Firmicutes", node = 63, ln.offset = 1.05, lab.offset = 1.1, mark.node = FALSE, lwd = 2, cex = 2)
arc.cladelabels(text = "Actinobacteria", node = 99, ln.offset = 1.05, lab.offset = 1.1, mark.node = FALSE, lwd = 2, cex = 2)
arc.cladelabels(text = "Proteobacteria", node = 104, ln.offset = 1.05, lab.offset = 1.1, mark.node = FALSE, lwd = 2, cex = 2)
dev.off()

# need to reset the graphics to get the next plot to fit properly
graphics.off()


# update the trait phenogram so that
# its just the trait (not CIs) and so that
# different clades (phyla) are coloured differently

plotTree(tree,pts=F,node.numbers=T)

# so node 104 is proteobacteria, 63 is actinobacteria, 68 is firmicutes

# paint the clades in different states
paintedtree<-paintSubTree(tree,node=104,state="2")
paintedtree<-paintSubTree(paintedtree,node=99,state="3")
paintedtree<-paintSubTree(paintedtree,node=63,state="4")

# now let's plot using plotSimmap to ensure
# that the correct branches were painted

cols<-c("black","orange","red", "blue"); names(cols)<-1:4
plotSimmap(paintedtree,cols,pts=F,lwd=3,node.numbers=T)

# good, now plot the traitgram
svg("results/Tpk_phenogram_nolabs.svg", width = 10, height = 10)
par(xaxt="n",yaxt="n")
phenogram(paintedtree,x=Tpk_matrix,colors=cols,ftype = "off", ylab = "", xlab = "", lwd = 4)
par(xaxt="s",yaxt="s")
axis(1, xaxp = c(0, 3223.57, 4), lwd = 2, labels = FALSE)
axis(2, labels = FALSE, lwd = 2)
dev.off()


## now conduct some actual tests of phylogenetic heritability of Tpk

# Important to note that the trait needs to be normally distributed for lambda and K - so check normality of Tpk

ggplot(trait_frame, aes(x = Tpk)) + geom_histogram()
ggplot(trait_frame, aes(x = Tpk)) + geom_density()

# log transform the data for normality

trait_frame$logTpk <- log(trait_frame$Tpk)
trait_frame$recip_pk <- 1/(trait_frame$Tpk)

Tpk_matrix <- data.matrix(trait_frame)[,5]

# The first test is Blomberg’s K, which compares the variance of PICs to what we would expect under a Brownian motion model.
# K = 1 means that relatives resemble one another as much as we should expect under BM;
# K < 1 means that there is less “phylogenetic signal” than expected under BM, 
# while K > 1 means that there is more. 
# A significant p-value returned from phylosig tells you that there is significant phylogenetic signal - 
# that is, close relatives are more similar than random pairs of species.

# λ measures the similarity of the covariances among species to the covariances expected under Brownian motion; 
# whereas K might be more usefully thought of as a measure of the partitioning of variance. 
# If K>1 then variance tends to be among clades; while if K<1 then variance is within clades (with BM as reference). 
# The variance on K for a given process is quite large.

# Phylogenetic heritabilities between 0 and 1 reflect deviations from random evolution

phylosig(tree,Tpk_matrix,test=TRUE)

# Secondly we use pagels λ
# λ close to zero, means phylogenetic signal equivalent to that expected if the data arose on a star phylogeny 
# (that is, no phylogenetic signal). λ = 1 corresponds to a Brownian motion model; 0 < λ < 1 is in between.

phylosig(tree,Tpk_matrix,method="lambda",test=TRUE)

# Finally we can test more complicated models of trait evolution

# brownian motion model
fitBM<-fitContinuous(tree,Tpk_matrix)
# 'early-burst' model (EB) in which character change tends to be concentrated towards the base of the tree
fitEB<-fitContinuous(tree,Tpk_matrix,model="EB")
# Ornstein-Uhlenbeck model (OU), which is used to model trait evolution with the tendency towards a central value - 
# such as under constant stabilizing selection.
fitOU<-fitContinuous(tree,Tpk_matrix,model="OU")

# pagels lambda
fitlambda <- fitContinuous(tree,Tpk_matrix,model="lambda")

# check out the results

fitBM
fitEB
fitOU
fitlambda

aic.vals<-setNames(c(fitBM$opt$aicc,fitOU$opt$aicc,fitEB$opt$aicc),
                   c("BM","OU","EB"))
aic.vals

aic.w(aic.vals)
# Brownian motion model appears to best fit the data.
