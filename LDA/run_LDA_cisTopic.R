# Run LDA with collapsed gibbs sampler (Griffiths & Steyvers (2004)) for provided dataset and parameters in cisTopic: https://github.com/aertslab/cisTopic
# 2 datasets: 1) simulated melanoma dataset provided by cisTopic (100 cells, 112701 regions); 2) 10x 5kPBMC dataset provided by 10x and also used in cisTopic (5335 cells, 97997 regions)
# Extracted from the runModels function in cisTopic package

# the number of topics can be customized. if length(topic)>1, it is suggested to us SNOW for parallelization.

library(lda) # core lda function
library(parallel) 
library(doSNOW) # parallelisation for running model for different number of topics simultaneously
library(plyr) #import llpy from plyr
library(Matrix)

# default parameter setting in cisTopic (The selected parameters are adapted from Griffiths & Steyvers (2004))
topic=c(2, 10, 20, 30, 40, 50)
seed=123
iterations = 500
burnin = 250
alpha = 50
beta=0.1

#####################################################################
######   Build LDA model for the similated melanoma dataset   #######
#####################################################################

# load data
cellList <- readRDS(file = 'data_simulated/cellList_simulated.rds')
regionList <- readRDS(file = 'data_simulated/regionList_simulated.rds')

# set parameters for lda model for simulated melanoma dataset
## For this dataset, #Topic=10 is the best model which gives the highest log likelihood
topic=c(2, 5:15, 20, 25) # set the number of topics
seed=987
nCores=length(topic)
burnin = 120
iterations = 150

# prepare for parallelisation
cl <- makeCluster(nCores)
registerDoSNOW(cl)
clusterEvalQ(cl, library(lda)) # load package lda for each progress
clusterExport(cl, c("cellList", "topic", "regionList", "iterations", "burnin", "alpha", "beta"), envir=environment())
clusterSetRNGStream(cl, seed)


# build lda model on simulated dataset -using parallelisation
set.seed(seed)
system.time(
  models_lda_simulated <- parLapply(cl,topic, function(t) lda.collapsed.gibbs.sampler(documents = cellList, K = t, vocab = regionList, num.iterations = iterations, alpha = alpha/t, eta = beta, burnin = burnin, compute.log.likelihood = TRUE)[-1])
)


# build lda model on simulated dataset - if length(topic)==1, no need to parallelise
#set.seed(seed)
#models_lda_simulated_melanoma <- suppressWarnings(llply(.data=topic, .fun = function(t) lda.collapsed.gibbs.sampler(cellList, K = t, regionList, num.iterations=iterations, alpha=alpha/t, eta=beta, compute.log.likelihood = TRUE, burnin=burnin)[-1], .progress = progress_text(char = ".")))
                                                        
print(length(models_lda_simulated))
print(attributes(models_lda_simulated[[1]]))
                                               
saveRDS(models_lda_simulated, file = 'data_simulated/models_lda_sim.rds')



######################################################################
###########   Build LDA model for 10x 5kPBMC dataset   ###############
######################################################################

#load data
cellList <- readRDS(file = 'data_10x/cellList_10x5kpbmc.rds')
regionList <- readRDS(file = 'data_10x/regionList_10x5kpbmc.rds')

# set parameters for lda model for 10x 5kPBMC dataset
# For this dataset, likelihood increases monotonically for the selected number of topics; Topic number=30 best balances model complexity and log likelihood;                                                
topic=c(2, 5, 10, 15, 20, 25, 30, 35, 40) # set the number of topics
seed=987
nCores=length(topic)
burnin = 120
iterations = 150

# prepare for parallelisation                                                        
cl <- makeCluster(nCores)
registerDoSNOW(cl)
clusterEvalQ(cl, library(lda)) # load package lda for each progress
clusterExport(cl, c("cellList", "topic", "regionList", "iterations", "burnin", "alpha", "beta"), envir=environment())
clusterSetRNGStream(cl, seed)                                                        
                                                        
# build lda model on 10xPBMC dataset -using parallelisation
set.seed(seed)
system.time(
  models_lda_10xPBMC <- parLapply(cl,topic, function(t) lda.collapsed.gibbs.sampler(documents = cellList, K = t, vocab = regionList, num.iterations = iterations, alpha = alpha/t, eta = beta, burnin = burnin, compute.log.likelihood = TRUE)[-1])
)
                                                        
                                                        
# build lda model for 10kPBMC dataset - if length(topic)==1, no need to parallelise
#set.seed(seed)                                                       
#models_lda_pbmc <- suppressWarnings(llply(.data=topic, .fun = function(t) lda.collapsed.gibbs.sampler(cellList, K = t, regionList, num.iterations=iterations, alpha=alpha/t, eta=beta, compute.log.likelihood = TRUE, burnin=burnin)[-1], .progress = progress_text(char = ".")))
                                     
print(length(models_lda_10xPBMC))
print(attributes(models_lda_10xPBMC[[1]]))

saveRDS(models_lda_10xPBMC, file = 'data_10x/models_lda_10x.rds')