#### Mehrnoosh Oghbaie
### 07/26/2019
## Following functions are run in this document
# 1. Build networks
# 2. Build mapping matrix includes correction
# 3. Tune parameters
# 4. Leave-One-Out Cross validation
## * functions that take relatively long time to run
## ** functions that take couple of days to run

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

# install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","ggthemes","magrittr","igraph","sqldf","diffusr",
                   "readxl","ggplot2", "stringr","XML","RCurl","xml2","gsubfn","httr", "rvest",
                   "plotly","AnnotationDbi","GO.db", "R6","Matrix", "ggforce","wordspace")

bioconductor.packages <- c("PSICQUIC","biomaRt","graphite")
source("Functions/functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)
source("Functions/build_networks.R")
source("Functions/tune_parameters.R")
source("Functions/build_interface_db.R")
## 1. Building the adjacency and bipartite mapping networks
### 1.1 Import geneMania and Bioplex interactions
### 1.2 Import ComplexPortal database
### 1.3 Download 

net = Networks$new()$
  getGeneManiaNetworks()$
  extractPubmedComplexPortal()$
  getDirectInteraction()$
  saveOldProteinInteraction()$
  getPhysicalInetractionfromDB()$
  getIdentifierMapping()$
  buildOldDirectGeneGraph()$
  buildOldDirectProteinNet()$
  buildColocalizationNet()$
  net$buildPredictionNet()$
  buildGeneticNet()$
  buildCoexpressionNet()$
  buildPathwayNet()$
  buildxLinkNet()$
  unionGeneList()$
  unionProteinList()$
  buildDiseaseSimilarityNet()$
  buildGoBipartiteMappingNet()$
  saveLists2DB()

#Saving network object
save(net, file="backup/networksBackup.RData")

## 2. Tuning parameters
### 2.1 Mapping gene/disease - level networks to protein-level using bipartite mapping network
### 2.2 Integrating multiple-kernels to a monoplex network using Rige-regression in 10 folds each with 10% missing
### 2.3 Make different combination of networks with average parameters from previous step
### 2.4 Leave-One-Out Cross validation for different combination of network
### 2.5 Build a Multiplex with the network built in previous step and disease similarity

tune = Tune$new()$
  loadBackupGeneProtein()$
  mapGeneDisease2Protein()$
  makeComplexPortalDB()$
  tuneLayers()$
  takeAverageParameters()$
  buildCombinationLayers()$
  leaveOneOutCVMonoplex()$
  #buildMultiplex()$
  #leaveOneOutCVMultiplex()$
  visualizeResult()

save(tune, file="backup/tuneBackup.RData")

## 3. Saving interface databases in InterfaceDB object
interfaceDB = InterfaceDB$new()$
  loadDB()
interfaceDB$directInteraction()

save(interfaceDB, file="../interface/db/interfacedbBackup.RData") 
