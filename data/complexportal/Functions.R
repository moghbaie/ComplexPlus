# Mehrnoosh Oghbaie
# 08/30/2018
# Repository for all the functions


######################################################################################
# Either download or install the required library from CRAN or bioconductor
#######################################################################################

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  if (length(bioconductor.packages) > 0) {
    source("http://bioconductor.org/biocLite.R")
  }
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      biocLite(p) 
      library(p, character.only=T)
    }
  }
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p) 	
      library(p, character.only=T)  	
    }	
  }
}

###########################################################################################################################
# Writing a function that get the url of each file (Protein Complex) and integrate interaction data to one csv file
###########################################################################################################################

Convert2CSV <- function(filename,dir){
  complex <- strsplit(basename(filename),".xml")[[1]]
  data <- xmlParse(filename)
  xml_data <- xmlToList(data)
  
  ## reading the interaction list
  
  complexID <- as.list(xml_data[["entry"]][["interactionList"]][["interaction"]][["xref"]])
  complex_name <- complexID[grepl("CPX",complexID)]$secondaryRef[["id"]]
  binding <- as.list(xml_data[["entry"]][["interactionList"]][["interaction"]][["inferredInteractionList"]])
  
  if (length(binding)!=0){
    list <- data.frame(matrix(ncol=2,nrow=0))
    for (bind in binding){
      interact1 <- as.integer(bind[1]$participant$participantFeatureRef)
      interact2 <- as.integer(bind[2]$participant$participantFeatureRef)
      interact <- c(interact1,interact2)
      list <- rbind(list,interact)}
    colnames(list) <- c("id_A","id_B")
    
    ## reading the interactor list
    nodes <-as.list(xml_data[["entry"]][["interactorList"]])
    ref_list <- data.frame(matrix(ncol=3,nrow=0))
    for(node in nodes){
      ref <- as.integer(unlist(node$.attrs["id"]))
      protein <- node$names$shortLabel
      protein_fullname <- ifelse(is.null(node$names$fullName),NA, node$names$fullName)
      nodel <- cbind(as.integer(ref),protein,protein_fullname)
      ref_list <- rbind(ref_list,unname(nodel))
    }
    colnames(ref_list) <- c("ref","protein","protein_fullname")
    ## reading the mapping between interactors and featured interactors in interaction list
    links <- as.list(xml_data[["entry"]][["interactionList"]][["interaction"]][["participantList"]])
    link_list <- data.frame(matrix(ncol=2, nrow=0))
    for(link in links){
      intRef <-as.integer(link$interactorRef)
      featurelist <- link$featureList
      feat <-c()
      if(!is.null(featurelist)){
        for(feature in featurelist){
          feat <-rbind(feat,as.integer(feature$.attrs["id"]))
        }
        linkRef <- cbind(feat,intRef)
        link_list <- rbind(link_list, linkRef)
      }
    }
    colnames(link_list) <- c("id","ref")
    
    ## matching and merging the data
    link_list$protein <- ref_list$protein[match(link_list$ref,ref_list$ref)]
    link_list$protein_fullname <- ref_list$protein_fullname[match(link_list$ref,ref_list$ref)]
    list$protein_A <- link_list$protein[match(list$id_A,link_list$id)]
    list$protein_B <- link_list$protein[match(list$id_B,link_list$id)]
    list$protein_fullname_A <- link_list$protein_fullname[match(list$id_A,link_list$id)]
    list$protein_fullname_B <- link_list$protein_fullname[match(list$id_B,link_list$id)]
    list$complexID <- complex_name
    ## Saving each complex separately
    write.csv(list, paste0(dir,complex,".csv"))
    
    
  }else{print(complex_name)}
  ##return the list of interactors for each complex
  return(list)
}


