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
      BiocManager::install(p) 
      library(p, character.only=T)
    }
  }
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p,character.only = T, dependencies = T ) 	
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
    ref_list <- data.frame(matrix(ncol=2,nrow=0))
    for(node in nodes){
      ref <- as.integer(unlist(node$.attrs["id"]))
      protein <- node$names$shortLabel
      nodel <- cbind(as.integer(ref),protein)
      ref_list <- rbind(ref_list,nodel)
    }
    colnames(ref_list) <- c("ref","protein")
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
    list$protein_A <- link_list$protein[match(list$id_A,link_list$id)]
    list$protein_B <- link_list$protein[match(list$id_B,link_list$id)]
    ## Saving each complex separately
    write.csv(list, paste0(dir,complex,".csv"))
    list$complex <- complex
    
  }
  ##return the list of interactors for each complex
  return(list)
}


##########################################################################################################################
## Copying file

my.file.copy <- function(from, to) {
  #if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  if(!file.exists(to)){
    file.copy(from,to)
  }
}

###########################################################################################################################
## PSCIQUIC related functions
###########################################################################################################################


#########################################################################################################
# A function that take pubmed, imex, doi and mint IDs from one column and create separate column for them

fillSource <- function(dx){
  if(grepl("\\|",dx["publicationID"])){
    x <- strsplit(unlist(dx["publicationID"]), "\\|")[[1]]
  }else{
    x<- dx["publicationID"]
  }
  if(any(grepl("pubmed",unlist(x)))){
    dx["pubmed"] <- gsub("pubmed:","",x[grepl("pubmed",x)])
  } 
  if(any(grepl("imex",x))){
    dx["imex"] <- gsub("imex:","",x[grepl("imex",x)])
  }
  if(any(grepl("doi",x))){
    dx["doi"] <- gsub("doi:","",x[grepl("doi",x)])
  } 
  if(any(grepl("mint",x))){
    dx["mint"] <- gsub("mint:","",x[grepl("mint",x)])[1]
  }
  return(dx)
}
###################################################################################################################
# Get content of the parenthesis
fillMethod <- function(dx){
  x <- unlist(as.character(dx[["detectionMethod"]]))
  return(paste0(unlist(regmatches(x, gregexpr("(?<=\\().*?(?=\\))", x, perl=T))[[1]]), collapse="|"))
}

###################################################################################################################
# Get confidence score value from PSCIQUIC output
# score:
# intact-miscore:
# mentha-score:
# author-confidence:Z-score
# author score:
fillconfidenceScore <- function(dx){
  if(grepl("\\|",dx["confidenceScore"])){
    x <- strsplit(unlist(dx["confidenceScore"]), "\\|")[[1]]
  }else{
    x<- dx["confidenceScore"]
  }
  if(any(grepl("intact-miscore:",unlist(x)))){
    dx["confidence"] <- gsub("intact-miscore:","",x[grepl("intact-miscore:",x)])
  } 
  else if(any(grepl("mentha-score:",x))){
    dx["confidence"] <- gsub("mentha-score:","",x[grepl("mentha-score:",x)])
  }
  else if(any(grepl("score:",x)&!grepl("author score:",x))){
    dx["confidence"] <- gsub("score:","",x[grepl("score:",x)&!grepl("author score:",x)])
  } 
  return(dx)
}

####################################################################################################################
# Convert different gene ID types to gene name
####################################################################################################################

transfer2Gene <- function(x){
  if(grepl("entrez gene/locuslink:", unlist(x[c("A")]))){
    x["A"] <- entrez_geneName_mapping_hs[match(as.integer(gsub("entrez gene/locuslink:","",unlist(x[c("A")]))),entrez_geneName_mapping_hs[,"entrezgene"]),"external_gene_name"]
    x["B"] <- entrez_geneName_mapping_hs[match(as.integer(gsub("entrez gene/locuslink:","",unlist(x[c("B")]))),entrez_geneName_mapping_hs[,"entrezgene"]),"external_gene_name"]
  } else if(grepl("innatedb:",unlist(x[c("A")]))){
    x["A"] <- ensembl_geneName_mapping_hs[match(gsub("ensembl:","",unlist(x[c("altA")])),ensembl_geneName_mapping_hs[,"ensembl_gene_id"]),"external_gene_name"]
    x["B"] <- ensembl_geneName_mapping_hs[match(gsub("ensembl:","",unlist(x[c("altB")])),ensembl_geneName_mapping_hs[,"ensembl_gene_id"]),"external_gene_name"]
  } else if(grepl("uniprotkb:",unlist(x[c("A")]))){
    x["A"] <- protein_geneName_mapping_hs[match(gsub("uniprotkb:","",unlist(x[c("A")])),protein_geneName_mapping_hs[,"uniprot_gn"]),"external_gene_name"]
    x["B"] <- protein_geneName_mapping_hs[match(gsub("uniprotkb:","",unlist(x[c("B")])),protein_geneName_mapping_hs[,"uniprot_gn"]),"external_gene_name"]
  }
  return(x)
}

#####################################################################################################################
getGeneName <- function(x,col="altA") {
  lst <- strsplit(as.character(x[[col]]),"\\|")[[1]]
  uniprots <- gsub("uniprotkb:", "",lst[grepl("uniprotkb:",lst)])
  
  genes <- gsub("entrez gene/locuslink:", "",lst[grepl("entrez gene/locuslink:",lst)])
  y <- genes
  if(length(uniprots)>0){
    protein_gene <- biomaRt::getBM(attributes=c('external_gene_name', 'uniprot_gn_id'), 
                                   filters = 'uniprot_gn_id', 
                                   values =  uniprots, 
                                   mart = mart)
    y <-  c(protein_gene[["external_gene_name"]],y)
  }
  return( paste0(unique(y),collapse="|"))
}

x <- tbl_direct4[1,]

apply(tbl_direct4[1:10,] ,1, getGeneName, col="altA")
######################################################################################################################
### Merging link in a network 

Merge_parallel <- function(x){
  y <- x %>%
    dplyr::group_by(Gene_A,Gene_B) %>%
    dplyr::summarise(score_max = max(Weight, na.rm=T), Source = paste(unique(Source[Source!=0]), collapse="|"))
  return(y)
} 

##########################################################################################################
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# # check if the HGNC symbols that we will use to generate the networks
####  exist in Biomart.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
check.HGNC.symbols <- function(GeneNames){
  
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  ensemble_HGNC <- getBM(attributes=c('hgnc_symbol'), 
                         filters='hgnc_symbol', values=GeneNames, mart=ensembl)
  
  return(ensemble_HGNC$hgnc_symbol)
}

##############################################################################################################
#### Get pathway data using graphite
##############################################################################################################
GetPathwayData <- function (db.name, db.specie){
  db <- pathways(db.specie, db.name)
  OneDb.HGNC <- data.frame()
  
  for (i in 1:length(db)) {
    current.pathway <- db[[i]]
    current.pathway.info <- try(data.frame(src = current.pathway@protEdges$src,
                                           dest = current.pathway@protEdges$dest,
                                           direction = current.pathway@protEdges$direction,
                                           type = current.pathway@protEdges$type,
                                           name = current.pathway@title),silent=TRUE) 
    OneDb.HGNC <- rbind(OneDb.HGNC,current.pathway.info)
  }
  
  OneDb.HGNC<-OneDb.HGNC[complete.cases(OneDb.HGNC),] #remove NA rows
  
  #transform the directed edge into undirected edge
  OneDb.HGNC[OneDb.HGNC == "directed"] <- "undirected" 
  
  #create a sorted labelled per edge
  OneDb.HGNC$name1 <- pmin(as.vector(OneDb.HGNC$src), as.vector(OneDb.HGNC$dest))
  OneDb.HGNC$name2 <- pmax(as.vector(OneDb.HGNC$src), as.vector(OneDb.HGNC$dest))
  OneDb.HGNC$inter.label <- paste(OneDb.HGNC$name1,OneDb.HGNC$name2, sep="_")
  
  #select the unique inter.label = unique edge
  OneDb.HGNC <- subset(OneDb.HGNC,!(duplicated(OneDb.HGNC$inter.label)))
  
  # select the 2 columns of interest
  NetworkRaw<-as.matrix(OneDb.HGNC[,c(1,2)])
  colnames(NetworkRaw) <- c("SymbolA", "SymbolB")
  return(NetworkRaw)
  
}


################################################################################################################
#' Build an igraph network	
#' Take df with 2 columns, format "ncol" (one line per interaction, the 2 interactors are separated by a tab)
#' returns a simplified network as an igraph object.	
build.network <- function (raw.network) {	
  net <- igraph::graph.edgelist(raw.network, directed=F)	
  # simplify
  net <- igraph::simplify(net,remove.multiple = TRUE, remove.loops = TRUE)
  
  # We additionally delete isolated nodes.
  isolates <- which(degree(net, mode = c("all")) == 0) - 1
  net <- delete.vertices(net, names(isolates))
  
  return (net) 	
}
################################################################################################################
### Normalize a sparse matrix

normalizeM<- function(X){
  rowsum <- Matrix::rowSums(X, na.rm = TRUE)
  colsum <- Matrix::colSums(X, na.rm = TRUE)
  maxDegree <- pmax(rowsum,colsum)
  D12 <- Diagonal(x=sqrt(maxDegree))
  InD12 <- solve(D12)
  XP_normalize <- InD12%*%X%*%InD12
  return(XP_normalize)
}
#################################################################################################################
### Conjugate gradient

# A is a symmetric matrix of nxn
# b and x are matrix of nx1 or just a vector of n

conjgrad = function(A,b,x0){
  r <- b - A %*% x0
  p=r
  rsold <- t(r)%*%r
  threshold <- 1e-10
  i=1
  rsnew = t(r)%*% r
  while((i<=dim(b)[1] && rsnew[1,1]>threshold)){
    Ap = A%*%p
    alpha = rsold/(t(p)%*% Ap)
    x0 = x0 +alpha[1,1]*p
    r = r-alpha[1,1]*Ap
    rsnew = t(r)%*% r
    p <- r +(rsnew/rsold)[1,1]*p;
    rsold = rsnew
    i <- i+1;
  }
  return(x0)
}

##################################################################################################################
#####' convert uniprot ID to gene name
getGeneName <- function(x){
  list <- getBM(attributes=c('uniprot_gn', 'external_gene_name'), 
                filters = 'uniprot_gn', 
                values = strsplit(unlist(unname(x['Identifiers'])),"\\|")[[1]], 
                mart = mart)
  return(paste(unlist(unname(list["external_gene_name"])),collapse="|"))
}

#################################################################################################################
#' Enrich list of uniprot ID based complex protein distribution
enrichComplex <- function(significant_uniprot,Complex_protein,geneName=1){
  ### list of annotated proteins
  k <- length(significant_uniprot)
  
  ### list of proteins in human complexes
  N <-length(unique(unlist(Complex_protein$uniprotID)))
  
  selected_complex_count <- Complex_protein %>%
    filter(uniprotID %in% significant_uniprot) %>%
    group_by(ComplexName) %>%
    summarise(count=as.numeric(as.character(count[1])),
              selected_count= as.numeric(n()),
              GeneRatio=selected_count/count,
              complete_members=Identifers[1],
              Identifiers=paste(uniprotID, collapse="|")) %>%
    mutate(pValue = phyper(q=as.numeric(selected_count),m=as.numeric(count), n =N-as.numeric(count), k=k , lower.tail=FALSE))
  
  
  selected_complex_count$p.adj <- p.adjust(selected_complex_count$pValue, "BH")
  selected_complex_count_test <- selected_complex_count %>% filter(p.adj<5e-02&as.numeric(count)>1) %>% arrange(p.adj,desc(selected_count))
  if(geneName==1){
    mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    selected_complex_count_test$geneName <- apply(selected_complex_count_test,1, getGeneName)
  }
  
  return(selected_complex_count_test)
}

#######################################################################################
## Build gene and proteins networks
# Create empty columns
createEmptyCol <- function(tbl){
  tbl$pubmed <- NA
  tbl$imex <- NA
  tbl$doi <- NA
  tbl$mint <- NA
  return(tbl)
}


convertEnsembl2Swisspro <- function(ensembl_gene_id) {
  mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  return(biomaRt::getBM(attributes=c('ensembl_gene_id', 'uniprotswissprot'), 
                        filters = 'ensembl_gene_id', 
                        values =as.character(ensembl_gene_id), 
                        mart = mart))
}


convertUniprot2Swisspro <- function(uniprot_gn) {
  mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  return(biomaRt::getBM(attributes=c('uniprot_gn_id', 'uniprotswissprot'), 
                        filters = 'uniprot_gn_id', 
                        values =as.character(uniprot_gn), 
                        mart = mart))
}

convertSwisspro2ensembl <- function(uniprot_gn) {
  mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  return(biomaRt::getBM(attributes=c('uniprotswissprot', 'ensembl_gene_id'), 
                        filters = 'uniprotswissprot', 
                        values =as.character(uniprot_gn), 
                        mart = mart))
}


getPubmed <- function(x) {
  xx <- strsplit(as.character(x[["Publication.Identifier.s."]]), "\\|")[[1]]
  xxx <- xx[grepl("pubmed:",xx)]
  xxxx <- as.integer(gsub("pubmed:","",xxx))
  return(xxxx)
}

getWeight <- function(x) {
  xx <- strsplit(as.character(x[["Confidence.value.s."]]), "\\|")[[1]]
  xxx <- xx[grepl("intact-miscore:",xx)]
  xxxx <- as.numeric(gsub("intact-miscore:","",xxx))
  return(xxxx)
}

get_Gene_Name <- function(x,identifier_mappings){
  name <- identifier_mappings %>% dplyr::filter(Preferred_Name==x, Source=="Gene Name") %>% dplyr::select(Name)
  return(as.character(name$Name))
}
#################################################################################################################
### Normalize and map to target network

normalizeConvert<- function(M_Gene2Protein,X){
  XP <- t(M_Gene2Protein)%*% X %*% M_Gene2Protein
  rowsum <- Matrix::rowSums(XP, na.rm = TRUE)
  colsum <- Matrix::colSums(XP, na.rm = TRUE)
  maxDegree <- pmax(rowsum,colsum)
  D12 <- Matrix::Diagonal(x=sqrt(maxDegree))
  InD12 <- solve(D12)
  XP_normalize <- InD12%*%XP%*%InD12
  return(XP_normalize)
}

### Normalize
normalizeM<- function(X){
  rowsum <- Matrix::rowSums(X, na.rm = TRUE)
  colsum <- Matrix::colSums(X, na.rm = TRUE)
  maxDegree <- pmax(rowsum,colsum)
  D12 <- Matrix::Diagonal(x=sqrt(maxDegree))
  InD12 <- solve(D12)
  XP_normalize <- InD12%*%X%*%InD12
  return(XP_normalize)
}

### Fast Correlation
sparse.cor4 <- function(x){
  n <- nrow(x)
  cMeans <- Matrix::colMeans(x)
  covmat <- (as.matrix(Matrix::crossprod(x)) - n*Matrix::tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(Matrix::diag(covmat)) 
  cormat <- covmat/Matrix::tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

#########################
webScrapperOrpha <- function(x) {
  extract <- unique(c(ifelse(nchar(gsub(".*\\[(.*)\\].*", "\\1", strsplit(as.character(x[["Disease"]]) ,"\\|")[[1]]))>15,"",gsub(".*\\[(.*)\\].*", "\\1", strsplit(as.character(x[["Disease"]]) ,"\\|")[[1]])),
                      regmatches(x[["Disease"]], gregexpr("(?<=\\().*?(?=\\))", x[["Disease"]], perl=T))[[1]]))
  return(extract[grepl("^(EFO|Orpha|ORPHA)", extract)])
}



#####################
Random_Walk_Restart <- function(Network_Matrix, r,prox_vector ){
  
  ### We define the threshold and the number maximum of iterations for the randon walker.
  Threeshold <- 1e-10
  
  
  ### We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1
  
  restart_vector <-  prox_vector
  
  while(residue >= Threeshold){
    
    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(Network_Matrix %*% prox_vector) + r*restart_vector
    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
    iter <- iter + 1; 
  }
  return(prox_vector) 
}