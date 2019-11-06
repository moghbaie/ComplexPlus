

#### Conjugate Gradient
conjgrad = function(A,b,x0){
  r <- b - A %*% x0
  p=r
  rsold <- t(r)%*%r
  threshold <- 1e-16
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

#########################################################################################
### Run random restart on edges
Random_Walk_Restart <- function(A, r,b ){
  
  ### We define the threshold and the number maximum of iterations for the randon walker.
  Threeshold <- 1e-4
  NetworkSize <- ncol(A)
  
  ### We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1
  
  #### We define the prox_vector(The vector we will move after the first RW iteration. We start from The seed. We have to take in account
  #### that the walker with restart in some of the Seed genes, depending on the score we gave in that file).
  
  prox_vector  <- b
  restart_vector <-  prox_vector
  
  while(residue >= Threeshold){
    
    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(A %*% prox_vector) + r*restart_vector
    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
    iter <- iter + 1; 
  }
  return(prox_vector) 
}
#########################################################################################
### Get protein membership 
getProteinComplexMembership <- function(result,Complex_protein,protein_list ){
  dz <- Complex_protein %>%
    dplyr::filter(uniprotID %in% toupper(result[["Protein"]]))%>%
    dplyr::group_by(uniprotID) %>% dplyr::summarise(paste0(ComplexName, collapse="|"))
  
  residue <- result[["Protein"]][!result[["Protein"]]%in% dz[["uniprotID"]]]
  dz2 <- data.frame(cbind(residue,rep("",length(residue))))
  colnames(dz) <- colnames(dz2) <- c("uniprotID","Complex")
  dz <- rbind(data.frame(dz),dz2)
  result[["Membership"]] <- as.character(dz[["Complex"]])[match(as.character(result[["Protein"]]),dz$uniprotID)]
  result[["count"]] <- ifelse(is.na(result[["Membership"]])|result[["Membership"]]=="", 0 ,stringr::str_count(result[["Membership"]], pattern = "\\|")+1)
  
  return(result)
}

#########################################################################################

#idmapping <- read.delim("F:/Network_Analysis/Interface/hm/identifier_mappings.txt")
#idmapping <- idmapping %>% mutate(Preferred_Name = as.character(Preferred_Name),
#                                  Name = as.character(Name),
 #                                 Source = as.character(Source))
#save(idmapping, file="db/idmapping.RData")
#load(file="db/idmapping.RData")

UniProt2Gene <- function(x){
  preferred<- idmapping%>% 
    dplyr::filter(Name == x) %>% 
    dplyr::select(Preferred_Name)
  Gene_Name <- idmapping %>% 
    dplyr::filter(Preferred_Name %in% preferred$Preferred_Name )%>% 
    dplyr::filter(Source=="Gene Name")%>%
    dplyr::select(Name)
  return(as.character(Gene_Name$Name[1]))
}
