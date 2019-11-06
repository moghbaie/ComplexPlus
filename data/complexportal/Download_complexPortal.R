# Mehrnoosh Oghbaie
# 10/05/2018
# reading files in ftp server and converting interaction data to table


## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

source("Functions.R")
CRAN.packages <- c("XML","RCurl","xml2","httr","readxl","readr")
bioconductor.packages <- c()
#source("Functions/Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)


#################################################################################
## homo_sapiens 
## setting the url

local_url <- "F:/ComplexPortal/xml/"


#url <- "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/psi25/Homo_sapiens/"

## getting list of files on ftp folder
#filenames <- getURL(url, ftp.use.epsv = FALSE,dirlistonly = TRUE) 
#filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")

filenames_homo <- list.files(paste0(local_url,"H.sapiens/"), full.names = TRUE)

## making empty dataframe to save the total interaction data from each species
Total_list_homosapien <- data.frame(matrix(ncol=7,nrow=0))
colnames(Total_list_homosapien) <- c("id_A", "id_B", "protein_A", "protein_B","protein_fullname_A", "protein_fullname_B" ,"complex")

control_hm <- data.frame(matrix(NA, ncol=2,nrow=length(filenames_homo)))
for(i in 1:length(filenames_homo)){
  print(i)
  file <- filenames_homo[i]
  dir <- "F:/ComplexPortal/csv/H.sapiens/"
  list <- Convert2CSV(file,dir)
  Total_list_homosapien <- rbind(Total_list_homosapien,list)
  control_hm[i,1] <- list$complexID[1]
  control_hm[i,2] <- dim(list)[1]
}


write.csv(unique(unlist(Total_list_homosapien[,c("protein_A","protein_B")])),"protein_list.csv")

uniprot_hm <- read_excel("uniprotID/uniprot_hm.xlsx")

Total_list_homosapien$uniprot_A <- NA
Total_list_homosapien$uniprot_B <- NA

Total_list_homosapien$uniprot_B[grepl("-pro",Total_list_homosapien$protein_B)] <- 
  apply(Total_list_homosapien[grepl("-pro",Total_list_homosapien$protein_B),], 1, function(x) 
    toupper(strsplit(x[["protein_B"]],"-pro")[[1]][1]))

Total_list_homosapien$uniprot_A[grepl("-pro",Total_list_homosapien$protein_A)] <- 
  apply(Total_list_homosapien[grepl("-pro",Total_list_homosapien$protein_A),], 1, function(x) 
    toupper(strsplit(x[["protein_A"]],"-pro")[[1]][1]))

Total_list_homosapien$uniprot_B[!grepl("-pro",Total_list_homosapien$protein_B)] <- uniprot_hm$Entry[match(as.character(Total_list_homosapien$protein_B[!grepl("-pro",Total_list_homosapien$protein_B)]),tolower(uniprot_hm[["yourlist"]]))]  
Total_list_homosapien$uniprot_A[!grepl("-pro",Total_list_homosapien$protein_A)] <- uniprot_hm$Entry[match(as.character(Total_list_homosapien$protein_A[!grepl("-pro",Total_list_homosapien$protein_A)]),tolower(uniprot_hm[["yourlist"]]))]  


Total_list_homosapien$testA <- apply(Total_list_homosapien, 1, function(x) strsplit(x[["protein_A"]],"-")[[1]][1])
Total_list_homosapien$testB <- apply(Total_list_homosapien, 1, function(x) strsplit(x[["protein_B"]],"-")[[1]][1])

Total_list_homosapien$uniprot_A[is.na(Total_list_homosapien$uniprot_A)][grep('^[a-z][0-9][a-z0-9{4}]*$',as.character(Total_list_homosapien$testA[is.na(Total_list_homosapien$uniprot_A)]))] <- 
  Total_list_homosapien$testA[is.na(Total_list_homosapien$uniprot_A)][grep('^[a-z][0-9][a-z0-9{4}]*$',as.character(Total_list_homosapien$testA[is.na(Total_list_homosapien$uniprot_A)]))]

Total_list_homosapien$uniprot_B[is.na(Total_list_homosapien$uniprot_B)][grep('^[a-z][0-9][a-z0-9{4}]*$',as.character(Total_list_homosapien$testB[is.na(Total_list_homosapien$uniprot_B)]))] <- 
  Total_list_homosapien$testB[is.na(Total_list_homosapien$uniprot_B)][grep('^[a-z][0-9][a-z0-9{4}]*$',as.character(Total_list_homosapien$testB[is.na(Total_list_homosapien$uniprot_B)]))]

#View(Total_list_homosapien %>%filter(is.na(uniprot_A) | is.na(uniprot_B)))
#Total_list_homosapien_complete <- Total_list_homosapien%>%filter(!is.na(uniprot_A)&!is.na(uniprot_B))

write.csv(Total_list_homosapien ,"Total_list_homosapien.csv")
#################################################################################
## Saccharaomyces Cerevisiae
## setting the url

url <- "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/psi25/Saccharomyces_cerevisiae/"

filenames_sc <- list.files(paste0(local_url,"S.cerevisiae/"), full.names = TRUE)

## getting list of files on ftp folder
filenames <- getURL(url, ftp.use.epsv = FALSE,dirlistonly = TRUE) 
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")

## making empty dataframe to save the total interaction data from each species
Total_list_saccharomyces_cerevisiae <- data.frame(matrix(ncol=5,nrow=0))
colnames(Total_list_saccharomyces_cerevisiae) <- c("id_A", "id_B", "protein_A", "protein_B", "complex")

for(i in 1:length(filenames_sc)){
  file <- filenames_sc[i]
  dir <- "F:/ComplexPortal/csv/S.cerevisiae/"
  list <- Convert2CSV(file,dir)
  Total_list_saccharomyces_cerevisiae <- rbind(Total_list_saccharomyces_cerevisiae,list)
}


write.csv(unique(unlist(Total_list_saccharomyces_cerevisiae[,c("protein_A","protein_B")])),"protein_list_sc.csv")

uniprot_sc <- read_excel("uniprotID/uniprot_sc.xlsx")

Total_list_saccharomyces_cerevisiae$uniprot_A <- NA
Total_list_saccharomyces_cerevisiae$uniprot_B <- NA

#Total_list_saccharomyces_cerevisiae$uniprot_B[grepl("-pro",Total_list_saccharomyces_cerevisiae$protein_B)] <- 
#  apply(Total_list_saccharomyces_cerevisiae[grepl("-pro",Total_list_saccharomyces_cerevisiae$protein_B),], 1, function(x) 
#    toupper(strsplit(x[["protein_B"]],"-pro")[[1]][1]))

#Total_list_saccharomyces_cerevisiae$uniprot_A[grepl("-pro",Total_list_saccharomyces_cerevisiae$protein_A)] <- 
#  apply(Total_list_saccharomyces_cerevisiae[grepl("-pro",Total_list_saccharomyces_cerevisiae$protein_A),], 1, function(x) 
#    toupper(strsplit(x[["protein_A"]],"-pro")[[1]][1]))

Total_list_saccharomyces_cerevisiae$uniprot_B[!grepl("-pro",Total_list_saccharomyces_cerevisiae$protein_B)] <- 
  uniprot_sc$Entry[match(as.character(Total_list_saccharomyces_cerevisiae$protein_B[!grepl("-pro",Total_list_saccharomyces_cerevisiae$protein_B)]),tolower(uniprot_sc[["yourlist"]]))] 

Total_list_saccharomyces_cerevisiae$uniprot_A[!grepl("-pro",Total_list_saccharomyces_cerevisiae$protein_A)] <- 
  uniprot_sc$Entry[match(as.character(Total_list_saccharomyces_cerevisiae$protein_A[!grepl("-pro",Total_list_saccharomyces_cerevisiae$protein_A)]),tolower(uniprot_sc[["yourlist"]]))]  


#View(Total_list_saccharomyces_cerevisiae %>%filter(is.na(uniprot_A) | is.na(uniprot_B)))
#Total_list_saccharomyces_cerevisiae <- Total_list_saccharomyces_cerevisiae%>%filter(!is.na(uniprot_A)&!is.na(uniprot_B))

write.csv(Total_list_saccharomyces_cerevisiae ,"Total_list_saccharomyces_cerevisiae.csv")

#################################################################################
## Mus_musculus
## setting the url

filenames_mm<- list.files(paste0(local_url,"M.musculus/"), full.names = TRUE)

url <- "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/psi25/Mus_musculus/"

## getting list of files on ftp folder
filenames <- getURL(url, ftp.use.epsv = FALSE,dirlistonly = TRUE) 
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")

## making empty dataframe to save the total interaction data from each species
Total_list_mus_musculus <- data.frame(matrix(ncol=5,nrow=0))
colnames(Total_list_mus_musculus) <- c("id_A", "id_B", "protein_A", "protein_B", "complex")


for(i in 1:length(filenames_mm)){
  file <- filenames_mm[i]
  dir <- "F:/Network_Analysis/ComplexPortal/mus_musculus/"
  list <- Convert2CSV(file,dir)
  Total_list_mus_musculus <- rbind(Total_list_mus_musculus,list)
}


uniprot_mm <- read_excel("uniprotID/uniprot_mm.xlsx")

Total_list_mus_musculus$uniprot_A <- NA
Total_list_mus_musculus$uniprot_B <- NA

Total_list_mus_musculus$uniprot_B[grepl("-pro",Total_list_mus_musculus$protein_B)] <- 
  apply(Total_list_mus_musculus[grepl("-pro",Total_list_mus_musculus$protein_B),], 1, function(x) 
    toupper(strsplit(x[["protein_B"]],"-pro")[[1]][1]))

Total_list_mus_musculus$uniprot_A[grepl("-pro",Total_list_mus_musculus$protein_A)] <- 
  apply(Total_list_mus_musculus[grepl("-pro",Total_list_mus_musculus$protein_A),], 1, function(x) 
    toupper(strsplit(x[["protein_A"]],"-pro")[[1]][1]))

Total_list_mus_musculus$uniprot_B[!grepl("-pro",Total_list_mus_musculus$protein_B)] <- uniprot_mm$Entry[match(as.character(Total_list_mus_musculus$protein_B[!grepl("-pro",Total_list_mus_musculus$protein_B)]),tolower(uniprot_mm[["yourlist"]]))]  
Total_list_mus_musculus$uniprot_A[!grepl("-pro",Total_list_mus_musculus$protein_A)] <- uniprot_mm$Entry[match(as.character(Total_list_mus_musculus$protein_A[!grepl("-pro",Total_list_mus_musculus$protein_A)]),tolower(uniprot_mm[["yourlist"]]))]  


Total_list_mus_musculus$testA <- apply(Total_list_mus_musculus, 1, function(x) strsplit(x[["protein_A"]],"-")[[1]][1])
Total_list_mus_musculus$testB <- apply(Total_list_mus_musculus, 1, function(x) strsplit(x[["protein_B"]],"-")[[1]][1])

Total_list_mus_musculus$uniprot_A[is.na(Total_list_mus_musculus$uniprot_A)][grep('^[a-z][0-9][a-z0-9{4}]*$',as.character(Total_list_mus_musculus$testA[is.na(Total_list_mus_musculus$uniprot_A)]))] <- 
  Total_list_mus_musculus$testA[is.na(Total_list_mus_musculus$uniprot_A)][grep('^[a-z][0-9][a-z0-9{4}]*$',as.character(Total_list_mus_musculus$testA[is.na(Total_list_mus_musculus$uniprot_A)]))]

Total_list_mus_musculus$uniprot_B[is.na(Total_list_mus_musculus$uniprot_B)][grep('^[a-z][0-9][a-z0-9{4}]*$',as.character(Total_list_mus_musculus$testB[is.na(Total_list_mus_musculus$uniprot_B)]))] <- 
  Total_list_mus_musculus$testB[is.na(Total_list_mus_musculus$uniprot_B)][grep('^[a-z][0-9][a-z0-9{4}]*$',as.character(Total_list_mus_musculus$testB[is.na(Total_list_mus_musculus$uniprot_B)]))]

#Total_list_mus_musculus <- Total_list_mus_musculus%>%filter(!is.na(uniprot_A)&!is.na(uniprot_B))

write.csv(Total_list_mus_musculus ,"Total_list_mus_musculus.csv")

#################################################################################
## E.coli
## setting the url

filenames_ec<- list.files(paste0(local_url,"E.coli/"), full.names = TRUE)

url <- "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/psi25/Mus_musculus/"

## getting list of files on ftp folder
filenames <- getURL(url, ftp.use.epsv = FALSE,dirlistonly = TRUE) 
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")

## making empty dataframe to save the total interaction data from each species
Total_list_e_coli <- data.frame(matrix(ncol=5,nrow=0))
colnames(Total_list_e_coli) <- c("id_A", "id_B", "protein_A", "protein_B", "complex")


for(i in 1:length(filenames_ec)){
  file <- filenames_ec[i]
  dir <- "F:/ComplexPortal/csv/E.coli/"
  list <- Convert2CSV(file,dir)
  Total_list_e_coli <- rbind(Total_list_e_coli,list)
}

write.csv(unique(unlist(Total_list_e_coli[,c("protein_A","protein_B")])),"protein_list_ec.csv")

uniprot_ec <- read_excel("uniprotID/uniprot_ec.xlsx")

Total_list_e_coli$uniprot_A <- NA
Total_list_e_coli$uniprot_B <- NA


Total_list_e_coli$uniprot_B[!grepl("-pro",Total_list_e_coli$protein_B)] <- uniprot_ec$Entry[match(as.character(Total_list_e_coli$protein_B[!grepl("-pro",Total_list_e_coli$protein_B)]),tolower(uniprot_ec[["yourlist"]]))]  
Total_list_e_coli$uniprot_A[!grepl("-pro",Total_list_e_coli$protein_A)] <- uniprot_ec$Entry[match(as.character(Total_list_e_coli$protein_A[!grepl("-pro",Total_list_e_coli$protein_A)]),tolower(uniprot_ec[["yourlist"]]))]  


Total_list_e_coli$testA <- apply(Total_list_e_coli, 1, function(x) strsplit(x[["protein_A"]],"-")[[1]][1])
Total_list_e_coli$testB <- apply(Total_list_e_coli, 1, function(x) strsplit(x[["protein_B"]],"-")[[1]][1])


write.csv(Total_list_e_coli ,"Total_list_e_coli.csv")

