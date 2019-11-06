rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

# install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","magrittr","igraph","sqldf","diffusr", "dnet", "plotly")
bioconductor.packages <- c("biomaRt")
source("../Functions/Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

X9606_WHOLE_ORGANISM_integrated <- read_delim("9606-WHOLE_ORGANISM-integrated.txt", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)


getGeneName <- function(uniprotID) {
  mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  return(biomaRt::getBM(attributes=c('uniprot_gn', 'external_gene_name'), 
                        filters = 'uniprot_gn', 
                        values =as.character(uniprotID), 
                        mart = mart))
}



HUMAN_9606_idmapping <- read_delim("F:/Network_Analysis/Mapping/HUMAN_9606_idmapping.dat", header=FALSE)
X9606_WHOLE_ORGANISM_integrated$uniprotID <- HUMAN_9606_idmapping$V1[match(X9606_WHOLE_ORGANISM_integrated$string_external_id, HUMAN_9606_idmapping$V3)]

gene_list <- getGeneName(X9606_WHOLE_ORGANISM_integrated$uniprotID)
X9606_WHOLE_ORGANISM_integrated$geneName <- gene_list$external_gene_name[match(X9606_WHOLE_ORGANISM_integrated$uniprotID,gene_list$uniprot_gn)]

X9606_WHOLE_ORGANISM_integrated$checked <- ifelse(X9606_WHOLE_ORGANISM_integrated$uniprotID %in% colnames(tuneDF),1,0)

dt <- X9606_WHOLE_ORGANISM_integrated%>% filter(checked==1)
dt2 <- dt%>%group_by(geneName)%>% summarise(abundance_total =sum(abundance))


tuneDF_list <- getGeneName(colnames(tuneDF))
tuneDF_list$abundance <- X9606_WHOLE_ORGANISM_integrated$abundance[match(as.character(tuneDF_list$uniprot_gn),as.character(X9606_WHOLE_ORGANISM_integrated$uniprotID))]

tuneDF_list2 <-merge(x=tuneDF_list,y=dt2,by.x="external_gene_name",by.y = "geneName",all.x=TRUE)

mind <- min(tuneDF_list2$abundance,na.rm=T)
sdd <- sd(tuneDF_list2$abundance,na.rm=T)

min(tuneDF_list2$abundance[!tuneDF_list2$abundance>mind+2*sdd&tuneDF_list2$abundance!=0],na.rm=T)

tuneDF_list2$abundance[is.na(tuneDF_list2$abundance)] <- 0.001
tuneDF_list2$ratio <- tuneDF_list2$abundance/tuneDF_list2$abundance_total

##### Import Gene2protein mapping
load("F:/Network_Analysis/Mapping/graph/Gene2Protein.RData")
M_Gene2Protein <- M_Gene2Protein[,colnames(M_Gene2Protein)!=""]
dim(M_Gene2Protein)
M_Gene2Protein <- M_Gene2Protein[sort(rownames(M_Gene2Protein)),sort(colnames(M_Gene2Protein))]
