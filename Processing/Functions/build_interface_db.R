# Mehrnoosh Oghbaie
# 10/22/2019
# Build db used in interface

########################################################################################

InterfaceDB <- R6Class("InterfaceDB",  list(
  tuneBackup = "backup/tuneBackup.RData", # tune backup address
  idmap = "../data/mapping/db/geneMania/idmapping.txt",
  ChEBI = "../data/complexportal/ChEBI_Results.tsv",
  tuning = "backup/tuning_Layers.RData",
  edges_link = "..//data/complexportal/Total_list_homosapien.csv",
  db_hm = NA,
  Total_edges = NA, 
  idmapping = NA,
  protein_list = NA,
  ChEBI_Results = NA,
  M_Disease2Protein = NA,
  mim_morbid = NA,
  gM1 = NA,
  uniprot_ubiquitin_H_sapiens = NA,
  A = NA, 
  p0 = NA,
  complex_protein = NA,
  M_Pathway.normalized = NA,
  M_xl_graph_corrected.normalized = NA,
  PPI_layer = NA,
  # Load gene/protein backup
  
  loadDB = function(){
    file_names = c(self$tuneBackup,
                   self$tuning)
    lapply(file_names,load,.GlobalEnv)
    
    ##  idmapping
    self$idmapping <- read.delim("E:/ProteinComplexEntrant/data/mapping/db/geneMania/identifier_mappings.txt")
    
    ##  db_hm
    db <- tune$complexPortalDB
    db[["cc"]] <- apply(db,1, function(x) ifelse(x[["Go.CC"]] =="",0,str_count(x[["Go.CC"]], "\\|")+1))
    self$db_hm <- db
    
    ## Total edges
    Total_edges <- Total_list_homosapien <- read.csv(self$edges_link, row.names=1, stringsAsFactors=FALSE)
    self$Total_edges <- Total_edges
  
    ##  protein_list
    ensembl=useMart("ensembl")
    mart = useDataset("hsapiens_gene_ensembl",mart=ensembl)
    View(listAttributes(mart))
    protein_gene <- biomaRt::getBM(attributes=c('external_gene_name', 'uniprotswissprot'), 
                   filters = 'uniprotswissprot', 
                   values =tune$proteinList, 
                   mart = mart)
    self$protein_list <- protein_gene
    
    ## ChEBI_Results
    ChEBI_Results <- read_delim(self$ChEBI, 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
    self$ChEBI_Results <- ChEBI_Results
    
    ## M_Disease2Protein
    load(tune$disease2ProteinBipartiteMapping)
    self$M_Disease2Protein <- M_Disease2Protein
    
    ##  mim_morbid
  
    mim_morbid <- biomaRt::getBM(attributes=c('mim_morbid_accession', 'mim_morbid_description'), 
                                 filters = 'mim_morbid_accession', 
                                 values =colnames(self$M_Disease2Protein), 
                                 mart = mart)
    self$mim_morbid <- mim_morbid
    
    ##  gM1
    load(tune$go2ProteinBipartiteMapping)
    g <- graph_from_data_frame(as_edgelist(g_GO2Protein, names = TRUE)[E(g_GO2Protein)$category=="C",], directed=FALSE)
    go.cc <- as_edgelist(g_GO2Protein, names = TRUE)[E(g_GO2Protein)$category=="C",][,1]
    gM <- as_adjacency_matrix(g)
    self$gM1 <- gM[rownames(gM) %in% unique(go.cc), !colnames(gM) %in% unique( go.cc)]
    
    ##  uniprot_ubiquitin_H_sapiens
    uniprot_ubiquitin_H_sapiens = read_excel(tune$ubiquitinProteinHomo)
    self$uniprot_ubiquitin_H_sapiens = uniprot_ubiquitin_H_sapiens %>% dplyr::filter(grepl("ubiquitin|Ubiquitin",`Protein names`))
   
    ##  A 
    A <- Matrix(0, ncol = dim(tuning$Pathway)[2], nrow = dim(tuning$Pathway)[1])
    for(name in names(tuning)){
      A <- A+ tune$tuneList$Average_parameters[name,"params_avg"]*tuning[[name]]
    }
    self$A <- A
    
    ##  p0
    p0 <- Matrix(0, nrow=dim(self$A)[1], ncol=1)
    rownames(p0) <- rownames(A)
    self$p0 <- p0
    
    load(file = self$tuning)
    
    ##  M_Pathway.normalized
    self$M_Pathway.normalized <- tuning$Pathway
    
    ##  M_xl_graph_corrected.normalized
    self$M_xl_graph_corrected.normalized <- tuning$CrossLink
    
    ## PPI_layer
    #  M_gene_PPI_graph.normalized2
    #  M_protein_PPI_graph.normalized
    PPI_layer <- tune$tuneList$Average_parameters["Physical_Protein","params_avg"]* tuning$Physical_Protein +
      tune$tuneList$Average_parameters["Physical_G","params_avg"]* tuning$Physical_G
    self$PPI_layer <- PPI_layer
    
    ## Complex_protein 
    Complex_protein2 <- data.frame(matrix(0, ncol=6, nrow=0))
    colnames(Complex_protein2) <- c("Complex_ID", "ComplexName", "uniprotID", "count","Identifiers", "geneName") 
    for(i in 1:dim(self$db_hm)[1]){
      dbf <- self$db_hm[i,]
      uniprot <- gsub( "\\([0-9]\\)","",strsplit(dbf$Identifiers,"\\|")[[1]])
      
      geneUniprot <- biomaRt::getBM(attributes=c('external_gene_name', 'uniprotswissprot'), 
                                    filters = 'uniprotswissprot', 
                                    values =uniprot, 
                                    mart = mart)
      dt <- data.frame(cbind(rep(dbf$X.Complex.ac, dbf$count),
                             rep(dbf$Recommended.name, dbf$count),
                             gsub( "\\([0-9]\\)","",strsplit(dbf$Identifiers,"\\|")[[1]]),
                             rep(dbf$count, dbf$count),
                             paste0(gsub( "\\([0-9]\\)","",strsplit(dbf$Identifiers,"\\|")[[1]]),collapse="|"),
                             geneUniprot$external_gene_name[match(uniprot, geneUniprot$uniprotswissprot)]
      ))
      colnames(dt) <- c("Complex_ID", "ComplexName", "uniprotID", "count","Identifiers", "geneName")
      Complex_protein2 <- rbind(Complex_protein2, dt)
    }
    
    self$complex_protein <- Complex_protein2
    
    invisible(self)
  },
  
  directInteraction = function(){
    ## Downloading direct interaction from PSICQUIC
    psicquic = PSICQUIC()
    DB = providers(psicquic)
    tbl_direct = PSICQUIC::interactions(psicquic, species="9606",type="direct interaction", provider=DB,speciesExclusive= TRUE)
    #tbl_direct = addGeneInfo(IDMapper("9606"),tbl_direct) 
    tbl_direct = createEmptyCol(tbl_direct)
    tbl_direct2 = data.frame(t(apply(tbl_direct, 1, fillSource)))
    tbl_direct2[["detectionMethod"]] = apply(tbl_direct2,1,fillMethod)
    tbl_direct2$confidence = NA
    tbl_direct3 = data.frame(t(apply(tbl_direct2, 1, fillconfidenceScore)))
    tbl_direct3$pubmed = as.integer(as.character(tbl_direct3$pubmed))
    tbl_direct4 <- tbl_direct3
    tbl_direct4[is.na(tbl_direct4$year),"year"] <- apply(tbl_direct4,1, function(x) 
      as.integer(gsub("[\\(\\)]", "", regmatches(as.character(x[["firstAuthor"]]), gregexpr("\\(.*?\\)", as.character(x[["firstAuthor"]])))[[1]])))
    tbl_direct4 <- data.frame(tbl_direct4)
      
    tbl_direct4[["A.name"]] <- apply(tbl_direct4 ,1, getGeneName, col="altA")
    tbl_direct4[["B.name"]] <- apply(tbl_direct4 ,1, getGeneName, col="altB")
    
    ### Save the final table in a Rdata file
    
    save(tbl_direct4 , file="../interface/db/directinteraction.RData")
    invisible(self)
  }
))

for(i in 1:360){
  print(i)
  m <- 1+(i-1)*1000
  n <- 1000*i
  print(m)
  print(n)
  tbl_direct4[m:n,"A.name"] <- apply(tbl_direct4[m:n,] ,1, getGeneName, col="altA")
  tbl_direct4[m:n,"B.name"] <- apply(tbl_direct4[m:n,] ,1, getGeneName, col="altB")
}


