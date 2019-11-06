# Mehrnoosh Oghbaie
# 08/30/2019
# Tuning parameters


Tune <- R6Class("Tune",  list(
  geneManiaInput = "../data/input_geneMania/networks.txt", # List of publication from GeneMania
  geneManiaNetworks = NA,
  pubmed.complex = NA,
  complexPortalHS = "../data/complexportal/records/H.sapiens.tsv", # ComplexPortal database directory
  interactionComplexPortalHomo = "../data/complexportal/Total_list_homosapien.csv",
  ubiquitinProteinHomo = "../data/ubiquitin/uniprot_ubiquitin_HS.xlsx",
  geneBackup = "../data/layers/Gene_Backup.RData", # Final gene_level graphs
  proteinBackup = "../data/layers/Protein_Backup.RData", # Final protein_level graphs
  normalizedLayersBackup = "../data/layers/Normalized_Layers_Backup.RData",
  diseaseSimilarityGraph = "../data/layers/disease/graph/disease_similarity.RData",
  protein2GeneBipartiteMapping = "../data/mapping/graph/protein2Gene.RData",
  disease2ProteinBipartiteMapping = "../data/mapping/graph/Disease2Protein.RData",
  go2ProteinBipartiteMapping = "../data/mapping/graph/GO2Protein.RData",
  # Merged listof genes, proteins, diseases and go
  geneList = NA,
  proteinList = NA,
  diseaseList = NA,
  goList = NA,
  tuneList = NA,
  ubiquitinList = NA,
  complexPortalDB = NA,
  interactionHSDB = NA,
  Complex_protein_edege = NA ,
  mixedNetworks = list(),
  tuning = list(),
  totalResult = NA,
  ppiNet = NA,
  # Load gene/protein backup
  loadBackupGeneProtein = function(){
    print("loadBackupGeneProtein")
    file_names = c(self$geneBackup,
                   self$proteinBackup,
                   self$diseaseSimilarityGraph,
                   self$protein2GeneBipartiteMapping,
                   self$disease2ProteinBipartiteMapping,
                   self$go2ProteinBipartiteMapping)
    lapply(file_names,load,.GlobalEnv)
    self$proteinList = protein_list2
    self$geneList = gene_list
    invisible(self)
  },
  mapGeneDisease2Protein = function(){
    print("mapGeneDisease2Protein")
    #self$loadBackupGeneProtein();
    adj1 <- get.adjacency(protein2Gene, attr="Weight")
    M_protein2Gene <- adj1[rownames(adj1) %in% protein_list2[order(protein_list2)], colnames(adj1) %in% gene_list[order(gene_list)]]
    M_protein2Gene <- M_protein2Gene[,order(colnames(M_protein2Gene))]
    
    ## Map gene-level network using mapping matrix and normalize them
    M_gene_PPI_graph <- get.adjacency(gene_PPI_graph, attr="Weight")
    M_gene_PPI_graph <- M_gene_PPI_graph[order(rownames(M_gene_PPI_graph)), order(colnames(M_gene_PPI_graph))]
    M_gene_PPI_graph.normalized <- normalizeConvert(t(M_protein2Gene),M_gene_PPI_graph)
    extra_proteins <- protein_list2[!protein_list2 %in% colnames(M_gene_PPI_graph.normalized)]
    gene_PPI_graph2 <- graph_from_adjacency_matrix(M_gene_PPI_graph.normalized, weighted=TRUE)
    gene_PPI_graph2 <- add.vertices(gene_PPI_graph2, nv=length(extra_proteins),attr = list(name=extra_proteins))
    M_gene_PPI_graph.normalized2 <- get.adjacency(gene_PPI_graph2, attr="weight")
    M_gene_PPI_graph.normalized2 <- M_gene_PPI_graph.normalized2[order(rownames(M_gene_PPI_graph.normalized2)), order(colnames(M_gene_PPI_graph.normalized2))]
    
    M_coexpression_graph <- get.adjacency(coexpression_graph, attr="Weight")
    M_coexpression_graph <- M_coexpression_graph[order(rownames(M_coexpression_graph)), order(colnames(M_coexpression_graph))]
    M_coexpression_graph.normalized <- normalizeConvert(t(M_protein2Gene),M_coexpression_graph)
    coexpression_graph2 <- graph_from_adjacency_matrix(M_coexpression_graph.normalized, weighted=TRUE)
    coexpression_graph2 <- add.vertices(coexpression_graph2, nv=length(extra_proteins),attr = list(name=extra_proteins))
    M_coexpression_graph.normalized2 <- get.adjacency(coexpression_graph2, attr="weight")
    M_coexpression_graph.normalized2 <- M_coexpression_graph.normalized2[order(rownames(M_coexpression_graph.normalized2)), order(colnames(M_coexpression_graph.normalized2))]
    
    M_genetic_graph <- get.adjacency(genetic_graph, attr="Weight")
    M_genetic_graph <- M_genetic_graph[order(rownames(M_genetic_graph)), order(colnames(M_genetic_graph))]
    M_genetic_graph.normalized <- normalizeConvert(t(M_protein2Gene),M_genetic_graph)
    genetic_graph2 <- graph_from_adjacency_matrix(M_genetic_graph.normalized, weighted=TRUE)
    genetic_graph2 <- add.vertices(genetic_graph2, nv=length(extra_proteins),attr = list(name=extra_proteins))
    M_genetic_graph.normalized2 <- get.adjacency(genetic_graph2, attr="weight")
    M_genetic_graph.normalized2 <- M_genetic_graph.normalized2[order(rownames(M_genetic_graph.normalized2)), order(colnames(M_genetic_graph.normalized2))]
    
    M_prediction_graph <- get.adjacency(prediction_graph, attr="Weight")
    M_prediction_graph <- M_prediction_graph[order(rownames(M_prediction_graph)), order(colnames(M_prediction_graph))]
    M_prediction_graph.normalized <- normalizeConvert(t(M_protein2Gene),M_prediction_graph)
    prediction_graph2 <- graph_from_adjacency_matrix(M_prediction_graph.normalized, weighted=TRUE)
    prediction_graph2 <- add.vertices(prediction_graph2, nv=length(extra_proteins),attr = list(name=extra_proteins))
    M_prediction_graph.normalized2 <- get.adjacency(prediction_graph2, attr="weight")
    M_prediction_graph.normalized2 <- M_prediction_graph.normalized2[order(rownames(M_prediction_graph.normalized2)), order(colnames(M_prediction_graph.normalized2))]
    
    M_colocalization_graph <- get.adjacency(colocalization_graph, attr="Weight")
    M_colocalization_graph <- M_colocalization_graph[order(rownames(M_colocalization_graph)), order(colnames(M_colocalization_graph))]
    M_colocalization_graph.normalized <- normalizeConvert(t(M_protein2Gene),M_colocalization_graph)
    colocalization_graph2 <- graph_from_adjacency_matrix(M_colocalization_graph.normalized, weighted=TRUE)
    colocalization_graph2 <- add.vertices(colocalization_graph2, nv=length(extra_proteins),attr = list(name=extra_proteins))
    M_colocalization_graph.normalized2 <- get.adjacency(colocalization_graph2, attr="weight")
    M_colocalization_graph.normalized2 <- M_colocalization_graph.normalized2[order(rownames(M_colocalization_graph.normalized2)), order(colnames(M_colocalization_graph.normalized2))]
    
    ## Normalize protein level networks
    M_Pathway <- get.adjacency(Pathway)
    M_Pathway <- M_Pathway[order(rownames(M_Pathway)), order(colnames(M_Pathway))]
    M_Pathway.normalized <- normalizeM(M_Pathway)
    
    M_protein_PPI_graph <- get.adjacency(protein_PPI_graph, attr="Weight")
    M_protein_PPI_graph <- M_protein_PPI_graph[order(rownames(M_protein_PPI_graph)), order(colnames(M_protein_PPI_graph))]
    M_protein_PPI_graph.normalized <- normalizeM(M_protein_PPI_graph)
    
    M_xl_graph <- get.adjacency(xl_graph)
    M_xl_graph <- M_xl_graph[rownames(M_xl_graph) %in% rownames(M_protein_PPI_graph),colnames(M_xl_graph) %in% colnames(M_protein_PPI_graph)]
    M_xl_graph.normalized <- normalizeM(M_xl_graph)
  
    # Normalize disease similarity network
    ## Importing Protein to Disease mapping
    M_Disease2Protein <- M_Disease2Protein[rownames(M_Disease2Protein)!="",]
    M_Disease2Protein <- M_Disease2Protein[sort(rownames(M_Disease2Protein)),sort(colnames(M_Disease2Protein))]

   # Disease <- graph.data.frame(DiseaseSimilarity, directed = FALSE)
    adj2 <- get.adjacency(disease_similarity_graph,sparse = TRUE) 
    AdjacencyMatrix.Disease <- adj2[rownames(adj2) %in% as.character(uniqueDisease$unique.unlist.DiseaseSimilarity..), colnames(adj2) %in% as.character(uniqueDisease$unique.unlist.DiseaseSimilarity..)]
    AdjacencyMatrix.Disease<- AdjacencyMatrix.Disease[order(rownames(AdjacencyMatrix.Disease)), order(colnames(AdjacencyMatrix.Disease))]
    M_Disease.similarity.normalized <- normalizeConvert(t(M_Disease2Protein),AdjacencyMatrix.Disease)
    
    Complex_protein_edege <- Matrix(0, nrow=dim(M_Pathway.normalized)[1], ncol=dim(M_Pathway.normalized)[2])
    colnames(Complex_protein_edege) <- colnames(M_Pathway.normalized)
    rownames(Complex_protein_edege) <- rownames(M_Pathway.normalized)
    
    Total_list_homosapien <- read_csv(self$interactionComplexPortalHomo)
    self$interactionHSDB <- Total_list_homosapien
    df <- unique(Total_list_homosapien[,c("uniprot_A", "uniprot_B")]) %>% filter(uniprot_A %in% colnames(Complex_protein_edege)& uniprot_B %in% colnames(Complex_protein_edege))
    
    for(i in 1:dim(df)[1]){
      Complex_protein_edege[df[["uniprot_A"]][i],df[["uniprot_B"]][i]]<- 1
      Complex_protein_edege[df[["uniprot_B"]][i],df[["uniprot_A"]][i]] <- 1
    }
    self$Complex_protein_edege <- Complex_protein_edege

    tuning <- list()
    tuning[["Pathway"]] <- M_Pathway.normalized
    tuning[["Physical_Protein"]] <- M_protein_PPI_graph.normalized 
    tuning[["Coexpression_G"]] <- M_coexpression_graph.normalized2
    tuning[["Colocalization_G"]] <- M_colocalization_graph.normalized2
    tuning[["Genetic_G"]] <- M_genetic_graph.normalized2
    tuning[["Physical_G"]] <- M_gene_PPI_graph.normalized2
    tuning[["Prediction_G"]] <-  M_prediction_graph.normalized2
    tuning[["CrossLink"]] <- M_xl_graph.normalized
    tuning[["DiseaseSimilarity"]] <- M_Disease.similarity.normalized
    
    
    save(M_gene_PPI_graph.normalized2,
         M_coexpression_graph.normalized2,
         M_genetic_graph.normalized2,
         M_prediction_graph.normalized2,
         M_colocalization_graph.normalized2,
         M_Pathway.normalized,
         M_protein_PPI_graph.normalized,
         M_xl_graph.normalized,
         M_Disease.similarity.normalized, 
         Complex_protein_edege,
         file = self$normalizedLayersBackup)
    
    save(tuning ,file="backup/tuning_Layers.RData")
    invisible(self)
  },
  makeComplexPortalDB = function(){
    print("makeComplexPortalDB")
    ### Finding out the cell component from complexportal
    H.sapiens <- read.delim(self$complexPortalHS)
    db <- H.sapiens %>% dplyr::select(X.Complex.ac,
                                      Recommended.name , 
                                      Description ,
                                      Identifiers..and.stoichiometry..of.molecules.in.complex,
                                      Complex.assembly,
                                      Go.Annotations,
                                      Disease) %>% dplyr::rename(Identifiers= Identifiers..and.stoichiometry..of.molecules.in.complex )%>%
      mutate(Complex.assembly =toupper(gsub("\\.|-","",Complex.assembly)),
             X.Complex.ac = as.character(X.Complex.ac),
             Identifiers = as.character(Identifiers))
    
    for( i in db$X.Complex.ac[grepl("CPX",db$Identifiers)]){
      j = gsub("\\s*\\([^\\)]+\\)","",strsplit(db$Identifiers[db$X.Complex.ac==i] ,"\\|")[[1]][1])
      db$Identifiers[db$X.Complex.ac==i] <- paste0(c(db$Identifiers[db$X.Complex.ac==j],strsplit(db$Identifiers[db$X.Complex.ac==i] ,"\\|")[[1]][2]), collapse="|")
    }
    
    cc = unique(as_edgelist(g_GO2Protein, names = TRUE)[,1][E(g_GO2Protein)$category=="C"])
    db$Go.CC <- NA
    for(i in 1:dim(db)[1]){
      db$Go.CC[i] <- paste0(strsplit(as.character(db[["Go.Annotations"]][i]) ,"\\|")[[1]][grepl(paste0(cc, collapse="|"),                                                                                           strsplit(as.character(db[["Go.Annotations"]][i]) ,"\\|")[[1]])],collapse="|")
    }
    
    db$count <- NA
    db$count <- str_count(as.character(db$Identifiers), "\\|") +1
    
    ### Get ubiquitin list
    uniprot_ubiquitin_H_sapiens = read_excel(self$ubiquitinProteinHomo)
    uniprot_ubiquitin_H_sapiens = uniprot_ubiquitin_H_sapiens %>% dplyr::filter(grepl("ubiquitin|Ubiquitin",`Protein names`))
    
    dt <- data.frame(uniprotID=protein_list2)
    dt$ubiquitin <- ifelse(as.character(dt$uniprotID) %in% uniprot_ubiquitin_H_sapiens$Entry,1,0)
    ubiquitin_list <- as.character(dt$uniprotID[dt$ubiquitin==1])
    
    ##### make matrix of go to proteins
    g <- graph_from_data_frame(as_edgelist(g_GO2Protein, names = TRUE)[E(g_GO2Protein)$category=="C",], directed=FALSE)
    go.cc <- as_edgelist(g_GO2Protein, names = TRUE)[E(g_GO2Protein)$category=="C",][,1]
    gM <- as_adjacency_matrix(g)
    gM1 <- gM[rownames(gM) %in% unique(go.cc), !colnames(gM) %in% unique( go.cc)]
    
    db$checksum <- NA
    db$component_count <- NA
    db$filter <- NA
    missing <- c()
    for(i in 1:dim(db)[1]){
      go_cc <- gsub("\\s*\\([^\\)]+\\)","",strsplit(db$Go.CC[i] ,"\\|")[[1]])
      Identifiers <- gsub("\\s*\\([^\\)]+\\)","",strsplit(as.character(db$Identifiers[i]) ,"\\|")[[1]])
      Identifiers <-unlist(lapply(Identifiers,function(x) strsplit(x,"-")[[1]][1]))
      Identifiers <- Identifiers[!grepl("CHEBI:|_9606",Identifiers)]
      protein_interest <- NA
      if(!is.null(dim(gM1[rownames(gM1) %in% go_cc,]))){
        colsum <- Matrix::colSums(gM1[rownames(gM1) %in% go_cc,])
        protein_interest <- names(colsum[colsum!=0])
        db$component_count[i] <- length(unique(c(protein_interest,ubiquitin_list)))
        # Check if all of the complex members exists in the define cell component
        db$filter[i] <- ifelse(length(Identifiers[Identifiers %in% c(protein_interest,ubiquitin_list)])>=length(Identifiers),1,0)
        # Check if all of the complex member exists in our networks
        db$checksum[i] <- ifelse(length(Identifiers[Identifiers %in% protein_list2])>=length(Identifiers),1,0)
        missing <- c(missing,Identifiers[!Identifiers %in% protein_list2])
      }
      if(length(gM1[rownames(gM1) %in% go_cc,])>0& is.null(dim(gM1[rownames(gM1) %in% go_cc,])))
      {
        gM2 <- gM1[rownames(gM1) %in% go_cc,]
        protein_interest <- names(gM2==1)
        db$component_count[i] <- length(unique(c(protein_interest,ubiquitin_list)))
        db$filter[i] <- ifelse(length(Identifiers[Identifiers %in% c(protein_interest,ubiquitin_list)])>=length(Identifiers),1,0)
        db$checksum[i] <- ifelse(length(Identifiers[Identifiers %in% protein_list2])>=length(Identifiers),1,0)
        missing <- c(missing,Identifiers[!Identifiers %in% protein_list2])
      }
    }
    Total_list_homosapien <- self$interactionHSDB
    dk <- Total_list_homosapien[,c("protein_A","protein_B","uniprot_A", "uniprot_B","complexID")] %>% filter(uniprot_A %in% colnames(self$Complex_protein_edege)& uniprot_B %in% colnames(self$Complex_protein_edege))
    
    db$number.links <- 0
    for(i in 1:dim(db)[1]){
      db$number.links[i] <-  dim(dk[dk$complexID == db[["X.Complex.ac"]][i],])[1]
    }
    self$complexPortalDB = db
    invisible(self)
    
  }, 
  tuneLayers = function(){
   print("tuneLayers")
    load(file="backup/tuning_Layers.RData")
    #self$makeComplexPortalDB();
    db <- self$complexPortalDB
    db_more2link <- db %>% dplyr::filter(number.links>2)

    folds <- rep_len(1:10, length(db_more2link$X.Complex.ac))
    folds <- sample(folds, length(db_more2link$X.Complex.ac))
    flds = lapply(1:10,function(x) db_more2link$X.Complex.ac[folds==x])
    Total_list_homosapien <- self$interactionHSDB
    dk <- Total_list_homosapien[,c("protein_A","protein_B","uniprot_A", "uniprot_B","complexID")] %>% filter(uniprot_A %in% colnames(self$Complex_protein_edege)& uniprot_B %in% colnames(self$Complex_protein_edege))
    
    tuneList = list()
    for(j in 1:10){
      dk2 = dk %>% dplyr::filter(!complexID %in% db_more2link$X.Complex.ac[folds==j]) %>%dplyr::select(uniprot_A, uniprot_B)
      tuneList[[paste0("Missing run ",j)]] = paste0(db_more2link$X.Complex.ac[folds==j], collapse =",") 
      M_Pathway.normalized <- tuning[["Pathway"]]
      Complex_protein_edege2 <- Matrix(0, nrow=dim(M_Pathway.normalized)[1], ncol=dim(M_Pathway.normalized)[2])
      colnames(Complex_protein_edege2) <- colnames(M_Pathway.normalized)
      rownames(Complex_protein_edege2) <- rownames(M_Pathway.normalized)
      for(i in 1:dim(dk2)[1]){
        Complex_protein_edege2[dk2[["uniprot_A"]][i],dk2[["uniprot_B"]][i]]<- 1
        Complex_protein_edege2[dk2[["uniprot_B"]][i],dk2[["uniprot_A"]][i]] <- 1
      }
      df <- data.frame(matrix(NA, nrow=length(names(tuning)), ncol=6))
      colnames(df) <- c("layer", "S", "denom","nom", "prim","prim1")
      
      for(i in 1: 9){
        name <- names(tuning)[i]
        X <- tuning[[name]]
        S <- sum(diag(t(X)*X))
        link <- Complex_protein_edege2
        om <- X
        size = dim(om)[1]*dim(om)[1]
        dim(om) <-  c(size,1)
        dim(link) <- c(size,1)
        nom <- t(om)%*%link
        denom <- crossprod(om)
        prim <- nom[1,1]/denom[1,1]
        prim1 <- (nom+(1/9)*S)/(denom+S)
        df[i,"layer"] <- name
        df[["S"]][df$layer==name ] <- S
        df[i, "denom"] <- denom[1,1]
        df[i, "nom"] <- nom[1,1]
        df[i, "prim"] <- prim
        df[i, "prim1"] <- prim1
      }
      
      df$params <- df$prim/sum(df$prim, na.rm=T)
      df$params1 <- df$prim1/sum(df$prim1, na.rm=T)
      tuneList[[paste0("Parameters ",j)]]<- df
      
    }
    self$tuneList = tuneList
    invisible(self)
  }, 
  takeAverageParameters = function(){
    print("takeAverageParameters")
    d <- lapply(1:10,function(x) self$tuneList[[names(self$tuneList)[grepl("Parameters",names(self$tuneList))][x]]][,"params"])
    colMeans(do.call(rbind,d))
    
    d1 <- lapply(1:10,function(x) self$tuneList[[names(self$tuneList)[grepl("Parameters",names(self$tuneList))][x]]][,"params1"])
    colMeans(do.call(rbind,d1))
    
    dTune <- data.frame(cbind(
      colMeans(do.call(rbind,d)),
      100*apply(do.call(rbind,d), 2, sd)/colMeans(do.call(rbind,d)),
      colMeans(do.call(rbind,d1)),
      100*apply(do.call(rbind,d1), 2, sd)/colMeans(do.call(rbind,d1))
    ))
    colnames(dTune) <- c("params_avg", "coeff_var%", "params1_avg", "coeff_var1%")
    rownames(dTune) <- self$tuneList$`Parameters 1`$layer
    self$tuneList[["Average_parameters"]] <- dTune
    
    invisible(self)
  },
  buildCombinationLayers = function(){
    print("buildCombinationLayers")
    load(file="backup/tuning_Layers.RData")
  order=list(c(2),
             c(2,1),
             c(2,1,8),
             c(2,1,8,9),
             c(2,1,8,6),
             c(2,1,8,6,3),
             c(2,1,8,6,3,4),
             c(2,1,8,6,3:5),
             c(2,1,8,6,3:5,7),
             c(2,1,8,3:7,9))
  mixedNetworks <- list()
  for(i in 1:length(order)){
    nets <- tuning[names(tuning)[order[i][[1]]]]
    params1 <- self$tuneList$Average_parameters[["params1_avg"]][order[i][[1]]]
    params1 <- params1/sum(params1)
    n <- dim(tuning[[names(tuning)[1]]])[1]
    mixed1 <- Matrix(0, nrow = n, ncol= n)
    if(length(order[i][[1]])==1){
      num1 <- order[i][[1]][1]
      mixed1 <-  params1*tuning[[names(tuning)[num1]]]
    }else{
      k <- 1
      for(num1 in order[i][[1]]){
        mixed1 <- mixed1 + params1[k]*tuning[[names(tuning)[num1]]]
        k <- k+1
      }
    }
   
    mixed1S <- Matrix::Diagonal(dim(mixed1)[1])- 0.3*mixed1
    mixedNetworks[[paste0(c("params1",paste0(names(tuning)[order[i][[1]]], collapse="|")), collapse ="|")]] <- mixed1S
    
    params <- self$tuneList$Average_parameters[["params_avg"]][order[i][[1]]]
    params <- params/sum(params)
    mixed <- Matrix(0, nrow = n, ncol= n)
    if(length(order[i][[1]])==1){
      num <- order[i][[1]][1]
      mixed <-  params*tuning[[names(tuning)[num]]]
    }else{
      k <- 1
      for(num in order[i][[1]]){
        mixed <- mixed + params[k]*tuning[[names(tuning)[num]]]
        k <- k+1
      }
    }
    
    mixedS <- Diagonal(dim(mixed)[1])- 0.3*mixed
    mixedNetworks[[paste0(c("params",paste0(names(tuning)[order[i][[1]]], collapse="|")), collapse ="|")]] <- mixedS
  
    }
    save(mixedNetworks, file="backup/mixedNetworks.RData")
    invisible(self)
  },
  leaveOneOutCVMonoplex = function(){
    print("leaveOneOutCVMonoplex")
    load(file="backup/mixedNetworks.RData")
    self$ubiquitinList = read_excel(self$ubiquitinProteinHomo)
   
    db <- self$complexPortalDB
    for(name in names(mixedNetworks)){
    db[[name]] <- NA
    }
  
  Total_result <- data.frame(matrix(NA, nrow=0, ncol= length(colnames(db))+1))
  colnames(Total_result) <- c(colnames(db),"Protein")
  db_sel <- db %>% filter(count>2&checksum==1)
  for(j in 1:dim(db_sel)[1]){
    print(j)
    ds <- db_sel[j,]
    db_sel[j, "Identifiers"]
    complex_protein <- gsub("\\s*\\([^\\)]+\\)","",strsplit(as.character(db_sel$Identifiers[j]) ,"\\|")[[1]])
    complex_protein <- complex_protein[complex_protein %in% protein_list2]
    
    if(length(complex_protein)>1){
      res <- data.frame(cbind(ds[rep(seq_len(nrow(ds)), each=length(complex_protein)),], 
                              Protein =complex_protein))
      colnames(res) <- colnames(Total_result)
    }
    
    selected_protein <- c()
    go_id <- gsub("\\s*\\([^\\)]+\\)","",strsplit(as.character(res$Go.CC[1]) ,"\\|")[[1]])
    go_map <- data.frame(as_edgelist(g_GO2Protein, names = TRUE))
    go_map$X1 <- as.character(go_map$X1)
    go_map$X2 <- as.character(go_map$X2)
    
    selected_protein <- unique(c(unique(go_map$X2[go_map$X1 %in% go_id]),self$ubiquitinList[["Entry"]]))
    
      for(h in 1:length(names(mixedNetworks))){
        nam <- names(mixedNetworks)[h]
        AS <- mixedNetworks[[nam]]
        if(sum(complex_protein %in% selected_protein)==length(complex_protein)){
          AS1 <- AS[rownames(AS) %in% selected_protein, colnames(AS) %in% selected_protein]
        }else{
          AS1 <- AS
        }
       
        for(y in 1:length(complex_protein)){
          seed <- complex_protein[!complex_protein %in% complex_protein[y]]
          size <- dim(AS1)[1]
          p0 <- matrix(0,nrow=size, ncol=1)
          rownames(p0) <- colnames(AS1)
          p1=p0
          p1[rownames(p1) %in% seed,] <- 1.0/length(seed)
          b <- p1*0.7
          x0 <- p1
          run  <- conjgrad(AS1,b,x0)
          result <- data.frame(rownames(run),matrix(run))
          result <- result[order(-result[["matrix.run."]]),]
          result <- result[!result[["rownames.run."]] %in% seed,]
          rownames(result) <- NULL
          if(length(as.numeric(rownames(result[result[["rownames.run."]]==as.character(res[['Protein']][y]),])))>0){
            res[y,nam] <- as.numeric(rownames(result[as.character(result[["rownames.run."]])==as.character(res[["Protein"]][y]),]))
          }
        }
      }
    Total_result <- rbind(Total_result, res)
    print(dim(Total_result)[2])
  }
  self$totalResult <- Total_result
  self$ppiNet <- mixedNetworks[[names(mixedNetworks)[grepl("params1",names(mixedNetworks))&grepl("Prediction_G",names(mixedNetworks))&!grepl("DiseaseSimilarity",names(mixedNetworks))]]]

  save(Total_result, file= "backup/total_result.RData")
    invisible(self)
  },
  visualizeResult = function(){
    Total_result <- self$totalResult
    dd <- Total_result[colnames(Total_result)[grepl("params1",colnames(Total_result))]]
    colnames(dd) <- c("Physic_P",
                      "Physic_P+Pathway_P",
                      "Physic_P+Pathway_P+CrossLink",
                      "Physic_P+Pathway_P+CrossLink+DiseaseSim", 
                      "Physic_P+Pathway_P+CrossLink+Physic_G",
                      "Physic_P+Pathway_P+CrossLink+Physic_G+Coex_G",
                      "Physic_P+Pathway_P+CrossLink+Physic_G+Coex_G+Coloc_G", 
                      "Physic_P+Pathway_P+CrossLink+Physic_G+Coex_G+Coloc_G+Genetic_G",
                      "Physic_P+Pathway_P+CrossLink+Coex_G+Coloc_G+Genetic_G+Physic_G+Prediction_G",
                      "Physic_P+Pathway_P+CrossLink+Coex_G+Coloc_G+Genetic_G+Physic_G+Prediction_G+DiseaseSim")
    dt <- melt(dd)
    
    file.pdf = "image/LOOCV_with-prier_RR.pdf"
    pdf(file.pdf,width=8,height=10) 
    p1 <- ggplot(dt, aes(value)) +
      stat_ecdf(aes(color = variable)) +
      scale_fill_brewer(palette = "Dark2")+theme_bw()+
      theme(legend.direction = "vertical", legend.position = "bottom", legend.box = "vertical")+
      facet_zoom(xlim =c(0,1000), ylim=c(0.82,0.95))+
      ggtitle("Leave-One-Out Cross Validation") +
      xlab("Rank") + ylab("ecdf")
    print(p1)
    dev.off()
    
    dd2 <- Total_result[colnames(Total_result)[(!grepl("params1",colnames(Total_result)))&grepl("params",colnames(Total_result))]]
    colnames(dd2) <- colnames(dd)
    dt2 <- melt(dd2)
    
    file.pdf = "image/LOOCV_without-prier_RR.pdf"
    pdf(file.pdf,width=8,height=10) 
    p2 <- ggplot(dt2, aes(value)) +
      stat_ecdf(aes(color = variable)) +
      scale_fill_brewer(palette = "Dark2")+theme_bw()+
      theme(legend.direction = "vertical", legend.position = "bottom", legend.box = "vertical")+
      facet_zoom(xlim =c(0,1000), ylim=c(0.82,0.95))+
      ggtitle("Leave-One-Out Cross Validation") +
      xlab("Rank") + ylab("ecdf")
    print(p2)
    dev.off()
    
    x <- dd[[colnames(dd)[length(names(dd))]]]
    y <- dd2[[colnames(dd2)[length(names(dd2))]]]
    decdf <- function(x, y)  {ecdf(x)(1:length(x)) - ecdf(y)(1:length(y))}
    diff <- decdf(x, y)
    diff <- diff[diff!=0]
    df <- data.frame(diff) 
    

    ggplot(df, aes(x=diff, fill="#FF6666")) + 
      geom_histogram(aes(y=..density..),alpha=0.5, colour="#FF6666", fill="#FF6666")+
      theme_classic()+
      geom_vline(aes(xintercept=mean(diff)),
                  color="blue", linetype="dashed", size=1)+
      geom_vline(aes(xintercept=0),
                 color="green", size=1)
    
    dt <- data.frame(matrix(NA, nrow=0, ncol=9))
    colnames(dt) <- c(colnames(tune$tuneList[["Parameters 1"]]),"run")
    for (name in names(tune$tuneList)[grepl("Parameters",names(tune$tuneList))]){
      dt <- rbind(dt,cbind(tune$tuneList[[name]],run= rep(name,9)))
    }
    ggplot(dt, aes(fill=layer, y=prim, x=run)) + 
      geom_bar(position="fill", stat="identity")+
      theme(axis.text.x=element_text(angle=-90))
    
  }
  
  ))





