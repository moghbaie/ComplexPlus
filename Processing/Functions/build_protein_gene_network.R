# Mehrnoosh Oghbaie
# 08/30/2018
# Build gene and protein level networkS

Networks <- R6Class("Networks",  list(
  geneManiaInput = "../data/input_geneMania/networks.txt", # List of publication from GeneMania
  geneManiaNetworks = NA,
  pubmed.complex = NA,
  complexPortalHS = "../data/complexportal/records/H.sapiens.tsv", # ComplexPortal database directory
  identifierMappingGeneFile = "../data/mapping/db/genemania/identifier_mappings.txt",# Identifier gene mapping from GeneMania
  totalDirectInteractionDB = "../data/layers/physical_interaction/db/total_direct.RData", # Direct protein interactions downloaded from PSIQUIC
  oldDirectProtein = "../data/layers/physical_interaction/db/old_direct_protein.RData", # Old direct protein interactions used in GeneMnaia, BioPlex and ComplexpOrtal
  # Interactions from GeneMania and Bioplex have been saved in sql database 
  physicalInteractionDB = "../data/layers/physical_interaction/db/Physical_Interaction.sqlite", #physical interactions from GeneMania and BioPlex
  colocalizationDB = "../data/layers/co-localization/db/Co-localization.sqlite",
  predictionDB = "../data/layers/predicted/db/Predicted.sqlite",
  coexpressionDB = "../data/layers/co-expression/db/Co-expression.sqlite",
  geneticDB = "../data/layers/genetic_interaction/db/Genetic_Interaction.sqlite", # Genetic interaction database from GeneMania
  # Cross-link data from xLinkDB homo-sapiens
  xLinkDB = "../data/layers/xlink/db/xLinkDB_Homo.RData",
  # Initial graph
  oldDirectGene = "../data/layers/physical_interaction/db/old_direct_gene.RData", # Old direct gene interactions used in GeneMnaia, BioPlex and ComplexpOrtal,
  oldGeneGraph = "../data/layers/physical_interaction/graph/old_gene_PPI_complexportal.RData",
  oldProteinGraph = "../data/layers/physical_interaction/graph/old_protein_PPI_complexportal.RData",
  oldColocalizationGraph = "../data/layers/co-localization/graph/old_colocalization.RData",
  oldPredictionGraph = "../data/layers/predicted/graph/old_prediction.RData",
  oldGeneticGraph = "../data/layers/genetic_interaction/graph/old_genetic_interaction.RData",
  oldCoexpressionGraph = "../data/layers/co-expression/graph/old_Coexpression.RData",
  newProteinPathwayGraph = "../data/layers/pathway/graph/new_pathway_protein.RData",
  xLinkGraph = "../data/layers/xLink/graph/xlinkDB.RData",
  # Merged list
  geneList = NA,
  proteinList = NA,
  identifierMappingGene = NA,
  identifierMappingsProtein = NA,
  humanProteinAbundanceFile = "../data/mapping/db/9606-WHOLE_ORGANISM-integrated.txt", # Abundance ratio file from PaxDB
  identifierMappingsProteinFile = "../data/mapping/db/HUMAN_9606_idmapping.dat",
  geneBackup = "../data/layers/Gene_Backup.RData", # Final gene_level graphs
  proteinBackup = "../data/layers/Protein_Backup.RData", # Final protein_level graphs
  protein2GeneBipartiteMapping = "../data/mapping/graph/protein2Gene.RData",
  ## Functions that were done before: 
  # 1. Bioplex and GeneMania are merged.
  # 2. Result of each cathegory is saved in Layers under their directory
  getGeneManiaNetworks = function(){
    # GeneMania data repository http://genemania.org/data/current/Homo_sapiens/, date = 2017-03-14
    # BioPlex data repository https://bioplex.hms.harvard.edu/downloadInteractions.php, only data from the first three studies
    # Reading GeneMania database references table (included BioPlex-version 3.1)
    self$geneManiaNetworks = read.csv(self$geneManiaInput, sep="")
    invisible(self)
  },
  extractPubmedComplexPortal = function(){
    # Reading homo sapien complextab from ComplexPortal website
    # Getting the latest Homo Sapiens ComplexTab from https://www.ebi.ac.uk/complexportal/complex/organisms
    H.sapiens = read.delim(self$complexPortalHS)
    self$pubmed.complex = unique(unlist(apply(H.sapiens,1, function(y) lapply(strsplit(as.character(y[["Cross.references"]]),"\\|")[[1]][grepl("pubmed:",
                                                                                                                                           strsplit(as.character(y[["Cross.references"]]),"\\|")[[1]])],
                                                                          function(x) gsub("\\s*\\([^\\)]+\\)","",strsplit(x,"\\:")[[1]][2]))
                                          )
                                    ))
    invisible(self)
  },
  getDirectInteraction = function(){
    ################################################################################################################################
    #Download direct molecular interactions from PSIQUIC
    psicquic = PSICQUIC()
    DB = providers(psicquic)
    #################################################################################################################################
    tbl_direct = PSICQUIC::interactions(psicquic, species="9606",type="direct interaction", provider=DB,speciesExclusive= TRUE)
    tbl_direct = addGeneInfo(IDMapper("9606"),tbl_direct) 
    tbl_direct = createEmptyCol(tbl_direct)
    lst = strsplit(tbl_direct$publicationID, "\\|")
    tbl_direct2 = data.frame(t(apply(tbl_direct, 1, fillSource)))
    tbl_direct2$detectionMethod = apply(tbl_direct2,1,fillMethod)
    tbl_direct2$confidence = NA
    colnames(tbl_direct3) = colnames(tbl_direct2)
    tbl_direct3 = data.frame(t(apply(tbl_direct2, 1, fillconfidenceScore)))
    tbl_direct3$pubmed = as.integer(as.character(tbl_direct3$pubmed))
    tbl_direct4 <- tbl_direct3
    tbl_direct4[is.na(tbl_direct4$year),"year"] <- apply(tbl_direct4[is.na(tbl_direct4$year),],1, function(x) 
      as.integer(gsub("[\\(\\)]", "", regmatches(as.character(x[["firstAuthor"]]), gregexpr("\\(.*?\\)", as.character(x[["firstAuthor"]])))[[1]])))
    
  #  tbl_direct4 = merge(tbl_direct3, join_list, by.x="pubmed", by.y="pmid", all.x=TRUE)
    
    ################################################################################################################################
    ### Save the final table in a Rdata file
    save(tbl_direct4 , file=self$totalDirectInteractionDB)
    invisible(self)
  },
  saveOldProteinInteraction = function(){
    # Save the old protein interactions 
    load(self$totalDirectInteractionDB)
    # Get pubmed ID that exists in GeneMania
    networks_Physic = self$geneManiaNetworks$Pubmed_ID[self$geneManiaNetworks$Network_Group_Name=="Physical Interactions"&self$geneManiaNetworks$Pubmed_ID!=0]
   
     tbl_direct4[is.na(tbl_direct4$year),"year"] <- apply(tbl_direct4[is.na(tbl_direct4$year),],1, function(x) 
      as.integer(gsub("[\\(\\)]", "", regmatches(as.character(x[["firstAuthor"]]), gregexpr("\\(.*?\\)", as.character(x[["firstAuthor"]])))[[1]])))
    
     ## also adding exosome publication c("15231747","17412707","26496610","11719186", "29298432", "26496610")
    old_direct_protein = tbl_direct4 %>%
      filter( grepl("uniprotkb",A)) %>%
      filter(pubmed %in% unique(c(networks_Physic,"15231747","17412707","26496610","11719186", "29298432", "26496610",self$pubmed.complex))) %>%
      dplyr::select(A,B, confidence, pubmed,detectionMethod, firstAuthor,year) %>%
      distinct()
    
    #unique(old_direct_protein$pubmed)[unique(old_direct_protein$pubmed) %in% self$pubmed.complex]
    save(old_direct_protein, file=self$oldDirectProtein)
    invisible(self)
  }, 
  getPhysicalInetractionfromDB = function(){
    load(self$oldDirectProtein)
    protein_pmid = unique(old_direct_protein[["pubmed"]])
    #networks_Physic = self$geneManiaNetworks$Pubmed_ID[self$geneManiaNetworks$Network_Group_Name=="Physical Interactions"&self$geneManiaNetworks$Pubmed_ID!=0]
    ppi = dbConnect(SQLite(), dbname=self$physicalInteractionDB)
    dbListTables(ppi)
    dbListFields(ppi,"Physical_Interaction_records")
    ppi_list = dbGetQuery(ppi, "SELECT Gene_A, Gene_B, Weight, Source FROM Physical_Interaction_records")
    dbDisconnect(ppi)
    old_direct_gene = ppi_list[!as.integer(ppi_list$Source) %in% protein_pmid,]
    save(old_direct_gene, file=self$oldDirectGene)
    invisible(self)
  },
  getIdentifierMapping = function(){
    # Identifier mapping file came from GeneMania. Perhaps, I should write with biomaRT
    self$identifierMappingGene = read.delim(self$identifierMappingGeneFile)
    invisible(self)
  },
  buildOldDirectGeneGraph = function(){
    ## Build the old direct gene network
    load(self$oldDirectGene)
    old_direct_gene2 = Merge_parallel(old_direct_gene)
    gene_PPI_graph = graph_from_data_frame(old_direct_gene2[c("Gene_A","Gene_B")], directed = FALSE)
    E(gene_PPI_graph)$Weight = as.numeric(as.character(old_direct_gene2$score_max))
    E(gene_PPI_graph)$Source = old_direct_gene2$Source
    save(gene_PPI_graph, file=self$oldGeneGraph)
    invisible(self)
  },
  buildOldDirectProteinNet = function(){
    # Build old direct protein network
    load(self$oldDirectProtein)
    old_direct_protein2 = old_direct_protein %>% 
      mutate("Gene_A"=as.character(gsub("uniprotkb:","",A)), "Gene_B"=as.character(gsub("uniprotkb:","",B)), "Weight"=as.numeric(as.character(confidence)), "Source"=pubmed) %>%
      dplyr::select(Gene_A, Gene_B, Weight, Source) %>%
      filter(!grepl("refseq:",Gene_B)&!grepl("refseq:",Gene_A)& !is.na(Weight)) %>% mutate(Gene_B=gsub("-","",Gene_B),Gene_A=gsub("-","",Gene_A))
    #old_direct_protein2 = rbind(old_direct_protein2, self$exosomeLinks)
    
    ## Merge multiple weight
    protein_PPI_new_list <- convertUniprot2Swisspro(unique(unlist(old_direct_protein2[,c("Gene_A","Gene_B")])))
    protein_PPI_new_list <- protein_PPI_new_list[protein_PPI_new_list$uniprotswissprot!="",]
    old_direct_protein2$Gene_A <- protein_PPI_new_list$uniprotswissprot[match(old_direct_protein2$Gene_A,protein_PPI_new_list$uniprot_gn)]
    old_direct_protein2$Gene_B <- protein_PPI_new_list$uniprotswissprot[match(old_direct_protein2$Gene_B,protein_PPI_new_list$uniprot_gn)]
    old_direct_protein2 <- old_direct_protein2%>% filter(!is.na(Gene_A)&!is.na(Gene_B))
    old_direct_protein3 <- Merge_parallel(old_direct_protein2)
    
    protein_PPI_graph <- graph_from_data_frame(old_direct_protein3 [c("Gene_A","Gene_B")], directed = FALSE)
    E(protein_PPI_graph)$Weight <- old_direct_protein3$score_max
    E(protein_PPI_graph)$Source <- as.character(old_direct_protein3$Source)
    
    save(protein_PPI_graph, file=self$oldProteinGraph)
    invisible(self)
  },
  buildColocalizationNet = function(){
    ## Build old co-localization network from GeneMania
    col = dbConnect(SQLite(), dbname=self$colocalizationDB)
    dbListTables(col)
    dbListFields(col,"Co_localization_records")
    col_list <- dbGetQuery(col, "SELECT Gene_A, Gene_B, Weight, Source FROM Co_localization_records")
    dbDisconnect(col)
    #### Build graph from both old colocalization data
    ## Merge multiple weight
    col_list2 <- Merge_parallel(col_list)
    colocalization_graph <- graph_from_data_frame(col_list2[c("Gene_A","Gene_B")], directed = FALSE)
    E(colocalization_graph)$Weight <- col_list2$score_max
    E(colocalization_graph)$Source <- col_list2$Source
    save(colocalization_graph, file = self$oldColocalizationGraph)
    invisible(self)
  },
  buildPredictionNet = function(){
    # Reading prediction records from Sqlite
    db = dbConnect(SQLite(), dbname=self$predictionDB)
    dbListFields(db, "Predicted_records")
    prediction_list <- dbGetQuery(db, "SELECT Gene_A, Gene_B, Weight, Source FROM Predicted_records")
    dbDisconnect(db)
    
    prediction_list$Weight <- as.numeric(as.character(prediction_list$Weight))
    prediction_list <- prediction_list[complete.cases(prediction_list[c("Gene_A","Gene_B")]),]
    ## Merge multiple weight
    prediction_list2 <- Merge_parallel(prediction_list)
    
    prediction_graph <- graph_from_data_frame(prediction_list2[c("Gene_A","Gene_B")], directed = FALSE)
    E(prediction_graph)$Weight <- prediction_list2$score_max
    E(prediction_graph)$Source <- prediction_list2$Source

    save(prediction_graph, file=self$oldPredictionGraph)
    invisible(self)
  },
  buildGeneticNet = function(){
    
    gen = dbConnect(SQLite(), dbname=self$geneticDB)
    dbListTables(gen)
    dbListFields(gen, "Genetic_Interaction_records")
    Genetic_list <- dbGetQuery(gen, "SELECT Gene_A, Gene_B, Weight, Source FROM Genetic_Interaction_records")
    dbDisconnect(gen)
    
    ## Merge multiple weight
    Genetic_list2 <- Merge_parallel(Genetic_list)
    genetic_graph <- graph_from_data_frame(Genetic_list2[c("Gene_A","Gene_B")], directed = FALSE)
    E(genetic_graph)$Weight <- Genetic_list2$score_max
    E(genetic_graph )$Source <-Genetic_list2$Source
    save(genetic_graph, file=self$oldGeneticGraph)
    invisible(self)
  }, 
  buildCoexpressionNet = function(){
    co = dbConnect(SQLite(), dbname=self$coexpressionDB)
    dbListTables(co)
    dbListFields(co, "Co_expression_records")
    co_exp <- dbGetQuery(co, "SELECT Gene_A, Gene_B, Weight, Source FROM Co_expression_records")
    dbDisconnect(co)
    ## Merge multiple weight
    co_exp2 <- Merge_parallel(co_exp)
    coexpression_graph <- graph_from_data_frame(co_exp2[c("Gene_A","Gene_B")], directed = FALSE)
    E(coexpression_graph )$Weight <- co_exp2$score_max
    E(coexpression_graph)$Source <- co_exp2$Source
    save(coexpression_graph, file=self$oldCoexpressionGraph)
    invisible(self)
  },
  buildPathwayNet = function(){
    ##################################################################################################
    # Build pathway graph from nci, reactome and panther
    pathwayDatabases() 
    # We create Kegg interaction Network using the GetPathwayData funciton. 
    ################################################################################################################
    nciRaw <- GetPathwayData ("nci","hsapiens")
    pantherRaw <- GetPathwayData ("panther","hsapiens")
    reactomeRaw <- GetPathwayData ("reactome","hsapiens")
    
    protein_Pathway_list <- convertUniprot2Swisspro(unique(c(unlist(nciRaw),unlist(pantherRaw ),unlist(reactomeRaw))))
    protein_Pathway_list <- protein_Pathway_list[protein_Pathway_list$uniprotswissprot!="",]
    
    nciRaw[,1] <- protein_Pathway_list$uniprotswissprot[match(as.character(nciRaw[,1]),protein_Pathway_list$uniprot_gn)]
    nciRaw[,2]<- protein_Pathway_list$uniprotswissprot[match(as.character(nciRaw[,2]),protein_Pathway_list$uniprot_gn)]
    nciRaw <- nciRaw[!is.na(nciRaw[,1])& !is.na(nciRaw[,2]),]
    
    pantherRaw[,1] <- protein_Pathway_list$uniprotswissprot[match(as.character(pantherRaw[,1]),protein_Pathway_list$uniprot_gn)]
    pantherRaw[,2]<- protein_Pathway_list$uniprotswissprot[match(as.character(pantherRaw[,2]),protein_Pathway_list$uniprot_gn)]
    pantherRaw <- pantherRaw[!is.na(pantherRaw[,1])& !is.na(pantherRaw[,2]),]
    
    reactomeRaw[,1] <- protein_Pathway_list$uniprotswissprot[match(as.character(reactomeRaw[,1]),protein_Pathway_list$uniprot_gn)]
    reactomeRaw[,2]<- protein_Pathway_list$uniprotswissprot[match(as.character(reactomeRaw[,2]),protein_Pathway_list$uniprot_gn)]
    reactomeRaw <- reactomeRaw[!is.na(reactomeRaw[,1])& !is.na(reactomeRaw[,2]),]
    
    
    # Build network.
    nci <- build.network(nciRaw)
    panther <- build.network(pantherRaw)
    reactome <- build.network(reactomeRaw)
    
    PathwayRaw <- graph.union(reactome, nci,panther)
    Pathway <- igraph::simplify(PathwayRaw)
    
    save(Pathway, file = self$newProteinPathwayGraph)
    invisible(self)
  },
  buildxLinkNet = function(){
    # Download xLink interaction data from http://xlinkdb.gs.washington.edu/xlinkdb/
    load(file=self$xLinkDB)
    xl_df_hm$distance <- as.numeric(as.character(xl_df_hm$distance))
    xl_net <- xl_df_hm %>% 
      dplyr::select(accessionA, accessionB, distance) %>% 
      dplyr::filter(!is.na(distance)) %>%
      group_by(accessionA,accessionB)%>%
      dplyr::summarize(distance = mean(distance))%>%
      dplyr::filter(distance < 25.0) %>%
      distinct()
    xl_graph <- graph_from_data_frame(xl_net[c("accessionA", "accessionB")], directed = FALSE)
    save(xl_graph, file = self$xLinkGraph)
    invisible(self)
  },
  unionGeneList = function(){
    file_names = c(self$oldGeneGraph,
                   self$oldColocalizationGraph, 
                   self$oldPredictionGraph,
                   self$oldGeneticGraph,
                   self$oldCoexpressionGraph )
    lapply(file_names,load,.GlobalEnv)
    gene_list = unique(c(as_ids(V(gene_PPI_graph)), as_ids(V(coexpression_graph)), as_ids(V(genetic_graph)), as_ids(V(prediction_graph)), as_ids(V(colocalization_graph))))
    
    gene_PPI_graph = add.vertices(gene_PPI_graph, nv=length(gene_list[!gene_list %in% as_ids(V(gene_PPI_graph))]),
                                   attr = list(name=gene_list[!gene_list %in% as_ids(V(gene_PPI_graph))]))
    
    coexpression_graph = add.vertices(coexpression_graph, nv=length(gene_list[!gene_list %in% as_ids(V(coexpression_graph))]),
                                       attr = list(name=gene_list[!gene_list %in% as_ids(V(coexpression_graph))]))
    
    genetic_graph = add.vertices(genetic_graph, nv=length(gene_list[!gene_list %in% as_ids(V(genetic_graph))]),
                                  attr = list(name=gene_list[!gene_list %in% as_ids(V(genetic_graph))]))
    
    prediction_graph = add.vertices(prediction_graph, nv=length(gene_list[!gene_list %in% as_ids(V(prediction_graph))]),
                                     attr = list(name=gene_list[!gene_list %in% as_ids(V(prediction_graph))]))
    
    colocalization_graph = add.vertices(colocalization_graph, nv=length(gene_list[!gene_list %in% as_ids(V(colocalization_graph))]),
                                         attr = list(name=gene_list[!gene_list %in% as_ids(V(colocalization_graph))]))
    self$geneList = gene_list
    save(gene_list, gene_PPI_graph,coexpression_graph,genetic_graph, prediction_graph, colocalization_graph, 
         file=self$geneBackup)
    invisible(self)
    },
  unionProteinList = function(){
    
    file_names = c(self$oldProteinGraph,
                   self$newProteinPathwayGraph,
                   self$xLinkGraph)
    lapply(file_names,load,.GlobalEnv)
    
    self$proteinList<- unique(c(as_ids(V(Pathway)), as_ids(V(protein_PPI_graph)),as_ids(V(xl_graph))))
    ####################################################################################################
    ### Union the proteins from protein-level network and those from gene-level that exists abundance rate    
    ## The file was downloaded from https://pax-db.org/download
    X9606_WHOLE_ORGANISM_integrated <- read_delim(self$humanProteinAbundanceFile, 
                                                  "\t", escape_double = FALSE, trim_ws = TRUE)
    
    HUMAN_9606_idmapping <- read_delim(self$identifierMappingsProteinFile, 
                                       "\t", escape_double = FALSE, col_names = FALSE, 
                                       trim_ws = TRUE)
    self$identifierMappingsProtein <- HUMAN_9606_idmapping
    X9606_WHOLE_ORGANISM_integrated$uniprotID <- HUMAN_9606_idmapping$X1[match(X9606_WHOLE_ORGANISM_integrated$string_external_id, HUMAN_9606_idmapping$X3)]
    swiss_uniprot_list <- convertUniprot2Swisspro(X9606_WHOLE_ORGANISM_integrated$uniprotID)
    X9606_WHOLE_ORGANISM_integrated$uniprotswissprot <- swiss_uniprot_list$uniprotswissprot[match(X9606_WHOLE_ORGANISM_integrated$uniprotID, swiss_uniprot_list$uniprot_gn)]
    
    ## Write a convert list from ensembl gene to swissprot
    ensembl_swisspro_list <- convertEnsembl2Swisspro(self$geneList)
    ensembl_swisspro_list <- ensembl_swisspro_list[ensembl_swisspro_list$uniprotswissprot!="",]
    
    ### Find genes with abundance value
    ## map it first to gene_list
    total_list_gene <- merge(x=X9606_WHOLE_ORGANISM_integrated,y=ensembl_swisspro_list,by="uniprotswissprot",all.y=TRUE)
    total_list_gene1 <- total_list_gene%>% dplyr::select(uniprotswissprot,ensembl_gene_id,abundance)
    total_list_gene1$abundance[is.na(total_list_gene1$abundance)] <- 0
    
    ## Convert ProteinList to ensemble gene 
    swisspro_ensemble_protein_list <- convertSwisspro2ensembl(self$proteinList)
    ## merge the two list
    total_list_protein <- merge(x=X9606_WHOLE_ORGANISM_integrated,y=swisspro_ensemble_protein_list,by="uniprotswissprot",all.y=TRUE)
    total_list_protein1 <- total_list_protein %>% dplyr::select(uniprotswissprot,ensembl_gene_id,abundance)
    total_list_protein1$abundance[is.na(total_list_protein1$abundance)] <- 0
    
    total_protein_gene_list <- unique(rbind(total_list_protein1,total_list_gene1))
    total_protein_gene_list$abundance[is.na(total_protein_gene_list$abundance)]<- 0
    total_abundance <- total_protein_gene_list %>% distinct() %>% dplyr::group_by(ensembl_gene_id) %>% dplyr::summarise(abundance_total =sum(abundance,na.rm=T), count=n())
    total_protein_gene_list2 <- merge(x=total_protein_gene_list,y=total_abundance,by="ensembl_gene_id",all.x=TRUE)
    total_protein_gene_list2$ratio <- total_protein_gene_list2$abundance/total_protein_gene_list2$abundance_total
    total_protein_gene_list2$ratio[total_protein_gene_list2$count==1] <- 1
    total_protein_gene_list2$ratio[is.na(total_protein_gene_list2$ratio)] <- 1/total_protein_gene_list2$count[is.na(total_protein_gene_list2$ratio)]
    
    ## Final protein list
    protein_list2 <- unique(total_protein_gene_list$uniprotswissprot)
    self$proteinList <- protein_list2
    
    common_list <- total_protein_gene_list2 %>% filter(ensembl_gene_id %in% self$geneList)
    different_list <- total_protein_gene_list2 %>% filter(!ensembl_gene_id %in% self$geneList)
    
    protein2Gene <- graph_from_data_frame(common_list[common_list$ratio!=0,c("uniprotswissprot","ensembl_gene_id")], directed = FALSE)
    E(protein2Gene)$Weight <- common_list$ratio[common_list$ratio!=0]
    
    protein2Gene <- add.vertices(protein2Gene, nv=length(self$geneList[!self$geneList %in% common_list$ensembl_gene_id]),
                                 attr = list(name=self$geneList[!self$geneList %in% common_list$ensembl_gene_id]))
    protein2Gene <- add.vertices(protein2Gene, nv=length(self$proteinList[!self$proteinList %in% common_list$ensembl_gene_id]),
                                 attr = list(name=self$proteinList[!self$proteinList %in% common_list$ensembl_gene_id]))
    
    save(protein2Gene, file=self$protein2GeneBipartiteMapping)
    
    
    Pathway <- add.vertices(Pathway, nv=length(protein_list2[!protein_list2 %in% as_ids(V(Pathway))]),
                            attr = list(name=protein_list2[!protein_list2 %in% as_ids(V(Pathway))]))
    
    protein_PPI_graph <- add.vertices(protein_PPI_graph, nv=length(protein_list2[!protein_list2 %in% as_ids(V(protein_PPI_graph))]),
                                      attr = list(name=protein_list2[!protein_list2 %in% as_ids(V(protein_PPI_graph))]))
    
    xl_graph <- add.vertices(xl_graph, nv=length(protein_list2[!protein_list2 %in% as_ids(V(xl_graph))]),
                             attr = list(name=protein_list2[!protein_list2 %in% as_ids(V(xl_graph))]))
    
    save(protein_list2, Pathway, protein_PPI_graph, xl_graph,
         file=self$proteinBackup)
    
    invisible(self)
  }
  )
)






