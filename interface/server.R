


source("functions.R", local=TRUE)
load(file="db/interfacedbBackup.RData")

ChEBI_Results <- interfaceDB$ChEBI_Results
idmapping <- interfaceDB$idmapping
db_hm <- interfaceDB$db_hm
Total_edges <- interfaceDB$Total_edges
gM1 <- interfaceDB$gM1
Complex_protein <- interfaceDB$complex_protein
p0 <- interfaceDB$p0
A <- interfaceDB$A
uniprot_ubiquitin_H_sapiens <- interfaceDB$uniprot_ubiquitin_H_sapiens
protein_list <- interfaceDB$protein_list
PPI_layer <- interfaceDB$PPI_layer
M_xl_graph_corrected.normalized <- interfaceDB$M_xl_graph_corrected.normalized
M_Pathway.normalized <- interfaceDB$M_Pathway.normalized
M_Disease2Protein <- interfaceDB$M_Disease2Protein
mim_morbid <- interfaceDB$mim_morbid
goa <- interfaceDB$goa

rm(list=c("interfaceDB"))

load(file="db/directinteraction.RData")
directInteraction <- tbl_direct4
rm(list=c("tbl_direct4"))


server <- function(input, output,session) {
  
  hintjs(session, options = list("hintButtonLabel"="Hope this hint was helpful"),
         events = list("onhintclose"=I('alert("Wasn\'t that hint helpful")')))
  
  observeEvent(input$help,
               introjs(session, options = list("nextLabel"="Onwards and Upwards",
                                               "prevLabel"="Did you forget something?",
                                               "skipLabel"="Don't be a quitter"),
                       events = list(onbeforechange = readCallback("switchTabs")))
  )
  
  d <- reactive({
    dist <- switch(input$species,
                   hm = hm,
                   sc = sc)
  })
  
  
  output$value <- renderText({ input$num })
  
  BuildDatabaseC <- reactive({
    species <- input$species
    if(species=="hm"){
      df <- db_hm[db_hm$number.links!=0,]
    } 
    else if(species=="sc"){
    }
    
    # Pruning homo_sapiens
    databaseC <- df %>% 
      dplyr::select(X.Complex.ac,
                    Recommended.name,
                    Description, 
                    Identifiers) %>%
      dplyr::rename("Complex_ID"=X.Complex.ac,"Name"=Recommended.name,
                    "Members"=Identifiers)%>%
      dplyr::mutate(Members=gsub("\\(.)*","",Members))
    return(databaseC)
  })
  # Filter data based on selections
  output$table <- DT::renderDataTable(DT::datatable({
    data <-  BuildDatabaseC()
    return(data[,c("Complex_ID","Name")])
  }), 
  selection = 'single' )
  
  #Calibrate = reactive({
  #    s = input$table_rows_selected[1]
  #    df = BuildDatabaseC()
  #    if(length(s)){
  #      name = unname(df$Name[s])
  #    }else{
  #      name = unname(df$Name[461])
  #    }
  #    return(df[s,"Name"])
  #  })
  
  output$CalibrateOptions <- renderUI({
    selectizeInput(
      inputId ="calibrate",
      label = "Add other proteins to complex and recalculate:",
      choices = protein_list[["external_gene_name"]],
      multiple = TRUE, 
      selected = NULL,  
      options = list(maxItems = 30,'plugins' = list('remove_button'))
    )
  })
  
  #updateSelectizeInput(session, "calibrate",  choices = data, server = TRUE)
  
  output$Description <- renderUI({ 
    s = input$table_rows_selected[1]
    df = BuildDatabaseC()
    if(length(s)){
      Description = unname(df$Description[s])
      desc <- paste("<h3>Function</h3>","<div>",Description,"</div>")
      HTML(desc)
    }
    #else{
    #  Description = unname(df$Description[461])
    #}
    
  })
  
  callback = "function(table) {
      table.on('click.dt', 'tr', function() {
            table.$('tr.selected').removeClass('selected');
            $(this).toggleClass('selected');            
        Shiny.onInputChange('rows',
                            table.rows('.selected').data()[0][0]);
      });
    }"
  
  callback = "function(table3) {
    table.on('click.dt', 'tr', function() {
      table.$('tr.selected').removeClass('selected');
      $(this).toggleClass('selected');            
    Shiny.onInputChange('rows',
            table.rows('.selected').data()[0][0]);
});
  }"
  
#  UploadMappling <- reactive({
#    species <- input$species
#    idmapping <- read.delim(paste0(species,"/identifier_mappings.txt"))
#    idmapping <- idmapping %>% mutate(Preferred_Name = as.character(Preferred_Name),
#                                      Name = as.character(Name),
#                                      Source = as.character(Source))
#    return(idmapping)
#  })
  
  convertComplex <- reactive({
    s = input$table_rows_selected[1]
    if(length(s)){
      lst = unlist(strsplit(BuildDatabaseC()$Members[s], "[|]"))
      #idmapping <- UploadMappling
      members <- unlist(lapply(lst, UniProt2Gene))
      members <- members[!is.na(members )]
      return(members)
    }
    #else{
    #  lst = unlist(strsplit(BuildDatabaseC()$Members[461], "[|]"))
    #}
  })
  
  output$complex = renderPrint({
    s = input$table_rows_selected[1]
    if(length(s)){
      name = BuildDatabaseC()$Name[s]
      members <- convertComplex()
      cat(paste0("Complex ",name," is selected.\n The memebers are:"))
      cat(members, sep = ", ")
    }
    #else{
    #  name = BuildDatabaseC()$Name[461]
    #}
  })
  
  makeInitialGraph <- reactive({
    s = input$table_rows_selected[1]
    if(length(s)){
      Complex_ID = BuildDatabaseC()$Complex_ID[s]
      Members <- db_hm$Identifiers[db_hm$X.Complex.ac==Complex_ID]
      nodes=unique(strsplit(gsub("\\(.)*","",Members),"\\|")[[1]])
      edges <- Total_edges %>% dplyr::filter(complexID == Complex_ID)
      chebi_select <- ChEBI_Results[ChEBI_Results$ID %in%  unique(nodes[grepl("CHEBI:",nodes)]),]
      chemi_complex <- unique(c(edges$protein_A[is.na(edges$uniprot_A)],edges$protein_B[is.na(edges$uniprot_B)]))
      
      if(length(chebi_select[["NAME"]])==1){
        edges$uniprot_B[is.na(edges$uniprot_B)] <- chebi_select[["ID"]]
        edges$uniprot_A[is.na(edges$uniprot_A)] <- chebi_select[["ID"]]
      }
      
      if(length(chebi_select[["NAME"]])>1){
        edges$uniprot_B[is.na(edges$uniprot_B)] <- unlist(lapply(edges$protein_B[is.na(edges$uniprot_B)], function(x) chebi_select[["ID"]][grepl(x, chebi_select[["NAME"]])]))
        edges$uniprot_A[is.na(edges$uniprot_A)] <- unlist(lapply(edges$protein_A[is.na(edges$uniprot_A)], function(x) chebi_select[["ID"]][grepl(x, chebi_select[["NAME"]])]))
      }
      
      edges[grepl("pro",edges$protein_A),"uniprot_A"] <- edges[grepl("pro",edges$protein_A),"protein_A"] 
      edges[grepl("pro",edges$protein_B),"uniprot_B"] <- edges[grepl("pro",edges$protein_B),"protein_B"] 
      
      edges2 <- edges%>% dplyr::select("uniprot_A", "uniprot_B") %>% mutate(uniprot_A = toupper(uniprot_A),uniprot_B = toupper(uniprot_B)) %>% distinct()
      
      nodes_geneName <- protein_list[["external_gene_name"]][match(nodes,protein_list[["uniprotswissprot"]])]
      nodesDS <- data.frame(nodes, nodes_geneName)
      node_uniprot <- lapply(as.character(nodesDS$nodes), function(x) strsplit(x,"-")[[1]][1])
      nodes3 <- data.frame(id = nodes,
                           label = ifelse(!is.na(nodesDS$nodes_geneName),as.character(nodesDS$nodes_geneName),as.character(nodesDS$nodes)),
                           shape = ifelse( grepl("-PRO_",nodesDS$nodes),"triangle",
                                           ifelse(grepl("CHEBI:",node_uniprot),"diamond","dot")), 
                           color = list(background = "linen", 
                                        border = "black",
                                        highlight = "yellow"),
                           font = list(size=12, bold= TRUE),
                           shadow = list(enabled = TRUE, size = 12))
      
      edges3 <- data.frame(from = edges2[["uniprot_A"]], to = edges2[["uniprot_B"]])
    
      return(list( nodes3,edges3))
    }
    #else{
    #  Complex_ID = BuildDatabaseC()$Complex_ID[461]
    #}
    
  })
  
  output$complexMembers = renderVisNetwork({
    s = input$table_rows_selected[1]
    if(length(s)){
      list_nodes_edges <- makeInitialGraph()
      list_edges <- list_nodes_edges[[2]]
      list_edges$smooth <- FALSE
      
      visNetwork(list_nodes_edges[[1]], list_edges, height = "1000px", width = "100%",
                 main = list(text = paste0(BuildDatabaseC()$Name[s]), style = "color:black;font-size:18px;text-align:center;"))%>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visPhysics(solver = "barnesHut",maxVelocity = 25, minVelocity = 20)%>%
        visLayout(randomSeed = 123)
    }
  })
  
  makeDiseasePlot = reactive({
    s = input$table_rows_selected[1]
    if(length(s)){
      Complex_ID = BuildDatabaseC()$Complex_ID[s]
      Members <- db_hm$Identifiers[db_hm$X.Complex.ac==Complex_ID]
      nodes=as.character(unique(strsplit(gsub("\\(.)*","",Members),"\\|")[[1]]))
      node_uniprot <- lapply(as.character(nodes), function(x) strsplit(x,"-")[[1]][1])
      
      if(sum(M_Disease2Protein[rownames(M_Disease2Protein) %in% node_uniprot,])){
        colsum <- Matrix::colSums(M_Disease2Protein[rownames(M_Disease2Protein) %in% node_uniprot,])
        selected_disease <- M_Disease2Protein[rownames(M_Disease2Protein) %in% node_uniprot, colnames(M_Disease2Protein) %in% names(colsum[colsum!=0])]
        if(class(selected_disease)=="dgCMatrix"){
          uniprot_disease <- rownames(selected_disease)
          
        } else{uniprot_disease <- names(selected_disease)}
        selected_disease <- matrix(selected_disease, ncol= length(colsum[colsum!=0]))
        colnames(selected_disease) <- lapply(mim_morbid[["mim_morbid_description"]][ match(names(colsum[colsum!=0]),mim_morbid[["mim_morbid_accession"]])], function(x) strsplit(x, ";;")[[1]][1])
        rownames(selected_disease) <- as.character(protein_list[["external_gene_name"]])[match(as.character(uniprot_disease),as.character(protein_list[["uniprotswissprot"]]))]
        
        p <- ggplot(reshape2::melt(selected_disease), aes(x=as.character(Var1), y=as.character(Var2), fill = value)) + 
          geom_tile(colour = "black") + 
          scale_fill_gradient(low = "white", high = "#00CC99", limits=c(0,1))+
          #coord_flip()+ 
          coord_equal()+
          theme_minimal()+
          ylab("Disease") + xlab("Gene Name")+
          theme(axis.text.x = element_text(angle=90,size = 11,face="bold"), 
                axis.text.y = element_text(size = 11, face="bold"),
                axis.title=element_text(size=14,face="bold"),
                legend.position = "none",
                plot.title = element_text(size = 18, face = "bold"))+ ggtitle("Gene disease relationship")
        return(p)
      }
    }
    #else{
    #  Complex_ID = BuildDatabaseC()$Complex_ID[461]
    #}
    
  })
  
  output$DiseasePlot = renderPlot({
    s = input$table_rows_selected[1]
    if(length(s)){
      p <- makeDiseasePlot()
      p
    }
  },width=1200)
  heights = reactive({
    s= input$table_rows_selected[1]
    heigh = ifelse(exists("s"), ifelse(stringr::str_count(db_hm$Identifiers[s], pattern = "\\|")>4,200+(stringr::str_count(db_hm$Identifiers[s], pattern = "\\|")-4)*50,200),200)
    return(heights)
  })
  
  
  output$Schematics = renderPlot({
    s = input$table_rows_selected[1]
    if(length(s)){
      Complex_ID = BuildDatabaseC()$Complex_ID[s]
      Members <- db_hm$Identifiers[db_hm$X.Complex.ac==Complex_ID]
      nodes= unique(strsplit(gsub("\\(.)*","",Members),"\\|")[[1]])
      nodes_swissprot <- nodes[nodes %in% protein_list[["uniprotswissprot"]]]
      
      list_selected <- paste0(nodes_swissprot,collapse=" ")
      prot_data <- drawProteins::get_features(list_selected)
      prot_data <- drawProteins::feature_to_dataframe(prot_data)
      
      p <- drawProteins::draw_canvas(prot_data)
      p <- drawProteins::draw_chains(p, prot_data)
      p <- drawProteins::draw_domains(p, prot_data)
      p <- drawProteins::draw_repeat(p, prot_data)
      p <- drawProteins::draw_motif(p, prot_data)
      #p <- drawProteins::draw_recept_dom(p, prot_data)
      p <- drawProteins::draw_phospho(p, prot_data, size = 2)
      # background and y-axis
      p <- p + theme_bw(base_size = 16) + # white backgnd & change text size
        theme(panel.grid.minor=element_blank(),
              panel.grid.major=element_blank()) +
        theme(axis.ticks = element_blank(),
              axis.text.y = element_blank()) +
        theme(panel.border = element_blank())
      # add titles
      rel_subtitle <- paste0("circles = phosphorylation sites\n")
      p <- p + labs(title = paste0("Schematic of ", db_hm$Recommended.name[db_hm$X.Complex.ac==Complex_ID]," proteins"),
                    subtitle = rel_subtitle)
      # move legend to top
      p <- p + theme(legend.position="top") + labs(fill="")
      p
    }
    #else{
    #  Complex_ID = BuildDatabaseC()$Complex_ID[461]
    #}
    #Complex_ID = BuildDatabaseC()$Complex_ID[s]
    
  },  height= reactive({
    s= input$table_rows_selected[1]
    ifelse(length(s),ifelse(stringr::str_count(db_hm$Identifiers[s],"\\|")>4,300+(stringr::str_count(db_hm$Identifiers[s],"\\|")-4)*40,300),300)}))
  
  getCalibrate = eventReactive(input$go,{
    if(length(input$calibrate)>0){
      
    }
    
    return(protein_list$uniprotswissprot[as.character(protein_list$external_gene_name) %in% as.character(input$calibrate)] )
  })
  
  
  getCellComponent <- function(s){
    go_cc <- gsub("\\s*\\([^\\)]+\\)","",strsplit(db_hm$Go.CC[s] ,"\\|")[[1]])
    protein_interest <- NA
    protein_filter <- NA
    if(!is.null(dim(gM1[rownames(gM1) %in% go_cc,]))){colsum <- Matrix::colSums(gM1[rownames(gM1) %in% go_cc,])
    protein_interest <- names(colsum[colsum!=0])
    protein_filter <- unique(c(protein_interest,uniprot_ubiquitin_H_sapiens$Entry))
    }
    if(length(gM1[rownames(gM1) %in% go_cc,])>0 & is.null(dim(gM1[rownames(gM1) %in% go_cc,])))
    {
      gM2 <- gM1[rownames(gM1) %in% go_cc,]
      protein_interest <- names(gM2==1)
      protein_filter <- unique(c(protein_interest,uniprot_ubiquitin_H_sapiens$Entry))
    }
    return(protein_filter)
  }
  
  RunRRW <- eventReactive(input$go, {
    s = input$table_rows_selected[1]
    if(length(s)){
      complex = unlist(strsplit(BuildDatabaseC()$Members[s], "[|]"))
      if(db_hm$cc[s]==1){
        interactors <- getCellComponent(s)
        AS <- A[rownames(A) %in% interactors, colnames(A) %in% interactors]
        AS1 <- Diagonal(dim(AS)[1])-0.3*AS
        ps1 <- p0[rownames(p0)%in% interactors, ]
        ps0 <- Matrix(ps1,ncol=1)
        rownames(ps0) <- names(ps1)
      }else{
        AS <- A
        AS1 <- Diagonal(dim(AS)[1])-0.3*AS
        ps0 <- p0
      }
      p1=ps0
      p1[rownames(p1) %in% unique(c(complex,getCalibrate())),] <- 1.0/length(unique(c(complex,getCalibrate())))
      b <- p1*0.7
      x0 <- p1
      run1  <- conjgrad(AS1,b,x0)
      result <- data.frame(Protein= rownames(run1),score = unname(run1[,1]))
      result[["geneName"]] <- protein_list$external_gene_name[match(result$Protein,protein_list$uniprotswissprot)]
      result$input <- ifelse(result$Protein %in% c(complex,getCalibrate()),1,0)
      
      result <- getProteinComplexMembership(result,Complex_protein,protein_list )
      result <- result[order(result$input,result$score),]
      result <- result[order(-result$score),]
      result <- result[!as.character(result$Protein) %in% as.character(complex),]
      
      Gos <- unlist(lapply(strsplit(as.character(db_hm[s,"Go.Annotations"]),"\\|")[[1]], function(x) gsub("\\(.*\\)","",x)))
      
      go1 <- unique(as.character(goa$GO[as.character(goa$GO) %in% Gos]))
      
      goa_sel <- goa[as.character(goa$GO) %in% Gos,]
      
      dim(goa_sel)
      length(unique(goa_sel$uniprotID))
      
      goa_sel <- goa %>% 
        dplyr::filter(GO %in% Gos)%>%
        dplyr::group_by(uniprotID) %>%
        dplyr::summarize(GO = paste0(GO,collapse="|"), 
                         InterPro.ID = paste0(unique(unlist(strsplit(as.character(ID), "\\|"))), collapse="|"), 
                         Description = paste0(unique(unlist(strsplit(as.character(Description), "\\,"))), collapse="|"),
                         GO.type = paste0(unique(GO.type) , collapse = "|"))
      
      if(dim(goa_sel)[1]>0){
        result <- cbind(result,goa_sel[match(as.character(result$Protein),goa_sel$uniprotID), c("GO","InterPro.ID","Description", "GO.type")])
      }
      
      rownames(result) <- NULL
      
      return(result)
    }
  },ignoreNULL = FALSE)
  
  getResult <- eventReactive(input$go, {
    RunRRW()
  },ignoreNULL = FALSE)
  
  
  makeResultNetwork <- eventReactive(input$go,{
    s = input$table_rows_selected[1]
    if(length(s)){
      complex = unlist(strsplit(BuildDatabaseC()$Members[s], "[|]"))
      result = getResult()
      result <- result[!as.character(result$Protein) %in% as.character(complex),]
      nodes3 <- makeInitialGraph()[[1]]
      edges3 <- makeInitialGraph()[[2]]
      top_proteins <- unique(as.character(result$Protein[1:input$num]),getCalibrate())
      result_protein <- unique(c(top_proteins,complex))
      
      Pathway_layer <- M_Pathway.normalized[rownames(M_Pathway.normalized) %in% result_protein ,colnames(M_Pathway.normalized) %in% result_protein ]
      g_Pathway <- graph_from_adjacency_matrix(Pathway_layer, weighted = T, mode = c("upper"))
      edge_Pathway <- data.frame(cbind(get.edgelist(g_Pathway) , round(E(g_Pathway)$weight,3)))
      colnames(edge_Pathway) <- c("from","to","value")
      edge_Pathway[["type"]] <- "Pathway"
      
      PPI_layer <- PPI_layer[rownames(PPI_layer) %in% result_protein ,colnames(PPI_layer) %in% result_protein ]
      g_PPI <- graph_from_adjacency_matrix(PPI_layer, weighted = T, mode = c("upper"))
      edge_PPI <- data.frame(cbind(get.edgelist(g_PPI) , round( E(g_PPI)$weight,3)))
      colnames(edge_PPI) <- c("from","to","value")
      edge_PPI[["type"]] <- "Physical Interaction"
      
      XL_layer <- M_xl_graph_corrected.normalized[rownames(M_xl_graph_corrected.normalized) %in% result_protein, colnames(M_xl_graph_corrected.normalized) %in% result_protein]
      g_XL_layer <- graph_from_adjacency_matrix(XL_layer, mode = c("upper"))
      edge_XL <- data.frame(cbind(get.edgelist(g_XL_layer)))
      colnames(edge_XL) <- c("from","to")
      
      if(dim(edge_XL)[1]>0){
        edge_XL[["value"]] <- 0.1
        edge_XL[["type"]] <- "XL"
        result_edges <- rbind(edge_Pathway, 
                              edge_PPI, 
                              edge_XL)
      }else{
        result_edges <- rbind(edge_Pathway, 
                              edge_PPI)
      }
      result_edges <- result_edges %>% dplyr::filter(!(from %in% complex & to %in% complex))
      edges3[["value"]] <- 0.2
      edges3[["type"]] <- "ComplexPortal"
      result_geneName <- protein_list$external_gene_name[match(top_proteins, protein_list$uniprotswissprot)]
      result_nodesDS <- data.frame(top_proteins, result_geneName)
      colnames(result_nodesDS) <- c("result_protein","result_geneName")
      result_uniprot <- lapply(as.character(result_nodesDS$result_protein), function(x) strsplit(x,"-")[[1]][1])
      
      top_interpro <- as.character(result$InterPro.ID[1:input$num])
      
      result_nodes <- data.frame(id = top_proteins,
                                 label = ifelse(!is.na(result_nodesDS$result_geneName),as.character(result_nodesDS$result_geneName),as.character(result_nodesDS$result_protein)),
                                 shape = ifelse( result_nodesDS$result_protein %in% uniprot_ubiquitin_H_sapiens[["Entry"]],"triangle","dot"), 
                                 color = list(background = ifelse(is.na(top_interpro),"#FF9966","red"), 
                                              border = "black",
                                              highlight = "yellow"),
                                 font = list(size=22, bold= TRUE),
                                 shadow = list(enabled = TRUE, size = 12)) 
      nodes3$font.size <- 22
      nodes3$size <- 25
      result1 <- result[as.character(result$Protein) %in% as.character(result_nodes$id),]
      diff <- -log(min(result1[["score"]][!result1$Protein%in% c(getCalibrate(),Complex)])) + log(max(result1[["score"]][!result1$Protein %in% c(getCalibrate(),Complex)]))
      result_nodes$size<- 12+13*(log(result1[["score"]])-log(min(result1[["score"]][!result1$Protein %in% c(getCalibrate(),Complex)])))/diff
      nodes <- rbind(nodes3,result_nodes)
      nodes[as.character(nodes$id) %in% getCalibrate(), "color.background"] <- "linen"
      nodes$size[as.character(nodes$id) %in% getCalibrate()] <- 25
      
      edges <- rbind(edges3,result_edges)
      save(nodes, edges,file="test.RData")
      return(list(nodes, edges))
    }
  })
  
  
  output$resultNetwork = renderVisNetwork({
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    s = input$table_rows_selected[1]
    if(length(s)){
      result_nodes_edges <- makeResultNetwork()
      nodes <- result_nodes_edges[[1]]
      edges <- result_nodes_edges[[2]]
      edges$smooth <- FALSE
      edges$dashes <- ifelse(edges$type=="Pathway",TRUE,FALSE)
      edges$color <- ifelse(edges$type=="ComplexPortal", "#C1CDCD",ifelse(edges$type=="Pathway", "#CD5B45", ifelse(edges$type=="Physical Interaction","#458B74","deepskyblue")))
      
      lnodes <- data.frame(label = c("Existing members","Closest interactors","interactors w domain" ), 
                           color.background = c("linen","#FF9966", "red"))
      
 
      # edges data.frame for legend
      ledges <- data.frame(color = c("#C1CDCD", "#CD5B45","green","deepskyblue"),
                           label = c("ComplexPortal","Pathway","Physical Interaction","XL"))
      
      visNetwork(nodes, edges, height = "400px", width = "100%",
                 main = list(text = paste0(BuildDatabaseC()$Name[s]), style = "color:black;font-size:18px;text-align:center;"))%>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE )%>%
        visLayout(randomSeed = 123) %>%
        visPhysics(solver = "hierarchicalRepulsion", maxVelocity = 25, minVelocity = 20)%>%
        visLegend(addEdges = ledges, addNodes = lnodes, width = 0.1, position ="left", useGroups = FALSE)%>%
        visInteraction(navigationButtons = TRUE)
    }
  })
  

  returnInteractions <- reactive({
    s = input$table_rows_selected[1]
    complex = unlist(strsplit(BuildDatabaseC()$Members[s], "[|]"))
    result = getResult()
    top_proteins <- as.character(result$Protein[1:input$num])
    result_protein <- unique(c(top_proteins,complex))
    dt <- directInteraction %>% 
      dplyr::filter(grepl(paste0(result_protein,collapse="|"),as.character(A))|grepl(paste0(result_protein,collapse="|"),as.character(B))) %>% 
      dplyr::select(A, B, A.name, B.name, pubmed, firstAuthor, confidence, provider, interactionID, imex, doi, mint) %>%
      dplyr::filter(as.character(A) != as.character(B))

    return(data.frame(dt))
  })
  
  output$table3 <- DT::renderDataTable(DT::datatable({
    dt <- returnInteractions() 

    
    return(dt)
  }), selection = 'single')
  
  
  output$plot <- renderPlotly({
    result = getResult()
    result <- result[result$input== 0,]
    result <- result[order(-result$score),]
    #result$Protein <- factor(result$Protein, levels = result$Score)
    res = result[1:input$num,]
    write.csv(res,"res.csv")
    res$geneName <- factor(res$geneName,levels=res$geneName[order(-res$score)])
    plot_ly(res) %>%
      add_trace(y = ~geneName, x = ~score, name = "Proteins score before calibration",type = "bar",orientation = 'h',marker = list(color = ~ifelse(count>0,"#00CC99","#FF9966")),
                hoverinfo = "text",
                text = ~paste(geneName,ifelse(Membership!="", paste("is Member of ",Membership),paste(" is not a member of known complex")))) %>%
      layout(title = "Proteins score", 
             yaxis = list(autorange="reversed"),
             xaxis = list(side = 'left',  showgrid = FALSE, zeroline = FALSE)
             #,
             #xaxis2 = list(side = 'right', overlaying = "x", showgrid = FALSE, zeroline = FALSE)
      )
    
  })
  
  
  #  observe({ 
  #    dt <- directInteraction %>% 
  #      dplyr::filter(grepl(paste0(uniprot,collapse="|"),as.character(A))|grepl(paste0(uniprot,collapse="|"),as.character(B))) %>% 
  #      dplyr::filter(as.character(A.name) != as.character(B.name))
  #    s3 = input$table3_rows_selected[1]
  #    url <-  a(href=paste0('https://www.ncbi.nlm.nih.gov/pubmed/?term=', dt$pubmed[s3]),dt$pubmed[s3],target="_parent")
  #  })
  
  output$frame <- renderUI({
    dt <- returnInteractions()
    s3 = input$table3_rows_selected[1]
    url <-  paste0('https://www.ncbi.nlm.nih.gov/pubmed/?term=',unname(dt$pubmed[s3]))
    my_test <- tags$iframe(src=url, height=450, width=1100)
    my_test
  })
  
  
  # Generate a summary of the data ----
  output$summary <- renderPrint({
    summary(d())
  })
  
  # Download file
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      write.csv(getResult(), file)
    }
  )
  
}