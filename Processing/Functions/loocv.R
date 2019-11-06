### Mehrnoosh Oghbaie
### 03/27/2019

rm(list=ls())
## setting working directory
setwd(
  paste0(dirname(rstudioapi::getActiveDocumentContext()$path))
)

# install the required packages.
CRAN.packages <- c("readr","data.table","reshape2","dplyr","ggthemes","magrittr","igraph","sqldf","diffusr", "dnet","readxl","ggplot2","AnnotationDbi")
bioconductor.packages <- c("biomaRt")
source("../Functions/Functions.R")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

#########################################################################################################################
## making sure that protein IDs are swissprot
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

convertUniprot2Swissprot <- function(uniprot_gn) {
  mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  return(biomaRt::getBM(attributes=c('uniprot_gn', 'uniprotswissprot'), 
                        filters = 'uniprot_gn', 
                        values =as.character(uniprot_gn), 
                        mart = mart))
}



#save(M_Pathway.normalized,
#     M_protein_PPI_graph.normalized,	
#     M_coexpression_graph.normalized2,
#     M_colocalization_graph.normalized2,	
#     M_genetic_graph.normalized2,
#     M_gene_PPI_graph.normalized2, 
#     M_prediction_graph.normalized2, 
#     AdjacencyMatrix.Disease.normalized,	
#     Function.normalized, 
#     M_xl_graph_corrected.normalized,
#     Complex_protein_edege,
#     file="F:/Network_Analysis/Python3/backup3.RData")

load(file="F:/Network_Analysis/Python3/backup3.RData")

tuning <- list()
tuning[["PPI_Protein"]] <- M_protein_PPI_graph.normalized 
tuning[["XL"]] <- M_xl_graph_corrected.normalized
tuning[["Pathway"]] <- M_Pathway.normalized 
tuning[["Coex_G"]] <- M_coexpression_graph.normalized2
tuning[["Coloc_G"]] <- M_colocalization_graph.normalized2
tuning[["Genetic_G"]] <- M_genetic_graph.normalized2
tuning[["PPI_G"]] <- M_gene_PPI_graph.normalized2
tuning[["Predict_G"]] <-  M_prediction_graph.normalized2
tuning[["Disease"]] <- AdjacencyMatrix.Disease.normalized

###
View(db_hm)
complex_link <- db_hm$X.Complex.ac[db_hm$number.links>2]

########################################################################################################################################
### Build networks of complex protein edge

Total_list_homosapien <- read_csv("F:/ComplexPortal/Total_list_homosapien.csv")


CV <- list()
for(j in 1:10){
  df <- data.frame(matrix(NA, nrow=length(names(tuning)), ncol=4))
  colnames(df) <- c("layer", "denom","nom", "prim")
  
  for(i in 1: length(names(tuning))){
    name <- names(tuning)[i]
    X <- tuning[[name]]
    
    ##################################################################################
    ##### Randomly get 10% of the complexes out 
    Complex_protein_edege <- Matrix(0, nrow=dim(M_Pathway.normalized)[1], ncol=dim(M_Pathway.normalized)[2])
    colnames(Complex_protein_edege) <- colnames(M_Pathway.normalized)
    rownames(Complex_protein_edege) <- rownames(M_Pathway.normalized)
    LOO_complexID <- sample(complex_link, 16, replace = FALSE)
    CV[[paste0("LOO_Complex",j)]] <- LOO_complexID
    dk <- unique(Total_list_homosapien[!as.character(Total_list_homosapien$complexID) %in% LOO_complexID,c("uniprot_A", "uniprot_B")]) %>% filter(uniprot_A %in% colnames(Complex_protein_edege)& uniprot_B %in% colnames(Complex_protein_edege))
    
    for(k in 1:dim(dk)[1]){
      Complex_protein_edege[dk[["uniprot_A"]][k],dk[["uniprot_B"]][k]] <- 1
      Complex_protein_edege[dk[["uniprot_B"]][k],dk[["uniprot_A"]][k]] <- 1
    }
    
    link <- Complex_protein_edege
    om <- X
    size = dim(om)[1]*dim(om)[1]
    dim(om) <-  c(size,1)
    dim(link) <- c(size,1)
    nom <- t(om)%*%link
    denom <- crossprod(om)
    prim <- nom[1,1]/denom[1,1]
    df[i,"layer"] <- name
    df[i, "denom"] <- denom[1,1]
    df[i, "nom"] <- nom[1,1]
    df[i, "prim"] <- prim
  }
 
  df$params <- df$prim/sum(df$prim, na.rm=T)
  CV[[paste0("parameter.run",j)]] <- df
}

save(CV, file="10_fold_CV.RData")


##########################################################################
total_paramater <- data.frame(matrix(NA, nr=0,nc=6))
for( i in 1:10){
  df <- CV[[paste0("parameter.run",i)]]
  df[["run"]] <- paste0("run",i)
  total_paramater <- rbind(total_paramater,df)
}

p1 <- ggplot( total_paramater, aes( x = run, y = params, fill = layer ) ) + 
  geom_bar( stat = "identity", position = "stack" ) +
  coord_flip() +
  scale_y_reverse()+
  scale_fill_brewer( palette = "YlGnBu" ) +
  theme_minimal() + 
  theme( axis.title.y = element_blank(),
                           legend.position = "left" )


number_link <- list()
for(l in names(CV)[grepl("LOO_Complex",names(CV))]){
  Complex_protein_edege <- Matrix(0, nrow=dim(M_Pathway.normalized)[1], ncol=dim(M_Pathway.normalized)[2])
  colnames(Complex_protein_edege) <- colnames(M_Pathway.normalized)
  rownames(Complex_protein_edege) <- rownames(M_Pathway.normalized)
  
  dk <- unique(Total_list_homosapien[!as.character(Total_list_homosapien$complexID) %in% CV[[l]],c("uniprot_A", "uniprot_B")]) %>% filter(uniprot_A %in% colnames(Complex_protein_edege)& uniprot_B %in% colnames(Complex_protein_edege))
  
  for(k in 1:dim(dk)[1]){
    Complex_protein_edege[dk[["uniprot_A"]][k],dk[["uniprot_B"]][k]] <- 1
    Complex_protein_edege[dk[["uniprot_B"]][k],dk[["uniprot_A"]][k]] <- 1
  }
  number_link[[l]] <- sum(Complex_protein_edege)
}

names(number_link)
data <- data.frame(name = names(number_link),
                   number.link = unlist(unname(number_link)))

p2 <- ggplot(data, aes(x=name, y=number.link)) +
  geom_bar(stat="identity",fill="grey") +
  #scale_y_reverse()+
  coord_flip()+
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())
  
gridExtra::grid.arrange(p1, p2, nrow = 1,widths = c(3, 1))

##########################################################################
####Validating the first run

goa_human <- read.delim("F:/Network_Analysis/Mapping/goa/goa_human.txt", header=FALSE)


for(i in 2:10){
  result <- data.frame(matrix(NA, nr=0, nc= 17))
  colnames(result) <- c("X.Complex.ac", "Recommended.name","Identifiers" , "Complex.assembly", 
                        "Go.CC", "number.links","count", "Protein","num.existed.prot", 
                        "Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Rank8")
  
  
  params <- CV[[paste0("parameter.run",i)]]
  #db_hm <- db_hm %>% filter(!as.character("X.Complex.ac") %in% CV[[paste0("LOO_Complex",i)]])
  
  A1 <- M_protein_PPI_graph.normalized
  
  A2 <- M_Pathway.normalized
  
  A3 <- (params[params$layer == "PPI_Protein","prim"] * M_protein_PPI_graph.normalized + 
           params[params$layer == "XL","prim"] * M_xl_graph_corrected.normalized)/sum(params[params$layer %in% c("PPI_Protein","XL"),"prim"])
  
  A4 <- (params[params$layer == "PPI_Protein","prim"] * M_protein_PPI_graph.normalized + 
           params[params$layer == "Pathway","prim"] * M_Pathway.normalized)/sum(params[params$layer %in% c("PPI_Protein","Pathway"),"prim"])
  
  A5 <- (params[params$layer == "PPI_Protein","prim"] * M_protein_PPI_graph.normalized + 
           params[params$layer == "Pathway","prim"] * M_Pathway.normalized +
           params[params$layer == "XL","prim"] * M_xl_graph_corrected.normalized)/sum(params[params$layer %in% c("PPI_Protein","Pathway","XL"),"prim"])
  
  A6 <- (params[params$layer == "PPI_Protein","prim"] * M_protein_PPI_graph.normalized + 
           params[params$layer == "Pathway","prim"] * M_Pathway.normalized +
           params[params$layer == "XL","prim"] * M_xl_graph_corrected.normalized+
           params[params$layer == "Coex_G","prim"] * M_coexpression_graph.normalized2+
           params[params$layer == "Coloc_G","prim"] * M_colocalization_graph.normalized2+
           params[params$layer == "Genetic_G","prim"] * M_genetic_graph.normalized2+
           params[params$layer == "PPI_G","prim"] * M_gene_PPI_graph.normalized2+
           params[params$layer == "Predict_G","prim"] * M_prediction_graph.normalized2)/sum(params[params$layer %in% c("PPI_Protein","Pathway","XL","Coloc_G","Coex_G","Genetic_G","PPI_G","Predict_G"),"prim"])
  
  A7 <- (params[params$layer == "PPI_Protein","prim"] * M_protein_PPI_graph.normalized + 
           params[params$layer == "Pathway","prim"] * M_Pathway.normalized +
           params[params$layer == "XL","prim"] * M_xl_graph_corrected.normalized+
           params[params$layer == "Disease","prim"] * AdjacencyMatrix.Disease.normalized)/sum(params[params$layer %in% c("PPI_Protein","Pathway","XL","Disease"),"prim"])
  
  A8 <- (params[params$layer == "PPI_Protein","prim"] * M_protein_PPI_graph.normalized + 
           params[params$layer == "Pathway","prim"] * M_Pathway.normalized +
           params[params$layer == "XL","prim"] * M_xl_graph_corrected.normalized+
           params[params$layer == "Coex_G","prim"] * M_coexpression_graph.normalized2+
           params[params$layer == "Coloc_G","prim"] * M_colocalization_graph.normalized2+
           params[params$layer == "Genetic_G","prim"] * M_genetic_graph.normalized2+
           params[params$layer == "PPI_G","prim"] * M_gene_PPI_graph.normalized2+
           params[params$layer == "Predict_G","prim"] * M_prediction_graph.normalized2+
           params[params$layer == "Disease","prim"] * AdjacencyMatrix.Disease.normalized)/sum(params[params$layer %in% c("PPI_Protein","Pathway","XL","Coloc_G","Coex_G","Genetic_G","PPI_G","Predict_G","Disease"),"prim"])
  
  AS1<- Diagonal(dim(A1)[1])- 0.3*A1
  AS2<- Diagonal(dim(A2)[1])- 0.3*A2
  AS3<- Diagonal(dim(A3)[1])- 0.3*A3
  AS4<- Diagonal(dim(A4)[1])- 0.3*A4
  AS5<- Diagonal(dim(A5)[1])- 0.3*A5
  AS6<- Diagonal(dim(A6)[1])- 0.3*A6
  AS7<- Diagonal(dim(A7)[1])- 0.3*A7
  AS8<- Diagonal(dim(A8)[1])- 0.3*A8
  
  
  db_hm_sel <- db_hm %>% filter(count>2)
  
  for(j in 1:dim(db_hm_sel)[1]){
    
    ds <- db_hm_sel[j,c("X.Complex.ac", "Recommended.name","Identifiers" , "Complex.assembly", "Go.CC", "number.links","count")]
    db_hm_sel[j, "Identifiers"]
    complex_protein <- gsub("\\s*\\([^\\)]+\\)","",strsplit(as.character(db_hm_sel$Identifiers[j]) ,"\\|")[[1]])
    complex_protein <- complex_protein[complex_protein %in% protein_list2]
    
    if(length(complex_protein)>0){
      res <- data.frame(cbind(ds[rep(seq_len(nrow(ds)), each=length(complex_protein)),], 
                              Protein =complex_protein,
                              num.existed.prot = rep(length(complex_protein),length(complex_protein)),
                              Rank1 = rep(NA,length(complex_protein)),
                              Rank2 = rep(NA,length(complex_protein)),
                              Rank3 = rep(NA,length(complex_protein)),
                              Rank4 = rep(NA,length(complex_protein)),
                              Rank5 = rep(NA,length(complex_protein)),
                              Rank6 = rep(NA,length(complex_protein)),
                              Rank7 = rep(NA,length(complex_protein)),
                              Rank8 = rep(NA,length(complex_protein))
      ))
      
      selected_protein <- c()
      
      if(db_hm_sel[j,"cc"]==1){
        go_protein_human <- goa_human %>% mutate(uniprotID=V2, GO.CC = V5 )%>% 
          dplyr::select(uniprotID,GO.CC) %>% 
          dplyr::filter(GO.CC %in% gsub("\\s*\\([^\\)]+\\)","",strsplit(as.character(db_hm_sel$Go.CC[j]) ,"\\|")[[1]])) %>% 
          dplyr::select(uniprotID) %>% 
          dplyr::mutate(uniprotID = as.character(uniprotID)) %>%
          distinct() %>% .$uniprotID
        selected_protein <- unique(c(go_protein_human,uniprot_ubiquitin_H_sapiens[["Entry"]]))
      }
      
      if(length(selected_protein)>0){
        AD1 <- A1[rownames(A1) %in% selected_protein, colnames(A1) %in% selected_protein]
        AD2 <- A2[rownames(A2) %in% selected_protein, colnames(A2) %in% selected_protein]
        AD3 <- A3[rownames(A3) %in% selected_protein, colnames(A3) %in% selected_protein]
        AD4 <- A4[rownames(A4) %in% selected_protein, colnames(A4) %in% selected_protein]
        AD5 <- A5[rownames(A5) %in% selected_protein, colnames(A5) %in% selected_protein]
        AD6 <- A6[rownames(A6) %in% selected_protein, colnames(A6) %in% selected_protein]
        AD7 <- A7[rownames(A7) %in% selected_protein, colnames(A7) %in% selected_protein]
        AD8 <- A1[rownames(A8) %in% selected_protein, colnames(A8) %in% selected_protein]
        
        AD1<- Diagonal(dim(AD1)[1])- 0.3*AD1
        AD2<- Diagonal(dim(AD2)[1])- 0.3*AD2
        AD3<- Diagonal(dim(AD3)[1])- 0.3*AD3
        AD4<- Diagonal(dim(AD4)[1])- 0.3*AD4
        AD5<- Diagonal(dim(AD5)[1])- 0.3*AD5
        AD6<- Diagonal(dim(AD6)[1])- 0.3*AD6
        AD7<- Diagonal(dim(AD7)[1])- 0.3*AD7
        AD8<- Diagonal(dim(AD8)[1])- 0.3*AD8
      }
      for(y in 1:length(complex_protein)){
        
        if(length(selected_protein)==0){
          seed <- complex_protein[!complex_protein %in% complex_protein[y]]
          size <- dim(AS1)[1]
          p0 <- matrix(0,nrow=size, ncol=1)
          rownames(p0) <- colnames(AS1)
          p1=p0
          p1[rownames(p1) %in% seed,] <- 1.0/length(seed)
          b <- p1*0.7
          x0 <- p1
          run1  <- conjgrad(AS1,b,x0)
          run2  <- conjgrad(AS2,b,x0)
          run3  <- conjgrad(AS3,b,x0)
          run4  <- conjgrad(AS4,b,x0)
          run5  <- conjgrad(AS5,b,x0)
          run6  <- conjgrad(AS6,b,x0)
          run7  <- conjgrad(AS7,b,x0)
          run8  <- conjgrad(AS8,b,x0)}
        
        if(length(selected_protein)>0){
          seed <- complex_protein[!complex_protein %in% complex_protein[y]]
          size <- dim(AD1)[1]
          p0 <- matrix(0,nrow=size, ncol=1)
          rownames(p0) <- colnames(AD1)
          p1=p0
          p1[rownames(p1) %in% seed,] <- 1.0/length(seed)
          b <- p1*0.7
          x0 <- p1
          run1  <- conjgrad(AD1,b,x0)
          run2  <- conjgrad(AD2,b,x0)
          run3  <- conjgrad(AD3,b,x0)
          run4  <- conjgrad(AD4,b,x0)
          run5  <- conjgrad(AD5,b,x0)
          run6  <- conjgrad(AD6,b,x0)
          run7  <- conjgrad(AD7,b,x0)
          run8  <- conjgrad(AD8,b,x0)
        }
        
        result1 <- data.frame(rownames(run1),matrix(run1))
        result1 <- result1[order(-result1$matrix.run1.),]
        result1 <- result1[!result1$rownames.run1. %in% seed,]
        rownames(result1) <- NULL
        if(length(as.numeric(rownames(result1[result1$rownames.run1.==as.character(res$Protein[y]),])))>0){
          res$Rank1[y] <- as.numeric(rownames(result1[result1$rownames.run1.==as.character(res$Protein[y]),]))
        }
        result2 <- data.frame(rownames(run2),matrix(run2))
        result2 <- result2[order(-result2$matrix.run2.),]
        result2 <- result2[!result2$rownames.run2. %in% seed,]
        rownames(result2) <- NULL
        if(length(as.numeric(rownames(result2[result2$rownames.run2.==as.character(res$Protein[y]),])))>0){
          res$Rank2[y] <- as.numeric(rownames(result2[result2$rownames.run2.==as.character(res$Protein[y]),]))
        }
        result3 <- data.frame(rownames(run3),matrix(run3))
        result3 <- result3[order(-result3$matrix.run3.),]
        result3 <- result3[!result3$rownames.run3. %in% seed,]
        rownames(result3) <- NULL
        if(length(as.numeric(rownames(result3[result3$rownames.run3.==as.character(res$Protein[y]),])))>0){
          res$Rank3[y] <- as.numeric(rownames(result3[result3$rownames.run3.==as.character(res$Protein[y]),]))
        }
        
        result4 <- data.frame(rownames(run4),matrix(run4))
        result4 <- result4[order(-result4$matrix.run4.),]
        result4 <- result4[!result4$rownames.run4. %in% seed,]
        rownames(result4) <- NULL
        if(length(as.numeric(rownames(result4[result4$rownames.run4.==as.character(res$Protein[y]),])))>0){
          res$Rank4[y] <- as.numeric(rownames(result4[result4$rownames.run4.==as.character(res$Protein[y]),]))
        }
        
        result5 <- data.frame(rownames(run5),matrix(run5))
        result5 <- result5[order(-result5$matrix.run5.),]
        result5 <- result5[!result5$rownames.run5. %in% seed,]
        rownames(result5) <- NULL
        if(length(as.numeric(rownames(result5[result5$rownames.run5.==as.character(res$Protein[y]),])))>0){
          res$Rank5[y] <- as.numeric(rownames(result5[result5$rownames.run5.==as.character(res$Protein[y]),]))
        }
        
        result6 <- data.frame(rownames(run6),matrix(run6))
        result6 <- result6[order(-result6$matrix.run6.),]
        result6 <- result6[!result6$rownames.run6. %in% seed,]
        rownames(result6) <- NULL
        if(length(as.numeric(rownames(result6[result6$rownames.run6.==as.character(res$Protein[y]),])))>0){
          res$Rank6[y] <- as.numeric(rownames(result6[result6$rownames.run6.==as.character(res$Protein[y]),]))
        }
        
        result7 <- data.frame(rownames(run7),matrix(run7))
        result7 <- result7[order(-result7$matrix.run7.),]
        result7 <- result7[!result7$rownames.run7. %in% seed,]
        rownames(result7) <- NULL
        if(length(as.numeric(rownames(result7[result7$rownames.run7.==as.character(res$Protein[y]),])))>0){
          res$Rank7[y] <- as.numeric(rownames(result7[result7$rownames.run7.==as.character(res$Protein[y]),]))
        }
        
        result8 <- data.frame(rownames(run8),matrix(run8))
        result8 <- result8[order(-result8$matrix.run8.),]
        result8 <- result8[!result8$rownames.run8. %in% seed,]
        rownames(result8) <- NULL
        if(length(as.numeric(rownames(result8[result8$rownames.run8.==as.character(res$Protein[y]),])))>0){
          res$Rank8[y] <- as.numeric(rownames(result8[result8$rownames.run8.==as.character(res$Protein[y]),]))
        }
        
      }
      
    }
    result <- rbind(result,res)
  }
  
  CV[[paste0("result_run",i)]] <- data.frame(result)
  save(CV, file="10_fold_CV.RData")
}


load("10_fold_CV.RData")
##### Visualize each run rank

png("../CV/CV_image/CV.png", width=1200, height = 1500)
par(mfrow=c(2,1))
for(i in 1:10){
  dA <- CV[[paste0("result_run",i)]]
  plot(ecdf(dA$Rank1), main=paste0("LOOCV - run",i),cex=1.3,cex.main=2,cex.lab=1.3,
       col="gray", xlab="Rank of the Leave out Protein", ylab="")
  lines(ecdf(dA$Rank2), col="navy",cex=0.9)
  lines(ecdf(dA$Rank3), col="purple", cex=0.8)
  lines(ecdf(dA$Rank4), col="red", cex=1.4)
  lines(ecdf(dA$Rank5), col="blue", cex=0.7)
  lines(ecdf(dA$Rank6), col="green",cex=0.9)
  lines(ecdf(dA$Rank7), col="orange",cex=0.8)
  lines(ecdf(dA$Rank8), cex=0.7)
  legend("bottomright" , 
         title="Cumulative ranking using different layers",
         c("PPI", "Pathway", "PPI + XL", "PPI + Pathway", "PPI + Pathway + XL",
           "PPI + Pathway + XL + Gene(all)","PPI + Pathway + XL + Disease",
           "PPI + Pathway + XL + Gene + Disease"),
         col=c("gray", "navy", "purple", "red", "blue", "green", "orange", "black"),
         pch=19, cex=1)
}

dev.off()



png("../CV/CV_image/CV_complex.png", width=1500, height = 1000)
par(mfrow=c(2,2))
for(i in 1:10){
  dA <- CV[[paste0("result_run",i)]]
  LOO_Complex <- CV[[paste0("LOO_Complex",i)]]
  dAA <- dA %>% group_by(X.Complex.ac, Recommended.name, count,num.existed.prot) %>% summarize(Rank1= mean(Rank1, na.rm=T),
                                                                Rank2= mean(Rank2, na.rm=T),
                                                                Rank3= mean(Rank3, na.rm=T),
                                                                Rank4= mean(Rank4, na.rm=T),
                                                                Rank5= mean(Rank5, na.rm=T),
                                                                Rank6= mean(Rank6, na.rm=T),
                                                                Rank7= mean(Rank7, na.rm=T),
                                                                Rank8= mean(Rank8, na.rm=T))

  dAA1 <- dAA %>% filter(! X.Complex.ac %in% LOO_Complex)
  dAA2 <- dAA %>% filter(X.Complex.ac %in% LOO_Complex)
  plot(dAA1$count,dAA1$Rank6, pch=19, cex=1.3,cex.main=1.7,cex.lab=1.3,
       xlab="Size", ylab="Rank", main="Rank vs. sum of unique members of complex")
  points(dAA2$count,dAA2$Rank6, pch=19, col="red", cex=1.9) 
  legend("topright" , 
         title="training/test complex",
         c(paste("training rank ave = ",round(mean(dAA1$Rank6, na.rm=T),2)),
           paste0("testing rank ave = ",round(mean(dAA2$Rank6, na.rm=T),2))),
         col=c("black","red"),
         pch=19, cex=1.3)
  
  
  }

dev.off()


png("../CV/CV_image/CV_complex_log.png", width=1600, height = 1000)
par(mfrow=c(2,2))
for(i in 1:10){
  dA <- CV[[paste0("result_run",i)]]
  LOO_Complex <- CV[[paste0("LOO_Complex",i)]]
  dAA <- dA %>% group_by(X.Complex.ac, Recommended.name, count,num.existed.prot) %>% summarize(Rank1= mean(Rank1, na.rm=T),
                                                                                               Rank2= mean(Rank2, na.rm=T),
                                                                                               Rank3= mean(Rank3, na.rm=T),
                                                                                               Rank4= mean(Rank4, na.rm=T),
                                                                                               Rank5= mean(Rank5, na.rm=T),
                                                                                               Rank6= mean(Rank6, na.rm=T),
                                                                                               Rank7= mean(Rank7, na.rm=T),
                                                                                               Rank8= mean(Rank8, na.rm=T))
  
  dAA1 <- dAA %>% filter(! X.Complex.ac %in% LOO_Complex)
  dAA2 <- dAA %>% filter(X.Complex.ac %in% LOO_Complex)
  plot(dAA1$count,log(dAA1$Rank6), pch=19, cex=1.3,cex.main=1.7,cex.lab=1.3,
       xlab="Size", ylab="Rank", main="Log rank/Sum of unique members of complex")
  points(dAA2$count,log(dAA2$Rank6), pch=19, col="red", cex=1.9) 
  legend("topright" , 
         title="training and test complex",
         c(paste("training rank ave = ",round(mean(dAA1$Rank6, na.rm=T),2)),
           paste0("testing rank ave = ",round(mean(dAA2$Rank6, na.rm=T),2))),
         col=c("black","red"),
         pch=19, cex=1.3)
}

dev.off()
