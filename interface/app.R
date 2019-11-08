
CRAN <- c("shiny","shinyjs","shinythemes","shinycssloaders","rintrojs","ggplot2","stringr",
                   "Matrix","igraph","visNetwork","DT","dplyr", "magrittr","plotly","reshape2")
bioconductor <- c("drawProteins")

lapply(c(CRAN,bioconductor), library, character.only = TRUE)
#setRepositories(ind = 1,addURLs = c(BioC = "https://bioconductor.org/packages/3.8/bioc"))
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path)))

source("server.R", local=TRUE)
source("ui.R", local=TRUE)

#############################################################################################################################
# Create Shiny app ----
shinyApp(ui, server)