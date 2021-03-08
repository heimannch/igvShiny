library(shiny)
library(igvShiny)
library(dplyr)
library(GenomicAlignments)
#----------------------------------------------------------------------------------------------------
# we need a local directory to write files - for instance, a vcf file representing a genomic
# region of interest.  we then tell shiny about that directory, so that shiny's built-in http server
# can serve up files we write there, ultimately consumed by igv.js
if(!dir.exists("tracks"))
  dir.create("tracks")
addResourcePath("tracks", "tracks")
#----------------------------------------------------------------------------------------------------
f <- system.file(package="igvShiny", "extdata", "gwas.RData")
stopifnot(file.exists(f))
# tbl.gwas <- iatlas.api.client::query_germline_gwas_results(datasets = "TCGA") %>% 
#   dplyr::select(SNPS = snp_rsid, cptid = snp_name, CHR_ID = snp_chr, CHR_POS = snp_bp, 'P.VALUE'= p_value, feature_display)
ns.sep <- "."
#----------------------------------------------------------------------------------------------------
igv_ui = function(id){
  
  ns <- NS(id)
  
  shinyUI(fluidPage(
    
    sidebarLayout(
      sidebarPanel(
        actionButton(ns("searchButton"), "Search"),
        textInput(ns("roi"), label=""),
        shiny::uiOutput(ns("features")),
        shiny::radioButtons(ns("feature_action"), "Select or Exclude?", choices = c("Select", "Exclude"), selected = "Exclude"),
        shiny::uiOutput(ns("search_snp")),
        actionButton(ns("addGwasTrackButton"), "Add GWAS Track"),
        div(style="background-color: white; width: 200px; height:30px; padding-left: 5px;
            margin-top: 10px; border: 1px solid blue;",
            htmlOutput(ns("chromLocDisplay"))),
        hr(),
        width=2
      ),
      mainPanel(
        igvShinyOutput(ns('igvShiny_0')),
        width=10
      )
      ) # sidebarLayout
  ))
}
#----------------------------------------------------------------------------------------------------
igv_server <-  function(input, output, session) {
  
  ns <- session$ns
  
  tbl.gwas <- reactive({
    iatlas.api.client::query_germline_gwas_results(datasets = "TCGA") %>% 
    dplyr::select(SNPS = snp_rsid, cptid = snp_name, CHR_ID = snp_chr, CHR_POS = snp_bp, 'P.VALUE'= p_value, feature_display)
  })
  
  #generating list of features
  immune_feat <- reactive({
    # tbl.gwas %>%
    #   dplyr::select(feature_display, category) %>%
    #   dplyr::group_by(category) %>%
    #   tidyr::nest(data = c(feature_display))%>%
    #   dplyr::mutate(data = purrr::map(data, tibble::deframe)) %>%
    #   tibble::deframe()
    tbl.gwas()$feature_display
  })
  
  output$features <- renderUI({
    shiny::selectizeInput(ns('immunefeature'), "Select Immune Feature(s)",
                          choices = immune_feat(),
                          selected =c("MHC2 21978456", "Th17 cells"),
                          multiple = TRUE)
  })
  
  #updating GWAS table with selections or exclusion of features
  gwas_df <- reactive({
    req(input$immunefeature)
    if(input$feature_action == "Exclude") tbl.gwas() %>% dplyr::filter(!(feature_display %in% input$immunefeature))
    else tbl.gwas() %>% dplyr::filter(feature_display %in% input$immunefeature) 
  })
  
  #generating the SNP id list
  output$search_snp <- renderUI({
    shiny::req(gwas_df())
    snp_options <- (gwas_df() %>% dplyr::filter(!is.na(SNPS)))$SNPS
    shiny::selectInput(ns("snp_int"), "Click on the plot or search for a SNP id:",
                       choices = c("", snp_options))
  })
  
  observeEvent(input$searchButton, {
    printf("--- search")
    searchString = isolate(input$roi)
    if(nchar(searchString) > 0)
      showGenomicRegion(session, id=session$ns("igvShiny_0"), searchString)
  })
  observeEvent(input$trackClick, {
    printf("--- trackclick event")
    x <- input$trackClick
    print(x)
  })
  # observeEvent(input[["igv-trackClick"]], {
  #   printf("--- igv-trackClick event")
  #   x <- input[["igv-trackClick"]]
  #   print(x)
  # })
  
  observeEvent(input$trackClick, {
    printf("--- igv-trackClick popup")
    x <- input$trackClick
 
    attribute.name.positions <- grep("name", names(x))
    attribute.value.positions <- grep("value", names(x))
    attribute.names <- as.character(x)[attribute.name.positions]
    attribute.values <- as.character(x)[attribute.value.positions]
    tbl <- data.frame(name=attribute.names,
                      value=attribute.values,
                      stringsAsFactors=FALSE)
    dialogContent <- renderTable(tbl)
    html <- HTML(dialogContent())
    showModal(modalDialog(html, easyClose=TRUE))
  })
  
  shiny::observeEvent(input$igvReady, {
    containerID <- input$igvReady
    printf("igv ready, %s", containerID)
    loadGwasTrack(session, id=session$ns("igvShiny_0"), trackName="GWAS", tbl=gwas_df(), deleteTracksOfSameName=FALSE)
  })
  
  observeEvent(input$addGwasTrackButton, {
    print(input$igvReadyTest)
    #showGenomicRegion(session, id=session$ns("igvShiny_0"), "chr19:45,248,108-45,564,645")
    loadGwasTrack(session, id=session$ns("igvShiny_0"), trackName="GWAS", tbl=gwas_df(),ymin = 5, ymax = 1+max(-log10(gwas_df()$P.VALUE)), deleteTracksOfSameName=FALSE)
  })
  observeEvent(input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]], {
    newLoc <- input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]]
    #observeEvent(input$genomicRegionChanged, {
    #newLoc <- input$genomicRegionChanged
    printf("new chromLocString: %s", newLoc)
    output$chromLocDisplay <- renderText({newLoc})
  })
  
  genomes <- c("hg38", "hg19", "mm10", "tair10", "rhos")
  loci <- c("chr5:88,466,402-89,135,305", "MEF2C", "Mef2c", "1:7,432,931-7,440,395", "NC_007494.2:370,757-378,078")
  i <- 2
  
  output$igvShiny_0 <- renderIgvShiny({
    igvShiny(list(
      genomeName="hg19",
      initialLocus="all", #loci[i],
      displayMode="SQUISHED"
    ))
  })
  
  
  
}

server <- function(input, output, session){
  callModule(igv_server, "igv")
}

ui <- fluidPage(
  igv_ui(id="igv")
)

runApp(shinyApp(ui = ui, server = server), port=9833)