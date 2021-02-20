library(shiny)
library(igvShiny)
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
tbl.gwas <- tbl.gwas <- iatlas.api.client::query_germline_gwas_results(datasets = "TCGA") %>% 
  dplyr::select(SNPS = snp_rsid, cptid = snp_name, CHR_ID = snp_chr, CHR_POS = snp_bp, 'P.VALUE'= p_value, feature_display)
ns.sep <- "."
#----------------------------------------------------------------------------------------------------
igv_ui = function(id){
  
  ns <- NS(id)
  printf("namespace: '%s'", ns("foo"))
  
  shinyUI(fluidPage(
    
    sidebarLayout(
      sidebarPanel(
        actionButton(ns("searchButton"), "Search"),
        textInput(ns("roi"), label=""),
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
  
  observeEvent(input$searchButton, {
    printf("--- search")
    searchString = isolate(input$roi)
    if(nchar(searchString) > 0)
      showGenomicRegion(session, id=session$ns("igvShiny_0"), searchString)
  })
  observeEvent(input$addGwasTrackButton, {
    printf("---- addGWASTrack")
    printf("current working directory: %s", getwd())
    #showGenomicRegion(session, id=session$ns("igvShiny_0"), "chr19:45,248,108-45,564,645")
    loadGwasTrack(session, id=session$ns("igvShiny_0"), trackName="GWAS", tbl=tbl.gwas,ymin = 6, ymax = 15, deleteTracksOfSameName=FALSE)
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
  
  output$igvShiny_0 <- renderIgvShiny(
    igvShiny(list(
      genomeName=genomes[i],
      initialLocus="all", #loci[i],
      displayMode="SQUISHED"
    ))
  )
}

server <- function(input, output, session){
  callModule(igv_server, "igv")
}

ui <- fluidPage(
  igv_ui(id="igv")
)

runApp(shinyApp(ui = ui, server = server), port=9833)