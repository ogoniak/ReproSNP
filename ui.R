#----------------------------------------------------------------------------------------------------
# USER INTERFACE
#----------------------------------------------------------------------------------------------------
ui <- dashboardPage(title = "ReproSNP",
                    
  #----------------------------------------------------------------------------------------------------
  # EXTERNAL WEBLINKS IN HEADER
  #----------------------------------------------------------------------------------------------------
  dashboardHeader(
    title = span(img(src = "ReproSNP_Logo.svg", width = 200)),
    tags$li(a(href = "https://www.ukbiobank.ac.uk/",
              img(src = "UKB_Logo.svg",
                  title = "UK Biobank", height = "30px"),
              style = "padding-top:10px; padding-bottom:10px;"),
            class = "dropdown"),
    tags$li(a(href = "https://www.medizin.uni-muenster.de/busch-group/home.html",
              img(src = "Busch_Logo.svg",
                  title = "Busch Group", height = "30px"),
              style = "padding-top:10px; padding-bottom:10px;"),
            class = "dropdown"),
    tags$li(a(href = "https://www.medizin.uni-muenster.de/imi/das-institut.html",
              img(src = "IMI_Logo.svg",
                  title = "IMI", height = "30px"),
              style = "padding-top:10px; padding-bottom:10px;"),
            class = "dropdown")
  ),
  
  #----------------------------------------------------------------------------------------------------
  # SIDEBAR
  #----------------------------------------------------------------------------------------------------
  dashboardSidebar(
    use_theme(mytheme),
    
    sidebarMenu(
      
      useShinyjs(),
      br(),
      actionBttn(inputId = "submit_button", label = "Submit", style = "simple", class = "btn-success", width = "100%", size = "lg", no_outline = TRUE),
      br(),
      
      #----------------------------------------------------------------------------------------------------
      # SNP SELECTION
      #----------------------------------------------------------------------------------------------------
      menuItem("Locus", tabName = "locus", startExpanded = TRUE,
               selectizeInput(
                 inputId   = "gene_symbol",
                 label     = "Gene symbol (RefSeq)",
                 selected  = "LITAF",
                 choices   = NULL
               ),
               textInput("display_coordinate_string", label = "Genomic Range:", value = NULL),
               radioButtons(inputId = "input_genomic_version", label = "Genome version", choiceNames = c("GRCh38/hg38", "GRCh37/hg19"), choiceValues = c("hg38","hg19"), selected = "hg38"),
               br()
      ),
      
      menuItem("Variant Filter", tabname = "varfilter", startExpanded = FALSE,
               radioButtons(inputId = "includeAF0", label = "Variants with allele frequency = 0", choiceNames = c("Include", "Exclude"), choiceValues = c(TRUE, FALSE), selected = FALSE),
               br()
      ),
      
      menuItem("Cohort I", tabName = "cohort1", startExpanded = FALSE,
               checkboxGroupInput("cohort1_Sex", "Sex", unique(phenotable$Sex), selected = "Male"),
               checkboxGroupInput("cohort1_Generation", "Generation", c("Silent Generation", "Baby Boomer", "Generation X"), selected = unique(phenotable$Generation)),
               sliderInput(inputId = "cohort1_Number_of_children", label = "Number of children", min = 0, max = 3, step = 1, value = c(1, 3), ticks = FALSE),
               radioButtons(inputId = "cohort1_include_NA", label = "Reported number of children", choiceNames = c("Include NA", "Exclude NA"), choiceValues = c(-1, -2), selected = -2),
               checkboxGroupInput("cohort1_Reported_infertility", "Reported infertility", unique(phenotable$Reported_infertility), selected = FALSE),
               br()
      ),
      
      menuItem("Cohort II", tabName = "cohort2", startExpanded = FALSE,
               checkboxGroupInput("cohort2_Sex", "Sex", unique(phenotable$Sex), selected = "Male"),
               checkboxGroupInput("cohort2_Generation", "Generation", c("Silent Generation", "Baby Boomer", "Generation X"), selected = unique(phenotable$Generation)),
               sliderInput(inputId = "cohort2_Number_of_children", label = "Number of children", min = 0, max = 3, step = 1, value = c(0, 0), ticks = FALSE),
               radioButtons(inputId = "cohort2_include_NA", label = "Reported number of children", choiceNames = c("Include NA", "Exclude NA"), choiceValues = c(-1, -2), selected = -2),
               checkboxGroupInput("cohort2_Reported_infertility", "Reported infertility", unique(phenotable$Reported_infertility), selected = c(FALSE, TRUE)),
               br()
      )
    )
  ),
  
  #----------------------------------------------------------------------------------------------------
  # BODY WITH TABLES AND PLOTS
  #----------------------------------------------------------------------------------------------------
  dashboardBody(
    
    tags$head(
      tags$style(HTML("
                      /* center the notification */
                      .shiny-notification {
                        position: fixed;
                        top: 30%;
                        left: 50%;
                        transform: translate(-50%, -50%);
                        z-index: 9999;
                      }
                    "))
    ),
    
    fluidRow(style = "margin-top: 20px;margin-left: 20px; margin-right: 20px; margin-bottom: 20px;",
             tabsetPanel(
               
               tabPanel(h5("Variants"),
                        br(),
                        box(width = "100%",
                            dataTableOutput("mergedTable"),
                            downloadButton("VariantTable_download", "")
                        )
               )
               
             )
    )
  )
                    
)
