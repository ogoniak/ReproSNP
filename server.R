#----------------------------------------------------------------------------------------------------
# SERVER
#----------------------------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  ###### use coordinate string as display field
  not_edit <- TRUE
  if (not_edit) {
    toggleState("display_coordinate_string")
  }
  
  ###### Reactive Values
  R <- reactiveValues()
  
  ###### SERVER-SIDE SELECTION
  updateSelectizeInput(
    session,
    inputId = "gene_symbol",
    label = "Gene symbol (RefSeq)",
    choices = refseq_name2,
    selected = "LITAF",
    server = TRUE
  )
  
  ###### INPUTS
  observeEvent(input$includeAF0, {
    R$includeAF0 <- input$includeAF0
  })
  
  observeEvent(input$input_genomic_version, {
    R$input_genomic_version <- input$input_genomic_version
    if (R$input_genomic_version == "hg19") {
      R$coordinate_string <- liftover_function(coordinates = R$coordinate_string, liftover = "hg19ToHg38")
      updateTextInput(session, "display_coordinate_string", value = R$coordinate_string)
    } else if (R$input_genomic_version == "hg38") {
      R$coordinate_string <- liftover_function(coordinates = R$coordinate_string, liftover = "hg38ToHg19")
      updateTextInput(session, "display_coordinate_string", value = R$coordinate_string)
    }
  })
  
  observeEvent(input$gene_symbol, {
    R$gene_symbol <- input$gene_symbol
    R$RefSeq_subset <- RefSeq[RefSeq$name2 == R$gene_symbol, ]
    R$chrom <- gsub(pattern = "chr", replacement = "", R$RefSeq_subset$chrom[1])
    R$txStart <- min(R$RefSeq_subset$txStart)
    R$txEnd <- max(R$RefSeq_subset$txEnd)
    R$gene_coordinates <- paste0(R$chrom, ":", R$txStart, "-", R$txEnd)
    R$coordinate_string <- R$gene_coordinates
    updateTextInput(session, "display_coordinate_string", value = R$coordinate_string)
    R$exonStarts <- as.numeric(unlist(strsplit(R$RefSeq_subset$exonStarts, split = ",")))
    R$exonEnds <- as.numeric(unlist(strsplit(R$RefSeq_subset$exonEnds, split = ",")))
  })
  
  output$gene_symbol <- renderText(R$gene_symbol)
  output$coordinate_string <- renderText(R$coordinate_string)
  
  output$chrom <- renderText({ R$chrom })
  output$start <- renderText({ R$start })
  output$end <- renderText({ R$end })
  
  
  ###### SUBMIT CALCULATION
  observeEvent(input$submit_button, {
    
    R$chrom <- unlist(str_split(string = R$coordinate_string, pattern = ":"))[1]
    R$coordinates <- unlist(str_split(string = R$coordinate_string, pattern = ":"))[2]
    R$start <- as.numeric(unlist(str_split(string = R$coordinates, pattern = "-"))[1])
    R$end <- as.numeric(unlist(str_split(string = R$coordinates, pattern = "-"))[2])
    
    ### delete temporary files
    commandline("rm /path/to/TMP/temporary*")
    
    ### FILTER COHORTS
    filter_cohort(
      phenotable = phenotable, dir = tmp_dir, cohort_name = "1",
      cohort_Sex = input$cohort1_Sex, cohort_Generation = input$cohort1_Generation,
      cohort_Number_of_children = input$cohort1_Number_of_children,
      cohort_Reported_infertility = input$cohort1_Reported_infertility,
      cohort_include_NA = input$cohort1_include_NA
    )
    filter_cohort(
      phenotable = phenotable, dir = tmp_dir, cohort_name = "2",
      cohort_Sex = input$cohort2_Sex, cohort_Generation = input$cohort2_Generation,
      cohort_Number_of_children = input$cohort2_Number_of_children,
      cohort_Reported_infertility = input$cohort2_Reported_infertility,
      cohort_include_NA = input$cohort2_include_NA
    )
    
    if ((R$end - R$start) > 2500000 | (R$end - R$start) < 0) {#2305000
      shinyalert(
        "Invalid genomic range",
        paste0(
          "Please submit a range of minimum length 0 and maximum length 2500000. Submitted genomic range was: ",
          (R$end - R$start)
        ),
        type = "error"
      )
    } else {
      
      progressSweetAlert(session = session, id = "myprogress", title = "Starting Calculation ...", display_pct = F, value = 0, striped = T, status = "primary")
      
      Sys.sleep(0.5)
      updateProgressBar(session = session, id = "myprogress", value = 10, title = "Filtering Variants ...")
      
      ### RUN BACKEND TOOLS
      bgen2vcf(chr = R$chrom, start = R$start, end = R$end, cohorts = 2)
      
      updateProgressBar(session = session, id = "myprogress", value = 40, title = "Creating Table ...")
      
      ### CREATE TABLE
      R$mergedTable <- mergeTables(tmp_dir = tmp_dir, n_cohorts = 2)
      R$mergedTable <- R$mergedTable[R$mergedTable$C1_OBS_CT != 0 | R$mergedTable$C2_OBS_CT != 0, ]  # filter
      R$mergedTable <- add_AF(R$mergedTable)
      R$mergedTable <- fisher(R$mergedTable)
      R$mergedTable <- R$mergedTable %>% mutate(dbSNP = apply(R$mergedTable["ID"], 1, link_database, url = dbsnp_url))
      R$mergedTable <- R$mergedTable %>% mutate(ClinVar = apply(R$mergedTable["ID"], 1, link_database, url = clinvar_url))
      R$mergedTable$POS <- sapply(strsplit(R$mergedTable$ID, split = ":"), function(x) x[[2]])
      
      if (R$includeAF0 == FALSE) {
        R$mergedTable <- R$mergedTable[R$mergedTable$C1_ALT_AF > 0 | R$mergedTable$C2_ALT_AF > 0, ]
      }
      
      updateProgressBar(session = session, id = "myprogress", value = 60, title = "Annotating ...")
      
      ### ANNOTATION
      R$dbsnp_annotation <- apply(R$mergedTable, 1, annotate_dbsnp, dbsnp = useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp"))
      R$dbsnp_annotation <- as.data.frame(R$dbsnp_annotation)
      colnames <- rownames(R$dbsnp_annotation)
      R$dbsnp_annotation <- transpose(R$dbsnp_annotation)
      names(R$dbsnp_annotation) <- colnames
      R$mergedTable <- cbind(R$mergedTable, R$dbsnp_annotation)
      
      ### RENAME AND REORDER
      new_names <- c(
        "CHROM", "ID", "REF", "ALT", "C1_ALT_count", "C1_count", "C2_ALT_count", "C2_count",
        "Cohort1 GT00", "Cohort1 GT01", "Cohort1 GT11", "C1_0", "C1_1", "C1_missing",
        "Cohort2 GT00", "Cohort2 GT01", "Cohort2 GT11", "C2_0", "C2_1", "C2_missing",
        "Cohort1 AAF", "Cohort2 AAF", "P_FisherValue", "P Fisher", "P Bonferroni", "dbSNP_Query", "ClinVar Query",
        "POS",
        "rsID", "dbSNP", "Clinical Significance", "dbSNP MAF", "dbSNP_MA"
      )
      colnames(R$mergedTable) <- new_names
      new_order <- c(
        "CHROM", "POS", "REF", "ALT", "dbSNP", "ClinVar Query", "Cohort1 AAF", "Cohort2 AAF",
        "dbSNP MAF", "Clinical Significance",
        "Cohort1 GT00", "Cohort1 GT01", "Cohort1 GT11", "Cohort2 GT00", "Cohort2 GT01", "Cohort2 GT11", "P Fisher", "P Bonferroni",
        "C1_0", "C1_1", "C1_missing", "C2_0", "C2_1", "C2_missing",  # 23
        "C1_ALT_count", "C1_count", "C2_ALT_count", "C2_count", "P_FisherValue",
        "rsID", "dbSNP_Query", "dbSNP_MA", "ID"
      )
      
      R$mergedTable <- R$mergedTable[, new_order]
      
      R$mergedTable[R$mergedTable == "NA"] <- NA
      
      updateProgressBar(session = session, id = "myprogress", value = 90, title = "Rendering DataTable ...")
      
      ### RENDER DATATABLE
      output$mergedTable <- renderDT({
        R$mergedTable
      },
      filter = "top",
      options = list(
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = list(list(extend = 'colvis', columns = 16:32)),
        columnDefs = list(list(visible = FALSE, targets = 16:32))
      ),
      server = FALSE,
      escape = FALSE,
      rownames = FALSE,
      extensions = c("Buttons"))
      
      # Downloads
      output$VariantTable_download <- downloadHandler(
        filename = "VariantTable.tsv",
        content = function(file) { write.table(R$mergedTable, file, row.names = FALSE, sep = "\t", quote = FALSE) }
      )
      
      updateProgressBar(session = session, id = "myprogress", value = 100)
      Sys.sleep(1)
      
      closeSweetAlert(session = session)
      sendSweetAlert(session = session, title = "Calculation completed", type = "success")
      
    }
    
  })
  
}
