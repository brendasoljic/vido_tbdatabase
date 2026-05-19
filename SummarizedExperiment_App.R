library(readxl)
library(dplyr)
library(SummarizedExperiment)
library(S4Vectors)
library(shiny)
library(DT)
library(pheatmap)
library(bslib)
library(plotly)

#new test line

# Helper to safely pick a column or fallback
safe_pick <- function(x, fallback = NA_character_) {
  if (length(x) == 0 || is.na(x)) fallback else x
}

load_experiment_se <- function(file, sheet, experiment_name,
                               species = "M. tuberculosis H37Rv",
                               omics_type = "Transcriptomics") {
  
  df <- read_excel(file, sheet = sheet)
  
  # Clean column names
  names(df) <- trimws(names(df))
  names(df) <- gsub("\\s+", "_", names(df))
  
  # ---- Detect gene columns ----
  gene_id_col   <- names(df)[grepl("Gene$|Rv|Primary_Target", names(df), ignore.case = TRUE)][1]
  gene_name_col <- names(df)[grepl("Gene_Name|Common_Name|Gene_Symbol", names(df), ignore.case = TRUE)][1]
  
  gene_id_col   <- safe_pick(gene_id_col, fallback = names(df)[1])
  gene_name_col <- safe_pick(gene_name_col, fallback = gene_id_col)
  
  # ---- MANUAL ASSAY COLUMN DETECTION ----
  if (experiment_name == "Shee Transcriptome") {
    assay_cols <- c("moxi2x_GSM1829746", "moxi4x_GSM1829747", "moxi8x_GSM1829748")
    assay_cols <- assay_cols[assay_cols %in% names(df)]  # keep only existing
  } else if (experiment_name == "Schubert Proteome") {
    assay_cols <- names(df)[grepl("^expr_", names(df))]
  } else {
    # Vilcheze uses log2 fold change columns
    assay_cols <- names(df)[grepl("log2", names(df), ignore.case = TRUE)]
  }
  
  if (length(assay_cols) == 0) {
    message(paste("No assay columns detected for", experiment_name))
  }
  
  # ---- Build colData ----
  if (length(assay_cols) == 0) {
    colData <- DataFrame(
      condition = character(0),
      experiment = character(0)
    )
  } else {
    colData <- DataFrame(
      condition = assay_cols,
      experiment = rep(experiment_name, length(assay_cols))
    )
  }
  rownames(colData) <- assay_cols
  
  # ---- Detect padj columns ----
  padj_cols <- names(df)[grepl("padj|FDR|adj", names(df), ignore.case = TRUE)]
  if (length(padj_cols) == 0) padj_cols <- character(0)
  
  sd_cols <- names(df)[grepl("^sd_", names(df), ignore.case = TRUE)]
  if (length(sd_cols) == 0) sd_cols <- character(0)
  
  # ---- Build rowData ----
  rowData <- DataFrame(
    gene_id   = df[[gene_id_col]],
    gene_name = df[[gene_name_col]]
  )
  
  missing <- is.na(rowData$gene_id) | rowData$gene_id == ""
  rowData$gene_id[missing] <- paste0("unknown_", seq_len(sum(missing)))
  
  dups <- duplicated(rowData$gene_id)
  if (any(dups)) {
    rowData$gene_id[dups] <- paste0(rowData$gene_id[dups], "_dup", seq_len(sum(dups)))
  }
  
  rownames(rowData) <- rowData$gene_id
  
  # ---- Build assay matrix ----
  if (length(assay_cols) > 0) {
    assay_mat <- as.matrix(df[, assay_cols])
    assay_mat <- suppressWarnings(apply(assay_mat, 2, as.numeric))
    rownames(assay_mat) <- rowData$gene_id
  } else {
    assay_mat <- matrix(nrow = nrow(rowData), ncol = 0)
    rownames(assay_mat) <- rowData$gene_id
    colnames(assay_mat) <- character(0)
  }
  
  # Precompute first_col safely
  if (length(assay_cols) > 0) {
    first_col <- assay_cols[1]
  } else {
    first_col <- NA_character_
  }
  
  # ---- Build DEG table ----
    
  deg_df <- df %>%
    transmute(
      gene_id   = .data[[gene_id_col]],
      gene_name = .data[[gene_name_col]],
      across(all_of(assay_cols), .names = "fc_{col}"),
      if (length(padj_cols) > 0) across(all_of(padj_cols), .names = "{col}"),
      if (length(sd_cols) > 0) across(all_of(sd_cols), .names = "{col}"),
      deg_direction = if (!is.na(first_col)) {
        case_when(
          .data[[first_col]] > 0 ~ "Up",
          .data[[first_col]] < 0 ~ "Down",
          TRUE ~ "No change"
        )
      } else {
        NA_character_
      },
      experiment = experiment_name
    )

  
  
  # ---- Experiment summary ----
  experiment_summary <- data.frame(
    experiment_name = experiment_name,
    omics_type = omics_type,
    species = species,
    n_genes = nrow(rowData),
    n_deg = nrow(deg_df),
    n_samples = ncol(assay_mat),
    stringsAsFactors = FALSE
  )
  
  # ---- Build SummarizedExperiment ----
  se <- SummarizedExperiment(
    assays = list(log2FC = assay_mat),
    rowData = rowData,
    colData = colData,
    metadata = list(
      deg_df = deg_df,
      summary = experiment_summary,
      sd_values = if (length(sd_cols) > 0) df[, sd_cols, drop = FALSE] else NULL
    )
  )
  
  return(se)
}

combine_se <- function(se_list) {
  
  all_genes <- Reduce(union, lapply(se_list, function(se) rownames(se)))
  
  master_rowData_raw <- rowData(se_list[[1]])
  master_rowData <- DataFrame(
    gene_id = all_genes,
    gene_name = master_rowData_raw$gene_name[match(all_genes, rownames(master_rowData_raw))]
  )
  rownames(master_rowData) <- master_rowData$gene_id
  
  aligned <- lapply(se_list, function(se) {
    current_assay <- assay(se)
    current_genes <- rownames(se)
    
    missing <- setdiff(all_genes, current_genes)
    
    if (length(missing) > 0) {
      empty_mat <- matrix(NA, nrow = length(missing), ncol = ncol(current_assay))
      rownames(empty_mat) <- missing
      colnames(empty_mat) <- colnames(current_assay)
      new_assay <- rbind(current_assay, empty_mat)
    } else {
      new_assay <- current_assay
    }
    
    new_assay <- new_assay[all_genes, , drop = FALSE]
    
    SummarizedExperiment(
      assays = list(log2FC = new_assay),
      rowData = DataFrame(dummy = rep(NA, length(all_genes))),
      colData = colData(se),
      metadata = metadata(se)
    )
  })
  
  combined <- do.call(cbind, aligned)
  rowData(combined) <- master_rowData
  
  combined_metadata <- list(
    deg_df = bind_rows(lapply(se_list, function(se) metadata(se)$deg_df)),
    summary = bind_rows(lapply(se_list, function(se) metadata(se)$summary))
  )
  
  metadata(combined) <- combined_metadata
  
  return(combined)
}

# ---- Load datasets ----
se_Vilcheze <- load_experiment_se(
  file = "data/Transcriptome_Vilcheze.xlsx",
  sheet = "Table S2c log2 fold change",
  experiment_name = "Vilcheze Transcriptome"
)

se_Shee <- load_experiment_se(
  file = "data/Transcriptome_Shee.xlsx",
  sheet = "Sheet1",
  experiment_name = "Shee Transcriptome"
)

se_Schubert <- load_experiment_se(
  file = "data/Proteome_Schubert_2015.xlsx",
  sheet = "Mtb absolute per condition",
  experiment_name = "Schubert Proteome",
  omics_type = "Proteomics"
)

se_combined <- combine_se(list(se_Vilcheze, se_Shee, se_Schubert))

# =========================
# Shiny App
# =========================

ui <- page_navbar(
  title = "MtB Experiment Dashboard",
  theme = bs_theme(bootswatch = "flatly"),
  
  nav_panel(
    "Overview",
    fluidPage(
      br(),
      
      fluidRow(
        column(
          width = 4,
          card(
            card_header("Total Experiments"),
            card_body(
              h2(
                nrow(
                  metadata(se_combined)$summary
                )
              )
            )
          )
        ),
        
        column(
          width = 4,
          card(
            card_header("Transcriptomic Experiments"),
            card_body(
              h2(
                sum(
                  metadata(se_combined)$summary$omics_type ==
                    "Transcriptomics"
                )
              )
            )
          )
        ),
        
        column(
          width = 4,
          card(
            card_header("Proteomic Experiments"),
            card_body(
              h2(
                sum(
                  metadata(se_combined)$summary$omics_type ==
                    "Proteomics"
                )
              )
            )
          )
        )
      ),
      
      br(),
      
      h3("Experiment Summaries"),
      
      h4("Transcriptomics Experiments"),
      DTOutput("transcript_table"),
      
      br(), br(),
      
      h4("Proteomics Experiments"),
      DTOutput("proteomics_table")
    )
  )
  ,
  
  nav_panel(
    "Gene Browser",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "selected_gene",
            "Search gene by ID:",
            choices = rownames(rowData(se_combined)),
            selected = rownames(rowData(se_combined))[1]
          )
        ),
        mainPanel(
          h3("Gene Information"),
          tableOutput("gene_info"),
          br(),
          
          h3("Transcriptomic Expression Values"),
          tableOutput("gene_expression_tx"),
          br(),
          
          h3("Proteomic Expression Values"),
          tableOutput("gene_expression_px")
        )
      )
    )
  ),
  
  nav_panel(
    "Heatmap",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "heatmap_gene",
            "Choose a gene:",
            choices = rownames(rowData(se_combined)),
            selected = rownames(rowData(se_combined))[1]
          ),
          numericInput("n_genes", "Nearby genes:", 20, min = 5, max = 100)
        ),
        mainPanel(
          plotlyOutput("heatmap_plot", height = "600px")
        )
      )
    )
  ),
  
  nav_panel(
    "DEGs",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "deg_filter",
            "Filter DEG direction:",
            choices = c("All", "Up", "Down")
          )
        ),
        mainPanel(
          DTOutput("deg_table")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$transcript_table <- renderDT({
    df <- metadata(se_combined)$summary
    df <- df[df$omics_type == "Transcriptomics", ]
    datatable(df, rownames = FALSE)
  })
  
  output$proteomics_table <- renderDT({
    df <- metadata(se_combined)$summary
    df <- df[df$omics_type == "Proteomics", ]
    datatable(df, rownames = FALSE)
  })
  
  output$experiment_table <- renderDT({
    datatable(metadata(se_combined)$summary, rownames = FALSE)
  })
  
  output$gene_info <- renderTable({
    rowData(se_combined)[input$selected_gene, , drop = FALSE]
  })
  
  output$gene_expression_tx <- renderTable({
    
    # Identify transcriptomic experiments
    tx_experiments <- metadata(se_combined)$summary$experiment_name[
      metadata(se_combined)$summary$omics_type == "Transcriptomics"
    ]
    
    tx_cols <- rownames(colData(se_combined))[colData(se_combined)$experiment %in% tx_experiments]
    
    if (length(tx_cols) == 0) {
      return(data.frame(Message = "No transcriptomic data available"))
    }
    
    conds <- tx_cols
    log2fc_vals <- as.numeric(assay(se_combined)[input$selected_gene, conds])
    
    deg_df <- metadata(se_combined)$deg_df
    
    # Compute padj values for transcriptomics
    pvals <- vapply(conds, function(cn) {
      
      exp_name <- as.character(colData(se_combined)[cn, "experiment"])
      suffix   <- sub("^log2_Fold_change_", "", cn)
      padj_col <- paste0("padj_", suffix)
      
      row_idx <- which(
        deg_df$gene_id == input$selected_gene &
          deg_df$experiment == exp_name
      )
      
      if (length(row_idx) == 1 && padj_col %in% colnames(deg_df)) {
        as.numeric(deg_df[row_idx, padj_col])
      } else {
        NA_real_
      }
    }, numeric(1))
    
    pvals_fmt <- ifelse(
      is.na(pvals),
      NA,
      formatC(pvals, format = "e", digits = 2)
    )
    
    data.frame(
      condition = conds,
      log2FC    = log2fc_vals,
      "P-Value" = pvals_fmt,
      check.names = FALSE
    )
  })
  
  output$gene_expression_px <- renderTable({
    
    # Identify proteomic experiments
    px_experiments <- metadata(se_combined)$summary$experiment_name[
      metadata(se_combined)$summary$omics_type == "Proteomics"
    ]
    
    px_cols <- rownames(colData(se_combined))[colData(se_combined)$experiment %in% px_experiments]
    
    if (length(px_cols) == 0) {
      return(data.frame(Message = "No proteomic data available"))
    }
    
    expr_vals <- as.numeric(assay(se_combined)[input$selected_gene, px_cols])
    
    data.frame(
      condition = px_cols,
      expression = expr_vals,
      check.names = FALSE
    )
  })
  
  
  output$heatmap_plot <- renderPlotly({
    
    gene_index <- match(input$heatmap_gene, rownames(assay(se_combined)))
    if (is.na(gene_index)) return(NULL)
    
    start <- max(1, gene_index - floor(input$n_genes / 2))
    end <- min(nrow(assay(se_combined)), gene_index + floor(input$n_genes / 2))
    
    # Keep only transcriptomic experiments
    transcript_cols <- rownames(colData(se_combined))[colData(se_combined)$experiment %in%
                                                        metadata(se_combined)$summary$experiment_name[
                                                          metadata(se_combined)$summary$omics_type == "Transcriptomics"
                                                        ]]
    
    mat <- assay(se_combined)[start:end, transcript_cols, drop = FALSE]
    
    
    df_long <- reshape2::melt(mat)
    colnames(df_long) <- c("Gene", "Condition", "Value")
    
    gene_info <- as.data.frame(rowData(se_combined))
    
    df_long$gene_name <- gene_info$gene_name[
      match(df_long$Gene, gene_info$gene_id)
    ]
    
    df_long$pmid <- "35975988"
    df_long$doi <- "10.1128/AAC.00592-22"
    df_long$publication <- 
      "Moxifloxacin-Mediated Killing of Mycobacterium tuberculosis Involves Respiratory Downshift, Reductive Stress, and Accumulation of Reactive Oxygen Species"
    
    # Compute p-values for transcriptomic experiments only
    deg_df <- metadata(se_combined)$deg_df
    
    df_long$p_value <- vapply(seq_len(nrow(df_long)), function(i) {
      gene <- df_long$Gene[i]
      cond <- df_long$Condition[i]
      
      exp_name <- as.character(colData(se_combined)[cond, "experiment"])
      
      # Identify the padj column for this condition
      suffix   <- sub("^log2_Fold_change_", "", cond)
      padj_col <- paste0("padj_", suffix)
      
      row_idx <- which(
        deg_df$gene_id == gene &
          deg_df$experiment == exp_name
      )
      
      if (length(row_idx) == 1 && padj_col %in% colnames(deg_df)) {
        as.numeric(deg_df[row_idx, padj_col])
      } else {
        NA_real_
      }
    }, numeric(1))
    
    df_long$p_value_fmt <- ifelse(
      is.na(df_long$p_value),
      "",
      paste0("P-value: ", formatC(df_long$p_value, format = "e", digits = 2))
    )
    
    df_long$hover_text <- paste0(
      "<b>", df_long$Gene, "</b>",
      "<br>", df_long$gene_name,
      "<br><br>",
      "Condition: ", df_long$Condition,
      "<br>log2FC: ", round(df_long$Value, 3),
      ifelse(df_long$p_value_fmt == "", "", paste0("<br>", df_long$p_value_fmt)),
      "<br><br>",
      "PMID: ", df_long$pmid,
      "<br>DOI: ", df_long$doi,
      "<br><br>",
      "Click for publication details"
    )
    
    
    p <- plot_ly(
      data = df_long,
      x = ~Condition,
      y = ~Gene,
      z = ~Value,
      type = "heatmap",
      colors = colorRamp(c("blue", "white", "red")),
      text = ~hover_text,
      hoverinfo = "text",
      source = "heatmap_click"
    )
    
    p <- event_register(p, "plotly_click")
    p
  })
  
  observeEvent(
    event_data("plotly_click", source = "heatmap_click"),
    {
      click <- event_data("plotly_click", source = "heatmap_click")
      
      if (!is.null(click)) {
        
        selected_gene <- click$y
        gene_info <- as.data.frame(rowData(se_combined))
        
        showModal(
          modalDialog(
            title = paste("Gene Details:", selected_gene),
            
            tags$strong("Gene Name:"),
            br(),
            gene_info$gene_name[
              match(selected_gene, gene_info$gene_id)
            ],
            
            br(), br(),
            
            tags$strong("Publication:"),
            br(),
            "Moxifloxacin-Mediated Killing of Mycobacterium tuberculosis Involves Respiratory Downshift, Reductive Stress, and Accumulation of Reactive Oxygen Species",
            
            br(), br(),
            
            tags$a(
              href = "https://pubmed.ncbi.nlm.nih.gov/35975988/",
              target = "_blank",
              "Open PubMed Page"
            ),
            
            br(), br(),
            
            tags$a(
              href = "https://doi.org/10.1128/AAC.00592-22",
              target = "_blank",
              "Open DOI Page"
            ),
            
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      }
    }
  )
  
  
  
  output$deg_table <- renderDT({
    df <- metadata(se_combined)$deg_df
    if (input$deg_filter != "All") {
      df <- df %>% filter(deg_direction == input$deg_filter)
    }
    datatable(df, rownames = FALSE)
  })
}

shinyApp(ui, server)