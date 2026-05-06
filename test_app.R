library(readxl)
library(dplyr)
library(SummarizedExperiment)
library(S4Vectors)
library(shiny)
library(DT)
library(pheatmap)
library(bslib)
library(plotly)


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
  
  # ---- ASSAY COLUMN DETECTION ----
  
  if (omics_type == "Transcriptomics") {
    
    if (experiment_name == "Shee Transcriptome") {
      assay_cols <- c("moxi2x_GSM1829746", "moxi4x_GSM1829747", "moxi8x_GSM1829748")
      assay_cols <- assay_cols[assay_cols %in% names(df)]
    } else {
      assay_cols <- names(df)[grepl("log2", names(df), ignore.case = TRUE)]
    }
    
  } else if (omics_type == "Proteomics") {
    
    # Expression columns renamed to expr_*
    assay_cols <- names(df)[grepl("^expr_", names(df))]
    
  } else {
    
    assay_cols <- character(0)
    
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
      across(all_of(padj_cols), .names = "padj_{col}"),
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
      summary = experiment_summary
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
  file = "C:/Users/mayac/OneDrive/Documents/vido_tbdatabase/Transcriptome_Vilcheze.xlsx",
  sheet = "Table S2c log2 fold change",
  experiment_name = "Vilcheze Transcriptome"
)

se_Shee <- load_experiment_se(
  file = "C:/Users/mayac/OneDrive/Documents/vido_tbdatabase/Transcriptome_Shee.xlsx",
  sheet = "Sheet1",
  experiment_name = "Shee Transcriptome"
)

# Combine transcriptomics only
se_transcriptome_combined <- combine_se(list(se_Vilcheze, se_Shee))

# Proteomics
se_Schubert <- load_experiment_se(
  file = "C:/Users/mayac/OneDrive/Documents/vido_tbdatabase/Proteome_Schubert.xlsx",
  omics_type = "Proteomics",
  sheet = "Mtb absolute per condition",
  experiment_name = "Schubert Proteome"
)

# se_proteome_combined <- combine_se(list(se_Proteomics1))

# For now, create an empty proteomics SE so the app runs
se_proteome_combined <- SummarizedExperiment(
  assays = list(log2FC = matrix(nrow = 0, ncol = 0)),
  rowData = DataFrame(),
  colData = DataFrame(),
  metadata = list(
    summary = data.frame(),
    deg_df = data.frame()
  )
)

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
      h3("Transcriptomics"),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "selected_gene_tx",
            "Search gene (Transcriptomics):",
            choices = rownames(rowData(se_transcriptome_combined)),
            selected = rownames(rowData(se_transcriptome_combined))[1]
          )
        ),
        mainPanel(
          h4("Gene Information"),
          tableOutput("gene_info_tx"),
          br(),
          h4("Expression Values"),
          tableOutput("gene_expression_tx")
        )
      ),
      
      br(), br(),
      
      h3("Proteomics"),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "selected_gene_px",
            "Search gene (Proteomics):",
            choices = rownames(rowData(se_proteome_combined)),
            selected = NULL
          )
        ),
        mainPanel(
          h4("Gene Information"),
          tableOutput("gene_info_px"),
          br(),
          h4("Expression Values"),
          tableOutput("gene_expression_px")
        )
      )
    )
  )
  ,
  
  nav_panel(
    "Heatmap",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "heatmap_gene",
            "Choose a gene:",
            choices = rownames(rowData(se_transcriptome_combined)),
            selected = rownames(rowData(se_transcriptome_combined))[1]
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
  
  output$gene_expression <- renderTable({
    data.frame(
      condition = colnames(assay(se_combined)),
      log2FC = as.numeric(assay(se_combined)[input$selected_gene, ])
    )
  })
  
  
  output$heatmap_plot <- renderPlotly({
    gene_index <- match(input$heatmap_gene, rownames(assay(se_transcriptome_combined)))
    if (is.na(gene_index)) return(NULL)
    
    start <- max(1, gene_index - floor(input$n_genes / 2))
    end <- min(nrow(assay(se_transcriptome_combined)), gene_index + floor(input$n_genes / 2))
    
    mat <- assay(se_transcriptome_combined)[start:end, , drop = FALSE]
    
    df_long <- reshape2::melt(mat)
    colnames(df_long) <- c("Gene", "Condition", "Value")
    
    plot_ly(
      data = df_long,
      x = ~Condition,
      y = ~Gene,
      z = ~Value,
      type = "heatmap",
      colors = colorRamp(c("blue", "white", "red")),
      hovertemplate = paste(
        "Gene: %{y}<br>",
        "Condition: %{x}<br>",
        "log2FC: %{z}<extra></extra>"
      )
    )
  })
  
  
  output$deg_table <- renderDT({
    df <- bind_rows(
      metadata(se_transcriptome_combined)$deg_df,
      metadata(se_proteome_combined)$deg_df
    )
    
    if (input$deg_filter != "All") {
      df <- df %>% filter(deg_direction == input$deg_filter)
    }
    
    datatable(df, rownames = FALSE)
  })
  
  # Transcriptomics
  output$gene_info_tx <- renderTable({
    rowData(se_transcriptome_combined)[input$selected_gene_tx, , drop = FALSE]
  })
  
  output$gene_expression_tx <- renderTable({
    data.frame(
      condition = colnames(assay(se_transcriptome_combined)),
      log2FC = as.numeric(assay(se_transcriptome_combined)[input$selected_gene_tx, ])
    )
  })
  
  # Proteomics
  output$gene_info_px <- renderTable({
    if (is.null(input$selected_gene_px)) return(NULL)
    rowData(se_proteome_combined)[input$selected_gene_px, , drop = FALSE]
  })
  
  output$gene_expression_px <- renderTable({
    if (is.null(input$selected_gene_px)) return(NULL)
    data.frame(
      condition = colnames(assay(se_proteome_combined)),
      log2FC = as.numeric(assay(se_proteome_combined)[input$selected_gene_px, ])
    )
  })
  
}

shinyApp(ui, server)
