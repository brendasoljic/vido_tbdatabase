library(readxl)
library(dplyr)
library(SummarizedExperiment)
library(S4Vectors)
library(shiny)
library(DT)
library(pheatmap)
library(bslib)
library(plotly)
library(reshape2)

# Helper to avoid crashes with new datasets
safe_pick <- function(x, fallback = NA_character_) {
  if (length(x) == 0 || is.na(x)) fallback else x
}

# ---- Load experiments function ----
load_experiment_se <- function(file, sheet, experiment_name,
                               species = "M. tuberculosis H37Rv",
                               omics_type = "Transcriptomics",
                               pmid = NA_character_,
                               pmcid = NA_character_,
                               doi = NA_character_,
                               publication = NA_character_) {
  
  # ---- Reading columns correctly ----
  # Will need work as we get new datasets, and could be standardized
  
  df <- read_excel(file, sheet = sheet)
  
  names(df) <- trimws(names(df))
  names(df) <- gsub("\\s+", "_", names(df))
  
  gene_id_col   <- names(df)[grepl("Gene$|Rv|Primary_Target", names(df), ignore.case = TRUE)][1]
  gene_name_col <- names(df)[grepl("Gene_Name|Common_Name|Gene_Symbol", names(df), ignore.case = TRUE)][1]
  
  gene_id_col   <- safe_pick(gene_id_col, fallback = names(df)[1])
  gene_name_col <- safe_pick(gene_name_col, fallback = gene_id_col)
  
  if (experiment_name == "Shee Transcriptome") {
    
    assay_cols <- c("moxi2x_GSM1829746", "moxi4x_GSM1829747", "moxi8x_GSM1829748")
    assay_cols <- assay_cols[assay_cols %in% names(df)]
    
  } else if (experiment_name == "Schubert Proteome") {
    
    assay_cols <- names(df)[grepl("^expr_", names(df))]
    
  } else if (experiment_name == "Bei Transcriptome") {
    
    # Absolute transcriptome (log2 RPKM)
    assay_cols <- "Mean_logRPKM"
    
  } else {
    
    assay_cols <- names(df)[grepl("log2", names(df), ignore.case = TRUE)]
  }
  
  
  if (length(assay_cols) == 0) {
    message(paste("No assay columns detected for", experiment_name))
  }
  
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
  
  padj_cols <- names(df)[grepl("padj|FDR|adj", names(df), ignore.case = TRUE)]
  if (length(padj_cols) == 0) padj_cols <- character(0)
  
  sd_cols <- names(df)[grepl("^sd_", names(df), ignore.case = TRUE)]
  if (length(sd_cols) == 0) sd_cols <- character(0)
  
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
  
  if (length(assay_cols) > 0) {
    assay_mat <- as.matrix(df[, assay_cols])
    assay_mat <- suppressWarnings(apply(assay_mat, 2, as.numeric))
    rownames(assay_mat) <- rowData$gene_id
  } else {
    assay_mat <- matrix(nrow = nrow(rowData), ncol = 0)
    rownames(assay_mat) <- rowData$gene_id
    colnames(assay_mat) <- character(0)
  }
  
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
  
  if (length(padj_cols) > 0) {
    padj_df <- df %>%
      select(all_of(padj_cols))
    
    deg_df <- bind_cols(deg_df, padj_df)
  }
  
  if (length(sd_cols) > 0) {
    sd_df <- df %>%
      select(all_of(sd_cols))
    
    deg_df <- bind_cols(deg_df, sd_df)
  }
  
  # ---- Experiment summary ----
  
  experiment_summary <- data.frame(
    experiment_name = experiment_name,
    omics_type = omics_type,
    species = species,
    n_genes = nrow(rowData),
    n_deg = nrow(deg_df),
    n_samples = ncol(assay_mat),
    pmid = pmid,
    pmcid = pmcid,
    doi = doi,
    publication = publication,
    pubmed_link = if (!is.na(pmid) && pmid != "") {
      paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")
    } else {
      NA_character_
    },
    pmc_link = if (!is.na(pmcid) && pmcid != "") {
      paste0("https://pmc.ncbi.nlm.nih.gov/articles/", pmcid, "/")
    } else {
      NA_character_
    },
    doi_link = if (!is.na(doi) && doi != "") {
      paste0("https://doi.org/", doi)
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
  
  se <- SummarizedExperiment(
    assays = list(log2FC = assay_mat),
    rowData = rowData,
    colData = colData,
    metadata = list(
      deg_df = deg_df,
      summary = experiment_summary,
      sd_values = if (length(sd_cols) > 0) {
        df[, sd_cols, drop = FALSE]
      } else {
        NULL
      }
    )
  )
  
  return(se)
}

combine_se <- function(se_list) {
  
  all_genes <- Reduce(
    union,
    lapply(se_list, function(se) {
      rownames(se)
    })
  )
  
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
    deg_df = bind_rows(
      lapply(se_list, function(se) {
        metadata(se)$deg_df
      })
    ),
    summary = bind_rows(
      lapply(se_list, function(se) {
        metadata(se)$summary
      })
    )
  )
  
  metadata(combined) <- combined_metadata
  
  return(combined)
}

# ---- Load datasets ----

se_Vilcheze <- load_experiment_se(
  file = "data/Transcriptome_Vilcheze.xlsx",
  sheet = "Table S2c log2 fold change",
  experiment_name = "Vilcheze Transcriptome",
  omics_type = "Transcriptomics",
  doi = "10.3389/fimmu.2022.909904",
  pmcid = "PMC9283954",
  publication = "Transcriptional profiling of Mycobacterium tuberculosis"
)

se_Shee <- load_experiment_se(
  file = "data/Transcriptome_Shee.xlsx",
  sheet = "Sheet1",
  experiment_name = "Shee Transcriptome",
  omics_type = "Transcriptomics",
  pmid = "35975988",
  doi = "10.1128/AAC.00592-22",
  publication = "Moxifloxacin-Mediated Killing of Mycobacterium tuberculosis Involves Respiratory Downshift, Reductive Stress, and Accumulation of Reactive Oxygen Species"
)

se_Schubert <- load_experiment_se(
  file = "data/Proteome_Schubert_2015.xlsx",
  sheet = "Mtb absolute per condition",
  experiment_name = "Schubert Proteome",
  omics_type = "Proteomics",
  pmid = "26094805",
  doi = "10.1016/j.chom.2015.06.001",
  publication = "Schubert Proteome 2015"
)

se_Bei <- load_experiment_se(
  file = "data/Transcriptome_Bei_2024.xlsx",
  sheet = "log2FC of 5th and 95th",
  experiment_name = "Bei Transcriptome",
  omics_type = "Transcriptomics",
  pmid = "38600064",
  pmcid = "PMC11006872",
  doi = "10.1038/s41467-024-47410-5",
  publication = "Genetically encoded transcriptional plasticity underlies stress adaptation in Mycobacterium tuberculosis"
)


se_combined <- combine_se(list(se_Vilcheze, se_Shee, se_Schubert, se_Bei))

# ---- Shiny App ----

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
              h2(nrow(metadata(se_combined)$summary))
            )
          )
        ),
        column(
          width = 4,
          card(
            card_header("Transcriptomic Experiments"),
            card_body(
              h2(sum(metadata(se_combined)$summary$omics_type == "Transcriptomics"))
            )
          )
        ),
        column(
          width = 4,
          card(
            card_header("Proteomic Experiments"),
            card_body(
              h2(sum(metadata(se_combined)$summary$omics_type == "Proteomics"))
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
  ),
  
  nav_panel(
    "Gene Browser",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          checkboxGroupInput(
            "selected_datasets",
            "Choose datasets:",
            choices = unique(colData(se_combined)$experiment),
            selected = unique(colData(se_combined)$experiment)
          ),
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
          
          h3("Expression Values"),
          tableOutput("gene_expression_all")
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
          checkboxGroupInput(
            "heatmap_datasets",
            "Choose transcriptomic datasets:",
            choices = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Transcriptomics"
            ],
            selected = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Transcriptomics"
            ]
          ),
          selectizeInput(
            "heatmap_gene",
            "Choose a gene:",
            choices = rownames(rowData(se_combined)),
            selected = rownames(rowData(se_combined))[1]
          ),
          selectizeInput(
            "heatmap_conditions",
            "Search/select conditions:",
            choices = rownames(colData(se_combined))[
              colData(se_combined)$experiment %in%
                metadata(se_combined)$summary$experiment_name[
                  metadata(se_combined)$summary$omics_type == "Transcriptomics"
                ]
            ],
            selected = rownames(colData(se_combined))[
              colData(se_combined)$experiment %in%
                metadata(se_combined)$summary$experiment_name[
                  metadata(se_combined)$summary$omics_type == "Transcriptomics"
                ]
            ],
            multiple = TRUE
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
    "Scatterplot",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "scatter_gene",
            "Choose a gene:",
            choices = rownames(rowData(se_combined)),
            selected = rownames(rowData(se_combined))[1]
          ),
          
          h4("Transcriptomic Datasets (X‑axis)"),
          checkboxGroupInput(
            "scatter_transcript",
            "Select transcriptomic datasets:",
            choices = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Transcriptomics"
            ],
            selected = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Transcriptomics"
            ]
          ),
          
          h4("Proteomic Datasets (Y‑axis)"),
          checkboxGroupInput(
            "scatter_protein",
            "Select proteomic datasets:",
            choices = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Proteomics"
            ],
            selected = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Proteomics"
            ]
          )
        ),
        
        mainPanel(
          h3("Transcriptome vs Proteome Scatterplot"),
          plotlyOutput("scatter_plot", height = "600px")
        )
      )
    )
  ),
  
  nav_panel(
    "DEGs",
    fluidPage(
      br(),
      
      fluidRow(
        column(
          width = 4,
          card(
            card_header("Total DEGs"),
            card_body(
              h2(textOutput("total_degs"))
            )
          )
        ),
        column(
          width = 4,
          card(
            card_header("Upregulated"),
            card_body(
              h2(textOutput("up_degs"))
            )
          )
        ),
        column(
          width = 4,
          card(
            card_header("Downregulated"),
            card_body(
              h2(textOutput("down_degs"))
            )
          )
        )
      ),
      
      br(),
      
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "deg_filter",
            "Filter DEG direction:",
            choices = c("All", "Up", "Down")
          ),
          selectInput(
            "deg_experiment_filter",
            "Filter experiment:",
            choices = c("All", unique(metadata(se_combined)$deg_df$experiment))
          )
        ),
        mainPanel(
          h3("Differentially Expressed Genes"),
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
  
  output$gene_info <- renderTable({
    rowData(se_combined)[input$selected_gene, , drop = FALSE]
  })
  
  output$gene_expression_all <- renderTable({
    
    conds <- rownames(colData(se_combined))[
      colData(se_combined)$experiment %in% input$selected_datasets
    ]
    
    if (length(conds) == 0) {
      return(data.frame(Message = "No data available for selected datasets"))
    }
    
    log2fc_vals <- as.numeric(assay(se_combined)[input$selected_gene, conds])
    
    deg_df <- metadata(se_combined)$deg_df
    
    pvals <- vapply(conds, function(cn) {
      
      exp_name <- as.character(colData(se_combined)[cn, "experiment"])
      
      suffix <- sub("^log2_Fold_change_", "", cn)
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
      Condition = conds,
      Expression = log2fc_vals,
      `P-Value` = pvals_fmt,
      check.names = FALSE
    )
  })
  
  output$heatmap_plot <- renderPlotly({
    
    gene_index <- match(input$heatmap_gene, rownames(assay(se_combined)))
    if (is.na(gene_index)) return(NULL)
    
    start <- max(1, gene_index - floor(input$n_genes / 2))
    end <- min(nrow(assay(se_combined)), gene_index + floor(input$n_genes / 2))
    
    transcript_cols <- input$heatmap_conditions
    
    transcript_cols <- transcript_cols[
      colData(se_combined)[transcript_cols, "experiment"] %in% input$heatmap_datasets
    ]
    
    if (length(transcript_cols) == 0) {
      return(NULL)
    }
    
    mat <- assay(se_combined)[start:end, transcript_cols, drop = FALSE]
    
    df_long <- reshape2::melt(mat)
    colnames(df_long) <- c("Gene", "Condition", "Value")
    
    gene_info <- as.data.frame(rowData(se_combined))
    
    df_long$gene_name <- gene_info$gene_name[
      match(df_long$Gene, gene_info$gene_id)
    ]
    
    summary_df <- metadata(se_combined)$summary
    
    df_long$experiment <- as.character(
      colData(se_combined)[df_long$Condition, "experiment"]
    )
    
    df_long$publication <- summary_df$publication[
      match(df_long$experiment, summary_df$experiment_name)
    ]
    
    df_long$pmid <- summary_df$pmid[
      match(df_long$experiment, summary_df$experiment_name)
    ]
    
    df_long$pmcid <- summary_df$pmcid[
      match(df_long$experiment, summary_df$experiment_name)
    ]
    
    df_long$doi <- summary_df$doi[
      match(df_long$experiment, summary_df$experiment_name)
    ]
    
    deg_df <- metadata(se_combined)$deg_df
    
    df_long$p_value <- vapply(seq_len(nrow(df_long)), function(i) {
      gene <- df_long$Gene[i]
      cond <- df_long$Condition[i]
      
      exp_name <- as.character(colData(se_combined)[cond, "experiment"])
      
      suffix <- sub("^log2_Fold_change_", "", cond)
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
        selected_condition <- click$x
        
        gene_info <- as.data.frame(rowData(se_combined))
        deg_df <- metadata(se_combined)$deg_df
        
        exp_name <- as.character(colData(se_combined)[selected_condition, "experiment"])
        
        summary_df <- metadata(se_combined)$summary
        pub_row <- summary_df[summary_df$experiment_name == exp_name, ]
        
        suffix <- sub("^log2_Fold_change_", "", selected_condition)
        padj_col <- paste0("padj_", suffix)
        
        row_idx <- which(
          deg_df$gene_id == selected_gene &
            deg_df$experiment == exp_name
        )
        
        if (length(row_idx) == 1 && padj_col %in% colnames(deg_df)) {
          p_value <- as.numeric(deg_df[row_idx, padj_col])
          p_value_fmt <- formatC(p_value, format = "e", digits = 2)
        } else {
          p_value_fmt <- "N/A"
        }
        
        showModal(
          modalDialog(
            title = paste("Gene Details:", selected_gene),
            
            tags$strong("Gene Name:"),
            br(),
            gene_info$gene_name[
              match(selected_gene, gene_info$gene_id)
            ],
            
            br(), br(),
            
            tags$strong("Condition:"),
            br(),
            selected_condition,
            
            br(), br(),
            
            tags$strong("P-Value:"),
            br(),
            p_value_fmt,
            
            br(), br(),
            
            tags$strong("Publication:"),
            br(),
            pub_row$publication[1],
            
            br(), br(),
            
            if (!is.na(pub_row$pubmed_link[1])) {
              tags$div(
                tags$a(
                  href = pub_row$pubmed_link[1],
                  target = "_blank",
                  "Open PubMed Page"
                )
              )
            },
            
            if (!is.na(pub_row$pmc_link[1])) {
              tags$div(
                tags$a(
                  href = pub_row$pmc_link[1],
                  target = "_blank",
                  "Open PMC Page"
                )
              )
            },
            
            if (!is.na(pub_row$doi_link[1])) {
              tags$div(
                tags$a(
                  href = pub_row$doi_link[1],
                  target = "_blank",
                  "Open DOI Page"
                )
              )
            },
            
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      }
    }
  )
  
  output$scatter_plot <- renderPlotly({
    
    gene <- input$scatter_gene
    if (is.null(gene) || gene == "") return(NULL)
    
    # ---- Identify selected transcriptomic conditions ----
    transcript_cols <- rownames(colData(se_combined))[
      colData(se_combined)$experiment %in% input$scatter_transcript
    ]
    
    if (length(transcript_cols) == 0) return(NULL)
    
    # ---- Extract transcriptomic log2FC values ----
    x_vals <- as.numeric(assay(se_combined)[gene, transcript_cols])
    
    # ---- Extract proteomic absolute abundance (Schubert only) ----
    proteomic_cols <- rownames(colData(se_combined))[
      colData(se_combined)$experiment %in% input$scatter_protein
    ]
    
    if (length(proteomic_cols) == 0) return(NULL)
    
    # Schubert is absolute abundance → take the mean if multiple columns
    y_val <- mean(as.numeric(assay(se_combined)[gene, proteomic_cols]), na.rm = TRUE)
    
    # ---- Build dataframe: one point per transcriptomic condition ----
    df <- data.frame(
      transcript_condition = transcript_cols,
      transcript_exp = colData(se_combined)[transcript_cols, "experiment"],
      x = x_vals,
      y = y_val
    )
    
    df$hover <- paste0(
      "<b>Gene: ", gene, "</b>",
      "<br>Transcriptomic condition: ", df$transcript_condition,
      "<br>Experiment: ", df$transcript_exp,
      "<br><br>log2FC (Transcriptome): ", round(df$x, 3),
      "<br>Protein abundance (absolute): ", round(df$y, 3)
    )
    
    plot_ly(
      df,
      x = ~x,
      y = ~y,
      type = "scatter",
      mode = "markers",
      text = ~hover,
      hoverinfo = "text",
      marker = list(size = 14, color = "firebrick")
    ) %>%
      layout(
        xaxis = list(title = "Transcriptomic log2FC"),
        yaxis = list(title = "Proteomic abundance"),
        title = paste("Transcriptome vs Proteome for", gene)
      )
  })
  
  
  
  output$total_degs <- renderText({
    nrow(metadata(se_combined)$deg_df)
  })
  
  output$up_degs <- renderText({
    sum(metadata(se_combined)$deg_df$deg_direction == "Up", na.rm = TRUE)
  })
  
  output$down_degs <- renderText({
    sum(metadata(se_combined)$deg_df$deg_direction == "Down", na.rm = TRUE)
  })
  
  output$deg_table <- renderDT({
    df <- metadata(se_combined)$deg_df
    
    if (input$deg_filter != "All") {
      df <- df %>% filter(deg_direction == input$deg_filter)
    }
    
    if (input$deg_experiment_filter != "All") {
      df <- df %>% filter(experiment == input$deg_experiment_filter)
    }
    
    datatable(
      df,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE
      )
    )
  })
}

shinyApp(ui, server)