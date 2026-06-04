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
                               expression_scale = "Relative",
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
    expression_scale = expression_scale,
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
  expression_scale = "Relative",
  doi = "10.3389/fimmu.2022.909904",
  pmcid = "PMC9283954",
  publication = "Transcriptional profiling of Mycobacterium tuberculosis"
)

se_Shee <- load_experiment_se(
  file = "data/Transcriptome_Shee.xlsx",
  sheet = "Sheet1",
  experiment_name = "Shee Transcriptome",
  omics_type = "Transcriptomics",
  expression_scale = "Relative",
  pmid = "35975988",
  doi = "10.1128/AAC.00592-22",
  publication = "Moxifloxacin-Mediated Killing of Mycobacterium tuberculosis Involves Respiratory Downshift, Reductive Stress, and Accumulation of Reactive Oxygen Species"
)

se_Schubert <- load_experiment_se(
  file = "data/Proteome_Schubert_2015.xlsx",
  sheet = "Mtb absolute per condition",
  experiment_name = "Schubert Proteome",
  omics_type = "Proteomics",
  expression_scale = "Absolute",
  pmid = "26094805",
  doi = "10.1016/j.chom.2015.06.001",
  publication = "Schubert Proteome 2015"
)

se_Bei <- load_experiment_se(
  file = "data/Transcriptome_Bei_2024.xlsx",
  sheet = "log2FC of 5th and 95th",
  experiment_name = "Bei Transcriptome",
  omics_type = "Transcriptomics",
  expression_scale = "Absolute",
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
          accordion(
            id = "gb_sections",
            
            accordion_panel(
              title = "Gene Information",
              tableOutput("gene_info")
            ),
            
            accordion_panel(
              title = "Relative Expression Comparison",
              
              fluidRow(
                column(
                  width = 4,
                  selectizeInput(
                    "gb_condA",
                    "Condition A:",
                    choices = rownames(colData(se_combined))[
                      colData(se_combined)$experiment %in%
                        metadata(se_combined)$summary$experiment_name[
                          metadata(se_combined)$summary$omics_type == "Transcriptomics" &
                            metadata(se_combined)$summary$expression_scale == "Relative"
                        ]
                    ]
                  )
                ),
                column(
                  width = 4,
                  selectizeInput(
                    "gb_condB",
                    "Condition B:",
                    choices = rownames(colData(se_combined))[
                      colData(se_combined)$experiment %in%
                        metadata(se_combined)$summary$experiment_name[
                          metadata(se_combined)$summary$omics_type == "Transcriptomics" &
                            metadata(se_combined)$summary$expression_scale == "Relative"
                        ]
                    ]
                  )
                )
              ),
              
              tableOutput("gb_relcomp")
            ),
            
            accordion_panel(
              title = "Expression Values by Experiment",
              tableOutput("gene_expression_all")
            )
          )
          
          
        )
      )
    )
  ),
  
  nav_panel(
    "Relative Comparison",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "relcomp_genes",
            "Select one or more genes:",
            choices = rownames(rowData(se_combined)),
            selected = rownames(rowData(se_combined))[1],
            multiple = TRUE
          ),
          
          h4("Condition A"),
          selectizeInput(
            "relcomp_condA",
            "Choose Condition A:",
            choices = rownames(colData(se_combined))[
              colData(se_combined)$experiment %in%
                metadata(se_combined)$summary$experiment_name[
                  metadata(se_combined)$summary$omics_type == "Transcriptomics" &
                    metadata(se_combined)$summary$expression_scale == "Relative"
                ]
            ]
          ),
          
          h4("Condition B"),
          selectizeInput(
            "relcomp_condB",
            "Choose Condition B:",
            choices = rownames(colData(se_combined))[
              colData(se_combined)$experiment %in%
                metadata(se_combined)$summary$experiment_name[
                  metadata(se_combined)$summary$omics_type == "Transcriptomics" &
                    metadata(se_combined)$summary$expression_scale == "Relative"
                ]
            ]
          )
        ),
        
        mainPanel(
          h3("Relative Expression Comparison"),
          tableOutput("relcomp_table")
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
              metadata(se_combined)$summary$omics_type == "Transcriptomics" &
                metadata(se_combined)$summary$expression_scale == "Relative"
            ],
            selected = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Transcriptomics" &
                metadata(se_combined)$summary$expression_scale == "Relative"
            ]
          ),
          
          uiOutput("heatmap_condition_selector"),
          
          # --- Gene selection mode toggle ---
          radioButtons(
            "heatmap_gene_mode",
            "Gene selection method:",
            choices = c("Individual genes" = "single",
                        "Nearby gene window" = "window"),
            selected = "single"
          ),
          
          # --- Dynamic UI for gene selection ---
          uiOutput("heatmap_gene_selector")
          
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
            "scatter_genes",
            "Choose one or more genes:",
            choices = rownames(rowData(se_combined)),
            selected = rownames(rowData(se_combined))[1],
            multiple = TRUE
          ),
          
          h4("Transcriptomic Datasets (X‑axis)"),
          checkboxGroupInput(
            "scatter_transcript",
            "Select transcriptomic datasets:",
            choices = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Transcriptomics" &
                metadata(se_combined)$summary$expression_scale == "Absolute"
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
              metadata(se_combined)$summary$omics_type == "Proteomics" &
                metadata(se_combined)$summary$expression_scale == "Absolute"
            ],
            selected = metadata(se_combined)$summary$experiment_name[
              metadata(se_combined)$summary$omics_type == "Proteomics"
            ]
          ),
          selectInput(
            "scatter_scale_x",
            "X‑axis scaling:",
            choices = c("Linear" = "linear", "Log2" = "log2"),
            selected = "linear"
          ),
          
          selectInput(
            "scatter_scale_y",
            "Y‑axis scaling:",
            choices = c("Linear" = "linear", "Log2" = "log2"),
            selected = "linear"
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
          numericInput(
            "deg_magnitude",
            "Minimum |log2FC| magnitude:",
            value = 0,
            min = 0,
            max = 20,
            step = 0.1
          ),
          
          selectizeInput(
            "deg_experiments",
            "Select datasets:",
            choices = unique(metadata(se_combined)$deg_df$experiment),
            multiple = TRUE,
            selected = unique(metadata(se_combined)$deg_df$experiment)
          ),
          
          uiOutput("deg_condition_selector")
          
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
  
  output$gb_relcomp <- renderTable({
    
    gene <- input$selected_gene
    condA <- input$gb_condA
    condB <- input$gb_condB
    
    if (is.null(gene) || is.null(condA) || is.null(condB)) {
      return(data.frame(Message = "Select two conditions to compare."))
    }
    
    # Extract log2FC values
    valA <- as.numeric(assay(se_combined)[gene, condA])
    valB <- as.numeric(assay(se_combined)[gene, condB])
    
    # Compute relative log2FC
    rel <- valA - valB
    
    data.frame(
      Metric = c("Condition A log2FC", "Condition B log2FC", "Relative log2FC (A - B)"),
      Value = c(valA, valB, rel),
      check.names = FALSE
    )
  })
  
  
  output$gene_expression_all <- renderUI({
    
    gene <- input$selected_gene
    
    # Identify selected conditions
    conds <- rownames(colData(se_combined))[
      colData(se_combined)$experiment %in% input$selected_datasets
    ]
    
    if (length(conds) == 0) {
      return(tags$div("No data available for selected datasets"))
    }
    
    # Extract values
    log2fc_vals <- as.numeric(assay(se_combined)[gene, conds])
    
    # Metadata
    deg_df <- metadata(se_combined)$deg_df
    summary_df <- metadata(se_combined)$summary
    
    # Compute p-values
    pvals <- vapply(conds, function(cn) {
      exp_name <- as.character(colData(se_combined)[cn, "experiment"])
      suffix <- sub("^log2_Fold_change_", "", cn)
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
    
    pvals_fmt <- ifelse(
      is.na(pvals),
      NA,
      formatC(pvals, format = "e", digits = 2)
    )
    
    # Build master df
    df <- data.frame(
      Condition = conds,
      Expression = log2fc_vals,
      `P-Value` = pvals_fmt,
      Experiment = colData(se_combined)[conds, "experiment"],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    # Split by experiment
    df_split <- split(df, df$Experiment)
    
    # Build UI output
    ui_list <- lapply(names(df_split), function(exp_name) {
      
      exp_row <- summary_df[summary_df$experiment_name == exp_name, ]
      
      subtitle <- paste0(
        exp_name, " – ",
        exp_row$expression_scale, " ",
        exp_row$omics_type
      )
      
      tagList(
        tags$h4(subtitle),
        tableOutput(paste0("tbl_", exp_name)),
        tags$br()
      )
    })
    
    # Render each table separately
    for (exp_name in names(df_split)) {
      local({
        nm <- exp_name
        df_exp <- df_split[[nm]]
        
        # Identify metadata columns (anything not Condition/Expression/Experiment)
        meta_cols <- setdiff(colnames(df_exp), c("Condition", "Expression", "Experiment"))
        
        # Keep only metadata columns that have at least one non-NA value
        keep_meta <- meta_cols[sapply(meta_cols, function(col) any(!is.na(df_exp[[col]])))]
        
        # Build final column order
        final_cols <- c("Condition", "Expression", keep_meta)
        
        output[[paste0("tbl_", nm)]] <- renderTable({
          df_exp[, final_cols, drop = FALSE]
        })
        
      })
    }
    
    do.call(tagList, ui_list)
  })
  
  output$relcomp_table <- renderTable({
    
    genes <- input$relcomp_genes
    condA <- input$relcomp_condA
    condB <- input$relcomp_condB
    
    if (is.null(genes) || is.null(condA) || is.null(condB)) {
      return(data.frame(Message = "Please select genes and both conditions."))
    }
    
    # Extract log2FC values
    valsA <- as.numeric(assay(se_combined)[genes, condA])
    valsB <- as.numeric(assay(se_combined)[genes, condB])
    
    # Compute relative log2FC
    rel_vals <- valsA - valsB
    
    # Build output table
    data.frame(
      Gene = genes,
      Condition_A = valsA,
      Condition_B = valsB,
      Relative_log2FC = rel_vals,
      check.names = FALSE
    )
  })
  
  
  output$heatmap_plot <- renderPlotly({
    
    # --- MULTI‑GENE SUPPORT ---
    genes <- input$heatmap_gene
    if (is.null(genes) || length(genes) == 0) return(NULL)
    
    # Convert gene IDs to row indices
    gene_index <- match(genes, rownames(assay(se_combined)))
    gene_index <- gene_index[!is.na(gene_index)]
    if (length(gene_index) == 0) return(NULL)
    
    # --- RANGE SELECTION ---
    if (input$heatmap_gene_mode == "single") {
      
      # Show exactly the selected genes
      start <- min(gene_index)
      end   <- max(gene_index)
      
    } else {
      
      # Window mode (center gene + window size)
      nwin <- input$n_genes
      if (is.null(nwin) || is.na(nwin) || nwin < 1) return(NULL)
      
      center <- gene_index[1]   # first selected gene is the center
      half <- floor(nwin / 2)
      
      start <- max(1, center - half)
      end   <- min(nrow(assay(se_combined)), center + half)
    }
    
    # --- CONDITION FILTERING ---
    transcript_cols <- input$heatmap_conditions
    transcript_cols <- transcript_cols[
      colData(se_combined)[transcript_cols, "experiment"] %in% input$heatmap_datasets
    ]
    if (length(transcript_cols) == 0) return(NULL)
    
    # --- EXTRACT MATRIX ---
    mat <- assay(se_combined)[start:end, transcript_cols, drop = FALSE]
    
    # --- LONG FORMAT ---
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
    
    # --- P‑VALUES ---
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
    
    # --- HOVER TEXT ---
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
    
    # --- PLOT ---
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
  
  
  output$heatmap_condition_selector <- renderUI({
    
    exps <- input$heatmap_datasets
    
    if (is.null(exps) || length(exps) == 0) {
      return(tags$div("Select one or more datasets to choose conditions."))
    }
    
    # Find all conditions belonging to selected experiments
    conds <- rownames(colData(se_combined))[
      colData(se_combined)$experiment %in% exps
    ]
    
    checkboxGroupInput(
      "heatmap_conditions",
      "Select conditions:",
      choices = conds,
      selected = conds
    )
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
  
  output$heatmap_gene_selector <- renderUI({
    
    if (input$heatmap_gene_mode == "single") {
      
      selectizeInput(
        "heatmap_gene",
        "Choose one or more genes:",
        choices = rownames(rowData(se_combined)),
        selected = rownames(rowData(se_combined))[1],
        multiple = TRUE,
        options = list(
          placeholder = "Type gene names…",
          maxOptions = 5000
        )
      )
      
      
      
    } else {
      
      # Window mode
      tagList(
        selectizeInput(
          "heatmap_gene",
          "Choose a gene:",
          choices = rownames(rowData(se_combined)),
          selected = rownames(rowData(se_combined))[1],
          options = list(
            placeholder = "Type a gene name…",
            maxOptions = 5000
          )
        )
        ,
        numericInput(
          "n_genes",
          "Nearby genes:",
          value = 20,
          min = 5,
          max = 100
        )
      )
    }
  })
  
  
  output$scatter_plot <- renderPlotly({
    
    selected_genes <- input$scatter_genes
    if (is.null(selected_genes) || length(selected_genes) == 0) return(NULL)
    
    # ---- Identify Bei (absolute transcriptome) ----
    bei_cols <- rownames(colData(se_combined))[
      colData(se_combined)$experiment == "Bei Transcriptome"
    ]
    
    # ---- Identify Schubert (absolute proteome) ----
    schubert_cols <- rownames(colData(se_combined))[
      colData(se_combined)$experiment == "Schubert Proteome"
    ]
    
    if (length(bei_cols) == 0 || length(schubert_cols) == 0) return(NULL)
    
    # ---- Extract absolute transcript abundance ----
    transcript_vals <- as.numeric(assay(se_combined)[, bei_cols])
    transcript_vals <- 2 ^ transcript_vals   # convert log2 → absolute
    
    # ---- Extract absolute protein abundance ----
    protein_vals <- as.numeric(assay(se_combined)[, schubert_cols])
    
    # ---- Build dataframe: one point per gene ----
    df <- data.frame(
      gene_id = rownames(se_combined),
      transcript = transcript_vals,
      protein = protein_vals,
      stringsAsFactors = FALSE
    )
    
    # ---- Independent axis scaling ----
    
    # X‑axis
    if (input$scatter_scale_x == "log2") {
      df <- df[df$transcript > 0, ]   # remove invalid values
      df$transcript_scaled <- log2(df$transcript)
    } else {
      df$transcript_scaled <- df$transcript
    }
    
    # Y‑axis
    if (input$scatter_scale_y == "log2") {
      df <- df[df$protein > 0, ]      # remove invalid values
      df$protein_scaled <- log2(df$protein)
    } else {
      df$protein_scaled <- df$protein
    }
    
    
    
    # ---- Mark selected genes ----
    df$selected <- df$gene_id %in% selected_genes
    
    # ---- Hover text ----
    df$hover <- paste0(
      "<b>", df$gene_id, "</b>",
      "<br>Transcript (absolute): ", round(df$transcript, 2),
      "<br>Protein (absolute): ", round(df$protein, 2)
    )
    
    # ---- Plot ----
    plot_ly() %>%
      
      # Background points (all genes)
      add_markers(
        data = df[!df$selected, ],
        x = ~transcript,
        y = ~protein,
        marker = list(size = 6, color = "lightgrey"),
        hoverinfo = "text",
        text = ~hover,
        name = "All genes"
      ) %>%
      
      # Highlighted points (selected genes)
      add_markers(
        data = df[df$selected, ],
        x = ~transcript_scaled,
        y = ~protein_scaled,
        marker = list(size = 12, color = "firebrick"),
        hoverinfo = "text",
        text = ~hover,
        name = "Selected genes"
      ) %>%
      
      # Add text labels for selected genes
      add_text(
        data = df[df$selected, ],
        x = ~transcript_scaled,
        y = ~protein_scaled,
        text = ~gene_id,
        textposition = "top center",
        textfont = list(size = 14, color = "black"),
        showlegend = FALSE,
        hoverinfo = "none"
      ) %>%
      
      layout(
        xaxis = list(title = "Absolute Transcript Abundance (Bei)"),
        yaxis = list(title = "Absolute Protein Abundance (Schubert)"),
        title = "Absolute Transcriptome vs Proteome"
      )
  })
  
  output$deg_condition_selector <- renderUI({
    
    exps <- input$deg_experiments
    
    if (is.null(exps) || length(exps) == 0) {
      return(tags$div("Select one or more datasets to choose conditions."))
    }
    
    # Find all conditions belonging to selected experiments
    conds <- rownames(colData(se_combined))[
      colData(se_combined)$experiment %in% exps
    ]
    
    checkboxGroupInput(
      "deg_conditions",
      "Select conditions:",
      choices = conds,
      selected = conds
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
    
    deg_df <- metadata(se_combined)$deg_df
    
    # ---- Filter rows by selected datasets ----
    if (!is.null(input$deg_experiments) && length(input$deg_experiments) > 0) {
      deg_df <- deg_df[deg_df$experiment %in% input$deg_experiments, ]
    }
    
    # ---- Determine selected conditions ----
    conds <- input$deg_conditions
    
    # Build fc column names
    fc_cols <- paste0("fc_", conds)
    fc_cols <- fc_cols[fc_cols %in% colnames(deg_df)]
    
    # ---- Direction filtering across selected conditions ----
    if (input$deg_filter != "All") {
      
      if (input$deg_filter == "Up") {
        keep_dir <- apply(deg_df[, fc_cols, drop = FALSE], 1, function(x) any(x > 0, na.rm = TRUE))
      } else if (input$deg_filter == "Down") {
        keep_dir <- apply(deg_df[, fc_cols, drop = FALSE], 1, function(x) any(x < 0, na.rm = TRUE))
      } else {
        keep_dir <- apply(deg_df[, fc_cols, drop = FALSE], 1, function(x) any(x == 0, na.rm = TRUE))
      }
      
      deg_df <- deg_df[keep_dir, ]
    }
    
    # ---- Magnitude filtering across selected conditions ----
    mag <- input$deg_magnitude
    
    if (!is.null(mag) && mag > 0) {
      keep_mag <- apply(deg_df[, fc_cols, drop = FALSE], 1, function(x) any(abs(x) >= mag, na.rm = TRUE))
      deg_df <- deg_df[keep_mag, ]
    }
    
    # ---- Determine which columns to keep ----
    base_cols <- c("gene_id", "gene_name", "deg_direction")
    
    # padj columns
    padj_cols <- paste0("padj_", sub("^log2_Fold_change_", "", conds))
    padj_cols <- padj_cols[padj_cols %in% colnames(deg_df)]
    
    # sd columns
    sd_cols <- paste0("sd_", sub("^log2_Fold_change_", "", conds))
    sd_cols <- sd_cols[sd_cols %in% colnames(deg_df)]
    
    keep_cols <- c(base_cols, fc_cols, padj_cols, sd_cols)
    keep_cols <- keep_cols[keep_cols %in% colnames(deg_df)]
    
    df_out <- deg_df[, keep_cols, drop = FALSE]
    
    datatable(df_out, rownames = FALSE)
  })
  
  
}

shinyApp(ui, server)