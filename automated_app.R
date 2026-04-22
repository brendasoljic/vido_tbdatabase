
load_experiment <- function(file, sheet, 
                            experiment_name,
                            species = "M. tuberculosis H37Rv",
                            omics_type = "Transcriptomics") {
  
  library(readxl)
  library(dplyr)
  
  # ---- Load data ----
  df <- read_excel(file, sheet = sheet)
  
  # Helper: safely choose a column name or fallback
  safe_pick <- function(x, fallback = NA_character_) {
    if (length(x) == 0 || is.na(x)) fallback else x
  }
  
  # ---- Detect gene columns ----
  gene_id_col <- names(df)[grepl("Gene$|Rv|Primary Target", names(df), ignore.case = TRUE)][1]
  gene_name_col <- names(df)[grepl("Gene Name|Common Name", names(df), ignore.case = TRUE)][1]
  gene_symbol_col <- names(df)[grepl("Gene Symbol", names(df), ignore.case = TRUE)][1]
  
  # Apply safe fallback logic
  gene_id_col     <- safe_pick(gene_id_col,     fallback = names(df)[1])   # fallback to first column
  gene_name_col   <- safe_pick(gene_name_col,   fallback = gene_id_col)    # fallback to gene_id
  gene_symbol_col <- safe_pick(gene_symbol_col, fallback = gene_id_col)    # fallback to gene_id
  
  # ---- Detect assay columns (log2 fold change) ----
  assay_cols <- names(df)[grepl("log2", names(df), ignore.case = TRUE)]
  
  # If none found, create an empty vector (Shiny app will handle gracefully)
  if (length(assay_cols) == 0) {
    assay_cols <- character(0)
  }
  
  # ---- Detect padj columns ----
  padj_cols <- names(df)[grepl("padj|FDR|adj", names(df), ignore.case = TRUE)]
  
  if (length(padj_cols) == 0) {
    padj_cols <- character(0)
  }
  
  
  # ---- Build genes_df ----
  genes_df <- df %>%
    transmute(
      gene_id = .data[[gene_id_col]],
      orf = .data[[gene_id_col]],
      gene_name = .data[[gene_name_col]],
      gene_symbol = .data[[gene_symbol_col]]
    )
  
  # ---- Build expression matrix ----
  expr_mat <- df %>%
    select(all_of(assay_cols)) %>%
    as.matrix()
  
  rownames(expr_mat) <- genes_df$gene_id
  expr_mat <- apply(expr_mat, 2, as.numeric)
  expr_mat <- as.matrix(expr_mat)
  
  # ---- Build DEG table ----
  deg_df <- df %>%
    transmute(
      gene_id = .data[[gene_id_col]],
      gene_name = .data[[gene_name_col]],
      gene_symbol = .data[[gene_symbol_col]],
      across(all_of(assay_cols), .names = "fc_{col}"),
      across(all_of(padj_cols), .names = "padj_{col}"),
      deg_direction = case_when(
        !!sym(assay_cols[1]) > 0 ~ "Up",
        !!sym(assay_cols[1]) < 0 ~ "Down",
        TRUE ~ "No change"
      )
    )
  
  # ---- Build experiment summary ----
  experiment_summary <- data.frame(
    experiment_name = experiment_name,
    omics_type = omics_type,
    species = species,
    n_genes = nrow(genes_df),
    n_deg = nrow(deg_df),
    n_samples = ncol(expr_mat),
    stringsAsFactors = FALSE
  )
  
  # ---- Return everything ----
  list(
    genes_df = genes_df,
    expr_mat = expr_mat,
    deg_df = deg_df,
    experiment_summary = experiment_summary
  )
}

Shee <- load_experiment(file = "C:/Users/mayac/OneDrive/Documents/VIDO/Transcriptome_Shee.xlsx", sheet = "Sheet1", experiment_name = "Shee")

print(Shee)

