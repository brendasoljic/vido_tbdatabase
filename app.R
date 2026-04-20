library(shiny)
library(DT)
library(dplyr)
library(pheatmap)
library(readxl)
library(bslib)

# =========================
# 1. LOAD THE REAL DATA
# =========================

# Your Excel file should be here:
# vido_tbdatabase/data/Transcriptome_Moxifloxacin.xlsx

file_path <- "data/Transcriptome_Moxifloxacin.xlsx"

# Skip the first 2 rows because row 1 is the description, row 2 is blank
raw_expr <- read_excel(file_path, sheet = "Sheet1", skip = 2)
raw_deg  <- read_excel(file_path, sheet = "Sheet2", skip = 2)

# Clean and standardize column names for the main expression sheet
genes_df <- raw_expr %>%
  transmute(
    gene_id = `Primary Target`,
    orf = ORF,
    gene_name = `Common Name of Primary Target`,
    gene_symbol = `Gene Symbol`
  )

# Build the expression matrix
# These are the 3 moxifloxacin treatment columns
expr_mat <- raw_expr %>%
  select(moxi2x_GSM1829746, moxi4x_GSM1829747, moxi8x_GSM1829748) %>%
  as.matrix()

# Use gene IDs as row names
rownames(expr_mat) <- genes_df$gene_id

# Force numeric just in case Excel imports something oddly
expr_mat <- apply(expr_mat, 2, as.numeric)
expr_mat <- as.matrix(expr_mat)
rownames(expr_mat) <- genes_df$gene_id
colnames(expr_mat) <- c("Moxi 2x", "Moxi 4x", "Moxi 8x")

# Clean DEG table
deg_df <- raw_deg %>%
  transmute(
    gene_id = RvId,
    gene_name = `Common Name of Primary Target`,
    gene_symbol = `Gene Symbol`,
    moxi2x = moxi2x_GSM1829746,
    moxi4x = moxi4x_GSM1829747,
    moxi8x = moxi8x_GSM1829748,
    deg_direction = DEG
  )

# Simple experiment summary table
experiment_summary <- data.frame(
  experiment_name = "Moxifloxacin transcriptome",
  omics_type = "Transcriptomics",
  species = "M. tuberculosis H37Rv",
  treatment = "Moxifloxacin",
  n_genes = nrow(genes_df),
  n_deg = nrow(deg_df),
  n_samples = ncol(expr_mat),
  stringsAsFactors = FALSE
)

# Summary counts for the dashboard
total_experiments <- 1
transcriptomics_experiments <- 1
proteomics_experiments <- 0

# =========================
# 2. USER INTERFACE (UI)
# =========================

ui <- page_navbar(
  title = "MtB Experiment Dashboard",
  theme = bs_theme(bootswatch = "flatly"),
  
  # ---- Overview tab ----
  nav_panel(
    "Overview",
    fluidPage(
      br(),
      fluidRow(
        column(
          width = 4,
          card(
            card_header("Total Experiments"),
            card_body(h2(total_experiments))
          )
        ),
        column(
          width = 4,
          card(
            card_header("Transcriptomic Experiments"),
            card_body(h2(transcriptomics_experiments))
          )
        ),
        column(
          width = 4,
          card(
            card_header("Proteomic Experiments"),
            card_body(h2(proteomics_experiments))
          )
        )
      ),
      br(),
      h3("Loaded Experiment"),
      DTOutput("experiment_table")
    )
  ),
  
  # ---- Gene browser tab ----
  nav_panel(
    "Gene Browser",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "selected_gene",
            "Search gene by Rv ID:",
            choices = genes_df$gene_id,
            selected = genes_df$gene_id[1],
            options = list(placeholder = "Type a gene ID...")
          ),
          checkboxInput("show_symbol_only", "Show genes with gene symbols only", value = FALSE)
        ),
        mainPanel(
          h3("Gene Table"),
          DTOutput("gene_table"),
          br(),
          h3("Selected Gene Information"),
          tableOutput("gene_info"),
          br(),
          h3("Selected Gene Expression Values"),
          tableOutput("gene_expression_values")
        )
      )
    )
  ),
  
  # ---- Heatmap tab ----
  nav_panel(
    "Heatmap",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectizeInput(
            "heatmap_gene",
            "Choose a gene:",
            choices = genes_df$gene_id,
            selected = genes_df$gene_id[1]
          ),
          selectInput(
            "heatmap_mode",
            "Heatmap mode:",
            choices = c("Selected gene only", "Selected gene + nearby genes", "Top variable genes"),
            selected = "Selected gene + nearby genes"
          ),
          numericInput(
            "n_genes",
            "Number of genes to show:",
            value = 20,
            min = 2,
            max = 100
          )
        ),
        mainPanel(
          h3("Expression Heatmap"),
          plotOutput("heatmap_plot", height = "600px")
        )
      )
    )
  ),
  
  # ---- DEG tab ----
  nav_panel(
    "DEGs",
    fluidPage(
      br(),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "deg_filter",
            "Filter DEG direction:",
            choices = c("All", "Up", "Down"),
            selected = "All"
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

# =========================
# 3. SERVER LOGIC
# =========================

server <- function(input, output, session) {
  
  # ---- Overview table ----
  output$experiment_table <- renderDT({
    datatable(
      experiment_summary,
      rownames = FALSE,
      options = list(dom = "t", scrollX = TRUE)
    )
  })
  
  # ---- Filtered gene table ----
  filtered_genes <- reactive({
    df <- genes_df
    
    if (isTRUE(input$show_symbol_only)) {
      df <- df %>% filter(!is.na(gene_symbol), gene_symbol != "")
    }
    
    df
  })
  
  output$gene_table <- renderDT({
    datatable(
      filtered_genes(),
      selection = "single",
      filter = "top",
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  # Update gene selector when user clicks a row in the gene table
  observeEvent(input$gene_table_rows_selected, {
    row_index <- input$gene_table_rows_selected
    current_table <- filtered_genes()
    
    if (length(row_index) == 1) {
      selected_id <- current_table$gene_id[row_index]
      
      updateSelectizeInput(
        session,
        "selected_gene",
        selected = selected_id
      )
      
      updateSelectizeInput(
        session,
        "heatmap_gene",
        selected = selected_id
      )
    }
  })
  
  # ---- Selected gene info ----
  selected_gene_info <- reactive({
    req(input$selected_gene)
    genes_df %>% filter(gene_id == input$selected_gene)
  })
  
  output$gene_info <- renderTable({
    selected_gene_info()
  })
  
  # ---- Selected gene expression values ----
  output$gene_expression_values <- renderTable({
    req(input$selected_gene)
    
    data.frame(
      condition = colnames(expr_mat),
      log2_fold_change = as.numeric(expr_mat[input$selected_gene, ])
    )
  })
  
  # ---- Heatmap data ----
  heatmap_data <- reactive({
    req(input$heatmap_gene)
    
    if (input$heatmap_mode == "Selected gene only") {
      mat <- expr_mat[input$heatmap_gene, , drop = FALSE]
      
    } else if (input$heatmap_mode == "Selected gene + nearby genes") {
      gene_index <- match(input$heatmap_gene, rownames(expr_mat))
      
      start <- max(1, gene_index - floor(input$n_genes / 2))
      end <- min(nrow(expr_mat), gene_index + floor(input$n_genes / 2))
      
      mat <- expr_mat[start:end, , drop = FALSE]
      
    } else {
      vars <- apply(expr_mat, 1, var, na.rm = TRUE)
      top_ids <- names(sort(vars, decreasing = TRUE))[1:min(input$n_genes, length(vars))]
      mat <- expr_mat[top_ids, , drop = FALSE]
    }
    
    mat
  })
  
  output$heatmap_plot <- renderPlot({
    mat <- heatmap_data()
    
    # pheatmap does not like a true 1-row heatmap very much,
    # so duplicate the row very slightly if only one gene is shown
    if (nrow(mat) == 1) {
      mat <- rbind(mat, mat + 1e-6)
      rownames(mat)[2] <- paste0(rownames(mat)[1], "_copy")
    }
    
    pheatmap(
      mat,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = paste("Heatmap:", input$heatmap_gene),
      fontsize_row = 8,
      fontsize_col = 12
    )
  })
  
  # ---- DEG filter ----
  filtered_deg <- reactive({
    df <- deg_df
    
    if (input$deg_filter != "All") {
      df <- df %>% filter(deg_direction == input$deg_filter)
    }
    
    df
  })
  
  output$deg_table <- renderDT({
    datatable(
      filtered_deg(),
      filter = "top",
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
}

# =========================
# 4. RUN THE APP
# =========================
shinyApp(ui = ui, server = server)