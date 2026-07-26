# =============================================================================
# 04_experiment_inputs.R - Enter per-taxon cell size
# =============================================================================
# A small Shiny app for the per-strain cell-size input that anchors the carbon
# per cell (and therefore absolute growth):
#
#   Cell size:
#     Enter each taxon's rod width + length (um) and a carbon density; writes
#         data/taxon_cell_sizes.csv   (per-taxon carbon per cell)
#     used by 05_oxygen_fits.R (via config's cell_carbon_of()).
#
# NOTE: This app does NOT touch data/Ninoc.csv. The inoculation table is
#       maintained separately and must not be overwritten here.
#
# RUN:  RStudio -> open this file -> "Run App"   |   Terminal -> Rscript 04_experiment_inputs.R
# =============================================================================

# ---- Locate the project -----------------------------------------------------
.script_dir <- local({
  d <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path))
      dirname(rstudioapi::getActiveDocumentContext()$path) else NA_character_
  }, error = function(e) NA_character_)
  if (length(d) == 0 || is.na(d) || !nzchar(d))
    d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (length(d) == 0 || is.na(d) || !nzchar(d)) {
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(a)) d <- dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE))
  }
  if (length(d) == 0 || is.na(d) || !nzchar(d)) d <- getwd()
  normalizePath(d, mustWork = FALSE)
})
base_dir   <- dirname(.script_dir)
data_dir   <- file.path(base_dir, "data")
tables_dir <- file.path(base_dir, "results", "tables")
for (d in c(data_dir, tables_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
cell_csv  <- file.path(data_dir, "taxon_cell_sizes.csv")

# ---- Packages ---------------------------------------------------------------
need <- c("shiny", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) stop("Install first: install.packages(c(",
                           paste(sprintf('"%s"', miss), collapse = ", "), "))")
library(shiny)

rod_volume <- function(width, length) {
  rad <- width / 2; cyl <- max(length - width, 0)
  pi * rad^2 * cyl + (4/3) * pi * rad^3
}
read_ok <- function(p) if (file.exists(p)) tryCatch(readr::read_csv(p, show_col_types = FALSE), error = function(e) NULL) else NULL

# ---- Which taxa? ------------------------------------------------------------
taxa <- NULL
long_csv <- file.path(tables_dir, "Oxygen_All_Long.csv")
lg <- read_ok(long_csv)
if (!is.null(lg) && "Taxon" %in% names(lg)) taxa <- sort(unique(as.character(lg$Taxon)))
if (is.null(taxa)) {
  cc <- read_ok(file.path(data_dir, "Cell_Counts.csv"))
  if (!is.null(cc) && "Taxon" %in% names(cc)) taxa <- sort(unique(as.character(cc$Taxon)))
}
if (is.null(taxa)) taxa <- c("Bacillus", "Yersinia")

# ---- Prefill helpers --------------------------------------------------------
DENSITY_DEFAULT <- 100; WIDTH_DEFAULT <- 0.65; LENGTH_DEFAULT <- 2.25
cell_prev <- read_ok(cell_csv)
cell_prev_of <- function(tx, col, default) {
  if (!is.null(cell_prev) && all(c("Taxon", col) %in% names(cell_prev))) {
    v <- suppressWarnings(as.numeric(cell_prev[[col]][as.character(cell_prev$Taxon) == tx]))
    v <- v[is.finite(v)]; if (length(v)) return(v[1])
  }
  default
}

# ---- UI ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Experiment inputs - cell size, per taxon"),
  sidebarLayout(
    sidebarPanel(width = 5,
      helpText("Rod width + length (um) per taxon; density converts volume to carbon per cell."),
      numericInput("density", "Carbon density (fg C / um^3):", value = DENSITY_DEFAULT, min = 1, step = 10),
      tags$hr(), strong("Cell dimensions (um):"), uiOutput("dim_inputs"), tags$hr(),
      actionButton("save_cell", "Save taxon_cell_sizes.csv", class = "btn-primary"),
      tags$br(), tags$br(), verbatimTextOutput("cell_status")),
    mainPanel(width = 7, h4("Preview"), tableOutput("cell_preview"))
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {
  output$dim_inputs <- renderUI({
    lapply(taxa, function(tx) fluidRow(column(4, tags$b(tx)),
      column(4, numericInput(paste0("w_", make.names(tx)), "width",
                             value = cell_prev_of(tx, "cell_width_um", WIDTH_DEFAULT), min = 0.01, step = 0.05)),
      column(4, numericInput(paste0("l_", make.names(tx)), "length",
                             value = cell_prev_of(tx, "cell_length_um", LENGTH_DEFAULT), min = 0.01, step = 0.05))))
  })
  cell_table <- reactive({
    dens <- if (is.null(input$density) || !is.finite(input$density)) DENSITY_DEFAULT else input$density
    rows <- lapply(taxa, function(tx) {
      w <- input[[paste0("w_", make.names(tx))]]; l <- input[[paste0("l_", make.names(tx))]]
      if (is.null(w) || is.null(l) || !is.finite(w) || !is.finite(l) || w <= 0 || l <= 0) return(NULL)
      v <- rod_volume(w, l)
      data.frame(Taxon = tx, cell_width_um = w, cell_length_um = l, cell_volume_um3 = round(v, 4),
                 carbon_density_fg_per_um3 = dens, cell_carbon_fg = round(v * dens, 4), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  })
  output$cell_preview <- renderTable({ cell_table() }, digits = 3)
  observeEvent(input$save_cell, {
    df <- cell_table()
    if (is.null(df) || nrow(df) == 0) { output$cell_status <- renderText("Nothing to save."); return() }
    ok <- tryCatch({ readr::write_csv(df, cell_csv); TRUE },
                   error = function(e) { output$cell_status <- renderText(paste("FAILED:", conditionMessage(e))); FALSE })
    if (isTRUE(ok)) output$cell_status <- renderText(paste0("Saved ", nrow(df), " taxa to:\n", cell_csv,
                                                           "\nRe-run 05_oxygen_fits.R to apply."))
  })
}

if (interactive()) {
  shinyApp(ui, server)
} else {
  message("Launching experiment-inputs app at http://127.0.0.1:7799 ...")
  runApp(shinyApp(ui, server), host = "127.0.0.1", port = 7799, launch.browser = TRUE)
}
