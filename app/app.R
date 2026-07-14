library(shiny)
library(DT)

project_dir <- normalizePath(
  if (file.exists(file.path("scripts", "analysis_functions.R"))) "." else "..",
  mustWork = TRUE
)
source(file.path(project_dir, "scripts", "analysis_functions.R"))

gene_choices <- available_genes("total", file.path(project_dir, analysis_data_dir))
dataset_choices <- names(analysis_datasets(file.path(project_dir, analysis_data_dir)))
default_gene <- if ("CD82" %in% gene_choices) "CD82" else gene_choices[[1]]
default_target <- if ("MYC" %in% gene_choices) "MYC" else gene_choices[[1]]

named_labels <- function(x) {
  stats::setNames(names(x), vapply(x, function(item) item[["label"]] %||% item[["title"]] %||% "", character(1)))
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background: #f7f8fa; color: #1f2933; }
      .container-fluid { max-width: 1480px; }
      .app-shell { display: grid; grid-template-columns: 280px minmax(0, 1fr); gap: 18px; padding: 18px 8px 24px; }
      .sidebar-panel, .main-panel { background: #ffffff; border: 1px solid #d9dee7; border-radius: 6px; }
      .sidebar-panel { padding: 14px; align-self: start; position: sticky; top: 12px; }
      .main-panel { padding: 12px 14px 16px; min-width: 0; }
      h1 { font-size: 20px; margin: 0 0 14px; font-weight: 650; }
      h2 { font-size: 15px; margin: 12px 0 10px; font-weight: 650; }
      .form-group { margin-bottom: 12px; }
      .control-label { font-size: 12px; font-weight: 650; color: #344054; margin-bottom: 4px; }
      .selectize-input, .form-control { min-height: 34px; border-radius: 4px; border-color: #cbd5e1; font-size: 13px; }
      .btn { border-radius: 4px; font-weight: 600; }
      .btn-primary { background: #1f6feb; border-color: #1f6feb; }
      .nav-tabs { border-bottom-color: #d9dee7; margin-bottom: 14px; }
      .nav-tabs > li > a { border-radius: 4px 4px 0 0; color: #334155; font-weight: 600; }
      .plot-wrap { background: #ffffff; border: 1px solid #e5e7eb; border-radius: 4px; padding: 8px; margin-bottom: 12px; overflow-x: auto; }
      .status-text { color: #667085; font-size: 12px; margin: 4px 0 10px; }
      .download-row { display: flex; gap: 8px; flex-wrap: wrap; margin: 8px 0 12px; }
      .export-panel { display: grid; grid-template-columns: repeat(4, minmax(100px, 1fr)); gap: 10px; padding: 10px; margin: 8px 0 12px; background: #f8fafc; border: 1px solid #e5e7eb; border-radius: 4px; }
      .export-panel .form-group { margin-bottom: 0; }
      .export-actions { display: flex; gap: 8px; align-items: end; flex-wrap: wrap; }
      table.dataTable { font-size: 12px; }
      @media (max-width: 900px) {
        .app-shell { grid-template-columns: 1fr; }
        .export-panel { grid-template-columns: repeat(2, minmax(100px, 1fr)); }
        .sidebar-panel { position: static; }
      }
    "))
  ),
  div(
    class = "app-shell",
    div(
      class = "sidebar-panel",
      h1("DLBCL Public Dataset Explorer"),
      selectizeInput("gene", "Gene", choices = NULL, selected = default_gene, options = list(placeholder = "Type a gene symbol")),
      selectInput("dataset", "Dataset", choices = dataset_choices, selected = "total"),
      h2("Survival"),
      selectInput("surv_grouping", "Expression split", choices = default_grouping_methods, selected = "median"),
      selectInput("surv_outcome", "Outcome", choices = names(survival_outcomes), selected = "PFS"),
      selectInput("surv_stratification", "Stratification", choices = named_labels(survival_stratifications), selected = "basic"),
      actionButton("run_survival", "Run survival", class = "btn-primary"),
      h2("Boxplot"),
      selectInput("box_comparison", "Comparison", choices = named_labels(boxplot_comparisons), selected = "MYC_RNA"),
      actionButton("run_boxplot", "Run boxplot", class = "btn-primary"),
      h2("Correlation"),
      selectizeInput("target_gene", "Target gene", choices = NULL, selected = default_target, options = list(placeholder = "Type a target gene")),
      actionButton("run_correlation", "Run correlation", class = "btn-primary")
    ),
    div(
      class = "main-panel",
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          "Survival",
          div(class = "status-text", textOutput("survival_status", inline = TRUE)),
          div(class = "plot-wrap", plotOutput("survival_plot", height = "560px")),
          div(class = "plot-wrap", plotOutput("survival_risk_table", height = "260px")),
          div(
            class = "export-panel",
            selectInput("surv_plot_format", "Plot format", choices = c("png", "pdf"), selected = "png"),
            numericInput("surv_plot_width", "KM width (in)", value = 6.5, min = 3, max = 20, step = 0.25),
            numericInput("surv_plot_height", "KM height (in)", value = 5, min = 3, max = 20, step = 0.25),
            numericInput("surv_plot_dpi", "PNG DPI", value = 300, min = 72, max = 600, step = 25),
            numericInput("surv_table_width", "Risk width (in)", value = 6.5, min = 3, max = 20, step = 0.25),
            numericInput("surv_table_height", "Risk height (in)", value = 2.2, min = 1.5, max = 20, step = 0.25),
            div(class = "export-actions", downloadButton("download_survival_plot", "KM plot"), downloadButton("download_survival_risk", "Risk table"))
          ),
          div(class = "download-row", downloadButton("download_survival_summary", "Summary CSV"), downloadButton("download_survival_samples", "Samples CSV")),
          DTOutput("survival_summary")
        ),
        tabPanel(
          "Boxplots",
          div(class = "status-text", textOutput("boxplot_status", inline = TRUE)),
          div(class = "plot-wrap", plotOutput("boxplot_plot", height = "560px")),
          div(
            class = "export-panel",
            selectInput("boxplot_format", "Plot format", choices = c("png", "pdf"), selected = "png"),
            numericInput("boxplot_width", "Width (in)", value = 4.7, min = 3, max = 20, step = 0.25),
            numericInput("boxplot_height", "Height (in)", value = 4.2, min = 3, max = 20, step = 0.25),
            numericInput("boxplot_dpi", "PNG DPI", value = 300, min = 72, max = 600, step = 25),
            div(class = "export-actions", downloadButton("download_boxplot_plot", "Plot"))
          ),
          div(class = "download-row", downloadButton("download_boxplot_summary", "Summary CSV"), downloadButton("download_boxplot_pairwise", "Pairwise CSV")),
          DTOutput("boxplot_summary")
        ),
        tabPanel(
          "Correlations",
          div(class = "status-text", textOutput("correlation_status", inline = TRUE)),
          div(class = "plot-wrap", plotOutput("correlation_plot", height = "560px")),
          div(
            class = "export-panel",
            selectInput("correlation_format", "Plot format", choices = c("png", "pdf"), selected = "png"),
            numericInput("correlation_width", "Width (in)", value = 5.2, min = 3, max = 20, step = 0.25),
            numericInput("correlation_height", "Height (in)", value = 4.6, min = 3, max = 20, step = 0.25),
            numericInput("correlation_dpi", "PNG DPI", value = 300, min = 72, max = 600, step = 25),
            div(class = "export-actions", downloadButton("download_correlation_plot", "Plot"))
          ),
          div(class = "download-row", downloadButton("download_correlation_summary", "Summary CSV"), downloadButton("download_correlation_values", "Values CSV")),
          DTOutput("correlation_summary")
        ),
        tabPanel(
          "Samples",
          tabsetPanel(
            tabPanel("Survival", DTOutput("survival_samples")),
            tabPanel("Boxplot", DTOutput("boxplot_samples")),
            tabPanel("Correlation", DTOutput("correlation_values"))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  updateSelectizeInput(session, "gene", choices = gene_choices, selected = default_gene, server = TRUE)
  updateSelectizeInput(session, "target_gene", choices = gene_choices, selected = default_target, server = TRUE)

  survival_result <- eventReactive(input$run_survival, {
    req(input$gene, nzchar(input$gene))
    withProgress(message = "Running survival analysis", value = 0.2, {
      run_survival_analysis(
        dataset_name = input$dataset,
        gene = toupper(input$gene),
        grouping = input$surv_grouping,
        outcome_name = input$surv_outcome,
        analysis_type = input$surv_stratification,
        data_dir = file.path(project_dir, analysis_data_dir)
      )
    })
  }, ignoreNULL = FALSE)

  boxplot_result <- eventReactive(input$run_boxplot, {
    req(input$gene, nzchar(input$gene))
    withProgress(message = "Running boxplot analysis", value = 0.2, {
      run_boxplot_analysis(
        dataset_name = input$dataset,
        gene = toupper(input$gene),
        comparison_name = input$box_comparison,
        data_dir = file.path(project_dir, analysis_data_dir)
      )
    })
  }, ignoreNULL = FALSE)

  correlation_result <- eventReactive(input$run_correlation, {
    req(input$gene, nzchar(input$gene), input$target_gene, nzchar(input$target_gene))
    withProgress(message = "Running correlation analysis", value = 0.2, {
      run_correlation_analysis(
        dataset_name = input$dataset,
        gene = toupper(input$gene),
        target_gene = toupper(input$target_gene),
        data_dir = file.path(project_dir, analysis_data_dir)
      )
    })
  }, ignoreNULL = FALSE)

  valid_result <- function(result) {
    validate(need(isTRUE(result$valid), result$reason %||% "Analysis could not be completed."))
  }

  write_plot_file <- function(plot, file, format, width, height, dpi) {
    if (format == "pdf") {
      ggsave(file, plot, device = "pdf", width = width, height = height)
    } else {
      ggsave(file, plot, device = "png", width = width, height = height, dpi = dpi)
    }
  }

  plot_filename <- function(prefix, format, ...) {
    paste(safe_name(paste(prefix, ..., sep = "_")), format, sep = ".")
  }

  output$survival_status <- renderText({
    result <- survival_result()
    if (isTRUE(result$valid)) {
      paste("Log-rank", format_pvalue(result$pvalue), "| samples:", nrow(result$sample_groups))
    } else {
      result$reason
    }
  })

  output$survival_plot <- renderPlot({
    result <- survival_result()
    valid_result(result)
    result$plot
  }, res = 120)

  output$survival_risk_table <- renderPlot({
    result <- survival_result()
    valid_result(result)
    result$risk_table
  }, res = 120)

  output$survival_summary <- renderDT({
    result <- survival_result()
    valid_result(result)
    datatable(result$summary, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE))
  })

  output$survival_samples <- renderDT({
    result <- survival_result()
    valid_result(result)
    datatable(result$sample_groups, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE))
  })

  observeEvent(survival_result(), {
    result <- survival_result()
    if (isTRUE(result$valid)) {
      updateNumericInput(session, "surv_table_height", value = result$risk_table_export_height)
    }
  })

  output$boxplot_status <- renderText({
    result <- boxplot_result()
    if (isTRUE(result$valid)) {
      paste(result$test$method, format_pvalue(result$test$pvalue), "| samples:", nrow(result$sample_groups))
    } else {
      result$reason
    }
  })

  output$boxplot_plot <- renderPlot({
    result <- boxplot_result()
    valid_result(result)
    result$plot
  }, res = 120)

  output$boxplot_summary <- renderDT({
    result <- boxplot_result()
    valid_result(result)
    datatable(result$summary, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE))
  })

  output$boxplot_samples <- renderDT({
    result <- boxplot_result()
    valid_result(result)
    datatable(result$sample_groups, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE))
  })

  output$correlation_status <- renderText({
    result <- correlation_result()
    if (isTRUE(result$valid)) {
      spearman <- result$tests[result$tests$method == "Spearman", , drop = FALSE]
      paste("Spearman rho =", format_r(spearman$estimate[[1]]), format_pvalue(spearman$pvalue[[1]]), "| samples:", nrow(result$sample_values))
    } else {
      result$reason
    }
  })

  output$correlation_plot <- renderPlot({
    result <- correlation_result()
    valid_result(result)
    result$plot
  }, res = 120)

  output$correlation_summary <- renderDT({
    result <- correlation_result()
    valid_result(result)
    datatable(result$summary, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE))
  })

  output$correlation_values <- renderDT({
    result <- correlation_result()
    valid_result(result)
    datatable(result$sample_values, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE))
  })

  output$download_survival_summary <- downloadHandler(
    filename = function() paste("survival_summary", input$dataset, toupper(input$gene), input$surv_stratification, input$surv_grouping, input$surv_outcome, "csv", sep = "."),
    content = function(file) write_csv(survival_result()$summary, file)
  )
  output$download_survival_plot <- downloadHandler(
    filename = function() {
      plot_filename("survival_km", input$surv_plot_format, input$dataset, toupper(input$gene), input$surv_stratification, input$surv_grouping, input$surv_outcome)
    },
    content = function(file) {
      result <- survival_result()
      valid_result(result)
      write_plot_file(result$plot, file, input$surv_plot_format, input$surv_plot_width, input$surv_plot_height, input$surv_plot_dpi)
    }
  )
  output$download_survival_risk <- downloadHandler(
    filename = function() {
      plot_filename("survival_risk_table", input$surv_plot_format, input$dataset, toupper(input$gene), input$surv_stratification, input$surv_grouping, input$surv_outcome)
    },
    content = function(file) {
      result <- survival_result()
      valid_result(result)
      write_plot_file(result$risk_table, file, input$surv_plot_format, input$surv_table_width, input$surv_table_height, input$surv_plot_dpi)
    }
  )
  output$download_survival_samples <- downloadHandler(
    filename = function() paste("survival_samples", input$dataset, toupper(input$gene), input$surv_stratification, input$surv_grouping, input$surv_outcome, "csv", sep = "."),
    content = function(file) write_csv(survival_result()$sample_groups, file)
  )
  observeEvent(boxplot_result(), {
    result <- boxplot_result()
    if (isTRUE(result$valid)) {
      updateNumericInput(session, "boxplot_width", value = result$plot_width)
      updateNumericInput(session, "boxplot_height", value = result$plot_height)
    }
  })
  output$download_boxplot_summary <- downloadHandler(
    filename = function() paste("boxplot_summary", input$dataset, toupper(input$gene), input$box_comparison, "csv", sep = "."),
    content = function(file) write_csv(boxplot_result()$summary, file)
  )
  output$download_boxplot_plot <- downloadHandler(
    filename = function() {
      plot_filename("boxplot", input$boxplot_format, input$dataset, toupper(input$gene), input$box_comparison)
    },
    content = function(file) {
      result <- boxplot_result()
      valid_result(result)
      write_plot_file(result$plot, file, input$boxplot_format, input$boxplot_width, input$boxplot_height, input$boxplot_dpi)
    }
  )
  output$download_boxplot_pairwise <- downloadHandler(
    filename = function() paste("boxplot_pairwise", input$dataset, toupper(input$gene), input$box_comparison, "csv", sep = "."),
    content = function(file) write_csv(boxplot_result()$pairwise, file)
  )
  output$download_correlation_summary <- downloadHandler(
    filename = function() paste("correlation_summary", input$dataset, toupper(input$gene), toupper(input$target_gene), "csv", sep = "."),
    content = function(file) write_csv(correlation_result()$summary, file)
  )
  output$download_correlation_plot <- downloadHandler(
    filename = function() {
      plot_filename("correlation", input$correlation_format, input$dataset, toupper(input$gene), toupper(input$target_gene))
    },
    content = function(file) {
      result <- correlation_result()
      valid_result(result)
      write_plot_file(result$plot, file, input$correlation_format, input$correlation_width, input$correlation_height, input$correlation_dpi)
    }
  )
  output$download_correlation_values <- downloadHandler(
    filename = function() paste("correlation_values", input$dataset, toupper(input$gene), toupper(input$target_gene), "csv", sep = "."),
    content = function(file) write_csv(correlation_result()$sample_values, file)
  )
}

shinyApp(ui, server)
