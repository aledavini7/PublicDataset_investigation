# Shared analysis helpers for the curated Sha/REMoDLB dataset.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survival)
  library(survminer)
  library(tibble)
  library(tidyr)
})

analysis_data_dir <- "sha/curated_gene_symbols"
default_genes <- c("CXCR4", "CXCL12", "CD82")
default_grouping_methods <- c("median", "quartile")
default_target_genes <- c("MYC", "BCL2", "BCL6", "CXCR4", "CXCL12", "CD82")
drop_values <- c("", "n/a", "NA", "NaN")

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

analysis_datasets <- function(data_dir = analysis_data_dir) {
  list(
    total = list(
      expression = file.path(data_dir, "sha_total_expression.rds"),
      clinical = file.path(data_dir, "clin_sha_tot.rds")
    ),
    rchop = list(
      expression = file.path(data_dir, "sha_rchop_expression.rds"),
      clinical = file.path(data_dir, "clin_sha_rchop.rds")
    ),
    rbchop = list(
      expression = file.path(data_dir, "sha_rbchop_expression.rds"),
      clinical = file.path(data_dir, "clin_sha_rbchop.rds")
    )
  )
}

survival_outcomes <- list(
  OS = list(time = "Mesi_OS", event = "Evento_OS", xlab = "Time (months)", break_time_by = 12, axis_step = 6, axis_padding = 2, axis_left_padding = 5),
  PFS = list(time = "Mesi_PFS", event = "Evento_PFS", xlab = "Time (months)", break_time_by = 12, axis_step = 6, axis_padding = 2, axis_left_padding = 5)
)

survival_stratifications <- list(
  basic = list(type = "expression_only", label = "Gene expression"),
  myc_expression = list(type = "computed", column = "MYC_expression_group", label = "MYC expression median"),
  myc_rna = list(type = "metadata", column = "MYC_RNA", label = "MYC RNA", order = c("average", "high")),
  bcl2_rna = list(type = "metadata", column = "BCL2_RNA", label = "BCL2 RNA", order = c("average", "high")),
  myc_rearrangement = list(
    type = "metadata",
    column = "MYC_rearrangement",
    label = "MYC rearrangement",
    order = c("not-rearranged", "rearranged"),
    labels = c("not-rearranged" = "WT", "rearranged" = "REARRANGED")
  ),
  bcl2_rearrangement = list(
    type = "metadata",
    column = "BCL2_rearrangement",
    label = "BCL2 rearrangement",
    order = c("not-rearranged", "rearranged"),
    labels = c("not-rearranged" = "WT", "rearranged" = "REARRANGED")
  ),
  bcl6_rearrangement = list(
    type = "metadata",
    column = "BCL6_rearrangement",
    label = "BCL6 rearrangement",
    order = c("not-rearranged", "rearranged"),
    labels = c("not-rearranged" = "WT", "rearranged" = "REARRANGED")
  ),
  hit_status = list(
    type = "metadata",
    column = "Hit_rearrangement",
    label = "Hit rearrangement",
    order = c("MYC-normal", "single-hit", "double-hit", "MYC-rearranged NOS")
  ),
  coo_class = list(type = "metadata", column = "COO_class", label = "COO class", order = c("GCB", "ABC", "UNC")),
  expressor_rna = list(
    type = "metadata",
    column = "expressor_RNA",
    label = "RNA double-expressor",
    order = c("non-double-expressor", "double-expressor"),
    labels = c("non-double-expressor" = "NON_DOUBLE", "double-expressor" = "DOUBLE")
  ),
  expressor_ihc = list(
    type = "metadata",
    column = "expressor_IHC",
    label = "IHC double-expressor",
    order = c("non-double-expressor", "double-expressor"),
    labels = c("non-double-expressor" = "NON_DOUBLE", "double-expressor" = "DOUBLE")
  ),
  myc_ihc = list(
    type = "metadata",
    column = "MYC_IHC",
    label = "MYC IHC",
    order = c("negative", "borderline negative", "positive"),
    labels = c("negative" = "NEGATIVE", "borderline negative" = "BORDERLINE_NEGATIVE", "positive" = "POSITIVE")
  ),
  bcl2_ihc = list(
    type = "metadata",
    column = "BCL2_IHC",
    label = "BCL2 IHC",
    order = c("negative", "positive"),
    labels = c("negative" = "NEGATIVE", "positive" = "POSITIVE")
  )
)

boxplot_comparisons <- list(
  MYC_RNA = list(column = "MYC_RNA", title = "MYC high vs average", order = c("average", "high")),
  MYC_rearrangement = list(
    column = "MYC_rearrangement",
    title = "MYC rearranged vs WT",
    order = c("not-rearranged", "rearranged"),
    labels = c("not-rearranged" = "WT", "rearranged" = "REARRANGED")
  ),
  BCL2_RNA = list(column = "BCL2_RNA", title = "BCL2 high vs average", order = c("average", "high")),
  BCL2_rearrangement = list(
    column = "BCL2_rearrangement",
    title = "BCL2 rearranged vs WT",
    order = c("not-rearranged", "rearranged"),
    labels = c("not-rearranged" = "WT", "rearranged" = "REARRANGED")
  ),
  BCL6_rearrangement = list(
    column = "BCL6_rearrangement",
    title = "BCL6 rearranged vs WT",
    order = c("not-rearranged", "rearranged"),
    labels = c("not-rearranged" = "WT", "rearranged" = "REARRANGED")
  ),
  Hit_rearrangement = list(
    column = "Hit_rearrangement",
    title = "Hit rearrangement status",
    order = c("MYC-normal", "single-hit", "double-hit", "MYC-rearranged NOS")
  ),
  COO_class = list(column = "COO_class", title = "COO class", order = c("GCB", "ABC", "UNC")),
  expressor_RNA = list(
    column = "expressor_RNA",
    title = "RNA double-expressor status",
    order = c("non-double-expressor", "double-expressor")
  ),
  expressor_IHC = list(
    column = "expressor_IHC",
    title = "IHC double-expressor status",
    order = c("non-double-expressor", "double-expressor")
  ),
  MYC_IHC = list(column = "MYC_IHC", title = "MYC IHC status", order = c("negative", "borderline negative", "positive")),
  BCL2_IHC = list(column = "BCL2_IHC", title = "BCL2 IHC status", order = c("negative", "positive"))
)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("^_|_$", "", x)
}

format_pvalue <- function(pvalue) {
  if (is.na(pvalue)) {
    return("p = NA")
  }
  if (pvalue < 0.0001) {
    "p < 0.0001"
  } else {
    paste0("p = ", signif(pvalue, 2))
  }
}

pvalue_stars <- function(pvalue) {
  if (is.na(pvalue)) {
    "ns"
  } else if (pvalue <= 0.0001) {
    "****"
  } else if (pvalue <= 0.001) {
    "***"
  } else if (pvalue <= 0.01) {
    "**"
  } else if (pvalue <= 0.05) {
    "*"
  } else {
    "ns"
  }
}

format_r <- function(r_value) {
  if (is.na(r_value)) {
    return("NA")
  }
  sprintf("%.2f", r_value)
}

load_dataset <- function(dataset_name, data_dir = analysis_data_dir) {
  datasets <- analysis_datasets(data_dir)
  if (!dataset_name %in% names(datasets)) {
    stop("Unknown dataset: ", dataset_name)
  }

  paths <- datasets[[dataset_name]]
  expression <- readRDS(paths$expression)
  clinical <- readRDS(paths$clinical)
  clinical$Sample_ID <- as.character(clinical$Sample_ID)

  if (!all(colnames(expression) %in% clinical$Sample_ID)) {
    missing <- setdiff(colnames(expression), clinical$Sample_ID)
    stop("Clinical metadata missing sample(s): ", paste(missing, collapse = ", "))
  }

  clinical <- clinical[match(colnames(expression), clinical$Sample_ID), , drop = FALSE]

  if (!identical(colnames(expression), clinical$Sample_ID)) {
    stop("Could not align clinical metadata to expression columns.")
  }

  list(expression = expression, clinical = clinical)
}

available_genes <- function(dataset_name = "total", data_dir = analysis_data_dir) {
  rownames(load_dataset(dataset_name, data_dir)$expression)
}

make_output_dirs <- function(output_dir, gene, dataset_name, analysis_name) {
  paths <- list(
    plots = file.path(output_dir, "plots", gene, dataset_name, analysis_name),
    tables = file.path(output_dir, "tables", gene, dataset_name, analysis_name)
  )
  dir.create(paths$plots, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$tables, showWarnings = FALSE, recursive = TRUE)
  paths
}

make_gene_data <- function(expression, clinical, gene) {
  if (!gene %in% rownames(expression)) {
    stop("Gene not found in expression matrix: ", gene)
  }
  if (!"MYC" %in% rownames(expression)) {
    stop("MYC not found in expression matrix.")
  }

  gene_expression <- as.numeric(expression[gene, clinical$Sample_ID])
  myc_expression <- as.numeric(expression["MYC", clinical$Sample_ID])
  myc_cutoff <- median(myc_expression, na.rm = TRUE)

  clinical %>%
    mutate(
      gene = gene,
      expression = gene_expression,
      MYC_expression = myc_expression,
      MYC_expression_cutoff = myc_cutoff,
      MYC_expression_group = if_else(MYC_expression >= myc_cutoff, "MYC_HIGH", "MYC_LOW")
    )
}

add_expression_group <- function(data, method) {
  if (method == "median") {
    cutoff <- median(data$expression, na.rm = TRUE)
    data %>%
      mutate(
        grouping = method,
        cutoff_low = cutoff,
        cutoff_high = cutoff,
        expression_group = if_else(expression >= cutoff, "HIGH", "LOW")
      )
  } else if (method == "quartile") {
    q1 <- as.numeric(quantile(data$expression, 0.25, na.rm = TRUE, names = FALSE))
    q3 <- as.numeric(quantile(data$expression, 0.75, na.rm = TRUE, names = FALSE))
    data %>%
      mutate(
        grouping = method,
        cutoff_low = q1,
        cutoff_high = q3,
        expression_group = case_when(
          expression <= q1 ~ "LOW",
          expression >= q3 ~ "HIGH",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(expression_group))
  } else {
    stop("Unknown grouping method: ", method)
  }
}

make_survival_data <- function(grouped, analysis_type, min_modifier_n = 5) {
  if (!analysis_type %in% names(survival_stratifications)) {
    stop("Unknown survival stratification: ", analysis_type)
  }

  stratification <- survival_stratifications[[analysis_type]]

  if (stratification$type == "expression_only") {
    return(
      grouped %>%
        mutate(
          modifier = "none",
          modifier_group = "none",
          modifier_group_raw = "none",
          survival_group = expression_group
        )
    )
  }

  if (!stratification$column %in% colnames(grouped)) {
    stop("Clinical/derived column not found: ", stratification$column)
  }

  grouped %>%
    mutate(modifier_group_raw = as.character(.data[[stratification$column]])) %>%
    filter(
      !is.na(modifier_group_raw),
      !tolower(modifier_group_raw) %in% tolower(drop_values)
    ) %>%
    mutate(
      modifier_group = if ("labels" %in% names(stratification)) {
        dplyr::recode(modifier_group_raw, !!!stratification$labels, .default = toupper(modifier_group_raw))
      } else {
        toupper(modifier_group_raw)
      },
      modifier_order = if ("order" %in% names(stratification)) {
        match(modifier_group_raw, stratification$order)
      } else {
        match(modifier_group_raw, unique(modifier_group_raw))
      }
    ) %>%
    filter(!is.na(modifier_order)) %>%
    add_count(modifier_group, name = "modifier_n") %>%
    filter(modifier_n >= min_modifier_n) %>%
    arrange(modifier_order, modifier_group_raw, expression_group) %>%
    mutate(
      modifier = stratification$column,
      modifier_group = factor(modifier_group, levels = unique(modifier_group)),
      survival_group = paste(as.character(modifier_group), expression_group, sep = "_")
    ) %>%
    select(-modifier_order, -modifier_n)
}

fit_survival <- function(data, outcome) {
  analysis_data <- data %>%
    filter(
      !is.na(.data[[outcome$time]]),
      !is.na(.data[[outcome$event]]),
      !is.na(survival_group)
    ) %>%
    mutate(survival_group = factor(survival_group, levels = unique(survival_group)))

  if (n_distinct(analysis_data$survival_group) < 2) {
    return(NULL)
  }

  surv_formula <- as.formula(paste0("Surv(", outcome$time, ", ", outcome$event, ") ~ survival_group"))
  fit <- survfit(surv_formula, data = analysis_data)
  fit$call$formula <- surv_formula
  logrank <- survdiff(surv_formula, data = analysis_data)
  pvalue <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)

  list(data = analysis_data, formula = surv_formula, fit = fit, pvalue = pvalue)
}

summarise_survival_fit <- function(fit_result, dataset_name, gene, grouping, outcome_name, analysis_type) {
  data <- fit_result$data
  outcome_def <- survival_outcomes[[outcome_name]]

  group_counts <- data %>%
    count(survival_group, name = "n_group") %>%
    mutate(survival_group = as.character(survival_group))

  event_counts <- data %>%
    group_by(survival_group) %>%
    summarise(events = sum(.data[[outcome_def$event]] == 1), .groups = "drop") %>%
    mutate(survival_group = as.character(survival_group))

  med <- surv_median(fit_result$fit)
  median_stats <- setNames(med$median, sub("^survival_group=", "", med$strata))

  tibble(
    dataset = dataset_name,
    gene = gene,
    analysis_type = analysis_type,
    grouping = grouping,
    outcome = outcome_name,
    pvalue_logrank = fit_result$pvalue,
    n_total = nrow(data),
    cutoff_low = unique(data$cutoff_low),
    cutoff_high = unique(data$cutoff_high),
    myc_expression_cutoff = unique(data$MYC_expression_cutoff),
    max_time = max(data[[outcome_def$time]], na.rm = TRUE)
  ) %>%
    crossing(survival_group = as.character(levels(data$survival_group))) %>%
    left_join(group_counts, by = "survival_group") %>%
    left_join(event_counts, by = "survival_group") %>%
    mutate(median_survival = unname(median_stats[survival_group]))
}

survival_palette <- function(n) {
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442")
  if (n <= length(base)) {
    base[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, palette = "Dark 3")
  }
}

survival_group_palette <- function(groups) {
  groups <- as.character(groups)
  base_low <- "#0072B2"
  base_high <- "#D55E00"
  modifier_low <- c("#56B4E9", "#0072B2", "#8EC9F0", "#004C99", "#78A9CF", "#2B6F9E")
  modifier_high <- c("#E69F00", "#D55E00", "#F0C36A", "#A64000", "#E07B39", "#B86B00")

  if (setequal(groups, c("LOW", "HIGH"))) {
    return(unname(c(LOW = base_low, HIGH = base_high)[groups]))
  }

  expression_part <- ifelse(grepl("_HIGH$", groups), "HIGH", ifelse(grepl("_LOW$", groups), "LOW", NA_character_))
  modifier_part <- sub("_(HIGH|LOW)$", "", groups)
  modifiers <- unique(modifier_part[!is.na(expression_part)])

  if (length(modifiers) > 0 && all(!is.na(expression_part))) {
    colors <- vapply(seq_along(groups), function(i) {
      modifier_index <- match(modifier_part[[i]], modifiers)
      if (expression_part[[i]] == "LOW") {
        modifier_low[[((modifier_index - 1) %% length(modifier_low)) + 1]]
      } else {
        modifier_high[[((modifier_index - 1) %% length(modifier_high)) + 1]]
      }
    }, character(1))
    return(unname(colors))
  }

  survival_palette(length(groups))
}

survival_axis_config <- function(outcome, max_time) {
  axis_step <- outcome$axis_step %||% 6
  axis_padding <- outcome$axis_padding %||% 2
  axis_left_padding <- outcome$axis_left_padding %||% 5
  x_axis_max <- ceiling((max_time + axis_padding) / axis_step) * axis_step
  x_axis_min <- -min(axis_left_padding, x_axis_max * 0.065)
  list(
    x_axis_min = x_axis_min,
    x_axis_max = x_axis_max,
    breaks = seq(0, x_axis_max, by = outcome$break_time_by %||% 12),
    xlab = outcome$xlab %||% "Time (months)"
  )
}

make_survival_plot <- function(fit_result, dataset_name, gene, grouping, outcome_name, analysis_type) {
  levels_group <- levels(fit_result$data$survival_group)
  n_groups <- length(levels_group)
  outcome <- survival_outcomes[[outcome_name]]
  max_time <- max(fit_result$data[[outcome$time]], na.rm = TRUE)
  axis <- survival_axis_config(outcome, max_time)
  risk_table_height <- min(0.45, max(0.24, 0.08 + 0.055 * n_groups))
  risk_table_font_size <- dplyr::case_when(
    n_groups <= 2 ~ 4.8,
    n_groups <= 4 ~ 4.0,
    TRUE ~ 3.4
  )

  plot <- suppressMessages(suppressWarnings(ggsurvplot(
    fit_result$fit,
    data = fit_result$data,
    conf.int = FALSE,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.title = "",
    risk.table.height = risk_table_height,
    risk.table.fontsize = risk_table_font_size,
    xlim = c(axis$x_axis_min, axis$x_axis_max),
    break.time.by = outcome$break_time_by %||% 12,
    xlab = axis$xlab,
    ylab = paste0(outcome_name, " probability"),
    title = paste(gene, outcome_name, dataset_name, analysis_type, grouping, sep = " - "),
    legend.title = "Group",
    legend.labs = levels_group,
    palette = survival_group_palette(levels_group),
    ggtheme = theme_bw(base_size = 11) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)),
    risk.table.y.text = FALSE,
    tables.theme = theme_bw(base_size = 10)
  )))

  suppressMessages(suppressWarnings({
    plot$plot <- plot$plot +
      scale_x_continuous(limits = c(axis$x_axis_min, axis$x_axis_max), breaks = axis$breaks, expand = expansion(mult = c(0.005, 0.015))) +
      scale_y_continuous(limits = c(0, 1.03), expand = expansion(mult = c(0.015, 0.035))) +
      coord_cartesian(xlim = c(axis$x_axis_min, axis$x_axis_max), ylim = c(0, 1.03), clip = "off") +
      theme(plot.title = element_text(size = 11, face = "bold"))
    plot$table <- plot$table +
      scale_x_continuous(limits = c(axis$x_axis_min, axis$x_axis_max), breaks = axis$breaks, expand = expansion(mult = c(0.005, 0.015))) +
      coord_cartesian(xlim = c(axis$x_axis_min, axis$x_axis_max), clip = "off")
  }))

  list(
    plot = plot$plot,
    risk_table = plot$table,
    risk_table_export_height = max(2, 1.1 + 0.42 * n_groups)
  )
}

survival_sample_groups <- function(fit_result, dataset_name, gene, grouping, outcome_name, analysis_type) {
  fit_result$data %>%
    transmute(
      dataset = dataset_name,
      gene = gene,
      analysis_type = analysis_type,
      grouping = grouping,
      outcome = outcome_name,
      Sample_ID,
      expression,
      expression_group,
      MYC_expression,
      MYC_expression_group,
      MYC_RNA,
      BCL2_RNA,
      MYC_rearrangement,
      BCL2_rearrangement,
      BCL6_rearrangement,
      Hit_rearrangement,
      COO_class,
      expressor_RNA,
      expressor_IHC,
      MYC_IHC,
      BCL2_IHC,
      modifier,
      modifier_group,
      modifier_group_raw,
      survival_group,
      cutoff_low,
      cutoff_high,
      MYC_expression_cutoff
    )
}

save_survival_outputs <- function(plots, summary_rows, sample_groups, dataset_name, gene, grouping, outcome_name, analysis_type, output_dir) {
  paths <- make_output_dirs(output_dir, gene, dataset_name, analysis_type)
  filename_base <- safe_name(paste(dataset_name, gene, analysis_type, grouping, outcome_name, sep = "_"))

  plot_pdf <- file.path(paths$plots, paste0(filename_base, ".pdf"))
  plot_png <- file.path(paths$plots, paste0(filename_base, ".png"))
  risk_table_pdf <- file.path(paths$plots, paste0(filename_base, "_risk_table.pdf"))
  risk_table_png <- file.path(paths$plots, paste0(filename_base, "_risk_table.png"))
  summary_csv <- file.path(paths$tables, paste0(filename_base, "_summary.csv"))
  sample_groups_csv <- file.path(paths$tables, paste0(filename_base, "_sample_groups.csv"))

  suppressMessages(suppressWarnings({
    ggsave(plot_pdf, plots$plot, device = "pdf", width = 6.5, height = 5)
    ggsave(plot_png, plots$plot, device = "png", width = 6.5, height = 5, dpi = 300)
    ggsave(risk_table_pdf, plots$risk_table, device = "pdf", width = 6.5, height = plots$risk_table_export_height)
    ggsave(risk_table_png, plots$risk_table, device = "png", width = 6.5, height = plots$risk_table_export_height, dpi = 300)
  }))
  write_csv(summary_rows, summary_csv)
  write_csv(sample_groups, sample_groups_csv)

  c(
    plot_pdf = plot_pdf,
    plot_png = plot_png,
    risk_table_pdf = risk_table_pdf,
    risk_table_png = risk_table_png,
    summary_csv = summary_csv,
    sample_groups_csv = sample_groups_csv
  )
}

run_survival_analysis <- function(dataset_name, gene, grouping, outcome_name, analysis_type,
                                  data_dir = analysis_data_dir, min_modifier_n = 5,
                                  save_outputs = FALSE, output_dir = "results/gene_survival") {
  dataset <- load_dataset(dataset_name, data_dir)
  gene_data <- make_gene_data(dataset$expression, dataset$clinical, gene)
  grouped <- add_expression_group(gene_data, grouping)
  analysis_data <- make_survival_data(grouped, analysis_type, min_modifier_n)
  fit_result <- fit_survival(analysis_data, survival_outcomes[[outcome_name]])

  if (is.null(fit_result)) {
    return(list(valid = FALSE, reason = "Fewer than two survival groups after filtering"))
  }

  summary_rows <- summarise_survival_fit(fit_result, dataset_name, gene, grouping, outcome_name, analysis_type)
  sample_groups <- survival_sample_groups(fit_result, dataset_name, gene, grouping, outcome_name, analysis_type)
  plots <- make_survival_plot(fit_result, dataset_name, gene, grouping, outcome_name, analysis_type)
  output_paths <- NULL

  if (save_outputs) {
    output_paths <- save_survival_outputs(plots, summary_rows, sample_groups, dataset_name, gene, grouping, outcome_name, analysis_type, output_dir)
  }

  list(
    valid = TRUE,
    data = fit_result$data,
    fit = fit_result$fit,
    pvalue = fit_result$pvalue,
    summary = summary_rows,
    sample_groups = sample_groups,
    plot = plots$plot,
    risk_table = plots$risk_table,
    risk_table_export_height = plots$risk_table_export_height,
    output_paths = output_paths
  )
}

make_boxplot_data <- function(expression, clinical, gene, comparison) {
  if (!gene %in% rownames(expression)) {
    stop("Gene not found in expression matrix: ", gene)
  }
  if (!comparison$column %in% colnames(clinical)) {
    stop("Clinical column not found: ", comparison$column)
  }

  clinical %>%
    transmute(
      Sample_ID,
      expression = as.numeric(expression[gene, Sample_ID]),
      group_raw = as.character(.data[[comparison$column]])
    ) %>%
    filter(
      !is.na(expression),
      !is.na(group_raw),
      !tolower(group_raw) %in% tolower(drop_values)
    ) %>%
    mutate(
      group = if ("labels" %in% names(comparison)) {
        dplyr::recode(group_raw, !!!comparison$labels, .default = group_raw)
      } else {
        toupper(group_raw)
      },
      group_order = if ("order" %in% names(comparison)) {
        match(group_raw, comparison$order)
      } else {
        match(group_raw, unique(group_raw))
      }
    ) %>%
    arrange(group_order, group_raw) %>%
    mutate(group = factor(group, levels = unique(group)))
}

validate_boxplot_data <- function(plot_data, min_group_n = 5) {
  group_counts <- plot_data %>%
    count(group, name = "n_group") %>%
    mutate(group = as.character(group))

  usable_groups <- group_counts %>%
    filter(n_group >= min_group_n)

  if (nrow(usable_groups) < 2) {
    return(list(valid = FALSE, reason = "Fewer than two groups with enough samples."))
  }

  plot_data <- plot_data %>%
    filter(as.character(group) %in% usable_groups$group) %>%
    mutate(group = factor(as.character(group), levels = usable_groups$group))

  list(valid = TRUE, data = plot_data, group_counts = usable_groups)
}

run_association_test <- function(plot_data) {
  n_groups <- n_distinct(plot_data$group)
  if (n_groups == 2) {
    test <- wilcox.test(expression ~ group, data = plot_data, exact = FALSE)
    list(method = "Wilcoxon rank-sum test", pvalue = unname(test$p.value))
  } else {
    test <- kruskal.test(expression ~ group, data = plot_data)
    list(method = "Kruskal-Wallis test", pvalue = unname(test$p.value))
  }
}

run_pairwise_tests <- function(plot_data, dataset_name, gene, comparison_name) {
  groups <- levels(plot_data$group)
  pairs <- t(combn(groups, 2))

  lapply(seq_len(nrow(pairs)), function(i) {
    group_1 <- pairs[i, 1]
    group_2 <- pairs[i, 2]
    test_data <- plot_data %>%
      filter(group %in% c(group_1, group_2)) %>%
      mutate(group = droplevels(group))
    test <- wilcox.test(expression ~ group, data = test_data, exact = FALSE)

    tibble(
      dataset = dataset_name,
      gene = gene,
      comparison = comparison_name,
      group_1 = group_1,
      group_2 = group_2,
      n_group_1 = sum(test_data$group == group_1),
      n_group_2 = sum(test_data$group == group_2),
      test_method = "Pairwise Wilcoxon rank-sum test",
      pvalue = unname(test$p.value)
    )
  }) %>%
    bind_rows() %>%
    mutate(
      pvalue_adj_method = "BH",
      pvalue_adj = p.adjust(pvalue, method = "BH"),
      pvalue_label = vapply(pvalue_adj, format_pvalue, character(1)),
      pvalue_stars = vapply(pvalue_adj, pvalue_stars, character(1))
    )
}

summarise_boxplot_groups <- function(plot_data, dataset_name, gene, comparison_name, test_result) {
  plot_data %>%
    group_by(group) %>%
    summarise(
      n_group = n(),
      median_expression = median(expression),
      mean_expression = mean(expression),
      sd_expression = sd(expression),
      q1_expression = quantile(expression, 0.25, names = FALSE),
      q3_expression = quantile(expression, 0.75, names = FALSE),
      min_expression = min(expression),
      max_expression = max(expression),
      .groups = "drop"
    ) %>%
    mutate(
      dataset = dataset_name,
      gene = gene,
      comparison = comparison_name,
      test_method = test_result$method,
      pvalue = test_result$pvalue,
      pvalue_label = format_pvalue(test_result$pvalue)
    ) %>%
    select(dataset, gene, comparison, group, n_group, everything())
}

boxplot_palette <- function(n) {
  base <- c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF", "#E69F00", "#56B4E9", "#D55E00")
  if (n <= length(base)) {
    base[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, palette = "Dark 3")
  }
}

make_boxplot <- function(plot_data, dataset_name, gene, comparison_name, comparison, pairwise_results) {
  group_counts <- plot_data %>%
    count(group, name = "n_group") %>%
    mutate(
      group = as.character(group),
      x_label = paste0(toupper(group), "\n(N=", n_group, ")")
    )

  plot_data <- plot_data %>%
    mutate(
      group_chr = as.character(group),
      x_label = group_counts$x_label[match(group_chr, group_counts$group)],
      x_label = factor(x_label, levels = group_counts$x_label)
    )

  n_groups <- nrow(group_counts)
  n_pairwise <- nrow(pairwise_results)
  y_min <- min(plot_data$expression, na.rm = TRUE)
  y_max <- max(plot_data$expression, na.rm = TRUE)
  y_range <- max(y_max - y_min, 0.5)
  bracket_start <- y_max + 0.08 * y_range
  bracket_step <- 0.11 * y_range
  y_upper <- y_max + (0.28 + 0.11 * max(n_pairwise - 1, 0)) * y_range

  p <- ggplot(plot_data, aes(x = x_label, y = expression, color = x_label)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.45, fill = NA) +
    geom_jitter(width = 0.18, size = 1.1, alpha = 0.65) +
    scale_color_manual(values = boxplot_palette(n_groups), guide = "none") +
    labs(
      title = comparison$title,
      subtitle = toupper(dataset_name),
      x = NULL,
      y = paste0(gene, " expression level")
    ) +
    coord_cartesian(ylim = c(y_min - 0.03 * y_range, y_upper), clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "italic", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 10),
      plot.margin = margin(8, 14, 8, 8)
    )

  if (n_pairwise > 0) {
    pairwise_plot <- pairwise_results %>%
      mutate(
        x_start = match(group_1, group_counts$group),
        x_end = match(group_2, group_counts$group),
        bracket_y = bracket_start + (row_number() - 1) * bracket_step,
        text_y = bracket_y + 0.045 * y_range
      )

    for (i in seq_len(nrow(pairwise_plot))) {
      row <- pairwise_plot[i, ]
      p <- p +
        annotate("segment", x = row$x_start, xend = row$x_end, y = row$bracket_y, yend = row$bracket_y, linewidth = 0.35) +
        annotate("segment", x = row$x_start, xend = row$x_start, y = row$bracket_y - 0.03 * y_range, yend = row$bracket_y, linewidth = 0.35) +
        annotate("segment", x = row$x_end, xend = row$x_end, y = row$bracket_y - 0.03 * y_range, yend = row$bracket_y, linewidth = 0.35) +
        annotate("text", x = (row$x_start + row$x_end) / 2, y = row$text_y, label = row$pvalue_stars, size = 3.2)
    }
  }

  list(
    plot = p,
    width = max(4, 1.15 * n_groups + 2.4),
    height = max(4.2, 4.2 + 0.18 * max(n_pairwise - 1, 0))
  )
}

boxplot_sample_groups <- function(plot_data, dataset_name, gene, comparison_name) {
  plot_data %>%
    transmute(
      dataset = dataset_name,
      gene = gene,
      comparison = comparison_name,
      Sample_ID,
      expression,
      group = as.character(group)
    )
}

save_boxplot_outputs <- function(plot_result, summary_rows, pairwise_rows, sample_groups, dataset_name, gene, comparison_name, output_dir) {
  paths <- make_output_dirs(output_dir, gene, dataset_name, comparison_name)
  filename_base <- safe_name(paste(dataset_name, gene, comparison_name, sep = "_"))
  plot_pdf <- file.path(paths$plots, paste0(filename_base, ".pdf"))
  plot_png <- file.path(paths$plots, paste0(filename_base, ".png"))
  summary_csv <- file.path(paths$tables, paste0(filename_base, "_summary.csv"))
  pairwise_csv <- file.path(paths$tables, paste0(filename_base, "_pairwise.csv"))
  sample_groups_csv <- file.path(paths$tables, paste0(filename_base, "_sample_groups.csv"))

  ggsave(plot_pdf, plot_result$plot, device = "pdf", width = plot_result$width, height = plot_result$height)
  ggsave(plot_png, plot_result$plot, device = "png", width = plot_result$width, height = plot_result$height, dpi = 300)
  write_csv(summary_rows, summary_csv)
  write_csv(pairwise_rows, pairwise_csv)
  write_csv(sample_groups, sample_groups_csv)

  c(plot_pdf = plot_pdf, plot_png = plot_png, summary_csv = summary_csv, pairwise_csv = pairwise_csv, sample_groups_csv = sample_groups_csv)
}

run_boxplot_analysis <- function(dataset_name, gene, comparison_name,
                                 data_dir = analysis_data_dir, min_group_n = 5,
                                 save_outputs = FALSE, output_dir = "results/gene_boxplots") {
  if (!comparison_name %in% names(boxplot_comparisons)) {
    stop("Unknown boxplot comparison: ", comparison_name)
  }

  dataset <- load_dataset(dataset_name, data_dir)
  comparison <- boxplot_comparisons[[comparison_name]]
  plot_data <- make_boxplot_data(dataset$expression, dataset$clinical, gene, comparison)
  validation <- validate_boxplot_data(plot_data, min_group_n)

  if (!validation$valid) {
    return(list(valid = FALSE, reason = validation$reason))
  }

  plot_data <- validation$data
  test_result <- run_association_test(plot_data)
  pairwise_rows <- run_pairwise_tests(plot_data, dataset_name, gene, comparison_name)
  summary_rows <- summarise_boxplot_groups(plot_data, dataset_name, gene, comparison_name, test_result)
  plot_result <- make_boxplot(plot_data, dataset_name, gene, comparison_name, comparison, pairwise_rows)
  sample_groups <- boxplot_sample_groups(plot_data, dataset_name, gene, comparison_name)
  output_paths <- NULL

  if (save_outputs) {
    output_paths <- save_boxplot_outputs(plot_result, summary_rows, pairwise_rows, sample_groups, dataset_name, gene, comparison_name, output_dir)
  }

  list(
    valid = TRUE,
    data = plot_data,
    test = test_result,
    summary = summary_rows,
    pairwise = pairwise_rows,
    sample_groups = sample_groups,
    plot = plot_result$plot,
    plot_width = plot_result$width,
    plot_height = plot_result$height,
    output_paths = output_paths
  )
}

make_correlation_data <- function(expression, clinical, gene, target_gene) {
  missing_genes <- setdiff(c(gene, target_gene), rownames(expression))
  if (length(missing_genes) > 0) {
    stop("Gene(s) not found in expression matrix: ", paste(missing_genes, collapse = ", "))
  }

  clinical %>%
    transmute(
      Sample_ID,
      gene_expression = as.numeric(expression[gene, Sample_ID]),
      target_expression = as.numeric(expression[target_gene, Sample_ID])
    ) %>%
    filter(!is.na(gene_expression), !is.na(target_expression))
}

run_correlation_tests <- function(correlation_data) {
  pearson <- cor.test(correlation_data$gene_expression, correlation_data$target_expression, method = "pearson")
  spearman <- cor.test(correlation_data$gene_expression, correlation_data$target_expression, method = "spearman", exact = FALSE)

  tibble(
    method = c("Pearson", "Spearman"),
    estimate = c(unname(pearson$estimate), unname(spearman$estimate)),
    pvalue = c(unname(pearson$p.value), unname(spearman$p.value)),
    pvalue_label = vapply(c(unname(pearson$p.value), unname(spearman$p.value)), format_pvalue, character(1))
  )
}

summarise_correlation <- function(correlation_data, test_rows, dataset_name, gene, target_gene) {
  test_rows %>%
    mutate(
      dataset = dataset_name,
      gene = gene,
      target_gene = target_gene,
      n_samples = nrow(correlation_data),
      gene_median = median(correlation_data$gene_expression),
      target_median = median(correlation_data$target_expression),
      gene_mean = mean(correlation_data$gene_expression),
      target_mean = mean(correlation_data$target_expression),
      gene_sd = sd(correlation_data$gene_expression),
      target_sd = sd(correlation_data$target_expression)
    ) %>%
    select(dataset, gene, target_gene, n_samples, method, estimate, pvalue, pvalue_label, everything())
}

make_correlation_plot <- function(correlation_data, test_rows, dataset_name, gene, target_gene) {
  spearman <- test_rows %>% filter(method == "Spearman") %>% slice(1)
  pearson <- test_rows %>% filter(method == "Pearson") %>% slice(1)
  label <- paste0(
    "Spearman rho = ", format_r(spearman$estimate), ", ", format_pvalue(spearman$pvalue),
    " | Pearson r = ", format_r(pearson$estimate), ", ", format_pvalue(pearson$pvalue)
  )

  x_min <- min(correlation_data$target_expression)
  x_max <- max(correlation_data$target_expression)
  y_min <- min(correlation_data$gene_expression)
  y_max <- max(correlation_data$gene_expression)
  x_range <- max(x_max - x_min, 0.5)
  y_range <- max(y_max - y_min, 0.5)

  ggplot(correlation_data, aes(x = target_expression, y = gene_expression)) +
    geom_point(color = "#0072B2", size = 1.5, alpha = 0.65) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#D55E00", fill = "#D55E00", alpha = 0.16, linewidth = 0.7) +
    labs(
      title = paste(gene, "vs", target_gene),
      subtitle = paste(toupper(dataset_name), label, sep = "\n"),
      x = paste0(target_gene, " expression level"),
      y = paste0(gene, " expression level")
    ) +
    coord_cartesian(
      xlim = c(x_min - 0.04 * x_range, x_max + 0.04 * x_range),
      ylim = c(y_min - 0.04 * y_range, y_max + 0.08 * y_range),
      clip = "off"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5, lineheight = 1.1),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10),
      plot.margin = margin(10, 12, 8, 8)
    )
}

correlation_sample_values <- function(correlation_data, dataset_name, gene, target_gene) {
  correlation_data %>%
    transmute(
      dataset = dataset_name,
      gene = gene,
      target_gene = target_gene,
      Sample_ID,
      gene_expression,
      target_expression
    )
}

save_correlation_outputs <- function(plot, summary_rows, sample_values, dataset_name, gene, target_gene, output_dir) {
  paths <- make_output_dirs(output_dir, gene, dataset_name, target_gene)
  filename_base <- safe_name(paste(dataset_name, gene, target_gene, sep = "_"))
  plot_pdf <- file.path(paths$plots, paste0(filename_base, ".pdf"))
  plot_png <- file.path(paths$plots, paste0(filename_base, ".png"))
  summary_csv <- file.path(paths$tables, paste0(filename_base, "_summary.csv"))
  sample_values_csv <- file.path(paths$tables, paste0(filename_base, "_sample_values.csv"))

  ggsave(plot_pdf, plot, device = "pdf", width = 5.2, height = 4.6)
  ggsave(plot_png, plot, device = "png", width = 5.2, height = 4.6, dpi = 300)
  write_csv(summary_rows, summary_csv)
  write_csv(sample_values, sample_values_csv)

  c(plot_pdf = plot_pdf, plot_png = plot_png, summary_csv = summary_csv, sample_values_csv = sample_values_csv)
}

run_correlation_analysis <- function(dataset_name, gene, target_gene,
                                     data_dir = analysis_data_dir,
                                     save_outputs = FALSE, output_dir = "results/gene_correlations") {
  dataset <- load_dataset(dataset_name, data_dir)
  correlation_data <- make_correlation_data(dataset$expression, dataset$clinical, gene, target_gene)

  if (nrow(correlation_data) < 3) {
    return(list(valid = FALSE, reason = "Fewer than three complete samples."))
  }

  test_rows <- run_correlation_tests(correlation_data)
  summary_rows <- summarise_correlation(correlation_data, test_rows, dataset_name, gene, target_gene)
  plot <- make_correlation_plot(correlation_data, test_rows, dataset_name, gene, target_gene)
  sample_values <- correlation_sample_values(correlation_data, dataset_name, gene, target_gene)
  output_paths <- NULL

  if (save_outputs) {
    output_paths <- save_correlation_outputs(plot, summary_rows, sample_values, dataset_name, gene, target_gene, output_dir)
  }

  list(
    valid = TRUE,
    data = correlation_data,
    tests = test_rows,
    summary = summary_rows,
    sample_values = sample_values,
    plot = plot,
    output_paths = output_paths
  )
}

.app_cache <- new.env(parent = emptyenv())

app_dataset_registry <- function(project_dir = ".") {
  project_dir <- normalizePath(project_dir, mustWork = FALSE)
  cache_key <- paste("registry", project_dir, sep = "::")
  if (exists(cache_key, envir = .app_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .app_cache, inherits = FALSE))
  }

  sha_dir <- file.path(project_dir, analysis_data_dir)
  lenz_dir <- file.path(project_dir, "curated_datasets", "lenz")
  lenz_info_path <- file.path(lenz_dir, "dataset_info.rds")
  lenz_info <- if (file.exists(lenz_info_path)) readRDS(lenz_info_path) else NULL

  registry <- list()

  if (all(file.exists(unlist(analysis_datasets(sha_dir))))) {
    registry$sha_remodlb <- list(
      id = "sha_remodlb",
      label = "Sha/REMoDLB",
      reference = "Sha et al. / REMoDLB",
      summary = "928 patients; R-CHOP and R-B-CHOP arms; OS and PFS; COO, MYC/BCL2/BCL6 rearrangements, double-hit and double-expressor annotations.",
      data_type = "Microarray expression, curated gene symbols",
      loader = "sha",
      data_dir = sha_dir,
      cohorts = list(
        total = list(label = "Total", n = 928),
        rchop = list(label = "R-CHOP", n = 469),
        rbchop = list(label = "R-B-CHOP", n = 459)
      ),
      outcomes = survival_outcomes,
      stratifications = survival_stratifications,
      boxplot_comparisons = boxplot_comparisons,
      default_cohort = "total",
      default_outcome = "PFS",
      default_stratification = "basic",
      default_boxplot = "MYC_RNA"
    )
  }

  if (all(file.exists(file.path(lenz_dir, c("expression.rds", "clinical.rds"))))) {
    lenz_cohorts <- if (!is.null(lenz_info$cohorts)) lenz_info$cohorts else list(
      total = list(label = "Total", n = NA_integer_),
      chop = list(label = "CHOP", filter_column = "Treatment", filter_value = "CHOP", n = NA_integer_),
      rchop = list(label = "R-CHOP", filter_column = "Treatment", filter_value = "R-CHOP", n = NA_integer_)
    )
    registry$lenz_gse10846 <- list(
      id = "lenz_gse10846",
      label = lenz_info$display_name %||% "Lenz GSE10846",
      reference = lenz_info$reference_label %||% "Lenz et al. GSE10846",
      summary = paste0(
        lenz_info$n_samples %||% 414, " patients; CHOP and R-CHOP cohorts; OS only; ",
        "COO, stage, gender, and treatment annotations. PFS and rearrangement metadata are not available in the inspected public files."
      ),
      data_type = lenz_info$expression_type %||% "Microarray expression, gene symbols",
      loader = "single_expression_clinical",
      data_dir = lenz_dir,
      expression = file.path(lenz_dir, "expression.rds"),
      clinical = file.path(lenz_dir, "clinical.rds"),
      cohorts = lenz_cohorts,
      outcomes = list(OS = list(time = "OS_years", event = "OS_event", label = "Overall survival", xlab = "Time (years)", break_time_by = 1, axis_step = 1, axis_padding = 0.5, axis_left_padding = 0.4)),
      stratifications = list(
        basic = list(type = "expression_only", label = "Gene expression"),
        myc_expression = list(type = "computed", column = "MYC_expression_group", label = "MYC expression median"),
        treatment = list(type = "metadata", column = "Treatment", label = "Treatment", order = c("CHOP", "R-CHOP"), available_in = "total"),
        coo_class = list(type = "metadata", column = "COO_class", label = "COO class", order = c("GCB", "ABC", "UNC")),
        stage = list(type = "metadata", column = "Stage", label = "Stage", order = c("1", "2", "3", "4")),
        gender = list(type = "metadata", column = "Gender", label = "Gender", order = c("female", "male"))
      ),
      boxplot_comparisons = list(
        Treatment = list(column = "Treatment", title = "Treatment group", order = c("CHOP", "R-CHOP"), available_in = "total"),
        COO_class = list(column = "COO_class", title = "COO class", order = c("GCB", "ABC", "UNC")),
        Stage = list(column = "Stage", title = "Ann Arbor stage", order = c("1", "2", "3", "4")),
        Gender = list(column = "Gender", title = "Gender", order = c("female", "male"))
      ),
      default_cohort = "total",
      default_outcome = "OS",
      default_stratification = "basic",
      default_boxplot = "COO_class"
    )
  }

  assign(cache_key, registry, envir = .app_cache)
  registry
}

app_dataset_choices <- function(project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  stats::setNames(names(registry), vapply(registry, `[[`, character(1), "label"))
}

app_cohort_choices <- function(dataset_id, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  dataset <- registry[[dataset_id]]
  stats::setNames(names(dataset$cohorts), vapply(dataset$cohorts, function(x) paste0(x$label, if (!is.null(x$n) && !is.na(x$n)) paste0(" (n=", x$n, ")") else ""), character(1)))
}

load_app_dataset <- function(dataset_id, cohort_id, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  if (!dataset_id %in% names(registry)) {
    stop("Unknown app dataset: ", dataset_id)
  }
  dataset <- registry[[dataset_id]]

  if (identical(dataset$loader, "sha")) {
    loaded <- load_dataset(cohort_id, dataset$data_dir)
  } else {
    expression <- readRDS(dataset$expression)
    clinical <- readRDS(dataset$clinical)
    clinical$Sample_ID <- as.character(clinical$Sample_ID)
    cohort <- dataset$cohorts[[cohort_id]]
    if (is.null(cohort)) {
      stop("Unknown cohort: ", cohort_id)
    }
    if (!is.null(cohort$filter_column) && !is.null(cohort$filter_value)) {
      clinical <- clinical %>%
        filter(.data[[cohort$filter_column]] %in% cohort$filter_value)
      expression <- expression[, clinical$Sample_ID, drop = FALSE]
    }
    clinical <- clinical[match(colnames(expression), clinical$Sample_ID), , drop = FALSE]
    loaded <- list(expression = expression, clinical = clinical)
  }

  if (!identical(colnames(loaded$expression), loaded$clinical$Sample_ID)) {
    stop("Could not align clinical metadata to expression columns for ", dataset_id, " / ", cohort_id)
  }

  loaded$dataset_label <- dataset$label
  loaded$cohort_label <- dataset$cohorts[[cohort_id]]$label
  loaded
}

load_app_clinical <- function(dataset_id, cohort_id, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  if (!dataset_id %in% names(registry)) {
    stop("Unknown app dataset: ", dataset_id)
  }
  dataset <- registry[[dataset_id]]
  cohort <- dataset$cohorts[[cohort_id]]
  if (is.null(cohort)) {
    stop("Unknown cohort: ", cohort_id)
  }

  if (identical(dataset$loader, "sha")) {
    datasets <- analysis_datasets(dataset$data_dir)
    clinical <- readRDS(datasets[[cohort_id]]$clinical)
  } else {
    clinical <- readRDS(dataset$clinical)
    clinical$Sample_ID <- as.character(clinical$Sample_ID)
    if (!is.null(cohort$filter_column) && !is.null(cohort$filter_value)) {
      clinical <- clinical %>%
        filter(.data[[cohort$filter_column]] %in% cohort$filter_value)
    }
  }

  clinical
}

app_available_genes <- function(dataset_id, cohort_id = NULL, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  cohort_id <- cohort_id %||% registry[[dataset_id]]$default_cohort
  cache_key <- paste("genes", normalizePath(project_dir, mustWork = FALSE), dataset_id, sep = "::")
  if (exists(cache_key, envir = .app_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .app_cache, inherits = FALSE))
  }

  dataset <- registry[[dataset_id]]
  expression_path <- if (identical(dataset$loader, "sha")) {
    file.path(dataset$data_dir, "sha_total_expression.rds")
  } else {
    dataset$expression
  }
  genes <- rownames(readRDS(expression_path))
  assign(cache_key, genes, envir = .app_cache)
  genes
}

available_config_items <- function(items, clinical, cohort_id, genes = NULL) {
  keep <- vapply(items, function(item) {
    if (!is.null(item$available_in) && !cohort_id %in% item$available_in) {
      return(FALSE)
    }
    if (identical(item$type, "expression_only")) {
      return(TRUE)
    }
    if (identical(item$type, "computed")) {
      return(is.null(genes) || "MYC" %in% genes)
    }
    column <- item$column
    if (is.null(column) || !column %in% colnames(clinical)) {
      return(FALSE)
    }
    values <- as.character(clinical[[column]])
    values <- values[!is.na(values) & !tolower(values) %in% tolower(drop_values)]
    dplyr::n_distinct(values) >= 2
  }, logical(1))
  items[keep]
}

app_outcomes <- function(dataset_id, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  registry[[dataset_id]]$outcomes
}

app_survival_stratifications <- function(dataset_id, cohort_id, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  clinical <- load_app_clinical(dataset_id, cohort_id, project_dir)
  genes <- app_available_genes(dataset_id, cohort_id, project_dir)
  available_config_items(registry[[dataset_id]]$stratifications, clinical, cohort_id, genes)
}

app_boxplot_comparisons <- function(dataset_id, cohort_id, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  clinical <- load_app_clinical(dataset_id, cohort_id, project_dir)
  genes <- app_available_genes(dataset_id, cohort_id, project_dir)
  available_config_items(registry[[dataset_id]]$boxplot_comparisons, clinical, cohort_id, genes)
}

app_dataset_summary <- function(dataset_id, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  dataset <- registry[[dataset_id]]
  cohort_labels <- paste(vapply(dataset$cohorts, `[[`, character(1), "label"), collapse = ", ")
  outcome_labels <- paste(names(dataset$outcomes), collapse = ", ")
  list(
    title = dataset$label,
    reference = dataset$reference,
    summary = dataset$summary,
    data_type = dataset$data_type,
    cohorts = cohort_labels,
    outcomes = outcome_labels
  )
}

make_gene_data_app <- function(expression, clinical, gene) {
  if (!gene %in% rownames(expression)) {
    stop("Gene not found in expression matrix: ", gene)
  }
  gene_expression <- as.numeric(expression[gene, clinical$Sample_ID])
  data <- clinical %>%
    mutate(
      gene = gene,
      expression = gene_expression
    )

  if ("MYC" %in% rownames(expression)) {
    myc_expression <- as.numeric(expression["MYC", clinical$Sample_ID])
    myc_cutoff <- median(myc_expression, na.rm = TRUE)
    data <- data %>%
      mutate(
        MYC_expression = myc_expression,
        MYC_expression_cutoff = myc_cutoff,
        MYC_expression_group = if_else(MYC_expression >= myc_cutoff, "MYC_HIGH", "MYC_LOW")
      )
  }

  data
}

make_survival_data_app <- function(grouped, analysis_type, stratifications, min_modifier_n = 5) {
  if (!analysis_type %in% names(stratifications)) {
    stop("Unknown survival stratification: ", analysis_type)
  }

  stratification <- stratifications[[analysis_type]]

  if (stratification$type == "expression_only") {
    return(
      grouped %>%
        mutate(
          modifier = "none",
          modifier_group = "none",
          modifier_group_raw = "none",
          survival_group = expression_group
        )
    )
  }

  if (!stratification$column %in% colnames(grouped)) {
    stop("Clinical/derived column not found: ", stratification$column)
  }

  grouped %>%
    mutate(modifier_group_raw = as.character(.data[[stratification$column]])) %>%
    filter(
      !is.na(modifier_group_raw),
      !tolower(modifier_group_raw) %in% tolower(drop_values)
    ) %>%
    mutate(
      modifier_group = if ("labels" %in% names(stratification)) {
        dplyr::recode(modifier_group_raw, !!!stratification$labels, .default = toupper(modifier_group_raw))
      } else {
        toupper(modifier_group_raw)
      },
      modifier_order = if ("order" %in% names(stratification)) {
        match(modifier_group_raw, stratification$order)
      } else {
        match(modifier_group_raw, unique(modifier_group_raw))
      }
    ) %>%
    filter(!is.na(modifier_order)) %>%
    add_count(modifier_group, name = "modifier_n") %>%
    filter(modifier_n >= min_modifier_n) %>%
    arrange(modifier_order, modifier_group_raw, expression_group) %>%
    mutate(
      modifier = stratification$column,
      modifier_group = factor(modifier_group, levels = unique(modifier_group)),
      survival_group = paste(as.character(modifier_group), expression_group, sep = "_")
    ) %>%
    select(-modifier_order, -modifier_n)
}

summarise_survival_fit_app <- function(fit_result, dataset_label, cohort_label, gene, grouping, outcome_name, analysis_type, outcomes) {
  data <- fit_result$data
  outcome_def <- outcomes[[outcome_name]]

  group_counts <- data %>%
    count(survival_group, name = "n_group") %>%
    mutate(survival_group = as.character(survival_group))

  event_counts <- data %>%
    group_by(survival_group) %>%
    summarise(events = sum(.data[[outcome_def$event]] == 1), .groups = "drop") %>%
    mutate(survival_group = as.character(survival_group))

  med <- surv_median(fit_result$fit)
  median_stats <- setNames(med$median, sub("^survival_group=", "", med$strata))

  tibble(
    dataset = dataset_label,
    cohort = cohort_label,
    gene = gene,
    analysis_type = analysis_type,
    grouping = grouping,
    outcome = outcome_name,
    pvalue_logrank = fit_result$pvalue,
    n_total = nrow(data),
    cutoff_low = unique(data$cutoff_low),
    cutoff_high = unique(data$cutoff_high),
    myc_expression_cutoff = if ("MYC_expression_cutoff" %in% names(data)) unique(data$MYC_expression_cutoff) else NA_real_,
    max_time = max(data[[outcome_def$time]], na.rm = TRUE)
  ) %>%
    crossing(survival_group = as.character(levels(data$survival_group))) %>%
    left_join(group_counts, by = "survival_group") %>%
    left_join(event_counts, by = "survival_group") %>%
    mutate(median_survival = unname(median_stats[survival_group]))
}

make_survival_plot_app <- function(fit_result, dataset_label, cohort_label, gene, grouping, outcome_name, analysis_type, outcomes) {
  levels_group <- levels(fit_result$data$survival_group)
  n_groups <- length(levels_group)
  outcome <- outcomes[[outcome_name]]
  max_time <- max(fit_result$data[[outcome$time]], na.rm = TRUE)
  axis <- survival_axis_config(outcome, max_time)
  risk_table_height <- min(0.45, max(0.24, 0.08 + 0.055 * n_groups))
  risk_table_font_size <- dplyr::case_when(
    n_groups <= 2 ~ 4.8,
    n_groups <= 4 ~ 4.0,
    TRUE ~ 3.4
  )

  plot <- suppressMessages(suppressWarnings(ggsurvplot(
    fit_result$fit,
    data = fit_result$data,
    conf.int = FALSE,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.title = "",
    risk.table.height = risk_table_height,
    risk.table.fontsize = risk_table_font_size,
    xlim = c(axis$x_axis_min, axis$x_axis_max),
    break.time.by = outcome$break_time_by %||% 12,
    xlab = axis$xlab,
    ylab = paste0(outcome_name, " probability"),
    title = paste(gene, outcome_name, dataset_label, cohort_label, analysis_type, grouping, sep = " - "),
    legend.title = "Group",
    legend.labs = levels_group,
    palette = survival_group_palette(levels_group),
    ggtheme = theme_bw(base_size = 11) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)),
    risk.table.y.text = FALSE,
    tables.theme = theme_bw(base_size = 10)
  )))

  suppressMessages(suppressWarnings({
    plot$plot <- plot$plot +
      scale_x_continuous(limits = c(axis$x_axis_min, axis$x_axis_max), breaks = axis$breaks, expand = expansion(mult = c(0.005, 0.015))) +
      scale_y_continuous(limits = c(0, 1.03), expand = expansion(mult = c(0.015, 0.035))) +
      coord_cartesian(xlim = c(axis$x_axis_min, axis$x_axis_max), ylim = c(0, 1.03), clip = "off") +
      theme(plot.title = element_text(size = 11, face = "bold"))
    plot$table <- plot$table +
      scale_x_continuous(limits = c(axis$x_axis_min, axis$x_axis_max), breaks = axis$breaks, expand = expansion(mult = c(0.005, 0.015))) +
      coord_cartesian(xlim = c(axis$x_axis_min, axis$x_axis_max), clip = "off")
  }))

  list(
    plot = plot$plot,
    risk_table = plot$table,
    risk_table_export_height = max(2, 1.1 + 0.42 * n_groups)
  )
}

survival_sample_groups_app <- function(fit_result, dataset_label, cohort_label, gene, grouping, outcome_name, analysis_type) {
  data <- fit_result$data %>%
    mutate(
      dataset = dataset_label,
      cohort = cohort_label,
      gene = gene,
      analysis_type = analysis_type,
      grouping = grouping,
      outcome = outcome_name
    )
  preferred <- c(
    "dataset", "cohort", "gene", "analysis_type", "grouping", "outcome",
    "Sample_ID", "expression", "expression_group", "MYC_expression",
    "MYC_expression_group", "modifier", "modifier_group", "modifier_group_raw",
    "survival_group", "cutoff_low", "cutoff_high", "MYC_expression_cutoff"
  )
  data %>% select(any_of(preferred), everything())
}

run_app_survival_analysis <- function(dataset_id, cohort_id, gene, grouping, outcome_name, analysis_type,
                                      project_dir = ".", min_modifier_n = 5) {
  registry <- app_dataset_registry(project_dir)
  dataset_config <- registry[[dataset_id]]
  dataset <- load_app_dataset(dataset_id, cohort_id, project_dir)
  stratifications <- app_survival_stratifications(dataset_id, cohort_id, project_dir)
  outcomes <- dataset_config$outcomes

  gene_data <- make_gene_data_app(dataset$expression, dataset$clinical, gene)
  grouped <- add_expression_group(gene_data, grouping)
  analysis_data <- make_survival_data_app(grouped, analysis_type, stratifications, min_modifier_n)
  fit_result <- fit_survival(analysis_data, outcomes[[outcome_name]])

  if (is.null(fit_result)) {
    return(list(valid = FALSE, reason = "Fewer than two survival groups after filtering"))
  }

  summary_rows <- summarise_survival_fit_app(fit_result, dataset_config$label, dataset$cohort_label, gene, grouping, outcome_name, analysis_type, outcomes)
  sample_groups <- survival_sample_groups_app(fit_result, dataset_config$label, dataset$cohort_label, gene, grouping, outcome_name, analysis_type)
  plots <- make_survival_plot_app(fit_result, dataset_config$label, dataset$cohort_label, gene, grouping, outcome_name, analysis_type, outcomes)

  list(
    valid = TRUE,
    data = fit_result$data,
    fit = fit_result$fit,
    pvalue = fit_result$pvalue,
    summary = summary_rows,
    sample_groups = sample_groups,
    plot = plots$plot,
    risk_table = plots$risk_table,
    risk_table_export_height = plots$risk_table_export_height
  )
}

run_app_boxplot_analysis <- function(dataset_id, cohort_id, gene, comparison_name,
                                     project_dir = ".", min_group_n = 5) {
  registry <- app_dataset_registry(project_dir)
  dataset_config <- registry[[dataset_id]]
  comparisons <- app_boxplot_comparisons(dataset_id, cohort_id, project_dir)
  if (!comparison_name %in% names(comparisons)) {
    stop("Unknown boxplot comparison: ", comparison_name)
  }

  dataset <- load_app_dataset(dataset_id, cohort_id, project_dir)
  cohort_label <- dataset$cohort_label
  comparison <- comparisons[[comparison_name]]
  plot_data <- make_boxplot_data(dataset$expression, dataset$clinical, gene, comparison)
  validation <- validate_boxplot_data(plot_data, min_group_n)

  if (!validation$valid) {
    return(list(valid = FALSE, reason = validation$reason))
  }

  plot_data <- validation$data
  test_result <- run_association_test(plot_data)
  comparison_label <- paste(cohort_label, comparison_name, sep = "_")
  pairwise_rows <- run_pairwise_tests(plot_data, dataset_config$label, gene, comparison_label) %>%
    mutate(cohort = cohort_label, .after = dataset)
  summary_rows <- summarise_boxplot_groups(plot_data, dataset_config$label, gene, comparison_label, test_result) %>%
    mutate(cohort = cohort_label, .after = dataset)
  plot_result <- make_boxplot(plot_data, paste(dataset_config$label, cohort_label, sep = " - "), gene, comparison_name, comparison, pairwise_rows)
  sample_groups <- boxplot_sample_groups(plot_data, dataset_config$label, gene, comparison_label) %>%
    mutate(cohort = cohort_label, .after = dataset)

  list(
    valid = TRUE,
    data = plot_data,
    test = test_result,
    summary = summary_rows,
    pairwise = pairwise_rows,
    sample_groups = sample_groups,
    plot = plot_result$plot,
    plot_width = plot_result$width,
    plot_height = plot_result$height
  )
}

run_app_correlation_analysis <- function(dataset_id, cohort_id, gene, target_gene, project_dir = ".") {
  registry <- app_dataset_registry(project_dir)
  dataset_config <- registry[[dataset_id]]
  dataset <- load_app_dataset(dataset_id, cohort_id, project_dir)
  correlation_data <- make_correlation_data(dataset$expression, dataset$clinical, gene, target_gene)

  if (nrow(correlation_data) < 3) {
    return(list(valid = FALSE, reason = "Fewer than three complete samples."))
  }

  test_rows <- run_correlation_tests(correlation_data)
  dataset_label <- paste(dataset_config$label, dataset$cohort_label, sep = " - ")
  summary_rows <- summarise_correlation(correlation_data, test_rows, dataset_label, gene, target_gene)
  plot <- make_correlation_plot(correlation_data, test_rows, dataset_label, gene, target_gene)
  sample_values <- correlation_sample_values(correlation_data, dataset_label, gene, target_gene)

  list(
    valid = TRUE,
    data = correlation_data,
    tests = test_rows,
    summary = summary_rows,
    sample_values = sample_values,
    plot = plot
  )
}
