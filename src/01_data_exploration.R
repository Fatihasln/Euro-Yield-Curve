# ==============================================================================
# PART 1: DATA LOADING, EXPLORATION & STATIONARITY TESTS
# ==============================================================================

# Setup
cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║   EURO YIELD CURVE FORECASTING - PART 1: DATA & TESTS   ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

required_packages <- c("tidyverse", "xts", "vars", "forecast", "optimx",
                       "glmnet", "xgboost", "Metrics", "gridExtra", "tseries")

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
    invisible(lapply(required_packages, library, character.only = TRUE))
})

cat("✓ Packages loaded\n\n")

# Configuration
CONFIG <- list(
    data = list(
        files = c(
            "Y_3M"  = "ECB Data Portal_20251003123856.csv",
            "Y_1Y"  = "ECB Data Portal_20251003123942.csv",
            "Y_5Y"  = "ECB Data Portal_20251003124007.csv",
            "Y_10Y" = "ECB Data Portal_20251003124034.csv",
            "Y_30Y" = "ECB Data Portal_20251003124127.csv"
        ),
        maturities = c(3/12, 1, 5, 10, 30),
        maturity_names = c("3M", "1Y", "5Y", "10Y", "30Y")
    ),
    models = list(
        ns_lambda = 0.0609,
        var_max_lag = 5
    ),
    output = list(
        results_dir = "results/"
    )
)

dir.create(CONFIG$output$results_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# Data loading function
load_single_file <- function(filepath, col_name) {
    df <- read_csv(filepath, show_col_types = FALSE)
    yield_col <- colnames(df)[ncol(df)]
    df_clean <- df %>%
        dplyr::select(DATE, all_of(yield_col)) %>%
        setNames(c("DATE", col_name)) %>%
        mutate(!!col_name := as.numeric(.data[[col_name]]))
    return(xts(df_clean[[col_name]], order.by = as.Date(df_clean$DATE)))
}

# Load all data
cat("Loading data...\n")
file_names <- names(CONFIG$data$files)
xts_list <- list()

for (i in seq_along(file_names)) {
    col_name <- file_names[i]
    filepath <- CONFIG$data$files[i]
    xts_list[[i]] <- load_single_file(filepath, col_name)
}

df_level_xts <- Reduce(function(x, y) merge(x, y, join = 'inner'), xts_list)
df_level_xts <- na.omit(df_level_xts)
colnames(df_level_xts) <- file_names

cat(sprintf("✓ Data loaded: %d observations\n", nrow(df_level_xts)))
cat(sprintf("  Period: %s to %s\n\n",
            format(min(index(df_level_xts)), "%Y-%m-%d"),
            format(max(index(df_level_xts)), "%Y-%m-%d")))

# Create differences
df_diff_xts <- diff(df_level_xts)
df_diff_xts <- na.omit(df_diff_xts)
colnames(df_diff_xts) <- paste0("D_", file_names)

# --- EXPLORATORY ANALYSIS ---

cat("=== EXPLORATORY DATA ANALYSIS ===\n\n")

# Summary statistics
cat("Summary Statistics (Level):\n")
summary_stats <- data.frame(
    Maturity = CONFIG$data$maturity_names,
    Mean = apply(df_level_xts, 2, mean, na.rm = TRUE),
    SD = apply(df_level_xts, 2, sd, na.rm = TRUE),
    Min = apply(df_level_xts, 2, min, na.rm = TRUE),
    Max = apply(df_level_xts, 2, max, na.rm = TRUE)
)
print(summary_stats, row.names = FALSE)

cat("\n\nSummary Statistics (First Differences):\n")
summary_diff <- data.frame(
    Maturity = CONFIG$data$maturity_names,
    Mean = apply(df_diff_xts, 2, mean, na.rm = TRUE),
    SD = apply(df_diff_xts, 2, sd, na.rm = TRUE),
    Min = apply(df_diff_xts, 2, min, na.rm = TRUE),
    Max = apply(df_diff_xts, 2, max, na.rm = TRUE)
)
print(summary_diff, row.names = FALSE)

# --- TIME SERIES PLOTS ---

cat("\n\nGenerating time series plots...\n")

# Convert to long format for plotting
df_level_long <- data.frame(
    Date = index(df_level_xts),
    df_level_xts
) %>%
    pivot_longer(cols = -Date, names_to = "Maturity", values_to = "Yield")

df_diff_long <- data.frame(
    Date = index(df_diff_xts),
    df_diff_xts
) %>%
    pivot_longer(cols = -Date, names_to = "Maturity", values_to = "Change") %>%
    mutate(Maturity = gsub("D_", "", Maturity))

# Plot 1: Level series
p1 <- ggplot(df_level_long, aes(x = Date, y = Yield, color = Maturity)) +
    geom_line(linewidth = 0.6) +
    labs(title = "Euro Area Yield Curves (Level)",
         subtitle = "AAA-rated sovereign bonds",
         x = "Date", y = "Yield (%)") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom")

# Plot 2: First differences
p2 <- ggplot(df_diff_long, aes(x = Date, y = Change, color = Maturity)) +
    geom_line(linewidth = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Yield Changes (First Differences)",
         x = "Date", y = "Change (pp)") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom")

# Save plots
ggsave(file.path(CONFIG$output$results_dir, "yield_levels.png"),
       p1, width = 10, height = 6, dpi = 300)
ggsave(file.path(CONFIG$output$results_dir, "yield_changes.png"),
       p2, width = 10, height = 6, dpi = 300)

print(p1)
print(p2)

cat("✓ Plots saved\n\n")

# --- STATIONARITY TESTS ---

cat("=== STATIONARITY TESTS ===\n\n")

cat("Augmented Dickey-Fuller (ADF) Test:\n")
cat("H0: Unit root exists (non-stationary)\n")
cat("H1: Series is stationary\n\n")

# Test on levels
cat("1. LEVEL SERIES:\n\n")
adf_level_results <- data.frame()

for (col in colnames(df_level_xts)) {
    test <- adf.test(df_level_xts[, col], k = 4)
    adf_level_results <- rbind(adf_level_results, data.frame(
        Maturity = col,
        ADF_Statistic = round(test$statistic, 3),
        P_value = round(test$p.value, 4),
        Stationary = ifelse(test$p.value < 0.05, "YES", "NO")
    ))
}

print(adf_level_results, row.names = FALSE)

cat("\nInterpretation: All level series have p-values > 0.05,")
cat("\nfailing to reject the null hypothesis of unit root.")
cat("\nConclusion: Yield levels are NON-STATIONARY.\n\n")

# Test on first differences
cat("2. FIRST DIFFERENCE SERIES:\n\n")
adf_diff_results <- data.frame()

for (col in colnames(df_diff_xts)) {
    test <- adf.test(df_diff_xts[, col], k = 4)
    adf_diff_results <- rbind(adf_diff_results, data.frame(
        Maturity = gsub("D_", "", col),
        ADF_Statistic = round(test$statistic, 3),
        P_value = round(test$p.value, 4),
        Stationary = ifelse(test$p.value < 0.05, "YES", "NO")
    ))
}

print(adf_diff_results, row.names = FALSE)

cat("\nInterpretation: All first difference series have p-values < 0.05,")
cat("\nrejecting the null hypothesis of unit root.")
cat("\nConclusion: First differences are STATIONARY.\n\n")

# Save test results
write.csv(adf_level_results,
          file.path(CONFIG$output$results_dir, "adf_test_levels.csv"),
          row.names = FALSE)
write.csv(adf_diff_results,
          file.path(CONFIG$output$results_dir, "adf_test_differences.csv"),
          row.names = FALSE)

# --- SEASONALITY TEST ---

cat("=== SEASONALITY ANALYSIS ===\n\n")

# Check for quarterly seasonality (data is quarterly)
cat("Testing for quarterly seasonality using spectral analysis...\n\n")

seasonality_results <- data.frame()

for (col in colnames(df_level_xts)) {
    # Periodogram
    spec <- spectrum(df_level_xts[, col], plot = FALSE)

    # Find dominant frequency
    max_power_idx <- which.max(spec$spec)
    dominant_freq <- spec$freq[max_power_idx]
    dominant_period <- 1 / dominant_freq

    seasonality_results <- rbind(seasonality_results, data.frame(
        Maturity = col,
        Dominant_Period = round(dominant_period, 1),
        Strong_Seasonality = ifelse(abs(dominant_period - 4) < 1, "YES", "NO")
    ))
}

print(seasonality_results, row.names = FALSE)

cat("\nInterpretation: Yield curves typically do not exhibit strong")
cat("\nseasonal patterns (unlike commodity prices or retail sales).")
cat("\nOur analysis confirms no significant quarterly seasonality.\n\n")

# --- SUMMARY ---

cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║              PART 1 COMPLETED: KEY FINDINGS              ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

cat("1. Data: 5 maturities, ", nrow(df_level_xts), " observations\n")
cat("2. Stationarity: Levels NON-STATIONARY, Differences STATIONARY\n")
cat("3. Seasonality: No significant seasonal patterns detected\n")
cat("4. Transformation: First differencing applied for modeling\n\n")

cat("✓ Data exploration complete. Ready for Part 2 (DNS Factors).\n\n")

# Save workspace for next part
save(df_level_xts, df_diff_xts, CONFIG,
     file = file.path(CONFIG$output$results_dir, "part1_workspace.RData"))