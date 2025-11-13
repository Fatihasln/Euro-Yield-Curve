# ==============================================================================
# PART 2: DYNAMIC NELSON-SIEGEL FACTOR EXTRACTION
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║   EURO YIELD CURVE FORECASTING - PART 2: DNS FACTORS    ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# Load workspace from Part 1
load(file.path(CONFIG$output$results_dir, "part1_workspace.RData"))

cat("Workspace loaded from Part 1\n\n")

# --- NELSON-SIEGEL FUNCTIONS ---

nelson_siegel_basis <- function(tau, lambda = 0.0609) {
    # Nelson-Siegel basis functions
    # tau: maturity (in years)
    # lambda: decay parameter (Diebold & Li 2006: 0.0609)

    L <- 1.0  # Level (constant)
    S <- (1.0 - exp(-lambda * tau)) / (lambda * tau)  # Slope
    C <- S - exp(-lambda * tau)  # Curvature

    return(data.frame(L = L, S = S, C = C))
}

fit_nelson_siegel_robust <- function(yields, maturities, lambda = 0.0609) {
    # Fit Nelson-Siegel model to a cross-section of yields
    # Returns: list(success, params)

    if (length(yields) < 3 || any(is.na(yields))) {
        return(list(success = FALSE, params = c(NA, NA, NA)))
    }

    # Objective: minimize sum of squared errors
    objective <- function(params) {
        basis <- nelson_siegel_basis(maturities, lambda)
        fitted <- params[1] * basis$L + params[2] * basis$S + params[3] * basis$C
        return(sum((yields - fitted)^2))
    }

    # Initial guess
    initial_guess <- c(
        mean(yields, na.rm = TRUE),  # Beta1 (Level)
        yields[length(yields)] - yields[1],  # Beta2 (Slope)
        0  # Beta3 (Curvature)
    )

    # Optimization
    result <- tryCatch({
        opt <- optimx(
            par = initial_guess,
            fn = objective,
            method = "Nelder-Mead",
            control = list(abstol = 1e-6, maxit = 500)
        )

        if (!is.null(opt) && nrow(opt) > 0) {
            params <- as.numeric(opt[1, 1:3])

            # Validate parameters
            if (!any(is.na(params)) && all(is.finite(params))) {
                # Check fit quality
                basis <- nelson_siegel_basis(maturities, lambda)
                fitted <- params[1] * basis$L + params[2] * basis$S + params[3] * basis$C
                rmse <- sqrt(mean((yields - fitted)^2))

                # Accept if RMSE < 1% (reasonable fit)
                if (rmse < 1.0) {
                    return(list(success = TRUE, params = params, rmse = rmse))
                }
            }
        }
        return(list(success = FALSE, params = c(NA, NA, NA)))
    }, error = function(e) {
        return(list(success = FALSE, params = c(NA, NA, NA)))
    })

    return(result)
}

extract_dns_factors <- function(df_level, maturities, lambda) {
    # Extract DNS factors for all time periods

    n_rows <- nrow(df_level)
    factors_matrix <- matrix(NA, nrow = n_rows, ncol = 3)
    success_count <- 0
    fit_quality <- numeric(n_rows)

    cat("Extracting DNS factors (Nelson-Siegel decomposition)...\n")
    cat("This may take several minutes...\n\n")

    pb <- txtProgressBar(min = 0, max = n_rows, style = 3)

    for (i in 1:n_rows) {
        row_data <- as.numeric(df_level[i, ])
        result <- fit_nelson_siegel_robust(row_data, maturities, lambda)

        if (result$success) {
            factors_matrix[i, ] <- result$params
            fit_quality[i] <- result$rmse
            success_count <- success_count + 1
        }

        if (i %% 500 == 0) setTxtProgressBar(pb, i)
    }

    close(pb)

    # Create dataframe
    factors_df <- data.frame(factors_matrix)
    colnames(factors_df) <- c("Beta1_L", "Beta2_S", "Beta3_C")

    # Convert to xts
    factors_xts <- xts(factors_df, order.by = index(df_level))
    factors_clean <- na.omit(factors_xts)

    # Report statistics
    success_rate <- 100 * success_count / n_rows
    avg_rmse <- mean(fit_quality[fit_quality > 0], na.rm = TRUE)

    cat(sprintf("\n\n✓ DNS Factor Extraction Complete:\n"))
    cat(sprintf("  - Success rate: %.1f%% (%d/%d)\n", success_rate, success_count, n_rows))
    cat(sprintf("  - Average fit RMSE: %.4f%%\n", avg_rmse))
    cat(sprintf("  - Usable observations: %d\n\n", nrow(factors_clean)))

    if (success_rate < 90) {
        warning("DNS extraction success rate below 90%. Consider checking data quality.")
    }

    return(factors_clean)
}

# --- EXTRACT FACTORS ---

cat("=== DYNAMIC NELSON-SIEGEL FACTOR EXTRACTION ===\n\n")

cat("Model specification:\n")
cat("Y(τ) = β₁ · 1 + β₂ · [(1-exp(-λτ))/(λτ)] + β₃ · [(1-exp(-λτ))/(λτ) - exp(-λτ)]\n\n")
cat("where:\n")
cat("  β₁ = Level factor (long-term level)\n")
cat("  β₂ = Slope factor (short-term dynamics)\n")
cat("  β₃ = Curvature factor (medium-term hump)\n")
cat("  λ  = 0.0609 (Diebold & Li 2006)\n\n")

dns_factors_level <- extract_dns_factors(
    df_level_xts,
    CONFIG$data$maturities,
    CONFIG$models$ns_lambda
)

# Calculate first differences of factors
dns_factors_diff <- diff(dns_factors_level)
dns_factors_diff <- na.omit(dns_factors_diff)

cat(sprintf("Factor differences: %d observations\n\n", nrow(dns_factors_diff)))

# --- FACTOR ANALYSIS ---

cat("=== FACTOR CHARACTERISTICS ===\n\n")

# Summary statistics
cat("Factor Summary Statistics (Level):\n")
factor_summary <- data.frame(
    Factor = c("Beta1_L (Level)", "Beta2_S (Slope)", "Beta3_C (Curvature)"),
    Mean = apply(dns_factors_level, 2, mean, na.rm = TRUE),
    SD = apply(dns_factors_level, 2, sd, na.rm = TRUE),
    Min = apply(dns_factors_level, 2, min, na.rm = TRUE),
    Max = apply(dns_factors_level, 2, max, na.rm = TRUE)
)
print(factor_summary, row.names = FALSE)

# Factor correlations
cat("\n\nFactor Correlations:\n")
factor_cor <- cor(dns_factors_level)
print(round(factor_cor, 3))

cat("\n\nInterpretation:\n")
cat("- Level (β₁): Captures long-run expectations\n")
cat("- Slope (β₂): Captures short-term policy stance\n")
cat("- Curvature (β₃): Captures medium-term risk premia\n\n")

# --- FACTOR PLOTS ---

cat("Generating factor plots...\n")

# Convert factors to long format
dns_long <- data.frame(
    Date = index(dns_factors_level),
    dns_factors_level
) %>%
    pivot_longer(cols = -Date, names_to = "Factor", values_to = "Value") %>%
    mutate(Factor = recode(Factor,
                          Beta1_L = "Level (β₁)",
                          Beta2_S = "Slope (β₂)",
                          Beta3_C = "Curvature (β₃)"))

# Plot factors
p_factors <- ggplot(dns_long, aes(x = Date, y = Value, color = Factor)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~Factor, ncol = 1, scales = "free_y") +
    labs(title = "Nelson-Siegel Latent Factors",
         subtitle = "Extracted from Euro Area Yield Curve",
         x = "Date", y = "Factor Value") +
    theme_minimal() +
    theme(legend.position = "none")

ggsave(file.path(CONFIG$output$results_dir, "dns_factors.png"),
       p_factors, width = 10, height = 8, dpi = 300)

print(p_factors)

# Factor changes
dns_diff_long <- data.frame(
    Date = index(dns_factors_diff),
    dns_factors_diff
) %>%
    pivot_longer(cols = -Date, names_to = "Factor", values_to = "Change") %>%
    mutate(Factor = recode(Factor,
                          Beta1_L = "Level (β₁)",
                          Beta2_S = "Slope (β₂)",
                          Beta3_C = "Curvature (β₃)"))

p_factors_diff <- ggplot(dns_diff_long, aes(x = Date, y = Change, color = Factor)) +
    geom_line(linewidth = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~Factor, ncol = 1, scales = "free_y") +
    labs(title = "Factor Changes (First Differences)",
         x = "Date", y = "Change") +
    theme_minimal() +
    theme(legend.position = "none")

ggsave(file.path(CONFIG$output$results_dir, "dns_factors_changes.png"),
       p_factors_diff, width = 10, height = 8, dpi = 300)

print(p_factors_diff)

cat("\n✓ Factor plots saved\n\n")

# --- FIT QUALITY ASSESSMENT ---

cat("=== FIT QUALITY ASSESSMENT ===\n\n")

# Reconstruct yields for one sample date
sample_idx <- floor(nrow(df_level_xts) / 2)
sample_date <- index(df_level_xts)[sample_idx]
actual_yields <- as.numeric(df_level_xts[sample_idx, ])
factors <- as.numeric(dns_factors_level[sample_date, ])

reconstructed_yields <- sapply(CONFIG$data$maturities, function(tau) {
    basis <- nelson_siegel_basis(tau, CONFIG$models$ns_lambda)
    factors[1] * basis$L + factors[2] * basis$S + factors[3] * basis$C
})

cat(sprintf("Sample fit assessment (Date: %s):\n\n", sample_date))

fit_comparison <- data.frame(
    Maturity = CONFIG$data$maturity_names,
    Actual = actual_yields,
    Fitted = reconstructed_yields,
    Error = actual_yields - reconstructed_yields
)

print(fit_comparison, row.names = FALSE)

cat(sprintf("\nRMSE: %.4f%%\n", sqrt(mean(fit_comparison$Error^2))))

cat("\nInterpretation: Low RMSE indicates good fit quality.\n")
cat("Nelson-Siegel captures the overall shape of the yield curve.\n\n")

# --- SAVE WORKSPACE ---

cat("Saving workspace for Part 3...\n")

dns_factors <- list(level = dns_factors_level, diff = dns_factors_diff)

save(df_level_xts, df_diff_xts, dns_factors, CONFIG,
     file = file.path(CONFIG$output$results_dir, "part2_workspace.RData"))

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║           PART 2 COMPLETED: DNS FACTORS EXTRACTED        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

cat("Summary:\n")
cat("1. Factor extraction: ", nrow(dns_factors_level), " observations\n")
cat("2. Three factors: Level, Slope, Curvature\n")
cat("3. Fit quality: Validated and acceptable\n")
cat("4. Ready for forecasting (Part 3)\n\n")