# ==============================================================================
# PART 3: MODEL ESTIMATION & FORECASTING - FINAL CLEAN VERSION
# 4 Robust Models × 5 Maturities
# ==============================================================================

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║   EURO YIELD CURVE FORECASTING - PART 3: FORECASTING    ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

load(file.path(CONFIG$output$results_dir, "part2_workspace.RData"))
cat("Workspace loaded from Part 2\n\n")

# Use only robust models
MODEL_NAMES <- c("RW", "ARIMA", "VAR", "LASSO_VAR")

# Evaluation function
evaluate_forecasts <- function(actual, forecasts) {
    valid_idx <- !is.na(forecasts) & !is.na(actual) & 
                 is.finite(forecasts) & is.finite(actual)
    
    if (sum(valid_idx) < 10) {
        return(list(RMSFE = NA, MAFE = NA))
    }
    
    actual_clean <- actual[valid_idx]
    forecasts_clean <- forecasts[valid_idx]
    
    rmsfe <- sqrt(mean((actual_clean - forecasts_clean)^2))
    mafe <- mean(abs(actual_clean - forecasts_clean))
    
    return(list(RMSFE = rmsfe, MAFE = mafe))
}

# Forecasting function
forecast_single_maturity <- function(mat_idx, mat_name) {
    
    cat(sprintf("\n╔═══════════════════════════════════════╗\n"))
    cat(sprintf("║  FORECASTING: %-23s ║\n", mat_name))
    cat(sprintf("╚═══════════════════════════════════════╝\n\n"))
    
    TARGET_COL <- colnames(df_level_xts)[mat_idx]
    DF_TARGET_COL <- paste0("D_", TARGET_COL)
    
    TOTAL_OBS <- nrow(df_diff_xts)
    TEST_SIZE <- ceiling(TOTAL_OBS * 0.20)
    TRAIN_SIZE <- TOTAL_OBS - TEST_SIZE
    WINDOW_SIZE <- TRAIN_SIZE
    
    cat(sprintf("Train: %d obs | Test: %d obs\n\n", TRAIN_SIZE, TEST_SIZE))
    
    results <- data.frame(
        DATE = index(df_diff_xts)[(TRAIN_SIZE + 1):TOTAL_OBS],
        Actual = as.numeric(df_diff_xts[(TRAIN_SIZE + 1):TOTAL_OBS, DF_TARGET_COL])
    )
    
    for (model in MODEL_NAMES) {
        results[[paste0("Forecast_", model)]] <- NA
    }
    
    start_time <- Sys.time()
    
    for (t in 1:TEST_SIZE) {
        start_idx <- t
        end_idx <- WINDOW_SIZE + t - 1
        
        train_diff <- df_diff_xts[start_idx:end_idx, ]
        train_target <- as.numeric(train_diff[, DF_TARGET_COL])
        
        # 1. RANDOM WALK
        rw_forecast <- mean(train_target, na.rm = TRUE)
        results[t, "Forecast_RW"] <- rw_forecast
        
        # 2. ARIMA
        tryCatch({
            if (t == 1) {
                arima_fit <- auto.arima(train_target, d = 0, seasonal = FALSE)
                ARIMA_ORDER <<- arimaorder(arima_fit)
                cat(sprintf("ARIMA order: (%d,%d,%d)\n", 
                           ARIMA_ORDER[1], ARIMA_ORDER[2], ARIMA_ORDER[3]))
            } else {
                arima_fit <- Arima(train_target, order = ARIMA_ORDER)
            }
            arima_forecast <- forecast(arima_fit, h = 1)
            results[t, "Forecast_ARIMA"] <- as.numeric(arima_forecast$mean[1])
        }, error = function(e) {
            results[t, "Forecast_ARIMA"] <- rw_forecast
        })
        
        # 3. VAR
        tryCatch({
            if (nrow(train_diff) > CONFIG$models$var_max_lag) {
                var_select <- VARselect(train_diff, 
                                       lag.max = CONFIG$models$var_max_lag, 
                                       type = "const")
                VAR_LAGS <- var_select$selection["SC(n)"]
                var_model <- VAR(train_diff, p = VAR_LAGS, type = "const")
                var_forecast <- predict(var_model, n.ahead = 1)
                results[t, "Forecast_VAR"] <- var_forecast$fcst[[DF_TARGET_COL]][1, "fcst"]
            } else {
                results[t, "Forecast_VAR"] <- rw_forecast
            }
        }, error = function(e) {
            results[t, "Forecast_VAR"] <- rw_forecast
        })
        
        # 4. LASSO-VAR
        if (nrow(train_diff) > CONFIG$models$var_max_lag) {
            X_full <- embed(coredata(train_diff), CONFIG$models$var_max_lag + 1)
            target_pos <- which(colnames(train_diff) == DF_TARGET_COL)
            Y_lasso <- X_full[, target_pos]
            X_lasso <- X_full[, -(1:ncol(train_diff))]
            X_test <- matrix(tail(X_lasso, 1), nrow = 1)
            
            tryCatch({
                cv_lasso <- cv.glmnet(X_lasso, Y_lasso, alpha = 1, 
                                     family = "gaussian", nfolds = 5)
                lasso_fcst <- predict(cv_lasso, s = "lambda.min", newx = X_test)
                results[t, "Forecast_LASSO_VAR"] <- as.numeric(lasso_fcst[1])
            }, error = function(e) {
                results[t, "Forecast_LASSO_VAR"] <- rw_forecast
            })
        } else {
            results[t, "Forecast_LASSO_VAR"] <- rw_forecast
        }
        
        if (t %% 200 == 0) {
            cat(sprintf("  Progress: %d/%d (%.1f%%)\n", t, TEST_SIZE, 100 * t / TEST_SIZE))
        }
    }
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("✓ Completed in %.1f seconds\n", elapsed))
    
    return(results)
}

# --- RUN FOR ALL MATURITIES ---

cat("=== FORECASTING ALL MATURITIES ===\n")
cat(strrep("=", 60), "\n")

all_results <- list()

for (i in seq_along(CONFIG$data$maturity_names)) {
    mat_name <- CONFIG$data$maturity_names[i]
    all_results[[mat_name]] <- forecast_single_maturity(i, mat_name)
}

# --- PERFORMANCE EVALUATION ---

cat("\n\n=== PERFORMANCE EVALUATION ===\n\n")

performance_all <- data.frame()

for (mat_name in names(all_results)) {
    res <- all_results[[mat_name]]
    
    for (model in MODEL_NAMES) {
        metrics <- evaluate_forecasts(res$Actual, 
                                       res[[paste0("Forecast_", model)]])
        performance_all <- rbind(performance_all, data.frame(
            Maturity = mat_name,
            Model = model,
            RMSFE = metrics$RMSFE,
            MAFE = metrics$MAFE
        ))
    }
}

# Best model per maturity
best_models <- performance_all %>%
    group_by(Maturity) %>%
    slice_min(RMSFE, n = 1) %>%
    ungroup()

cat("BEST MODELS BY MATURITY:\n\n")
print(best_models %>% dplyr::select(Maturity, Model, RMSFE, MAFE), row.names = FALSE)

cat("\n\nFULL PERFORMANCE TABLE:\n\n")
perf_wide <- performance_all %>%
    dplyr::select(Maturity, Model, RMSFE) %>%
    pivot_wider(names_from = Model, values_from = RMSFE)

print(perf_wide, row.names = FALSE)

# --- VISUALIZATION ---

cat("\n\nGenerating plots...\n")

# Combined plot for all maturities
plot_list <- list()

for (mat_name in CONFIG$data$maturity_names) {
    res <- all_results[[mat_name]]
    
    res_long <- res %>%
        pivot_longer(cols = starts_with("Forecast_"),
                     names_to = "Model", values_to = "Forecast") %>%
        mutate(Model = gsub("Forecast_", "", Model)) %>%
        group_by(Model) %>%
        arrange(DATE) %>%
        mutate(Actual_Cum = cumsum(Actual),
               Forecast_Cum = cumsum(Forecast)) %>%
        ungroup()
    
    best_model_mat <- best_models$Model[best_models$Maturity == mat_name]
    
    p <- ggplot(res_long, aes(x = DATE)) +
        geom_line(aes(y = Actual_Cum), color = "black", linewidth = 0.8) +
        geom_line(aes(y = Forecast_Cum, color = Model), 
                  linewidth = 0.5, alpha = 0.7) +
        labs(title = sprintf("%s Yield", mat_name),
             subtitle = sprintf("Best: %s", best_model_mat),
             x = NULL, y = "Cumulative Change (pp)") +
        scale_color_brewer(palette = "Set1") +
        theme_minimal() +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold", size = 10))
    
    plot_list[[mat_name]] <- p
}

combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 2)

ggsave(file.path(CONFIG$output$results_dir, "forecasts_all_maturities.png"),
       combined_plot, width = 12, height = 10, dpi = 300)

cat("✓ Plots saved\n")

# --- SAVE RESULTS ---

cat("\nSaving results...\n")

write.csv(performance_all, 
          file.path(CONFIG$output$results_dir, "performance_all_maturities.csv"),
          row.names = FALSE)

write.csv(best_models,
          file.path(CONFIG$output$results_dir, "best_models_summary.csv"),
          row.names = FALSE)

for (mat_name in names(all_results)) {
    write.csv(all_results[[mat_name]],
              file.path(CONFIG$output$results_dir, 
                       sprintf("forecast_results_%s.csv", mat_name)),
              row.names = FALSE)
}

cat("✓ All results saved\n")

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║            PROJECT COMPLETED SUCCESSFULLY                ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

cat("SUMMARY:\n")
cat(sprintf("- Models: %d (RW, ARIMA, VAR, LASSO-VAR)\n", length(MODEL_NAMES)))
cat(sprintf("- Maturities: %d (3M, 1Y, 5Y, 10Y, 30Y)\n", length(CONFIG$data$maturity_names)))
cat(sprintf("- Total forecasts: %d\n", length(MODEL_NAMES) * length(CONFIG$data$maturity_names)))
cat("\nFiles saved in:", CONFIG$output$results_dir, "\n\n")