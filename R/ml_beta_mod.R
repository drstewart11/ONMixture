#' Run species-specific BART model for recruitment analysis
#'
#' @param count Data frame with count data including pond_name, year, site, YC_small, BS_small
#' @param mgmt Data frame with management data including pond_name, year, and stocking/removal variables
#' @param hab Data frame with habitat variables including pond_name, year, site, pH, wtemp, doxygen, wcond, veg, wdepth
#' @param species Character string: "YCHUB" or "BSHINER"
#' @param kfold Integer. Number of folds for cross-validation-based variable selection (default = 5)
#' @return List with BART model output, selected variables, and PDP plots
#' @export
run_species_bart <- function(count, mgmt, hab, species = "YCHUB", source_path, kfold = 5) {

  if (!(species %in% c("YCHUB", "BSHINER"))) stop("Species must be 'YCHUB' or 'BSHINER'")

  count$pname <- as.numeric(as.factor(count$pond_name))
  count$yr <- as.numeric(as.factor(count$year))
  mgmt$pname <- as.numeric(as.factor(mgmt$pond_name))
  mgmt$yr <- as.numeric(as.factor(mgmt$year))
  hab$pname <- as.numeric(as.factor(hab$pond_name))
  hab$yr <- as.numeric(as.factor(hab$year))

  count$y <- if (species == "YCHUB") as.numeric(count$YC_small) else as.numeric(count$BS_small)
  newdat <- data.table(pname = count$pname, year = count$yr, site = count$site, y = count$y)

  countdata <- newdat %>%
    group_by(pname, year, site) %>%
    summarise(y = mean(y, na.rm = TRUE), .groups = "drop")

  habdata <- data.table(
    pname = hab$pname, year = hab$yr, site = hab$site,
    pH = as.numeric(hab$pH), wtemp = as.numeric(hab$wtemp),
    wcond = as.numeric(hab$wcond), oxy = as.numeric(hab$doxygen),
    veg = as.numeric(hab$veg), wdepth = as.numeric(hab$wdepth)
  )

  data_merged <- merge(countdata, habdata, by = c("pname", "year", "site"))
  data_clean <- data_merged[!is.na(data_merged$y), ]
  imputed_data <- mice(data_clean, method = "pmm", m = 5, printFlag = FALSE)
  habdata_clean <- complete(imputed_data)

  newhab <- habdata_clean %>%
    group_by(pname, year, site) %>%
    summarise(
      veg = mean(veg), depth = mean(wdepth), temp = mean(wtemp),
      cond = mean(wcond), oxy = mean(oxy), ph = mean(pH), .groups = "drop"
    )

  aligned_data <- merge(
    data_clean[, c("pname", "year", "site", "y")],
    newhab, by = c("pname", "year", "site"), all = FALSE
  )

  # Define full predictor set
  if(species=="YCHUB"){
  full_vars <- c("depth", "veg", "temp", "ph")
  y <- as.numeric(aligned_data$y)
  }else if(species=="BSHINER"){
    full_vars <- c("depth", "temp", "oxy","veg", "ph")
    y <- as.numeric(aligned_data$y)
  }


  # Cross-validation variable selection
  best_rmse <- Inf
  best_subset <- NULL
  all_subsets <- unlist(lapply(1:length(full_vars), function(k) combn(full_vars, k, simplify = FALSE)), recursive = FALSE)

  cv_results <- lapply(all_subsets, function(subset_vars) {
    folds <- createFolds(y, k = kfold, list = TRUE)
    rmse_vals <- numeric(kfold)
    for (i in seq_along(folds)) {
      test_idx <- folds[[i]]
      train_idx <- setdiff(seq_along(y), test_idx)
      X_train <- aligned_data[train_idx, subset_vars, drop = FALSE]
      y_train <- y[train_idx]
      X_test <- aligned_data[test_idx, subset_vars, drop = FALSE]
      y_test <- y[test_idx]
      fit <- suppressMessages(suppressWarnings(
        gbart(x.train = X_train, y.train = y_train, x.test = X_test)
      ))
      preds <- rowMeans(fit$yhat.test)
      rmse_vals[i] <- sqrt(mean((y_test - preds)^2))
    }
    list(vars = subset_vars, rmse = mean(rmse_vals))
  })

  best_model <- cv_results[[which.min(sapply(cv_results, function(x) x$rmse))]]
  best_subset <- best_model$vars
  best_rmse <- best_model$rmse

  message("Selected predictors: ", paste(best_subset, collapse = ", "))

  X <- as.data.frame(aligned_data[, ..best_subset])
  bart_model <- gbart(x.train = X, y.train = y)
  varimp <- bart_model$varcount.mean

  generate_pdp <- function(var_name) {
    grid_points <- 500
    X_new <- as.data.frame(t(colMeans(X)))
    var_seq <- seq(min(X[[var_name]]), max(X[[var_name]]), length.out = grid_points)
    X_new <- X_new[rep(1, grid_points), ]
    X_new[[var_name]] <- rep(var_seq, length.out = nrow(X_new))
    y_pred_manual <- predict(bart_model, newdata = X_new)
    y_pred_means <- colMeans(y_pred_manual)
    plot(var_seq, y_pred_means, type = "l", col = "blue",
         xlab = var_name, ylab = "Predicted Effect on Recruitment",
         main = paste("Partial Dependence of", var_name, "on Recruitment"))
  }

  # Save PDP figure to file by species
  output_dir <- file.path(source_path)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_filename <- file.path(output_dir, paste0("BART_PDP_", species, ".tiff"))
  tiff(filename = output_filename, width = 8, height = 8, units = "in", res = 300)
  par(mfrow = c(2, 2))
  lapply(best_subset, generate_pdp)
  par(mfrow = c(1, 1))
  dev.off()

  return(list(model = bart_model,
              importance = varimp,
              predictors = best_subset,
              aligned_data = aligned_data,
              figure_file = output_filename))
}
