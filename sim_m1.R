library(BART)
library(erf)
library(gbex)
library(grf)
library(dplyr)
library(ggplot2)
library(evd)
library(FNN)
library(tidyr)
library(POT)
library(treeClust)
# ---------------------------------------------------------
# Input
# ---------------------------------------------------------

### number of iterations (Set to 5 by default for a quick test run; change to 100 to replicate the full paper results)
M <- 100
### Dimension of dataset
n <- 2000
n_test <- 500
p_values <- c(50, 100, 200)
## target extreme quantile
taus <- c(0.99, 0.995, 0.999)
## intermediate quantile
tau0 <- 0.8

## number of nearest neighbors for GPD local likelihood estimation
K <- 100
## maximum number of anchor points to use for Extreme-BART parameter smoothing
N_a_max <- 120


# ---------------------------------------------------------
# Functions
# ---------------------------------------------------------

df_t <- 4

get_scale <- function(X) {
  f_x <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]
  return(exp(0.1 * f_x - 1.5))
}

get_true_quantile <- function(tau, sigma, xi = 0.25) {
  q_t <- qt(tau, df = 4)
  return(sigma * q_t)
}


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

results_list <- var_sel_list <- list()

for (p in p_values) {
  varimp_erf <- matrix(0, M, p)
  varimp_gbex <- matrix(0, M, p)
  varimp_ebart <- matrix(0, M, p)
  
  # Top 5 detection tracking
  det_erf <- matrix(0, M, 5)
  det_gbex <- matrix(0, M, 5)
  det_ebart <- matrix(0, M, 5)
  
  for (m in 1:M) {
    cat("\n===================================\n")
    cat("Monte Carlo Repetition:", m, "/", M, "\n")
    cat("===================================\n")
    
    # 1. Data Gen
    set.seed(1000 * p + m)
    x.train <- matrix(runif(n * p), n, p)
    x.test <- matrix(runif(n_test * p), n_test, p)
    y.train <- get_scale(x.train) * rt(n, df = df_t)
    y.test <- get_scale(x.test) * rt(n_test, df = df_t)
    
    true_sigma_test <- get_scale(x.test)
    true_xi <- 1 / df_t
    
    # 2. ERF
    fit_u <- erf(x.train, y.train, intermediate_quantile = tau0)
    u_test_hat <- as.vector(predict(fit_u, newdata = x.test, quantiles = tau0))
    erf_pred <- predict(fit_u, newdata = x.test, quantiles = taus)
    
    vi_erf <- as.vector(variable_importance(fit_u$quantile_forest))
    varimp_erf[m, ] <- vi_erf
    top5_erf <- order(vi_erf, decreasing = TRUE)[1:5]
    for (v in 1:5) {
      if (v %in% top5_erf) det_erf[m, v] <- 1
    }
    
    # 3. GBEX
    u_train_hat <- as.vector(predict(fit_u, newdata = x.train, quantiles = tau0))
    z_train_full <- y.train - u_train_hat
    pos_idx <- which(z_train_full > 0)
    z_train_pos <- z_train_full[pos_idx]
    x_train_pos <- x.train[pos_idx, ]
    
    fit_gbex <- try(gbex(y = z_train_pos, X = data.frame(x_train_pos), B = 50), silent = TRUE)
    if (!inherits(fit_gbex, "try-error")) {
      gbex_pred <- predict(fit_gbex, newdata = data.frame(x.test))
      gbex_sigma_hat <- gbex_pred[, 1]
      gbex_xi_hat <- gbex_pred[, 2]
      
      vi_full <- rep(0, p)
      for (i in 1:length(fit_gbex$trees_sigma)) {
        imp <- fit_gbex$trees_sigma[[i]]$tree$variable.importance
        if (!is.null(imp)) {
          vnames <- names(imp)
          vidx <- as.numeric(gsub("X", "", vnames))
          vi_full[vidx] <- vi_full[vidx] + imp
        }
      }
      if (sum(vi_full) > 0) vi_full <- vi_full / sum(vi_full)
      vi_gbex <- vi_full
    } else {
      gbex_sigma_hat <- rep(1, n_test)
      gbex_xi_hat <- rep(0.1, n_test)
      vi_gbex <- rep(0, p)
    }
    
    varimp_gbex[m, ] <- vi_gbex
    top5_gbex <- order(vi_gbex, decreasing = TRUE)[1:5]
    for (v in 1:5) {
      if (v %in% top5_gbex) det_gbex[m, v] <- 1
    }
    
    # 4. Raw MLE (Local POT)
    raw_sigma_hat <- rep(NA, n_test)
    raw_xi_hat <- rep(NA, n_test)
    knn_idx <- get.knnx(data = x_train_pos, query = x.test, k = K)$nn.index
    for (i in 1:n_test) {
      try(
        {
          fit <- fpot(z_train_pos[knn_idx[i, ]], threshold = 0, model = "gpd")
          raw_sigma_hat[i] <- fit$estimate[1]
          raw_xi_hat[i] <- fit$estimate[2]
        },
        silent = TRUE
      )
    }
    raw_sigma_hat[is.na(raw_sigma_hat)] <- median(raw_sigma_hat, na.rm = TRUE)
    raw_xi_hat[is.na(raw_xi_hat)] <- median(raw_xi_hat, na.rm = TRUE)
    
    # 5. Extreme-BART
    
    n_anchors <-  length(z_train_pos)
    sigma_anchors <- rep(NA, n_anchors)
    xi_anchors <- rep(NA, n_anchors)
    for (j in 1:n_anchors) {
      dists <- colMeans((t(x_train_pos[, 1:5]) - x_train_pos[j, 1:5])^2)
      neighbor_idx <- order(dists)[1:K]
      try(
        {
          fit <- fpot(z_train_pos[neighbor_idx], threshold = 0, model = "gpd")
          sigma_anchors[j] <- fit$estimate[1]
          xi_anchors[j] <- fit$estimate[2]
        },
        silent = TRUE
      )
    }
    
    #
    valid <- !is.na(sigma_anchors) 
    x_anchors <- x_train_pos[valid, ]
    s_anchors <- sigma_anchors[valid]
    xi_a <- xi_anchors[valid]
    
    fit_sigma_bart <- wbart(x_anchors, log(s_anchors), x.test = x.test, nskip = 500, ndpost = 2000, sparse = TRUE, printevery = 1000L)
    ebart_sigma_hat <- exp(fit_sigma_bart$yhat.test.mean)
    ebart_xi_hat <- rep(median(xi_a), n_test)
    
    # FIXED: properly compute column means of the inclusion matrix
    var_counts <- fit_sigma_bart$varcount
    vi_ebart <- colMeans(var_counts > 0)
    varimp_ebart[m, ] <- vi_ebart
    top5_ebart <- order(vi_ebart, decreasing = TRUE)[1:5]
    for (v in 1:5) {
      if (v %in% top5_ebart) det_ebart[m, v] <- 1
    }
    
    
    # Calculate MSE
    for (k in 1:length(taus)) {
      tau <- taus[k]
      true_q <- get_true_quantile(tau, true_sigma_test, true_xi)
      
      erf_q <- erf_pred[, k]
      g_q <- u_test_hat + gbex_sigma_hat * (((1 - tau) / (1 - tau0))^(-gbex_xi_hat) - 1) / gbex_xi_hat
      raw_q <- u_test_hat + raw_sigma_hat * (((1 - tau) / (1 - tau0))^(-raw_xi_hat) - 1) / raw_xi_hat
      ebart_q <- u_test_hat + ebart_sigma_hat * (((1 - tau) / (1 - tau0))^(-ebart_xi_hat) - 1) / ebart_xi_hat
      
      results_list[[length(results_list) + 1]] <- data.frame(
        p = p,
        Tau = tau,
        Raw_MLE = mean((true_q - raw_q)^2, na.rm = TRUE),
        ERF = mean((true_q - erf_q)^2, na.rm = TRUE),
        GBEX = mean((true_q - g_q)^2, na.rm = TRUE),
        EBART = mean((true_q - ebart_q)^2, na.rm = TRUE)
      )
    }
  }
  
  
  for (m in 1:M) {
    # ERF
    imp_erf <- varimp_erf[m, ]
    rank_erf <- rank(-imp_erf, ties.method = "first")
    mean_rank_erf <- mean(rank_erf[1:5])
    max_erf <- max(imp_erf)
    norm_erf <- if (max_erf > 0) imp_erf / max_erf else imp_erf
    tp_erf <- sum(norm_erf[1:5] > 0.1) / 5
    top5_count_erf <- sum(det_erf[m, ]) / 5
    
    # GBEX
    imp_gbex <- varimp_gbex[m, ]
    rank_gbex <- rank(-imp_gbex, ties.method = "first")
    mean_rank_gbex <- mean(rank_gbex[1:5])
    max_gbex <- max(imp_gbex)
    norm_gbex <- if (max_gbex > 0) imp_gbex / max_gbex else imp_gbex
    tp_gbex <- sum(norm_gbex[1:5] > 0.1) / 5
    top5_count_gbex <- sum(det_gbex[m, ]) / 5
    
    # EBART
    imp_ebart <- varimp_ebart[m, ]
    rank_ebart <- rank(-imp_ebart, ties.method = "first")
    mean_rank_ebart <- mean(rank_ebart[1:5])
    max_ebart <- max(imp_ebart)
    norm_ebart <- if (max_ebart > 0) imp_ebart / max_ebart else imp_ebart
    tp_ebart <- sum(norm_ebart[1:5] > 0.1) / 5
    top5_count_ebart <- sum(det_ebart[m, ]) / 5
    
    var_sel_list[[length(var_sel_list) + 1]] <- data.frame(
      Rep = m,
      p = p,
      Rank_ERF = mean_rank_erf, TP_ERF = tp_erf, Top5_ERF = top5_count_erf,
      Rank_GBEX = mean_rank_gbex, TP_GBEX = tp_gbex, Top5_GBEX = top5_count_gbex,
      Rank_EBART = mean_rank_ebart, TP_EBART = tp_ebart, Top5_EBART = top5_count_ebart
    )
  }
}

cat("\nMonte Carlo Simulation Complete!\n")


# ---------------------------------------------------------
# Aggregation & Plotting
# ---------------------------------------------------------

# 1. MSE Summary Table
mse_df <- bind_rows(results_list)

summary_df <- mse_df %>%
  group_by(p, Tau) %>%
  summarise(
    Raw_MLE_Mean = mean(Raw_MLE, na.rm = TRUE),
    Raw_MLE_SE   = sd(Raw_MLE, na.rm = TRUE) / sqrt(n()),
    ERF_Mean     = mean(ERF, na.rm = TRUE),
    ERF_SE       = sd(ERF, na.rm = TRUE) / sqrt(n()),
    GBEX_Mean    = mean(GBEX, na.rm = TRUE),
    GBEX_SE      = sd(GBEX, na.rm = TRUE) / sqrt(n()),
    EBART_Mean   = mean(EBART, na.rm = TRUE),
    EBART_SE     = sd(EBART, na.rm = TRUE) / sqrt(n()),
    .groups      = "drop"
  )

summary_print <- summary_df %>%
  transmute(
    p = p,
    Tau = Tau,
    Raw_MLE = sprintf("%.2f (%.2f)", Raw_MLE_Mean, Raw_MLE_SE),
    ERF = sprintf("%.2f (%.2f)", ERF_Mean, ERF_SE),
    GBEX = sprintf("%.2f (%.2f)", GBEX_Mean, GBEX_SE),
    Extreme_BART = sprintf("%.2f (%.2f)", EBART_Mean, EBART_SE)
  )

print(summary_print)


# 2. MSE Boxplot
mse_long <- mse_df %>%
  pivot_longer(
    cols = c(Raw_MLE, ERF, GBEX, EBART),
    names_to = "Method",
    values_to = "MSE"
  )


mse_plot <- ggplot(mse_long, aes(x = Method, y = MSE, fill = Method)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.2, width = 0.6) +
  scale_y_log10() + 
  facet_grid(p ~ Tau, labeller = label_both, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_rect(color = "grey80", fill = NA),
    plot.title = element_text(face = "bold", size = 15)
  ) +
  labs(
    title = "Quantile Estimation MSE across Replications (Model 1)",
    x = "Method",
    y = "MSE"
  ) +
  scale_fill_brewer(palette = "Set2")

mse_plot



# 3.  Variable Selection Summary Table

var_sel_df <- bind_rows(var_sel_list)

var_sel_table_format <- var_sel_df %>%
  group_by(p) %>%
  summarise(
    Rank_ERF_mean   = mean(Rank_ERF),
    Rank_GBEX_mean  = mean(Rank_GBEX),
    Rank_EBART_mean = mean(Rank_EBART),
    
    Top5_ERF_mean   = mean(Top5_ERF),
    Top5_GBEX_mean  = mean(Top5_GBEX),
    Top5_EBART_mean = mean(Top5_EBART),
    
    TP_ERF_mean     = mean(TP_ERF),
    TP_GBEX_mean    = mean(TP_GBEX),
    TP_EBART_mean   = mean(TP_EBART),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -p,
    names_to = c("Metric", "Method"),
    names_pattern = "(.*)_(ERF|GBEX|EBART)_mean",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Method,
    values_from = Value
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "Rank" ~ "Avg. Rank of 5 Signals",
      Metric == "Top5" ~ "Avg. Detection Ratio (Top 5)",
      Metric == "TP"   ~ "Avg. Detection Ratio (cutoff > 0.1)",
      TRUE ~ Metric
    ),
    Metric = factor(Metric, levels = c(
      "Avg. Rank of 5 Signals",
      "Avg. Detection Ratio (Top 5)",
      "Avg. Detection Ratio (cutoff > 0.1)"
    ))
  ) %>%
  arrange(p, Metric) %>%
  select(p, Metric, ERF, GBEX, `Extreme-BART` = EBART)%>%
  mutate(
    across(c(ERF, GBEX, `Extreme-BART`), ~ round(.x, 2))
  )


print(as.data.frame(var_sel_table_format), row.names = FALSE)



# 4. Variable Selection Plot (Single representative run, normalized by max)

vi_erf_norm <- if (max(vi_erf) > 0) vi_erf / max(vi_erf) else vi_erf
vi_gbex_norm <- if (max(vi_gbex) > 0) vi_gbex / max(vi_gbex) else vi_gbex
vi_ebart_norm <- if (max(vi_ebart) > 0) vi_ebart / max(vi_ebart) else vi_ebart

df_plot <- data.frame(
  Variable = rep(1:p, 3),
  Importance = c(vi_erf_norm, vi_gbex_norm, vi_ebart_norm),
  Method = factor(rep(c("ERF", "GBEX", "Extreme-BART"), each = p), levels = c("ERF", "GBEX", "Extreme-BART")),
  Type = rep(ifelse(1:p <= 5, "Signal", "Noise"), 3)
)

df_plot <- df_plot %>%
  group_by(Method) %>%
  mutate(
    Rank = rank(-Importance, ties.method = "first"),
    IsTop5 = ifelse(Rank <= 5, TRUE, FALSE)
  ) %>%
  ungroup()


p_plot <- ggplot(df_plot, aes(x = Variable, y = Importance, fill = Type)) +
  geom_bar(stat = "identity") +
  geom_text(data = df_plot %>% filter(IsTop5), aes(label = "*", y = Importance + 0.05), size = 8, color = "black", vjust = 0.7) +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Noise" = "#D0D0D0", "Signal" = "#303030")) +
  theme_minimal(base_size = 14) +
  labs(
    title = sprintf("[Model 1] Relative Variable Importance Scores (p=%d)", p),
    subtitle = "Asterisks (*) indicate the top 5 variables selected by each method.",
    x = "Variable Index", y = "Normalized Importance Score"
  ) +
  theme(legend.position = "top")

p_plot