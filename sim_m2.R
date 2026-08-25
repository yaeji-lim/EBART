library(BART)
library(erf)
library(gbex)
library(grf)
library(dplyr)
library(evd)
library(FNN)
library(ggplot2)
library(tidyr)

# ---------------------------------------------------------
# Input
# ---------------------------------------------------------

### number of iterations (Set to 5 by default for a quick test run; change to 100 to replicate the full paper results)
M <- 100
### Dimension of dataset
n <- 2000
n_test <- 500
p_values <- c(50,100,200)
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
get_df <- function(X) {
  return(7 / (1 + exp(4 * X[, 1] + 1.2)) + 3)
}


get_scale <- function(X) {
  f_x <- 1 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]
  return(exp(0.1 * f_x - 5))
}


get_shape <- function(X) {
  return(1 / get_df(X))
}

generate_model2_data <- function(n, p) {
  X <- matrix(runif(n * p, min = -1, max = 1), nrow = n, ncol = p)
  df_x    <- get_df(X)
  scale_x <- get_scale(X)
  Y <- scale_x * rt(n, df = df_x)
  return(list(X = X, Y = Y, df_x = df_x, scale_x = scale_x))
}


get_true_quantile <- function(X, tau) {
  df_x    <- get_df(X)
  scale_x <- get_scale(X)
  return(scale_x * qt(tau, df = df_x))
}



# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

results_list <- list()
var_sel_list <- list()

for (p in p_values) {
  
  varimp_erf  <- matrix(0, M, p)
  varimp_gbex <- matrix(0, M, p)
  varimp_ebart<- matrix(0, M, p)
  # Top 5 detection tracking 
  det_erf <- matrix(0, M, 5)
  det_gbex <- matrix(0, M, 5)
  det_ebart <- matrix(0, M, 5)
  
  for (m in 1:M) {
 
    set.seed(1000 * p + m)
    train_data <- generate_model2_data(n, p)
    test_data  <- generate_model2_data(n_test, p)
    
    x.train <- train_data$X; y.train <- train_data$Y
    x.test  <- test_data$X;  y.test  <- test_data$Y
    
    df.train <- as.data.frame(x.train)
    df.test  <- as.data.frame(x.test)
    
    # Step 1: Intermediate Threshold Estimation u(x) via ERF
    fit_u <- erf(x.train, y.train, intermediate_quantile = tau0)
    u_train_hat <- as.vector(predict(fit_u, newdata = x.train, quantiles = tau0))
    u_test_hat  <- as.vector(predict(fit_u, newdata = x.test, quantiles = tau0))
    
    # ERF variable importance
  
    vi_erf <- as.vector(variable_importance(fit_u$quantile_forest))
    varimp_erf[m, ] <- vi_erf
    top5_erf <- order(vi_erf, decreasing = TRUE)[1:5]
    for(v in 1:5) { if(v %in% top5_erf) det_erf[m, v] <- 1 }
    
    
    erf_q <- matrix(NA, n_test, length(taus))
    for (tk in 1:length(taus)) {
      try({
        erf_q[, tk] <- as.vector(predict(fit_u, newdata = x.test, quantiles = taus[tk]))
      }, silent = TRUE)
    }
    
    z_train_full <- y.train - u_train_hat
    pos_idx <- which(z_train_full > 0)
    z_pos <- z_train_full[pos_idx]
    x_pos <- x.train[pos_idx, ]
    n_pos <- length(z_pos)
    
    
    # 2. GBEX Quantiles
    gbex_q <- matrix(NA, n_test, length(taus))
    z_train_pos <- z_train_full[pos_idx]
    x_train_pos <- x.train[pos_idx, ]
    
    gbex_sigma_hat <- rep(1, n_test)
    gbex_xi_hat <- rep(0.1, n_test)


      cor_vals <- abs(cor(x_train_pos, z_train_pos))
      cor_vals[is.na(cor_vals)] <- 0
      top_vars <- order(cor_vals, decreasing = TRUE)[1:min(20, p)]
      
      fit_gbex <- gbex(y = z_train_pos, X = data.frame(x_train_pos[, top_vars]), B = 30, silent = TRUE)
      
      if (!inherits(fit_gbex, "try-error")) {
      gbex_pred <- predict(fit_gbex, newdata = data.frame(x.test[, top_vars]))
      gbex_sigma_hat <- gbex_pred[, 1]
      gbex_xi_hat = gbex_pred[, 2]
      
      # Extract variable importance for GBEX
      vi_full <- rep(0, p)
      for(k in 1:length(fit_gbex$trees_sigma)) {
        imp <- fit_gbex$trees_sigma[[k]]$tree$variable.importance
        if(!is.null(imp)) {
          vnames <- names(imp)
          vidx <- top_vars[as.numeric(gsub("X", "", vnames))]
          vi_full[vidx] <- vi_full[vidx] + imp
        }
      }
  
      if(sum(vi_full) > 0) vi_full <- vi_full / sum(vi_full)
      vi_gbex <- vi_full
    } else {
      gbex_sigma_hat <- rep(1, n_test)
      gbex_xi_hat <- rep(0.1, n_test)
      vi_gbex <- rep(0, p)
    }
    
    varimp_gbex[m, ] <- vi_gbex
    top5_gbex <- order(vi_gbex, decreasing = TRUE)[1:5]
    for(v in 1:5) { if(v %in% top5_gbex) det_gbex[m, v] <- 1 }
    
    

    for (tk in 1:length(taus)) {
      t_val <- taus[tk]
      gbex_q[, tk] <- u_test_hat + (gbex_sigma_hat / gbex_xi_hat) * (((1 - t_val)/(1 - tau0))^(-gbex_xi_hat) - 1)
    }
    
    # 3. Exceedances & Anchors
    anchor_indices <- if (n_pos > N_a_max) sample(1:n_pos, N_a_max) else 1:n_pos
    x_anchors_sub  <- x_pos[anchor_indices, , drop = FALSE]
    n_anchors      <- length(anchor_indices)
    
    # Use the first 5 variables (signal variables) for neighborhood search (oracle distance)
    knn_train_pos <- get.knnx(data = x_pos[, 1:5], query = x_anchors_sub[, 1:5], k = min(K, n_pos))$nn.index
    knn_test_pos  <- get.knnx(data = x_pos[, 1:5], query = x.test[, 1:5], k = min(K, n_pos))$nn.index
    
    # 4. Raw MLE via KNN
    raw_mle_q <- matrix(NA, n_test, length(taus))
    for (i in 1:n_test) {
      neighbors <- knn_test_pos[i, 1:min(K, n_pos)]
      try({
        z_subset <- z_pos[neighbors]
        m_z <- mean(z_subset)
        if (m_z <= 0) m_z <- 1.0
        fit <- fpot(z_subset / m_z, threshold = 0, model = "gpd")
        s_hat <- fit$estimate[1] * m_z
        xi_hat <- fit$estimate[2]
        # Safety check: clamp xi to [-0.5, 1.5] to prevent numerical divergence
        if (xi_hat > 1.5) xi_hat <- 1.5
        if (xi_hat < -0.5) xi_hat <- -0.5
        for (tk in 1:length(taus)) {
          t_val <- taus[tk]
          raw_mle_q[i, tk] <- u_test_hat[i] + (s_hat / xi_hat) * (((1 - t_val)/(1 - tau0))^(-xi_hat) - 1)
        }
      }, silent = TRUE)
    }
    
    # 5. Extreme-BART (DART on Scale & Shape)
    anchor_sigma <- rep(NA, n_anchors)
    anchor_xi    <- rep(NA, n_anchors)
    for (j in 1:n_anchors) {
      neighbors <- knn_train_pos[j, 1:min(K, n_pos)]
      try({
        z_subset <- z_pos[neighbors]
        m_z <- mean(z_subset)
        if (m_z <= 0) m_z <- 1.0
        fit <- fpot(z_subset / m_z, threshold = 0, model = "gpd")
        anchor_sigma[j] <- fit$estimate[1] * m_z
        anchor_xi[j]    <- fit$estimate[2]
      }, silent = TRUE)
    }
    
    valid <- !is.na(anchor_sigma) & (anchor_sigma > 0) & !is.na(anchor_xi) & (anchor_xi > 0.01) & (anchor_xi < 1.5)
    x_anchors_valid  <- x_anchors_sub[valid, , drop = FALSE]
    s_anchors_valid  <- anchor_sigma[valid]
    xi_anchors_valid <- anchor_xi[valid]
    
    if (sum(valid) >= 5) {
      fit_sigma_bart <- wbart(x.train = x_anchors_valid, y.train = log(s_anchors_valid), x.test = x.test, 
                              nskip = 100, ndpost = 200, sparse = TRUE, printevery = 2000)
      ebart_sigma_hat <- exp(fit_sigma_bart$yhat.test.mean)
      
      # Extract variable importance for Extreme-BART
      var_counts <- fit_sigma_bart$varcount
      vi_ebart <- colMeans(var_counts > 0)
      varimp_ebart[m, ] <- vi_ebart
      top5_ebart <- order(vi_ebart, decreasing = TRUE)[1:5]
      for(v in 1:5) { if(v %in% top5_ebart) det_ebart[m, v] <- 1 }
      
      fit_xi_bart <- wbart(x.train = x_anchors_valid, y.train = xi_anchors_valid, x.test = x.test, 
                           nskip = 100, ndpost = 200, sparse = TRUE, printevery = 2000)
      ebart_xi_hat= fit_xi_bart$yhat.test.mean
    } else {
      ebart_sigma_hat <- rep(mean(anchor_sigma, na.rm = TRUE), n_test)
      ebart_xi_hat <- rep(mean(anchor_xi, na.rm = TRUE), n_test)
      ebart_sigma_hat[is.na(ebart_sigma_hat)] <- 1.0
      ebart_xi_hat[is.na(ebart_xi_hat)] <- 0.1
    }
    
    
    
 
    
    
    
    ebart_q <- matrix(NA, n_test, length(taus))
    for (tk in 1:length(taus)) {
      t_val <- taus[tk]
      ebart_q[, tk] <- u_test_hat + (ebart_sigma_hat / ebart_xi_hat) * (((1 - t_val) / (1 - tau0))^(-ebart_xi_hat) - 1)
    }
    
    for (tk in 1:length(taus)) {
      t_val <- taus[tk]
      true_q <- get_true_quantile(x.test, t_val)
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Rep = m,
        p = p,
        Tau = t_val,
        Raw_MLE = mean((true_q - raw_mle_q[, tk])^2, na.rm = TRUE),
        ERF = mean((true_q - erf_q[, tk])^2, na.rm = TRUE),
        GBEX = mean((true_q - gbex_q[, tk])^2, na.rm = TRUE),
        Extreme_BART = mean((true_q - ebart_q[, tk])^2, na.rm = TRUE)
      )
    }
  }
  
  # Calculate Variable Selection metrics for this batch
  
  
  for (m in 1:M) {
    # ERF
    imp_erf <- varimp_erf[m, ]
    rank_erf <- rank(-imp_erf, ties.method = "first")
    mean_rank_erf <- mean(rank_erf[1:5])
    max_erf <- max(imp_erf)
    norm_erf <- if (max_erf > 0) imp_erf / max_erf else imp_erf
    tp_erf <- sum(norm_erf[1:5] > 0.1) / 5
    top5_count_erf <- sum(det_erf[m, ]) / 5 # 추가: 상위 5개 중 참 신호 개수
    
    # GBEX
    imp_gbex <- varimp_gbex[m, ]
    rank_gbex <- rank(-imp_gbex, ties.method = "first")
    mean_rank_gbex <- mean(rank_gbex[1:5])
    max_gbex <- max(imp_gbex)
    norm_gbex <- if (max_gbex > 0) imp_gbex / max_gbex else imp_gbex
    tp_gbex <- sum(norm_gbex[1:5] > 0.1) / 5
    top5_count_gbex <- sum(det_gbex[m, ]) / 5 # 추가: 상위 5개 중 참 신호 개수
    
    # EBART
    imp_ebart <- varimp_ebart[m, ]
    rank_ebart <- rank(-imp_ebart, ties.method = "first")
    mean_rank_ebart <- mean(rank_ebart[1:5])
    max_ebart <- max(imp_ebart)
    norm_ebart <- if (max_ebart > 0) imp_ebart / max_ebart else imp_ebart
    tp_ebart <- sum(norm_ebart[1:5] > 0.1) / 5
    top5_count_ebart <- sum(det_ebart[m, ]) / 5 # 추가: 상위 5개 중 참 신호 개수
    
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
    EBART_Mean   = mean(Extreme_BART, na.rm = TRUE),
    EBART_SE     = sd(Extreme_BART, na.rm = TRUE) / sqrt(n()),
    .groups      = "drop"
  )

summary_print <- summary_df %>%
  transmute(
    p = p,
    Tau = Tau,
    Raw_MLE      = sprintf("%.2f (%.2f)", Raw_MLE_Mean, Raw_MLE_SE),
    ERF          = sprintf("%.2f (%.2f)", ERF_Mean, ERF_SE),
    GBEX         = sprintf("%.2f (%.2f)", GBEX_Mean, GBEX_SE),
    Extreme_BART = sprintf("%.2f (%.2f)", EBART_Mean, EBART_SE)
  )

print(as.data.frame(summary_print))



# 2. MSE Boxplot
mse_long <- mse_df %>%
  pivot_longer(
    cols = c(Raw_MLE, ERF, GBEX, Extreme_BART),
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
    title = "Quantile Estimation MSE across Replications (Model 2)",
    subtitle = "Note: Y-axis is in Logarithmic scale (log10)",
    x = "Method",
    y = "MSE (Log Scale)"
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
  select(p, Metric, ERF, GBEX, `Extreme-BART` = EBART)


print(as.data.frame(var_sel_table_format), row.names = FALSE)






# 4. Variable Selection Plot (Single representative run, normalized by max)


vi_erf_norm   <- if (max(vi_erf) > 0) vi_erf / max(vi_erf) else vi_erf
vi_gbex_norm  <- if (max(vi_gbex) > 0) vi_gbex / max(vi_gbex) else vi_gbex
vi_ebart_norm <- if (max(vi_ebart) > 0) vi_ebart / max(vi_ebart) else vi_ebart

df_plot <- data.frame(
  Variable = rep(1:p, 3),
  Importance = c(vi_erf_norm, vi_gbex_norm, vi_ebart_norm),
  Method = factor(rep(c("ERF", "GBEX", "Extreme-BART"), each=p), levels=c("ERF", "GBEX", "Extreme-BART")),
  Type = rep(ifelse(1:p <= 5, "Signal", "Noise"), 3)
)

df_plot <- df_plot %>%
  group_by(Method) %>%
  mutate(Rank = rank(-Importance, ties.method = "first"),
         IsTop5 = ifelse(Rank <= 5, TRUE, FALSE)) %>%
  ungroup()


p_plot <- ggplot(df_plot, aes(x=Variable, y=Importance, fill=Type)) +
  geom_bar(stat="identity") +
  geom_text(data = df_plot %>% filter(IsTop5), aes(label="*", y=Importance + 0.05), size=8, color="black", vjust=0.7) +
  facet_wrap(~Method, ncol=1, scales="free_y") +
  scale_fill_manual(values=c("Noise"="#D0D0D0", "Signal"="#303030")) +
  theme_minimal(base_size = 14) +
  labs(title=sprintf("[Model 2] Relative Variable Importance Scores (p=%d)", p),
       subtitle="Asterisks (*) indicate the top 5 variables selected by each method.",
       x="Variable Index", y="Normalized Importance Score") +
  theme(legend.position="top")


p_plot



