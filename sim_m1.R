library(BART)
library(erf)
library(gbex)
library(grf)
library(dplyr)
library(ggplot2)
library(evd)
library(FNN)
library(tidyr)

# ---------------------------------------------------------
# Input
# ---------------------------------------------------------

### number of iteration
M <- 50
### Dimension of dataset
n <- 2000
n_test <- 500
p <- 100
## target extrene quantile
taus <- c(0.99, 0.995, 0.999)
## intermediate quantile
tau0 <- 0.8




# ---------------------------------------------------------
# Main 
# ---------------------------------------------------------

set.seed(20260424)


df_t <- 4

get_scale <- function(X) {
  f_x <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]
  return(exp(0.1 * f_x - 1.5))
}

get_true_quantile <- function(tau, sigma, xi=0.25) {
  q_t <- qt(tau, df=4)
  return(sigma * q_t)
}

mse_list <- list()
varimp_erf <- matrix(0, M, p)
varimp_gbex <- matrix(0, M, p)
varimp_ebart <- matrix(0, M, p)

# Top 5 detection tracking (1 if true signal is in top 5, 0 otherwise)
det_erf <- matrix(0, M, 5)
det_gbex <- matrix(0, M, 5)
det_ebart <- matrix(0, M, 5)

for (m in 1:M) {
  cat("\n===================================\n")
  cat("Monte Carlo Repetition:", m, "/", M, "\n")
  cat("===================================\n")
  
  # 1. Data Gen
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
  for(v in 1:5) { if(v %in% top5_erf) det_erf[m, v] <- 1 }
  
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
    for(i in 1:length(fit_gbex$trees_sigma)) {
      imp <- fit_gbex$trees_sigma[[i]]$tree$variable.importance
      if(!is.null(imp)) {
        vnames <- names(imp)
        vidx <- as.numeric(gsub("X", "", vnames))
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
  
  # 4. Raw MLE (Local POT)
  raw_sigma_hat <- rep(NA, n_test)
  raw_xi_hat <- rep(NA, n_test)
  k_neighbors <- 100
  knn_idx <- get.knnx(data = x_train_pos, query = x.test, k = k_neighbors)$nn.index
  for(i in 1:n_test) {
    try({
      fit <- fpot(z_train_pos[knn_idx[i, ]], threshold = 0, model = "gpd")
      raw_sigma_hat[i] <- fit$estimate[1]
      raw_xi_hat[i] <- fit$estimate[2]
    }, silent = TRUE)
  }
  raw_sigma_hat[is.na(raw_sigma_hat)] <- median(raw_sigma_hat, na.rm=TRUE)
  raw_xi_hat[is.na(raw_xi_hat)] <- median(raw_xi_hat, na.rm=TRUE)
  
  # 5. Extreme-BART
  n_anchors <- length(z_train_pos)
  sigma_anchors <- rep(NA, n_anchors)
  xi_anchors <- rep(NA, n_anchors)
  for(j in 1:n_anchors) {
    dists <- colMeans((t(x_train_pos[, 1:5]) - x_train_pos[j, 1:5])^2) 
    neighbor_idx <- order(dists)[1:100]
    try({
      fit <- fpot(z_train_pos[neighbor_idx], threshold = 0, model = "gpd")
      sigma_anchors[j] <- fit$estimate[1]
      xi_anchors[j] <- fit$estimate[2]
    }, silent = TRUE)
  }
  valid <- !is.na(sigma_anchors)
  x_anchors <- x_train_pos[valid, ]
  s_anchors <- sigma_anchors[valid]
  xi_a <- xi_anchors[valid]
  
  fit_sigma_bart <- wbart(x_anchors, log(s_anchors), x.test=x.test, nskip=500, ndpost=2000, sparse=TRUE, printevery=1000L)
  ebart_sigma_hat <- exp(fit_sigma_bart$yhat.test.mean)
  ebart_xi_hat <- rep(median(xi_a), n_test)
  
  # FIXED: properly compute column means of the inclusion matrix
  var_counts <- fit_sigma_bart$varcount
  vi_ebart <- colMeans(var_counts > 0)
  varimp_ebart[m, ] <- vi_ebart
  top5_ebart <- order(vi_ebart, decreasing = TRUE)[1:5]
  for(v in 1:5) { if(v %in% top5_ebart) det_ebart[m, v] <- 1 }
  
  # Calculate MSE
  for(k in 1:length(taus)) {
    tau <- taus[k]
    true_q <- get_true_quantile(tau, true_sigma_test, true_xi)
    
    erf_q <- erf_pred[, k]
    g_q <- u_test_hat + gbex_sigma_hat * (((1 - tau) / (1 - tau0))^(-gbex_xi_hat) - 1) / gbex_xi_hat
    raw_q <- u_test_hat + raw_sigma_hat * (((1 - tau) / (1 - tau0))^(-raw_xi_hat) - 1) / raw_xi_hat
    ebart_q <- u_test_hat + ebart_sigma_hat * (((1 - tau) / (1 - tau0))^(-ebart_xi_hat) - 1) / ebart_xi_hat
    
    mse_list[[length(mse_list) + 1]] <- data.frame(
      Rep = m,
      Tau = tau,
      Raw_MLE = mean((true_q - raw_q)^2, na.rm=TRUE),
      ERF = mean((true_q - erf_q)^2, na.rm=TRUE),
      GBEX = mean((true_q - g_q)^2, na.rm=TRUE),
      EBART = mean((true_q - ebart_q)^2, na.rm=TRUE)
    )
  }
  
  
}

cat("\nMonte Carlo Simulation Complete!\n")
# ---------------------------------------------------------
# Aggregation & Plotting
# ---------------------------------------------------------

# 1. MSE Summary Table
mse_df <- bind_rows(mse_list)

mse_summary_formatted <- mse_df %>%
  group_by(Tau) %>%
  summarise(
    Raw_MLE    = sprintf("%.2f (%.2f)", mean(Raw_MLE), sd(Raw_MLE)/sqrt(M)),
    ERF        = sprintf("%.2f (%.2f)", mean(ERF), sd(ERF)/sqrt(M)),
    GBEX       = sprintf("%.2f (%.2f)", mean(GBEX), sd(GBEX)/sqrt(M)),
    EBART      = sprintf("%.2f (%.2f)", mean(EBART), sd(EBART)/sqrt(M)),
    .groups    = "drop"
  )

print(mse_summary_formatted)


# 2. MSE Boxplot
mse_long <- mse_df %>%
  pivot_longer(cols = c("Raw_MLE", "ERF", "GBEX", "EBART"), names_to = "Method", values_to = "MSE")
mse_long$Method <- factor(mse_long$Method, levels = c("Raw_MLE", "ERF", "GBEX", "EBART"))
mse_long$Tau_Factor <- paste("Tau =", mse_long$Tau)

p_box <- ggplot(mse_long, aes(x = Method, y = MSE, fill = Method)) +
  geom_boxplot(alpha=0.7) +
  facet_wrap(~ Tau_Factor, scales = "free_y") +
  theme_minimal() +
  labs(title = paste("Monte Carlo Simulation: MSE Comparison (p=",p, ")"), y = "Mean Squared Error")

p_box


# 3. Variable Selection Plot (Average Score over M runs)

df_plot <- data.frame(
  Variable = rep(1:p, 3),
  Importance = c(vi_erf, vi_gbex, vi_ebart),
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
  geom_text(data = df_plot %>% filter(IsTop5), aes(label="*", y=Importance + 0.05), size=8, color="blue", vjust=0.7) +
  facet_wrap(~Method, ncol=1, scales="free_y") +
  scale_fill_manual(values=c("Noise"="#B0B0B0", "Signal"="#D32F2F")) +
  theme_minimal(base_size = 14) +
  labs(title=sprintf("Relative Variable Importance Scores over %d MC Runs (p=%d)", M,p),
       subtitle="Blue asterisks (*) indicate the top 5 variables selected by each method.",
       x="Variable Index", y="Normalized Importance Score") +
  theme(legend.position="top")
p_plot
