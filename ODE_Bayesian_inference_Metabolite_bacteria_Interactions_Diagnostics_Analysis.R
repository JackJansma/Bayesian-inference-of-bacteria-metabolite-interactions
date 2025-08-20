# Libraries needed
library(deSolve)
library(rstan)
library(tidyr)
library(ggplot2)
library(dplyr)
library(posterior)
library(stringr)

#-------------------------------------------------------------------------------------
# Diagnostics and analysis
#-------------------------------------------------------------------------------------

# Check for divergent transitions and pathological behaviour
lapply(seq_along(fits), function(i) {
  cat("==== Diagnostics for model", i, "====\n")
  rstan::check_hmc_diagnostics(fits[[i]])
})

# Check traceplots for proper mixing of multiple chains
traceplots <- lapply(seq_along(fits), function(i) {
  traceplot(fits[[i]], pars = "Cij") +
    ggtitle(paste("Traceplot for model", i)) +
    scale_x_continuous(labels = scales::label_number(accuracy = 1),
                       breaks = scales::breaks_pretty(n = 2)) +
    scale_y_continuous(labels = scales::label_scientific(digits = 1))
})

# To view the trace plots of the first simulated community:
print(traceplots[[1]])

# To check the mean, credibility interval, Rhat and n_eff values
# Make sure to check sigma_ar and cij values
# If ar1 is set to FALSE, the parameter rho should be NA
# y0_sim is used in the generated quantities block in Rstan
params_list <- lapply(seq_along(fits[[1]]), function(i) {
  rstan::summary(fits[[i]], pars = c("sigma_ar", "Cij"))$summary
})

params_list[1]

# Summarize the fit of all parameters
fit_summaries <- lapply(fits, function(f) {
  rstan::summary(f)$summary
})

# Extract n_eff and Rhat. Ideally Rhat values should be <1.01 and n_eff should be > warm up iterations
diagnostics_df <- lapply(fits, function(f) {
  summ <- rstan::summary(f)$summary
  data.frame(
    parameter = rownames(summ),
    Rhat = summ[, "Rhat"],
    neff = summ[, "n_eff"]
  )
}) %>%
  bind_rows(.id = "fit_id")

# Creates a dataframe and plot displaying a histogram of all the relevant Rhat values
selected_diagnostics <- diagnostics_df %>%
  filter(str_detect(parameter, "Cij\\[") | parameter == "sigma_ar")

Rhat <- ggplot(selected_diagnostics %>% filter(is.finite(Rhat)), aes(x = Rhat)) +
  geom_histogram(binwidth = 0.0001, fill = "steelblue", color = "black") +
  scale_y_continuous(limits = c(0, 20)) +
  labs(x = "Rhat", y = "Number of Parameters") +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

# Creates a plot displaying a histogram of all the relevant n_eff values
N_eff <- ggplot(selected_diagnostics %>% filter(is.finite(n_eff)), aes(x = n_eff)) +
  geom_histogram(binwidth = 100, fill = "darkgreen", color = "black") +
  labs(x = "n_eff", y = "Number of Parameters") +
  scale_y_continuous(limits = c(0, 12)) +
  theme_minimal()+
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

# Creates a histogram displaying the mean and credibility interval of the posterior distributions
# per estimated parameter. The stanplot is a ggplot object and can be altered accordingly.
Density_plot <- stan_plot(fits[[1]], 
                          pars = c("Cij"), 
                          show_density = FALSE,
                          inc_warmup = FALSE,
                          point_est = "mean",
                          ci_level = 0.95,
                          show_outer_line = FALSE,
                          outer_level = 0.95,
                          fill_color= "red",
                          outline_color = "black",
                          est_color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "black") +
  theme(
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = "Interaction coefficient (cij)", y = "Metabolite-bacteria pair")

# Posterior predictive check to compare the observations with generated observations using 
# draws from the obtained posterior distribution
posterior <- rstan::extract(fits[[1]])
delta <- 10e-6

y_sim <- posterior$y_sim
y_sim_corrected <- y_sim - delta # Delta was added to prevent ODE solver issues
y_obs <- y_obs_all[1,,,]

# Dimensions
R <- dim(y_sim_corrected)[1]
N_reps <- dim(y_sim_corrected)[2]
T <- dim(y_sim_corrected)[3]
L <- dim(y_sim_corrected)[4]

# -----------------------------
# Posterior predictive check for the average of the replicates
# -----------------------------
y_sim_mean_reps <- apply(y_sim_corrected, c(1, 3, 4), mean)
y_sim_mean <- apply(y_sim_mean_reps, c(2, 3), mean)
y_sim_lower <- apply(y_sim_corrected, c(3, 4), quantile, 0.025)
y_sim_upper <- apply(y_sim_corrected, c(3, 4), quantile, 0.975)

y_obs_mean <- apply(y_obs, c(2, 3), mean)
y_obs_sd <- apply(y_obs, c(2, 3), sd)

T_max <- 4
Obs <- 4
ts <- seq(10e-6, T_max, length.out = T)

df_across <- expand.grid(
  time = ts,
  variable = 1:L
)

df_across$obs <- as.vector(y_obs_mean)
df_across$obs_sd <- as.vector(y_obs_sd)
df_across$mean <- as.vector(y_sim_mean)
df_across$lower <- as.vector(y_sim_lower)
df_across$upper <- as.vector(y_sim_upper)

df_metabolites <- df_across[!df_across$variable %in% c("Number corresponding to the bacteria"), ]
df_bacteria    <- df_across[df_across$variable %in% c("Number corresponding to the bacteria"), ]

ggplot(df_metabolites, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "Simulated 95% CI"), alpha = 0.3) +
  geom_line(aes(y = mean, color = "Simulated Mean"), size = 1) +
  geom_line(aes(y = obs, color = "Observed Mean"), linetype = "dashed", linewidth = 1.2) +
  geom_errorbar(aes(ymin = obs - obs_sd, ymax = obs + obs_sd, color = "Observed Mean"), width = 0.2) +
  facet_wrap(~ variable, scales = "free_y", ncol = 5) +
  scale_color_manual(name = "Lines", values = c("Observed Mean" = "black", "Simulated Mean" = "blue")) +
  scale_fill_manual(name = "Ribbon", values = c("Simulated 95% CI" = "blue")) +
  labs(
    y = "Metabolite levels (Arbitrary units)",
    x = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        title = element_text(size = 13),
        strip.text = element_text(size = 14))

# To create a plot for the comparison of the true cij values with the corresponding posterior distributions
posterior_draws <- as_draws_df(fits[[1]])

# Extract Cij parameter names from the posterior
cij_names <- names(posterior_draws)[grepl("^Cij\\[", names(posterior_draws))]

cij_draws_long <- posterior_draws %>%
  select(all_of(cij_names)) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

# Assign true Cij values
true_df <- data.frame(
  parameter = cij_names,
  true_value = as.vector(Cij_sim_list[[1]])
)

# Merge posterior samples with true values
cij_draws_merged <- left_join(cij_draws_long, true_df, by = "parameter")

# Extract row and column indices from the parameter names
cij_draws_merged <- cij_draws_merged %>%
  mutate(
    row_idx = as.integer(str_extract(parameter, "(?<=\\[)\\d+")),
    col_idx = as.integer(str_extract(parameter, "(?<=,)\\d+"))
  ) %>%
  arrange(row_idx, col_idx) %>%
  mutate(parameter = factor(parameter, levels = unique(parameter)))

# Plot
ggplot(cij_draws_merged, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(aes(xintercept = true_value), color = "red", 
             linetype = "dashed",
             linewidth = 1) +
  facet_wrap(~ parameter, scales = "free_x", ncol = 6) +
  theme_minimal() +
  labs(x = "Posterior distribution", y = "Density")+
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::label_number(digits = 1))

# To compare the mean/median absolute error between the true cij values and corresponding posterior distributions
# Initialize storage
n_sims <- 5

rmse_mat <- matrix(NA, nrow = n_sims, ncol = 1)
mae_mat <- matrix(NA, nrow = n_sims, ncol = 1)
cor_mat <- matrix(NA, nrow = n_sims, ncol = 1)

for (rep in 1:n_sims) {
  true_vals <- as.vector(true_Cij_list[[rep]])
  est_vals  <- as.vector(posterior_Cij_list[[rep]])
  
  rmse_mat[rep] <- sqrt(mean((true_vals - est_vals)^2))
  mae_mat[rep]  <- mean(abs(true_vals - est_vals))
  cor_mat[rep]  <- cor(true_vals, est_vals)
}

mae_df <- as.data.frame(mae_mat)
mae_df$Replicate <- 1:n_sims

mae_long <- mae_df %>%
  pivot_longer(-Replicate, values_to = "MAE")

# Plot with raw values and error bars
ggplot(mae_long, aes(x = MAE)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(
    x = "Absolute Error",
    y = "Number of comparisons"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(mae_long$MAE), by = 0.04)
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    title = element_text(size = 14, face = "bold"))

# To compare the credibility intervals between multiple communities
# Initialize storage
ci_width_mean_vec <- numeric(n_sims)
ci_width_sd_vec <- numeric(n_sims)

for (i in 1:n_sims) {
  fit_i <- fits[[i]]
  post_Cij <- rstan::extract(fit_i, pars = "Cij")$Cij
  
  # Calculate 2.5% and 97.5% quantiles
  lower <- apply(post_Cij, c(2,3), quantile, probs = 0.025)
  upper <- apply(post_Cij, c(2,3), quantile, probs = 0.975)
  
  # Interval width per parameter
  widths <- upper - lower
  flat_widths <- as.vector(widths)
  
  # Store mean and sd
  ci_width_mean_vec[i] <- mean(flat_widths)
  ci_width_sd_vec[i] <- sd(flat_widths)
}

# Combine with observation counts
ci_width_df <- data.frame(
  Observations = "...",
  Mean_CI_Width = ci_width_mean_vec,
  SD_CI_Width = ci_width_sd_vec
)

# Plot CI width vs number of observations
ci_plot <- ggplot(ci_width_df, aes(x = "...", y = Mean_CI_Width)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  geom_errorbar(aes(ymin = Mean_CI_Width - SD_CI_Width,
                    ymax = Mean_CI_Width + SD_CI_Width),
                width = 0.2,
                color = "black") +
  theme_minimal() +
  labs(x = "...",
       y = "Average CI Width") +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.line = element_line(color = "black"),
    panel.border = element_blank()
  )
ci_plot