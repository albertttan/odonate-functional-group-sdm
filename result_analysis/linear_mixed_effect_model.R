library(car)
library(mgcv)
library(broom)
library(broom.mixed)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(tweedie)
library(DHARMa)

# Load data
data <- read.csv("result_analysis/results_processed.csv")

# Transformations for models 
data$distance_adjusted <- ifelse(data$distance == 0, 0.001, data$distance)
data$distance_binary <- ifelse(data$distance == 0, 0, 1)
data$distance_scaled <- scale(data$distance)
data_nonzero <- data[data$distance > 0,]
# data$year <- as.factor(data$year) if `year` is treated as the factor, we get rank deficiency

### Mean Suitability Analysis ----
hist(data$range_mean_suitability) # check distribution

# Model with genus as fixed effect
model_suitability <- lm(range_mean_suitability ~ year + ssp + genus, data = data)
testResiduals(model_suitability)
testUniformity(model_suitability)
testDispersion(model_suitability)

vif(model_suitability) # Check for multicollinearity

# Extract model coefficients with confidence intervals
model_summary <- tidy(model_suitability, conf.int = TRUE)
print("Model Coefficients:")
print(model_summary)

# Calculate marginal means for each genus
emmeans_genus <- emmeans(model_suitability, ~ genus)
print("Marginal Means for each Genus (averaged over SSP):")
print(emmeans_genus)

# Marginal means plot
p0 <- ggplot(as.data.frame(emmeans_genus), 
       aes(x = genus, y = emmean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    theme_minimal() +
    coord_flip() +
    labs(title = "Mean Suitability by Genus",
         subtitle = "Averaged over all SSP scenarios",
         x = "Genus",
         y = "Mean Suitability") +
    theme(axis.text.y = element_text(face = "italic"))
print(p0)

# Calculate separate marginal means for each genus by SSP
emmeans_genus_ssp <- emmeans(model_suitability, ~ genus + ssp)
print("Marginal Means for each Genus by SSP:")
print(emmeans_genus_ssp)

# Plot of marginal means by SSP
p1 <- ggplot(as.data.frame(emmeans_genus_ssp), 
       aes(x = genus, y = emmean, color = ssp, group = ssp)) +
    geom_point() +
    # geom_line() + 
    theme_minimal() +
    coord_flip() +
    labs(title = "Mean Suitability by Genus and SSP",
         x = "Genus",
         y = "Mean Suitability") +
    theme(axis.text.y = element_text(face = "italic")) +
    scale_color_brewer(palette = "Set2")
print(p1) # @Albert, SSP5 is a universal improvement in mean suitability for all functional groups

# Heatmap-style visualization
p2 <- ggplot(as.data.frame(emmeans_genus_ssp), 
       aes(x = ssp, y = genus, fill = emmean)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_minimal() +
    labs(title = "Mean Suitability Heatmap",
         x = "SSP Scenario",
         y = "Genus",
         fill = "Mean\nSuitability") +
    theme(axis.text.y = element_text(face = "italic"))
print(p2)

### Range Size ----

hist(data$range_size) # check distribution

# Model with genus as fixed effect 
model_size <- lm(range_size ~ year + ssp + genus, data = data)
summary(model_size)

# Diagnostic checks 
testResiduals(model_size)
testUniformity(model_size)
testDispersion(model_size)
vif(model_size, type="predictor")

# Extract model coefficients with confidence intervals
model_summary_size <- tidy(model_size, conf.int = TRUE)
print("Model Coefficients (Range Size):")
print(model_summary_size)

# Calculate marginal means for each genus
emmeans_genus_size <- emmeans(model_size, ~ genus)
print("Marginal Means for each Genus (Range Size, averaged over SSP):")
print(emmeans_genus_size)

# Marginal means plot
p0_size <- ggplot(as.data.frame(emmeans_genus_size), 
       aes(x = genus, y = emmean)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    theme_minimal() +
    coord_flip() +
    labs(title = "Range Size by Genus",
         subtitle = "Averaged over all SSP scenarios",
         x = "Genus",
         y = "Range Size") +
    theme(axis.text.y = element_text(face = "italic"))
print(p0_size)

# Calculate separate marginal means for each genus by SSP
emmeans_genus_ssp_size <- emmeans(model_size, ~ genus + ssp)
print("Marginal Means for each Genus by SSP (Range Size):")
print(emmeans_genus_ssp_size)

# Convert to data frame for plotting
emmeans_genus_ssp_size_df <- as.data.frame(emmeans_genus_ssp_size)

# Debug print to check SSP levels
print("SSP levels in data:")
print(unique(emmeans_genus_ssp_size_df$ssp))

# Plot of marginal means by SSP (normal points)
p1_size <- ggplot(emmeans_genus_ssp_size_df, 
       aes(x = genus, y = emmean, color = ssp)) +
    geom_point(size = 2) +
    theme_minimal() +
    coord_flip() +
    labs(title = "Range Size by Genus and SSP",
         subtitle = "All scenarios",
         x = "Genus",
         y = "Range Size") +
    theme(axis.text.y = element_text(face = "italic")) +
    scale_color_brewer(palette = "Set2") +
    scale_x_discrete(limits = rev(levels(as.factor(emmeans_genus_ssp_size_df$genus))))
print(p1_size)

# Plot of marginal means by SSP (showing baseline has near-perfect overlap with SSP2)
p1_size_baseline <- ggplot(emmeans_genus_ssp_size_df, 
       aes(x = genus, y = emmean, color = ssp)) +
    geom_point(aes(size = ssp == "baseline")) +
    scale_size_manual(values = c(2, 4)) +
    theme_minimal() +
    coord_flip() +
    labs(title = "Range Size by Genus and SSP",
         subtitle = "Highlighting baseline scenario",
         x = "Genus",
         y = "Range Size") +
    theme(axis.text.y = element_text(face = "italic")) +
    scale_color_brewer(palette = "Set2") +
    scale_x_discrete(limits = rev(levels(as.factor(emmeans_genus_ssp_size_df$genus))))
print(p1_size_baseline)

# Heatmap-style visualization
p2_size <- ggplot(emmeans_genus_ssp_size_df, 
       aes(x = ssp, y = genus, fill = emmean)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_minimal() +
    labs(title = "Range Size Heatmap",
         subtitle = "Including baseline scenario",
         x = "SSP Scenario",
         y = "Genus",
         fill = "Range\nSize") +
    theme(axis.text.y = element_text(face = "italic"))
print(p2_size)

### Centroid Change Distance ----
hist(data$distance)

# simpler approach with log transformation, using non-zero distances
model_distance_log <- lm(log(distance) ~ year + ssp + genus, 
                        data = data_nonzero)
summary(model_distance_log)

# Diagnostic checks 
testResiduals(model_distance_log)
testUniformity(model_distance_log) # outlier test significant, but this is fine
testDispersion(model_distance_log)
vif(model_distance_log, type="predictor")

# Extract model coefficients with confidence intervals
model_summary_distance <- tidy(model_distance_log, conf.int = TRUE)
print("Model Coefficients (Centroid Change Distance):")
print(model_summary_distance)

# Calculate marginal means for each genus (on log scale)
emmeans_genus_distance <- emmeans(model_distance_log, ~ genus)
print("Marginal Means for each Genus (Distance, averaged over SSP):")
print(emmeans_genus_distance)

# Marginal means plot (back-transformed to original scale)
p0_distance <- ggplot(as.data.frame(emmeans_genus_distance), 
       aes(x = genus, y = exp(emmean))) +
    geom_point() +
    geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), width = 0.2) +
    theme_minimal() +
    coord_flip() +
    labs(title = "Centroid Change Distance by Genus",
         subtitle = "Averaged over all SSP scenarios",
         x = "Genus",
         y = "Distance (km)") +
    theme(axis.text.y = element_text(face = "italic"))
print(p0_distance)

# Calculate separate marginal means for each genus by SSP
emmeans_genus_ssp_distance <- emmeans(model_distance_log, ~ genus + ssp)
print("Marginal Means for each Genus by SSP (Distance):")
print(emmeans_genus_ssp_distance)

# Conversion for plotting
emmeans_genus_ssp_distance_df <- as.data.frame(emmeans_genus_ssp_distance)
emmeans_genus_ssp_distance_df$response <- exp(emmeans_genus_ssp_distance_df$emmean)
emmeans_genus_ssp_distance_df$lower.CL <- exp(emmeans_genus_ssp_distance_df$lower.CL)
emmeans_genus_ssp_distance_df$upper.CL <- exp(emmeans_genus_ssp_distance_df$upper.CL)

# Plot of marginal means by SSP (normal points)
p1_distance <- ggplot(emmeans_genus_ssp_distance_df, 
       aes(x = genus, y = response, color = ssp)) +
    geom_point(size = 2) +
    theme_minimal() +
    coord_flip() +
    labs(title = "Centroid Change Distance by Genus and SSP",
         subtitle = "All scenarios",
         x = "Genus",
         y = "Distance (km)") +
    theme(axis.text.y = element_text(face = "italic")) +
    scale_color_brewer(palette = "Set2") +
    scale_x_discrete(limits = rev(levels(as.factor(emmeans_genus_ssp_distance_df$genus))))
print(p1_distance)


# Heatmap-style visualization
p2_distance <- ggplot(emmeans_genus_ssp_distance_df, 
       aes(x = ssp, y = genus, fill = response)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_minimal() +
    labs(title = "Centroid Change Distance Heatmap",
         subtitle = "Including baseline scenario",
         x = "SSP Scenario",
         y = "Genus",
         fill = "Distance\n(km)") +
    theme(axis.text.y = element_text(face = "italic"))
print(p2_distance)

# Integrated interpretations: huge changes in range size, small changes in habitat suitability
# Species spreading out to pockets

### Permutation Importance Analysis ----
# Define variables to exclude
exclude_vars <- c("species", "year", "ssp", "range_mean_suitability", 
                  "range_size", "range_centroid_x", "range_centroid_y",
                  "eval_runtime", "eval_training_auc", "eval_testing_auc",
                  "output_threshold", "genus", "distance")

# Get environmental predictor variables
env_vars <- names(data)[!names(data) %in% exclude_vars]
print("Environmental predictor variables:")
print(env_vars)

# Create a data frame with just the environmental variables
env_data <- data[, env_vars]

# Calculate correlation matrix for environmental variables
cor_matrix <- cor(env_data, use = "pairwise.complete.obs")

# Heatmap of correlations
p_cor <- ggplot(reshape2::melt(cor_matrix), 
                aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, limit = c(-1,1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Correlation Matrix of Environmental Variables",
         subtitle = "High correlations may indicate redundant variables",
         x = "",
         y = "",
         fill = "Correlation")
print(p_cor)

# Basic stats to assess relative importance
# Create a long-format data frame for analysis
importance_long <- reshape2::melt(env_data, 
                                 id.vars = NULL,
                                 variable.name = "env_var",
                                 value.name = "importance")

# Create a summary data frame with statistical information
importance_summary <- data.frame(
    variable = levels(importance_long$env_var),
    mean = tapply(importance_long$importance, importance_long$env_var, mean, na.rm = TRUE),
    sd = tapply(importance_long$importance, importance_long$env_var, sd, na.rm = TRUE),
    se = tapply(importance_long$importance, importance_long$env_var, 
                function(x) sd(x, na.rm = TRUE)/sqrt(length(na.omit(x))))
)

# Plot with error bars
p_importance <- ggplot(importance_summary, 
                      aes(x = reorder(variable, mean), y = mean)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Permutation Importance of Environmental Variables",
         subtitle = "Error bars represent standard error",
         x = "Variable",
         y = "Mean Permutation Importance")
print(p_importance)


# Print summary statistics
print("Summary of permutation importance:")
print(importance_summary[order(-importance_summary$mean), ])


