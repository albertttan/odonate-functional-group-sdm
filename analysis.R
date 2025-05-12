library(car)
library(boot)
library(broom)
library(rlang)
library(DHARMa)
library(emmeans)
library(ggplot2)
library(stringr)
library(reshape2)



### Load data

set.seed(0)

baseline_data <- read.csv("output/model/baseline.csv")
future_data <- read.csv("output/model/future.csv")

baseline_data$genus <- sub("^(\\w+).*", "\\1", baseline_data$species)
future_data$genus <- sub("^(\\w+).*", "\\1", future_data$species)



### Response variables analysis

response_analysis <- function(data, response_var, unit) {
    response_sym <- ensym(response_var)
    response_str <- as_string(response_sym)
    
    
    ## Fit linear model with logged response variable
    
    model_formula <- as.formula(paste0("log(", response_str, "_change) ~ year + ssp * genus"))
    model <- lm(model_formula, data = data)
    model_tidy <- as.data.frame(tidy(model, conf.int = TRUE))
    
    
    ## Assess model assumptions
    
    hist(data[[paste0(as_string(response_sym), "_change")]])  # Visually assess distribution
    vif(model)  # Acceptable
    testResiduals(model)  # Not Acceptable
    
    
    ## Bootstrapping
    
    boot_fn <- function(data, indices) {
        d <- data[indices, ]
        mod <- lm(model_formula, data = d)
        return(coef(mod))
    }
    boot_model <- boot(data = data, statistic = boot_fn, R = 1000)
    
    
    ## Bootstrapping: Compare confidence intervals
    
    orig_ci <- setNames(
        model_tidy[c("term", "conf.low", "conf.high")], 
        c("variable", "orig_lower", "orig_upper")
    )
    boot_ci_raw <- sapply(seq_along(coef(model)), function(i) {
        ci <- boot.ci(boot_model, type = "perc", index = i)$percent
        if (is.null(ci)) c(NA, NA) else ci[4:5]
    })
    boot_ci <- data.frame(boot_lower = boot_ci_raw[1, ], boot_upper = boot_ci_raw[2, ])
    boot_ci_comparison <- cbind(orig_ci, boot_ci)
    boot_ci_comparison  # Acceptable
    
    
    ## Bootstrapping: Compare medians
    
    boot_medians <- apply(boot_model$t, 2, median)
    boot_summary <- data.frame(
        boot_median = boot_medians,
        difference = coef(model) - boot_medians
    )
    boot_summary  # Acceptable
    
    
    ## Convert log to percent in marginal means results
    
    emmeans_exp <- function(result) {
        result[c("emmean", "SE", "lower.CL", "upper.CL")] <- lapply(
            result[c("emmean", "SE", "lower.CL", "upper.CL")], ifelse(unit == "%", function(x) exp(x) * 100, exp)
        )
        return(result)
    }
    
    
    ## Calculate and visualize marginal means of genus
    
    emmeans_genus <- emmeans_exp(as.data.frame(emmeans(model, ~ genus)))
    emmeans_genus
    
    plot_genus <- ggplot(emmeans_genus, 
                         aes(x = genus, y = emmean)) +
        geom_point() +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
        theme_minimal() +
        coord_flip() +
        labs(x = "Genus", y = paste0(str_to_title(response_str), " Change (", unit, ")")) +
        theme(axis.text.y = element_text(face = "italic"))
    ggsave(paste0("output/analysis/images/", response_str, "_genus.png"), plot_genus, dpi = 500)
    
    
    ## Calculate and visualize marginal means of genus by SSP
    
    emmeans_genus_ssp <- emmeans_exp(as.data.frame(emmeans(model, ~ genus + ssp)))
    emmeans_genus_ssp
    
    plot_ssp <- ggplot(emmeans_genus_ssp, 
                       aes(x = genus, y = emmean, color = str_to_upper(ssp), group = ssp)) +
        geom_point() +
        theme_minimal() +
        coord_flip() +
        labs(x = "Genus", y = paste0(str_to_title(response_str), " Change (", unit, ")"), color = "SSP") +
        theme(axis.text.y = element_text(face = "italic")) +
        scale_color_brewer(palette = "Set2")
    ggsave(paste0("output/analysis/images/", response_str, "_ssp.png"), plot_ssp, dpi = 500)
    
    
    ## Heatmap visualization
    
    plot_heatmap <- ggplot(emmeans_genus_ssp, 
                           aes(x = str_to_upper(ssp), y = genus, fill = emmean)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "steelblue") +
        theme_minimal() +
        labs(x = "SSP Scenario", y = "Genus", fill = paste0(str_to_title(response_str), "\nChange (", unit, ")")) +
        theme(axis.text.y = element_text(face = "italic"))
    ggsave(paste0("output/analysis/images/", response_str, "_heatmap.png"), plot_heatmap, dpi = 500)
    
    
    ## Return everything useful
    
    return(list(
        model = model,
        boot_model = boot_model,
        boot_ci_comparison = boot_ci_comparison,
        boot_summary = boot_summary,
        emmeans_genus = emmeans_genus,
        emmeans_genus_ssp = emmeans_genus_ssp,
        plots = list(
            genus = plot_genus,
            ssp = plot_ssp,
            heatmap = plot_heatmap
        )
    ))
}

result_suitability <- response_analysis(future_data, suitability, "%")
result_range <- response_analysis(future_data, range, "%")
result_distance <- response_analysis(future_data, distance, "km")
save(result_suitability, result_range, result_distance, file="output/analysis/results.RData")



### Permutation importance analysis


## Select environmental predictor variables

exclude_vars <- c("species", "genus", "suitability", "range", "centroid_x", "centroid_y",
                  "runtime", "threshold", "training_auc", "testing_auc")
env_data <- baseline_data[, names(baseline_data)[!names(baseline_data) %in% exclude_vars]]
names(env_data)[names(env_data) == "slope"] <- "Slope"


## Calculate and visualize correlation matrix for environmental variables

cor_matrix <- cor(env_data, use = "pairwise.complete.obs")
cor_matrix

plot_cor <- ggplot(reshape2::melt(cor_matrix), 
                aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, limit = c(-1,1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "", fill = "Correlation")
ggsave(paste0("output/analysis/images/variable_correlation.png"), plot_cor, dpi = 500)


## Summarize and visualize permutation importance

importance_long <- melt(env_data, id.vars = NULL, variable.name = "env_var", value.name = "importance")
# If N/A values are removed instead of considered as 0:
# importance_summary <- data.frame(
#     variable = levels(importance_long$env_var),
#     mean = tapply(importance_long$importance, importance_long$env_var, mean, na.rm = TRUE),
#     sd = tapply(importance_long$importance, importance_long$env_var, sd, na.rm = TRUE),
#     se = tapply(importance_long$importance, importance_long$env_var, 
#                 function(x) sd(x, na.rm = TRUE)/sqrt(length(na.omit(x))))
# )
importance_summary <- data.frame(
    variable = levels(importance_long$env_var),
    mean = tapply(importance_long$importance, importance_long$env_var, function(x) mean(replace(x, is.na(x), 0))),
    sd = tapply(importance_long$importance, importance_long$env_var, function(x) sd(replace(x, is.na(x), 0))),
    se = tapply(importance_long$importance, importance_long$env_var, 
                function(x) sd(replace(x, is.na(x), 0)/sqrt(length(x))))
)
importance_summary[order(-importance_summary$mean), ]
save(importance_summary, file="output/analysis/variable_pi.RData")

plot_importance <- ggplot(importance_summary, 
                      aes(x = reorder(variable, mean), y = mean)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Variable", y = "Permutation Importance")
ggsave(paste0("output/analysis/images/variable_pi.png"), plot_importance, dpi = 500)
