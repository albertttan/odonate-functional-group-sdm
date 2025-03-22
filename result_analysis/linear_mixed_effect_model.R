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

data <- read.csv("result_analysis/results_processed.csv")
data$distance_adjusted <- ifelse(data$distance == 0, 0.001, data$distance)
data$distance_binary <- ifelse(data$distance == 0, 0, 1)
data$distance_scaled <- scale(data$distance)
data_nonzero <- data[data$distance > 0,]
# If `year` is treated as the factor, the following warning would be:
# `fixed-effect model matrix is rank deficient so dropping 1 column / coefficient`
# data$year <- as.factor(data$year)

# Mean Suitability
hist(data$range_mean_suitability)
model_suitability <- lmer(range_mean_suitability ~ year + ssp + (1 | genus), data=data)
summary(model_suitability)
testResiduals(model_suitability)
testUniformity(model_suitability)
testDispersion(model_suitability)
vif(model_suitability, type="predictor")

# Range Size




hist(data$range_size)
model_size <- lmer(range_size ~ year + ssp + (1 | genus), data = data)
summary(model_size)
testResiduals(model_size)
testUniformity(model_size)
testDispersion(model_size)
vif(model_size, type="predictor")



data$genus_treatment <- factor(data$genus)

contrasts(data$genus_treatment) <- contr.treatment(length(levels(data$genus_treatment)))

model_size <- lmer(range_size ~ year + ssp + genus_treatment + (1 | genus), data = data)

fixed_effects <- tidy(model_size, effects = "fixed", conf.int = TRUE)

print(fixed_effects)





model_size_2 <- lm(range_size ~ year + ssp + genus_sum, data = data)
summary(model_size_2)
testResiduals(model_size_2)

genus_estimates <- coef(model_size_2)$genus
genus_estimates$last_genus <- -rowSums(genus_estimates[, -1])
print(genus_estimates)

random_effects <- tidy(model_size, effects = "ran_vals", conf.int = TRUE)

random_effects <- random_effects %>%
    select(term, estimate, conf.low, conf.high) %>%
    mutate(effect = "random")

variance_components <- tidy(model_size, effects = "ran_pars") %>%
    mutate(effect = "variance")

lm_results <- tidy(model_size_2, conf.int = TRUE) %>%
    mutate(effect = "fixed (lm)")

fixed_effects <- fixed_effects %>%
    mutate(effect = "fixed") %>%
    select(term, estimate, conf.low, conf.high, effect)

all_effects <- bind_rows(fixed_effects, random_effects, variance_components)

print(all_effects)
print(lm_results)

ggplot(all_effects, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color = effect)) +
    geom_point() +
    geom_errorbar(width = 0.2) +
    theme_minimal() +
    coord_flip() +
    labs(
        x = "Term", y = "Estimate") +
    scale_color_manual(values = c("fixed" = "blue", "random" = "green", "variance" = "purple"))

ggplot(lm_results, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_point(color = "red") +
    geom_errorbar(width = 0.2) +
    theme_minimal() +
    coord_flip() +
    labs(title = "Linear Model (lm) Estimates",
         x = "Term", y = "Estimate")




# Centroid Change Distance
hist(data$distance)
model_distance <- glmmTMB(distance ~ year + ssp + (1 | genus), data = data_nonzero, family = tweedie(link="log"), ziformula = ~ 1,
                      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(model_distance)
testResiduals(model_distance)
testUniformity(model_distance)
testDispersion(model_distance)
vif(model_distance, type="predictor")

# Integrated interpretations: huge changes in range size, small changes in habitat suitability
# Species spreading out to pockets

