library(tigris)
library(sf)
library(spdep)
library(dplyr)
library(brms)
library(ggplot2)
library(scales)

#Intercept from on Schein et al. 
#Using beta binomial distr to account for variance from geographically 
#localied concentration mentioned in (Weisburd et al., 2018)

options(tigris_use_cache = TRUE)

crime_ptsd <- chicago_crime_ptsd %>%
  mutate(
    zip = as.character(zip),
    ptsd_cases = round(ptsd_rate * population),
    log_violent = log(violent_rate + 1),
    hardship_s = as.numeric(scale(hardship_index)),
    uninsured_s = as.numeric(scale(uninsured_rate)),
    income_s = as.numeric(scale(med_income)),
    trials = population
  )

zcta <- zctas(cb = TRUE, year = 2020)

zcta_chicago <- zcta %>%
  filter(ZCTA5CE20 %in% crime_ptsd$zip) %>%
  rename(zip = ZCTA5CE20)

spatial_data <- zcta_chicago %>%
  left_join(crime_ptsd, by = "zip") %>%
  filter(!is.na(ptsd_cases)) 

nb <- poly2nb(spatial_data)
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(W) <- spatial_data$zip
colnames(W) <- spatial_data$zip

spatial_data$zip <- factor(spatial_data$zip, levels = rownames(W))

priors <- c(
  prior(normal(-3,0.5), class = Intercept),
  prior(normal(0,1), class = b),
  prior(exponential(1), class = phi)
)

model_crime <- brm(
  ptsd_cases | trials(trials) ~ log_violent,
  data = spatial_data,
  family = beta_binomial(),
  prior = priors,
  chains = 4, iter = 2000, cores = 4,
  control = list(adapt_delta = 0.99)
)

model_controls <- brm(
  ptsd_cases | trials(trials) ~ 
    log_violent + hardship_s + uninsured_s + income_s,
  data = spatial_data,
  family = beta_binomial(),
  prior = priors,
  chains = 4, iter = 2000, cores = 4,
  control = list(adapt_delta = 0.99)
)

model_interaction <- brm(
  ptsd_cases | trials(trials) ~ 
    log_violent * hardship_s + uninsured_s,
  data = spatial_data,
  family = beta_binomial(),
  prior = priors,
  chains = 4, iter = 2000, cores = 4,
  control = list(adapt_delta = 0.99)
)

model_spline <- brm(
  ptsd_cases | trials(trials) ~ 
    s(log_violent) + hardship_s + uninsured_s,
  data = spatial_data,
  family = beta_binomial(),
  prior = priors,
  chains = 4, iter = 2000, cores = 4,
  control = list(adapt_delta = 0.99)
)

model_random_zip <- brm(
  ptsd_cases | trials(trials) ~ 
    log_violent + hardship_s + uninsured_s + (1|zip),
  data = spatial_data,
  family = beta_binomial(),
  prior = priors,
  chains = 4, iter = 2000, cores = 4,
  control = list(adapt_delta = 0.99)
)

model_car <- brm(
  ptsd_cases | trials(trials) ~ 
    log_violent + hardship_s + uninsured_s +
    car(W, gr = zip),
  data = spatial_data,
  data2 = list(W = W),
  family = beta_binomial(),
  prior = priors,
  chains = 4, iter = 2000, cores = 4,
  control = list(adapt_delta = 0.99)
)

spatial_data$ptsd_est_prob_ctrl <- fitted(model_controls)[,1] / spatial_data$trials

ggplot(spatial_data, aes(x = violent_rate, y = ptsd_est_prob_ctrl)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Posterior Estimated PTSD Probability (Controls Model)",
    x = "Violent Crime Rate",
    y = "Estimated PTSD Probability"
  )
