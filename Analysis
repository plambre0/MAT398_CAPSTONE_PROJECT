library(sf)
library(spdep)
library(dplyr)
library(brms)
library(ggplot2)
library(scales)
library(car)
library(hglm)

summary(spatial_data_scaled)
hist(spatial_data_scaled$med_income_s, breaks = 20)
hist(spatial_data_scaled$uninsured_rate_s, breaks = 20)
hist(spatial_data_scaled$pm25_s, breaks = 20)
hist(spatial_data_scaled$hs_grad_rate_s, breaks = 20)
hist(spatial_data_scaled$unemployment_rate_s, breaks = 20)
hist(spatial_data_scaled$hardship_index_s, breaks = 20)
hist(spatial_data_scaled$prop_white_s, breaks = 20)
hist(spatial_data_scaled$log_violent_s, breaks = 20)

cor(spatial_data_scaled %>% select(contains("_s")),
    use="complete.obs", method = "kendall")
pheatmap::pheatmap(cor(spatial_data_scaled %>% select(contains("_s"))))

no_contr_bin <- glm(cbind(ptsd_cases, population - ptsd_cases) ~ 
                      log_violent_s,
                    data = spatial_data_scaled,
                    family = binomial())
summary(no_contr_bin)

contr_bin <- glm(cbind(ptsd_cases, population - ptsd_cases) ~ 
                   log_violent_s + med_income_s + uninsured_rate_s + 
                   hs_grad_rate_s + unemployment_rate_s +
                   lat_s * long_s,
                 data = spatial_data_scaled,
                 family = binomial())
summary(contr_bin)
vif(contr_bin)

priors_informed <- c(
  set_prior("normal(-5, 1)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b"),
  set_prior("exponential(1)", class = "phi")
)

brm_no_contr <- brm(ptsd_cases | trials(population) ~ log_violent_s,
                    data = spatial_data_scaled,
                    family = beta_binomial(),
                    prior = priors_informed,
                    chains = 4, iter = 2000, cores = 4,
                    control = list(adapt_delta = 0.99))
plot(brm_no_contr)
summary(brm_no_contr)
nuts_params(brm_no_contr)
posterior_summary(brm_no_contr, variable = "b_log_violent_s")
loo(brm_no_contr)
pp_check(brm_no_contr, type = "dens_overlay")
pp_check(brm_no_contr, type = "stat")
bayes_R2(brm_no_contr)
exp(posterior_summary(brm_no_contr, "b_log_violent_s"))
conditional_effects(brm_no_contr, "log_violent_s")

brm_contr_no_spat <- brm(ptsd_cases | trials(population) ~ 
                           log_violent_s + med_income_s + uninsured_rate_s + 
                           hs_grad_rate_s + unemployment_rate_s,
                         data = spatial_data_scaled,
                         family = beta_binomial(),
                         prior = priors_informed,
                         chains = 4, iter = 2000, cores = 4,
                         control = list(adapt_delta = 0.99))
plot(brm_contr_no_spat)
summary(brm_contr_no_spat)
nuts_params(brm_contr_no_spat)
posterior_summary(brm_contr_no_spat, variable = "b_log_violent_s")
loo(brm_contr_no_spat)
pp_check(brm_contr_no_spat, type = "dens_overlay")
pp_check(brm_contr_no_spat, type = "stat")
bayes_R2(brm_contr_no_spat)
exp(posterior_summary(brm_contr_no_spat, "b_log_violent_s"))
conditional_effects(brm_contr_no_spat, "log_violent_s")

brm_contr_spat <- brm(
  ptsd_cases | trials(population) ~ log_violent_s + med_income_s +
    uninsured_rate_s + hs_grad_rate_s + unemployment_rate_s +
    gp(lat, long),
  data = spatial_data_scaled,
  family = beta_binomial(),
  chains = 2, iter = 5000, cores = 4,
  prior = priors_informed,
  control = list(adapt_delta = 0.999)
)
plot(brm_contr_spat)
summary(brm_contr_spat)
nuts_params(brm_contr_spat)
posterior_summary(brm_contr_spat, variable = "b_log_violent_s")
loo(brm_contr_spat)
pp_check(brm_contr_spat, type = "dens_overlay")
pp_check(brm_contr_spat, type = "stat")
bayes_R2(brm_contr_spat)
exp(posterior_summary(brm_contr_spat, "b_log_violent_s"))
conditional_effects(brm_contr_spat, "log_violent_s")

model_no_contr_spat <- brm(
  formula = ptsd_cases | trials(population) ~ log_violent_s + t2(lat_s, long_s),
  data = spatial_data_scaled,
  family = beta_binomial(),
  prior = priors_informed
)
plot(model_no_contr_spat)
summary(model_no_contr_spat)
nuts_params(model_no_contr_spat)
posterior_summary(model_no_contr_spat, variable = "b_log_violent_s")
loo(model_no_contr_spat)
pp_check(model_no_contr_spat, type = "dens_overlay")
pp_check(model_no_contr_spat, type = "stat")
bayes_R2(model_no_contr_spat)
exp(posterior_summary(model_no_contr_spat, "b_log_violent_s"))
conditional_effects(model_no_contr_spat, "log_violent_s")

model_car <- brm(
  formula = ptsd_cases | trials(population) ~ 
    log_violent_s + uninsured_rate_s + hs_grad_rate_s + 
    car(adj_matrix, gr = zip, type = "icar"), 
  data = spatial_data_scaled,
  data2 = list(adj_matrix = adj_matrix), 
  family = beta_binomial(),
  prior = priors_informed,
  chains = 4, 
  iter = 2000, 
  cores = 4,
  control = list(adapt_delta = 0.97,
                 max_treedepth = 15)
)
plot(model_car)
summary(model_car)
nuts_params(model_car)
posterior_summary(model_car, variable = "b_log_violent_s")
loo(model_car)
pp_check(model_car, type = "dens_overlay")
pp_check(model_car, type = "stat")
bayes_R2(model_car)
model_car_res <- residuals(contr_bin, type="pearson")
moran.test(model_car_res, )
exp(posterior_summary(model_car, "b_log_violent_s"))
conditional_effects(model_car, "log_violent_s")

loo_compare(
  loo(brm_no_contr),
  loo(brm_contr_no_spat),
  loo(brm_contr_spat),
  loo(model_no_contr_spat),
  loo(model_car)
)
