library(car)
library(glmmTMB)
library(tmap)
library(spdep)
library(brms)
library(tidybayes)
library(lavaan)
library(lavaanPlot)


lapply(names(chicago_scaled), function(col_name) {
  plot(density(chicago_scaled[[col_name]]), main = col_name)
})

chicago_scaled_transf <- chicago_scaled
chicago_scaled_transf[, c("loc_afford_idx",
                          "active_transport",
                          "poverty_rate",
                          "unemployment_rate",
                          "hs_grad_rate",
                          "food_insecurity_rate",
                          "violent_crime_rate",
                          "drug_crime_rate",
                          "ptsd_rate")] <- 
  apply(chicago_scaled[, c("loc_afford_idx",
                           "active_transport",
                           "poverty_rate",
                           "unemployment_rate",
                           "hs_grad_rate",
                           "food_insecurity_rate",
                           "violent_crime_rate",
                           "drug_crime_rate",
                           "ptsd_rate")], 2, 
        function(x){sign(x) * abs(x)^(1/3)})

cor(chicago_scaled)
cor(chicago_scaled_transf)


cor(chicago_scaled_transf[, c("ptsd_rate", "violent_crime_rate", "poverty_rate",
                          "hs_grad_rate", "white_pct", "foreign_born", 
                          "drug_crime_rate")])

chicago_prcomp <- prcomp(chicago_scaled_transf)
plot(chicago_prcomp)
chicago_fa <- factanal(chicago_scaled_transf[, c("loc_afford_idx",
                                                 "active_transport",
                                                 "foreign_born",
                                                 "white_pct",
                                                 "population",
                                                 "poverty_rate",
                                                 "house_cost_burden",
                                                 "unemployment_rate",
                                                 "hs_grad_rate",
                                                 "food_insecurity_rate",
                                                 "medicaid_coverage",
                                                 "violent_crime_rate",
                                                 "drug_crime_rate")], 
                       factors = 5)
chicago_model <- '
  Hardship =~ poverty_rate + medicaid_coverage + loc_afford_idx
  Safety   =~ violent_crime_rate + unemployment_rate + active_transport + white_pct
'
chicago_model_fit <- lavaan::cfa(chicago_model, data = chicago_scaled_transf)
lavaan::summary(chicago_model_fit)
lavaan::standardizedSolution(chicago_model_fit)
lavaan::fitMeasures(chicago_model_fit, c("cfi", "tli", "rmsea", "srmr", "aic", "bic"))

glm_no_contr <- glm(ptsd_rate ~ violent_crime_rate, 
                    data = chicago_scaled_transf, 
                    family = gaussian())
summary(glm_no_contr)
glm_contr <- glm(ptsd_rate ~ violent_crime_rate + poverty_rate + foreign_born + white_pct,
                 data = chicago_scaled_transf, family = gaussian())
summary(glm_contr)
vif(glm_contr)
avPlots(glm_contr)
par(mfrow = c(2,2))
plot(glm_contr)
dev.off()
qqnorm(residuals(glm_contr))
qqline(residuals(glm_contr))
shapiro.test(residuals(glm_contr))
library(lmtest)
bptest(glm_contr)
summary(influence.measures(glm_contr))
plot(cooks.distance(glm_contr), type = "h")
abline(h = 4/length(cooks.distance(glm_contr)), col = "red")

summary(glm(ptsd_rate ~ violent_crime_rate + poverty_rate + 
              foreign_born + white_pct, 
            family = gaussian(), 
            data = chicago_scaled_transf[-c(3, 4, 32), ]))

glm_no_contr <- glm(ptsd_rate ~ poverty_rate, 
                    family = gaussian(), 
                    data = chicago_scaled_transf)
avPlots(glm_contr)
par(mfrow = c(2,2))
plot(glm_contr)
dev.off()
qqnorm(residuals(glm_contr))
qqline(residuals(glm_contr))
shapiro.test(residuals(glm_contr))
bptest(glm_contr)
summary(influence.measures(glm_contr))

model_participation <- glm(chicago$pos_res ~ white_pct + poverty_rate, 
                           family = poisson(), 
                           offset = log(chicago$population), 
                           data = chicago_scaled_transf)
summary(model_participation)
avPlots(glm_contr)
par(mfrow = c(2,2))
plot(glm_contr)
dev.off()
qqnorm(residuals(glm_contr))
qqline(residuals(glm_contr))
shapiro.test(residuals(glm_contr))
bptest(glm_contr)
summary(influence.measures(glm_contr))

model_participation_model_nb <- glmmTMB(chicago$pos_res ~ white_pct + poverty_rate + offset(log(chicago$population)), 
                                        family = nbinom2, 
                                        data = chicago_scaled_transf)
summary(model_participation_model_nb)
library(performance)
library(DHARMa)
simulationOutput <- DHARMa::simulateResiduals(model_participation_model_nb)
plot(simulationOutput)
DHARMa::testDispersion(simulationOutput)
DHARMa::testZeroInflation(simulationOutput)
DHARMa::testUniformity(simulationOutput)

chicago$zip[as.numeric(names(sort(residuals(model_participation_model_nb, 
                                            type = "deviance"))))]

summary(glm(chicago$pos_res ~ violent_crime_rate + poverty_rate + foreign_born + white_pct,
            offset = log(chicago$population), data = chicago_scaled_transf, family = poisson()))

chicago_all <- chicago_scaled_transf
chicago_all$pos_res <- chicago$pos_res
chicago_all$population_raw <- chicago$population
chicago_all$zip <- chicago$zip
chicago_all$pos <- numFactor(zip_chicago$lat, zip_chicago$long)

model_spatial <- glmmTMB(pos_res ~ violent_crime_rate + poverty_rate + 
                           foreign_born + white_pct + 
                           offset(log(population_raw)) +
                           exp(pos + 0 | zip),
                         data = chicago_all,
                         family = nbinom2)

model_final_clean <- glmmTMB(
  pos_res ~ violent_crime_rate + poverty_rate + foreign_born + white_pct + 
    offset(log(population_raw)) + 
    (1 | zip),
  data = chicago_all,
  family = nbinom2
)
simulationOutput <- DHARMa::simulateResiduals(model_final_clean)
plot(simulationOutput)
summary(model_final_clean)
ranef(model_final_clean)
VarCorr(model_final_clean)
AIC(model_participation_model_nb, model_final_clean)
BIC(model_participation_model_nb, model_final_clean)

W <- nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(W) <- chicago_all$zip

priors <- c(
  prior(normal(-6, 2), class = "Intercept"),
  prior(normal(0.3, 1), class = "b", coef = "violent_crime_rate"),
  prior(normal(-0.2, 1), class = "b", coef = "poverty_rate"),
  prior(normal(0.1, 1), class = "b", coef = "white_pct"),
  prior(normal(0.1, 1), class = "b", coef = "foreign_born"),
  prior(exponential(1), class = "sdcar"),
  prior(beta(2, 2), class = "rhocar")
)

model_bayesian_car <- brm(
  pos_res ~ violent_crime_rate + poverty_rate + foreign_born + white_pct + 
    offset(log(population_raw)) + 
    car(W, gr = zip, type = "bym2"), 
  data = chicago_all,
  data2 = list(W = W),
  family = negbinomial(),
  chains = 4, 
  iter = 2000, 
  cores = 4
)

model_bayesian_spatial_final <- brm(
  pos_res ~ violent_crime_rate + poverty_rate + foreign_born + white_pct + 
    offset(log(population_raw)) + 
    car(W, gr = zip, type = "bym2"), 
  data = chicago_all,
  data2 = list(W = W),
  family = negbinomial(),
  prior = priors,
  chains = 4, 
  iter = 4000,
  warmup = 2000,
  cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
pp_check(model_bayesian_spatial_final)
pp_check(model_bayesian_spatial_final, type = "stat")
loo(model_bayesian_spatial_final)
waic(model_bayesian_spatial_final)

spatial_draws <- as.matrix(model_bayesian_spatial_final, variable = "^rcar", regex = TRUE)
print(dim(spatial_draws)) 
prob_under <- colMeans(spatial_draws < 0)

loo_compare(
  loo(model_bayesian_car),
  loo(model_bayesian_spatial_final)
)

zip_levels <- levels(chicago_all$zip)
chicago_all$zip <- droplevels(chicago_all$zip)
zip_levels <- levels(chicago_all$zip)
print(length(zip_levels))
prob_df <- data.frame(zip = zip_levels, prob_under = prob_under)
zip_shapes_final <- merge(zip_shapes, prob_df, by = "zip")
head(prob_df[order(-prob_df$prob_under), ], 5)

mu <- predict(model_bayesian_spatial_final)[, 1]
chicago_all$resid <- (chicago_all$pos_res - mu) / sqrt(mu + (mu^2 / 1254756901)) 
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_result <- moran.test(chicago_all$resid, listw)
print(moran_result)
local_m <- localmoran(chicago_all$resid, listw)
chicago_all$quadrant <- "Insignificant"
m_res <- scale(chicago_all$resid)
lag_res <- scale(lag.listw(listw, chicago_all$resid))
chicago_all$quadrant[m_res > 0 & lag_res > 0 & local_m[,5] <= 0.05] <- "High-High"
chicago_all$quadrant[m_res < 0 & lag_res < 0 & local_m[,5] <= 0.05] <- "Low-Low"
chicago_all$quadrant[m_res > 0 & lag_res < 0 & local_m[,5] <= 0.05] <- "High-Low"
chicago_all$quadrant[m_res < 0 & lag_res > 0 & local_m[,5] <= 0.05] <- "Low-High"

zip_shapes_final <- merge(zip_shapes, chicago_all[, c("zip", "quadrant")], by = "zip")

tm_shape(zip_shapes_final) +
  tm_polygons("quadrant", 
              palette = c("High-High" = "red", "Low-Low" = "blue", 
                          "Low-High" = "lightblue", "High-Low" = "pink", 
                          "Insignificant" = "white"),
              title = "LISA Cluster Map (Residuals)")

model_sem <- '
  #Violence is predicted by SES
  violent_crime_rate ~ a*poverty_rate + foreign_born
  
  #PTSD is predicted by Violence and SES
  # We use the log-transformed or rate-adjusted PTSD screens
  ptsd_rate ~ b*violent_crime_rate + c_prime*poverty_rate + white_pct
  
  indirect := a * b
  total    := c_prime + (a * b)
'

fit_sem <- lavaan::sem(model_sem, data = chicago_all)
lavaan::summary(fit_sem, standardize = TRUE, rsquare = TRUE)

lavaanPlot(model = fit_sem, 
           labels = c(violent_crime_rate = "Violent Crime", 
                      poverty_rate = "Poverty Preve", 
                      white_pct = "Percent White", 
                      foreign_born = "Foreign Born",
                      pos_res = "PTSD Screens"),
           node_options = list(shape = "box", fontname = "Helvetica"), 
           edge_options = list(color = "grey"), 
           coefs = TRUE, 
           stand = TRUE, 
           stars = c("regress"))

bf_violence <- bf(violent_crime_rate ~ poverty_rate + foreign_born + car(W, gr = zip, type = "bym2"))
bf_ptsd     <- bf(pos_res ~ violent_crime_rate + poverty_rate + white_pct + 
                    offset(log(population_raw)) + car(W, gr = zip, type = "bym2"))

fit_spatial_sem <- brm(
  bf_violence + bf_ptsd + set_rescor(FALSE), 
  data = chicago_all,
  data2 = list(W = W),
  family = list(gaussian(), negbinomial()),
  prior = c(
    prior(normal(0, 1), class = "b", resp = "violentcrimerate"),
    prior(normal(0, 1), class = "b", resp = "posres")
  ),
  chains = 4, iter = 4000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

bf_violence <- bf(violent_crime_rate ~ poverty_rate + foreign_born + car(W, gr = zip, type = "bym2"))
bf_ptsd     <- bf(pos_res ~ violent_crime_rate + poverty_rate + white_pct + 
                    offset(log(population_raw)) + car(W, gr = zip, type = "bym2"))

priors_final <- c(
  # Predictors for Violence
  prior(normal(0, 1), class = "b", resp = "violentcrimerate"),
  #Predictors for PTSD
  prior(normal(0, 1), class = "b", resp = "posres"),
  #Spatial variance
  prior(exponential(1), class = "sdcar", resp = "violentcrimerate"),
  prior(exponential(1), class = "sdcar", resp = "posres")
)

fit_spatial_sem_final <- brm(
  bf_violence + bf_ptsd + set_rescor(FALSE), 
  data = chicago_all,
  data2 = list(W = W),
  family = list(gaussian(), negbinomial()),
  prior = priors_final,
  chains = 4, 
  iter = 4000, 
  warmup = 2000, 
  cores = 4,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)
