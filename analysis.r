library(car)
library(glmmTMB)
library(tmap)
library(spdep)
library(brms)
library(tidybayes)
library(lavaan)
library(lavaanPlot)
library(DHARMa)
library(performance)
library(lmtest)

################################################################################
#Preprocessing
################################################################################

# Density plots
lapply(names(chicago_scaled), function(col) {
  plot(density(chicago_scaled[[col]]), main = col)
})

# Cube-root transformation
vars_transform <- c("loc_afford_idx","active_transport","poverty_rate",
                    "unemployment_rate","hs_grad_rate","food_insecurity_rate",
                    "violent_crime_rate","drug_crime_rate","ptsd_rate")

chicago_transf <- chicago_scaled
chicago_transf[, vars_transform] <- apply(chicago_scaled[, vars_transform], 2, 
                                          function(x) sign(x) * abs(x)^(1/3))

# Correlations
cor(chicago_scaled)
cor(chicago_transf)

cor(chicago_transf[, c("ptsd_rate","violent_crime_rate","poverty_rate",
                       "hs_grad_rate","white_pct","foreign_born",
                       "drug_crime_rate")])

################################################################################
#PCA/FA
################################################################################

#PCA
pca_model <- prcomp(chicago_transf)
plot(pca_model)

#Factor Analysis
fa_model <- factanal(chicago_transf[, c(
  "loc_afford_idx","active_transport","foreign_born","white_pct","population",
  "poverty_rate","house_cost_burden","unemployment_rate","hs_grad_rate",
  "food_insecurity_rate","medicaid_coverage","violent_crime_rate",
  "drug_crime_rate")],
  factors = 5)

################################################################################
#CFA
################################################################################

cfa_spec <- '
  Hardship =~ poverty_rate + medicaid_coverage + loc_afford_idx
  Safety   =~ violent_crime_rate + unemployment_rate + active_transport + 
  white_pct
'

cfa_fit <- lavaan::cfa(cfa_spec, data = chicago_transf)
summary(cfa_fit)
standardizedSolution(cfa_fit)
fitMeasures(cfa_fit, c("cfi","tli","rmsea","srmr","aic","bic"))

################################################################################
#Linear Models ~ PTSD Rate
################################################################################

#No controls
lm_ptsd_crime <- glm(ptsd_rate ~ violent_crime_rate,
                     data = chicago_transf, family = gaussian())

#With controls
lm_ptsd_full <- glm(ptsd_rate ~ violent_crime_rate + poverty_rate +
                      foreign_born + white_pct,
                    data = chicago_transf, family = gaussian())

summary(lm_ptsd_full)
vif(lm_ptsd_full)

# Diagnostics
par(mfrow = c(2,2)); plot(lm_ptsd_full)
qqnorm(residuals(lm_ptsd_full)); qqline(residuals(lm_ptsd_full))
shapiro.test(residuals(lm_ptsd_full))
bptest(lm_ptsd_full)
summary(influence.measures(lm_ptsd_full))

cooksd <- cooks.distance(lm_ptsd_full)
plot(cooksd, type = "h"); abline(h = 4/length(cooksd), col = "red")

################################################################################
#Count Models ~ Pos Responses
################################################################################

# Poisson baseline
pois_ptsd_base <- glm(pos_res ~ white_pct + poverty_rate,
                      offset = log(population),
                      family = poisson(),
                      data = chicago_transf)

# Negative binomial
nb_ptsd_base <- glmmTMB(pos_res ~ white_pct + poverty_rate +
                          offset(log(population)),
                        family = nbinom2,
                        data = chicago_transf)

summary(nb_ptsd_base)

# DHARMa diagnostics
res_nb <- simulateResiduals(nb_ptsd_base)
plot(res_nb)
testDispersion(res_nb)
testZeroInflation(res_nb)
testUniformity(res_nb)

################################################################################
#Spatial Data Setup
################################################################################

chicago_all <- chicago_transf
chicago_all$pos_res <- chicago$pos_res
chicago_all$population_raw <- chicago$population
chicago_all$zip <- chicago$zip
chicago_all$pos <- numFactor(zip_chicago$lat, zip_chicago$long)

################################################################################
#Mixed & Spatial Models
################################################################################

#Random intercept model
nb_ptsd_mixed <- glmmTMB(
  pos_res ~ violent_crime_rate + poverty_rate + foreign_born + white_pct +
    offset(log(population_raw)) + (1 | zip),
  data = chicago_all,
  family = nbinom2
)

summary(nb_ptsd_mixed)

res_mixed <- simulateResiduals(nb_ptsd_mixed)
plot(res_mixed)

AIC(nb_ptsd_base, nb_ptsd_mixed)
BIC(nb_ptsd_base, nb_ptsd_mixed)


#Spatial Weights
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(W) <- chicago_all$zip
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

#Bayesian Spatial Models
priors_spatial <- c(
  prior(normal(-6, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sdcar"),
  prior(beta(2, 2), class = "rhocar")
)

bayes_spatial_nb <- brm(
  pos_res ~ violent_crime_rate + poverty_rate + foreign_born + white_pct +
    offset(log(population_raw)) + car(W, gr = zip, type = "bym2"),
  data = chicago_all,
  data2 = list(W = W),
  family = negbinomial(),
  prior = priors_spatial,
  chains = 4, iter = 4000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(bayes_spatial_nb)
pp_check(bayes_spatial_nb)
loo(bayes_spatial_nb)
waic(bayes_spatial_nb)


#Spatial Residual Analysis
mu <- predict(bayes_spatial_nb)[,1]
chicago_all$resid <- (chicago_all$pos_res - mu) /
  sqrt(mu + (mu^2 / 1254756901))

moran.test(chicago_all$resid, listw)
local_m <- localmoran(chicago_all$resid, listw)

#SEM (Mediation)
sem_spec <- '
  violent_crime_rate ~ a*poverty_rate + foreign_born
  ptsd_rate ~ b*violent_crime_rate + c_prime*poverty_rate + white_pct

  indirect := a*b
  total := c_prime + (a*b)
'

sem_fit <- lavaan::sem(sem_spec, data = chicago_all)
summary(sem_fit, standardized = TRUE, rsquare = TRUE)

lavaanPlot(sem_fit, coefs = TRUE, stand = TRUE)

#Bayesian Spatial SEM
# Model formulas
bf_violence <- bf(
  violent_crime_rate ~ poverty_rate + foreign_born +
    car(W, gr = zip, type = "bym2")
)

bf_ptsd <- bf(
  pos_res ~ violent_crime_rate + poverty_rate + white_pct +
    offset(log(population_raw)) +
    car(W, gr = zip, type = "bym2")
)

priors_sem <- c(
  #Violence
  prior(normal(0, 1), class = "b", resp = "violentcrimerate"),
  prior(exponential(1), class = "sdcar", resp = "violentcrimerate"),
  
  #PTSD
  prior(normal(0, 1), class = "b", resp = "posres"),
  prior(exponential(1), class = "sdcar", resp = "posres")
)

# Fit model
bayes_spatial_sem <- brm(
  bf_violence + bf_ptsd + set_rescor(FALSE),
  data = chicago_all,
  data2 = list(W = W),
  family = list(gaussian(), negbinomial()),
  prior = priors_sem,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(bayes_spatial_sem)
pp_check(bayes_spatial_sem)
loo(bayes_spatial_sem)
