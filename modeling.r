library(dplyr)
library(brms)
library(ggplot2)
library(brms)


crime_ptsd_clean <- chicago_crime_ptsd %>%
  mutate(
    ptsd_cases = round(ptsd_rate * population),
    log_violent = log(violent_rate) 
  )

crime_ptsd_final <- chicago_crime_ptsd %>%
  filter(!zip %in% c("60601", "60602", "60603", "60604")) %>%
  mutate(
    ptsd_cases = round(ptsd_rate * population),
    log_violent = log(violent_rate),
    hardship_s = scale(hardship_index)[,1],
    uninsured_s = scale(uninsured_rate)[,1],
    income_s = scale(med_income)[,1]
  )

#Intercept from on Schein et al.
#Using beta binomial distr to account for variance from geographically
#localied concentration mentioned in (Weisburd et al., 2018)
model_final <- brm(
  ptsd_cases | trials(population) ~ log_violent + hardship_s + uninsured_s,
  data = crime_ptsd_final,
  family = beta_binomial(),
  prior = c(
    prior(normal(-3, 0.5), class = Intercept), 
    prior(normal(0, 1), class = b),
    prior(exponential(1), class = phi)
  ),
  chains = 4, iter = 2000, cores = 4,
  control = list(adapt_delta = 0.99)
)

summary(model_final)
pp_check(model_final)

crime_ptsd_clean$ptsd_est <- fitted(model_final)[, 1]
ggplot(crime_ptsd_clean, aes(x = violent_rate, y = ptsd_est / population))
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))
  labs(
    title = "Corrected Posterior PTSD Probability",
    x = "Violent Crime Rate",
    y = "PTSD Probability (Estimated)"
  )

plot(conditional_effects(model_final), points = TRUE)
