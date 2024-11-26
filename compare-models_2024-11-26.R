#Compare a model based on taxa to one based on genes
#Taxa from Spencer et al Science 2021 - beneficial effect
#of Ruminococcacaeae, Clostridia, and Oscillospirales

library(tidyverse)
library(phyloseq)
library(yingtools2)
library(survival)
library(rms)


result <- ps %>%
  get.otu.melt() %>%
  select(sample, age10, cohort, ecog, tt_pfs_d, tt_os_d, pfs_event, event_os, Class, Order, Family, Species, numseqs) %>%
  group_by(sample) %>%
  summarise(
    rum = sum(numseqs[Family == "Ruminococcaceae"], na.rm=T),
    rum_lach = sum(numseqs[Family %in% c("Ruminococcaceae", "Lachnospiraceae")], na.rm=T),
    clostridiales = sum(numseqs[Order == "Clostridiales"], na.rm=T),
    clostridia = sum(numseqs[Class == "Clostridia"], na.rm=T),
    faecalibacterium_prausnitzii = sum(numseqs[Species == "Faecalibacterium_prausnitzii"], na.rm=T)
  ) %>%
  mutate(
    rum = replace_na(rum, 0),
    rum_lach = replace_na(rum_lach, 0),
    clostridiales = replace_na(clostridiales, 0),
    clostridia = replace_na(clostridia, 0),
    faecalibacterium_prausnitzii = replace_na(faecalibacterium_prausnitzii, 0)
  ) %>%
  left_join(tpt1_a, by = c("sample" = "WMS_SGPID"))


model1os <- coxph(Surv(tt_os_d, event_os)~log(rum), result)
model2os <- coxph(Surv(tt_os_d, event_os)~log(rum_lach), result)
model3os <- coxph(Surv(tt_os_d, event_os)~log(clostridiales), result)
model4os <- coxph(Surv(tt_os_d, event_os)~log(clostridia), result)
model5os <- coxph(Surv(tt_os_d, event_os)~log(eval(faecalibacterium_prausnitzii+1e-4)), result)
model6os <- coxph(Surv(tt_os_d, event_os)~log1p(pyruvate), result)





model1s <- coxph(Surv(tt_pfs_d, pfs_event)~log(rum), result)
model2s <- coxph(Surv(tt_pfs_d, pfs_event)~log(rum_lach), result)
model3s <- coxph(Surv(tt_pfs_d, pfs_event)~log(clostridiales), result)
model4s <- coxph(Surv(tt_pfs_d, pfs_event)~log(clostridia), result)
model5s <- coxph(Surv(tt_pfs_d, pfs_event)~log(eval(faecalibacterium_prausnitzii+1e-4)), result)
model6s <- coxph(Surv(tt_pfs_d, pfs_event)~log1p(pyruvate), result)


# Assuming you have a list of Cox models named model1, model2, ..., model6
models <- list(model1, model2, model3, model4, model5, model6)

# Initialize an empty data frame to store results
results_df <- data.frame(
  Model = character(),
  HR = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each model and extract HR, CI, and p-value
for (i in seq_along(models)) {
  model <- models[[i]]
  summary_model <- summary(model)

  # Extract the HR, CI, and p-value for the main variable of interest
  hr <- summary_model$coefficients[, "exp(coef)"]
  lower_ci <- summary_model$conf.int[, "lower .95"]
  upper_ci <- summary_model$conf.int[, "upper .95"]
  pvalue <- summary_model$coefficients[, "Pr(>|z|)"]

  # Add to results data frame
  results_df <- rbind(results_df, data.frame(
    Model = paste0("model", i),
    HR = hr,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    pvalue = pvalue
  ))
}

# Print the results
print(results_df)

ggplot(results_df) +
  geom_point(aes(y=fct_rev(fct_reorder(Model, HR)), x=HR)) +
  geom_errorbarh(aes(y=fct_rev(fct_reorder(Model, HR)), xmin=Lower_CI, xmax=Upper_CI), height=0.5)+
  geom_vline(xintercept=1, lty=2, col="gray70") +
  scale_x_continuous(trans="log") +
  theme_classic() +
  scale_y_discrete(labels=rev(c(
  "log(ACoA)",
  "log(Ruminococacceae + Lachnospiraceae)",
  "log(Clostridiales)",
  "log(Clostridia)",
  "log(Faecalibacterium prausnitizii + 1e-7)",
  "log(Ruminococcaceae)"
  ))) +
  labs(title = "HR (95% CI) univariable Cox -- our data",
       subtitle ="", y="Variable", x="HR")



model1 <- coxph(Surv(tt_os_d, event_os)~log(rum)+age10+factor(cohort)+ecog, result)
model2 <- coxph(Surv(tt_os_d, event_os)~log(rum_lach)+age10+factor(cohort)+ecog, result)
model3 <- coxph(Surv(tt_os_d, event_os)~log(clostridiales)+age10+factor(cohort)+ecog, result)
model4 <- coxph(Surv(tt_os_d, event_os)~log(clostridia)+age10+factor(cohort)+ecog, result)
model5 <- coxph(Surv(tt_os_d, event_os)~log(faecalibacterium_prausnitzii+0.000001)+age10+factor(cohort)+ecog, result)
model6 <- coxph(Surv(tt_os_d, event_os)~log(pyruvate)+age10+factor(cohort)+ecog, result)

result <- result[which(!is.na(result$ecog)),]

model1 <- coxph(Surv(tt_pfs_d, pfs_event)~log(rum)+age10+factor(cohort)+ecog, result)
model2 <- coxph(Surv(tt_pfs_d, pfs_event)~log(rum_lach)+age10+factor(cohort)+ecog, result)
model3 <- coxph(Surv(tt_pfs_d, pfs_event)~log(clostridiales)+age10+factor(cohort)+ecog, result)
model4 <- coxph(Surv(tt_pfs_d, pfs_event)~log(clostridia)+age10+factor(cohort)+ecog, result)
model5 <- coxph(Surv(tt_pfs_d, pfs_event)~log(faecalibacterium_prausnitzii+0.000001)+age10+factor(cohort)+ecog, result)
model6 <- coxph(Surv(tt_pfs_d, pfs_event)~log(pyruvate)+age10+factor(cohort)+ecog, result)

library(timeROC)

model1 <- coxph(Surv(tt_os_d, event_os)~log(rum), result)
model2 <- coxph(Surv(tt_os_d, event_os)~log(rum_lach), result)
model3 <- coxph(Surv(tt_os_d, event_os)~log(clostridiales), result)
model4 <- coxph(Surv(tt_os_d, event_os)~log(clostridia), result)
model5 <- coxph(Surv(tt_os_d, event_os)~log(faecalibacterium_prausnitzii+0.000001), result)
model6 <- coxph(Surv(tt_os_d, event_os)~log(pyruvate), result)
