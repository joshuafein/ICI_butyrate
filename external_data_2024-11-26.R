library(prodlim)
library(readxl)
library(tidyverse)
library(gtable)
library(gtsummary)
library(vegan)
library(compositions)


pyruvate <- read_xlsx("~/Library/CloudStorage/Box-Box/Joshua Fein (Internal)/ICI_butyrate/genelist.xlsx", sheet = 1)
glutarate <- read_xlsx("~/Library/CloudStorage/Box-Box/Joshua Fein (Internal)/ICI_butyrate/genelist.xlsx", sheet = 2)
fouraminobutyrate <- read_xlsx("~/Library/CloudStorage/Box-Box/Joshua Fein (Internal)/ICI_butyrate/genelist.xlsx", sheet = 3)
lysine <- read_xlsx("~/Library/CloudStorage/Box-Box/Joshua Fein (Internal)/ICI_butyrate/genelist.xlsx", sheet = 4)

l_pyruvate <- unlist(pyruvate[,1:10], use.names=F)
l_pyruvate <- l_pyruvate[which(!is.na(l_pyruvate))]


l_glutarate <- unlist(glutarate[,1:7], use.names=F)
l_glutarate <- l_glutarate[which(!is.na(l_glutarate))]

l_fouraminobutyrate <- unlist(fouraminobutyrate[,1:3], use.names=F)
l_fouraminobutyrate <- l_fouraminobutyrate[which(!is.na(l_fouraminobutyrate))]

l_lysine <- unlist(lysine[,1:8], use.names=F)
l_lysine <- l_lysine[which(!is.na(l_lysine))]


genelist <- data.frame(rbind(cbind(l_pyruvate, "pyruvate"),
                             cbind(l_glutarate, "gluatarate"),
                             cbind(l_fouraminobutyrate, "fouraminobutyrate"),
                             cbind(l_lysine, "lysine")))
genelist <- genelist %>%
  rename("gene_family_id" = "l_pyruvate" ,"gene_path_true" = "V2") %>%
  mutate(gene_family_id = as.numeric(gene_family_id))


result_df <- result_df %>%
  left_join(genelist)



############
wargo <- read.delim("combined_butyrate_952024.tsv")
wargo <- wargo %>%
  rename(sample_id = sampleid)
result_df$gene_family_id <- as.numeric(result_df$gene_family_id)

named_wargo <- wargo %>%
  left_join(result_df, by=c("Family"="gene_family_id"))

as_tibble(named_wargo)
named_wargo %>%
  group_by(gene_path_true, sample_id) %>%
  summarise(total_count = sum(Count)) %>%
  filter(!is.na(gene_path_true) & gene_path_true != "NA" & gene_path_true != "other") %>%
  ggplot() +
  geom_point(aes(x=fct_reorder(sample_id, total_count, sum), y=total_count, fill=gene_path_true), shape=21, size=4, alpha=0.6)+
  scale_fill_manual(values=colors)+
  scale_y_continuous(trans="sqrt")+
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position="bottom", axis.text.y=element_text(size=18), legend.text=element_text(size=18))

genename_wargo <- named_wargo %>%
  mutate(gene_name = if_else(grepl("^[^(\\[]", gene_info), # If string doesn't start with ( or [
                             sub("\\s*\\[.*", "", gene_info), # Remove starting from [
                             gene_info)) %>%
  mutate(gene_name = if_else(grepl("^\\(", gene_name), # If cleaned string starts with (
                             sub("^\\(.*?\\)\\s*", "", gene_name), # Remove the first ()
                             gene_name)) %>%
  mutate(gene_name = sub("\\s*\\(.*", "", gene_name)) %>% # Finally, remove starting from
  mutate(gene_name = sub("\\s*[\\[\\(].*", "", gene_name)) %>%
  group_by(sample_id, gene_path_true) %>%
  dplyr::filter(gene_path_true != "other" & !is.na(gene_path_true)) %>%
  dplyr::summarise(total_count = sum(Count)) %>%
  ungroup() %>%
  pivot_wider(id_cols= sample_id, names_from=gene_path_true, values_from=total_count)

w_meta <- read_csv("SraRunTable (1).txt")


w_meta_s <- w_meta %>%
  select(Run, pfs_d, pfsevent, fiber_cat, DSQfib, treatment, AGE, adv_substage, Antibiotics) %>%
  filter(Run %in% genename_wargo$sample_id)

pfs_wargo <- genename_wargo %>%
  left_join(w_meta_s, by=c("sample_id"="Run"))


w_meta_3 <- w_meta %>%
  filter(treatment %in% c("anti-PD1", "other ICB") &
           env_broad_scale == "Human" &
           !is.na(fiber_cat)) %>%
  select(Run, pfs_d, pfsevent, fiber_cat, DSQfib, treatment, AGE, adv_substage, Antibiotics)

pfs_wargo_2 <- genename_wargo %>%
  inner_join(w_meta_3, by=c("sample_id"="Run"))

plot(prodlim(Hist(pfs_d/30.44, pfsevent)~fiber_cat, pfs_wargo_2), xlim=c(0,60),
     confint=F, logrank=T, legend.cex=0.8, legend.x="bottomleft",
     legend.legend=c("Insufficient","Sufficient"))


  pfs_wargo_2$AGE[which(pfs_wargo_2$AGE == ">89")] <- 90
  pfs_wargo_2$AGE <- as.numeric(pfs_wargo_2$AGE)
#any ICB -- 93 patients

tbl_regression(coxph(Surv(pfs_d, pfsevent)~log1p(pyruvate), pfs_wargo_2), exponentiate = T) %>%
    modify_caption(caption="Progression-Free Survival") #p = 0.023, HR 0.79 [0.64, 0.97] for log increase

tbl_regression(coxph(Surv(pfs_d, pfsevent)~log1p(pyruvate)*Antibiotics, pfs_wargo_2), exponentiate = T) %>%
  modify_caption(caption="Progression-Free Survival") #p = 0.023, HR 0.79 [0.64, 0.97] for log increase


tbl_regression(coxph(Surv(pfs_d, pfsevent)~log1p(pyruvate) + fiber_cat, pfs_wargo_2), exponentiate = T) %>%
  modify_caption(caption="Progression-Free Survival") #p = 0.023, HR 0.79 [0.64, 0.97] for log increase

tbl_regression(coxph(Surv(pfs_d, pfsevent)~log1p(pyruvate)+DSQfib10, pfs_wargo_2), exponentiate = T) %>%
  modify_caption(caption="Progression-Free Survival") #p = 0.023, HR 0.79 [0.64, 0.97] for log increase


tbl_regression(coxph(Surv(pfs_d, pfsevent)~log1p(pyruvate)+adv_substage, pfs_wargo_2), exponentiate = T) %>%
  modify_caption(caption="Progression-Free Survival") #p = 0.072, HR 0.67 [0.67, 1.02] for log increase

boxplot(log1p(pfs_wargo_2$pyruvate)~pfs_wargo_2$fiber_cat)


tbl_regression(coxph(Surv(pfs_d, pfsevent)~DSQfib10, pfs_wargo_2), exponentiate = T) %>%
  modify_caption(caption="Progression-Free Survival") #p = 0.3, HR 0.95 [0.86, 1.04] for 1 DSQ increase

tbl_regression(coxph(Surv(pfs_d, pfsevent)~fiber_cat, pfs_wargo_2), exponentiate = T) %>%
  modify_caption(caption="Progression-Free Survival") #p = 0.2, HR 0.59 [0.25, 1.40]



library(prodlim)
library(survival)
library(gtsummary)
tbl_regression(coxph(Surv(pfs_d, pfsevent)~log1p(pyruvate), pfs_wargo), exponentiate = T) %>%
  modify_caption(caption="Progression-Free Survival") #p = 0.029, HR 0.79 [0.64, 0.98] for log increase

pfs_wargo <- pfs_wargo %>%
  mutate(med_p = cut(pyruvate, quantile(pyruvate, probs=seq(0,1,0.5)), c("below", "above"))) %>%
  mutate(med_p250 = cut(pyruvate, c(-Inf, 250, Inf)))

plot(prodlim(Hist(pfs_d/30.44, pfsevent)~med_p250, pfs_wargo), confint=T, marktime=T, xlim=c(0,48), axes=F, background.fg=NA,  atrisk.labelcol="black", atrisk.labels=levels(pfs_wargo$med_p),
     background.border=NA, atrisk=T, xlab="Months", ylab="Progression-Free Survival, %", legend=T, legend.x="topright", legend.title="", legend.cex=0.6, atrisk.at=seq(0,48, 12), atrisk.title="", atrisk.cex=0.8, lwd=1)
axis(side=1, pos=c(0,0), lwd=2, at=seq(0,48, 12))
axis(side=2, pos=0, lwd=2, at=seq(0,1, 0.25), labels=c("0", "25", "50", "75", "100"))
round(survdiff(Surv(pfs_d, pfsevent)~med_p, pfs_wargo)$pvalue,8) #p = 0.648

library(microViz)
library(viridis)
library(ggExtra)
library(tidyverse)


raw_mpa_phy <- readRDS("out_phyloseq80.RDS")

sampletable <- sample_data(raw_mpa_phy)
sampletable <- sampletable %>%
  as.data.frame()
sampletable$sample_id <- rownames(sampletable)
sampletable <- sampletable %>%
  as_tibble() %>%
  left_join(pfs_wargo, by = "sample_id") %>%
  select(-dummymetadata) %>%
  mutate(log_acoa = log(pyruvate))

sampletable <- sampletable %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

sample_data(raw_mpa_phy) <- sample_data(sampletable)


taxa_table <- as.data.frame(tax_table(raw_mpa_phy))
taxa_table <- taxa_table %>%
  mutate(Family = if_else(Order == "Bacteroidales", "Bacteroidales (order)", Family))
tax_table(raw_mpa_phy) <- tax_table(as.matrix(taxa_table))
raw_mpa_phy_icb <- subset_samples(raw_mpa_phy, treatment != "other systemic", )

common_min <- -0.9480444
common_max <- 7.9995015

primary_min <- 4.318380
primary_max <- 7.512117

#species pcoa
plot_ord_sp_wargo <- raw_mpa_phy_icb %>%
  filter(!is.infinite(log_acoa)) %>%
  tax_fix() %>%
  tax_transform(trans = "compositional", rank = "Species") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(fill = "log_acoa", size = 4, shape=21, alpha=0.6,
           plot_taxa = c(1:5)) +
  #scale_fill_viridis(option="B", begin = 0, end=1) +
  theme_classic()

library(scales)
pdf("figures/PCoA_species_acoa_wargo-25Pct24.pdf", height=6, width=6)
plot_ord_sp_wargo$data$scaled_fill <- rescale(plot_ord_sp_wargo$data$log_acoa, to = c(primary_min, primary_max), from = c(common_min, common_max))
plot_ord_sp_wargo + scale_fill_viridis(option="B", limits = c(primary_min, primary_max), oob = scales::squish)
dev.off()
summary(lm(plot_ord_sp_wargo$data$MDS1~plot_ord_sp_wargo$data$log_acoa))

ord_plot <- raw_mpa_phy_icb %>%
  tax_fix() %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ord_calc("PCA") %>%
  ord_plot(fill = "log_acoa", size = 4, shape=21, alpha=0.6,
           plot_taxa = 1:3) + #c("Lachnospiraceae", "Bacteroidaceae", "Ruminococcaceae")) +
  scale_fill_viridis(option="C", begin = 0, end=1, na.value="white") +
  theme_classic()




# Aitchinson distance
physeq_otu <- as(otu_table(raw_mpa_phy_icb), "matrix")
physeq_otu[physeq_otu == 0] <- 0.00001  # Replace zeros for log-transformation
clr_transformed <- apply(physeq_otu, 1, function(x) compositions::clr(x))
aitchison_dist <- vegdist(t(clr_transformed), method = "euclidean")
pcoa_results <- ordinate(raw_mpa_phy_icb, method = "PCoA", distance = aitchison_dist)
eigenvalues <- pcoa_results$values$Eigenvalues
variance_explained <- eigenvalues / sum(eigenvalues)

pdf("supp_6.16Sep2024.pdf", height=5, width=5)
ord_plot
dev.off()

View(ord_plot$data)

summary(with(ord_plot$data[which(ord_plot$data$log_acoa != -Inf & !is.na(ord_plot$data$log_acoa)
                                 & !is.nan(ord_plot$data$PC1)),], lm(PC1~log_acoa)))
#R^2 0.088
#p log_acoa 0.0001
#B log_acoa -0.04

##################



#Compare a model based on taxa to one based on genes
#Taxa from Spencer et al Science 2021 - beneficial effect
#of Ruminococcacaeae, Clostridia, and Oscillospirales

library(tidyverse)
library(phyloseq)
library(yingtools2)
library(survival)
library(rms)


result <- raw_mpa_phy_icb %>%
  get.otu.melt() %>%
  select(sample, pfs_d, pfsevent, pyruvate, Class, Order, Family, Genus, Species, numseqs) %>%
  group_by(sample, pfs_d, pfsevent, pyruvate) %>%
  summarise(
    rum = sum(numseqs[Family == "Ruminococcaceae"], na.rm=T),
    rum_lach = sum(numseqs[Family %in% c("Ruminococcaceae", "Lachnospiraceae")], na.rm=T),
    clostridiales = sum(numseqs[Order == "Clostridiales"], na.rm=T),
    osc = sum(numseqs[Order == "Oscillospirales"], na.rm=T),
    clostridia = sum(numseqs[Class == "Clostridia"], na.rm=T),
    faecalibacterium = sum(numseqs[Genus == "Faecalibacterium"], na.rm=T),
    faecalibacterium_prausnitzii = sum(numseqs[Species == "Faecalibacterium_prausnitzii"], na.rm=T)
  ) %>%
  mutate(
    rum = replace_na(rum, 0),
    osc = replace_na(osc, 0),
    rum_lach = replace_na(rum_lach, 0),
    clostridiales = replace_na(clostridiales, 0),
    clostridia = replace_na(clostridia, 0),
    faecalibacterium = replace_na(faecalibacterium, 0),
    faecalibacterium_prausnitzii = replace_na(faecalibacterium_prausnitzii, 0)
  )

model1 <- coxph(Surv(pfs_d, pfsevent)~log(rum), result)
model2 <- coxph(Surv(pfs_d, pfsevent)~log(rum_lach), result)
model3 <- coxph(Surv(pfs_d, pfsevent)~log(clostridiales), result)
model4 <- coxph(Surv(pfs_d, pfsevent)~log(clostridia), result)
model5 <- coxph(Surv(pfs_d, pfsevent)~log(eval(faecalibacterium_prausnitzii+1e-7)), result)
model6 <- coxph(Surv(pfs_d, pfsevent)~log1p(pyruvate), result)
model7 <- coxph(Surv(pfs_d, pfsevent)~log1p(faecalibacterium), result)
model8 <- coxph(Surv(pfs_d, pfsevent)~log1p(osc), result)


models <- list(model1, model2, model3, model4, model5, model6,
               model1s, model2s, model3s, model4s, model5s, model6s)

osmodels <- list(model1os, model2os, model3os, model4os, model5os, model6os)
# Initialize an empty data frame to store results
results_df <- data.frame(
  Model = character(),
  HR = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)


for (i in seq_along(models)) {
  model <- models[[i]]
  summary_model <- summary(model)


  hr <- summary_model$coefficients[, "exp(coef)"]
  lower_ci <- summary_model$conf.int[, "lower .95"]
  upper_ci <- summary_model$conf.int[, "upper .95"]
  pvalue <- summary_model$coefficients[, "Pr(>|z|)"]


  results_df <- rbind(results_df, data.frame(
    Model = paste0("model", i),
    HR = hr,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    pvalue = pvalue
  ))
}

results_df_os <- data.frame(
  Model = character(),
  HR = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)


for (i in seq_along(osmodels)) {
  model <- osmodels[[i]]
  summary_model <- summary(model)


  hr <- summary_model$coefficients[, "exp(coef)"]
  lower_ci <- summary_model$conf.int[, "lower .95"]
  upper_ci <- summary_model$conf.int[, "upper .95"]
  pvalue <- summary_model$coefficients[, "Pr(>|z|)"]


  results_df_os <- rbind(results_df_os, data.frame(
    Model = paste0("model", i),
    HR = hr,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    pvalue = pvalue
  ))
}


results_df$source = c(rep("wargo", 6), rep("seres",6))
results_df$modelname = rep(c("log(Ruminococcaceae)", "log(Ruminococcacaeae + Lachnospiraceae",
                             "log(Clostridiales)", "log(Clostridia)", "log(F prausnitizii + 1e-7)", "log(ACoA"),2)
pdf("acoa/figures/univar_wargo.pdf", height=4, width=8)
ggplot(results_df) +
  geom_point(aes(y=fct_rev(fct_reorder(modelname, HR)), x=HR)) +
  geom_errorbarh(aes(y=fct_rev(fct_reorder(modelname, HR)), xmin=Lower_CI, xmax=Upper_CI), height=0.5)+
  geom_vline(xintercept=1, lty=2, col="gray70") +
  scale_x_continuous(trans="log") +
  theme_classic() +
  facet_grid(.~source) +
  labs(title = "HR (95% CI) univariable Cox",
       subtitle ="", y="Variable", x="HR")
dev.off()

