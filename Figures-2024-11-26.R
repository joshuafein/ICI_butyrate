library(rms)
library(survival)
library(prodlim)
library(tidyverse)
library(yingtools2)
library(ape)
library(ggtree)
library(tidyverse)
library(phyloseq)
library(yingtools2)
library(microViz)
library(viridis)

# Cox model w RCS spline
model_spline_os <- coxph(Surv(tt_os_d, event_os)~rcs(pyruvate,4)+age10+factor(cohort)+ecog, tpt1_a)
ptemp <- termplot(model_spline_os, se=T, plot=F)
buterm <- ptemp$pyruvate
center <- buterm$y[which(buterm$y == median(buterm$y))]
ytemp <- buterm$y + outer(buterm$se, c(0, -1.96, 1.96), '*')
exp_ytemp <- exp(ytemp - center)

spline_data <- data.frame(buty = buterm$x, Estimate = exp_ytemp[,1],
                          Lower = exp_ytemp[,2], Upper = exp_ytemp[,3])

pdf("figures/hr_spline.pdf", height=6, width=6)
ggplot(spline_data, aes(x = buty)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = Estimate), color = "blue") +
  geom_rug(sides = "b") +  # Add rug plot at the bottom ('b') of the plot
  scale_y_log10() +  # Log scale for y-axis
  labs(x = "Pyruvate-pathway genes present", y = "Hazard ratio") +
  theme_minimal()
dev.off()

tpt1_a <- tpt1_a %>%
  mutate(med_p = cut(pyruvate, quantile(pyruvate, probs=seq(0,1,0.5)), c("below", "above")))

#Figure 2A - OS
pdf("figures/os.pdf", height=6, width=4)
plot(prodlim(Hist(tt_os_d/30.44, event_os)~med_p, tpt1_a), confint=T, marktime=T, xlim=c(0,48), axes=F, background.fg=NA,  atrisk.labelcol="black", atrisk.labels=levels(tpt1_a$med_p),
     background.border=NA, atrisk=T, xlab="Months", ylab="Overall Survival, %", legend=F, legend.x="topright", legend.title="", legend.cex=0.6, atrisk.at=seq(0,48, 12), atrisk.title="", atrisk.cex=0.8, lwd=1)
axis(side=1, pos=c(0,0), lwd=2, at=seq(0,48, 12))
axis(side=2, pos=0, lwd=2, at=seq(0,1, 0.25), labels=c("0", "25", "50", "75", "100"))
dev.off()
round(survdiff(Surv(tt_os_d, event_os)~med_p, tpt1_a)$pvalue,8) #p = 0.0080

d <- datadist(tpt1_a)
options(datadist = d)
cox_model <- cph(Surv(tt_os_d, event_os) ~ rcs(pyruvate, 4) + age10 + cohort + ecog, data = tpt1_a, x=TRUE, y=TRUE)
summary(cox_model)

cox_model <- coxph(Surv(tt_os_d, event_os) ~ log(total) + age10 + cohort + ecog, data = tpt1_a)
summary(cox_model)
cox_model <- coxph(Surv(tt_pfs_d, pfs_event) ~ log(total) + age10 + cohort + ecog, data = tpt1)
summary(cox_model)


cox_model <- coxph(Surv(tt_os_d, event_os) ~ log(pyruvate) + age10 + cohort + ecog, data = tpt1_a)
summary(cox_model)
cox_model <- coxph(Surv(tt_pfs_d, pfs_event) ~ log(pyruvate) + age10 + cohort + ecog, data = tpt1_a)
summary(cox_model)

#Figure 2B - PFS
pdf("figures/pfs.pdf", height=6, width=4)
plot(prodlim(Hist(tt_pfs_d/30.44, pfs_event)~med_p, tpt1_a), confint=T, marktime=T, xlim=c(0,48), axes=F, background.fg=NA,  atrisk.labelcol="black", atrisk.labels=levels(tpt1_a$med_p),
     background.border=NA, atrisk=T, xlab="Months", ylab="Progression-Free Survival, %", legend=T, legend.x="topright", legend.title="", legend.cex=0.6, atrisk.at=seq(0,48, 12), atrisk.title="", atrisk.cex=0.8, lwd=1)
axis(side=1, pos=c(0,0), lwd=2, at=seq(0,48, 12))
axis(side=2, pos=0, lwd=2, at=seq(0,1, 0.25), labels=c("0", "25", "50", "75", "100"))
dev.off()
round(survdiff(Surv(tt_pfs_d, pfs_event)~med_p, tpt1_a)$pvalue,8) #p = 0.070




ps <- serse_phylo %>% filter(sample %in% tpt1$WMS_SGPID)


sample_df <- data.frame(sample_data(ps))
merged_df <- left_join(sample_df, tpt1_a %>% select(id_int, `pyruvate`), by = "id_int")
merged_df <- left_join(merged_df, sdiv %>% rownames_to_column("WMS_SGPID") %>% left_join(ids %>% select(WMS_SGPID)), by=c("experiments" = "WMS_SGPID"))
merged_df$log_total <- log(merged_df$total)
merged_df$log_acoa <- log(merged_df$pyruvate)
new_sample_data <- sample_data(as(merged_df, "data.frame"))
rownames(new_sample_data) <- rownames(sample_data(ps))
sample_data(ps) <- new_sample_data



#species pcoa
pdf("figures/PCoA_species_acoa-25Oct2024.pdf", height=6, width=6)
ps %>%
  filter(SAMPID %in% tpt1$SAMPID) %>%
  tax_fix() %>%
  tax_transform(trans = "compositional", rank = "Species") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(fill = "log_acoa", size = 4, shape=21, alpha=0.6,
           plot_taxa = c(1:5)) +
  scale_fill_viridis(option="B", begin = 0, end=1) +
  theme_classic() +
  coord_cartesian(ylim=c(-7.5,9)) +
  scale_x_reverse()
dev.off()
plot_ord_sp <- ps %>%
  filter(SAMPID %in% tpt1$SAMPID) %>%
  tax_fix() %>%
  tax_transform(trans = "compositional", rank = "Species") %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(fill = "log_acoa", size = 4, shape=21, alpha=0.6,
           plot_taxa = c(1:5)) +
  scale_fill_viridis(option="B", begin = 0, end=1) +
  theme_classic()

summary(lm(plot_ord_sp$data$MDS1~eval(-plot_ord_sp$data$log_acoa)))



source("scripts/add_bacteroidales.R") #this script updates all Bacteroidales order to
###                                    list family as "Bacteroidales_order

pdf("figures/PCA_biplot_acoa.pdf", height=6, width=6)
ps2 %>%
  filter(SAMPID %in% tpt1$SAMPID) %>%
  tax_fix() %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ord_calc("PCA") %>%
  ord_plot(fill = "log_acoa", size = 4, shape=21, alpha=0.6,
           plot_taxa = c("Lachnospiraceae", "Bacteroidales_order", "Ruminococcaceae")) +
  scale_fill_viridis(option="B", begin = 0, end=1) +
  theme_classic()
dev.off()




compare <- ps %>% get.otu.melt()

compare_order <- compare %>%
  group_by(SAMPID, Order) %>%
  summarise(total_order = sum(pctseqs)) %>%
  left_join(tpt1_a %>% select(SAMPID, `pyruvate`))

compare_family <- compare %>%
  group_by(SAMPID, Family) %>%
  summarise(total_family = sum(pctseqs)) %>%
  left_join(tpt1_a %>% select(SAMPID, `pyruvate`))


#bacteroidaceae
compare_bac <- compare_order %>%
  filter(Order == "Bacteroidales")

pdf("figures/bacteroidetes.pdf", height=4, width=4)
ggplot(compare_bac) +
  geom_smooth(aes(x=pyruvate, y=total_order), method=lm, color="#68879b") +
  geom_point(aes(x=pyruvate, y=total_order),
             shape=21, size=3, fill="#68879b", alpha=0.8) +
  theme_classic()
dev.off()
cor.test(compare_bac$pyruvate, compare_bac$total_order)
#r2 = 0.2570, p = 5.6e-11


#ruminococcaceae
compare_rum <- compare_family %>%
  filter(Family == "Ruminococcaceae")

pdf("figures/ruminococcaceae.pdf", height=4, width=4)
ggplot(compare_rum) +
  geom_smooth(aes(x=pyruvate, y=total_family), method=lm, color="#e05f5f") +
  geom_point(aes(x=pyruvate, y=total_family),
             shape=21, size=3, fill="#e05f5f", alpha=0.8) +
  theme_classic()
dev.off()
cor.test(compare_rum$pyruvate, compare_rum$total_family)
#r2 = 0.227, p=1.04e-9



#lachnospiraceae
compare_lac <- compare_family %>%
  filter(Family == "Lachnospiraceae")
t
pdf("figures/lachnospiraceae.pdf", height=4, width=4)
ggplot(compare_lac) +
  geom_smooth(aes(x=pyruvate, y=total_family), method=lm, color="#e05f5f") +
  geom_point(aes(x=pyruvate, y=total_family),
             shape=21, size=3, fill="#e05f5f", alpha=0.8) +
  theme_classic()
dev.off()
cor.test(compare_lac$pyruvate, compare_lac$total_family)
#r2 = 0.177, p=1.1e-7


###### Supplemental figures


# Fit a Cox model with a restricted cubic spline
model_spline_os <- coxph(Surv(tt_os_d, event_os)~rcs(total,4)+age10+factor(cohort)+ecog, tpt1_a)
ptemp <- termplot(model_spline_os, se=T, plot=F)
buterm <- ptemp$total
center <- buterm$y[which(buterm$y == median(buterm$y))]
ytemp <- buterm$y + outer(buterm$se, c(0, -1.96, 1.96), '*')
exp_ytemp <- exp(ytemp - center)

# Create a spline_data frame for plotting
spline_data <- data.frame(buty = buterm$x, Estimate = exp_ytemp[,1],
                          Lower = exp_ytemp[,2], Upper = exp_ytemp[,3])

pdf("figures/supp_hr_spline_total.pdf", height=6, width=6)
ggplot(spline_data, aes(x = buty)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = Estimate), color = "blue") +
  geom_rug(sides = "b") +  # Add rug plot at the bottom ('b') of the plot
  scale_y_log10() +  # Log scale for y-axis
  labs(x = "Total butyrate-producing genes present", y = "Hazard ratio") +
  coord_cartesian(xlim=c(0,2000))+
  theme_minimal()
dev.off()


univ2 <- function(cox, k=1, l=1){
  coef <- exp(coef(cox))
  confint <- exp(confint(cox))
  p <- summary(cox)$coefficients[,5]
  out <- data.frame("HR"=round(coef,4), round(confint,4), "p"=round(p,9))
  q <- nrow(out)
  if(l== -1){out <- out[q,]
  }else{out <- out[k:l,]}
  return(out)

}



forest_os <- rbind(
univ2(coxph(Surv(tt_os_d, event_os)~log(total)+age10+cohort+ecog, tpt1_a)),
univ2(coxph(Surv(tt_os_d, event_os)~log(pyruvate)+age10+cohort+ecog, tpt1_a)),
univ2(coxph(Surv(tt_os_d, event_os)~log(lysine+0.0000001)+age10+cohort+ecog, tpt1_a)),
univ2(coxph(Surv(tt_os_d, event_os)~log(fouraminobutyrate+0.000001)+age10+cohort+ecog, tpt1_a)),
univ2(coxph(Surv(tt_os_d, event_os)~log_glutarate+age10+cohort+ecog, tpt1_a)))

pdf("forest_hr_pathway.pdf", height=2, width=4)
forest_os %>%
  rename_with(~ c("HR", "CI_low", "CI_high", "p"), .cols = everything()) %>%
  rownames_to_column("scfa") %>%
  ggplot() +
  geom_vline(xintercept=1, lty=2, size=0.25)+
  geom_point(aes(x=HR, y=fct_reorder(scfa, HR, .desc = T)))+
  geom_linerange(aes(xmin=CI_low, xmax=CI_high, y=fct_reorder(scfa, HR)), size=0.5) +
  scale_x_continuous(transform = "log", breaks=round(exp(seq(-2,1.5,0.5)),1)) +
  theme_classic() +
  coord_cartesian(xlim=exp(c(-2, 0.3)))
dev.off()

forest_pfs <- rbind(
  univ2(coxph(Surv(tt_pfs_d, pfs_event)~log(total)+age10+cohort+ecog, tpt1_a)),
  univ2(coxph(Surv(tt_pfs_d, pfs_event)~log(pyruvate)+age10+cohort+ecog, tpt1_a)),
  univ2(coxph(Surv(tt_pfs_d, pfs_event)~log(lysine+0.0000001)+age10+cohort+ecog, tpt1_a)),
  univ2(coxph(Surv(tt_pfs_d, pfs_event)~log(fouraminobutyrate+0.000001)+age10+cohort+ecog, tpt1_a)),
  univ2(coxph(Surv(tt_pfs_d, pfs_event)~log_glutarate+age10+cohort+ecog, tpt1_a)))

#pdf("forest_hr_pathway_pfs.pdf", height=2, width=4)
forest_pfs %>%
  rename_with(~ c("HR", "CI_low", "CI_high", "p"), .cols = everything()) %>%
  rownames_to_column("scfa") %>%
  ggplot() +
  geom_vline(xintercept=1, lty=2, size=0.25)+
  geom_point(aes(x=HR, y=fct_reorder(scfa, HR, .desc = T)))+
  geom_linerange(aes(xmin=CI_low, xmax=CI_high, y=fct_reorder(scfa, HR)), size=0.5) +
  scale_x_continuous(transform = "log", breaks=round(exp(seq(-2,1.5,0.5)),1)) +
  theme_classic() +
  coord_cartesian(xlim=exp(c(-2, 0.3)))
dev.off()



pdf("figures/but_pyr.pdf", height=4, width=4)
ggplot(scfa2) +
  geom_smooth(aes(x=log(pyruvate), y=Butyric.acid), method=lm, color="#e05f5f") +
  geom_point(aes(x=log(pyruvate), y=Butyric.acid),
             shape=21, size=3, fill="#e05f5f", alpha=0.8) +
  theme_classic()
cor.test(scfa2$Butyric.acid, log(scfa2$pyruvate))
#R2 = 0.0823, p = 0.022
dev.off()

pdf("figures/but_tot.pdf", height=4, width=4)
ggplot(scfa2) +
  geom_smooth(aes(x=log(total), y=Butyric.acid), method=lm, color="#e05f5f") +
  geom_point(aes(x=log(total), y=Butyric.acid),
             shape=21, size=3, fill="#e05f5f", alpha=0.8) +
  theme_classic()
dev.off()
cor.test(scfa2$Butyric.acid, log(scfa2$total))
#R2 = 0.035, p = 0.138



#species level
ps %>%
  filter(SAMPID %in% tpt1$SAMPID) %>%
  tax_fix() %>%
  tax_transform(trans = "compositional", rank = "Species") %>%
  dist_calc("bray") %>%
  ord_calc("PCoA") %>%
  ord_plot(fill = "Shannon", size = 4, shape=21, alpha=0.6) +
  scale_fill_viridis(option="B", begin = 0, end=1) +
  theme_classic()





# Wald test for spline terms
fit <- coxph(Surv(tt_os_d, event_os) ~ rcs(pyruvate, 4) + age10 + factor(cohort) + ecog, data = tpt1_a)
cov_matrix <- vcov(fit)
coef_spline <- coef(fit)[grep("rcs\\(pyruvate", names(coef(fit)))]
cov_spline <- cov_matrix[grep("rcs\\(pyruvate", rownames(cov_matrix)), grep("rcs\\(pyruvate", colnames(cov_matrix))]
wald_stat <- t(coef_spline) %*% solve(cov_spline) %*% coef_spline
p_value <- 1 - pchisq(wald_stat, length(coef_spline))
wald_stat
p_value
