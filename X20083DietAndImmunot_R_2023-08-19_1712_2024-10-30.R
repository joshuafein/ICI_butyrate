rm(list=ls())
graphics.off()

library(Hmisc)

data=read.csv('DATA_dietary_analysis.csv', na.strings = c("", "NA"))


library(dplyr)
library(gtsummary)
library(survminer)
library(forestplot)


label(data$age)="Age at time of dietary questionnaire (in years)"
label(data$sex)="Sex"
label(data$race)="Race"
label(data$latinx)="Ethnicity"
label(data$cancer_type)="Malignancy type"
label(data$varhistoyn_uc)="Variant histology? (i.e. anything other than urothelial carcinoma, NOS)"
label(data$visceral_mets)="Visceral metastases present? (per Bajorin risk factors: lung liver, bone or other non-lymph node metastasis)"
label(data$liver_mets)="Liver metastases present?"
label(data$non_cc_rcc_yn)="Non-Clear Cell Histology?"
label(data$imdc)="IMDC (International Metastatic Renal Cell Carcinoma Database Consortium) group?  (IMDC risk factors: Karnofsky performance status score < 80 Time from original diagnosis to initiation of targeted therapy < 1 year Hemoglobin less than the lower limit of normal Serum calcium greater than the upper limit of normal Neutrophil count greater than the upper limit of normal Platelet count greater than the upper limit of normal)"
label(data$io_regimen)="Planned immunotherapy regimen"
label(data$bmi)="BMI at time of questionnaire in kg/m2"
label(data$ecogps)="ECOG performance status at time of questionnaire"
label(data$pfs_mo)="Progression-free survival time interval in months"
label(data$pod_status)="Disease progression status"
label(data$os_mo)="Overall survival time interval in months"
label(data$os_status)="Overall survival status"

data$sex.factor = factor(data$sex,levels=c("1","2","3"))
data$race.factor = factor(data$race,levels=c("1","2","3","4","99","5"))
data$latinx.factor = factor(data$latinx,levels=c("1","2","3"))
data$cancer_type.factor = factor(data$cancer_type,levels=c("1","2","3","5","6","4"))
data$varhistoyn_uc.factor = factor(data$varhistoyn_uc,levels=c("1","0"))
data$visceral_mets.factor = factor(data$visceral_mets,levels=c("1","2","3"))
data$liver_mets.factor = factor(data$liver_mets,levels=c("0","1","99"))
data$non_cc_rcc_yn.factor = factor(data$non_cc_rcc_yn,levels=c("1","0"))
data$imdc.factor = factor(data$imdc,levels=c("0","1","2","3"))
data$io_regimen.factor = factor(data$io_regimen,levels=c("1","2","3","4","14","5","6","7","8","18","9","17","10","11","16","12","13","15","99"))
data$ecogps.factor = factor(data$ecogps,levels=c("0","1","2","3","4","5","6"))
data$kps.factor = factor(data$kps,levels=c("0","1","2","3","4","5","6","7","8","9","10","99"))
data$pod_status.factor = factor(data$pod_status,levels=c("0","1"))
data$os_status.factor = factor(data$os_status,levels=c("0","1"))

levels(data$sex.factor)=c("Male","Female","Other")
levels(data$race.factor)=c("White","Black","Asian","Native American","Other","Prefers not to disclose/Unknown")
levels(data$latinx.factor)=c("Latinx","Non-Latinx","Prefers not to disclose/Unknown")
levels(data$cancer_type.factor)=c("urothelial cancer of the bladder","urothelial cancer of the ureter","urothelial cancer of the renal pelvis","urothelial cancer of the ureter and renal pelvis","urothelial cancer of the urethra","renal cell carcinoma")
levels(data$varhistoyn_uc.factor)=c("Yes","No")
levels(data$visceral_mets.factor)=c("Yes","No","Unknown")
levels(data$liver_mets.factor)=c("Yes","No","Unknown")
levels(data$non_cc_rcc_yn.factor)=c("Yes","No")
levels(data$imdc.factor)=c("Favorable risk: None of the above risk factors present.","Intermediate risk: 1 or 2 of the above risk factors present.","Poor risk: 3 or more risk factors present.","Unknown")
levels(data$io_regimen.factor)=c("pembrolizumab 1st line metastatic","atezolizumab 1st line metastatic","pembrolizumab metastatic, 2nd line+","atezolizumab, metastatic, 2nd line+","nivolumab, metastatic, 2nd line+","maintenance avelumab","adjuvant nivolumab","adjuvant pembrolizumab","pembrolizumab + axitinib, metastatic, 1st line","pembrolizumab + axitinib, metastatic, 2nd line+","ipilimumab + nivolumab, metastatic, 1st line","ipilimumab + nivolumab, metastatic, 2nd line+","cabozantinib + nivolumab, metastatic, 1st line","avelumab + axitinib, metastatic, 1st line","pembrolizumab + lenvatinib, metastatic, 1st line","pembrolizumab + lenvatinib, metastatic, 2nd line+","pembrolizumab, NMIBC","investigational protocol regimen","unknown")
levels(data$ecogps.factor)=c("0","1","2","3","4","5","Unknown")
levels(data$kps.factor)=c("100","90","80","70","60","50","40","30","20","10","0","Unknown")
levels(data$pod_status.factor)=c("Censored","Progressed")
levels(data$os_status.factor)=c("Censored","Deceased")

data <- mutate(data,io_line = case_when(io_regimen==1~'PD-(L)1 1st line',
                                                io_regimen==2~'PD-(L)1 1st line',
                                        io_regimen==5~'Maintenance avelumab',
                                        io_regimen==3~'PD-(L)1 2nd line+',
                                        io_regimen==4~'PD-(L)1 2nd line+',
                                        io_regimen==14~'PD-(L)1 2nd line+',
                                        io_regimen==6~'Adjuvant PD-1',
                                        io_regimen==7~'Adjuvant PD-1',
                                        io_regimen==13~'Pembrolizumab for NMIBC',
                                        io_regimen==17~'Ipilimumab plus nivolumab 2nd line+',
                                        io_regimen==9~'1st line ipilimumab plus nivolumab',
                                        io_regimen==8~'1st line IO + TKI',
                                        io_regimen==10~'1st line IO + TKI',
                                        io_regimen==11~'1st line IO + TKI',
                                        io_regimen==16~'1st line IO + TKI',
                                        io_regimen==18~'IO + TKI 2nd line+',
                                        io_regimen==12~'IO + TKI 2nd line+',
                                        io_regimen==15~'Investigational regimen'))

data$io_line.factor <- as.factor(data$io_line)

data <- mutate(data,ecog_kpsConverted = case_when(ecogps==0~"0",
                                                  ecogps==1~"1",
                                                  ecogps==2~"2",
                                                  ecogps==3~"3",
                                                  ecogps==4~"4",
                                                  ecogps==5~"5",
                                                  ecogps==6&kps==0~"0",
                                                  ecogps==6&kps==1~"0",
                                                  ecogps==6&kps==2~"1",
                                                  ecogps==6&kps==3~"1",
                                                  ecogps==6&kps==4~"2",
                                                  ecogps==6&kps==5~"2",
                                                  ecogps==6&kps==6~"3",
                                                  ecogps==6&kps==7~"3",
                                                  ecogps==6&kps==8~"4",
                                                  ecogps==6&kps==9~"4",
                                                  ecogps==6&kps==10~"5",
                                                  ecogps==6&kps==99~"99"))

data$ecog_kpsConverted[data$ecog_kpsConverted==99] <- NA

data<-mutate(data,ecog_1to5 = case_when(ecog_kpsConverted<1~0,
                                        ecog_kpsConverted>0~1))

data<-mutate(data,ecog_1to5_final = case_when(ecog_1to5=="1"~1,
                                                ecog_1to5=="0"~0,
                                                is.na(ecog_1to5)==TRUE~1))

data[is.na(data$liver_mets.factor),"liver_mets.factor"]<-"No"

data$age.numeric <-as.numeric(data$age)
data <- mutate(data,ageby10 = age.numeric/10)

data <- mutate(data,cancer_RCCvsUC = ifelse(cancer_type==4,0,1))
data <- mutate(data,cancer_RCCvsUC_NAME = ifelse(cancer_type==4,"RCC","UC"))

data<-mutate(data,Impact.TMB.Imputed = case_when(is.na(Impact.TMB.Score)==FALSE~Impact.TMB.Score,
                                                                       is.na(Impact.TMB.Score)==TRUE&cancer_RCCvsUC==1~5.8,
                                                                       is.na(Impact.TMB.Score)==TRUE&cancer_RCCvsUC==0~3)) ##imputed median TMB of 5.8 for missing values for UC based on the 2017 bladder TCGA publication; imputed TMB =3 for RCC

mUCplusmRCdata=filter(data, io_regimen!=6 & io_regimen!=13 & io_regimen!=15 & io_regimen!=7 & io_regimen!=99)
mUCplusmRCdataFiber=filter(mUCplusmRCdata,((calor>=600 & calor<=4200 & sex.factor=="Male") 
                                           |(calor>=500 & calor<=3500 & sex.factor=="Female"))&(nblank<=70))


data <- mutate(data,  fiber_dichotomized = case_when(aofib>=15~1,
                                                  aofib<15~0))

data$fiber_dichotomized_factor = factor(data$fiber_dichotomized,levels=c("0","1"))
levels(data$fiber_dichotomized_factor)=c("Less than 15 g/d","15+ g/day")

label(data$ecog_kpsConverted)="ECOG performance status"
label(data$latinx.factor)="Ethnicity"
label(data$race.factor)="Race"
label(data$visceral_mets.factor)="Visceral metastases"
label(data$liver_mets.factor)="Liver metastases"
label(data$io_line)="Line of therapy"
label(data$sex.factor)="Sex"
label(data$age)="Age (in years)"
label(data$varhistoyn_uc)="Variant histology present"
label(data$io_regimen.factor)="Immunotherapy regimen"
label(data$aofib)="Fiber intake (g/day)"
label(data$bmi) = "BMI (kg/m2)"
label(data$imdc.factor) = "IMDC Risk Category"
label(data$non_cc_rcc_yn.factor) = "Non-clear cell histology"


UCdata=filter(data,cancer_type!=4)


label(data$Impact.TMB.Score)="TMB (mut/Mb)"
label(data$Impact.TMB.Imputed)="TMB (mut/Mb)"
label(UCdata$Impact.TMB.Score)="TMB (mut/Mb)"

UCdata$TMBhigherThan10<-ifelse(UCdata$Impact.TMB.Score>=10,1,0)

UCdata<-mutate(UCdata,TMBhigherThan10.Imputed = case_when(is.na(TMBhigherThan10)==FALSE~TMBhigherThan10,
                                                          is.na(TMBhigherThan10)==TRUE~0))


mUCdata=filter(UCdata,io_regimen!=6 & io_regimen!=13 & io_regimen!=15)


dataUCavelumab=filter(mUCdata,io_regimen==5)
dataUC_NOTavelumab=filter(mUCdata,io_regimen!=5)


data$calor<-as.numeric(data$calor)
sum(data$nblank > 70, na.rm = TRUE) 
sum(data$calor<600 & data$calor>4200 & data$sex.factor=="Male", na.rm = TRUE) 
sum(data$calor<500 & data$calor>3500 & data$sex.factor=="Female", na.rm = TRUE) 

mUCdataFiber=filter(mUCdata,((calor>=600 & calor<=4200 & sex.factor=="Male") 
                    |(calor>=500 & calor<=3500 & sex.factor=="Female"))&(nblank<=70))
dataUCavelumabFiber=filter(dataUCavelumab,((calor>=600 & calor<=4200 & sex.factor=="Male") 
                                 |(calor>=500 & calor<=3500 & sex.factor=="Female"))&(nblank<=70))
dataUC_NOTavelumabFiber=filter(dataUC_NOTavelumab,((calor>=600 & calor<=4200 & sex.factor=="Male") 
                                           |(calor>=500 & calor<=3500 & sex.factor=="Female"))&(nblank<=70))


RCdata=filter(data,cancer_type==4)

mRCdata=filter(RCdata,io_regimen!=7 & io_regimen!=15 & io_regimen!=99)
dataRCIpiNivo=filter(RCdata,io_regimen==9 | io_regimen==17)
dataRCioTKI=filter(RCdata,io_regimen==8 | io_regimen==18 | io_regimen==10 | io_regimen==11 | io_regimen==16 | io_regimen==12)


mRCdataFiber=filter(mRCdata,((calor>=600 & calor<=4200 & sex.factor=="Male") 
                             |(calor>=500 & calor<=3500 & sex.factor=="Female"))&(nblank<=70))
dataRCIpiNivoFiber=filter(dataRCIpiNivo,((calor>=600 & calor<=4200 & sex.factor=="Male") 
                                         |(calor>=500 & calor<=3500 & sex.factor=="Female"))&(nblank<=70))
dataRCioTKIFiber=filter(dataRCioTKI,((calor>=600 & calor<=4200 & sex.factor=="Male") 
                                     |(calor>=500 & calor<=3500 & sex.factor=="Female"))&(nblank<=70))

AllMetdata=filter(data,io_regimen!=6 & io_regimen!=13 & io_regimen!=7 & io_regimen!=15)


AllMetdataFiber=filter(AllMetdata,((calor>=600 & calor<=4200 & sex.factor=="Male") 
                                   |(calor>=500 & calor<=3500 & sex.factor=="Female"))&(nblank<=70))


tbl_summary(mUCdataFiber%>%select(age,sex.factor,race.factor,latinx.factor,
                             bmi,ecog_kpsConverted,varhistoyn_uc,
                             visceral_mets.factor,liver_mets.factor,
                             Impact.TMB.Score,
                             io_line)) %>% modify_caption("**Supplementary Table. Baseline characteristics of patients with mUC**")


mRCdata2 <- mRCdataFiber %>% select(age,sex.factor,race.factor,latinx.factor,
                               bmi,
                               ecog_kpsConverted,
                               non_cc_rcc_yn.factor,
                               imdc.factor,
                               io_line)

mRCdata2 %>% tbl_summary()


# Figure 1

  makeSpline <- function(df, survtime, survevent, var,
                         knots=knots, centerpoint = NULL,
                         title=paste("Spline with", knots, "knots"),
                         xlab = xlab, ylab=ylab) {
    
    require(survival)
    require(tidyverse)
    
    model_spline <- coxph(Surv(survtime, survevent)~pspline(var,knots), df) #
    ptemp <- termplot(model_spline, se=T, plot=F)
    
    splinem <- ptemp$var
    
    if(is.null(centerpoint)){centerpoint = median(var, na.rm=T)}
   
    if(centerpoint %in% splinem$x){
      c_row = which(splinem$x == centerpoint)
    } else {
      c_row = which(abs(splinem$x - centerpoint) == min(abs(splinem$x - centerpoint)))
    }
    
    center <- splinem$y[c_row]
    
    ytemp <- splinem$y + outer(splinem$se, c(0, -1.96, 1.96), '*')
    
    exp_ytemp <- exp(ytemp - center)
    
    # Create a spline_data frame for plotting
    spline_data <- data.frame(xvar = splinem$x, Estimate = exp_ytemp[,1],
                              Lower = exp_ytemp[,2], Upper = exp_ytemp[,3])
    
    ggplot(spline_data, aes(x = xvar)) +
      geom_hline(yintercept=1, lty=2)+
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.5) +
      geom_line(aes(y = Estimate), color = "blue") +
      geom_rug(sides = "b") +  # Add rug plot at the bottom ('b') of the plot
      scale_y_log10() +  # Log scale for y-axis
      scale_x_log10() +
      labs(x = xlab, y = ylab, title=title) +
      theme_classic()
  }
  


#Figure 1A
makeSpline(AllMetdataFiber, AllMetdataFiber$pfs_mo, AllMetdataFiber$pod_status, 
           AllMetdataFiber$aofib,
           xlab = "Daily average fiber intake (g/day)", ylab="Hazard ratio (95% CI)",
           title="Smoothing spline for progression-free survival",
           knots=2,
           centerpoint = 15)

#Figure 1B
makeSpline(AllMetdataFiber, AllMetdataFiber$os_mo, AllMetdataFiber$os_status, 
           AllMetdataFiber$aofib,
           xlab = "Daily average fiber intake (g/day)", ylab="Hazard ratio (95% CI)",
           title="Smoothing spline for overall survival",
           knots=2, 
           centerpoint = 15)


#Figure 1C
kmPFS<-survfit(Surv(pfs_mo, pod_status)~fiber_dichotomized,
               type="kaplan-meier", data=AllMetdataFiber)

ggsurvplot(kmPFS,risk.table=TRUE,
           xlab="Time (months)", 
           ylab="Progression-free survival",
                      pval = TRUE,
           legend=c(0.8,0.9),
           title="Patients with urothelial or renal cell cancer on ICB",
           legend.title="Fiber intake",
           legend.labs=c("<15 g/day","15+ g/day"),
           #surv.median.line = "hv",
                      conf.int=TRUE,
           break.x.by=6,
           axes.offset=FALSE,
           xlim=c(0,25),
           data = AllMetdataFiber)

#Figure 1D
kmOS<-survfit(Surv(os_mo, os_status)~fiber_dichotomized,
               type="kaplan-meier", data=AllMetdataFiber)

ggsurvplot(kmOS,risk.table=TRUE,
           xlab="Time (months)", 
           ylab="Overall survival",
           pval = TRUE,
           legend="none",
           title="Patients with urothelial or renal cell cancer on ICB",
           legend.title="Fiber intake",
           legend.labs=c("<15 g/day","15+ g/day"),
           #surv.median.line = "hv",
           conf.int=TRUE,
           break.x.by=6,
           axes.offset=FALSE,
           xlim=c(0,25),
           data = AllMetdataFiber)

# Figure 1E
cox_pfs_all_uv <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized,
                    AllMetdataFiber)
summary(cox_pfs_all_uv)

cox_pfs_RC <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized,
                   mRCdataFiber)
summary(cox_pfs_RC)

cox_pfs_RC_IpiNivo <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized,
                           dataRCIpiNivoFiber)
summary(cox_pfs_RC_IpiNivo)

cox_pfs_RC_TKI_IO <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized,
                          dataRCioTKIFiber)
summary(cox_pfs_RC_TKI_IO)

cox_pfs_UC <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized,
                           mUCdataFiber)
summary(cox_pfs_UC)

cox_pfs_UC_NOTavelumab <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized,
                   dataUC_NOTavelumab)
summary(cox_pfs_UC_NOTavelumab)

cox_pfs_UC_avelumab <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized,
                            dataUCavelumab)
summary(cox_pfs_UC_avelumab)


            
forestplot_data <- tibble::tibble(mean  = c(0.4409,0.3485, 0.4555, 0.3034, 0.6998, 0.8952, 0.8994),
                            lower = c(0.2582,0.1608, 0.1071, 0.09031, 0.3332, 0.3767, 0.1789),
                            upper = c(0.7528,0.7551, 1.937, 1.02, 1.47, 2.128, 4.521),
                            study = c("All","RCC", "RCC ipi-nivo", "RCC TKI-ICB",
                                      "UC", "UC non-maintenance ICB", "UC maintenance ICB"),
                            n = c("88","48","15","24","40","26","14"),
                            events = c("57","28","8","12","29","21","8"),
                            HR = c("0.44","0.35", "0.46", "0.30", "0.70", "0.90", "0.90"))

forestplot_data |>
  forestplot(labeltext = c(study, n, events, HR),
             clip = c(0.09, 2.5),
             xlog = TRUE,
             xlab = "Hazard ratio for PFS (95% CI)",
             xticks = c(log(0.1),log(0.25),log(0.5),log(1.0),log(2.5)),
             title = "Univariable progression-free survival subgroup analyis") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_header(study = c("Subgroup"),
                n = c("n"),
                events = c("Events"),
                HR = c("HR")) 


# Figure 1F
cox_pfs_all <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized + ageby10 + ecog_1to5_final,
                   AllMetdataFiber)
summary(cox_pfs_all)

cox_pfs_RC <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized + ageby10 + ecog_1to5_final,
                   mRCdataFiber)
summary(cox_pfs_RC)

cox_pfs_RC_IpiNivo <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized+ ageby10 + ecog_1to5_final,
                           dataRCIpiNivoFiber)
summary(cox_pfs_RC_IpiNivo)

cox_pfs_RC_TKI_IO <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized+ ageby10 + ecog_1to5_final,
                          dataRCioTKIFiber)
summary(cox_pfs_RC_TKI_IO)

cox_pfs_UC <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized+ ageby10 + ecog_1to5_final,
                   mUCdataFiber)
summary(cox_pfs_UC)

cox_pfs_UC_NOTavelumab <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized+ ageby10 + ecog_1to5_final,
                               dataUC_NOTavelumabFiber)
summary(cox_pfs_UC_NOTavelumab)

cox_pfs_UC_avelumab <-coxph(Surv(pfs_mo,pod_status)~ fiber_dichotomized+ ageby10 + ecog_1to5_final,
                            dataUCavelumabFiber)
summary(cox_pfs_UC_avelumab)

forestplot_data <- tibble::tibble(mean  = c(0.405, 0.3185, 0.6549),
                                  lower = c(0.2298, 0.1322, 0.3038),
                                  upper = c(0.7138, 0.7676, 1.412),
                                  study = c("All", "RCC", "UC"),
                                  n = c("88","48","40"),
                                  events = c("57","28","29"),
                                  HR = c("0.41","0.32","0.65"))

forestplot_data |>
  forestplot(labeltext = c(study, n, events, HR),
             clip = c(0.09, 2.5),
             xlog = TRUE,
             xlab = "Hazard ratio for PFS (95% CI)",
             xticks = c(log(0.1),log(0.25),log(0.5),log(1.0),log(2.5)),
             title = "Multivariable progression-free survival subgroup analyis") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_header(study = c("Subgroup"),
                n = c("n"),
                events = c("Events"),
                HR = c("HR"))

#Additional analyses
cox_pfs_allMet_mv <-coxph(Surv(pfs_mo,
                               pod_status)~ fiber_dichotomized + 
                            ageby10 + 
                            ecog_1to5_final +
                            cancer_RCCvsUC_NAME,
                          AllMetdataFiber)

summary(cox_pfs_allMet_mv)


cox_pfs_allMet_mv2 <-coxph(Surv(pfs_mo,
                                pod_status)~ fiber_dichotomized + 
                             calor +
                             ageby10 + 
                             ecog_1to5_final +
                             cancer_RCCvsUC_NAME,
                           AllMetdataFiber)

summary(cox_pfs_allMet_mv2)


cox_pfs_allMet_mv3 <-coxph(Surv(pfs_mo,
                                pod_status)~ fiber_dichotomized + 
                             Impact.TMB.Imputed +
                             ageby10 + 
                             ecog_1to5_final +
                             cancer_RCCvsUC_NAME,
                           AllMetdataFiber)

summary(cox_pfs_allMet_mv3)


#Extended Data Figure 2
kmPFS<-survfit(Surv(pfs_mo, pod_status)~fiber_dichotomized,
               type="kaplan-meier", data=mUCdataFiber)

ggsurvplot(kmPFS,risk.table=TRUE,
           xlab="Time (months)", 
           ylab="Progression-free survival",
           #pval = TRUE,
           legend="none",
           title="Patients with urothelial cancer on ICB",
           legend.title="Fiber intake",
           legend.labs=c("<15 g/day","15+ g/day"),
           #surv.median.line = "hv",
           conf.int=TRUE,
           break.x.by=6,
           axes.offset=FALSE,
           xlim=c(0,16),
           data = mUCdataFiber)

kmPFS<-survfit(Surv(pfs_mo, pod_status)~fiber_dichotomized,
               type="kaplan-meier", data=mRCdataFiber)

ggsurvplot(kmPFS,risk.table=TRUE,
           xlab="Time (months)", 
           ylab="Progression-free survival",
           #pval = TRUE,
           legend="none",
           title="Patients with renal cell cancer on ICB",
           legend.title="Fiber intake",
           legend.labs=c("<15 g/day","15+ g/day"),
           #surv.median.line = "hv",
           conf.int=TRUE,
           break.x.by=6,
           axes.offset=FALSE,
           #xlim=c(0,16),
           data = mRCdataFiber)
