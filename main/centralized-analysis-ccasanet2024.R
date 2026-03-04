#---------------------------------------------------------------------------------#
#      Functions defined for centralized survival analysis of ccasanet data       #
#---------------------------------------------------------------------------------#

# for each primary analysis, there are several versions
## include haiti or exclude haiti
## use .95 or .90 quantile of distribution of last observations for defining closing dates
## with imputation or without imputation

library(survival)
library(rms)
library(mice)
library(dplyr)
library(stringr)
library(rlang)

# data, load the cleaned version of hiv

# major treatment change #
## identify relevant date
regimenchange.date = function(data, quantile = .95){
  n = nrow(data)
  if(quantile == .95){
    tmp = sapply(1:n, function(x){
      if(!is.na(data$ARTchange_status.95[x])){
        if(data$ARTchange_status.95[x]==2){
          data$ARTchange_y_d[x]
        }
        else if(data$ARTchange_status.95[x]==1){
          data$death_y_d[x]
        }
        else{
          min(data$lastobs_d[x], data$siteclosedate.95[x])
        }
      }
      else{
        NA
      }
    }) %>% as.Date()
  }
  else{ # .90 quantile
    tmp = sapply(1:n, function(x){
      if(!is.na(data$ARTchange_status.90[x])){
        if(data$ARTchange_status.90[x]==2){
          data$ARTchange_y_d[x]
        }
        else if(data$ARTchange_status.90[x]==1){
          data$death_y_d[x]
        }
        else{
          min(data$lastobs_d[x], data$siteclosedate.90[x])
        }
      }
      else{
        NA
      }
    }) %>% as.Date()    
  }
  return(tmp)
}

## unadjusted analysis of major treatment change, plot of cumulative incidence
regimenchange.cuminc = function(data, quantile=.95, includeHaiti=FALSE){
  # remove initial ART regimen of "Other" type
  data = data %>% subset(initialARTtype !="Other")
  data$initialARTtype = factor(data$initialARTtype)
  #
  data$ARTchange_status.date = regimenchange.date(data = data, quantile = quantile)
  data$ARTchange_status.days = as.integer(data$ARTchange_status.date) - as.integer(data$recart_d)
  # some have censored dates prior to the date of receiving initial ART due to the definition of closing dates
  data = data %>% subset(ARTchange_status.days>=0)
  # define time on study in years
  data$ARTchange_status.years = data$ARTchange_status.days/365.25
  if(!includeHaiti){
    data = data %>% subset(site!="haiti")
    data$site = factor(data$site)
    colorset = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")[-4] 
  }
  else{ # include Haiti
    colorset = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")
  }
  if(quantile==.95){
    data$ARTchange_status = factor(data$ARTchange_status.95, levels = c(0,1,2), labels = c("right censored", "death", "major treatment change"))
  }
  else{
    data$ARTchange_status = factor(data$ARTchange_status.90, levels = c(0,1,2), labels = c("right censored", "death", "major treatment change"))
  }
  cic_trt = npsurv(Surv(time = ARTchange_status.years, event = ARTchange_status)~site, data = data)
  cic_trt_overall = npsurv(Surv(time = ARTchange_status.years, event = ARTchange_status)~1, data = data)
  # plot
  survplot(cic_trt_overall, state = "major treatment change", conf = "none",
           label.curves = F,
           lty = "dashed",
           col = "black",
           xlab = "Time from ART initiation (years)",
           ylab = "Cumulative incidence of major treatment change",
           n.risk = T,
           y.n.risk = -0.2, adj.n.risk = 1, sep.n.risk = 2,
           ylim = c(0,1), xlim = c(0,9),
           cex.xlab = 0.7, cex.ylab = 0.7)
  survplot(cic_trt, state = "major treatment change", conf = "none",
           label.curves = F,
           lty = "solid",
           col = colorset,
           xlab = "Time from ART initiation (years)",
           ylab = "Cumulative incidence of major treatment change", add = T)
  legend("topleft",
         col = c("black", colorset),
         lty = c("dashed", rep("solid", length(colorset))),
         legend = c("overall", levels(data$site)),
         bty = "n",
         cex = 0.7
  )
}

## adjusted analysis for major treatment change, with imputation or without imputation
regimenchange.cscox = function(data, quantile=.95, complete=FALSE, includeHaiti=FALSE, includeHonduras=FALSE, merge=FALSE){
  # remove initial ART regimen of "Other" type
  data = data %>% subset(initialARTtype !="Other")
  data$initialARTtype = factor(data$initialARTtype)
  #
  data$ARTchange_status.date = regimenchange.date(data = data, quantile = quantile)
  data$ARTchange_status.days = as.integer(data$ARTchange_status.date) - as.integer(data$recart_d)
  # some have censored dates prior to the date of receiving initial ART due to the definition of closing dates
  data = data %>% subset(ARTchange_status.days>=0)
  # define time on study in years
  data$ARTchange_status.years = data$ARTchange_status.days/365.25
  if(quantile==.95){
    data$ARTchange_status = factor(data$ARTchange_status.95, levels = c(0,1,2), labels = c("right censored", "death", "major treatment change"))
  }
  else{
    data$ARTchange_status = factor(data$ARTchange_status.90, levels = c(0,1,2), labels = c("right censored", "death", "major treatment change"))
  }
  if(!includeHaiti){
    data = data %>% subset(site!="haiti")
  }
  if(!includeHonduras){
    data = data %>% subset(site!="honduras")
  }
  data$site = factor(data$site) 
  if(complete){ # without imputation
    if(!merge){
      data = data %>% subset(!is.na(sex)  & !is.na(mode_new) & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log) & !is.na(ARTchange_status))     
      fit = coxph(Surv(ARTchange_status.years, ARTchange_status=="major treatment change")~strata(site)+sex+age+mode_new+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+calendaryear+initialARTtype, data = data)      
    }
    else{
      data = data %>% subset(!is.na(sex_mode)  & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log) & !is.na(ARTchange_status))     
      fit = coxph(Surv(ARTchange_status.years, ARTchange_status=="major treatment change")~strata(site)+sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+calendaryear+initialARTtype, data = data)      
    }
  }
  else{ # with imputation
    var_model = c("ARTchange_status.years", "ARTchange_status", "site", "sex", "age", "mode_new", "clinicalAIDS", "cd4_v_initial_sqrt", "rna_v_initial_log", "calendaryear", "initialARTtype")
    imputed_data = mice(data %>% select(all_of(var_model)) , m = 5, method = "pmm", seed = 1117, printFlag = F)    
    if(!merge){
      # fit
      fit <- with(imputed_data,
                  coxph(Surv(ARTchange_status.years, ARTchange_status=="major treatment change") ~ strata(site)+sex+age+mode_new+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+calendaryear+initialARTtype))      
    }
    else{
      # derive sex_mode
      imputed_data_d = complete(imputed_data, action = "long", include = TRUE)
      imputed_data_d$sex_mode = sapply(1:nrow(imputed_data_d), function(i) {
        ifelse(imputed_data_d$sex[i]=="Female", "Female",
               ifelse(imputed_data_d$mode_new=="MSM", "Male-MSM", "Male-Other"))
      }) %>% factor(levels = c("Female", "Male-MSM", "Male-Other"))
      imputed_data_d = as.mids(imputed_data_d)
      # fit
      fit <- with(imputed_data_d,
                  coxph(Surv(ARTchange_status.years, ARTchange_status=="major treatment change") ~ strata(site)+sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+calendaryear+initialARTtype))       
    }
    fit = summary(pool(fit))
  }
  return(list(trt=data, fit=fit))
} 

# virologic failure #
## identify relevant date
virologicfailure.date = function(data, quantile=0.95){
  n = nrow(data)
  if(quantile==0.95){
    tmp = sapply(1:n, function(x){
      if(!is.na(data$vlfailure_status.95[x])){
        if(data$vlfailure_status.95[x]==2){
          data$vlfailure_y_d[x]
        }
        else if(data$vlfailure_status.95[x]==1){
          data$vlfailure_y_d[x]
        }
        else{
          min(data$vlfailure_y_d[x], data$siteclosedate.95[x])
        }
      }
      else{
        NA
      }
    }) %>% as.Date()    
  }
  else{
    tmp = sapply(1:nrow(data), function(x){
      if(!is.na(data$vlfailure_status.90[x])){
        if(data$vlfailure_status.90[x]==2){
          data$vlfailure_y_d[x]
        }
        else if(data$vlfailure_status.90[x]==1){
          data$vlfailure_y_d[x]
        }
        else{
          min(data$vlfailure_y_d[x], data$siteclosedate.90[x])
        }
      }
      else{
        NA
      }
    }) %>% as.Date()
  }
  return(tmp)
}

## unadjusted analysis of virologic failure, plot of cumulative incidence 
virologicfailure.cuminc = function(data, quantile=.95, includeHaiti=FALSE){
  # remove initial ART regimen of "Other" type
  data = data %>% subset(initialARTtype !="Other")
  data$initialARTtype = factor(data$initialARTtype)
  #
  data$vlfailure_status.date = virologicfailure.date(data = data, quantile = quantile)
  data$vlfailure_status.days = as.integer(data$vlfailure_status.date) - as.integer(data$recart_d)
  # some have censored dates prior to the date of receiving initial ART due to the definition of closing dates
  data = data %>% subset(vlfailure_status.days>=0)
  # define time on study in years
  data$vlfailure_status.years = data$vlfailure_status.days/365.25
  if(!includeHaiti){
    data = data %>% subset(site!="haiti")
    data$site = factor(data$site)
    colorset = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")[-4] 
  }
  else{ # include Haiti
    colorset = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")
  }
  if(quantile==.95){
    data$vlfailure_status = factor(data$vlfailure_status.95, levels = c(0,1,2), labels = c("right censored", "death", "virologic failure"))
  }
  else{
    data$vlfailure_status = factor(data$vlfailure_status.90, levels = c(0,1,2), labels = c("right censored", "death", "virologic failure"))
  }
  cic_vl = npsurv(Surv(time = vlfailure_status.years, event = vlfailure_status)~site, data = data)
  cic_vl_overall = npsurv(Surv(time = vlfailure_status.years, event = vlfailure_status)~1, data = data)
  # plot
  survplot(cic_vl_overall, state = "virologic failure", conf = "none",
           label.curves = F,
           lty = "dashed",
           col = "black",
           xlab = "Time from ART initiation (years)",
           ylab = "Cumulative incidence of virologic failure",
           n.risk = T,
           y.n.risk = -0.2, adj.n.risk = 1, sep.n.risk = 2,
           ylim = c(0,ifelse(length(levels(data$site))==6, 0.8, 0.4)), xlim = c(0,9),
           cex.xlab = 0.7, cex.ylab = 0.7)
  survplot(cic_vl, state = "virologic failure", conf = "none",
           label.curves = F,
           lty = "solid",
           col = colorset,
           xlab = "Time from ART initiation (years)",
           ylab = "Cumulative incidence of virologic failure", add = T)
  legend("topleft",
         col = c("black", colorset),
         lty = c("dashed", rep("solid", length(colorset))),
         legend = c("overall", levels(data$site)),
         bty = "n",
         cex = 0.7
  )
}

virologicfailure.cscox = function(data, quantile=.95, complete=FALSE, includeHaiti=FALSE, includeHonduras=FALSE, merge=FALSE){
  # remove initial ART regimen of "Other" type
  data = data %>% subset(initialARTtype !="Other")
  data$initialARTtype = factor(data$initialARTtype)
  #
  data$vlfailure_status.date = virologicfailure.date(data = data, quantile = quantile)
  data$vlfailure_status.days = as.integer(data$vlfailure_status.date) - as.integer(data$recart_d)
  # some have censored dates prior to the date of receiving initial ART due to the definition of closing dates
  data = data %>% subset(vlfailure_status.days>=0)
  # define time on study in years
  data$vlfailure_status.years = data$vlfailure_status.days/365.25
  if(quantile==.95){
    data$vlfailure_status = factor(data$vlfailure_status.95, levels = c(0,1,2), labels = c("right censored", "death", "virologic failure"))
  }
  else{
    data$vlfailure_status = factor(data$vlfailure_status.90, levels = c(0,1,2), labels = c("right censored", "death", "virologic failure"))
  }
  if(!includeHaiti){
    data = data %>% subset(site!="haiti")
  }
  if(!includeHonduras){
    data = data %>% subset(site!="honduras")
  }
  data$site = factor(data$site)
  if(complete){ # without imputation
    if(!merge){
      data = data %>% subset(!is.na(sex)  & !is.na(mode_new) & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log) & !is.na(vlfailure_status))     
      fit = coxph(Surv(vlfailure_status.years, vlfailure_status=="virologic failure")~strata(site)+sex+age+mode_new+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+calendaryear+initialARTtype, data = data)      
    }
    else{
      data = data %>% subset(!is.na(sex_mode)  & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log) & !is.na(vlfailure_status))
      fit = coxph(Surv(vlfailure_status.years, vlfailure_status=="virologic failure")~strata(site)+sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+calendaryear+initialARTtype, data = data)      
    }
  }
  else{ # with imputation
    var_model = c("vlfailure_status.years", "vlfailure_status", "site", "sex", "age", "mode_new", "clinicalAIDS", "cd4_v_initial_sqrt", "rna_v_initial_log", "calendaryear", "initialARTtype")
    imputed_data = mice(data %>% select(all_of(var_model)) , m = 5, method = "pmm", seed = 1117, printFlag = F)
    if(!merge){
      fit <- with(imputed_data,
                  coxph(Surv(vlfailure_status.years, vlfailure_status=="virologic failure") ~ strata(site)+sex+age+mode_new+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+calendaryear+initialARTtype))      
    }
    else{
      # derive sex_mode
      imputed_data_d = complete(imputed_data, action = "long", include = TRUE)
      imputed_data_d$sex_mode = sapply(1:nrow(imputed_data_d), function(i) {
        ifelse(imputed_data_d$sex[i]=="Female", "Female",
               ifelse(imputed_data_d$mode_new=="MSM", "Male-MSM", "Male-Other"))
      }) %>% factor(levels = c("Female", "Male-MSM", "Male-Other"))
      imputed_data_d = as.mids(imputed_data_d)
      fit <- with(imputed_data_d,
                  coxph(Surv(vlfailure_status.years, vlfailure_status=="virologic failure") ~ strata(site)+sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+calendaryear+initialARTtype))
    }
    fit = summary(pool(fit))
  }
  return(list(vl=data, fit=fit)) 
}
