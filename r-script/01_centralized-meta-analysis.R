#------------------------------------------------------------------------------------------#
#   Centralized and Meta-analysis for CCASAnet Federated Survival Analysis Study          #
#------------------------------------------------------------------------------------------#
# This script performs:
# 1. Data import and preparation
# 2. Summary statistics by site
# 3. Centralized stratified Cox regression (complete case and multiple imputation)
# 4. Meta-analysis (complete case and multiple imputation)
#
# Analyses focus on 3 sites: Brazil, Chile, Mexico
# Outcomes: (1) virologic failure; (2) major regimen change
# Closing date defined as the .95 quantile of the distribution of last observations per site


# --- Load libraries and source files --- #

source("./main/loadlibrary-ccasanet.R")
source("./main/ccasanet2024.R")
source("./main/centralized-analysis-ccasanet2024.R")
source("./main/distributed-analysis-ccasanet.R")


# --- Data import and preparation --- #

# import raw data
datalist = importccasanet()
datalist = convertdate(datalist)

# unpack data
basic  = datalist$basic
art    = datalist$art
visit  = datalist$visit
follow = datalist$follow
ce     = datalist$ce
cd4    = datalist$cd4
rna    = datalist$rna
center = datalist$center

# derive variables
## probable route of infection
basic$mode_new = modecategory(basic$mode, basic$gender, basic$mode_oth)
## clinical AIDS status at ART initiation
basic$clinicalAIDS = clinicalaids(basic$aids_first_y, basic$aids_first_d, basic$recart_d)
## CD4 count at/near ART initiation (within 180 days before and 7 days after)
basic$cd4_v_initial = cd4atinitialART(cd4, basic)
basic$cd4_v_initial_sqrt = sqrt(basic$cd4_v_initial)
## viral load at/near ART initiation (within 180 days before ART initiation)
basic$rna_v_initial = vlatinitialART(rna, basic)
basic$rna_v_initial_log = log10(basic$rna_v_initial)
## initial ART type (COBI reclassified as non-PI)
art$pi_new = definenewpi(art)
basic$initialARTtype = identifyinitialARTtype(art, basic)
## calendar year and age at ART initiation
basic$calendaryear = as.integer(format(basic$recart_d, "%Y"))
basic$age = as.numeric(basic$recart_d - basic$birth_d) / 365.25
## last observation date and death date
basic$lastobs_d = getlastobs_d(follow, basic)
basic$death_y_d = follow$death_d[match(basic$patient_id, follow$patient_id)]

# identify outcomes
## virologic failure
vlfailure = identifyvlfailure(rna, basic)
basic$vlfailure_y   = vlfailure$vlfailure_y
basic$vlfailure_y_d = vlfailure$vlfailure_y_d
## major regimen change
regimenchange = identifyregimenchange(art, basic)
basic$ARTchange_y     = regimenchange$ARTchange_y
basic$ARTchange_y_d   = regimenchange$ARTchange_y_d
basic$ARTchange_y_oth = regimenchange$ARTchange_y_oth

# define site-specific closing dates
## .95 quantile of the distribution of last observations per site
closedate.95 = sapply(levels(factor(basic$site)), function(x) {
  quantile(basic$lastobs_d[basic$site==x], probs = .95, na.rm = TRUE) %>% as.Date()
})
## .90 quantile
closedate.90 = sapply(levels(factor(basic$site)), function(x) {
  quantile(basic$lastobs_d[basic$site==x], probs = .90, na.rm = TRUE) %>% as.Date()
})

# add site-specific closing dates to basic
basic$siteclosedate.95 = closedate.95[basic$site]
basic$siteclosedate.90 = closedate.90[basic$site]

# define outcome status variables (0 = right censored, 1 = death, 2 = event)
basic$ARTchange_status.95 = identifyregimenchange.status(basic, follow, closedate.95)
basic$ARTchange_status.90 = identifyregimenchange.status(basic, follow, closedate.90)
basic$vlfailure_status.95 = identifyvlfailure.status(basic, follow, closedate.95)
basic$vlfailure_status.90 = identifyvlfailure.status(basic, follow, closedate.90)

# sex + probable route of infection combined variable
basic$sex_mode = sapply(1:nrow(basic), function(i) {
  ifelse(is.na(basic$sex[i]), NA,
         ifelse(basic$sex[i]=="Female", "Female",
                ifelse(is.na(basic$mode_new[i]), NA,
                       ifelse(basic$mode_new[i]=="MSM", "Male-MSM", "Male-Other"))))
}) %>% factor(levels = c("Female", "Male-MSM", "Male-Other"))

# create analysis dataset
hiv = basic


# --- Variables of interest --- #

## major regimen change
trt.var.model = c("ARTchange_status.years", "ARTchange_status", "site", "sex_mode", "age",
                  "clinicalAIDS", "cd4_v_initial_sqrt", "rna_v_initial_log", "calendaryear", "initialARTtype")
## virologic failure
vl.var.model = c("vlfailure_status.years", "vlfailure_status", "site", "sex_mode", "age",
                 "clinicalAIDS", "cd4_v_initial_sqrt", "rna_v_initial_log", "calendaryear", "initialARTtype")

# exclude 'Other' initial ART type for all analyses
hiv = hiv %>% subset(initialARTtype!="Other")
hiv$initialARTtype = factor(hiv$initialARTtype, levels = c("NNRTI", "PI", "INSTI"))


# --- Summary tables --- #

# all sites (excluding Peru)
hiv %>% select(all_of(intersect(trt.var.model, vl.var.model))) %>%
  tbl_summary(by=site, missing = "ifany",
              label = list(sex_mode = "Sex+Probable Route of Infection", age="Age",
                           clinicalAIDS = "Clinical AIDS",
                           cd4_v_initial_sqrt = "CD4 count square root",
                           rna_v_initial_log = "Log10 viral load",
                           initialARTtype = "Initial ART",
                           calendaryear = "Calendar year"),
              type = list(calendaryear ~ "continuous"),
              missing_text = "NA") %>%
  add_overall() %>% bold_labels() %>% as_gt() %>%
  gt::tab_header(title = "Summary Table by site",
                 subtitle = "Excluding Peru")

# 3 sites only: Brazil, Chile, Mexico
hiv %>% select(all_of(intersect(trt.var.model, vl.var.model))) %>%
  subset(site=="brazil" | site=="chile" | site=="mexico") %>%
  mutate(site = factor(site, levels=c("brazil", "chile", "mexico"))) %>%
  tbl_summary(by=site, missing = "ifany",
              label = list(sex_mode = "Sex+Probable Route of Infection", age="Age",
                           clinicalAIDS = "Clinical AIDS",
                           cd4_v_initial_sqrt = "CD4 count square root",
                           rna_v_initial_log = "Log10 viral load",
                           initialARTtype = "Initial ART",
                           calendaryear = "Calendar year"),
              type = list(calendaryear ~ "continuous"),
              missing_text = "NA") %>%
  add_overall() %>% bold_labels() %>% as_gt() %>%
  gt::tab_header(title = "Summary Table by site",
                 subtitle = "Only keep 3 sites")


# --- Centralized analysis: complete case --- #

# use .95 quantile closing date; exclude Haiti; merge sex and route of infection
trt.95 = regimenchange.cscox(data = hiv, quantile = .95, complete = TRUE, includeHaiti = FALSE, merge = TRUE)
vl.95  = virologicfailure.cscox(data = hiv, quantile = .95, complete = TRUE, includeHaiti = FALSE, merge = TRUE)

# focus on 3 sites: Brazil, Chile, Mexico
trt.95.3sites = trt.95$trt %>%
  subset(site=="brazil" | site=="chile" | site=="mexico") %>%
  mutate(site = factor(site, levels=c("brazil", "chile", "mexico")))
vl.95.3sites = vl.95$vl %>%
  subset(site=="brazil" | site=="chile" | site=="mexico") %>%
  mutate(site = factor(site, levels=c("brazil", "chile", "mexico")))

## major regimen change
trt.95.3sites.fit = coxph(Surv(ARTchange_status.years, ARTchange_status=="major treatment change") ~
                            strata(site)+sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+
                            rna_v_initial_log+calendaryear+initialARTtype,
                          data = trt.95.3sites)
### sample size
data.frame(site = c("brazil", "chile", "mexico", "total"),
           sample.size = c(table(trt.95.3sites$site) %>% unname(), nrow(trt.95.3sites)),
           event.size  = c(sapply(levels(trt.95.3sites$site), function(x)
             trt.95.3sites %>% subset(site==x & ARTchange_status=="major treatment change") %>% nrow()) %>% unname(),
             trt.95.3sites %>% subset(ARTchange_status=="major treatment change") %>% nrow()))
### results
cbind(summary(trt.95.3sites.fit)$conf.int[,-2],
      std.error = sqrt(diag(trt.95.3sites.fit$var)),
      pvalue = summary(trt.95.3sites.fit)$coefficients[, "Pr(>|z|)"])

## virologic failure
vl.95.3sites.fit = coxph(Surv(vlfailure_status.years, vlfailure_status=="virologic failure") ~
                           strata(site)+sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+
                           rna_v_initial_log+calendaryear+initialARTtype,
                         data = vl.95.3sites)
### sample size
data.frame(site = c("brazil", "chile", "mexico", "total"),
           sample.size = c(table(vl.95.3sites$site) %>% unname(), nrow(vl.95.3sites)),
           event.size  = c(sapply(levels(vl.95.3sites$site), function(x)
             vl.95.3sites %>% subset(site==x & vlfailure_status=="virologic failure") %>% nrow()) %>% unname(),
             vl.95.3sites %>% subset(vlfailure_status=="virologic failure") %>% nrow()))
### results
cbind(summary(vl.95.3sites.fit)$conf.int[,-2],
      std.error = sqrt(diag(vl.95.3sites.fit$var)),
      pvalue = summary(vl.95.3sites.fit)$coefficients[, "Pr(>|z|)"])


# --- Centralized analysis: multiple imputation --- #

m = 10 # number of imputed datasets

## variables for imputation
trt.var.model.imp = c("ARTchange_status.years", "ARTchange_status", "site", "sex", "age",
                      "mode_new", "clinicalAIDS", "cd4_v_initial_sqrt", "rna_v_initial_log",
                      "calendaryear", "initialARTtype")
vl.var.model.imp  = c("vlfailure_status.years", "vlfailure_status", "site", "sex", "age",
                      "mode_new", "clinicalAIDS", "cd4_v_initial_sqrt", "rna_v_initial_log",
                      "calendaryear", "initialARTtype")

## prepare analysis datasets with .95 quantile closing date
trt.95.mi = regimenchange.cscox(data = hiv, quantile = .95, complete = FALSE, includeHaiti = FALSE, merge = TRUE)
vl.95.mi  = virologicfailure.cscox(data = hiv, quantile = .95, complete = FALSE, includeHaiti = FALSE, merge = TRUE)

## focus on 3 sites
trt.95.mi.3sites = trt.95.mi$trt %>%
  subset(site=="brazil" | site=="chile" | site=="mexico") %>%
  mutate(site = factor(site, levels=c("brazil", "chile", "mexico")))
vl.95.mi.3sites = vl.95.mi$vl %>%
  subset(site=="brazil" | site=="chile" | site=="mexico") %>%
  mutate(site = factor(site, levels=c("brazil", "chile", "mexico")))

## local multiple imputation by site (predictive mean matching)
### major regimen change
trt.95.mi.3sites.split = split(trt.95.mi.3sites, trt.95.mi.3sites$site)
trt.95.mi.3sites.split.mi = sapply(trt.95.mi.3sites.split, function(x)
  mice(x %>% select(all_of(trt.var.model.imp)), m=m, method = "pmm", seed = 0108, printFlag = FALSE),
  simplify = F)
trt.95.mi.3sites.split.mi = sapply(trt.95.mi.3sites.split.mi, function(x) {
  sapply(1:m, function(y) complete(x, action=y) %>%
           mutate(sex_mode = ifelse(sex=="Female", "Female",
                                   ifelse(mode_new=="MSM", "Male-MSM", "Male-Other")) %>%
                    factor(levels = c("Female", "Male-MSM", "Male-Other"))), simplify = F)
}, simplify = F)
## pool sites for each imputed dataset
trt.95.mi.3sites.split.mi.d = sapply(1:m, function(x) {
  rbind(trt.95.mi.3sites.split.mi[["brazil"]][[x]],
        trt.95.mi.3sites.split.mi[["chile"]][[x]],
        trt.95.mi.3sites.split.mi[["mexico"]][[x]])
}, simplify = F)
## fit Cox regression on each imputed dataset and pool via Rubin's rules
trt.95.mi.3sites.split.mi.fit = sapply(trt.95.mi.3sites.split.mi.d, function(x) {
  coxph(Surv(ARTchange_status.years, ARTchange_status=="major treatment change") ~
          strata(site)+sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+
          rna_v_initial_log+calendaryear+initialARTtype, data = x)
}, simplify = F)
trt.95.mi.3sites.split.mi.fit = sapply(trt.95.mi.3sites.split.mi.fit, function(x) {
  tmp = summary(x)$coefficients[,c(1:3, 5)]
  colnames(tmp) = c("estimate","hr","std.error","pvalue")
  ci.tmp = summary(x)$conf.int[,3:4]
  as.data.frame(cbind(tmp, ci.tmp))
}, simplify = F)
trt.95.mi.3sites.split.mi.fit = mi.pool.rubin(trt.95.mi.3sites.split.mi.fit,
                                               nrow(trt.95.mi.3sites.split.mi.d[[1]]))
### results
trt.95.mi.3sites.split.mi.fit %>%
  mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
         upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
  select(all_of(c("hr","lower.95","upper.95", "std.error", "p.value")))

### virologic failure
vl.95.mi.3sites.split = split(vl.95.mi.3sites, vl.95.mi.3sites$site)
vl.95.mi.3sites.split.mi = sapply(vl.95.mi.3sites.split, function(x)
  mice(x %>% select(all_of(vl.var.model.imp)), m=m, method = "pmm", seed = 0108, printFlag = FALSE),
  simplify = F)
vl.95.mi.3sites.split.mi = sapply(vl.95.mi.3sites.split.mi, function(x) {
  sapply(1:m, function(y) complete(x, action=y) %>%
           mutate(sex_mode = ifelse(sex=="Female", "Female",
                                   ifelse(mode_new=="MSM", "Male-MSM", "Male-Other")) %>%
                    factor(levels = c("Female", "Male-MSM", "Male-Other"))), simplify = F)
}, simplify = F)
## pool sites for each imputed dataset
vl.95.mi.3sites.split.mi.d = sapply(1:m, function(x) {
  rbind(vl.95.mi.3sites.split.mi[["brazil"]][[x]],
        vl.95.mi.3sites.split.mi[["chile"]][[x]],
        vl.95.mi.3sites.split.mi[["mexico"]][[x]])
}, simplify = F)
## fit Cox regression on each imputed dataset and pool via Rubin's rules
vl.95.mi.3sites.split.mi.fit = sapply(vl.95.mi.3sites.split.mi.d, function(x) {
  coxph(Surv(vlfailure_status.years, vlfailure_status=="virologic failure") ~
          strata(site)+sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+
          rna_v_initial_log+calendaryear+initialARTtype, data = x)
}, simplify = F)
vl.95.mi.3sites.split.mi.fit = sapply(vl.95.mi.3sites.split.mi.fit, function(x) {
  tmp = summary(x)$coefficients[,c(1:3,5)]
  colnames(tmp) = c("estimate","hr","std.error","pvalue")
  as.data.frame(tmp)
}, simplify = F)
vl.95.mi.3sites.split.mi.fit = mi.pool.rubin(vl.95.mi.3sites.split.mi.fit,
                                              nrow(vl.95.mi.3sites.split.mi.d[[1]]))
### results
vl.95.mi.3sites.split.mi.fit %>%
  mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
         upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
  select(all_of(c("hr","lower.95","upper.95", "std.error", "p.value")))


# --- Meta-analysis --- #

## complete case
## local estimates are obtained from *_initialize.json files generated during the
## simulated federated analysis (see 02_simulated-federated-analysis.R), which contain
## site-specific Cox regression estimates (bhat_i) and variances (Vhat_i)
### major regimen change
site_1_AD = fromJSON(file = "./regimenchange/chile/complete-case/chile_initialize.json")
site_2_AD = fromJSON(file = "./regimenchange/chile/complete-case/brazil_initialize.json")
site_3_AD = fromJSON(file = "./regimenchange/chile/complete-case/mexico_initialize.json")
trt.95.3sites.fit.meta = meta(list(site_1_AD, site_2_AD, site_3_AD))
trt.95.3sites.fit.meta %>%
  select(-all_of(c("estimate"))) %>%
  mutate(name=names(coefficients(trt.95.3sites.fit))) %>%
  column_to_rownames(var = "name") %>%
  select(all_of(c("hr", "lower.95", "upper.95", "std.error", "p.value")))

### virologic failure
site_1_AD = fromJSON(file = "./virologicfailure/chile/complete-case/chile_initialize.json")
site_2_AD = fromJSON(file = "./virologicfailure/chile/complete-case/brazil_initialize.json")
site_3_AD = fromJSON(file = "./virologicfailure/chile/complete-case/mexico_initialize.json")
vl.95.3sites.fit.meta = meta(list(site_1_AD, site_2_AD, site_3_AD))
vl.95.3sites.fit.meta %>%
  select(-all_of(c("estimate"))) %>%
  mutate(name=names(coefficients(vl.95.3sites.fit))) %>%
  column_to_rownames(var = "name") %>%
  select(all_of(c("hr", "lower.95", "upper.95", "std.error", "p.value")))

## multiple imputation
### major regimen change
trt.95.mi.3sites.split.mi.fit.meta = sapply(trt.95.mi.3sites.split.mi, function(x) {
  sapply(x, function(y)
    coxph(Surv(ARTchange_status.years, ARTchange_status=="major treatment change") ~
            sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+
            calendaryear+initialARTtype, data = y), simplify = F)
}, simplify = F)
trt.95.mi.3sites.split.mi.fit.meta = sapply(trt.95.mi.3sites.split.mi.fit.meta, function(x) {
  summary(pool(x))
}, simplify = F)
trt.95.mi.3sites.split.mi.fit.meta = meta.mi(trt.95.mi.3sites.split.mi.fit.meta)
trt.95.mi.3sites.split.mi.fit.meta %>%
  select(-all_of(c("estimate"))) %>%
  mutate(name=names(coefficients(trt.95.3sites.fit))) %>%
  column_to_rownames(var = "name") %>%
  select(all_of(c("hr", "lower.95", "upper.95", "std.error", "p.value")))

### virologic failure
vl.95.mi.3sites.split.mi.fit.meta = sapply(vl.95.mi.3sites.split.mi, function(x) {
  sapply(x, function(y)
    coxph(Surv(vlfailure_status.years, vlfailure_status=="virologic failure") ~
            sex_mode+age+clinicalAIDS+cd4_v_initial_sqrt+rna_v_initial_log+
            calendaryear+initialARTtype, data = y), simplify = F)
}, simplify = F)
vl.95.mi.3sites.split.mi.fit.meta = sapply(vl.95.mi.3sites.split.mi.fit.meta, function(x) {
  summary(pool(x))
}, simplify = F)
vl.95.mi.3sites.split.mi.fit.meta = meta.mi(vl.95.mi.3sites.split.mi.fit.meta)
vl.95.mi.3sites.split.mi.fit.meta %>%
  select(-all_of(c("estimate"))) %>%
  mutate(name=names(coefficients(vl.95.3sites.fit))) %>%
  column_to_rownames(var = "name") %>%
  select(all_of(c("hr", "lower.95", "upper.95", "std.error", "p.value")))
