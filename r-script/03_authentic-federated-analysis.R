#------------------------------------------------------------------------------------------#
#    Authentic Federated Analysis for CCASAnet Federated Survival Analysis Study          #
#------------------------------------------------------------------------------------------#
# This script implements the authentic federated analysis protocol.
# Unlike the simulated federated analysis, data are never pooled across sites.
# Each participating site (Brazil, Chile, Mexico) runs its own local code and uploads
# intermediate JSON files to a shared folder (e.g., OneDrive); only the lead site
# aggregates results. Chile is the lead site in all analyses.
#
# The pda package (Privacy-preserving Distributed Algorithms) implements the
# ODAC algorithm for distributed stratified Cox regression.
#
# Outcomes: (1) virologic failure; (2) major regimen change
# Analysis scenarios: complete case and multiple imputation (m imputed datasets)
#
# Communication protocol (for each outcome x scenario):
#
#   Step 1 [Chile]:         pda.setup()   -- generate control.json; upload to shared folder;
#                                            inform Brazil and Mexico to download
#   Step 2a [Brazil]:       pda.local()   -- initialize; upload brazil_initialize.json;
#                                            inform Chile to download
#   Step 2b [Mexico]:       pda.local()   -- initialize; upload mexico_initialize.json;
#                                            inform Chile to download
#   Step 2c [Chile]:        pda.local()   -- initialize after downloading brazil/mexico
#                                            _initialize.json; upload updated control.json;
#                                            inform Brazil and Mexico to download
#   Step 3a [Brazil]:       pda.derive()  -- derive; upload brazil_derive.json;
#                                            inform Chile to download
#   Step 3b [Mexico]:       pda.derive()  -- derive; upload mexico_derive.json;
#                                            inform Chile to download
#   Step 3c [Chile]:        pda.derive()  -- derive after downloading brazil/mexico
#                                            _derive.json; move to estimation
#   Step 4 [Chile]:         pda.estimate() -- aggregate and estimate
#   Step 5 [Chile]:         pda.summary()  -- extract and pool results
#                           dist.meta()    -- meta-analysis from local estimates
#
# See authentic-implementation/ for site-specific HTML protocols.


# --- Load libraries and source files --- #

source("./main/loadlibrary-ccasanet.R")
source("./main/ccasanet2024.R")
source("./main/distributed-analysis-ccasanet.R")
source("./main/authentic-distributed-analysis.R")


# --- Site-specific data preparation --- #
# Each site prepares its own local dataset independently.
# The setup object contains shared parameters (seed, number of imputations, variable lists).

## shared setup parameters
setup = list(
  seed             = 0108,
  no.of.imputation = 10,
  vl.imp.var  = c("vlfailure_status.years", "vlfailure_status", "sex", "age",
                  "mode_new", "clinicalAIDS", "cd4_v_initial_sqrt", "rna_v_initial_log",
                  "calendaryear", "initialARTtype"),
  trt.imp.var = c("ARTchange_status.years", "ARTchange_status", "sex", "age",
                  "mode_new", "clinicalAIDS", "cd4_v_initial_sqrt", "rna_v_initial_log",
                  "calendaryear", "initialARTtype")
)
m = setup[["no.of.imputation"]]

## each site loads its own data (site-specific files are not included in this repository)
## the code below illustrates the Chile site; Brazil and Mexico follow the same structure

# setwd("your-working-directory") # set to the site's local working directory

## Chile
chile.vl  = readRDS("chile_all_vl.rds")
chile.trt = readRDS("chile_all_trt.rds")
### complete case
chile.vl.cc  = chile.vl  %>% subset(!is.na(sex_mode) & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log))
chile.trt.cc = chile.trt %>% subset(!is.na(sex_mode) & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log))
### multiple imputation (local, per site)
chile.vl.mi  = mice.local.merge(
  mice(chile.vl  %>% select(all_of(setup[["vl.imp.var"]])),
       m=m, method = "pmm", seed = setup[["seed"]], printFlag = FALSE),
  m = m)
chile.trt.mi = mice.local.merge(
  mice(chile.trt %>% select(all_of(setup[["trt.imp.var"]])),
       m=m, method = "pmm", seed = setup[["seed"]], printFlag = FALSE),
  m = m)

## Brazil (run at the Brazil site)
brazil.vl  = readRDS("brazil_all_vl.rds")
brazil.trt = readRDS("brazil_all_trt.rds")
brazil.vl.cc  = brazil.vl  %>% subset(!is.na(sex_mode) & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log))
brazil.trt.cc = brazil.trt %>% subset(!is.na(sex_mode) & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log))
brazil.vl.mi  = mice.local.merge(
  mice(brazil.vl  %>% select(all_of(setup[["vl.imp.var"]])),
       m=m, method = "pmm", seed = setup[["seed"]], printFlag = FALSE),
  m = m)
brazil.trt.mi = mice.local.merge(
  mice(brazil.trt %>% select(all_of(setup[["trt.imp.var"]])),
       m=m, method = "pmm", seed = setup[["seed"]], printFlag = FALSE),
  m = m)

## Mexico (run at the Mexico site)
mexico.vl  = readRDS("mexico_all_vl.rds")
mexico.trt = readRDS("mexico_all_trt.rds")
mexico.vl.cc  = mexico.vl  %>% subset(!is.na(sex_mode) & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log))
mexico.trt.cc = mexico.trt %>% subset(!is.na(sex_mode) & !is.na(cd4_v_initial_sqrt) & !is.na(rna_v_initial_log))
mexico.vl.mi  = mice.local.merge(
  mice(mexico.vl  %>% select(all_of(setup[["vl.imp.var"]])),
       m=m, method = "pmm", seed = setup[["seed"]], printFlag = FALSE),
  m = m)
mexico.trt.mi = mice.local.merge(
  mice(mexico.trt %>% select(all_of(setup[["trt.imp.var"]])),
       m=m, method = "pmm", seed = setup[["seed"]], printFlag = FALSE),
  m = m)


# ======================================================================================== #
#                        VIROLOGIC FAILURE -- COMPLETE CASE                                #
# ======================================================================================== #

# Step 1 [Chile -- lead site]: define control and generate control.json
control <- list(project_name = 'virologic_failure_distributed_survival_analysis_ccasanet',
                step         = 'initialize',
                sites        = c("chile", "brazil", "mexico"),
                heterogeneity = TRUE,
                model        = 'ODAC',
                family       = 'cox',
                outcome      = "Surv(vlfailure_status.years, vlfailure_status==\"virologic failure\")",
                variables    = c('sex_mode', 'age', 'clinicalAIDS', 'cd4_v_initial_sqrt',
                                 'rna_v_initial_log', 'calendaryear', 'initialARTtype'),
                xlev         = list(sex_mode     = c('Female','Male-MSM', "Male-Other"),
                                    clinicalAIDS = c('Not clinical AIDS','Clinical AIDS'),
                                    initialARTtype = c('NNRTI','PI','INSTI')),
                optim_maxit  = 10000,
                init_method  = "meta",
                optim_method = "Nelder-Mead",
                lead_site    = "chile",
                upload_date  = as.character(Sys.time()))
set.seed(setup[["seed"]])
pda.setup(leadsite = "chile", control, dir = getwd(), m = 1, complete = TRUE, outcome = "vl")
# [Chile] upload control.json to shared folder; inform Brazil and Mexico to download

# Step 2 [Initialization]: Brazil and Mexico run first, then Chile
## Step 2a [Brazil]: download control.json; initialize; upload brazil_initialize.json
set.seed(setup[["seed"]])
pda.local(brazil.vl.cc, "brazil", getwd(), m = 1, complete = TRUE, outcome = "vl")
# [Brazil] upload brazil_initialize.json to shared folder; inform Chile to download

## Step 2b [Mexico]: download control.json; initialize; upload mexico_initialize.json
set.seed(setup[["seed"]])
pda.local(mexico.vl.cc, "mexico", getwd(), m = 1, complete = TRUE, outcome = "vl")
# [Mexico] upload mexico_initialize.json to shared folder; inform Chile to download

## Step 2c [Chile]: download brazil_initialize.json and mexico_initialize.json; initialize;
##                  upload updated control.json; inform Brazil and Mexico to download
set.seed(setup[["seed"]])
pda.local(chile.vl.cc, "chile", getwd(), m = 1, complete = TRUE, outcome = "vl")
# [Chile] upload updated control.json to shared folder; inform Brazil and Mexico to download

# Step 3 [Derivation]: Brazil and Mexico run first, then Chile
## Step 3a [Brazil]: download updated control.json; derive; upload brazil_derive.json
set.seed(setup[["seed"]])
pda.derive(brazil.vl.cc, "brazil", getwd(), m = 1, complete = TRUE, outcome = "vl")
# [Brazil] upload brazil_derive.json to shared folder; inform Chile to download

## Step 3b [Mexico]: download updated control.json; derive; upload mexico_derive.json
set.seed(setup[["seed"]])
pda.derive(mexico.vl.cc, "mexico", getwd(), m = 1, complete = TRUE, outcome = "vl")
# [Mexico] upload mexico_derive.json to shared folder; inform Chile to download

## Step 3c [Chile]: download brazil_derive.json and mexico_derive.json; derive
set.seed(setup[["seed"]])
pda.derive(chile.vl.cc, "chile", getwd(), m = 1, complete = TRUE, outcome = "vl")
# [Chile] move to estimation

# Step 4 [Chile -- lead site]: estimate
set.seed(setup[["seed"]])
pda.estimate(chile.vl.cc, "chile", getwd(), m = 1, complete = TRUE, outcome = "vl")

# Step 5 [Chile -- lead site]: extract results and save
vl.dist.cc = pda.summary("chile", c("brazil", "mexico"), getwd(), 1, complete = TRUE, outcome = "vl")
saveRDS(vl.dist.cc, file = "vl_dist_cc.rds")
vl.meta.cc = meta.mi(dist.meta(c("chile", "brazil", "mexico"), getwd(), 1, complete = TRUE, outcome = "vl"))
saveRDS(vl.meta.cc, file = "vl_meta_cc.rds")


# ======================================================================================== #
#                     VIROLOGIC FAILURE -- MULTIPLE IMPUTATION                            #
# ======================================================================================== #

# Step 1 [Chile -- lead site]: generate m control.json files (one per imputed dataset)
control <- list(project_name = 'virologic_failure_distributed_survival_analysis_ccasanet',
                step         = 'initialize',
                sites        = c("chile", "brazil", "mexico"),
                heterogeneity = TRUE,
                model        = 'ODAC',
                family       = 'cox',
                outcome      = "Surv(vlfailure_status.years, vlfailure_status==\"virologic failure\")",
                variables    = c('sex_mode', 'age', 'clinicalAIDS', 'cd4_v_initial_sqrt',
                                 'rna_v_initial_log', 'calendaryear', 'initialARTtype'),
                xlev         = list(sex_mode     = c('Female','Male-MSM', "Male-Other"),
                                    clinicalAIDS = c('Not clinical AIDS','Clinical AIDS'),
                                    initialARTtype = c('NNRTI','PI','INSTI')),
                optim_maxit  = 10000,
                init_method  = "meta",
                optim_method = "Nelder-Mead",
                lead_site    = "chile",
                upload_date  = as.character(Sys.time()))
set.seed(setup[["seed"]])
pda.setup(leadsite = "chile", control, dir = getwd(), m = m, outcome = "vl")
# [Chile] upload m control.json files to shared folder (one sub-folder per imputed dataset);
#         inform Brazil and Mexico to download

# Step 2 [Initialization]: Brazil and Mexico run first, then Chile
## Step 2a [Brazil]: download control.json files; initialize for each imputed dataset;
##                   upload brazil_initialize.json for each imputed dataset
set.seed(setup[["seed"]])
pda.local(brazil.vl.mi, "brazil", getwd(), m = m, outcome = "vl")
# [Brazil] upload brazil_initialize.json (x m) to shared folder; inform Chile to download

## Step 2b [Mexico]: download control.json files; initialize for each imputed dataset;
##                   upload mexico_initialize.json for each imputed dataset
set.seed(setup[["seed"]])
pda.local(mexico.vl.mi, "mexico", getwd(), m = m, outcome = "vl")
# [Mexico] upload mexico_initialize.json (x m) to shared folder; inform Chile to download

## Step 2c [Chile]: download brazil/mexico _initialize.json (x m); initialize;
##                  upload updated control.json (x m); inform Brazil and Mexico to download
set.seed(setup[["seed"]])
pda.local(chile.vl.mi, "chile", getwd(), m = m, outcome = "vl")
# [Chile] upload updated control.json (x m) to shared folder; inform Brazil and Mexico to download

# Step 3 [Derivation]: Brazil and Mexico run first, then Chile
## Step 3a [Brazil]: download updated control.json files; derive for each imputed dataset;
##                   upload brazil_derive.json for each imputed dataset
set.seed(setup[["seed"]])
pda.derive(brazil.vl.mi, "brazil", getwd(), m = m, outcome = "vl")
# [Brazil] upload brazil_derive.json (x m) to shared folder; inform Chile to download

## Step 3b [Mexico]: download updated control.json files; derive for each imputed dataset;
##                   upload mexico_derive.json for each imputed dataset
set.seed(setup[["seed"]])
pda.derive(mexico.vl.mi, "mexico", getwd(), m = m, outcome = "vl")
# [Mexico] upload mexico_derive.json (x m) to shared folder; inform Chile to download

## Step 3c [Chile]: download brazil/mexico _derive.json (x m); derive
set.seed(setup[["seed"]])
pda.derive(chile.vl.mi, "chile", getwd(), m = m, outcome = "vl")
# [Chile] move to estimation

# Step 4 [Chile -- lead site]: estimate for each imputed dataset
set.seed(setup[["seed"]])
pda.estimate(chile.vl.mi, "chile", getwd(), m = m, outcome = "vl")

# Step 5 [Chile -- lead site]: pool results across imputed datasets (Rubin's rules) and save
vl.dist.mi = pda.summary("chile", c("brazil", "mexico"), getwd(), m = m, outcome = "vl")
saveRDS(vl.dist.mi, file = "vl_dist_mi.rds")
vl.meta.mi = meta.mi(dist.meta(c("chile", "brazil", "mexico"), getwd(), m = m, outcome = "vl"))
saveRDS(vl.meta.mi, file = "vl_meta_mi.rds")


# ======================================================================================== #
#                      MAJOR REGIMEN CHANGE -- COMPLETE CASE                              #
# ======================================================================================== #

# Step 1 [Chile -- lead site]: define control and generate control.json
control <- list(project_name  = 'major_regimen_change_distributed_survival_analysis_ccasanet',
                step          = 'initialize',
                sites         = c("chile", "brazil", "mexico"),
                heterogeneity = TRUE,
                model         = 'ODAC',
                family        = 'cox',
                outcome       = "Surv(ARTchange_status.years, ARTchange_status==\"major treatment change\")",
                variables     = c('sex_mode', 'age', 'clinicalAIDS', 'cd4_v_initial_sqrt',
                                  'rna_v_initial_log', 'calendaryear', 'initialARTtype'),
                xlev          = list(sex_mode     = c('Female','Male-MSM', 'Male-Other'),
                                     clinicalAIDS = c('Not clinical AIDS','Clinical AIDS'),
                                     initialARTtype = c('NNRTI','PI','INSTI')),
                optim_maxit   = 10000,
                init_method   = "meta",
                optim_method  = "Nelder-Mead",
                lead_site     = "chile",
                upload_date   = as.character(Sys.time()))
set.seed(setup[["seed"]])
pda.setup(leadsite = "chile", control, dir = getwd(), m = 1, complete = TRUE, outcome = "trt")
# [Chile] upload control.json to shared folder; inform Brazil and Mexico to download

# Step 2 [Initialization]: Brazil and Mexico run first, then Chile
## Step 2a [Brazil]: download control.json; initialize; upload brazil_initialize.json
set.seed(setup[["seed"]])
pda.local(brazil.trt.cc, "brazil", getwd(), m = 1, complete = TRUE, outcome = "trt")
# [Brazil] upload brazil_initialize.json to shared folder; inform Chile to download

## Step 2b [Mexico]: download control.json; initialize; upload mexico_initialize.json
set.seed(setup[["seed"]])
pda.local(mexico.trt.cc, "mexico", getwd(), m = 1, complete = TRUE, outcome = "trt")
# [Mexico] upload mexico_initialize.json to shared folder; inform Chile to download

## Step 2c [Chile]: download brazil_initialize.json and mexico_initialize.json; initialize;
##                  upload updated control.json; inform Brazil and Mexico to download
set.seed(setup[["seed"]])
pda.local(chile.trt.cc, "chile", getwd(), m = 1, complete = TRUE, outcome = "trt")
# [Chile] upload updated control.json to shared folder; inform Brazil and Mexico to download

# Step 3 [Derivation]: Brazil and Mexico run first, then Chile
## Step 3a [Brazil]: download updated control.json; derive; upload brazil_derive.json
set.seed(setup[["seed"]])
pda.derive(brazil.trt.cc, "brazil", getwd(), m = 1, complete = TRUE, outcome = "trt")
# [Brazil] upload brazil_derive.json to shared folder; inform Chile to download

## Step 3b [Mexico]: download updated control.json; derive; upload mexico_derive.json
set.seed(setup[["seed"]])
pda.derive(mexico.trt.cc, "mexico", getwd(), m = 1, complete = TRUE, outcome = "trt")
# [Mexico] upload mexico_derive.json to shared folder; inform Chile to download

## Step 3c [Chile]: download brazil_derive.json and mexico_derive.json; derive
set.seed(setup[["seed"]])
pda.derive(chile.trt.cc, "chile", getwd(), m = 1, complete = TRUE, outcome = "trt")
# [Chile] move to estimation

# Step 4 [Chile -- lead site]: estimate
set.seed(setup[["seed"]])
pda.estimate(chile.trt.cc, "chile", getwd(), m = 1, complete = TRUE, outcome = "trt")

# Step 5 [Chile -- lead site]: extract results and save
trt.dist.cc = pda.summary("chile", c("brazil", "mexico"), getwd(), 1, complete = TRUE, outcome = "trt")
saveRDS(trt.dist.cc, file = "trt_dist_cc.rds")
trt.meta.cc = meta.mi(dist.meta(c("chile", "brazil", "mexico"), getwd(), 1, complete = TRUE, outcome = "trt"))
saveRDS(trt.meta.cc, file = "trt_meta_cc.rds")


# ======================================================================================== #
#                   MAJOR REGIMEN CHANGE -- MULTIPLE IMPUTATION                           #
# ======================================================================================== #

# Step 1 [Chile -- lead site]: generate m control.json files
control <- list(project_name  = 'major_regimen_change_distributed_survival_analysis_ccasanet',
                step          = 'initialize',
                sites         = c("chile", "brazil", "mexico"),
                heterogeneity = TRUE,
                model         = 'ODAC',
                family        = 'cox',
                outcome       = "Surv(ARTchange_status.years, ARTchange_status==\"major treatment change\")",
                variables     = c('sex_mode', 'age', 'clinicalAIDS', 'cd4_v_initial_sqrt',
                                  'rna_v_initial_log', 'calendaryear', 'initialARTtype'),
                xlev          = list(sex_mode     = c('Female','Male-MSM', 'Male-Other'),
                                     clinicalAIDS = c('Not clinical AIDS','Clinical AIDS'),
                                     initialARTtype = c('NNRTI','PI','INSTI')),
                optim_maxit   = 10000,
                init_method   = "meta",
                optim_method  = "Nelder-Mead",
                lead_site     = "chile",
                upload_date   = as.character(Sys.time()))
set.seed(setup[["seed"]])
pda.setup(leadsite = "chile", control, dir = getwd(), m = m, outcome = "trt")
# [Chile] upload m control.json files to shared folder (one sub-folder per imputed dataset);
#         inform Brazil and Mexico to download

# Step 2 [Initialization]: Brazil and Mexico run first, then Chile
## Step 2a [Brazil]: download control.json files; initialize for each imputed dataset;
##                   upload brazil_initialize.json for each imputed dataset
set.seed(setup[["seed"]])
pda.local(brazil.trt.mi, "brazil", getwd(), m = m, outcome = "trt")
# [Brazil] upload brazil_initialize.json (x m) to shared folder; inform Chile to download

## Step 2b [Mexico]: download control.json files; initialize for each imputed dataset;
##                   upload mexico_initialize.json for each imputed dataset
set.seed(setup[["seed"]])
pda.local(mexico.trt.mi, "mexico", getwd(), m = m, outcome = "trt")
# [Mexico] upload mexico_initialize.json (x m) to shared folder; inform Chile to download

## Step 2c [Chile]: download brazil/mexico _initialize.json (x m); initialize;
##                  upload updated control.json (x m); inform Brazil and Mexico to download
set.seed(setup[["seed"]])
pda.local(chile.trt.mi, "chile", getwd(), m = m, outcome = "trt")
# [Chile] upload updated control.json (x m) to shared folder; inform Brazil and Mexico to download

# Step 3 [Derivation]: Brazil and Mexico run first, then Chile
## Step 3a [Brazil]: download updated control.json files; derive for each imputed dataset;
##                   upload brazil_derive.json for each imputed dataset
set.seed(setup[["seed"]])
pda.derive(brazil.trt.mi, "brazil", getwd(), m = m, outcome = "trt")
# [Brazil] upload brazil_derive.json (x m) to shared folder; inform Chile to download

## Step 3b [Mexico]: download updated control.json files; derive for each imputed dataset;
##                   upload mexico_derive.json for each imputed dataset
set.seed(setup[["seed"]])
pda.derive(mexico.trt.mi, "mexico", getwd(), m = m, outcome = "trt")
# [Mexico] upload mexico_derive.json (x m) to shared folder; inform Chile to download

## Step 3c [Chile]: download brazil/mexico _derive.json (x m); derive
set.seed(setup[["seed"]])
pda.derive(chile.trt.mi, "chile", getwd(), m = m, outcome = "trt")
# [Chile] move to estimation

# Step 4 [Chile -- lead site]: estimate for each imputed dataset
set.seed(setup[["seed"]])
pda.estimate(chile.trt.mi, "chile", getwd(), m = m, outcome = "trt")

# Step 5 [Chile -- lead site]: pool results across imputed datasets (Rubin's rules) and save
trt.dist.mi = pda.summary("chile", c("brazil", "mexico"), getwd(), m = m, outcome = "trt")
saveRDS(trt.dist.mi, file = "trt_dist_mi.rds")
trt.meta.mi = meta.mi(dist.meta(c("chile", "brazil", "mexico"), getwd(), m = m, outcome = "trt"))
saveRDS(trt.meta.mi, file = "trt_meta_mi.rds")


# ======================================================================================== #
#                CUMULATIVE INCIDENCE FUNCTION -- FEDERATED ESTIMATION                    #
# ======================================================================================== #
# The Aalen-Johansen estimator for the CIF is computed locally at each site and
# aggregated at the lead site. pooled.time is shared by the lead site to align
# event times across sites before computing the overall federated CIF.

source("./main/federated-cumulative-incidence-function.R")

## Step 1 [All sites]: compute local CIF and extract local event times; upload to shared folder
### Brazil
brazil.trt.time = time_CIF(brazil.trt)
brazil.trt.cif  = obtain_CIF(obtain_counts_CIF(data = brazil.trt, outcome = "trt"))
brazil.vl.time  = time_CIF(brazil.vl)
brazil.vl.cif   = obtain_CIF(obtain_counts_CIF(data = brazil.vl, outcome = "vl"))
brazil = list(trt.time = brazil.trt.time,
              trt.cif  = brazil.trt.cif,
              vl.time  = brazil.vl.time,
              vl.cif   = brazil.vl.cif)
save(brazil, file = "brazil-federated-cumulative-incidence-function.RData")
# [Brazil] upload brazil-federated-cumulative-incidence-function.RData to shared folder

### Mexico
mexico.trt.time = time_CIF(mexico.trt)
mexico.trt.cif  = obtain_CIF(obtain_counts_CIF(data = mexico.trt, outcome = "trt"))
mexico.vl.time  = time_CIF(mexico.vl)
mexico.vl.cif   = obtain_CIF(obtain_counts_CIF(data = mexico.vl, outcome = "vl"))
mexico = list(trt.time = mexico.trt.time,
              trt.cif  = mexico.trt.cif,
              vl.time  = mexico.vl.time,
              vl.cif   = mexico.vl.cif)
save(mexico, file = "mexico-federated-cumulative-incidence-function.RData")
# [Mexico] upload mexico-federated-cumulative-incidence-function.RData to shared folder

### Chile
chile.trt.time = time_CIF(chile.trt)
chile.trt.cif  = obtain_CIF(obtain_counts_CIF(data = chile.trt, outcome = "trt"))
chile.vl.time  = time_CIF(chile.vl)
chile.vl.cif   = obtain_CIF(obtain_counts_CIF(data = chile.vl, outcome = "vl"))
chile = list(trt.time = chile.trt.time,
             trt.cif  = chile.trt.cif,
             vl.time  = chile.vl.time,
             vl.cif   = chile.vl.cif)
save(chile, file = "chile-federated-cumulative-incidence-function.RData")
# [Chile] upload chile-federated-cumulative-incidence-function.RData to shared folder

## Step 2 [Chile -- lead site]: download all sites' RData; merge local event times;
##                               broadcast pooled.time to shared folder
pooled.time = list(
  trt = merge_time_CIF(list(chile.trt.time, brazil.trt.time, mexico.trt.time)),
  vl  = merge_time_CIF(list(chile.vl.time,  brazil.vl.time,  mexico.vl.time))
)
save(pooled.time, file = "pooled-time.RData")
# [Chile] upload pooled-time.RData to shared folder; inform Brazil and Mexico to download

## Step 3 [All sites]: download pooled-time.RData; compute pooled-time CIF counts; upload
load("pooled-time.RData")

### Brazil
brazil.trt.pooled = obtain_counts_CIF(data = brazil.trt, outcome = "trt",
                                       pooled = TRUE, event.time = pooled.time$trt)
brazil.vl.pooled  = obtain_counts_CIF(data = brazil.vl,  outcome = "vl",
                                       pooled = TRUE, event.time = pooled.time$vl)
brazil.pooled = list(trt = brazil.trt.pooled, vl = brazil.vl.pooled)
save(brazil.pooled, file = "brazil-pooled-cumulative-incidence-function.RData")
# [Brazil] upload brazil-pooled-cumulative-incidence-function.RData to shared folder

### Mexico
mexico.trt.pooled = obtain_counts_CIF(data = mexico.trt, outcome = "trt",
                                       pooled = TRUE, event.time = pooled.time$trt)
mexico.vl.pooled  = obtain_counts_CIF(data = mexico.vl,  outcome = "vl",
                                       pooled = TRUE, event.time = pooled.time$vl)
mexico.pooled = list(trt = mexico.trt.pooled, vl = mexico.vl.pooled)
save(mexico.pooled, file = "mexico-pooled-cumulative-incidence-function.RData")
# [Mexico] upload mexico-pooled-cumulative-incidence-function.RData to shared folder

### Chile
chile.trt.pooled = obtain_counts_CIF(data = chile.trt, outcome = "trt",
                                      pooled = TRUE, event.time = pooled.time$trt)
chile.vl.pooled  = obtain_counts_CIF(data = chile.vl,  outcome = "vl",
                                      pooled = TRUE, event.time = pooled.time$vl)
chile.pooled = list(trt = chile.trt.pooled, vl = chile.vl.pooled)
save(chile.pooled, file = "chile-pooled-cumulative-incidence-function.RData")

## Step 4 [Chile -- lead site]: download all sites' pooled RData;
##                               compute overall federated CIF from aggregated counts
trt.pooled.cif = obtain_CIF(list(brazil.pooled$trt, chile.pooled$trt, mexico.pooled$trt),
                              pooled = TRUE)
vl.pooled.cif  = obtain_CIF(list(brazil.pooled$vl,  chile.pooled$vl,  mexico.pooled$vl),
                              pooled = TRUE)
