#------------------------------------------------------------------------------------------#
#     Simulated Federated Analysis for CCASAnet Federated Survival Analysis Study         #
#------------------------------------------------------------------------------------------#
# This script performs a simulated federated analysis in which all participating sites'
# data are available locally. The ODACH algorithm (via the pda package) is applied to
# replicate the communication protocol of a genuine federated analysis without actual
# data transfer across institutions.
#
# Analyses focus on 3 sites: Brazil, Chile, Mexico
# Outcomes: (1) virologic failure; (2) major regimen change
#
# Prerequisites: run 01_centralized-meta-analysis.R to prepare the following objects:
#   trt.95.3sites, vl.95.3sites             (complete-case analysis datasets)
#   trt.95.mi.3sites.split.mi.d             (list of m imputed datasets, major regimen change)
#   vl.95.mi.3sites.split.mi.d              (list of m imputed datasets, virologic failure)
#   trt.95.3sites.fit, vl.95.3sites.fit     (centralized Cox fits, for labeling)
#   m                                       (number of imputations)


# --- Load libraries and source files --- #

source("./main/loadlibrary-ccasanet.R")
source("./main/distributed-analysis-ccasanet.R")

mydir = getwd()


# --- Complete-case simulated federated analysis --- #

## major regimen change
### Chile as lead site
trt.95.3sites.chile.fit.dist = regimenchange.cscox.distributed(
  data      = list(trt.95.3sites),
  leadsite  = "chile",
  mydir     = mydir,
  method    = "Nelder-Mead",
  merge     = TRUE
)
trt.95.3sites.chile.fit.dist[[1]] %>%
  as.data.frame() %>%
  select(all_of(c("hr", "lower.95", "upper.95", "std.error", "p.value"))) %>%
  mutate(name=names(coefficients(trt.95.3sites.fit))) %>%
  column_to_rownames(var = "name")

## virologic failure
### Chile as lead site
vl.95.3sites.chile.fit.dist = virologicfailure.cscox.distributed(
  data      = list(vl.95.3sites),
  leadsite  = "chile",
  mydir     = mydir,
  method    = "Nelder-Mead",
  merge     = TRUE
)
vl.95.3sites.chile.fit.dist[[1]] %>%
  as.data.frame() %>%
  select(all_of(c("hr", "lower.95", "upper.95", "std.error", "p.value"))) %>%
  mutate(name=names(coefficients(vl.95.3sites.fit))) %>%
  column_to_rownames(var = "name")

# --- Multiple imputation simulated federated analysis --- #
# Local imputation is performed at each site (see 01_centralized-meta-analysis.R).
# The ODACH algorithm is applied to each imputed dataset; results are pooled via Rubin's rules.

## major regimen change (Chile as lead site)
trt.95.mi.3sites.split.mi.fit.dist = regimenchange.cscox.distributed(
  data      = trt.95.mi.3sites.split.mi.d,
  leadsite  = "chile",
  mydir     = mydir,
  method    = "Nelder-Mead",
  merge     = TRUE
)
mi.pool.rubin(trt.95.mi.3sites.split.mi.fit.dist,
              nrow(trt.95.mi.3sites.split.mi.d[[1]])) %>%
  as.data.frame() %>%
  mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
         upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
  select(all_of(c("hr", "lower.95", "upper.95", "std.error", "p.value"))) %>%
  mutate(name=names(coefficients(trt.95.3sites.fit))) %>%
  column_to_rownames(var = "name")

## virologic failure (Chile as lead site)
vl.95.mi.3sites.split.mi.fit.dist = virologicfailure.cscox.distributed(
  data      = vl.95.mi.3sites.split.mi.d,
  leadsite  = "chile",
  mydir     = mydir,
  method    = "Nelder-Mead",
  merge     = TRUE
)
mi.pool.rubin(vl.95.mi.3sites.split.mi.fit.dist,
              nrow(vl.95.mi.3sites.split.mi.d[[1]])) %>%
  as.data.frame() %>%
  mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
         upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
  select(all_of(c("hr", "lower.95", "upper.95", "std.error", "p.value"))) %>%
  mutate(name=names(coefficients(vl.95.3sites.fit))) %>%
  column_to_rownames(var = "name")


# --- Wasserstein distance: comparison to centralized analysis --- #
# Wasserstein distance measures proximity between methods in terms of
# point estimate and standard error jointly.

## complete case: major regimen change
trt.95.3sites.fit.mod = data.frame(
  estimate  = trt.95.3sites.fit$coefficients,
  std.error = sqrt(diag(trt.95.3sites.fit$var))
)
data.frame(
  meta             = cal.wd(trt.95.3sites.fit.mod, trt.95.3sites.fit.meta),
  distributed.chile = cal.wd(trt.95.3sites.fit.mod, trt.95.3sites.chile.fit.dist[[1]])
) %>%
  mutate(name=names(coefficients(trt.95.3sites.fit))) %>%
  column_to_rownames(var = "name")

## complete case: virologic failure
vl.95.3sites.fit.mod = data.frame(
  estimate  = vl.95.3sites.fit$coefficients,
  std.error = sqrt(diag(vl.95.3sites.fit$var))
)
data.frame(
  meta               = cal.wd(vl.95.3sites.fit.mod, vl.95.3sites.fit.meta),
  distributed.chile  = cal.wd(vl.95.3sites.fit.mod, vl.95.3sites.chile.fit.dist[[1]])
) %>%
  mutate(name=names(coefficients(vl.95.3sites.fit))) %>%
  column_to_rownames(var = "name")

## multiple imputation: major regimen change
data.frame(
  meta              = cal.wd(trt.95.mi.3sites.split.mi.fit,
                             trt.95.mi.3sites.split.mi.fit.meta),
  distributed.chile = cal.wd(trt.95.mi.3sites.split.mi.fit,
                             mi.pool.rubin(trt.95.mi.3sites.split.mi.fit.dist,
                                          nrow(trt.95.mi.3sites.split.mi.d[[1]])) %>%
                               as.data.frame())
) %>%
  mutate(name=names(coefficients(trt.95.3sites.fit))) %>%
  column_to_rownames(var = "name")

## multiple imputation: virologic failure
data.frame(
  meta              = cal.wd(vl.95.mi.3sites.split.mi.fit,
                             vl.95.mi.3sites.split.mi.fit.meta),
  distributed.chile = cal.wd(vl.95.mi.3sites.split.mi.fit,
                             mi.pool.rubin(vl.95.mi.3sites.split.mi.fit.dist,
                                          nrow(vl.95.mi.3sites.split.mi.d[[1]])) %>%
                               as.data.frame())
) %>%
  mutate(name=names(coefficients(vl.95.3sites.fit))) %>%
  column_to_rownames(var = "name")
