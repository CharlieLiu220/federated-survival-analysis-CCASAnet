#------------------------------------------------------------------------------------------#
#       Figures and Tables for CCASAnet Federated Survival Analysis Manuscript            #
#------------------------------------------------------------------------------------------#
# This script reproduces the main figures and tables reported in the manuscript.
#
# Prerequisites:
#   - Run 01_centralized-meta-analysis.R to produce:
#       trt.95.3sites, vl.95.3sites                (complete-case data)
#       trt.95.3sites.fit, vl.95.3sites.fit         (centralized Cox fits)
#       trt.95.3sites.fit.meta, vl.95.3sites.fit.meta  (complete-case meta-analysis)
#       trt.95.mi.3sites.split.mi.fit               (MI centralized, major regimen change)
#       vl.95.mi.3sites.split.mi.fit                (MI centralized, virologic failure)
#       trt.95.mi.3sites.split.mi.fit.meta          (MI meta-analysis, major regimen change)
#       vl.95.mi.3sites.split.mi.fit.meta           (MI meta-analysis, virologic failure)
#       trt.95.mi.3sites.split.mi.fit.dist          (MI simulated federated, major regimen change)
#       vl.95.mi.3sites.split.mi.fit.dist           (MI simulated federated, virologic failure)
#       trt.95.mi.3sites.split.mi.d                 (list of m imputed datasets, major regimen change)
#       vl.95.mi.3sites.split.mi.d                  (list of m imputed datasets, virologic failure)
#   - Run 02_simulated-federated-analysis.R to produce:
#       trt.95.3sites.chile.fit.dist, trt.95.3sites.brazil.fit.dist  (complete-case, dist.)
#       vl.95.3sites.chile.fit.dist, vl.95.3sites.brazil.fit.dist, vl.95.3sites.mexico.fit.dist
#   - Run 03_authentic-federated-analysis.R and have available (saved as .rds):
#       trt_dist_cc, trt_dist_mi   (authentic federated, major regimen change)
#       vl_dist_cc,  vl_dist_mi    (authentic federated, virologic failure)
#       trt.pooled.cif, vl.pooled.cif      (federated overall CIF)
#       chile, brazil, mexico              (local CIF objects from each site)


# --- Load libraries and source files --- #

source("./main/loadlibrary-ccasanet.R")
source("./main/distributed-analysis-ccasanet.R")
source("./main/federated-cumulative-incidence-function.R")

# load authentic federated analysis results (saved by 03_authentic-federated-analysis.R)
trt_dist_cc = readRDS("trt_dist_cc.rds")
trt_dist_mi = readRDS("trt_dist_mi.rds")
vl_dist_cc  = readRDS("vl_dist_cc.rds")
vl_dist_mi  = readRDS("vl_dist_mi.rds")


# ======================================================================================== #
#               FIGURE: Cumulative Incidence Function (CIF)                               #
# ======================================================================================== #
# Overlay of: pooled (centralized) CIF, federated overall CIF, and site-specific CIFs.
# Dotted black lines: Aalen-Johansen estimator from the centralized (pooled) data.
# Solid colored lines: federated Aalen-Johansen estimator computed without data sharing.

## color palette: Chile = #fdc086, Brazil = #beaed4, Mexico = #386cb0
color3 = c("#fdc086", "#beaed4", "#386cb0")

## CIF from centralized analysis (for comparison)
### major regimen change
trt.pooled.data = rbind(trt.95.3sites %>% subset(site=="chile"),
                        trt.95.3sites %>% subset(site=="brazil"),
                        trt.95.3sites %>% subset(site=="mexico")) %>%
  mutate(site = factor(site, levels = c("chile", "brazil", "mexico")))
trt_cif_overall = npsurv(Surv(time = ARTchange_status.years, event = ARTchange_status) ~ 1,
                          data = trt.pooled.data)
trt_cif_local   = npsurv(Surv(time = ARTchange_status.years, event = ARTchange_status) ~ site,
                          data = trt.pooled.data)
### virologic failure
vl.pooled.data = rbind(vl.95.3sites %>% subset(site=="chile"),
                       vl.95.3sites %>% subset(site=="brazil"),
                       vl.95.3sites %>% subset(site=="mexico")) %>%
  mutate(site = factor(site, levels = c("chile", "brazil", "mexico")))
vl_cif_overall = npsurv(Surv(time = vlfailure_status.years, event = vlfailure_status) ~ 1,
                         data = vl.pooled.data)
vl_cif_local   = npsurv(Surv(time = vlfailure_status.years, event = vlfailure_status) ~ site,
                         data = vl.pooled.data)

## Figure: CIF for major regimen change
par(xaxs = "i")
survplot(trt_cif_overall, state = "major treatment change", conf = "none",
         label.curves = F,
         lty = 5,
         col = "black",
         dots = TRUE,
         xlab = expression(bold("Time from ART initiation (years)")),
         ylab = expression(bold("Cumulative incidence of major treatment change")),
         n.risk = T,
         y.n.risk = -0.2, adj.n.risk = 0.5,
         sep.n.risk = 1, cex.n.risk = 0.8,
         xlim = c(0,9), ylim = c(0, 0.8),
         cex.xlab = 1, cex.ylab = 1)
survplot(trt_cif_local, state = "major treatment change", conf = "none",
         label.curves = F,
         lty = 5,
         col = "black",
         dots = TRUE,
         add = T)
lines(trt.pooled.cif$event.time, trt.pooled.cif$cuminc,
      type = "s", col = "#7fc97f", lty = 4, lwd = 2)
lines(chile$trt.cif$event.time, chile$trt.cif$cuminc,
      type = "s", col = "#fdc086", lty = 4, lwd = 2)
lines(brazil$trt.cif$event.time, brazil$trt.cif$cuminc,
      type = "s", col = "#beaed4", lty = 4, lwd = 2)
lines(mexico$trt.cif$event.time, mexico$trt.cif$cuminc,
      type = "s", col = "#386cb0", lty = 4, lwd = 2)
legend("topleft",
       legend = c("pooled", "distributed-overall", "distributed-Chile",
                  "distributed-Brazil", "distributed-Mexico"),
       col = c("black", "#7fc97f", color3),
       lty = c(5, 4, 4, 4, 4),
       lwd = 2,
       bty = "n",
       text.font = 2)
par(xaxs = "r")

## Figure: CIF for virologic failure
par(xaxs = "i")
survplot(vl_cif_overall, state = "virologic failure", conf = "none",
         label.curves = F,
         lty = 5,
         col = "black",
         dots = TRUE,
         xlab = expression(bold("Time from ART initiation (years)")),
         ylab = expression(bold("Cumulative incidence of virologic failure")),
         n.risk = T,
         y.n.risk = -0.1, adj.n.risk = 0.5,
         sep.n.risk = 1, cex.n.risk = 0.8,
         xlim = c(0,9), ylim = c(0, 0.4),
         cex.xlab = 1, cex.ylab = 1)
survplot(vl_cif_local, state = "virologic failure", conf = "none",
         label.curves = F,
         lty = 5,
         col = "black",
         dots = TRUE,
         add = T)
lines(vl.pooled.cif$event.time, vl.pooled.cif$cuminc,
      type = "s", col = "#7fc97f", lty = 4, lwd = 2)
lines(chile$vl.cif$event.time, chile$vl.cif$cuminc,
      type = "s", col = "#fdc086", lty = 4, lwd = 2)
lines(brazil$vl.cif$event.time, brazil$vl.cif$cuminc,
      type = "s", col = "#beaed4", lty = 4, lwd = 2)
lines(mexico$vl.cif$event.time, mexico$vl.cif$cuminc,
      type = "s", col = "#386cb0", lty = 4, lwd = 2)
legend("topleft",
       legend = c("pooled", "distributed-overall", "distributed-Chile",
                  "distributed-Brazil", "distributed-Mexico"),
       col = c("black", "#7fc97f", color3),
       lty = c(5, 4, 4, 4, 4),
       lwd = 2,
       bty = "n",
       text.font = 2)
par(xaxs = "r")


# ======================================================================================== #
#          FIGURE: Forest plot -- complete case, major regimen change                     #
# ======================================================================================== #

df.trt = rbind(
  data.frame(
    hr       = summary(trt.95.3sites.fit)$conf.int[,"exp(coef)"],
    lower.95 = summary(trt.95.3sites.fit)$conf.int[,"lower .95"],
    upper.95 = summary(trt.95.3sites.fit)$conf.int[,"upper .95"],
    p.value  = summary(trt.95.3sites.fit)$coefficients[, "Pr(>|z|)"]
  ),
  trt.95.3sites.chile.fit.dist[[1]] %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  trt_dist_cc$results[[1]]          %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  trt.95.3sites.brazil.fit.dist[[1]] %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  trt.95.3sites.fit.meta             %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value")))
)
df.trt$coefficient = rep(names(coefficients(trt.95.3sites.fit)), times = 5) %>%
  factor(levels = names(coefficients(trt.95.3sites.fit)),
         labels = c("male-MSM", "male-Other", "age", "clinical AIDS",
                    "CD4 count", "viral load", "calenday year", "PI", "INSTI"))
df.trt$method = rep(c("Centralized", "Distributed-Chile", "Distributed-Chile-Authentic",
                      "Distributed-Brazil", "Meta-analysis"),
                    each = length(coefficients(trt.95.3sites.fit)))
rownames(df.trt) = NULL
plot.forest(df = df.trt, outcome = "Major Regimen Change")


# ======================================================================================== #
#          FIGURE: Forest plot -- complete case, virologic failure                        #
# ======================================================================================== #

df.vl = rbind(
  data.frame(
    hr       = summary(vl.95.3sites.fit)$conf.int[,"exp(coef)"],
    lower.95 = summary(vl.95.3sites.fit)$conf.int[,"lower .95"],
    upper.95 = summary(vl.95.3sites.fit)$conf.int[,"upper .95"],
    p.value  = summary(vl.95.3sites.fit)$coefficients[, "Pr(>|z|)"]
  ),
  vl.95.3sites.chile.fit.dist[[1]]  %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  vl_dist_cc$results[[1]]           %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  vl.95.3sites.brazil.fit.dist[[1]] %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  vl.95.3sites.mexico.fit.dist[[1]] %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  vl.95.3sites.fit.meta             %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value")))
)
df.vl$coefficient = rep(names(coefficients(vl.95.3sites.fit)), times = 6) %>%
  factor(levels = names(coefficients(vl.95.3sites.fit)),
         labels = c("male-MSM", "male-Other", "age", "clinical AIDS",
                    "CD4 count", "viral load", "calenday year", "PI", "INSTI"))
df.vl$method = rep(c("Centralized", "Distributed-Chile", "Distributed-Chile-Authentic",
                     "Distributed-Brazil", "Distributed-Mexico", "Meta-analysis"),
                   each = length(coefficients(vl.95.3sites.fit)))
rownames(df.vl) = NULL
plot.forest(df = df.vl, outcome = "Virologic Failure")


# ======================================================================================== #
#       FIGURE: Forest plot -- multiple imputation, major regimen change                  #
# ======================================================================================== #

df.trt.mi = rbind(
  trt.95.mi.3sites.split.mi.fit %>%
    mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
           upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
    select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  mi.pool.rubin(trt.95.mi.3sites.split.mi.fit.dist,
                nrow(trt.95.mi.3sites.split.mi.d[[1]])) %>%
    mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
           upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
    select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  mi.pool.rubin(trt_dist_mi$results, sum(trt_dist_mi$site_size)) %>%
    mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
           upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
    select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  trt.95.mi.3sites.split.mi.fit.meta %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value")))
)
df.trt.mi$coefficient = rep(names(coefficients(trt.95.3sites.fit)), times = 4) %>%
  factor(levels = names(coefficients(trt.95.3sites.fit)),
         labels = c("male-MSM", "male-Other", "age", "clinical AIDS",
                    "CD4 count", "viral load", "calenday year", "PI", "INSTI"))
df.trt.mi$method = rep(c("Centralized", "Distributed-Chile", "Distributed-Chile-Authentic",
                         "Meta-analysis"),
                       each = length(coefficients(trt.95.3sites.fit)))
rownames(df.trt.mi) = NULL
plot.forest(df = df.trt.mi, outcome = "Major Regimen Change")


# ======================================================================================== #
#       FIGURE: Forest plot -- multiple imputation, virologic failure                     #
# ======================================================================================== #

df.vl.mi = rbind(
  vl.95.mi.3sites.split.mi.fit %>%
    mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
           upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
    select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  mi.pool.rubin(vl.95.mi.3sites.split.mi.fit.dist,
                nrow(vl.95.mi.3sites.split.mi.d[[1]])) %>%
    mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
           upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
    select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  mi.pool.rubin(vl_dist_mi$results, sum(vl_dist_mi$site_size)) %>%
    mutate(lower.95 = exp(estimate - qnorm(.975)*std.error),
           upper.95 = exp(estimate + qnorm(.975)*std.error)) %>%
    select(all_of(c("hr", "lower.95", "upper.95", "p.value"))),
  vl.95.mi.3sites.split.mi.fit.meta %>% select(all_of(c("hr", "lower.95", "upper.95", "p.value")))
)
df.vl.mi$coefficient = rep(names(coefficients(vl.95.3sites.fit)), times = 4) %>%
  factor(levels = names(coefficients(trt.95.3sites.fit)),
         labels = c("male-MSM", "male-Other", "age", "clinical AIDS",
                    "CD4 count", "viral load", "calenday year", "PI", "INSTI"))
df.vl.mi$method = rep(c("Centralized", "Distributed-Chile", "Distributed-Chile-Authentic",
                        "Meta-analysis"),
                      each = length(coefficients(vl.95.3sites.fit)))
rownames(df.vl.mi) = NULL
plot.forest(df = df.vl.mi, outcome = "Virologic Failure")


# ======================================================================================== #
#        TABLE: Wasserstein distance -- complete case                                     #
# ======================================================================================== #

## major regimen change
trt.95.3sites.fit.mod = data.frame(
  estimate  = trt.95.3sites.fit$coefficients,
  std.error = sqrt(diag(trt.95.3sites.fit$var))
)
data.frame(
  meta                     = cal.wd(trt.95.3sites.fit.mod, trt.95.3sites.fit.meta),
  distributed.chile        = cal.wd(trt.95.3sites.fit.mod, trt.95.3sites.chile.fit.dist[[1]]),
  distributed.chile.auth   = cal.wd(trt.95.3sites.fit.mod, trt_dist_cc$results[[1]]),
  distributed.brazil       = cal.wd(trt.95.3sites.fit.mod, trt.95.3sites.brazil.fit.dist[[1]])
) %>%
  mutate(name=names(coefficients(trt.95.3sites.fit))) %>%
  column_to_rownames(var = "name")

## virologic failure
vl.95.3sites.fit.mod = data.frame(
  estimate  = vl.95.3sites.fit$coefficients,
  std.error = sqrt(diag(vl.95.3sites.fit$var))
)
data.frame(
  meta                   = cal.wd(vl.95.3sites.fit.mod, vl.95.3sites.fit.meta),
  distributed.chile      = cal.wd(vl.95.3sites.fit.mod, vl.95.3sites.chile.fit.dist[[1]]),
  distributed.chile.auth = cal.wd(vl.95.3sites.fit.mod, vl_dist_cc$results[[1]]),
  distributed.brazil     = cal.wd(vl.95.3sites.fit.mod, vl.95.3sites.brazil.fit.dist[[1]]),
  distributed.mexico     = cal.wd(vl.95.3sites.fit.mod, vl.95.3sites.mexico.fit.dist[[1]])
) %>%
  mutate(name=names(coefficients(vl.95.3sites.fit))) %>%
  column_to_rownames(var = "name")


# ======================================================================================== #
#        TABLE: Wasserstein distance -- multiple imputation                               #
# ======================================================================================== #

## major regimen change
data.frame(
  meta                   = cal.wd(trt.95.mi.3sites.split.mi.fit,
                                  trt.95.mi.3sites.split.mi.fit.meta),
  distributed.chile      = cal.wd(trt.95.mi.3sites.split.mi.fit,
                                  mi.pool.rubin(trt.95.mi.3sites.split.mi.fit.dist,
                                               nrow(trt.95.mi.3sites.split.mi.d[[1]])) %>%
                                    as.data.frame()),
  distributed.chile.auth = cal.wd(trt.95.mi.3sites.split.mi.fit,
                                  mi.pool.rubin(trt_dist_mi$results, sum(trt_dist_mi$site_size)))
) %>%
  mutate(name=names(coefficients(trt.95.3sites.fit))) %>%
  column_to_rownames(var = "name")

## virologic failure
data.frame(
  meta                   = cal.wd(vl.95.mi.3sites.split.mi.fit,
                                  vl.95.mi.3sites.split.mi.fit.meta),
  distributed.chile      = cal.wd(vl.95.mi.3sites.split.mi.fit,
                                  mi.pool.rubin(vl.95.mi.3sites.split.mi.fit.dist,
                                               nrow(vl.95.mi.3sites.split.mi.d[[1]])) %>%
                                    as.data.frame()),
  distributed.chile.auth = cal.wd(vl.95.mi.3sites.split.mi.fit,
                                  mi.pool.rubin(vl_dist_mi$results, sum(vl_dist_mi$site_size)))
) %>%
  mutate(name=names(coefficients(vl.95.3sites.fit))) %>%
  column_to_rownames(var = "name")
