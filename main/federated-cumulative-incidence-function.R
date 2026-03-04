###################################################
# Plot Cumulative Incidence Function/Curve (CIF)
# in a federated manner
# Kaixing Liu, kaixing.liu@vanderbilt.edu
###################################################

# require(dplyr)
# require(stringr)
# require(survival)
# require(rms)
# require(Hmisc)
packages <- c("dplyr", "stringr", "survival", "rms", "Hmisc")

# install any that are not already installed
to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install)

# load all
lapply(packages, library, character.only = TRUE)

# obtain Aalen–Johansen Estimator for the CIF
## at each site
## and overall

obtain_counts_CIF = function(data, outcome, pooled = FALSE,
                      event.time = NULL){
  if(outcome=="trt"){
    status.list = c("right censored", "death", "major treatment change")
  }
  else if(outcome=="vl"){
    status.list = c("right censored", "death", "virologic failure")
  }
  else{
    stop("Outcome is either 'major regimen change' or 'virologic failure'.")
  }
  status = data[, colnames(data)[grepl("_status$",colnames(data))]]
  time = data[, colnames(data)[grepl(".years$",colnames(data))]]
  if(!pooled){
    event.time = sort(unique(time))
    tmp = sapply(event.time, function(t) {
      n.t = which(time>=t) %>% length();
      m.t = which(time==t & status==status.list[3]) %>% length();
      d.t = which(time==t & status==status.list[2]) %>% length();
      c(index.time = which(event.time==t),
        event.time = t,
        event.no = m.t,
        death.no = d.t,
        atrisk.no = n.t)
    }) %>% t() %>% as.data.frame()
  }
  else{
    if(!is.null(event.time)){
      tmp = sapply(event.time, function(t) {
        n.t = which(time>=t) %>% length();
        m.t = which(time==t & status==status.list[3]) %>% length();
        d.t = which(time==t & status==status.list[2]) %>% length();
        c(index.time = which(event.time==t),
          event.time = t,
          event.no = m.t,
          death.no = d.t,
          atrisk.no = n.t)
      }) %>% t()
    }
    else{
      stop("'event.time' must be provided for overall cumulative incidence function.")
    }
  }
  return(tmp)
}

obtain_CIF = function(counts_CIF, pooled=FALSE){
  tmp = counts_CIF
  if(pooled){
    header = tmp[[1]][,1:2]
    tmp = cbind(header, Reduce("+", 
                            sapply(tmp, function(x) x[,-c(1:2)],
                                   simplify = F))) %>% as.data.frame()
  }
  tmp$hazard = tmp$event.no/tmp$atrisk.no
  tmp$surv.prev = sapply(1:nrow(tmp), function(x) {
    if(x==1){
      1;
    }
    else{
      prod(sapply(1:(x-1), function(t) {
        1 - (tmp$event.no[t] +
               tmp$death.no[t])/tmp$atrisk.no[t]
      }, simplify = T))
    }
  })
  tmp$inc = tmp$hazard*tmp$surv.prev
  tmp$cuminc = sapply(1:nrow(tmp), function(x) {
    sum(tmp$inc[1:x])
  })
  return(tmp)
}

# obtain time at each site
time_CIF = function(data){
  time = data[, colnames(data)[grepl(".years$",colnames(data))]]
  return(sort(unique(time)))
}

# merge/pool time at lead site/central server
merge_time_CIF = function(timelist){
  return(sort(unique(unlist(timelist))))
}
