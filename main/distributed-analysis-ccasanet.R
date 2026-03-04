#---------------------------------------------------------------------------------#
#      Functions defined for distributed survival analysis of ccasanet data       #
#---------------------------------------------------------------------------------#

# load library #
library("rjson")
library(rlang)
library(metafor)

# override menu for non-interactive running
menu <- function(choices, title = NULL) { return(1) } ## this doesn't work tbh

# currently
# only consider 3 sites
# Brazil, Chile, Mexico
# with global or local imputation
# try with different lead sites

# major treatment change #

regimenchange.cscox.distributed = function(data, leadsite, mydir, method=NULL, merge=FALSE){
  results = vector("list", length(data)) # data should be a list of data frames
  for(i in 1:length(data)){
    dat.split = split(data[[i]], data[[i]]$site)
    sites = levels(data[[i]]$site)
    # create new directory
    newdir = paste(mydir, "/regimenchange", sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }
    newdir = paste(newdir, "/", leadsite, sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }    
    if(length(data)==1){
      newdir = paste(newdir, "/complete-case", sep = "")      
    }
    else{
      newdir = paste(newdir, "/", i, sep = "")      
    }
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }
    # set up
    control <- list(project_name = 'Federated_trt_global_mi_3sites',
                    step = 'initialize',
                    sites = sites,
                    heterogeneity = TRUE,
                    model = 'ODAC',
                    family = 'cox',
                    outcome = "Surv(ARTchange_status.years, ARTchange_status==\"major treatment change\")",
                    # xlev = list(sex=c('F', 'M')),  #levels of all categorical X's, with the first being the reference
                    optim_maxit = 10000,
                    init_method = "meta",
                    #optim_method = "Nelder-Mead",
                    #optim_method = "BFGS",
                    lead_site = leadsite,
                    upload_date = as.character(Sys.time()))
    if(merge){
      control[["variables"]] = c('sex_mode', 'age', 'clinicalAIDS', 'cd4_v_initial_sqrt', 'rna_v_initial_log', 'calendaryear', 'initialARTtype')
      control[["xlev"]] =  list(sex_mode=c('Female','Male-MSM', 'Male-Other'), clinicalAIDS=c('Not clinical AIDS','Clinical AIDS'), initialARTtype=c('NNRTI','PI','INSTI'))
    }
    else{
      control[["variables"]] = c('sex', 'age', 'mode_new', 'clinicalAIDS', 'cd4_v_initial_sqrt', 'rna_v_initial_log', 'calendaryear', 'initialARTtype')
      control[["xlev"]] = list(sex=c('Female','Male'), mode_new=c('Sexual','MSM','Other'), clinicalAIDS=c('Not clinical AIDS','Clinical AIDS'), initialARTtype=c('NNRTI','PI','INSTI'))
    }
    if(!is.null(method)){
      control$optim_method = method
    }
    # set up for lead site
    pda(site_id = leadsite, control = control, dir = newdir, upload_without_confirm = T)
    # initialize
    for(site in sites[-which(sites==leadsite)]){
      pda(site_id = site, ipdata = dat.split[[site]], dir=newdir, upload_without_confirm = T)      
    }
    ## initialize for lead site
    pda(site_id = leadsite, ipdata = dat.split[[leadsite]], dir=newdir, upload_without_confirm = T)
    # derivative
    for(site in sites[-which(sites==leadsite)]){
      pda(site_id = site, ipdata = dat.split[[site]], dir=newdir, upload_without_confirm = T)      
    }
    ##
    pda(site_id = leadsite, ipdata = dat.split[[leadsite]], dir=newdir, upload_without_confirm = T)
    # estimate
    pda(site_id = leadsite, ipdata = dat.split[[leadsite]], dir=newdir, upload_without_confirm = T)
    # results
    config <- getCloudConfig(site_id = leadsite, dir=newdir)
    fit.pda <- pdaGet(name = paste0(c(leadsite, "estimate"), collapse = "_"), config = config)
    control <- pdaGet('control', config)
    ##
    results.tmp = cbind(fit.pda$btilde, exp(fit.pda$btilde), sqrt(diag(solve(fit.pda$Htilde)/nrow(data[[i]]))))
    colnames(results.tmp) = c("estimate","hr", "std.error")
    results.tmp = as.data.frame(results.tmp)
    results.tmp$lower.95 = exp(results.tmp$estimate - qnorm(.975) * results.tmp$std.error)
    results.tmp$upper.95 = exp(results.tmp$estimate + qnorm(.975) * results.tmp$std.error)
    results.tmp$p.value = 2*pnorm(abs(fit.pda$btilde)/results.tmp$std.error,
                                     lower.tail = F)
    results[[i]] = results.tmp
  }
  return(results)
}

# virologic failure #

virologicfailure.cscox.distributed = function(data, leadsite, mydir, method=NULL, merge=FALSE){
  results = vector("list", length(data)) # data should be a list of data frames
  for(i in 1:length(data)){
    dat.split = split(data[[i]], data[[i]]$site)
    sites = levels(data[[i]]$site)
    # create new directory
    newdir = paste(mydir, "/virologicfailure", sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }
    newdir = paste(newdir, "/", leadsite, sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }    
    if(length(data)==1){
      newdir = paste(newdir, "/complete-case", sep = "")      
    }
    else{
      newdir = paste(newdir, "/", i, sep = "")      
    }
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }    
    # set up
    control <- list(project_name = 'Federated_vl_global_mi_3sites',
                    step = 'initialize',
                    sites = sites,
                    heterogeneity = TRUE,
                    model = 'ODAC',
                    family = 'cox',
                    outcome = "Surv(vlfailure_status.years, vlfailure_status==\"virologic failure\")",
                    # xlev = list(sex=c('F', 'M')),  #levels of all categorical X's, with the first being the reference
                    optim_maxit = 10000,
                    #optim_method = method,
                    init_method = "meta",
                    lead_site = leadsite,
                    upload_date = as.character(Sys.time()))
    if(merge){
      control[["variables"]] = c('sex_mode', 'age', 'clinicalAIDS', 'cd4_v_initial_sqrt', 'rna_v_initial_log', 'calendaryear', 'initialARTtype')
      control[["xlev"]] =  list(sex_mode=c('Female','Male-MSM', 'Male-Other'), clinicalAIDS=c('Not clinical AIDS','Clinical AIDS'), initialARTtype=c('NNRTI','PI','INSTI'))
    }
    else{
      control[["variables"]] = c('sex', 'age', 'mode_new', 'clinicalAIDS', 'cd4_v_initial_sqrt', 'rna_v_initial_log', 'calendaryear', 'initialARTtype')
      control[["xlev"]] = list(sex=c('Female','Male'), mode_new=c('Sexual','MSM','Other'), clinicalAIDS=c('Not clinical AIDS','Clinical AIDS'), initialARTtype=c('NNRTI','PI','INSTI'))
    }
    if(!is.null(method)){
      control$optim_method = method
    }    
    # set up for lead site
    pda(site_id = leadsite, control = control, dir = newdir, upload_without_confirm = T)
    # initialize
    for(site in sites[-which(sites==leadsite)]){
      pda(site_id = site, ipdata = dat.split[[site]], dir=newdir, upload_without_confirm = T)      
    }
    ## initialize for lead site
    pda(site_id = leadsite, ipdata = dat.split[[leadsite]], dir=newdir, upload_without_confirm = T)
    # derivative
    for(site in sites[-which(sites==leadsite)]){
      pda(site_id = site, ipdata = dat.split[[site]], dir=newdir, upload_without_confirm = T)      
    }
    ##
    pda(site_id = leadsite, ipdata = dat.split[[leadsite]], dir=newdir, upload_without_confirm = T)
    # estimate
    pda(site_id = leadsite, ipdata = dat.split[[leadsite]], dir=newdir, upload_without_confirm = T)
    # results
    config <- getCloudConfig(site_id = leadsite, dir=newdir)
    fit.pda <- pdaGet(name = paste0(c(leadsite, "estimate"), collapse = "_"), config = config)
    control <- pdaGet('control', config)
    ##
    results.tmp = cbind(fit.pda$btilde, exp(fit.pda$btilde), sqrt(diag(solve(fit.pda$Htilde)/nrow(data[[i]]))))
    colnames(results.tmp) = c("estimate","hr", "std.error")
    results.tmp = as.data.frame(results.tmp)
    results.tmp$lower.95 = exp(results.tmp$estimate - qnorm(.975) * results.tmp$std.error)
    results.tmp$upper.95 = exp(results.tmp$estimate + qnorm(.975) * results.tmp$std.error)
    results.tmp$p.value = 2*pnorm(abs(fit.pda$btilde)/results.tmp$std.error,
                                  lower.tail = F)
    results[[i]] = results.tmp
  }
  return(results)
}

# naive meta analysis
## based on cox ph fit in json file from each site
### <https://ppw.kuleuven.be/mesrg/documents/sr2015-course-text/workshop-r/Manual-metafor.pdf>
meta = function(fit.list, random=TRUE){
  p = length(fit.list[[1]]$bhat_i)
  fit.meta = fit.meta.se = c()
  for (i in 1:p){ # for each covariate
    if(random){
      meta_results = rma.uni(yi = sapply(fit.list, function(x) x$bhat_i[i], simplify = T),
                             vi = sapply(fit.list, function(x) x$Vhat_i[i], simplify = T))      
    }
    else{
      meta_results = rma.uni(yi = sapply(fit.list, function(x) x$bhat_i[i], simplify = T),
                             vi = sapply(fit.list, function(x) x$Vhat_i[i], simplify = T),
                             method = "FE")
    }
    fit.meta[i] = meta_results$beta
    fit.meta.se[i] = meta_results$se
  }
  # compute lower & upper bound
  upper.95 = fit.meta + qnorm(.975)*fit.meta.se
  lower.95 = fit.meta - qnorm(.975)*fit.meta.se
  pvalue = 2*pnorm(abs(fit.meta/fit.meta.se), lower.tail = F)
  return(data.frame(estimate=fit.meta,
                    std.error=fit.meta.se,
                    hr=exp(fit.meta),
                    lower.95=exp(lower.95),
                    upper.95=exp(upper.95),
                    p.value = pvalue))
}

meta.mi = function(fit.list, random=TRUE){
  p = nrow(fit.list[[1]])
  fit.meta = fit.meta.se = c()
  for (i in 1:p){ # for each covariate
    if(random){
      meta_results = rma.uni(yi = sapply(fit.list, function(x) x[i, "estimate"], simplify = T),
                             vi = sapply(fit.list, function(x) x[i, "std.error"]^2, simplify = T))      
    }
    else{
      meta_results = rma.uni(yi = sapply(fit.list, function(x) x[i, "estimate"], simplify = T),
                             vi = sapply(fit.list, function(x) x[i, "std.error"]^2, simplify = T),
                             method = "FE")
    }
    fit.meta[i] = meta_results$beta
    fit.meta.se[i] = meta_results$se
  }
  # compute lower & upper bound
  upper.95 = fit.meta + qnorm(.975)*fit.meta.se
  lower.95 = fit.meta - qnorm(.975)*fit.meta.se
  pvalue = 2*pnorm(abs(fit.meta/fit.meta.se), lower.tail = F)
  return(data.frame(estimate=fit.meta,
                    std.error=fit.meta.se,
                    hr=exp(fit.meta),
                    lower.95=exp(lower.95),
                    upper.95=exp(upper.95),
                    p.value = pvalue))
}


# pool results of complete models via multiple imputation using Rubin's rule #
## Multiple Imputation for Nonresponse in Surveys, p76
## <https://bookdown.org/mwheymans/bookmi/rubins-rules.html#degrees-of-freedom-and-p-values> use adjusted degree of freedom

require(dplyr)

mi.pool.rubin = function(results, samplesize){
  # results is a list
  m = length(results)
  # coefficient estimate
  Qbar = results %>% sapply(function(x) x$estimate) %>% rowMeans()
  Ubar = results %>% sapply(function(x) (x$std.error)^2) %>% rowMeans()
  B = results %>% sapply(function(x) (x$estimate-Qbar)%*%t(x$estimate-Qbar), simplify = F)
  B = Reduce("+", B)/(m-1)
  # variance estimate
  Tm = Ubar + (1+(1/m))*diag(B)
  # p-value
  pvalue = sapply(1:length(Qbar), function(x) {
    #rm.tmp = (1+(1/m))*B[x,x]/Ubar[x]
    rm.tmp = (1+(1/m))*B[x,x]/Tm[x] 
    nu.tmp = (m-1)*(1+(1/rm.tmp))^2
    nu.tmp.obs = ((samplesize-length(Qbar)+1)/(samplesize-length(Qbar)+3))*(samplesize-length(Qbar))*(1-(rm.tmp/(1+rm.tmp)))
    nu.tmp = (nu.tmp*nu.tmp.obs)/(nu.tmp+nu.tmp.obs)
    p.tmp = 2*pt(abs(Qbar[x])/sqrt(Tm[x]), df=nu.tmp, lower.tail = F)
    c(df=nu.tmp, pvalue=p.tmp)
  })
  return(data.frame(estimate = Qbar,
                    hr = exp(Qbar),
                    std.error = sqrt(Tm),
                    #df = pvalue["df", ],  
                    p.value = pvalue["pvalue", ]))
}


# calculate Wasserstein distance #
## compare different methods in terms of the recovery of the centralized analysis

cal.wd = function(benchmark, tobecompared){
  return(sqrt((benchmark$estimate-tobecompared$estimate)^2+
                (benchmark$std.error-tobecompared$std.error)^2))
}

# show imputed values for a specific variable among observed complete data #

show.imputation = function(mids, data=NULL, var, overall=TRUE){
  if(class(mids)=="mids"){
    dat = complete(mids, action=0)
    tmp = complete(mids, action="long")
    index = is.na(dat[,var])
    tmp$imputed = rep(index, length(unique(tmp$.imp)))    
  }
  else{
    dat = data
    index = is.na(dat[,var])
    tmp = bind_rows(mids)
    tmp$.imp = rep(1:length(mids), each=nrow(dat))
    tmp$imputed = rep(index, length(mids))
  }
  if(sum(str_detect(colnames(dat), "ARTchange"))!=0){
    outcome = "major regimen change"
  }
  else{
    outcome = "virologic failure"
  }
  if(overall){
    if(is.numeric(dat[,var])){
      ggplot(tmp, aes(x=.imp, y=!!sym(var), color=imputed)) +
        geom_jitter(alpha = 0.5) +
        labs(x="index of imputed dataset",
             y="value",
             title = paste(outcome, var, "overall", sep = ", ")) +
        theme_minimal()
    }
  }
  else{
    if(is.numeric(dat[,var])){
      ggplot(tmp, aes(x=.imp, y=!!sym(var), color=imputed)) +
        geom_jitter(alpha = 0.5) +
        labs(x="index of imputed dataset",
             y="value",
             title = paste(outcome, var, "by site", sep = ", ")) +
        theme_minimal() +
        facet_wrap(~site, ncol = 1)
    }    
  }
}

# forest plot for comparison #

## auxiliary function ##
### for way larger upper endpoints of CI
# piecewise_trans <- trans_new(
#   name = "piecewise",
#   transform = function(x) ifelse(x < 20, x, sqrt(x) + 8),
#   inverse = function(x) ifelse(x < 20, x, (x - 8)^2)
# )

plot.forest = function(df, outcome){
  if(length(unique(df$method))==4 & length(which(unique(df$method)=="Distributed-Brazil")) > 0){
    # custom_colors = c("Centralized" = "#CFAE70", "Meta-analysis" = "#E0D5C0", 
    #                   #"Distributed" = "#FEEEBA"
    #                   "Distributed-Chile" = "#B49248", "Distributed-Brazil" = "#ECB748")
    custom_colors = c("Centralized" = "#7fc97f", "Meta-analysis" = "#beaed4",
                      "Distributed-Chile" = "#fdc086", "Distributed-Brazil" = "#ffff99")
  }
  else if(length(unique(df$method))==4 & length(which(unique(df$method)=="Distributed-Chile-Authentic")) > 0){
    custom_colors = c("Centralized" = "#7fc97f", "Meta-analysis" = "#beaed4",
                      "Distributed-Chile" = "#fdc086", "Distributed-Chile-Authentic" = "#FEEEBA")    
  }
  else if(length(unique(df$method))==3){
    custom_colors = c("Centralized" = "#7fc97f", "Meta-analysis" = "#beaed4",
                      "Distributed-Chile" = "#fdc086")
  }
  else if(length(unique(df$method))==5 & !(length(which(unique(df$method)=="Distributed-Chile-Authentic")) > 0) ){
    # custom_colors = c("Centralized" = "#CFAE70", "Meta-analysis" = "#E0D5C0", 
    #                   #"Distributed" = "#FEEEBA"
    #                   "Distributed-Chile" = "#B49248", "Distributed-Brazil" = "#ECB748",
    #                   "Distributed-Mexico" = "#946E24")
    custom_colors = c("Centralized" = "#7fc97f", "Meta-analysis" = "#beaed4",
                      "Distributed-Chile" = "#fdc086", "Distributed-Brazil" = "#ffff99",
                      "Distributed-Mexico" = "#386cb0")    
  }
  else if(length(unique(df$method))==5 & length(which(unique(df$method)=="Distributed-Chile-Authentic")) > 0 ){
    custom_colors = c("Centralized" = "#7fc97f", "Meta-analysis" = "#beaed4",
                      "Distributed-Chile" = "#fdc086", "Distributed-Chile-Authentic" = "#FEEEBA",
                      "Distributed-Brazil" = "#ffff99")     
  }
  else{
    custom_colors = c("Centralized" = "#7fc97f", "Meta-analysis" = "#beaed4",
                      "Distributed-Chile" = "#fdc086",
                      "Distributed-Chile-Authentic" = "#FEEEBA",
                      "Distributed-Brazil" = "#ffff99",
                      "Distributed-Mexico" = "#386cb0")
  }
  if(length(unique(df$method))==3){
    ggplot(df, aes(x = coefficient, y = hr, color = method, shape = method)) +
      geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Plot coefficients as points
      geom_errorbar(aes(ymin = lower.95, ymax = upper.95), width = 0.2, 
                    position = position_dodge(width = 0.5)) + # Vertical error bars
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") + # Horizontal reference line at 1
      scale_color_manual(values = custom_colors) +
      theme_minimal() +
      labs(x = "coefficient", y = "HR", color = "Method", shape = "Method",
           title = outcome,
           subtitle = "Forest plot of Hazard Ratio (HR) from different methods") +
      theme_minimal() +
      theme(
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(face = "bold", size = 14),
        legend.position = "top"
        # ,panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
      )    
  }
  else if(length(unique(df$method))==4){
    ggplot(df, aes(x = coefficient, y = hr, color = method, shape = method)) +
      geom_point(size = 3, position = position_dodge(width = 0.5*(length(unique(df$method)))/(length(unique(df$method))-1)) ) +  # Plot coefficients as points
      geom_errorbar(aes(ymin = lower.95, ymax = upper.95), width = 0.2*(length(unique(df$method)))/(length(unique(df$method))-1), 
                    position = position_dodge(width = 0.5*(length(unique(df$method)))/(length(unique(df$method))-1))) + # Vertical error bars
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") + # Horizontal reference line at 1
      scale_color_manual(values = custom_colors) +
      theme_minimal() +
      labs(x = "coefficient", y = "HR", color = "Method", shape = "Method",
           title = outcome,
           subtitle = "Forest plot of Hazard Ratio (HR) from different methods") +
      theme_minimal() +
      theme(
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 14*(length(unique(df$method))-1)/(length(unique(df$method))), face = "bold"),
        legend.title = element_text(face = "bold", size = 14*(length(unique(df$method))-1)/(length(unique(df$method)))),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(face = "bold", size = 14),
        legend.position = "top"
        # ,panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
      )
  }
  else if(length(unique(df$method))==5){
    ggplot(df, aes(x = coefficient, y = hr, color = method, shape = method)) +
      geom_point(size = 3, position = position_dodge(width = 0.8)) +  # Plot coefficients as points
      geom_errorbar(aes(ymin = lower.95, ymax = upper.95), width = 0.2, position = position_dodge(width = 0.8)) + # Vertical error bars
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") + # Horizontal reference line at 1
      scale_color_manual(values = custom_colors) +
      theme_minimal() +
      labs(x = "coefficient", y = "HR", color = "Method", shape = "Method",
           title = outcome,
           subtitle = "Forest plot of Hazard Ratio (HR) from different methods") +
      theme_minimal() +
      theme(
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10.5, face = "bold"),
        legend.title = element_text(face = "bold", size = 10.5),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(face = "bold", size = 14),
        legend.position = "top"
        # ,panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
      ) +
      guides(
        color = guide_legend(nrow = 2),
        shape = guide_legend(nrow = 2)
      )
  }
  else{
    ggplot(df, aes(x = coefficient, y = hr, color = method, shape = method)) +
      geom_point(size = 3, position = position_dodge(width = 0.8)) +  # Plot coefficients as points
      geom_errorbar(aes(ymin = lower.95, ymax = upper.95), width = 0.2, position = position_dodge(width = 0.8)) + # Vertical error bars
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") + # Horizontal reference line at 1
      scale_color_manual(values = custom_colors) +
      theme_minimal() +
      labs(x = "coefficient", y = "HR", color = "Method", shape = "Method",
           title = outcome,
           subtitle = "Forest plot of Hazard Ratio (HR) from different methods") +
      theme_minimal() +
      theme(
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10.5, face = "bold"),
        legend.title = element_text(face = "bold", size = 10.5),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(face = "bold", size = 14),
        legend.position = "top"
        # ,panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
      ) +
      guides(
        color = guide_legend(nrow = 2),
        shape = guide_legend(nrow = 2)
      )    
  }
}
