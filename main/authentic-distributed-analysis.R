#---------------------------------------------------------------------------------#
# Functions defined for AUTHENTIC distributed survival analysis on CCASAnet data  #
#---------------------------------------------------------------------------------#

require(pda)
require(dplyr)
require(rjson)
require(metafor)
require(mice)


# local imputation #
# dat.list is the object returned from mice
mice.local.merge = function(dat.list, m){
  return(sapply(1:m, function(x) complete(dat.list, action = x) %>% 
                  mutate(sex_mode = ifelse(sex=="Female", "Female",
                                           ifelse(mode_new=="MSM", "Male-MSM", "Male-Other")) %>%
                           factor(levels = c("Female", "Male-MSM", "Male-Other"))), simplify = F))
}



# set up #
pda.setup = function(leadsite, control, dir, m, complete = FALSE, outcome){
  if(complete){
    message(sprintf("📢 We will conduct distributed stratified Cox regression for the outcome '%s' based on complete cases.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))    
  }
  else{
    message(sprintf("📢 We will conduct distributed stratified Cox regression for the outcome '%s' with multiple imputation.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))
  }  
  ctrl_tmp = control
  file_path <- ifelse(complete, 
                      sprintf("federated_learning_experiments_with_sites/%s/complete-case", outcome),
                      sprintf("federated_learning_experiments_with_sites/%s", outcome))  
  newdir = paste(dir, "/", outcome, sep = "")
  if(!dir.exists(newdir)){
    dir.create(newdir)
  }  
  if(!complete){
    for(i in 1:m){
      newdir = paste(dir, "/", outcome, "/", i, sep = "")
      ctrl_tmp[["project_name"]] = paste(control[["project_name"]], "_", i, sep = "")
      if(!dir.exists(newdir)){
        dir.create(newdir)
      }    
      pda(site_id = leadsite, control = ctrl_tmp, dir = newdir, upload_without_confirm = T) 
      message(sprintf("✅  [%s]  You can now upload control.json  regarding the %d-th imputed  data set to the OneDrive folder, %s/%d." , 
                      leadsite, i, file_path, i))
    }    
  }
  else{
    newdir = paste(dir, "/", outcome, "/", "complete-case", sep = "")
    ctrl_tmp[["project_name"]] = paste(control[["project_name"]], "_", "complete-case", sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }    
    pda(site_id = leadsite, control = ctrl_tmp, dir = newdir, upload_without_confirm = T)      
    message(sprintf("✅  [%s]  You can now upload control.json to the OneDrive folder, %s." , leadsite, file_path) )
  }
  message(sprintf("✅ [%s] Please inform participating sites to download  control.json from OneDrive." , leadsite ))
}

# initialization #
# dat.list is a list of imputed datasets
pda.local = function(dat.list, localsite, dir, m, complete = FALSE, outcome){
  json_name <- sprintf("%s_initialize.json", localsite)
  file_path <- ifelse(complete, 
                      sprintf("federated_learning_experiments_with_sites/%s/complete-case", outcome),
                      sprintf("federated_learning_experiments_with_sites/%s", outcome))
  if(complete){
    message(sprintf("📢 We are conducting distributed stratified Cox regression for the outcome '%s' based on complete cases.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))    
  }
  else{
    message(sprintf("📢 We are conducting distributed stratified Cox regression for the outcome '%s' with multiple imputation.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))
  }
  if(!complete){
    if(length(dat.list)!=m){
      stop("Check the number of imputed datasets.")
    }
    for(i in 1:m){
      newdir = paste(dir, "/", outcome, "/", i, sep = "")
      if(!dir.exists(newdir)){
        dir.create(newdir)
      }    
      pda(site_id = localsite, ipdata = dat.list[[i]], dir = newdir, upload_without_confirm = T)
      if(localsite != "chile"){
        message(sprintf(
          "✅ [%s] You can now upload %s regarding the %i-th imputed data set to the OneDrive folder , %s/%i.",
          localsite, json_name, i, file_path, i))
      }
      else{
        message(sprintf(
          "✅ [%s] You can now upload the  updated control.json  regarding the %i-th imputed data set to the OneDrive folder, %s.",
          localsite, i, file_path))
      }
    }
  }
  else{
    newdir = paste(dir, "/", outcome, "/", "complete-case", sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }    
    pda(site_id = localsite, ipdata = dat.list, dir = newdir, upload_without_confirm = T)
    if(localsite!="chile"){
      message(sprintf(
        "✅ [%s] You can now upload %s to the OneDrive folder, %s.",
        localsite, json_name, file_path))
    }
    else{
      message(sprintf(
        "✅ [%s] You can now upload the  updated control.json to the OneDrive folder, %s.",
        localsite, file_path))
    }
  }
  if(localsite!="chile"){
    message(sprintf("✅ [%s] Please inform the lead site to download  %s from OneDrive." ,  localsite, json_name))    
  }
  else{
    message(sprintf("✅ [%s] Please inform the participating sites to download updated control.json from OneDrive." ,  localsite ))   
  }
}

# derivation #
# dat.list is a list of imputed datasets
pda.derive = function(dat.list, localsite, dir, m, complete = FALSE, outcome){
  if(complete){
    message(sprintf("📢 We are conducting distributed stratified Cox regression for the outcome '%s' based on complete cases.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))    
  }
  else{
    message(sprintf("📢 We are conducting distributed stratified Cox regression for the outcome '%s' with multiple imputation.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))
  }  
  json_name <- sprintf("%s_derive.json", localsite)
  file_path <- ifelse(complete, 
                      sprintf("federated_learning_experiments_with_sites/%s/complete-case", outcome),
                      sprintf("federated_learning_experiments_with_sites/%s", outcome))  
  if(!complete){
    if(length(dat.list)!=m){
      stop("Check the number of imputed datasets.")
    }
    for(i in 1:m){
      newdir = paste(dir, "/", outcome, "/", i, sep = "")
      if(!dir.exists(newdir)){
        dir.create(newdir)
      }    
      pda(site_id = localsite, ipdata = dat.list[[i]], dir = newdir, upload_without_confirm = T)
      if(localsite != "chile"){
        message(sprintf(
          "✅ [%s] You can now upload %s regarding the %i-th imputed data set to the OneDrive folder , %s/%i",
          localsite, json_name, i, file_path, i))
      }
      else{
        message(sprintf(
          "✅ [%s] You can now move to the 'Estimation' stage.",
          localsite))
      }
    }
  }
  else{
    newdir = paste(dir, "/", outcome, "/", "complete-case", sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }    
    pda(site_id = localsite, ipdata = dat.list, dir = newdir, upload_without_confirm = T)
    if(localsite!="chile"){
      message(sprintf(
        "✅ [%s] You can now upload %s to the OneDrive folder, %s",
        localsite, json_name, file_path))
    }
    else{
      message(sprintf(
        "✅ [%s] You can now move to the 'Estimation' stage.",
        localsite))
    }
  }
  if(localsite!="chile"){
    message(sprintf("✅ [%s] Please inform the lead site to download  %s from OneDrive." ,  localsite, json_name))    
  }
}

# estimation #
pda.estimate = function(dat.list, leadsite, dir, m, complete = FALSE, outcome){
  if(complete){
    message(sprintf("📢 We are conducting distributed stratified Cox regression for the outcome '%s' based on complete cases.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))    
  }
  else{
    message(sprintf("📢 We are conducting distributed stratified Cox regression for the outcome '%s' with multiple imputation.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))
  }  
  if(!complete){
    if(length(dat.list)!=m){
      stop("Check the number of imputed datasets.")
    }
    for(i in 1:m){
      newdir = paste(dir, "/", outcome, "/", i, sep = "")
      if(!dir.exists(newdir)){
        dir.create(newdir)
      }    
      pda(site_id = leadsite, ipdata = dat.list[[i]], dir = newdir, upload_without_confirm = T)
      message(sprintf("✅ Done for the %d-th imputed data set!", i))
    }    
  }
  else{
    newdir = paste(dir, "/", outcome, "/", "complete-case", sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }    
    pda(site_id = leadsite, ipdata = dat.list, dir = newdir, upload_without_confirm = T)      
  }
  message(sprintf(
    "✅ [%s] You can now move to the 'Results' stage.",
    leadsite))  
}

# summary results of distributed analysis #
pda.summary = function(leadsite, localsites, dir, m, complete = FALSE, outcome){
  if(complete){
    message(sprintf("📢 We are conducting distributed stratified Cox regression for the outcome '%s' based on complete cases.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))    
  }
  else{
    message(sprintf("📢 We are conducting distributed stratified Cox regression for the outcome '%s' with multiple imputation.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))
  }  
  file_path <- sprintf("federated_learning_experiments_with_sites/%s", outcome)      
  results = vector("list", m)
  size = vector("numeric", length = 1+length(localsites))
  if(!complete){
    for(i in 1:m){
      newdir = paste(dir, "/", outcome, "/", i, sep = "")
      if(!dir.exists(newdir)){
        dir.create(newdir)
      }    
      config <- getCloudConfig(site_id = leadsite, dir=newdir)
      fit.pda <- pdaGet(name = paste0(c(leadsite, "estimate"), collapse = "_"), config = config)
      for(j in 1:length(size)){
        if(j==1){
          size[[j]] = fit.pda$site_size
        }
        else{
          tmp = fromJSON(file = paste(newdir, "/",localsites[j-1], "_initialize.json", sep = ""))
          size[[j]] = tmp$site_size
        }
      }
      ##
      results.tmp = cbind(fit.pda$btilde, exp(fit.pda$btilde), sqrt(diag(solve(fit.pda$Htilde)/sum(size))))
      colnames(results.tmp) = c("estimate","hr", "std.error")
      results.tmp = as.data.frame(results.tmp)
      results.tmp$lower.95 = exp(results.tmp$estimate - qnorm(.975) * results.tmp$std.error)
      results.tmp$upper.95 = exp(results.tmp$estimate + qnorm(.975) * results.tmp$std.error)
      results.tmp$p.value = 2*pnorm(abs(fit.pda$btilde)/results.tmp$std.error,
                                    lower.tail = F)
      results[[i]] = results.tmp    
    }    
  }
  else{
    newdir = paste(dir, "/", outcome, "/", "complete-case", sep = "")
    if(!dir.exists(newdir)){
      dir.create(newdir)
    }    
    config <- getCloudConfig(site_id = leadsite, dir=newdir)
    fit.pda <- pdaGet(name = paste0(c(leadsite, "estimate"), collapse = "_"), config = config)
    for(j in 1:length(size)){
      if(j==1){
        size[[j]] = fit.pda$site_size
      }
      else{
        tmp = fromJSON(file = paste(newdir, "/",localsites[j-1], "_initialize.json", sep = ""))
        size[[j]] = tmp$site_size
      }
    }
    ##
    results.tmp = cbind(fit.pda$btilde, exp(fit.pda$btilde), sqrt(diag(solve(fit.pda$Htilde)/sum(size))))
    colnames(results.tmp) = c("estimate","hr", "std.error")
    results.tmp = as.data.frame(results.tmp)
    results.tmp$lower.95 = exp(results.tmp$estimate - qnorm(.975) * results.tmp$std.error)
    results.tmp$upper.95 = exp(results.tmp$estimate + qnorm(.975) * results.tmp$std.error)
    results.tmp$p.value = 2*pnorm(abs(fit.pda$btilde)/results.tmp$std.error,
                                  lower.tail = F)
    results[[1]] = results.tmp        
  }
  message(sprintf(
    "✅ [%s] Done! You can upload the results for the distributed  analysis in this scenario to the OneDrive folder, %s" ,
    leadsite, file_path))  
  return(list(results = results,
              site_size = size))
}

# meta with multiple imputation #
dist.meta = function(sites, dir, m, complete = FALSE, outcome){
  if(complete){
    message(sprintf("📢 We are conducting meta-analysis for the outcome '%s' based on complete cases.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))    
  }
  else{
    message(sprintf("📢 We are conducting meta-analysis for the outcome '%s' with multiple imputation.", ifelse(outcome=="vl", "Virological Failure", "Major Regimen Change")))
  }  
  file_path <- sprintf("federated_learning_experiments_with_sites/%s", outcome)    
  results = vector("list", length(sites))
  if(!complete){
    for(j in 1:length(sites)){
      tmp = vector("list", m)
      for(i in 1:m){
        newdir = paste(dir, "/", outcome, "/", i, sep = "")
        if(!dir.exists(newdir)){
          dir.create(newdir)
        }
        result_tmp = fromJSON(file = paste(newdir, "/", sites[j], "_initialize.json", sep = "")) 
        tmp[[i]] = data.frame(estimate = result_tmp$bhat_i,
                              std.error = sqrt(result_tmp$Vhat_i))
      }
      results[[j]] = mi.pool.rubin(tmp, result_tmp$site_size)
    }    
  }
  else{
    for(j in 1:length(sites)){
        newdir = paste(dir, "/", outcome, "/", "complete-case", sep = "")
        if(!dir.exists(newdir)){
          dir.create(newdir)
        }
        result_tmp = fromJSON(file = paste(newdir, "/", sites[j], "_initialize.json", sep = "")) 
        tmp = data.frame(estimate = result_tmp$bhat_i,
                              std.error = sqrt(result_tmp$Vhat_i))
        results[[j]] = tmp
      }
  }
  message(sprintf(
    "✅ [%s] Done! You can upload the results for the meta-analysis in this scenario to the OneDrive folder, %s" ,
    sites[1], file_path))  
  return(results)
}
