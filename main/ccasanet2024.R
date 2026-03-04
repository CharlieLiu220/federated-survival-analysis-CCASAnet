#--------------------------------------------------#
# Functions defined for analysis of ccasanet data  #
#--------------------------------------------------#

# global objects #
# <https://www.aidsmap.com/about-hiv/types-antiretroviral-medications> #
NNRTI = c("nnrti", "nnrti1", "nnrti2")
NRTI = c("nrti")
PI = c("pi")
PI_new = c("pi_new")
BPI = c("rtv_drug")
INSTI = c("ii1", "ii2")
OTHER = c("t20", "ccr5")



# import data all together #
importccasanet = function(changedir=FALSE){
  if(changedir==TRUE){
    basic = read.csv("../ccasanet_database_20241007/basic.csv")
    art = read.csv("../ccasanet_database_20241007/art.csv")
    visit = read.csv("../ccasanet_database_20241007/visit.csv")
    follow = read.csv("../ccasanet_database_20241007/follow.csv")
    ce = read.csv("../ccasanet_database_20241007/ce.csv")
    cd4 = read.csv("../ccasanet_database_20241007/lab_cd4.csv")
    rna = read.csv("../ccasanet_database_20241007/lab_rna.csv")
    center = read.csv("../Data-Protected/center.csv")    
  }
  else{
    basic = read.csv("./ccasanet_database_20241007/basic.csv")
    art = read.csv("./ccasanet_database_20241007/art.csv")
    visit = read.csv("./ccasanet_database_20241007/visit.csv")
    follow = read.csv("./ccasanet_database_20241007/follow.csv")
    ce = read.csv("./ccasanet_database_20241007/ce.csv")
    cd4 = read.csv("./ccasanet_database_20241007/lab_cd4.csv")
    rna = read.csv("./ccasanet_database_20241007/lab_rna.csv")
    center = read.csv("./Data-Protected/center.csv")    
  }
  return(list(basic=basic, art=art, visit=visit, follow=follow, ce=ce,
              cd4 = cd4, rna = rna, center = center))
}

# convert date variable from character type to Date type #
convertdate = function(datalist){
  basic = datalist$basic
  art = datalist$art
  visit = datalist$visit
  cd4 = datalist$cd4
  rna = datalist$rna
  center = datalist$center
  follow = datalist$follow
  ce = datalist$ce
  ## basic
  date_v_basic = colnames(basic)[str_detect(colnames(basic), "_d$")]
  for(v in date_v_basic){
    basic[,v] = as.Date(basic[,v], "%Y-%m-%d")
  }
  ## art
  art$art_sd = as.Date(art$art_sd, "%Y-%m-%d")
  art$art_ed = as.Date(art$art_ed,"%Y-%m-%d")
  ## visit
  visit$visit_d = as.Date(visit$visit_d, "%Y-%m-%d")
  ## cd4
  cd4$cd4_d = as.Date(cd4$cd4_d, "%Y-%m-%d")
  ## rna
  rna$rna_d = as.Date(rna$rna_d, "%Y-%m-%d")
  ## follow
  follow$lastobs_d = as.Date(follow$lastobs_d, "%Y-%m-%d")
  follow$death_d = as.Date(follow$death_d, "%Y-%m-%d")
  follow$l_alive_d = as.Date(follow$l_alive_d, "%Y-%m-%d")
  ## ce
  ce$ce_d = as.Date(ce$ce_d, "%Y-%m-%d")
  ## center
  center$close_d = as.Date(center$close_d, "%Y-%m-%d")
  center$open_d = as.Date(center$open_d, "%Y-%m-%d")
  center$lastenrol_d = as.Date(center$lastenrol_d, "%Y-%m-%d")
  return(list(basic=basic, art=art, visit=visit, follow=follow, ce=ce,
              cd4 = cd4, rna = rna, center = center))
}

# probable route of infection #
modecategory = function(mode, gender, mode_oth){
  n = length(mode)
  for(i in 1:n){
    if((mode[i]=="Homosexual contact" & gender[i]==1) | (mode[i]=="Bisexual" & gender[i]==1)){
      mode[i] = "MSM"
    }
    else if(mode[i]=="Generic Sexual" | mode[i]=="Heterosexual contact" | 
            (mode[i]=="Homosexual contact" & gender[i]!=1) | (mode[i]=="Bisexual" & gender[i]!=1) |
            (mode[i]=="Other (specify in mode_oth)" & mode_oth[i]=="Sexual")){
      # homosexual + unknown gender = sexual
      mode[i] = "Sexual"
    }
    else if(mode[i]=="Unknown" | mode[i]==""){
      mode[i] = "Unknown"
    }
    else{
      mode[i] = "Other"
    }
  }
  mode = factor(mode, levels = c("Sexual", "MSM", "Other", "Unknown"))
  return(mode)
}

# determine clinical aids status at ART initiation #
clinicalaids = function(aids_first_y, aids_first_d, recart_d){
  n = length(aids_first_y)
  aids_first = vector(length = n)
  for(i in 1:n){
    if(aids_first_y[i]==0){
      aids_first[i] = "Not clinical AIDS"
    }
    else{
      if(aids_first_d[i]-recart_d[i] <= 30){
        aids_first[i] = "Clinical AIDS"
      }
      else{
        aids_first[i] = "Not clinical AIDS"
      }
    }
  }
  aids_first = factor(aids_first, levels = c("Not clinical AIDS", "Clinical AIDS"))
  return(aids_first)
}

# find CD4 count at/near ART initiation #
cd4atinitialART = function(cd4, basic){
  n = nrow(basic)
  cd4count = vector(length = n)
  patients = unique(cd4$patient_id)
  for(i in 1:n){
    if(!(basic$patient_id[i] %in% patients)){
      cd4count[i] = NA
    }
    else{
      tmp = cd4 %>% subset(patient_id == basic$patient_id[i]) %>% arrange(cd4_d)
      if(nrow(tmp)==0){
        cd4count[i] = NA
      }
      else{
        daysdiff = as.numeric(tmp$cd4_d - basic$recart_d[i])
        closedays_index = order(abs(daysdiff))
        l = length(closedays_index)
        index = 1
        while(index < l | index == l){
          if(daysdiff[closedays_index[index]] >= -180 & daysdiff[closedays_index[index]] <= 7){
            cd4count[i] = tmp$cd4_v[closedays_index[index]]
            break
          }
          else{
            index <- index + 1
          }
        }
      }
      if(index == l+1){
        cd4count[i] = NA
      }      
    }
  }
  return(cd4count)
}

# find viral load at/near ART initiation #
vlatinitialART = function(rna, basic){
  n = nrow(basic)
  viralload = vector(length = n)
  patients = unique(rna$patient_id)
  for(i in 1:n){
    if(!(basic$patient_id[i] %in% patients)){
      viralload[i] = NA
    }
    else{
      tmp = rna %>% subset(patient_id == basic$patient_id[i]) %>%
        arrange(rna_d)
      if(nrow(tmp)==0){
        viralload[i] = NA
      }
      else{
        daysdiff = as.numeric(tmp$rna_d - basic$recart_d[i])
        closedays_index = order(abs(daysdiff))
        l = length(closedays_index)
        index = 1
        while(index < l | index == l){
          if(daysdiff[closedays_index[index]] >= -180 & daysdiff[closedays_index[index]] <=0){
            viralload[i] = tmp$rna_v[closedays_index[index]]
            break
          }
          else{
            index = index + 1
          }
        }
        if(index == l+1){
          viralload[i] = NA
        }
      }
    }
  }
  return(viralload)
} 

# COBI is not a PI but the data regard it as a PI #
definenewpi = function(art){
  pi_new = sapply(1:nrow(art), function(x) {
    if(str_detect(art$art_id[x], "COBI")){
      return(art$pi[x]-1)
    }
    else{
      return(art$pi[x])
    }
  })
  return(pi_new)
}

# identify initial ART type (more than HAART or non-HAART) #
## make sure patients selected from basic.csv all have available data in art.csv
identifyinitialARTtype = function(art, basic, factor=T){
  n = nrow(basic)
  initialART = vector(length = n)
  for(i in 1:n){
    tmp = art %>% subset(patient_id == basic$patient_id[i] &
                           art_sd == basic$recart_d[i])
    if(sum(tmp[, NNRTI])>0 & sum(tmp[, c(PI_new, INSTI, OTHER)])==0 & tmp[,BPI]=="No"){
      initialART[i] = "NNRTI"
    }
    else if(sum(tmp[, INSTI])>0 & sum(tmp[, c(NNRTI, PI_new, OTHER)])==0 & tmp[,BPI]=="No"){
      initialART[i] = "INSTI"
    }
    else if((sum(tmp[, PI_new])>0 | tmp[,BPI]=="Yes") & sum(tmp[, c(NNRTI, INSTI, OTHER)])==0){
      initialART[i] = "PI"
    }
    else{
      if(sum(tmp[, OTHER])>0){
        initialART[i] = "Other-t20/ccr5"
      }
      else if(sum(tmp[, NRTI])==3){
        initialART[i] = "Other-3NRTIs"
      }
      else{
        initialART[i] = "Other"
      }
    }
  }
  if(factor){
    initialART = factor(initialART, levels = c("NNRTI","PI", "INSTI","Other-t20/ccr5", "Other-3NRTIs","Other"),
                        labels = c("NNRTI", "PI", "INSTI", "Other", "Other","Other"))    
  }
  return(initialART)
}

# get date of the last observation for each patient in the basic #
getlastobs_d = function(follow, basic){
  index = sapply(basic$patient_id, function(x) which(follow$patient_id==x))
  #return(follow[index, "lastobs_d"])
  return(follow[index, "l_alive_d"])
}

# decide whether a measurement of viral load drops below 200 #
vldropbelow200 = function(rna_v, rna_l){
  # both rna_v and rna_l are vectors
  tmp = sapply(1:length(rna_v), function(x) {
    if( rna_v[x] >= 0){
      return(rna_v[x] < 200);
    }
    else{
      if(!is.na(rna_l[x])){
        ifelse(rna_l[x] < 200, return(T), return(F));
      }
      else{
        # if rna_l is unavailable, then interpret its negative value as smaller than 200
        return(T);
      }
    }    
  })
  return(tmp)
}

# auxiliary function, gettwostraightmorethan200 #
gettwostraightmorethan200 = function(rna_v, vldropbelow200){
  # rna_v here is a vector
  index = which(vldropbelow200==T)[1]
  n = length(rna_v)
  indicator = FALSE
  j = index
  while(!indicator & j< (n-1)){
    indicator = (rna_v[j+1] > 200) & (rna_v[j+2] > 200)
    j = j+1
  }
  if(indicator==T){
    return(list(y=TRUE, dindex=j+1))
  }
  else{
    return(list(y=FALSE, dindex=NA))
  }
}

# identify whether a patient has experienced virological failure within his/her follow up #
identifyvlfailure = function(rna, basic){
  n = nrow(basic)
  vlfailure_y = vector(length = n)
  vlfailure_y_d = vector(length = n)
  for(i in 1:n){
    recart_d_tmp = basic$recart_d[i]
    tmp = rna %>% subset(patient_id == basic$patient_id[i] &
                           rna_d >= recart_d_tmp) %>% arrange(rna_d)
    if(nrow(tmp)==0){
      vlfailure_y[i] = NA
      vlfailure_y_d[i] = NA
    }
    else{
      tmp$rna_v_below_200_y = vldropbelow200(tmp$rna_v, tmp$rna_l)
      tmp = tmp %>% subset(!is.na(rna_v_below_200_y))
      if(nrow(tmp)==0){
        # all missing indicators of below 200; should not be a problem now as all negative values without lower detection limit are regarded as <200
        vlfailure_y[i] = NA
        vlfailure_y_d[i] = NA
      }
      else{
        if(sum(tmp$rna_v_below_200_y)>0){
          below_200_y_d = tmp$rna_d[which(tmp$rna_v_below_200_y==T)[1]]
          tmp = tmp %>% subset(rna_d >= below_200_y_d)
          indicator_1000 = indicator_two200 = FALSE
          # criterion 3
          if(sum(tmp$rna_v > 1000)>0){
            indicator_1000 = TRUE
            indicator_1000_d = tmp$rna_d[which(tmp$rna_v>1000)[1]]
          }
          # criterion 2
          if(gettwostraightmorethan200(tmp$rna_v, tmp$rna_v_below_200_y)$y){
            indicator_two200 = TRUE
            indicator_two200_d = tmp$rna_d[gettwostraightmorethan200(tmp$rna_v, tmp$rna_v_below_200_y)$dindex]        
          }
          if(indicator_1000 & !indicator_two200){
            vlfailure_y[i] = indicator_1000
            vlfailure_y_d[i] = indicator_1000_d
          }
          else if(!indicator_1000 & indicator_two200){
            vlfailure_y[i] = indicator_two200
            vlfailure_y_d[i] = indicator_two200_d
          }
          else if(indicator_1000 & indicator_two200){
            vlfailure_y[i] = indicator_1000 & indicator_two200
            vlfailure_y_d[i] = min(indicator_1000_d, indicator_two200_d)
          }
          else{
            vlfailure_y[i] = FALSE
            vlfailure_y_d[i] = tail(tmp$rna_d, 1)
          }
        }
        else{
          tmp_new = tmp %>% subset(as.numeric(tmp$rna_d-recart_d_tmp) >= 180) # 30 days each month
          if(nrow(tmp_new)!=0){
            # criterion 1
            vlfailure_y[i] = TRUE
            vlfailure_y_d[i] = tail(tmp_new$rna_d, 1)
          }
          else{
            vlfailure_y[i] = FALSE
            vlfailure_y_d[i] = tail(tmp$rna_d, 1)
          }
        }
      }
    }      
  }
  return(data.frame(vlfailure_y = vlfailure_y,
                    vlfailure_y_d = as.Date(vlfailure_y_d)))
}

# identify whether a patient has experienced major treatment change #
## NNRTI <--> PI; PI <--> INSTI; NNRTI <--> INSTI
### auxiliary function: categorize the treatment regimen

ARTtype = function(art){
  n = nrow(art)
  ARTtype = vector(length = n)
  for(i in 1:n){
    if(sum(art[i, NNRTI])>0 & sum(art[i, c(PI_new, INSTI, OTHER)])==0 & art[i,BPI]=="No"){
      ARTtype[i] = "NNRTI"
    }
    else if(sum(art[i, INSTI])>0 & sum(art[i, c(NNRTI, PI_new, OTHER)])==0 & art[i,BPI]=="No"){
      ARTtype[i] = "INSTI"
    }
    else if((sum(art[i, PI_new])>0 | art[i,BPI]=="Yes") & sum(art[i, c(NNRTI, INSTI, OTHER)])==0){
      ARTtype[i] = "PI"
    }
    else{
      if(sum(art[i, OTHER])>0){
        ARTtype[i] = "Other-t20/ccr5"
      }
      else if(sum(art[i, NRTI])==3){
        ARTtype[i] = "Other-3NRTIs"
      }
      else{
        ARTtype[i] = "Other"
      }      
    }    
  }
  return(ARTtype)
}

anychangeofregimen = function(ARTtypehistory){
  n = length(ARTtypehistory)
  if(n==1){
    return(list(indicator=FALSE, index=n, direction = NA))
  }
  else{
    tocheck = sapply(1:(n-1), function(x) ARTtypehistory[x]!=ARTtypehistory[x+1], simplify = T)
    if(sum(tocheck)==0){
      return(list(indicator=FALSE, index = n, direction = NA))
    }
    else{
      # change involving Other is also regarded as a major treatment change
      firstindex = which(tocheck==TRUE)[1] 
      return(list(indicator=TRUE,
                  index = firstindex,
                  direction = paste(ARTtypehistory[firstindex],
                                    ARTtypehistory[firstindex+1], sep = "->")))
      #indicator=FALSE
      #j = 1
      #trueindex = which(tocheck==T)
      #while(!indicator & j <= length(trueindex)){
      #  if(ARTtypehistory[trueindex[j]]!="Other" & ARTtypehistory[trueindex[j]+1]!="Other"){
      #    indicator=TRUE
      #  }
      #  else{
      #    j = j+1
      #  }
      #}
      #if(indicator){
      #  return(list(indicator=indicator, index=trueindex[j]+1, 
      #              direction=paste(ARTtypehistory[trueindex[j]],
      #                              ARTtypehistory[trueindex[j]+1], sep = "->")))        
      #}
      #else{
      #  return(list(indicator=indicator, index=n, direction = NA))
      #}
    }
  }
}

identifyregimenchange = function(art, basic){
  n = nrow(basic)
  regimenchange_y = vector(length = n)
  regimenchange_y_d = vector(length = n)
  regimenchange_y_oth = vector(length = n)
  for(i in 1:n){
    tmp = art %>% subset(patient_id == basic$patient_id[i])
    tmp$ARTtype = ARTtype(tmp)
    changeregimen = anychangeofregimen(tmp$ARTtype)
    if(changeregimen$indicator){
      regimenchange_y[i] = T
      regimenchange_y_d[i] = tmp$art_ed[changeregimen$index]
      regimenchange_y_oth[i] = changeregimen$direction
    }
    else{
      regimenchange_y[i] = F
      regimenchange_y_d[i] = tmp$art_sd[changeregimen$index]
      regimenchange_y_oth[i] = changeregimen$direction
    }
  }
  return(data.frame(ARTchange_y = regimenchange_y,
                    ARTchange_y_d = as.Date(regimenchange_y_d),
                    ARTchange_y_oth = regimenchange_y_oth))
}

# identify composite outcome #

identifycomposite = function(basic){
  n = nrow(basic)
  composite_y = vector(length = n)
  composite_y_d = vector(length = n)
  for(i in 1:n){
    if(is.na(basic$ARTchange_y[i]) | is.na(basic$vlfailure_y[i])){
      composite_y[i] = NA
      composite_y_d[i] = NA
    }
    else if(basic$ARTchange_y[i] & basic$vlfailure_y[i]){
      composite_y[i] = TRUE
      composite_y_d[i] = min(basic$ARTchange_y_d[i], basic$vlfailure_y_d[i])
    }
    else if(basic$ARTchange_y[i] & !basic$vlfailure_y[i]){
      composite_y[i] = TRUE
      composite_y_d[i] = basic$ARTchange_y_d[i]      
    }
    else if(!basic$ARTchange_y[i] & basic$vlfailure_y[i]){
      composite_y[i] = TRUE
      composite_y_d[i] = basic$vlfailure_y_d[i] 
    }
    else{
      composite_y[i] = FALSE
      composite_y_d[i] = min(basic$ARTchange_y_d[i], basic$vlfailure_y_d[i])
    }
  }
  return(data.frame(composite_y = composite_y,
                    composite_y_d = as.Date(composite_y_d)))
}

# identify regimen change status #

identifyregimenchange.status = function(basic, follow, closedate){
  n = nrow(basic)
  status = vector(length = n)
  for(i in 1:n){
    site_tmp = basic$site[i]
    lastobs_tmp = follow$lastobs_d[follow$patient_id==basic$patient_id[i]]
    death_y_tmp = follow$death_y[follow$patient_id==basic$patient_id[i]]
    death_d_tmp = follow$death_d[follow$patient_id==basic$patient_id[i]]
    if(!is.na(basic$ARTchange_y[i])){
      if(basic$ARTchange_y[i]){ # event
        if(basic$ARTchange_y_d[i]<=closedate[site_tmp]){
          # not censored by the defined closing date
          status[i] = 2
        }
        else{
          # censored by the closing date
          status[i] = 0
        }
      }
      else{ # not event
        if(death_y_tmp){ # death
          if(death_d_tmp<=closedate[site_tmp]){ # not censored
            status[i] = 1
          }
          else{
            status[i] = 0
          }
        }
        else{
          status[i] = 0
        }
      }
    }
    else{
      status[i] = NA
    }
  }
  return(status)
}

# identify virological failure status #

identifyvlfailure.status = function(basic, follow, closedate){
  n = nrow(basic)
  status = vector(length = n)
  for(i in 1:n){
    site_tmp = basic$site[i]
    lastobs_tmp = follow$lastobs_d[follow$patient_id==basic$patient_id[i]]
    death_y_tmp = follow$death_y[follow$patient_id==basic$patient_id[i]]
    death_d_tmp = follow$death_d[follow$patient_id==basic$patient_id[i]]
    if(!is.na(basic$vlfailure_y[i])){
      if(basic$vlfailure_y[i]){ # event
        if(basic$vlfailure_y_d[i]<=closedate[site_tmp]){
          # not censored by the defined closing date
          status[i] = 2
        }
        else{
          # censored by the closing date
          status[i] = 0
        }
      }
      else{ # not event
        if(death_y_tmp){ # death
          if(death_d_tmp<=closedate[site_tmp]){ # not censored
            status[i] = 1
          }
          else{
            status[i] = 0
          }
        }
        else{
          status[i] = 0
        }
      }
    }
    else{
      status[i] = NA
    }
  }
  return(status)
}

# identify composite outcome status #

identifycomposite.status = function(basic, follow, closedate){
  n = nrow(basic)
  status = vector(length = n)
  for(i in 1:n){
    site_tmp = basic$site[i]
    lastobs_tmp = follow$lastobs_d[follow$patient_id==basic$patient_id[i]]
    death_y_tmp = follow$death_y[follow$patient_id==basic$patient_id[i]]
    death_d_tmp = follow$death_d[follow$patient_id==basic$patient_id[i]]
    if(!is.na(basic$composite_y[i])){
      if(basic$composite_y[i]){ # event
        if(basic$composite_y_d[i]<=closedate[site_tmp]){
          # not censored by the defined closing date
          status[i] = 2
        }
        else{
          # censored by the closing date
          status[i] = 0
        }
      }
      else{ # not event
        if(death_y_tmp){ # death
          if(death_d_tmp<=closedate[site_tmp]){ # not censored
            status[i] = 1
          }
          else{
            status[i] = 0
          }
        }
        else{
          status[i] = 0
        }
      }
    }
    else{
      status[i] = NA
    }
  }
  return(status)
}



