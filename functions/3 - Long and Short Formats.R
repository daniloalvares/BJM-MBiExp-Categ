format_change <- function(data_curr, data_prev){
  
  # Data by biomarker
  data_M <- data_prev[which(!is.na(data_prev$m_spike_serum)),]
  data_F <- data_prev[which(!is.na(data_prev$flc_serum)),]
  
  # Unique ID
  u_ID <- unique(data_curr$patientid)
  
  # Start and stop positions of longitudinal measurements per patient
  start_M <- stop_M <- rep(NA,length(u_ID))
  start_F <- stop_F <- rep(NA,length(u_ID))
  aux_M <- aux_F <- 0
  
  # Variables for short format
  sex <- race <- ecog <- iss <- chain_status <- rep(NA,length(u_ID))
  age <- albumin <- beta2_microglobulin <- creatinine <- hemoglobin <- rep(NA,length(u_ID))
  ldh <- lymphocyte <- neutrophil <- platelet <- rep(NA,length(u_ID))
  immunoglobulin_a <- immunoglobulin_g <- immunoglobulin_m <- rep(NA,length(u_ID))
  
  for(i in 1:length(u_ID)){
    pos <- which(data_curr$patientid==u_ID[i])
    pos_M <- which(data_M$patientid==u_ID[i])
    pos_F <- which(data_F$patientid==u_ID[i])
    
    # M-spike
    if(length(pos_M)>0){
      start_M[i] <- aux_M + 1
      stop_M[i] <- aux_M + length(pos_M)
      aux_M <- stop_M[i]
    }
    
    # FLC
    if(length(pos_F)>0){
      start_F[i] <- aux_F + 1
      stop_F[i] <- aux_F + length(pos_F)
      aux_F <- stop_F[i]
    }
    
    # Categorical baseline variables
    sex[i] <- data_curr$gender[pos][1]
    race[i] <- data_curr$race_ethnicity[pos][1]
    ecog[i] <- data_curr$ecoggrp[pos][1]
    iss[i] <- data_curr$iss_stage[pos][1]
    chain_status[i] <- data_curr$chain_status[pos][1]
    
    # Continuous baseline variables
    age[i] <- data_curr$age[pos][1]
    albumin[i] <- data_curr$albumin_bl[pos][1]
    beta2_microglobulin[i] <- data_curr$beta2_microglobulin_bl[pos][1]
    creatinine[i] <- data_curr$creatinine_bl[pos][1]
    hemoglobin[i] <- data_curr$hemoglobin_bl[pos][1]
    ldh[i] <- data_curr$ldh_bl[pos][1]
    lymphocyte[i] <- data_curr$lymphocyte_bl[pos][1]
    neutrophil[i] <- data_curr$neutrophil_count_bl[pos][1]
    platelet[i] <- data_curr$platelet_count_bl[pos][1]
    immunoglobulin_a[i] <- data_curr$immunoglobulin_a_bl[pos][1]
    immunoglobulin_g[i] <- data_curr$immunoglobulin_g_bl[pos][1]
    immunoglobulin_m[i] <- data_curr$immunoglobulin_m_bl[pos][1]
  }
  
  # Reference: Male
  sex <- ifelse(sex=="Female",1,0)
  
  # Reference: Non-Hispanic White
  race_0 <- ifelse(race=="Non-Hispanic White",1,0)
  race_1 <- ifelse(race=="Non-Hispanic Black",1,0)
  race_2 <- ifelse(race=="Hispanic or Latino (any)" | race=="Non-Hispanic Asian" | race=="Other Non-Hispanic",1,0)
  race_3 <- ifelse(race=="Not reported",1,0)
  
  # Reference: 0
  ecog_0 <- ifelse(ecog=="0",1,0)
  ecog_1 <- ifelse(ecog=="1",1,0)
  ecog_2 <- ifelse(ecog=="2" | ecog=="3+",1,0)
  ecog_3 <- ifelse(ecog=="Not reported",1,0)
  
  # Reference: Stage I
  iss_0 <- ifelse(iss=="Stage I",1,0)
  iss_1 <- ifelse(iss=="Stage II",1,0)
  iss_2 <- ifelse(iss=="Stage III",1,0)
  iss_3 <- ifelse(iss=="Unknown/not documented",1,0)
  
  Short <- data.frame(patientid=data_curr[!duplicated(data_curr$patientid),"patientid"],
                      linename=data_curr[!duplicated(data_curr$patientid),"linename"],
                      ID1=data_prev[!duplicated(data_prev$patientid),"ID1"],
                      sex=sex, race_0=race_0, race_1=race_1, race_2=race_2, race_3=race_3, 
                      ecog_0=ecog_0, ecog_1=ecog_1, ecog_2=ecog_2, ecog_3=ecog_3, 
                      iss_0=iss_0, iss_1=iss_1, iss_2=iss_2, iss_3=iss_3, chain_status=chain_status,
                      age=age, albumin=albumin, beta2_microglobulin=beta2_microglobulin,
                      creatinine=creatinine, hemoglobin=hemoglobin, ldh=ldh, lymphocyte=lymphocyte,
                      neutrophil=neutrophil, platelet=platelet, immunoglobulin_a=immunoglobulin_a, 
                      immunoglobulin_g=immunoglobulin_g, immunoglobulin_m=immunoglobulin_m,
                      start_M=start_M, stop_M=stop_M, start_F=start_F, stop_F=stop_F, 
                      duration=data_prev[!duplicated(data_prev$patientid),"surv"]/365.25)
  
  # Carfilzomib and Pomalidomide
  carfil <- ifelse(str_detect(Short$linename,"Carfilzomib"),1,0)
  pomali <- ifelse(str_detect(Short$linename,"Pomalidomide"),2,0)
  Short$drugcat <- carfil + pomali
  
  Short$drugcat[Short$drugcat==0] <- "Other"
  Short$drugcat[Short$drugcat==1] <- "Carfilzomib"
  Short$drugcat[Short$drugcat==2] <- "Pomalidomide"
  #Short$drugcat[Short$drugcat==3] <- "Carfilzomib+Pomalidomide"
  Short$drugcat[Short$drugcat==3] <- "Other"
  
  # Reference: Other
  Short$drug_0 <- ifelse(Short$drugcat=="Other",1,0)
  Short$drug_1 <- ifelse(Short$drugcat=="Carfilzomib",1,0)
  Short$drug_2 <- ifelse(Short$drugcat=="Pomalidomide",1,0)
  # Short$drug_3 <- ifelse(Short$drugcat=="Carfilzomib+Pomalidomide",1,0)
  
  
  # Longitudinal information
  M_Spike <- data.frame(patientid=data_M$patientid, ID1=data_M$ID1, y=data_M$m_spike_serum, time=as.numeric(data_M$timefromlot)/365.25)
  FLC <- data.frame(patientid=data_F$patientid, ID1=data_F$ID1, y=data_F$flc_serum, time=as.numeric(data_F$timefromlot)/365.25)
  

  return( list(Short=Short, Long=list(M_Spike=M_Spike, FLC=FLC)) )
  
}


datalot2 <- format_change(data_curr=datalot2_curr, data_prev=datalot1_prev)
datalot3 <- format_change(data_curr=datalot3_curr, data_prev=datalot2_prev)
datalot4 <- format_change(data_curr=datalot4_curr, data_prev=datalot3_prev)


rm(list=setdiff(ls(), c("datalot2")))
