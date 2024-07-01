vars <- c("age","albumin_bl","beta2_microglobulin_bl","creatinine_bl","hemoglobin_bl","ldh_bl","lymphocyte_bl",
          "neutrophil_count_bl","platelet_count_bl","immunoglobulin_a_bl","immunoglobulin_g_bl","immunoglobulin_m_bl")


# LOG SCALE
datalot1_prev[,vars] <- log(datalot1_prev[,vars] + 0.1)
datalot2_curr[,vars] <- log(datalot2_curr[,vars] + 0.1)
datalot2_prev[,vars] <- log(datalot2_prev[,vars] + 0.1)
datalot3_curr[,vars] <- log(datalot3_curr[,vars] + 0.1)
datalot3_prev[,vars] <- log(datalot3_prev[,vars] + 0.1)
datalot4_curr[,vars] <- log(datalot4_curr[,vars] + 0.1)


for(vv in vars){
  # Mean and standard deviation
  mean_bl <- mean(datalot1_prev[!duplicated(datalot1_prev$patientid),vv], na.rm = TRUE)
  sd_bl <- sd(datalot1_prev[!duplicated(datalot1_prev$patientid),vv], na.rm = TRUE)
  
  # STANDARDISATION
  datalot1_prev[,vv] <- (datalot1_prev[,vv] - mean_bl) / sd_bl
  datalot2_prev[,vv] <- (datalot2_prev[,vv] - mean_bl) / sd_bl
  datalot2_curr[,vv] <- (datalot2_curr[,vv] - mean_bl) / sd_bl
  datalot3_prev[,vv] <- (datalot3_prev[,vv] - mean_bl) / sd_bl
  datalot3_curr[,vv] <- (datalot3_curr[,vv] - mean_bl) / sd_bl
  datalot4_curr[,vv] <- (datalot4_curr[,vv] - mean_bl) / sd_bl
  
  # IMPUTATION
  datalot1_prev[is.na(datalot1_prev[,vv]),vv] <- 0
  datalot2_prev[is.na(datalot2_prev[,vv]),vv] <- 0
  datalot2_curr[is.na(datalot2_curr[,vv]),vv] <- 0
  datalot3_prev[is.na(datalot3_prev[,vv]),vv] <- 0
  datalot3_curr[is.na(datalot3_curr[,vv]),vv] <- 0
  datalot4_curr[is.na(datalot4_curr[,vv]),vv] <- 0
}


rm(list=setdiff(ls(), c("datalot1_prev","datalot2_prev","datalot3_prev","datalot2_curr","datalot3_curr","datalot4_curr")))
