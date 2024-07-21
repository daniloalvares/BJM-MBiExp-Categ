# Function for MAP
MAP_fc <- function(x){
    
  lim.inf <- min(x)-1
  lim.sup <- max(x)+1
  s <- density(x,from=lim.inf,to=lim.sup,bw=0.2)
  n <- length(s$y)
  v1 <- s$y[1:(n-2)]
  v2 <- s$y[2:(n-1)]
  v3 <- s$y[3:n]
  ix <- 1+which((v1<v2)&(v2>v3))
  out <- s$x[which(s$y==max(s$y))]
  
  return( out )
    
}


# Fitting the joint model
fit_jm <- function(data, iter=5000, warmup=1000, chains=3){
  
  model <- cmdstan_model(JM_Biexp_Cat_FIT)

  vars <- c("sex","race_1","race_2","race_3","ecog_1","ecog_2","ecog_3","iss_1","iss_2","iss_3","age","albumin","beta2_microglobulin",
            "creatinine","hemoglobin","ldh","lymphocyte","neutrophil","platelet","immunoglobulin_a","immunoglobulin_g","immunoglobulin_m","duration")
  X <- cbind(1, data$Short[,vars])
  
  # Model data
  N_M <- nrow(data$Long$M_Spike)
  N_F <- nrow(data$Long$FLC)
  n_M <- length(unique(data$Long$M_Spike$patientid))
  n_F <- length(unique(data$Long$FLC$patientid))
  y_M <- data$Long$M_Spike$y
  y_F <- data$Long$FLC$y
  times_M <- data$Long$M_Spike$time
  times_F <- data$Long$FLC$time
  n <- nrow(data$Short)
  nbeta <- ncol(X)
  ID_M1 <- (data$Long$M_Spike %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_F1 <- (data$Long$FLC %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_M2 <- data$Long$M_Spike$ID1
  ID_F2 <- data$Long$FLC$ID1
  z <- apply(cbind(data$Short$drug_1,data$Short$drug_2,data$Short$drug_0),1,which.max)

  # Setting initial values
  inits <- list(beta_raw=matrix(0,nrow=ncol(X),ncol=2), alpha_raw=matrix(0,nrow=6,ncol=2),
                theta_M=c(0,0,0), theta_F=c(0,0,0), sigma2_M=1, sigma2_F=1,
                Omega_M=diag(rep(1,3)), Omega_F=diag(rep(1,3)),
                bi_M=matrix(0,nrow=n_M,ncol=3), bi_F=matrix(0,nrow=n_F,ncol=3))
  init_fun <- rep(list(inits), chains)
 
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
    
  i.time <- Sys.time()
  fit <- suppressMessages(suppressWarnings( model$sample(data = list(N_M=N_M, N_F=N_F, n_M=n_M, n_F=n_F, n=n, 
                                                                     nbeta=nbeta, y_M=log(y_M+1), y_F=log(y_F+1),
                                                                     ID_M1=ID_M1, ID_M2=ID_M2, ID_F1=ID_F1, ID_F2=ID_F2, 
                                                                     times_M=times_M, times_F=times_F, z=z, X=X),
                                                         init            = init_fun,
                                                         iter_sampling   = iter,
                                                         iter_warmup     = warmup,                 
                                                         chains          = chains,
                                                         thin            = 1,
                                                         seed            = 123,
                                                         parallel_chains = getOption("mc.cores",chains)) ))
  e.time <- Sys.time()
  
  return( list(fit=fit, time=e.time-i.time) )
  
}


# Generated quantities from the fitted joint model
gq_jm <- function(fit, data){
  
  model <- cmdstan_model(JM_Biexp_Cat_VI)

  vars <- c("sex","race_1","race_2","race_3","ecog_1","ecog_2","ecog_3","iss_1","iss_2","iss_3","age","albumin","beta2_microglobulin",
            "creatinine","hemoglobin","ldh","lymphocyte","neutrophil","platelet","immunoglobulin_a","immunoglobulin_g","immunoglobulin_m","duration")
  X <- cbind(1, data$Short[,vars])
  
  # Model data
  N_M <- nrow(data$Long$M_Spike)
  N_F <- nrow(data$Long$FLC)
  n_M <- length(unique(data$Long$M_Spike$patientid))
  n_F <- length(unique(data$Long$FLC$patientid))
  y_M <- data$Long$M_Spike$y
  y_F <- data$Long$FLC$y
  times_M <- data$Long$M_Spike$time
  times_F <- data$Long$FLC$time
  n <- nrow(data$Short)
  nbeta <- ncol(X)
  ID_M1 <- (data$Long$M_Spike %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_F1 <- (data$Long$FLC %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_M2 <- data$Long$M_Spike$ID1
  ID_F2 <- data$Long$FLC$ID1
  z <- apply(cbind(data$Short$drug_1,data$Short$drug_2,data$Short$drug_0),1,which.max)
    
  i.time <- Sys.time()
  fit <- suppressMessages(suppressWarnings( model$generate_quantities(fit,
                                                                      data = list(N_M=N_M, N_F=N_F, n_M=n_M, n_F=n_F, n=n, 
                                                                                  nbeta=nbeta, y_M=log(y_M+1), y_F=log(y_F+1),
                                                                                  ID_M1=ID_M1, ID_M2=ID_M2, ID_F1=ID_F1, ID_F2=ID_F2, 
                                                                                  times_M=times_M, times_F=times_F, z=z, X=X,
                                                                                  IMP_B_M=1:n, IMP_B_F=1:n,
                                                                                  IMP_G_M=1:n, IMP_G_F=1:n,
                                                                                  IMP_D_M=1:n, IMP_D_F=1:n),
                                                                      seed = 123) ))
  e.time <- Sys.time()
  
  # WAIC and LOO
  ww <- waic(fit$draws("log_lik"))$estimates[3,1]
  ll <- loo(fit$draws("log_lik"))$estimates[3,1]

  return( list(fit=fit, time=e.time-i.time, waic=ww, loo=ll) )
  
}


# Performance metrics for the joint model
# - Longitudinal submodels: Individual weighted residuals (IWRES)
# - Categorical submodel: Accuracy, Precision, Recall, and F1-score
metrics_jm <- function(fit1, fit2, data){

  # LONGITUDINAL SUBMODELS
  # Estimated trajectory for M-spike and FLC from the bi-exponential submodels
  est_M <- apply(fit2$draws("nonlinpred_M"), 3, mean)
  est_F <- apply(fit2$draws("nonlinpred_F"), 3, mean)

  # Estimated residual standard deviation for M-spike and FLC from the bi-exponential submodels
  sig_M <- apply(sqrt(fit1$draws("sigma2_M")), 3, mean)
  sig_F <- apply(sqrt(fit1$draws("sigma2_F")), 3, mean)

  # Individual weighted residuals (IWRES)
  iwres_M <- (log(data$Long$M_Spike$y+1) - est_M) / sig_M
  iwres_F <- (log(data$Long$FLC$y+1) - est_F) / sig_F

  dta <- data.frame(time = c(data$Long$M_Spike$time,data$Long$FLC$time),
                  col = c(data$Long$M_Spike$patientid,data$Long$FLC$patientid),
                  biom = c(rep("M-spike",nrow(data$Long$M_Spike)),rep("FLC",nrow(data$Long$FLC))),
                  res = c(iwres_M,iwres_F))
  dta$biom <- factor(dta$biom, levels = c("M-spike","FLC"))

  print(ggplot(data = dta, aes(x=time, y=res)) + geom_point(aes(x=time, y=res), size=0.5, color="black") +
    geom_hline(yintercept=0, linetype="dashed", color="red") + xlab("Time (in years)") + ylab("IWRES") + theme_bw() +
    theme(legend.position="none", panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
    facet_wrap( ~ biom, ncol = 2) + scale_x_continuous(breaks = 0:6))

  # CATEGORICAL SUBMODEL
  z <- apply(cbind(data$Short$drug_1,data$Short$drug_2,data$Short$drug_0),1,which.max)
  n <- length(z)
  z_pred <- rep(NA,n)
  for(i in 1:n){
    cat1 <- MAP_fc(as.vector(fit2$draws("probs")[,,i]))
    cat2 <- MAP_fc(as.vector(fit2$draws("probs")[,,n+i]))
    cat3 <- MAP_fc(as.vector(fit2$draws("probs")[,,2*n+i]))
    z_pred[i] <- which.max(c(cat1, cat2, cat3))
  }

  # Accuracy, Precision, Recall, and F1-score
  met <- crf_evaluation(pred=z_pred, obs=z)
  print(met$overall[c(1,2,3,5)])

}


# Permutation-based variable importance
vi_jm <- function(fit1, fit2, data, runs=50, criterion="WAIC"){
  
  # WAIC and LOO from the fitted joint model
  ww <- fit2$waic
  ll <- fit2$loo

  # Variable positions for permutations
  i.change <- 1 + c(1,2,5,8,11:23)
  e.change <- 1 + c(1,4,7,10,11:23)

  model <- cmdstan_model(JM_Biexp_Cat_VI)

  vars <- c("sex","race_1","race_2","race_3","ecog_1","ecog_2","ecog_3","iss_1","iss_2","iss_3","age","albumin","beta2_microglobulin",
            "creatinine","hemoglobin","ldh","lymphocyte","neutrophil","platelet","immunoglobulin_a","immunoglobulin_g","immunoglobulin_m","duration")
  X <- cbind(1, data$Short[,vars])
  
  # Model data
  N_M <- nrow(data$Long$M_Spike)
  N_F <- nrow(data$Long$FLC)
  n_M <- length(unique(data$Long$M_Spike$patientid))
  n_F <- length(unique(data$Long$FLC$patientid))
  y_M <- data$Long$M_Spike$y
  y_F <- data$Long$FLC$y
  times_M <- data$Long$M_Spike$time
  times_F <- data$Long$FLC$time
  n <- nrow(data$Short)
  nbeta <- ncol(X)
  ID_M1 <- (data$Long$M_Spike %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_F1 <- (data$Long$FLC %>% group_by(patientid) %>% dplyr::mutate(ID = cur_group_id()))$ID
  ID_M2 <- data$Long$M_Spike$ID1
  ID_F2 <- data$Long$FLC$ID1
  z <- apply(cbind(data$Short$drug_1,data$Short$drug_2,data$Short$drug_0),1,which.max)
    
  waics <- matrix(NA, nrow=length(i.change)+6, ncol=runs)
  loos <- matrix(NA, nrow=length(i.change)+6, ncol=runs)

  for(j in 1:(length(i.change)+6)){
  
    print(j)
    for(k in 1:runs){
      Xpred <- X
      IMP_B_M <- IMP_B_F <- IMP_G_M <- IMP_G_F <- IMP_D_M <- IMP_D_F <- 1:n
      if(j <= length(i.change)){
        # Baseline covariates
        Xpred[,i.change[j]:e.change[j]] <- Xpred[sample(1:n, size = n, replace = F),i.change[j]:e.change[j]]
      }else{
        # M-spike baseline
        if(j == (length(i.change)+1)){ IMP_B_M <- sample(1:n, size = n, replace = F) }
        # FLC baseline
        if(j == (length(i.change)+2)){ IMP_B_F <- sample(1:n, size = n, replace = F) }
      
        # M-spike growth rate
        if(j == (length(i.change)+3)){ IMP_G_M <- sample(1:n, size = n, replace = F) }
        # FLC growth rate
        if(j == (length(i.change)+4)){ IMP_G_F <- sample(1:n, size = n, replace = F) }
      
        # M-spike decay rate
        if(j == (length(i.change)+5)){ IMP_D_M <- sample(1:n, size = n, replace = F) }
        # FLC decay rate
        if(j == (length(i.change)+6)){ IMP_D_F <- sample(1:n, size = n, replace = F) }
      }
    
      fit <- suppressMessages(suppressWarnings( model$generate_quantities(fit1,
                                                                          data = list(N_M=N_M, N_F=N_F, n_M=n_M, n_F=n_F, n=n, 
                                                                                      nbeta=nbeta, y_M=log(y_M+1), y_F=log(y_F+1),
                                                                                      ID_M1=ID_M1, ID_M2=ID_M2, ID_F1=ID_F1, ID_F2=ID_F2, 
                                                                                      times_M=times_M, times_F=times_F, z=z, X=Xpred,
                                                                                      IMP_B_M=IMP_B_M, IMP_B_F=IMP_B_F,
                                                                                      IMP_G_M=IMP_G_M, IMP_G_F=IMP_G_F,
                                                                                      IMP_D_M=IMP_D_M, IMP_D_F=IMP_D_F),
                                                                          seed = 123) ))
      if(criterion == "WAIC"){ 
         waics[j,k] <- waic(fit$draws("log_lik"))$estimates[3,1] 
      }else{
         loos[j,k] <- loo(fit$draws("log_lik"))$estimates[3,1]
      }
    }
  }

  varnames <- c("Sex","Ethnicity","ECOG","ISS","Age","Albumin","B2M","Creatinine","Hemoglobin","LDH","Lymphocyte","Neutrophil","Platelet",
                "IgA","IgG","IgM","RVd duration","Baseline M-spike","Growth M-spike","Decay M-spike","Baseline FLC","Growth FLC","Decay FLC")

  if(criterion == "WAIC"){
     dta <- data.frame(name=varnames, Importance=apply(waics-ww,1,mean), min=pmax(0,apply(waics-ww,1,min)), max=apply(waics-ww,1,max))
  }else{
     dta <- data.frame(name=varnames, Importance=apply(loos-ll,1,mean), min=pmax(0,apply(loos-ll,1,min)), max=apply(loos-ll,1,max))
  }
  dta <- transform(dta, variable=reorder(name, Importance))


  print(dta %>% mutate(name = fct_reorder(name, Importance)) %>%
          ggplot(aes(x=name, y=Importance)) + theme_bw() + geom_bar(stat="identity", fill="gray") +
          geom_errorbar(aes(x=name, ymin=min, ymax=max), width=0.3, colour="black", linewidth=0.4) +
          theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.title.y=element_blank()) + coord_flip())

  return( dta )
  
}
