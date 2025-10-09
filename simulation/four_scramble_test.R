.libPaths("~/rlibs")
source("./sim_neuro_data.R")
library(jsonlite)
library(R.matlab)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RNifti)
library(cmdstanr)
library(gptoolsStan)
library(stringr)
library(future.apply)
library(abind)
library(schoolmath)
library(PTHfftSURROGATES)

fishZ <- function(rho) {
  rhoZ <- .5 * log((1 + rho) / (1 - rho))
  return(rhoZ)
}


bootstrap <- function(data, n) {
  resampled_data <- lapply(1:n, function(i) {
    resample <- sample(data, replace = TRUE)
    # Perform desired operations on the resampled data, e.g., compute a statistic
    # and return the result
  })
  return(resampled_data)
}

coords <- df_long %>% 
  dplyr::select(x,y) %>% 
  unique() %>% 
  as.matrix()


est_cor <- matrix(NA, nrow = 53, ncol = 31)
est_zcor <- est_cor
for (c in 1:nrow(coords)) {
  x <- coords[c, 1]; y <- coords[c, 2]
  sub_cor <- matrix(NA, nrow = length(subjects), ncol = length(subjects))  
  for (i in 1:length(subjects)) {
    for (j in i:length(subjects)) {
      sub_cor[i,j] <- as.numeric(cor(scale(sub_data[[i]][x,y,1:50]), 
                                     scale(sub_data[[j]][x,y,1:50])))
      sub_cor[j,i] <- sub_cor[i,j]
      
    }
  }
  diag(sub_cor) <- 0
  rho <- median(sub_cor, na.rm = TRUE)
  est_cor[x, y] <- rho
  est_zcor[x, y] <- fishZ(rho)
  print(c)
  
}

M <- 1000
boot_samps <- bootstrap(subjects, M)
boot_zcor <- array(NA, dim = c(53, 31, M))

#43, 16; 838
#20, 24; 1239
#3, 13; 639
# cs <- which(coords[,1] %in% 8:18 & coords[,2] == 13)

plan(multisession, workers = min(12, 120))

# for (m in 1:M) {
res <- future_lapply(1:M,
                     
   function(m) {
     # cs <- which(coords[,1] %in% 8:18 & coords[,2] == 13)
     tryCatch({ 
        # boot_subs <- boot_samps[[m]]
       perm_est_zcor <- matrix(NA, nrow = 53, ncol = 31)
       library(schoolmath)
        for (c in 1:nrow(coords)) {
          x <- coords[c, 1]; y <- coords[c, 2]
          sub_cor <- matrix(NA, nrow = length(subjects), 
                            ncol = length(subjects))  
          for (i in subjects) {
            for (j in subjects) {
              ts1 <- as.numeric(scale(sub_data[[i]][x,y,1:50]))
              ts2 <- as.numeric(scale(sub_data[[j]][x,y,1:50]))
              # ts1 <- data.table::shift(ts1, 
              #                   n = sample(length(ts1), 1), 
              #                   type = "cyclic")
              ts1 <- phase_scramble2(ts2)
              sub_cor[i,j] <- as.numeric(cor(ts1,
                                             ts2))
              sub_cor[j,i] <- sub_cor[i,j]
              
            }
          }
          diag(sub_cor) <- 0
          rho <- median(sub_cor, na.rm = TRUE)
          perm_est_zcor[x, y] <- fishZ(rho)
          print(c)
          
        }
        print(m)
        # boot_zcor[,,m] <- 
        perm_est_zcor
     },
     error = function(e) {
       message(sprintf("Error in iteration %d: %s", m, e$message))
       return(NULL)   # or NULL, or some sentinel value
     })
}, future.stdout = TRUE, future.seed = TRUE)
# }

sig_map <- matrix(NA, nrow = 53, ncol = 31)
sim_arr <- simplify2array(res)
for (i in 1:53) {
  for (j in 1:31) {
    sig_map[i,j] <- mean(abs(sim_arr[i,j,]) > est_zcor[i,j], na.rm = TRUE)
    # sig_map[i,j] <- ifelse(sig_map[i,j] < .000000001, 1, 0)
  }
}


saveRDS(sig_map, "./sim_sig_map.rds")


































