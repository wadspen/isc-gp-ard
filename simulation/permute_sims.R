source("./simulation_functions.R")
library(future.apply)

args <- commandArgs()


pmethod <- args[6]
nsubjects <- as.numeric(args[7])

K <- 2000
M <- 500
set.seed(21)
seeds <- sample(2^31-1, M, replace = FALSE)

plan(multisession, workers = min(M, 120))  # Windows-friendly
res <- future_lapply(1:M, future.seed = TRUE,  
          FUN = function(m) {
            
            sims_test <- sim_fmri_data(nsubs = nsubjects, seed = seeds[m])
            hdr_df <- sims_test[[2]] %>% 
              filter(time <= 150)
            
            fun_dist <- get_est_fun_dist(hdr_df)
            
            cor_rd_df <- make_cor_rd_df(hdr_df)
            zcor_ts <- est_zcor(cor_rd_df)
            
            perm_zcors <- get_perm_zcors(cor_rd_df, K = K, perm = pmethod)
            act_res <- get_act_res(perm_zcors, hdr_df, zcor_ts)
            all_res <- act_res %>% 
              left_join(fun_dist, by = "voxel") %>% 
              mutate(m = m, method = pmethod)
             
            all_res
          
        })


saveRDS(res, paste0("./sim_res/perm_", pmethod, "_", nsubjects, ".rds"))
