source("./simulation_functions.R")
library(future.apply)

args <- commandArgs()



nsubjects <- as.numeric(args[6])
pmethod <- "blme"


mod <- cmdstan_model("../stan_models/blme.stan")
M <- 300
set.seed(21)
seeds <- sample(2^31-1, M, replace = FALSE)
m <- 4
plan(multisession, workers = min(M, 100))  # Windows-friendly
res <- future_lapply(1:M, future.seed = TRUE,  
          FUN = function(m) {
            
            # svox <- expand.grid(x = 7:8, y = 11:12) %>%
            #   mutate(voxel = paste(x, y, sep = "_"))
            sims_test <- sim_fmri_data(nsubs = nsubjects, seed = seeds[m])
            hdr_df <- sims_test[[2]] %>% 
              filter(time <= 150)#, voxel %in% svox$voxel)
            rm(sims_test)
            gc()
            fun_dist <- get_est_fun_dist(hdr_df)
            cor_rd_df <- make_cor_rd_df(hdr_df)
            
            stan_data <- get_blme_stan_data(cor_rd_df, nsubjects)
            
            start <- Sys.time()
            fit <- mod$sample(data = stan_data,
                              chains = 1,
                              iter_warmup = 1000,
                              iter_sampling = 1000)
            end <- Sys.time()
            fit_time <- as.numeric(difftime(end, start, units = "min"))
            draws <- fit$draws(format = "df")
            rm(fit); rm(stan_data)
            gc()
            un_vox = unique(cor_rd_df$voxel)
            
            par_res <- get_blme_par_res(draws, un_vox)
            
            all_res <- par_res %>% 
              left_join(fun_dist, by = "voxel") %>% 
              mutate(m = m, method = pmethod, fittime = fit_time)
            
            saveRDS(all_res, paste0("./sim_res/blme/", 
                                    pmethod, "_", nsubjects, "_iter", m, ".rds"))
            
            rm(all_res)
            gc()
            NULL
            
          
        })



