source("./simulation_functions.R")


nsubjects <- 10
pmethod <- "shift"
K <- 40
M <- 500
set.seed(21)
seeds <- sample(2^31-1, M, replace = FALSE)

mod <- cmdstan_model("../stan_models/spat_nngp_lm_hs_sb_ar1.stan", 
                     include_paths = gptools_include_path())


plan(multisession, workers = min(M, 110))  # Windows-friendly
res <- future_lapply(1:M,  
          function(m) {
            
            testvox <- expand.grid(x = 9:12, y = 7:9) %>% 
              mutate(voxel = paste(x, y, sep = "_"))
            sims_test <- sim_fmri_data(nsubs = nsubjects, seed = seeds[m])
            hdr_df <- sims_test[[2]] %>% 
              filter(time <= 150, voxel %in% testvox$voxel)
            
             
            stan_data <- get_stan_data(hdr_df)
            
            
            fit <- mod$sample(data = stan_data, 
                              chains = 1, 
                              # parallel_chains = 4,
                              iter_warmup = 1000, 
                              iter_sampling = 1000)
            
            un_vox <- unique(hdr_df$voxel)
            draws <- fit$draws(format = "df")
            
            predf_dist <- get_sgp_fun_dist(draws, hdr_df, un_vox)
            param_sums <- get_gp_param_sums(draws, un_vox, 
                                            nsubjects, ntime = stan_data$N)
            
            
            all_res <- get_gp_act_res(hdr_df, param_sums) %>% 
              mutate(m = m, method = "sgp")
            
          
        })



