source("./simulation_functions.R")
library(parallel)
library(doParallel)
library(doMC)

args <- commandArgs()
nsubjects <- as.numeric(args[6])
pmethod <- "sgp_ar1"
K <- 40
M <- 300
set.seed(21)
seeds <- sample(2^31-1, M, replace = FALSE)

mod <- cmdstan_model("../stan_models/spat_nngp_lm_hs_sb_ar1.stan", 
                     include_paths = gptools_include_path())

m <- 4
n.cores <- 110
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)



test <- foreach(m = 1:M, .errorhandling = "pass",
        .packages = c("tidyr", "dplyr", "stringr", "neuRosim", "MASS",
                      "purrr", "cmdstanr", "abind", "schoolmath",
                      "PTHfftSURROGATES", "gptoolsStan")

) %dopar% {
  
  # testvox <- expand.grid(x = 26:27, y = 15:16) %>% 
  #   mutate(voxel = paste(x, y, sep = "_"))
  sims_test <- sim_fmri_data(nsubs = nsubjects, seed = seeds[m])
  hdr_df <- sims_test[[2]] %>% 
    filter(time <= 60)#, voxel %in% testvox$voxel)
  
  rm(sims_test)
  gc()
  stan_data <- get_stan_data(hdr_df)
  
  start <- Sys.time()
  fit <- mod$sample(data = stan_data, 
                    chains = 1, 
                    # parallel_chains = 4,
                    iter_warmup = 1000, 
                    iter_sampling = 1000,
		    refresh = 2)
  end <- Sys.time()

  
  fit_time <- as.numeric(difftime(end, start, units = "min"))
  un_vox <- unique(hdr_df$voxel)
  draws <- fit$draws(format = "df")
  rm(fit);
  gc()
  predf_dist <- get_sgp_fun_dist(draws, hdr_df, un_vox)
  param_sums <- get_gp_param_sums(draws, un_vox, 
                                  nsubjects, ntime = stan_data$N)
  
  
  
  all_res <- get_gp_act_res(hdr_df, param_sums) %>% 
    mutate(m = m, method = pmethod, fittime = fit_time) %>%
    left_join(predf_dist, by = "voxel")
  
  saveRDS(all_res, paste0("./sim_res/sgp/", pmethod, "_", nsubjects, "_iter", m, ".rds")) 
  #rm(all_res)
  #gc()
  NULL
    
  }


print(test)


