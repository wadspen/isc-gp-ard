source("./simulation_functions.R")
library(future.apply)
library(parallel)
library(doParallel)
library(doMC)
#n.cores <- detectCores()
n.cores <- 110
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

args <- commandArgs()
nul <- as.numeric(args[6])
sigma_rho <- as.numeric(args[7])
nsubjects <- as.numeric(args[8])
eff <- as.numeric(args[9])
K <- 40
M <- 150
set.seed(21)
seeds <- sample(2^31-1, M, replace = FALSE)

mod <- cmdstan_model("../stan_models/spat_nngp_lm_hs2.stan", 
                     include_paths = gptools_include_path())

warm <- 2000
samp <- 3000
pmethod <- "sgp"
m <- 4
# plan(multisession, workers = min(M, 100))  # Windows-friendly
# res <- future_lapply(M:1,  
#           function(m) {


foreach(m = 1:M,
                    .packages = c("cmdstanr", "dplyr", "neuRosim", "tidyr",
                                  "MASS", "purrr", "abind", "schoolmath",
                                  "PTHfftSURROGATES", "gptoolsStan", "stringr")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %dopar% {
                      
            sims_test <- sim_fmri_data(nsubs = nsubjects, dims = c(11,11), 
                                       eff_size = eff, stim = "natural",
                                       seed = seeds[m])
            hdr_df <- sims_test[[2]] %>% 
              filter(time <= time_lim)#, voxel %in% testvox$voxel)
            
            rm(sims_test)
	          gc()
            stan_data <- get_stan_data(hdr_df, nu_l = nul, 
                                       sigma_rho_n = sigma_rho)
            
            start <- Sys.time()
            fit <- mod$sample(data = stan_data, 
                              chains = 1, 
                              # parallel_chains = 4,
                              iter_warmup = warm, 
                              iter_sampling = samp)
            end <- Sys.time()
	          fit_time <- as.numeric(difftime(end, start, units = "min"))
            un_vox <- unique(hdr_df$voxel)
            draws <- fit$draws(format = "df")
            rm(fit); rm(st_data)
	          gc()
            predf_dist <- get_sgp_fun_dist(draws, hdr_df, un_vox)
            param_sums <- get_gp_param_sums(draws, un_vox, 
                                            nsubjects, ntime = stan_data$N)
            
            
            all_res <- get_gp_act_res(hdr_df, param_sums) %>% 
              mutate(m = m, method = pmethod, fittime = fit_time) %>%
	            left_join(predf_dist, by = "voxel")
            
            all_res$nul <- nul
            all_res$sig_rho_p <- sigma_rho
            all_res$eff <- eff
            all_res$warmup <- warm
            all_res$sample <- samp
            all_res$time_lim <- time_lim
           
            saveRDS(all_res, paste0("./results/sgp/", 
                                    "nul", nul,
                                    "_srho", sigma_rho,
                                    "_nsub", nsubjects, 
                                    "_eff", eff,
                                    "_iter", m, ".rds")) 
            rm(all_res)
	          gc()
	          NULL
            
    }       
       # })



