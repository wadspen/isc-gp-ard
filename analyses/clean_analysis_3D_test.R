.libPaths("~/rlibs")
source("./analysis_functions.R")
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
task <- args[6]
proc_type <- args[7]
model <- args[8]
nut <- as.numeric(args[9])
nul <- as.numeric(args[10])
sigma_rho <- as.numeric(args[11])
print(proc_type)
print(model)
print(nut)
print(nul)
print(sigma_rho)
#proc_type <- "filter_gsr"
#task <- "checker"
#nut <- 1000
#nul <- 1000
#model <- "glm"
cmstn_file <- ifelse(model == "glm", "../stan_models/spat_nngp_glm.stan",
                     "../stan_models/spat_nngp_lm_hs.stan")

time_s = ifelse(task == "checker", 97, 288)
atlas_df <- get_atlas()
print("First checkpoint")
print(time_s)
all_bold <- read_bold_data(checker = str_detect(task, "checker"), filter = str_detect(proc_type, "filter"),
                           gsr = !str_detect(proc_type, "nogsr"), time_steps = time_s)
coords_int_full <- get_shared_coords(all_bold)

mod <- cmdstan_model(cmstn_file, 
                     include_paths = gptools_include_path())

save_dir <- paste0("./results/", task, "/", model, "/")
file_part <- paste0("nut", nut, 
                    "_nul", nul,
                    "_sig", sigma_rho)

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}
#proc_type <- "filter_gsr"
#task <- "checker"
#nut <- 1000
#nul <- 1000

warm <- 2000
samp <- 3000

data_loc <- paste("../../data_by_roi", task, proc_type, sep = "/")
file_names <- list.files(data_loc)
print("Second checkpoint")
options(future.debug = FALSE)
index <- 1




# plan(multisession, workers = min(length(file_names) + 4, 108))  # Windows-friendly
# res <- future_lapply(1:length(file_names),
#                     
#   function(index) {




foreach(m = M:1,
        .packages = c("cmdstanr", "neuRosim", "dplyr", "stringr", "tidyr", 
                      "MASS", "purrr", "abind", "schoolmath", 
                      "PTHfftSURROGATES", "gptoolsStan", "hrf", "abind")
        ,.errorhandling = "remove"
        ,.combine = rbind) %dopar% {
          
          file <- paste(data_loc, file_names[index], sep = "/")
          hdr_df_pb <- readRDS(file)
          ind <- unique(hdr_df_pb$roi)
          stan_data <- get_stan_data3D(hdr_df_pb, ind, nu_t = nut, nu_l = nul,
                                       sigma_rho_n = sigma_rho); stan_data[[1]]$C
          rm(hdr_df_pb)
          gc()
          
          if (model == "glm") {
            # library(hrf)
            EVs <- list(taskA=cbind(on=c(20,60,100,140), dr=rep(20,4)))
            TR <- 2.1
            nTime <- 97
            des <- make_design(EVs, nTime, TR)
            stan_data[[1]]$hrf <- as.numeric(scale(as.numeric(des$design)))
          }
          
          
          start <- Sys.time()
          fit <- mod$sample(data = stan_data[[1]], 
                            chains = 1, 
                            # parallel_chains = 4,
                            iter_warmup = warm, 
                            iter_sampling = samp)
          end <- Sys.time()
          fit_time <- as.numeric(difftime(end, start, units = "min"))
          
          
          res <- get_results(fit, stan_data)
          
          res$nul <- nul
          res$sig_rho_p <- sigma_rho
          res$index <- index
          res$fit_time <- fit_time
          res$proc_type <- proc_type
          
          saveRDS(res, paste0(save_dir, file_part,
                              "_index", index, 
                              "_times", time_s, 
                              "_pt_", proc_type, ".rds"))
          rm(res)
          gc()
          NULL
        }
# )

#save_file <- paste0(task, "_", proc_type, "_", model, "_nut", nut, "_nul", nul, ".rds")

#saveRDS(res, save_file)
