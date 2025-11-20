.libPaths("~/rlibs")
source("./analysis_functions.R")



mod <- cmdstan_model("../stan_models/spat_nngp_lm_hs2.stan", 
                     include_paths = gptools_include_path())

proc_type <- "filter_gsr"
task <- "checker"
nut <- 1000
nul <- 1000

data_loc <- paste("../../data_by_roi", task, proc_type, sep = "/")
file_names <- list.files(data_loc)

options(future.debug = FALSE)

plan(multisession, workers = min(length(file_names) + 4, 108))  # Windows-friendly
res <- future_lapply(1:length(file_names),
                     
  function(ind) {
    
    file <- paste(data_loc, file_names[ind], sep = "/")
    hdr_df_pb <- readRDS(file)
    stan_data <- get_stan_data(hdr_df_pb, ind, nu_t = nut, nu_l = nul)
    fit <- mod$sample(data = stan_data[[1]], 
                      chains = 1, 
                      # parallel_chains = 4,
                      iter_warmup = 1000, 
                      iter_sampling = 2000)
    
    
    res <- get_results(fit, stan_data)
    res
  }
)

save_file <- paste0(task, "_", proc_type, ".rds")

saveRDS(res, save_file)
