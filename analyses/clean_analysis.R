.libPaths("~/rlibs")
source("./analysis_functions.R")



mod <- cmdstan_model("../stan_models/spat_nngp_lm_hs2.stan", 
                     include_paths = gptools_include_path())


atlas_df <- get_atlas()
all_bold <- read_bold_data(subs = 6, checker = TRUE)
coords_int_full <- get_shared_coords(all_bold)
hdr_df_pb <- get_hdr_df_pb(all_bold, atlas_df)
#hdr_df_pb$roi <- 1
rois <- unique(hdr_df_pb$roi) 


options(future.debug = FALSE)

plan(multisession, workers = min(length(rois) + 4, 63))  # Windows-friendly
res <- future_lapply(rois,
                     
  function(ind) {
    stan_data <- get_stan_data(hdr_df_pb, ind)
    stan_data[[1]]$nul <- 1000
    fit <- mod$sample(data = stan_data[[1]], 
                      chains = 1, 
                      # parallel_chains = 4,
                      iter_warmup = 1000, 
                      iter_sampling = 2000)
    
    
    res <- get_results(fit, stan_data)
    res
  }
)



saveRDS(res, "./spat_lm_ucf_sb_checker.rds")
