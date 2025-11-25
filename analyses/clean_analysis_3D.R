.libPaths("~/rlibs")
source("./analysis_functions.R")

args <- commandArgs()
task <- args[6]
proc_type <- args[7]
model <- args[8]
nut <- as.numeric(args[9])
nul <- as.numeric(args[10])
print(proc_type)
print(model)
print(nut)
print(nul)
cmstn_file <- ifelse(model == "glm", "../stan_models/spat_nngp_glm.stan",
		     "../stan_models/spat_nngp_lm_hs2.stan")


atlas_df <- get_atlas()
all_bold <- read_bold_data(checker = str_detect(task, "checker"), filter = str_detect(proc_type, "filter"),
			   gsr = !str_detect(proc_type, "nogsr"))
coords_int_full <- get_shared_coords(all_bold)

mod <- cmdstan_model(cmstn_file, 
                     include_paths = gptools_include_path())

#proc_type <- "filter_gsr"
#task <- "checker"
#nut <- 1000
#nul <- 1000

data_loc <- paste("../../data_by_roi", task, proc_type, sep = "/")
file_names <- list.files(data_loc)

options(future.debug = FALSE)
index <- 1
plan(multisession, workers = min(length(file_names) + 4, 50))  # Windows-friendly
res <- future_lapply(1:length(file_names),
                    
  function(index) {
    library(hrf)
    EVs <- list(taskA=cbind(on=c(20,60,100,140), dr=rep(20,4)))
    TR <- 2.1
    nTime <- 97
    des <- make_design(EVs, nTime, TR) 
    file <- paste(data_loc, file_names[index], sep = "/")
    hdr_df_pb <- readRDS(file)
    ind <- unique(hdr_df_pb$roi)
    stan_data <- get_stan_data3D(hdr_df_pb, ind, nu_t = nut, nu_l = nul)
    rm(hdr_df_pb)
    gc()
    stan_data[[1]]$hrf <- as.numeric(scale(as.numeric(des$design)))
    fit <- mod$sample(data = stan_data[[1]], 
                      chains = 1, 
                      # parallel_chains = 4,
                      iter_warmup = 1000, 
                      iter_sampling = 1000)
    
    
    res <- get_results(fit, stan_data)
    saveRDS(res, paste0("./results/", task, "/", model, "/", proc_type,
		       	"/nut", nut, "/nul", nul, "/roi", ind, ".rds"))
    rm(res)
    gc()
 }
)

#save_file <- paste0(task, "_", proc_type, "_", model, "_nut", nut, "_nul", nul, ".rds")

#saveRDS(res, save_file)
