source("../simulation/simulation_functions.R")
source("./analysis_functions.R")
library(future.apply)
library(hrf)

args <- commandArgs()


pmethod <- args[6]
task <- args[7]
proc_type <- args[8]

K <- 2000

time_s = ifelse(task == "checker", 97, 288)
atlas_df <- get_atlas()
print("First checkpoint")
print(time_s)
all_bold <- read_bold_data(checker = str_detect(task, "checker"), filter = str_detect(proc_type, "filter"),
                           gsr = !str_detect(proc_type, "nogsr"), time_steps = time_s)
coords_int_full <- get_shared_coords(all_bold)

save_dir <- paste0("./results/", task, "/perm/", pmethod, "/") 
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

data_loc <- paste("../../data_by_roi", task, proc_type, sep = "/")
file_names <- list.files(data_loc)
options(future.debug = FALSE)

plan(multisession, workers = min(length(file_names) + 4, 108))  # Windows-friendly
res <- future_lapply(1:length(file_names), future.seed = TRUE,  
                     FUN = function(index) {
                       
                       file <- paste(data_loc, file_names[index], sep = "/")
                       hdr_df_pb <- readRDS(file)
                       
                       
                      
                       num_sub_vox <- length(unique(hdr_df_pb$voxel))* 
                         length(unique(hdr_df_pb$participant_id))
                       
                       hdr_df_pb$hdr_conv <- rep(hdr_conv, num_sub_vox)
                       
                       
                       hdr_df_pb <- hdr_df_pb %>% 
                         mutate(subject = as.numeric(str_extract(participant_id, "\\d+")))
                       
                       cor_rd_df <- make_cor_rd_df(hdr_df_pb)
                       zcor_ts <- est_zcor(cor_rd_df)
                       
                       perm_zcors <- get_perm_zcors(cor_rd_df, K = K, perm = pmethod)
                       perm_res <- get_perm_res(perm_zcors, hdr_df_pb, zcor_ts)
                       
                       
                       if (task == "checker") {
                         
                         
                         EVs <- list(taskA=cbind(on=c(20,60,100,140), dr=rep(20,4)))
                         TR <- 2.1
                         nTime <- 97
                         des <- make_design(EVs, nTime, TR)
                         hdr_conv <- as.numeric(scale(as.numeric(des$design)))
                         fun_dist <- get_est_fun_dist(hdr_df_pb)
                         
                         
                         all_res <- perm_res %>% 
                           left_join(fun_dist, by = "voxel") %>% 
                           mutate(method = pmethod) %>% 
                           mutate(roi = unique(hdr_df_pb$roi))
                       } else {
                         
                         all_res <- perm_res %>% 
                           mutate(method = pmethod) %>% 
                           mutate(roi = unique(hdr_df_pb$roi))
                           
                       }
                       
                       all_res
                       
                       
                     })



saveRDS(res, paste0(save_dir, proc_type, ".rds"))

