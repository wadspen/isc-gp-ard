.libPaths("~/rlibs")
library(jsonlite)
library(R.matlab)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RNifti)
# library(magick)
# library(hrf)
library(cmdstanr)
library(gptoolsStan)
library(stringr)
library(future.apply)
library(abind)
#library(INLA)

data_loc <- "../../dme_files/"


mod <- cmdstan_model("../stan_models/mc_nngp_hs.stan", 
                     include_paths = gptools_include_path())
subjects <- as.character(1:22)
subjects <- ifelse(nchar(subjects) == 1, paste0(0, subjects), subjects)

sub_names <- paste0("sub-", subjects)
fold2 <- "/ses-01/func/"
file <- "_ses-01_task-dme_run-01_bold/func_preproc/func_pp_filter_gsr_sm0.mni152.3mm.nii.gz"
mask_name <- "_ses-01_task-dme_run-01_bold/func_seg/wm_mask.nii"
all_bold <- list()
for (i in 1:length(subjects)) {
  sub <- paste0("sub-", subjects[i])
  file_path <- paste0(data_loc, sub, fold2, sub, file)
  mask_path <- paste0(data_loc, sub, fold2, sub, mask_name)
  if (subjects[i] == "17") {
    while (str_detect(file_path, "ses-01")) {
      file_path <- str_replace(file_path, "ses-01", "ses-03")
    }
  }
  bold <- readNifti(file_path)
  all_bold[[i]] <- bold[,,,1:288]
  print(subjects[i])
  
}



all_coords <- data.frame()
for (i in 1:length(all_bold)) {
  coords <- which(all_bold[[i]][,,,31] != 0, arr.ind = TRUE) %>%
    as.data.frame()
  coords$sub <- subjects[i]
  all_coords <- rbind(all_coords, coords)
}

full_coords <- all_coords %>% 
  mutate(voxel = paste(dim1, dim2, dim3, sep = "_")) %>% 
  group_by(voxel) %>% 
  summarise(n = n()) %>% 
  filter(n == length(subjects)) %>%
  unique()

coords_int_full <- all_coords %>%
  mutate(voxel = paste(dim1, dim2, dim3, sep = "_")) %>% 
  filter(voxel %in% unique(full_coords$voxel)) %>% 
  distinct(sub, dim1, dim2, dim3) %>%
  group_by(dim1, dim2, dim3) %>%
  summarise(n_groups = n_distinct(sub), .groups = "drop") %>%
  filter(n_groups == n_distinct(all_coords$sub)) %>% 
  filter(dim3 == 27) %>% 
  select(-n_groups)

colnames(coords_int_full) <- c("x", "y", "z")




options(future.debug = FALSE)


ab_df <- data.frame()
for (i in 1:length(all_bold)) {
  arr <- all_bold[[i]][,,27,]
  nx <- dim(arr)[1]
  ny <- dim(arr)[2]
  nt <- dim(arr)[3]
  
  df <- data.frame(participant_id = sub_names[i],
                   expand.grid(x = 1:nx, y = 1:ny),
                   matrix(aperm(arr, c(3,1,2)), nrow = nx*ny, ncol = nt, byrow = TRUE)
  )
  
  # Name time columns nicely
  colnames(df)[-(1:3)] <- paste0("t", 1:nt)
  ab_df <- rbind(ab_df, df)
  
}

hdr_df_ab <- ab_df %>% 
  pivot_longer(4:ncol(ab_df), names_to = "time", values_to = "hdr") %>% 
  mutate(time = as.numeric(str_replace(time, "t", ""))) %>% 
  filter(time < 51)

hdr_df_ab$hdr <- as.numeric(hdr_df_ab$hdr)

hdr_df_ab <- hdr_df_ab %>%
	group_by(participant_id, x, y) %>%
	mutate(hdr = as.numeric(scale(hdr, center = TRUE, scale = TRUE))) %>%
	ungroup() %>%
	filter(!is.na(hdr))
vs <- 25
param_means_all <- data.frame()
for (v in 1:3) {
  start <- (v - 1) * vs + 1
  end   <- v * vs
  
  hdr_df_pb <- hdr_df_ab %>% 
    filter(between(y, start, end)) %>% 
    mutate(voxel = paste(x, y, sep = "_"))


  
  


  coords_int <- coords_int_full %>% 
    mutate(voxel = paste(x, y, sep = "_")) %>% 
    filter(voxel %in% unique(hdr_df_pb$voxel))
    # filter(x %in% unique(hdr_df_pb$x) & y %in% unique(hdr_df_pb$y))
  
  voxels <- unique(coords_int$voxel)
  
  #coords_int <- coords_int[sample(nrow(coords_int), 36, replace = FALSE),]
  
  
  vox_ids <- ceiling(seq_along(voxels)/ 10)
  coords_int$vox_id <- vox_ids
  
  
  print(paste0("Iteration ", v, " of ", 3, " with ",
              dim(coords_int)[1], " voxels."))
  plan(multisession, workers = 120)  # Windows-friendly

# coords_int <- expand.grid(1:7, 1:5, 1)
res <- future_lapply(unique(vox_ids),
                     
                     function(ind) {
                       tryCatch({  
                         
                         # for (ind in 1:36) {
                         # point <- as.numeric(coords_int[ind, 1:3])
                         # xc <- point[1]; yc <- point[2]; zc <- point[3]
                         
         hdr_df <- hdr_df_pb %>% 
             # filter(x == xc, y == yc) %>%
             filter(voxel %in% coords_int$voxel[coords_int$vox_id == ind]) %>% 
             select(participant_id, hdr, time, voxel)
         
         n <- length(unique(hdr_df$time))
         x <- unique(hdr_df$time)
         time <- unique(hdr_df$time)
         un_vox <- unique(hdr_df$voxel)
                         
		   	 
# 		     y_dat <- hdr_df %>%
#   				 pivot_wider(time, names_from = participant_id, 
#   				             values_from = hdr) %>%
#   				 select(-time) %>%
#   				 t() %>%
#   				 as.matrix()
         
         y_dat_list <- list()
         for (i in 1:length(un_vox)) {
           y_dat <- hdr_df %>% 
             filter(voxel == un_vox[i]) %>% 
             # select(-time) %>% 
             # mutate(subject = paste("sub", subject, sep = "-")) %>% 
             pivot_wider(id_cols = time, names_from = participant_id, 
                         values_from = hdr) %>% 
             select(-time) %>% 
             t() %>% 
             as.matrix()
           # print(dim(y_dat))
           y_dat_list[[i]] <- y_dat
         }
         
         y_arr <- abind::abind(y_dat_list, along = 3)
			 
         edge_index <- rbind(1:(ncol(y_dat) - 1), 2:ncol(y_dat))
  			 
         stan_data <- list(N = n, S = nrow(y_dat), 
                           C = length(unique(hdr_df$voxel)),
                           y = y_arr,
                           x = matrix(unique(hdr_df$time), ncol = 1)*2,
                           dx = 1*2,
                           nf = n%/%2 + 1,
                           edge_index = edge_index,
                           mu_rho = log(1000),
                           sigma_rho = .1,
                           mu_sigma = 0,
                           sigma_sigma = .01,
                           mu_tau = 1,
                           sigma_tau = .1,
                           r = 1,
                           m = .5,
                           sigma = 3.5)
			 fit <- mod$sample(data = stan_data, chains = 1, 
			                   iter_warmup = 1000, iter_sampling = 1000)
			 
			 draws <- fit$draws(format = "df")
			 
			 
			 preds_df <- data.frame()
			 for (i in 1:length(un_vox)) {
			   drawszf <- draws %>% 
			     select(contains("f[") & contains(paste0(",",i,"]")))
			   
			   preds <- apply(drawszf, MARGIN = 2, FUN = mean)
			   upp <- apply(drawszf, MARGIN = 2, FUN = quantile, probs = .975)
			   # t <- 1:250
			   
			   # plot(preds~t, type = "l")
			   
			   pdf <- data.frame(pred = preds, time = unique(hdr_df$time), 
			                     voxel = un_vox[i])
			   
			   preds_df <- rbind(pdf, preds_df)

			 }
			 
			 rownames(preds_df) <- 1:nrow(preds_df)
			 
			 tau_sigs <- draws %>% 
			   select(contains("tau_sig")) %>% 
			   as.matrix() %>% 
			   as.numeric()
			 
			 lambdas <- draws %>% 
			   select(contains("lambda")) %>% 
			   pivot_longer(1:length(un_vox),
			                names_to = "param",
			                values_to = "post")
			 
			 
			 kappa_df <- data.frame(lambdas, 
			                        tausig = rep(tau_sigs, 
			                                     length(un_vox))) %>% 
			   rowwise() %>% 
			   mutate(kappa = 1/(1 + tausig^2 * post^2)) %>% 
			   mutate(param = str_replace(param, 
			                              "lambda", "kappa")) %>% 
			   group_by(param) %>% 
			   summarise(m = mean(kappa),
			             sd = sd(kappa))
			 
			 kappa_df$voxel <- un_vox
			 
			 rho_sum <- draws %>% 
			   select(contains("rho")) %>% 
			   pivot_longer(1:length(un_vox), 
			                names_to = "param", values_to = "post") %>% 
			   group_by(param) %>% 
			   summarise(m = mean(post, na.rm = TRUE),
			             sd = sd(post, na.rm = TRUE))
			 
			 
			 rho_sum$voxel <- un_vox
			 
			 
			 tau_sum <- draws %>% 
			   select(contains("tau")) %>% 
			   pivot_longer(1:length(un_vox), 
			                names_to = "param", values_to = "post") %>% 
			   group_by(param) %>% 
			   summarise(m = mean(post, na.rm = TRUE),
			             sd = sd(post, na.rm = TRUE))
			 
			 
			 tau_sum$voxel <- un_vox
			 
			 sigma_sum <- data.frame(param = "sigma", m = mean(draws$sigma, na.rm = TRUE), 
			                         sd = sd(draws$sigma, na.rm = TRUE), 
			                         voxel = un_vox)
			 
			 param_sums <- rbind(rho_sum, tau_sum, sigma_sum, kappa_df)
			 
			 fin_preds <- preds_df %>% 
			   left_join(param_sums, by = "voxel")
			 
			 fin_preds

			 
                         
                       },
                       error = function(e) {
                         message(sprintf("Error in iteration %d: %s", ind, e$message))
                         return(NULL)   # or NULL, or some sentinel value
                       })
                     })
    param_means <- do.call(rbind, res)
    
    param_means_all <- rbind(param_means_all, param_means)

}

saveRDS(param_means_all, "./test_means_hs.rds")
