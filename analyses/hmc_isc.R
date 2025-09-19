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
#library(INLA)

data_loc <- "../../dme_files/"


mod <- cmdstan_model("../stan_models/gtools_nn.stan", include_paths = gptools_include_path())
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
  mutate(time = as.numeric(str_replace(time, "t", "")))

hdr_df_ab$hdr <- as.numeric(hdr_df_ab$hdr)

hdr_df_ab <- hdr_df_ab %>%
	group_by(participant_id, x, y) %>%
	mutate(hdr = scale(hdr, center = TRUE, scale = TRUE)) %>%
	ungroup() %>%
	filter(!is.na(hdr))
n <- 25
param_means_all <- data.frame()
for (i in 1:3) {
  start <- (i - 1) * n + 1
  end   <- i * n
  
  hdr_df_pb <- hdr_df_ab %>% 
    filter(between(y, start, end))

 
  


coords_int <- coords_int_full %>% 
  filter(x %in% unique(hdr_df_pb$x) & y %in% unique(hdr_df_pb$y))

#coords_int <- coords_int[sample(nrow(coords_int), 36, replace = FALSE),]


print("what about Bob?")
plan(multisession, workers = 120)  # Windows-friendly

# coords_int <- expand.grid(1:7, 1:5, 1)
res <- future_lapply(seq_len(nrow(coords_int)),
                     
                     function(ind) {
                       tryCatch({  
                         
                         # for (ind in 1:36) {
                         point <- as.numeric(coords_int[ind, 1:3])
                         xc <- point[1]; yc <- point[2]; zc <- point[3]
                         
                         hdr_df <- hdr_df_pb %>% 
                           filter(x == xc, y == yc) %>% 
                           select(participant_id, hdr, time)
                         
		   	 hdr_df$hdr <- as.numeric(hdr_df$hdr)
		         y_dat <- hdr_df %>%
				 pivot_wider(time, names_from = participant_id, values_from = hdr) %>%
				 select(-time) %>%
				 t() %>%
				 as.matrix()
			 
			 time <- unique(hdr_df$time)
			 edge_index <- rbind(1:(ncol(y_dat) - 1), 2:ncol(y_dat))
		         stan_data <- list(N = length(unique(hdr_df$time)), S = nrow(y_dat), y = y_dat,
					   x = matrix(unique(hdr_df$time), ncol = 1),
					   edge_index = edge_index,
					   mu_rho = 0,
					   sigma_rho = .01,
					   mu_sigma = 0,
					   sigma_sigma = .1,
					   mu_tau = 1,
					   sigma_tau = .1,
					   sigma = .15)
			 fit <- mod$sample(data = stan_data, chains = 1, iter_warmup = 1000, iter_sampling = 1000)
			 draws <- fit$draws(format = "df")

			 drawsf <- draws %>%
				 select(contains("f["))

			 fitted_y = apply(drawsf, MARGIN = 2, FUN = mean)

                         sigma <- mean(draws$tau)
			 sigmasd <- sd(draws$tau)

			 theta1 <- mean(draws$rho)
			 theta1sd <- sd(draws$rho)

			 theta2 <- mean(draws$sigma)
			 theta2sd <- sd(draws$sigma)

                         #ldist <- result$marginals.hyperpar$`Range for i` %>%
                         #  as.data.frame()
                         #theta1 <- ldist$x[which.max(ldist$y)]
                         
                         #ldist <- result$marginals.hyperpar$`Stdev for i` %>%
                         #  as.data.frame()
                         #theta2 <- ldist$x[which.max(ldist$y)]
                         #message("Finished task ", ind)
			 data.frame(xc = xc, yc = yc, zc = zc, 
				    sigma = sigma, sigmasd = sigmasd,
                                    theta1 = theta1, theta1sd = theta1sd,
				    theta2 = theta2, theta2sd = theta2sd, 
                                    tot_size = mean(abs(fitted_y)),
                                    time = time, y = fitted_y)

			 
                         
                       },
                       error = function(e) {
                         message(sprintf("Error in iteration %d: %s", ind, e$message))
                         return(NULL)   # or NULL, or some sentinel value
                       })
                     })
    param_means <- do.call(rbind, res) %>% 
      mutate(coords = paste(xc, yc, sep = "_"))
    
    param_means_all <- rbind(param_means_all, param_means)

}

saveRDS(param_means_all, "hmc_par_means.rds")
