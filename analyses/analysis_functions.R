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


get_atlas <- function(atlas_loc = "../atlases/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_3mm.nii.gz") {
  atlas <- readNifti(atlas_loc)
  
  atlas_df <- as.data.frame(expand.grid(
    x = seq_len(dim(atlas)[1]),
    y = seq_len(dim(atlas)[2]),
    z = seq_len(dim(atlas)[3])
  ))
  
  # Add the values
  atlas_df$roi <- as.vector(atlas)
  return(atlas_df)
}

make_subjects <- function(subs = 1:22) {
  subjects <- as.character(subs)
  subjects <- ifelse(nchar(subjects) == 1, paste0(0, subjects), subjects)
  return(subjects)
}

read_bold_data <- function(data_loc = data_loc <- "../../dme_files/",
                         fold2 = "/ses-01/func/",
                         file = "_ses-01_task-dme_run-01_bold/func_preproc/func_pp_filter_gsr_sm0.mni152.3mm.nii.gz",
                         mask_name = "_ses-01_task-dme_run-01_bold/func_seg/wm_mask.nii",
                         subs = 1:22, 
                         checker = FALSE,
                         filter = TRUE,
                         gsr = TRUE,
			 time_steps = 97) {
  
  if (checker == TRUE) {
  	fold2 <- str_replace(fold2, "dme_run-01", "checker")
        file <- str_replace(file, "dme_run-01", "checker")
	mask_name <- str_replace(mask_name, "dme_run-01", "checker")
  }
  
  if (filter == FALSE) {file <- str_replace(file, "filter", "nofilt")}
  if (gsr == FALSE) {file <- str_replace(file, "gsr_", "")}
  subjects <- make_subjects(subs)
  sub_names <- paste0("sub-", subjects)
  
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
    all_bold[[i]] <- bold[,,,1:time_steps]
    
  }
  
  return(all_bold)
}


get_shared_coords <- function(all_bold, subs = 1:22, filtz = 27) {
  subjects <- make_subjects(subs)
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
    filter(dim3 %in% filtz) %>%
    dplyr::select(-n_groups)
  
  colnames(coords_int_full) <- c("x", "y", "z")
  return(coords_int_full)
}

get_hdr_df_pb <- function(all_bold, atlas_df, max_time = 50, subs = 1:22,
                          cross_sec = 27, coords_int_full) {
  coords_int_full <- coords_int_full %>% 
  subjects <- make_subjects(subs)
  sub_names <- paste0("sub-", subjects)
  nsubs <- length(subjects)
  
  ab_df <- data.frame()
  ab_ls <- list()
  for (i in 1:length(all_bold)) {
    
    
      arr <- all_bold[[i]][,,cross_sec,]
      if (length(cross_sec) == 1) {  
        nx <- dim(arr)[1]
        ny <- dim(arr)[2]
        nt <- dim(arr)[3]
        
        df <- data.frame(participant_id = sub_names[i],
                         expand.grid(x = 1:nx, y = 1:ny),
                         matrix(aperm(arr, c(3,1,2)), nrow = nx*ny, ncol = nt, 
                                byrow = TRUE))
        
        
        # Name time columns nicely
        colnames(df)[-(1:4)] <- paste0("t", 1:nt)
        ab_ls[[i]] <- df
        
      } else if (length(cross_sec) > 1) {
        nx <- dim(arr)[1]
        ny <- dim(arr)[2]
        nz <- dim(arr)[3]
        nt <- dim(arr)[4]
        
        df <- data.frame(participant_id = sub_names[i],
                         expand.grid(x = 1:nx, y = 1:ny, z = 1:nz),
                         matrix(aperm(arr, c(4,1,2,3)), nrow = nx*ny*nz, 
                                ncol = nt, 
                                byrow = TRUE))
                         
                         
        colnames(df)[-(1:4)] <- paste0("t", 1:nt)
        df <- df %>% 
          mutate(voxel = paste(x, y, z, sep = "_"))
        ab_ls[[i]] <- df
      }
      
    }
    
    ab_df <- do.call(rbind, ab_ls)
    start <- ifelse(length(cross_sec) == 1, 4, 5)
    hdr_df_ab <- ab_df %>% 
      pivot_longer(start:ncol(ab_df), names_to = "time", values_to = "hdr") %>% 
      mutate(time = as.numeric(str_replace(time, "t", ""))) %>% 
      filter(time <= max_time)
    
    hdr_df_ab$hdr <- as.numeric(hdr_df_ab$hdr)
    
    hdr_df_ab <- hdr_df_ab %>%
      group_by(participant_id, x, y) %>%
      mutate(hdr = as.numeric(scale(hdr, center = TRUE, scale = TRUE))) %>%
      ungroup() %>%
      filter(!is.na(hdr)) %>% 
      mutate(voxel = paste(x, y, sep = "_"))
    
    select_vox <- hdr_df_ab %>% 
      dplyr::select(participant_id, voxel) %>% 
      unique() %>% 
      group_by(voxel) %>% 
      summarise(n = n()) %>% 
      filter(n == nsubs)
    
    hdr_df_pb <- hdr_df_ab %>% 
      filter(voxel %in% select_vox$voxel) %>% 
      left_join(atlas_df %>% 
                  filter(z %in% cross_sec) %>% 
                  dplyr::select(-z), by = c("x", "y")) %>% 
      filter(roi != 0)
  
  return(hdr_df_pb)
}


get_stan_data <- function(hdr_df_pb, ind, sigma_rho_n) {
  edge_df <- hdr_df_pb %>%
    dplyr::select(x, y) %>% 
    unique() %>% 
    mutate(id = row_number(), voxel = paste(x, y, sep = "_"))
  
  
  edges <- edge_df %>%
    inner_join(edge_df, by = character(), suffix = c(".a", ".b")) %>%
    filter(id.a < id.b) %>%
    filter(abs(x.a - x.b) + abs(y.a - y.b) == 1) %>%
    dplyr::select(id.a, id.b)
  
  
  
  coords_int <- coords_int_full %>% 
    mutate(voxel = paste(x, y, sep = "_")) %>% 
    filter(voxel %in% unique(hdr_df_pb$voxel))
  
  voxels <- unique(coords_int$voxel)
  
  
  hdr_df <- hdr_df_pb %>%
    filter(roi == ind) %>% 
    dplyr::select(participant_id, hdr, time, x, y, voxel)
  
  edge_df <- hdr_df %>%
    dplyr::select(x, y) %>% 
    unique() %>% 
    mutate(id = row_number(), 
           voxel = paste(x, y, sep = "_"))
  
  
  edges <- edge_df %>%
    inner_join(edge_df, by = character(), suffix = c(".a", ".b")) %>%
    filter(id.a < id.b) %>%
    filter(abs(x.a - x.b) + abs(y.a - y.b) == 1) %>%
    dplyr::select(id.a, id.b)
  
  
  
  coords_int <- coords_int_full %>% 
    mutate(voxel = paste(x, y, sep = "_")) %>% 
    filter(voxel %in% unique(hdr_df$voxel))
  
  voxels <- unique(coords_int$voxel)
  
  n <- length(unique(hdr_df$time))
  x <- unique(hdr_df$time)
  time <- unique(hdr_df$time)
  un_vox <- edge_df$voxel
  
  
  
  
  y_dat_list <- list()
  for (i in 1:length(un_vox)) {
    y_dat <- hdr_df %>% 
      filter(voxel == un_vox[i]) %>% 
      pivot_wider(id_cols = time, names_from = participant_id, 
                  values_from = hdr) %>% 
      dplyr::select(-time) %>% 
      t() %>% 
      as.matrix()
    y_dat_list[[i]] <- y_dat
  }
  
  y_arr <- abind::abind(y_dat_list, along = 3)
  
  edge_index <- rbind(1:(ncol(y_dat) - 1), 2:ncol(y_dat))
  
  stan_data <- list(N = n, M = length(un_vox),
                    S = nrow(y_dat), 
                    n_edges = nrow(edges),
                    C = length(unique(hdr_df$voxel)),
                    y = y_arr,
                    x = matrix(unique(hdr_df$time), ncol = 1)*2,
                    dx = 1*2,
                    nf = n%/%2 + 1,
                    edge_index = edge_index,
                    mu_rho = log(10000),
                    sigma_rho = sigma_rho_n/((1 - sigma_rho_n)*sqrt(prod(dim(y_arr)[1:3]))),
                    mu_sigma = 0,
                    sigma_sigma = .01,
                    mu_tau = 1,
                    sigma_tau = .1,
                    r = 1,
                    nut = 1000,
                    nul = 1000,
                    m = 1,
                    node1 = as.integer(edges$id.a),
                    node2 = as.integer(edges$id.b))
  
  return(list(stan_data, hdr_df))
}



get_stan_data3D <- function(hdr_df_pb, ind, nu_t = 1000, nu_l = 1000, sigma_rho_n) {
  
  edge_df <- hdr_df_pb %>%
    dplyr::select(x, y, z) %>% 
    unique() %>% 
    mutate(id = row_number(),
           voxel = paste(x, y, z, sep = "_"))
  
  edges <- edge_df %>%
    inner_join(edge_df, by = character(), suffix = c(".a", ".b")) %>%
    filter(id.a < id.b) %>%
    filter(abs(x.a - x.b) + abs(y.a - y.b) + abs(z.a - z.b) == 1) %>%
    dplyr::select(id.a, id.b)
  
  
  
  coords_int <- coords_int_full %>% 
    mutate(voxel = paste(x, y, x, sep = "_")) %>% 
    filter(voxel %in% unique(hdr_df_pb$voxel))
  
  voxels <- unique(coords_int$voxel)
  
  
  hdr_df <- hdr_df_pb %>%
    filter(roi == ind) %>% 
    dplyr::select(participant_id, hdr, time, x, y, z, voxel)
  
  
  n <- length(unique(hdr_df$time))
  x <- unique(hdr_df$time)
  time <- unique(hdr_df$time)
  un_vox <- edge_df$voxel
  
  
  
  
  y_dat_list <- list()
  for (i in 1:length(un_vox)) {
    y_dat <- hdr_df %>% 
      filter(voxel == un_vox[i]) %>% 
      pivot_wider(id_cols = time, names_from = participant_id, 
                  values_from = hdr) %>% 
      dplyr::select(-time) %>% 
      t() %>% 
      as.matrix()
    y_dat_list[[i]] <- y_dat
  }
  
  y_arr <- abind::abind(y_dat_list, along = 3)
  if (dim(y_arr)[1] == 1) {y_arr <- y_arr[1,,]}
  edge_index <- rbind(1:(ncol(y_dat) - 1), 2:ncol(y_dat))
  
  stan_data <- list(N = n, M = length(un_vox),
                    S = nrow(y_dat), 
                    n_edges = nrow(edges),
                    C = length(unique(hdr_df$voxel)),
                    y = y_arr,
                    x = matrix(unique(hdr_df$time), ncol = 1)*2,
                    dx = 1*2,
                    nf = n%/%2 + 1,
                    edge_index = edge_index,
                    mu_rho = log(10000),
                    sigma_rho = sigma_rho_n/((1 - sigma_rho_n) * sqrt(prod(dim(y_arr)))),
                    mu_sigma = 0,
                    sigma_sigma = .01,
                    mu_tau = 1,
                    sigma_tau = .1,
                    r = 1,
                    nut = nu_t,
                    nul = nu_l,
                    m = 1,
                    node1 = as.integer(edges$id.a),
                    node2 = as.integer(edges$id.b))
  
  return(list(stan_data, hdr_df))
}




get_results <- function(fit, stan_data) {
  draws <- fit$draws(format = "df")
  hdr_df <- stan_data[[2]]
  un_vox <- unique(hdr_df$voxel)
  preds_df <- data.frame()
  nsub <- length(unique(hdr_df$participant_id))
  ntime <- length(unique(hdr_df$time))
  nvox <- length(unique(hdr_df$voxel))
  # for (i in 1:length(un_vox)) {
  #   drawszf <- draws %>% 
  #     dplyr::select(contains("pred_res") &
  #              contains(paste0("[",i,",")))
  #   
  #   preds <- apply(drawszf, MARGIN = 2, 
  #                  FUN = median)
  #   upp <- apply(drawszf, MARGIN = 2, 
  #                FUN = quantile, probs = .975)
  #   low <- apply(drawszf, MARGIN = 2, 
  #                FUN = quantile, probs = .025)
  #   # t <- 1:250
  #   
  #   # plot(preds~t, type = "l")
  #   
  #   pdf <- data.frame(pred = preds, 
  #                     upper = upp,
  #                     lower = low,
  #                     time = unique(hdr_df$time), 
  #                     voxel = un_vox[i])
  #   
  #   preds_df <- rbind(pdf, preds_df)
  #   
  # }
  
  # for (i in 1:length(un_vox)) {
  drawszf <- draws %>% 
    dplyr::select(contains("f["))
  
  preds <- apply(drawszf, MARGIN = 2, 
                 FUN = median)
  upp <- apply(drawszf, MARGIN = 2, 
               FUN = quantile, probs = .975)
  low <- apply(drawszf, MARGIN = 2, 
               FUN = quantile, probs = .025)
  # t <- 1:250
  
  # plot(preds~t, type = "l")
  
  normfs <- apply(drawszf, MARGIN = 1, FUN = function(x) {sum(x^2)})
  
  preds_df <- data.frame(pred = preds, 
                         upper = upp,
                         lower = low,
                         time = unique(hdr_df$time)) %>% 
    slice(rep(1:n(), each = length(un_vox)))
  
  # preds_df <- rbind(pdf, preds_df)
  
  # }
  preds_df$voxel <- rep(un_vox, each = ntime)
  
  rownames(preds_df) <- 1:nrow(preds_df)
  
  params <- draws %>% 
    dplyr::select(!contains("f[") & !contains("pred_res")) %>% 
    dplyr::select(contains("[")) %>% 
    mutate(draw = row_number()) %>% 
    unique()
  
  par_long <- params %>% 
    pivot_longer(-"draw", names_to = "param", values_to = "val") %>% 
    mutate(vox_num = as.numeric(str_extract(param, "\\d+"))) %>% 
    mutate(param = str_replace(param, "\\[\\d+\\]", ""))
  
  voxels <- data.frame(vox_num = 1:length(un_vox), 
                       voxel = un_vox)
  
  par_wide <- par_long %>% 
    left_join(voxels, by = "vox_num") %>% 
    pivot_wider(names_from = param, values_from = val) %>% 
    mutate(n = nsub*ntime*nvox) 
  
  par_wide$tau_sigma_rho <- rep(draws$tau_sigma_rho, length(un_vox))
  par_wide$normf <- rep(normfs, length(un_vox))
  
  param_sums <- par_wide %>% 
    rowwise() %>% 
    mutate(kappa = 1 / (1 + n * lambda^2 * tau_sigma_rho^2
                        * 1/(tau^2))) %>% 
    mutate(kappa2 = 1 / (1 + nsub * ntime * lambda^2 * tau_sigma_rho^2
                         * 1/(tau^2))) %>% 
    mutate(kappa3 = 1 / (1 + nsub * ntime * 4 * lambda^2 * tau_sigma_rho^2
                         * 1/(tau^2))) %>% 
    mutate(kappagp1 = 1 / (1 + normf * lambda^2 * tau_sigma_rho^2
                           * 1/(tau^2))) %>% 
    mutate(kappagp2 = 1 / (1 + nsub * normf * lambda^2 * tau_sigma_rho^2
                           * 1/(tau^2))) %>% 
    mutate(kappagp = 1 / (1 + nsub * nvox * normf * lambda^2 * tau_sigma_rho^2
                          * 1/(tau^2))) %>% 
    group_by(voxel) %>% 
    summarise(kappa = mean(kappa),
              kappa2 = mean(kappa2),
              kappa3 = mean(kappa3),
              kappagp = mean(kappagp),
              kappagp1 = mean(kappagp1),
              kappagp2 = mean(kappagp2),
              beta = mean(beta),
              upp_beta = quantile(beta, 
                                  probs = 0.975),
              low_beta = quantile(beta,
                                  probs = 0.025),
              lambda = mean(lambda), 
              tau = mean(tau),
              alpha = mean(alpha), 
              phi = mean(phi),
              n = unique(n),
              tau_sigma_rho = mean(tau_sigma_rho))
  
  
  
  
  fin_preds <- preds_df %>% 
    left_join(param_sums, by = "voxel")
  
  return(fin_preds)
}























