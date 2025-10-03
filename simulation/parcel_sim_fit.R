.libPaths("~/rlibs")
source("./sim_neuro_data.R")
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

a <- max(unique(hdr_sim$x))
b <- max(unique(hdr_sim$y))
parc_df <- expand.grid(x = 1:a, y = 1:b)

# Break into 3 bins along y, 4 bins along x
parc_df <- parc_df %>%
  mutate(
    xbin = cut(x, breaks = 6, labels = FALSE),
    ybin = cut(y, breaks = 5, labels = FALSE),
    roi = (ybin - 1) * 6 + xbin   # region index 1–12
  ) %>% 
  dplyr::select(-xbin, -ybin)

x_breaks <- quantile(parc_df$x, probs = seq(0, 1, length.out = 5))[-c(1,5)]
y_breaks <- quantile(parc_df$y, probs = seq(0, 1, length.out = 4))[-c(1,3)]

parc_df %>% 
  # filter(roi == 13) %>% 
ggplot( aes(x, y, fill = factor(roi))) +
  geom_tile(color = "grey90") +
  # geom_vline(xintercept = x_breaks + 0.5, color = "black") +
  # geom_hline(yintercept = y_breaks + 0.5, color = "black") +
  coord_equal() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Grid split into 12 regions (3x4)")

data_loc <- "../../dme_files/"


mod <- cmdstan_model("../stan_models/spat_nngp_hs.stan", 
                     include_paths = gptools_include_path())
hdr_df_ab <- hdr_sim %>%
  group_by(subject, x, y) %>%
  mutate(hdr = as.numeric(scale(hdr, center = TRUE, scale = TRUE))) %>%
  ungroup() %>%
  filter(!is.na(hdr)) %>% 
  filter(time < 51)

  
  hdr_df_pb <- hdr_df_ab %>% 
    mutate(voxel = paste(x, y, sep = "_")) %>% 
    left_join(parc_df, by = c("x", "y"))
  
  rois <- unique(hdr_df_pb$roi) 
  
  
  plan(multisession, workers = min(length(rois) + 4, 120))  # Windows-friendly
  
  function(ind) {
    tryCatch({  
  
  
      edge_df <- hdr_df_pb %>%
        dplyr::select(x, y) %>% 
        unique() %>% 
        mutate(id = row_number(), voxel = paste(x, y, sep = "_"))
      
      edges_x <- edge_df %>% 
        inner_join(edge_df, by = "x", suffix = c(".a", ".b")) %>%
        filter(id.a < id.b) %>%
        dplyr::select(id.a, id.b)
      
      edges_y <- edge_df %>%
        inner_join(edge_df, by = "y", suffix = c(".a", ".b")) %>%
        filter(id.a < id.b) %>%
        dplyr::select(id.a, id.b)
      
      edges <- bind_rows(edges_x, edges_y) %>% distinct()
      
      edges <- edge_df %>%
        inner_join(edge_df, by = character(), suffix = c(".a", ".b")) %>%
        filter(id.a < id.b) %>%
        filter(abs(x.a - x.b) + abs(y.a - y.b) == 1) %>%
        dplyr::select(id.a, id.b)
      
      
      
      coords_int <- coords %>% 
        as.data.frame() %>% 
        mutate(voxel = paste(x, y, sep = "_")) %>% 
        filter(voxel %in% unique(hdr_df_pb$voxel))
      
      voxels <- unique(coords_int$voxel)
      
      
      hdr_df <- hdr_df_pb %>%
        filter(roi == ind) %>% 
        dplyr::select(subject, hdr, time, x, y, voxel)
      
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
      
      
    
      
      voxels <- unique(coords_int$voxel)
      
      n <- length(unique(hdr_df$time))
      x <- unique(hdr_df$time)
      time <- unique(hdr_df$time)
      un_vox <- edge_df$voxel
      
      
      
      
      y_dat_list <- list()
      for (i in 1:length(un_vox)) {
        y_dat <- hdr_df %>% 
          filter(voxel == un_vox[i]) %>% 
          pivot_wider(id_cols = time, names_from = subject, 
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
                        sigma_rho = .1/sqrt(prod(dim(y_arr)[1:3])),
                        mu_sigma = 0,
                        sigma_sigma = .01,
                        mu_tau = 1,
                        sigma_tau = .1,
                        r = 1,
                        nut = 3000,
                        nul = 3,
                        m = 1,
                        node1 = as.integer(edges$id.a),
                        node2 = as.integer(edges$id.b))
      
       fit <- mod$sample(data = stan_data, 
                         chains = 1, 
                         # parallel_chains = 4,
                         iter_warmup = 1000, 
                         iter_sampling = 1000)
                           
       draws <- fit$draws(format = "df")
       preds_df <- data.frame()
       for (i in 1:length(un_vox)) {
         drawszf <- draws %>% 
           dplyr::select(contains("f[") & contains(paste0(",",i,"]")))
         
         preds <- apply(drawszf, MARGIN = 2, FUN = mean)
         upp <- apply(drawszf, MARGIN = 2, FUN = quantile, probs = .975)
         # t <- 1:250
         
         # plot(preds~t, type = "l")
         
         pdf <- data.frame(pred = preds, time = unique(hdr_df$time), 
                           voxel = un_vox[i])
         
         preds_df <- rbind(pdf, preds_df)
         
       }
       
       rownames(preds_df) <- 1:nrow(preds_df)
       
       params <- draws %>% 
         dplyr::select(!contains("f[")) %>% 
         dplyr::select(contains("[")) %>% 
         mutate(draw = row_number())
       
       par_long <- params %>% 
         pivot_longer(-"draw", names_to = "param", values_to = "val") %>% 
         mutate(vox_num = as.numeric(str_extract(param, "\\d+"))) %>% 
         mutate(param = str_replace(param, "\\[\\d+\\]", ""))
       
       voxels <- data.frame(vox_num = 1:length(un_vox), voxel = un_vox)
       
       par_wide <- par_long %>% 
         left_join(voxels, by = "vox_num") %>% 
         pivot_wider(names_from = param, values_from = val) %>% 
         mutate(n = prod(dim(y_arr)[1:3])) 
       
       par_wide$tau_sigma_rho <- rep(draws$tau_sigma_rho, length(un_vox))
       
       param_sums <- par_wide %>% 
         rowwise() %>% 
         mutate(kappa = 1 / (1 + n * lambda^2 * tau_sigma_rho^2 * phi^2
                             * 1/(tau^2))) %>% 
         group_by(voxel) %>% 
         summarise(kappa = mean(kappa),
                   sigma = mean(sigma),
                   rho = mean(rho),
                   lambda = mean(lambda), 
                   tau = mean(tau),
                   alpha = mean(alpha), 
                   phi = mean(phi),
                   n = unique(n),
                   tau_sigma_rho = mean(tau_sigma_rho))
       
       
       
       
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
  


saveRDS(param_means_all, "./spat_parc_sim.rds")
