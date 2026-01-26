# .libPaths("~/rlibs")
library(neuRosim)
# library(hrf)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(MASS)
library(tidyr)
library(purrr)
library(cmdstanr)
library(abind)
library(schoolmath)
library(PTHfftSURROGATES)
library(gptoolsStan)

is_stationary <- function(rho) {
  # Characteristic polynomial: 1 - rho1*z - rho2*z^2 - ...
  roots <- polyroot(c(1, -rho))
  all(Mod(roots) > 1)  # Stationary if all roots are outside unit circle
}

# Draw stationary AR(3) coefficients
draw_stationary_ar3 <- function() {
  repeat {
    rho <- mvrnorm(1, 
                   c(0.142,0.108,0.084),
                   diag(rep(.03, 3)))   # random candidate coefficients
    if (is_stationary(rho)) return(rho)
  }
}


my_shift <- function(ts, ns) {
  rsamp <- sample(unique(ns), 1)
  sts <- data.table::shift(ts, n = rsamp, type = "cyclic")
  return(sts)
}



sim_fmri_data <- function(nsubs = 22, dims = c(40, 22), seed = 16,
                          stim = "task", eff_size = 2.4) {
  set.seed(seed)
  nscan <- 250
  TR <- 2
  total.time <- nscan*TR
  nTime <- total.time/TR
  onsets.N1 <- c( 6.75, 15.75, 18.00, 27.00, 29.25, 31.50,
                  36.00, 42.75, 65.25, 74.25, 92.25, 112.50, 119.25,
                  123.75, 126.00, 137.25, 141.75, 144.00, 146.25, 155.25,
                  159.75, 162.00, 164.25, 204.75, 238.50)*TR
  onsets.N2 <- c(13.50, 40.50, 47.25, 56.25, 90.00, 94.50,
                 96.75, 135.00, 148.50, 184.50, 191.25, 202.50, 216.00,
                 234.00, 236.25)*TR#, 256.50, 261.00, 281.25, 290.25, 303.75,
  # 310.50, 319.50, 339.75, 342.00)*TR
  onsets.F1 <- c( 0.00, 2.25, 9.00, 11.25, 22.50, 45.00,
                  51.75, 60.75, 63.00, 76.50, 78.75, 85.50, 99.00,
                  101.25, 103.50, 117.00, 130.50, 150.75, 171.00, 189.00,
                  227.25)*TR#, 265.50, 283.50, 285.75, 288.00, 344.25)*TR
  onsets.F2 <- c(33.75, 49.50, 105.75, 153.00, 157.50, 168.75,
                 177.75, 180.00, 182.25, 198.00, 222.75, 240.75)*TR#, 254.25,
  # 267.75, 270.00, 274.40, 294.75, 299.25, 301.50, 315.00,
  # 317.25, 326.25, 333.00, 335.25, 337.50, 346.50)*TR
  if (stim == "task") {
    simple_onsets <- seq(3, 60, by = 16)*TR
    onsets <- list(simple_onsets)
    dur <- list(0)
  } else if (stim == "natural") {
    onsets <- list(onsets.N1, onsets.N2, onsets.F1, onsets.F2)
    dur <- list( 0, 0, 0)
    dur <- as.list(c(14, 6, 10, 24))
  }
  
  
  
  subjects <- 1:nsubs
  sub_data <- list()
  
  hdr_sim <- data.frame()
  for (i in 1:length(subjects)) {
    
    
    # region.1A.center <- c(8,8)
    region.1A.center <- dims/2
    region.1A.radius <- ceiling(rnorm(1, 2.5, .56))
    # region.1B.center <- c(28,16)
    # region.1B.radius <- ceiling(rnorm(1, 3.4, .5))
    # region.1C.center <- c(10,45,24)
    # region.1C.radius <- 3
    # region.2.center <- c(15,16,31)
    # region.2.radius <- 5
    # region.3.center <- c(12,16,13)
    # region.3.radius <- 5
    coord.regions <- list(region.1A.center)
    radius.regions <- c(region.1A.radius)
    onsets.regions <- list(onsets)
    # dur.regions <- list(dur)
    dur
    
    
    # region.1a.d <- list(160.46, 140.19, 200.16, 160.69)
    # region.1b.d <- list(140.51, 120.71, 160.55, 120.44)
    # region.1c.d <- list(120.53, 120.74, 140.02, 100.48)
    # region.2.d <- list( -0.24, 10.29, 80.18, 160.24)
    # region.3.d <- list(192.7, 50.04, 240.60, 50.83)
    
    # some quirk in the simulator make so that there are occasional issues if some effect sizes are 0. Hence the absolute value
    if (eff_size == 0) {
      effs <- abs(rnorm(length(onsets), eff_size, 1e-5))
    } else {
      effs <- abs(rnorm(length(onsets), eff_size, .34)) #change sample size to number of stimuli
    }
    
    region.1a.d <- as.list(effs)
    # region.1b.d <- as.list(rnorm(length(onsets), eff_size, .21))
    # region.1c.d <- as.list(rnorm(4, 1.5, .5))
    # region.2.d <- as.list(rnorm(4, 1.5, .5))
    # region.3.d <- as.list(rnorm(4, 1.5, .5))
    # effect <- list(region.1a.d)
    effect <- region.1a.d
    
    
    
    design <- simprepTemporal(
                              #regions=1,
                              onsets=onsets, durations=dur,
                              hrf="double-gamma", TR=TR, 
                              totaltime=total.time,
                              effectsize=effect)
    spatial <- simprepSpatial(regions=1,
                              coord=coord.regions, radius=radius.regions,
                              form="sphere", fading=.1)
    
    temp <- draw_stationary_ar3()
    sub_data[[i]] <- simVOLfmri(design=design, image=spatial,
                                base=0, SNR=.5, 
                                noise="mixture", 
                                type="rician",
                                rho.temp=temp,
                                rho.spat=abs(rnorm(1, 0.4,.1)),
                                w=c(0.05,0.1,0.01,0.09,0.05,0.7), 
                                dim=dims,
                                # template=baseline.bin,
                                spat="gaussRF")
    
    
    
    
    exp_design <- simprepTemporal(
                              #regions=1,
                              onsets=onsets, durations=dur,
                              hrf="double-gamma", TR=TR, 
                              totaltime=total.time,
                              effectsize=as.list(rep(eff_size, length(onsets))))
    
    signal_mask <- simVOLfmri(design = design, image = spatial,
                              base = 0, SNR = Inf, dim = dims, 
                              spat = "gaussRF")
    
    str_mask <- apply(signal_mask[,,], MARGIN = c(1,2), FUN = function(x) {
      mean(x)
    })
    
    mask <- apply(str_mask, MARGIN = 1:2, function(x) any(x != 0))
    
    arr <- sub_data[[i]]
    nx <- dim(arr)[1]
    ny <- dim(arr)[2]
    nt <- dim(arr)[3]
    
    df <- data.frame(subject = subjects[i],
                     expand.grid(x = 1:nx, y = 1:ny),
                     matrix(aperm(arr, c(3,1,2)), nrow = nx*ny, 
                            ncol = nt, byrow = TRUE)
    )
    
    mask_df <- data.frame(subject = subjects[i],
                          expand.grid(x = 1:nx, y = 1:ny),
                          active = matrix(aperm(mask, c(1,2)), nrow = nx*ny, 
                                          byrow = TRUE)
    )
    
    str_df <- data.frame(subject = subjects[i],
                         expand.grid(x = 1:nx, y = 1:ny),
                         msignal = matrix(aperm(str_mask, c(1,2)), 
                                          nrow = nx*ny, 
                                          byrow = TRUE)
    )
    
    # Name time columns nicely
    colnames(df)[-(1:3)] <- paste0("t", 1:nt)
    
    

    
    df_long <- df %>% 
      pivot_longer(4:ncol(df), names_to = "time", values_to = "hdr") %>% 
      mutate(time = as.numeric(str_replace(time, "t", ""))) %>% 
      left_join(mask_df, by = c("subject", "x", "y")) %>% 
      left_join(str_df, by = c("subject", "x", "y")) %>% 
      mutate(voxel = paste(x, y, sep = "_"))
    
    luv <- length(unique(df_long$voxel))
    
    spatial <- simprepSpatial(
      regions = 2,
      coord = list(c(1, 1, 1), c(1, 1, 1)),
      radius = c(1, 1),
      form = "sphere"
    )
    
    
    sim_hdr <- simVOLfmri(design=exp_design, image=spatial,
                          base=0, SNR=1, 
                          noise="none", 
                          type="rician",
                          rho.temp=0,
                          rho.spat=0,
                          w=c(1), 
                          dim=c(1,1,1),
                          # template=baseline.bin,
                          spat="gaussRF")
    
    hdr_convolved <- as.numeric(sim_hdr)
    hdr_convolved <- as.numeric(scale(hdr_convolved))
    df_long$hdr_conv <- rep(hdr_convolved, luv)
    hdr_sim <- rbind(hdr_sim, df_long)
    
    
    
    
    
    
  }
  
  
  # coords <- df_long %>% 
  #   dplyr::select(x,y) %>% 
  #   unique()
  
  
  return(list(sub_data, hdr_sim))
}

get_est_fun_dist <- function(hdr_df) {
  fun_dist <- hdr_df %>% 
    group_by(voxel, time) %>% 
    summarise(est_hdr = mean(hdr),
              hdr_conv = first(hdr_conv)) %>% 
    group_by(voxel) %>% 
    mutate(est_hdr = as.numeric(scale(est_hdr))) %>% 
    # dplyr::select(voxel, time, hdr_conv, est_hdr) %>% 
    unique() %>% 
    group_by(voxel) %>% 
    summarise(fun_dist = mean((est_hdr - hdr_conv)^2))
  
  return(fun_dist)
}


make_cor_rd_df <- function(hdr_df) {
  for_cor_df <- hdr_df %>% 
    dplyr::select(subject, voxel, time, hdr)

  df_pairs <- for_cor_df %>%
    rename(subject1 = subject, hdr1 = hdr) %>%
    inner_join(
      for_cor_df %>% rename(subject2 = subject, hdr2 = hdr),
      by = c("voxel", "time")
    ) %>% 
    filter(subject1 < subject2) %>% 
    arrange(voxel, subject1, subject2, time) %>% 
    mutate(pair = paste(subject1, subject2, sep = "_"))
  
  return(df_pairs)
}

est_zcor <- function(cor_rd_df, perm = FALSE) {
  tss <- cor_rd_df %>% 
    group_by(voxel, pair) %>% 
    summarise(corr = fishZ(cor(hdr1, hdr2))) %>% 
    group_by(voxel) %>% 
    summarise(pmcorr = median(corr))
  
  if (perm == FALSE) {
    tss$mcorr <- tss$pmcorr
    tss <- tss %>% 
      dplyr::select(-pmcorr)
  }
  
  
  return(tss)
}

perm_hdr <- function(cor_rd_df, method = "cyclic") {
  
  if (method == "cyclic") {
    hdr_perm <- cor_rd_df %>% 
      group_by(pair, voxel) %>% 
      mutate(n = n()) %>% 
      mutate(hdr2 = my_shift(hdr2, n))
  } else if (method == "shift") {
    hdr_perm <- cor_rd_df %>% 
      group_by(pair, voxel) %>% 
      mutate(hdr2 = phase_scramble2(hdr2))
  }
  
  return(hdr_perm)
}

get_perm_zcors <- function(cor_rd_df, K = 2000, perm = "cyclic") {
  all_perm_zcors <- list()
    for (i in 1:K) {
      hdr_perm <- perm_hdr(cor_rd_df, perm)
      perm_zcor <- est_zcor(hdr_perm, perm = TRUE)
      perm_zcor$rep <- i
      all_perm_zcors[[i]] <- perm_zcor
    }
  
  perm_zcors <- do.call(rbind, all_perm_zcors)
  return(perm_zcors)
}

get_act_res <- function(perm_zcors, hdr_df, zcor_ts) {
  sim_act <- perm_zcors %>% 
    left_join(zcor_ts, by = "voxel") %>% 
    mutate(dp = mcorr < pmcorr) %>% 
    group_by(voxel) %>% 
    summarise(pval = mean(dp)) %>% 
    ungroup() %>% 
    mutate(n = length(unique(voxel))) %>%
    mutate(active = ifelse(pval < .05/n, TRUE, FALSE))
  
  true_act <- hdr_df %>% 
    filter(time == 1) %>% 
    group_by(voxel) %>% 
    summarise(act_sum = sum(active))
  
  act_res <- true_act %>% 
    left_join(sim_act, by = "voxel")
  
  return(act_res)
}

get_perm_res <- function(perm_zcors, hdr_df, zcor_ts) {
  perm_res <- perm_zcors %>% 
    left_join(zcor_ts, by = "voxel") %>% 
    mutate(dp = mcorr < pmcorr) %>% 
    group_by(voxel) %>% 
    # summarise(pval = mean(dp)) %>% 
    mutate(pval = mean(dp)) %>% 
    ungroup() %>% 
    mutate(n = length(unique(voxel))) %>%
    mutate(active = ifelse(pval < .05/n, TRUE, FALSE))
  
  
  return(perm_res)
}

make_zcor_df <- function(df) {
  
  pairwise_corrs <- df %>%
    mutate(voxel = paste(x, y, sep = "_")) %>% 
    dplyr::select(subject, voxel, time, hdr) %>% 
    group_by(voxel) %>%
    summarise(
      corr_matrix = list({
        dat <- cur_data_all()
        
        # pivot so each subject is a column
        wide <- dat %>%
          pivot_wider(names_from = subject, values_from = hdr)
        
        # explicitly select only subject columns
        subject_cols <- setdiff(names(wide), c("time", "voxel"))
        mat <- as.data.frame(wide)[, subject_cols, drop = FALSE]
        
        # compute correlation matrix across subjects
        cor(mat, use = "pairwise.complete.obs")
      }),
      .groups = "drop"
    )
  
  
  isc_long <- map2_dfr(
    pairwise_corrs$voxel,
    pairwise_corrs$corr_matrix,
    function(loc_name, corr_mat) {
      as.data.frame(as.table(corr_mat)) %>%
        rename(n1 = Var1, n2 = Var2, corr = Freq) %>%
        mutate(voxel = loc_name)
    }
  )
  
  low_tr <- isc_long %>% 
    mutate(n1 = as.numeric(as.character(n1)),
           n2 = as.numeric(as.character(n2))) %>% 
    filter(n1 < n2) %>% 
    arrange(voxel, n1, n2) %>% 
    mutate(zcor = fishZ(corr)) %>% 
    dplyr::select(-corr)
  
  return(low_tr)
}






make_est_zcor <- function(sub_data, coords, dims = c(40, 22)) {
    subjects <- 1:length(sub_data)
    est_cor <- matrix(NA, nrow = dims[1], ncol = dims[2])
    est_zcor <- est_cor
    for (c in 1:nrow(coords)) {
      x <- coords[c, 1]; y <- coords[c, 2]
      sub_cor <- matrix(NA, nrow = length(subjects), ncol = length(subjects))  
      for (i in 1:length(subjects)) {
        for (j in i:length(subjects)) {
          sub_cor[i,j] <- as.numeric(cor(scale(sub_data[[i]][x,y,1:50]), 
                                         scale(sub_data[[j]][x,y,1:50])))
          sub_cor[j,i] <- sub_cor[i,j]
          
        }
      }
      diag(sub_cor) <- 0
      rho <- median(sub_cor, na.rm = TRUE)
      est_cor[x, y] <- rho
      est_zcor[x, y] <- fishZ(rho)
      
    }
    
    return(est_zcor)
}




cyclic_test <- function(sub_data, coords, dims = c(40, 22)) {
  subjects <- 1:length(sub_data)
  perm_est_zcor <- matrix(NA, nrow = dims[1], ncol = dims[2])
  library(schoolmath)
  for (c in 1:nrow(coords)) {
    x <- coords[c, 1]; y <- coords[c, 2]
    sub_cor <- matrix(NA, nrow = length(subjects), 
                      ncol = length(subjects))  
    for (i in subjects) {
      for (j in subjects) {
        ts1 <- as.numeric(scale(sub_data[[i]][x,y,1:50]))
        ts2 <- as.numeric(scale(sub_data[[j]][x,y,1:50]))
        ts1 <- data.table::shift(ts1,
                                 n = sample(length(ts1), 1),
                                 type = "cyclic")
        sub_cor[i,j] <- as.numeric(cor(ts1,
                                       ts2))
        sub_cor[j,i] <- sub_cor[i,j]
        
      }
    }
    diag(sub_cor) <- 0
    rho <- median(sub_cor, na.rm = TRUE)
    perm_est_zcor[x, y] <- fishZ(rho)
    
  }
  perm_est_zcor
}


rep_resamp_test <- function(resamp_func = NULL, est_zcor = NULL, 
                            sub_data, coords, dims = c(40, 22), 
                            K = 5000) {
  
  if (is.null(resamp_func)) {return()}
  
  corrs_list <- list()
  for (k in 1:K) {
    corrs_list[[k]] <- resamp_func(sub_data, coords, dims)
  }
  sim_arr <- simplify2array(corrs_list)
  sig_map <- matrix(NA, nrow = dims[1], ncol = dims[2])
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      sig_map[i,j] <- mean(abs(sim_arr[i,j,]) > est_zcor[i,j], na.rm = TRUE)
    }
  }
  
  df <- data.frame(
    y = rep(1:ncol(sig_map), each = nrow(sig_map)),
    x = rep(1:nrow(sig_map), times = ncol(sig_map)),
    pval = as.vector(sig_map)
  )
  
  return(df)
}






fishZ <- function(rho) {
  rhoZ <- .5 * log((1 + rho) / (1 - rho))
  return(rhoZ)
}






bootstrap <- function(data, n) {
  resampled_data <- lapply(1:n, function(i) {
    resample <- sample(data, replace = TRUE)
    # Perform desired operations on the resampled data, e.g., compute a statistic
    # and return the result
  })
  return(resampled_data)
}






################################################################
####################Sparse GP model functions###################
################################################################


get_stan_data <- function(hdr_df, nu_l = 1000, sigma_rho_n = 0.1) {
  hdr_df_pb <- hdr_df %>%
    group_by(subject, x, y) %>%
    mutate(hdr = as.numeric(scale(hdr, center = TRUE, scale = TRUE))) %>%
    ungroup() %>%
    filter(!is.na(hdr))
  
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
  
  
  hdr_df_red <- hdr_df_pb %>%
    dplyr::select(subject, hdr, time, x, y, voxel)
  
  edge_df <- hdr_df_red %>%
    dplyr::select(x, y) %>% 
    unique() %>% 
    mutate(id = row_number(), 
           voxel = paste(x, y, sep = "_"))
  
  
  edges <- edge_df %>%
    inner_join(edge_df, by = character(), suffix = c(".a", ".b")) %>%
    filter(id.a < id.b) %>%
    filter(abs(x.a - x.b) + abs(y.a - y.b) == 1) %>%
    dplyr::select(id.a, id.b)
  
  
  
  
  voxels <- unique(hdr_df_pb$voxel)
  
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
                    sigma_rho = sigma_rho_n/sqrt(prod(dim(y_arr)[1:3])),
                    mu_sigma = 0,
                    sigma_sigma = .01,
                    mu_tau = 1,
                    sigma_tau = .1,
                    r = 1,
                    nut = 1000,
                    nul = nu_l,
                    m = 1,
                    node1 = as.integer(edges$id.a),
                    node2 = as.integer(edges$id.b))
  
  return(stan_data)
}

get_sgp_fun_dist <- function(draws, hdr_df, un_vox) {
  predf_draws <- draws %>% 
    dplyr::select(contains("f["))
  
  predf <- predf_draws %>% 
    pivot_longer(everything(), names_to = "index", 
                 values_to = "draw") %>% 
    rowwise() %>% 
    extract(index, into = "time", regex = "(\\d+)") %>% 
    mutate(time = as.numeric(time))
    # rowwise() %>% 
    # extract(index, into = c("vox_num", "subject", "time"),
    #         regex = "\\[(\\d+),(\\d+),(\\d+)\\]", convert = TRUE)
  
  predf_sum <- predf %>% 
    # left_join(data.frame(voxel = un_vox, 
    #                      vox_num = 1:length(un_vox)),
    #           by = "vox_num") %>% 
    # group_by(voxel, time) %>% 
    group_by(time) %>% 
    summarise(pred = mean(draw),
              upper = quantile(draw, probs = 0.975),
              lower = quantile(draw, probs = 0.025))
  
  predf_dist <- predf_sum %>% 
    # group_by(voxel) %>% 
    # mutate(pred = as.numeric(scale(pred))) %>% 
    # left_join(hdr_df %>% 
    #             dplyr::select(voxel, time, hdr_conv), 
    #           by = c("voxel", "time")) %>% 
    left_join(hdr_df %>% 
                dplyr::select(time, hdr_conv), 
              by = c("time")) %>%
    dplyr::select(-upper, -lower) %>% 
    # group_by(voxel, time) %>% 
    group_by(time) %>% 
    summarise(pred = mean(pred),
              hdr_conv = mean(hdr_conv)
    ) %>% 
    # group_by(voxel) %>% 
    mutate(
      pred = as.numeric(scale(pred)),
      hdr_conv = as.numeric(scale(hdr_conv))) %>% 
    # group_by(voxel) %>% 
    summarise(fun_dist = mean((pred - hdr_conv)^2)
              
    )
  
  predf_dist <- data.frame(voxel = un_vox, fun_dist = predf_dist$fun_dist)
  
  mean_dist_tr <- hdr_df %>% 
    filter(active == TRUE) %>%
    group_by(subject) %>% 
    mutate(hdr = scale(hdr)) %>% 
    group_by(time) %>% 
    summarise(hdr = mean(hdr),
              hdr_conv = mean(hdr_conv)) %>% 
    ungroup() %>% 
    summarise(dist = mean(abs(scale(hdr) - scale(hdr_conv))^2))
  
  mean_dist <- hdr_df %>% 
    group_by(subject) %>% 
    mutate(hdr = scale(hdr)) %>% 
    group_by(time) %>% 
    summarise(hdr = mean(hdr),
              hdr_conv = mean(hdr_conv)) %>% 
    ungroup() %>% 
    summarise(dist = mean(abs(scale(hdr) - scale(hdr_conv))^2))
  
  predf_dist$mean_dist <- mean_dist$dist
  predf_dist$mean_tr <- mean_dist_tr$dist
  return(predf_dist)
}


get_gp_param_sums <- function(draws, un_vox, nsubjects, ntime) {
  nvox <- length(unique(un_vox))
  params <- draws %>% 
    dplyr::select(!contains("f") & !contains("pred")
                  & !contains("mean_res") & !contains("resid")) %>%  
    mutate(draw = row_number())
  
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
  
  par_long <- params %>% 
    pivot_longer(-"draw", names_to = "param", values_to = "val") %>% 
    mutate(vox_num = as.numeric(str_extract(param, "\\d+"))) %>% 
    mutate(param = str_replace(param, "\\[\\d+\\]", "")) 
  
  voxels <- data.frame(vox_num = 1:length(un_vox), voxel = un_vox)
  
  par_wide <- par_long %>% 
    left_join(voxels, by = "vox_num") %>% 
    pivot_wider(names_from = param, values_from = val) %>% 
    mutate(n = nsubjects*ntime*nvox) %>% #prod(dim(y_arr)[1:3])) %>% 
    filter(!is.na(voxel))
  
  par_wide$tau_sigma_rho <- rep(draws$tau_sigma_rho, length(un_vox))
  par_wide$rho <- rep(draws$rho, length(un_vox))
  par_wide$normf <- rep(normfs, length(un_vox))
  
  param_sums <- par_wide %>% 
    rowwise() %>% 
    mutate(kappa = 1 / (1 + n * lambda^2 * tau_sigma_rho^2 
                         * 1/(tau^2))) %>% 
    mutate(kappa2 = 1 / (1 + (nsubjects*ntime) * lambda^2 * tau_sigma_rho^2 
                         * 1/(tau^2))) %>% 
    mutate(kappa3 = 1 / (1 + (nsubjects*ntime*4) * lambda^2 * tau_sigma_rho^2 
                        * 1/(tau^2))) %>% 
    mutate(kappagp1 = 1 / (1 + normf * lambda^2 * tau_sigma_rho^2
                           * 1/(tau^2))) %>% 
    mutate(kappagp2 = 1 / (1 + nsubjects * normf * lambda^2 * tau_sigma_rho^2
                           * 1/(tau^2))) %>% 
    mutate(kappagp = 1 / (1 + nsubjects * nvox * normf * lambda^2 * tau_sigma_rho^2
                          * 1/(tau^2))) %>%
    group_by(voxel) %>% 
    summarise(kappa = mean(kappa),
              kappa2 = mean(kappa2),
              kappa3 = mean(kappa3),
              kappagp = mean(kappagp),
              kappagp1 = mean(kappagp1),
              kappagp2 = mean(kappagp2),
              rho = mean(rho),
              lambda = mean(lambda), 
              tau = mean(tau),
              alpha = mean(alpha),
              beta = mean(beta), 
              phi = mean(phi),
              n = unique(n),
              tau_sigma_rho = mean(tau_sigma_rho)) %>% 
    mutate(active = ifelse(kappa < .5, TRUE, FALSE))
  
  return(param_sums)
}



get_gp_act_res <- function(hdr_df, param_sums) {
  true_act <- hdr_df %>% 
    filter(time == 1) %>% 
    group_by(voxel) %>% 
    summarise(act_sum = sum(active))
  
  act_res <- true_act %>% 
    left_join(param_sums, by = "voxel")
  
  return(act_res)
}



get_blme_stan_data <- function(cor_rd_df, nsubjects) {
  pair_zcor_df <- cor_rd_df %>% 
    group_by(voxel, pair) %>% 
    summarise(zcor = fishZ(cor(hdr1, hdr2))) %>% 
    separate(pair, into = c("n1", "n2"), sep = "_", 
             convert = TRUE)
  
  nobs <- length(unique(cor_rd_df$pair))
  
  sub_space <- pair_zcor_df %>% 
    ungroup() %>% 
    dplyr::select(n1, n2) %>% 
    unique() %>% 
    as.matrix()
  
  W <- matrix(0, nrow = nobs, ncol = nobs)
  for (i in 1:nrow(sub_space)) {
    for (j in i:nrow(sub_space)) {
      nusubs <- length(unique(c(sub_space[i,], sub_space[j,])))
      W[i,j] <- ifelse(nusubs == 3, 1, 0)
      W[j,i] <- W[i,j]
    }
  }
  
  
  
  pz_wide <- pair_zcor_df %>% 
    mutate(comb = paste(n1, n2, sep = "_")) %>% 
    dplyr::select(zcor, voxel, comb, n1, n2) %>% 
    pivot_wider(names_from = voxel, values_from = zcor)
  
  ymat <- pz_wide %>% 
    dplyr::select(-comb, -n1, -n2) %>% 
    as.matrix()
  
  S <- nsubjects
  
  stan_data <- list(N = nrow(pair_zcor_df),
                    S = S,
                    R = ncol(ymat),
                    V = 1/2 * S * (S - 1), 
                    y = ymat,
                    subjects = 1:nsubjects,
                    n1_inds = pz_wide$n1,
                    n2_inds = pz_wide$n2,
                    n1s = pair_zcor_df$n1,
                    n2s = pair_zcor_df$n2,
                    ROIs = as.numeric(unique(factor(pair_zcor_df$voxel))),
                    rs = as.numeric(factor(pair_zcor_df$voxel)),
                    W = W,
                    I = diag(rep(1, ncol(W))))
  
  return(stan_data)
}


get_blme_par_res <- function(draws, un_vox) {
  par_res <- draws %>% 
    dplyr::select(contains("pi0")) %>% 
    pivot_longer(everything(), names_to = "param") %>% 
    mutate(vox_num = 
             as.numeric(str_extract(param, "(?<=\\[)\\d+(?=\\])"))) %>% 
    dplyr::select(-param) %>% 
    group_by(vox_num) %>% 
    summarise(upp = quantile(value, .975),
              low = quantile(value, .025),
              pred = mean(value)) %>%
    rowwise() %>% 
    mutate(active = low > 0) %>% 
    arrange(vox_num)
  
  par_res$voxel = un_vox
  
  return(par_res)
}







