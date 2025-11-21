atlas_df <- get_atlas() %>% 
  filter(roi != 0) %>% 
  mutate(voxel = paste(x, y, z, sep = "_"))


checkers <- c(TRUE, FALSE); checks <- c("checker", "dme")
filters <- c(TRUE, FALSE); filts <- c("filter", "nofilt")
gsrs <- c(TRUE, FALSE); gs <- c("gsr", "nogsr")


for (h in 1:2) {
  for (j in 1:2) {
    for (k in 1:2) {
      
        all_bold <- read_bold_data(checker = checkers[h], 
                                   filter = filters[j], 
                                   gsr = gsrs[k])
        coords_int_full <- get_shared_coords(all_bold, 
                                             filtz = 1:dim(all_bold[[1]])[3])
        
        coords_int <- coords_int_full %>% 
          mutate(voxel = paste(x, y, z, sep = "_"))
        ab_ls <- list()
        for (i in 1:length(all_bold)) {
          
          
          arr <- all_bold[[i]][,,,]
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
            df <- df %>% 
              filter(voxel %in% unique(coords_int$voxel))
            ab_ls[[i]] <- df
           print(i)
        }
        
        all_ab_df <- do.call(rbind, ab_ls) %>% 
          left_join(atlas_df, 
                    by = c("x", "y", "z", "voxel")) %>% 
          filter(!is.na(roi))
        
        
        ab_list <- all_ab_df %>% 
          group_by(roi) %>% 
          group_split()
        
        
        
        
        for (i in 1:length(ab_list)) {
          ab_df <- ab_list[[i]] %>% 
            dplyr::select(-voxel, -roi)
          
          roi <- unique(ab_list[[i]]$roi)
          
          hdr_df_ab <- ab_df %>%
            pivot_longer(start:ncol(ab_df), names_to = "time", values_to = "hdr") %>%
            mutate(time = as.numeric(str_replace(time, "t", "")))

          hdr_df_ab$hdr <- as.numeric(hdr_df_ab$hdr)

          hdr_df_ab <- hdr_df_ab %>%
            group_by(participant_id, x, y) %>%
            mutate(hdr = as.numeric(scale(hdr, center = TRUE, scale = TRUE))) %>%
            ungroup() %>%
            filter(!is.na(hdr)) %>%
            mutate(voxel = paste(x, y, z, sep = "_"))

          select_vox <- hdr_df_ab %>%
            dplyr::select(participant_id, voxel) %>%
            unique() %>%
            group_by(voxel) %>%
            summarise(n = n()) %>%
            filter(n == nsubs)

          hdr_df_pb <- hdr_df_ab %>%
            filter(voxel %in% select_vox$voxel) %>%
            left_join(atlas_df, by = c("x", "y", "z", "voxel")) %>%
            filter(roi != 0)
          
          file_name <- paste0("../../data_by_roi/", checks[h], "/",
                 paste(filts[j], gs[k], sep = "_"), "/roi_", roi,
                 ".rds")
          
          saveRDS(hdr_df_pb, file_name)
        
        }
        
    }
  }
}       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        