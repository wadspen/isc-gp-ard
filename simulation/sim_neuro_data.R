library(neuRosim)
# library(hrf)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

set.seed(14)
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
onsets <- list(onsets.N1, onsets.N2, onsets.F1, onsets.F2)
dur <- list(0, 0, 0, 0)








region.1A.center <- c(13,13)
region.1A.radius <- 4
region.1B.center <- c(40,18)
region.1B.radius <- 6
# region.1C.center <- c(10,45,24)
# region.1C.radius <- 3
# region.2.center <- c(15,16,31)
# region.2.radius <- 5
# region.3.center <- c(12,16,13)
# region.3.radius <- 5
coord.regions <- list(region.1A.center, region.1B.center)
radius.regions <- c(region.1A.radius,region.1B.radius)
onsets.regions <- list(onsets, onsets)
dur.regions <- list(dur, dur)


# region.1a.d <- list(160.46, 140.19, 200.16, 160.69)
# region.1b.d <- list(140.51, 120.71, 160.55, 120.44)
# region.1c.d <- list(120.53, 120.74, 140.02, 100.48)
# region.2.d <- list( -0.24, 10.29, 80.18, 160.24)
# region.3.d <- list(192.7, 50.04, 240.60, 50.83)

region.1a.d <- as.list(rnorm(4, .75, .5))
region.1b.d <- as.list(rnorm(4, .6, .5))
# region.1c.d <- as.list(rnorm(4, 1.5, .5))
# region.2.d <- as.list(rnorm(4, 1.5, .5))
# region.3.d <- as.list(rnorm(4, 1.5, .5))
effect <- list(region.1a.d,region.1b.d)





design <- simprepTemporal(regions=2,
                          onsets=onsets.regions, durations=dur.regions,
                          hrf="double-gamma", TR=TR, totaltime=total.time,
                          effectsize=effect)
spatial <- simprepSpatial(regions=2,
                          coord=coord.regions, radius=radius.regions,
                          form="sphere", fading=0.05)



subjects <- 1:22
sub_data <- list()

hdr_sim <- data.frame()
for (i in 1:length(subjects)) {
  sub_data[[i]] <- simVOLfmri(design=design, image=spatial,
                              base=0, SNR=.3, noise="spatial", type="rician",
                              rho.temp=c(0.142,0.108,0.084), rho.spat=0.4,
                              w=c(0.05,0.1,0.01,0.09,0.05,0.7), dim=c(53,63),
                              # template=baseline.bin,
                              spat="gaussRF")
  
  arr <- sub_data[[i]]
  nx <- dim(arr)[1]
  ny <- dim(arr)[2]
  nt <- dim(arr)[3]
  
  df <- data.frame(subject = subjects[i],
                   expand.grid(x = 1:nx, y = 1:ny),
                   matrix(aperm(arr, c(3,1,2)), nrow = nx*ny, 
                          ncol = nt, byrow = TRUE)
  )
  
  # Name time columns nicely
  colnames(df)[-(1:3)] <- paste0("t", 1:nt)
  
  df_long <- df %>% 
    pivot_longer(4:ncol(df), names_to = "time", values_to = "hdr") %>% 
    mutate(time = as.numeric(str_replace(time, "t", "")))
  
  hdr_sim <- rbind(hdr_sim, df_long)
  
  # Name time columns nicely
  colnames(df)[-(1:3)] <- paste0("t", 1:nt)
}

saveRDS(hdr_sim, "init_sims.rds")

coords <- df_long %>% 
  select(x,y) %>% 
  unique()






























