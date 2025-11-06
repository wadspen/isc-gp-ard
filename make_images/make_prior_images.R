library(GPBayes)
library(tidyr)
library(ggplot2)
library(dplyr)


dbessel <- function(theta) {
  m <- length(theta)
  dbes <- c()
  for (i in 1:m) {
    if (abs(theta[i]) == 0) {dbes[i] <- Inf}
    else {
      dbes[i] <- (1/pi) *BesselK(0, abs(theta[i]))
    }
  }
  return(dbes)
}



khs <- function(x, t) {
  num <- exp(-(t^2/(2 * x^2)))
  den1 <- sqrt(2*pi*x^2)
  den2 <- pi * (1 + x^2)
  
  return((2 * num) / (den1 * den2))
}


dhs <- function(theta) {
  d <- c()
  for (i in 1:length(theta)) {
    # print(i)
    d[i] <- integrate(khs, 0, Inf, t = theta[i])$value
  }
  return(d)
}


t <- seq(-7.999, 8, length.out = 2001)

distc <- data.frame(t) %>% 
  mutate(BE = dbessel(t),
         HS = dhs(t)) %>% 
  pivot_longer(2:3, names_to = "dist", values_to = "p") %>% 
  ggplot() +
  geom_line(aes(x = t, y = p, colour = dist, linetype = dist), 
            size = 1.1) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) +
  coord_cartesian(ylim = c(0, 1.1), xlim = c(-3,3)) +
  labs(x = expression(theta),
       y = expression(pi(theta))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_blank())

distt <- data.frame(t) %>% 
  mutate(BE = dbessel(t),
         HS = dhs(t)) %>% 
  pivot_longer(2:3, names_to = "dist", values_to = "p") %>% 
  ggplot() +
  geom_line(aes(x = t, y = p, colour = dist, linetype = dist), 
            size = 1.1) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) +
  coord_cartesian(xlim = c(.5,4), ylim = c(0,.45)) +
  labs(x = expression(theta),
       y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.position = c(.7,.74),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.width = unit(2.2, "cm"))

library(gridExtra)
grid.arrange(distc, distt, ncol = 2)

dkapn <- function(kappa, tau = .5) {
  one <- 1/sqrt(2 * pi)
  two <- 1/(tau * kappa^(3/2) * sqrt(1 - kappa))
  three <- exp(-(1 - kappa) / (2 * kappa * tau^2))
  
  return(one * two * three)
}

dkaphs <- function(kappa, tau = .5) {
  den1 <- sqrt(kappa * (1 - kappa))
  den2 <- (1 + kappa * (tau^2 - 1))
  return(tau / (pi * den1 * den2))
}

ks <- seq(1e-8, 1, length.out = 1001)
ykap <- dkaphs(ks, tau = .2)
ykapn <- dkapn(ks, tau = .2)
plot(ykap ~ ks, type = "l", ylim = c(0,2))
lines(ykapn ~ ks, col = "red")

kappa.5 <- data.frame(ks) %>% 
  mutate(tau = .5) %>%
  mutate(HS = dkaphs(ks, tau),
         BE = dkapn(ks, tau)) %>%
  # dplyr:;select(tau) %>% 
  pivot_longer(3:4, names_to = "dist", values_to = "p") %>% 
  # filter(ks < .50) %>% 
  ggplot() +
  geom_line(aes(x = ks, y = p, colour = dist, linetype = dist),
            size = 1.2) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) +
  labs(x = expression(kappa),
       y = expression(paste("p(", kappa, "|", tau, " = .5)"))) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.position = c(.82,.15),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.width = unit(2.2, "cm"))


kappa1 <- data.frame(ks) %>% 
  mutate(tau = 1) %>%
  mutate(HS = dkaphs(ks, tau),
         BE = dkapn(ks, tau)) %>%
  # dplyr:;select(tau) %>% 
  pivot_longer(3:4, names_to = "dist", values_to = "p") %>% 
  # filter(ks < .50) %>% 
  ggplot() +
  geom_line(aes(x = ks, y = p, colour = dist, linetype = dist),
            size = 1.2) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  scale_colour_manual(values = c("#E69F00", "#0072B2")) +
  labs(x = expression(kappa),
       y = expression(paste("p(", kappa, "|", tau, " = 1)"))) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_blank())

grid.arrange(kappa1, kappa.5, ncol = 2)

  











