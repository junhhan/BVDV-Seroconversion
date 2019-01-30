# ABC-SMC for estimating within-herd BVDV transmission parameter
rm(list=ls ())
library(ggplot2)

setwd("D:/Temp/RC")
farm <- "F04"

# Varibles for ABC-SMC
max_i <- 5000 # Number of particles
tuner <- 0.1
list_beta <- NULL # Particles of beta
list_rho <- NULL # Particles of rho
list_tau <- NULL # Particles of tau

# Observed value
n_hfr <- 40
n_sample <- 15
n_test01 <- 9
n_test02 <- 15

d_mate_start <- 250
d_mate_end <- 306

p_conc <- 0.0411

d_test01 <- 590
d_test02 <- 735

n_sim_pos01 <- 0
n_sim_pos02 <- 0

p_sim01 <- 0
p_sim02 <- 0

# Error setting
epsilon <- c(6, 5, 4, 3, 2, 1, 0)
#p_obs01 <- n_test01/n_sample
#p_obs02 <- n_test02/n_sample
#e_list <- c(.01, .05, .35, .65, .95, .99)
#
#dist01 <- (p_obs01*(1-p_obs01)/n_sample)^0.5
#dist02 <- (p_obs02*(1-p_obs02)/n_sample)^0.5
#
#epsilon <- NULL
#for (i in 1:length(e_list)) {
#  epsilon <- c(epsilon, (0.5 * (dist01*qnorm(1-e_list[i]/2, 0, 1))^2 + (dist02*qnorm(1-e_list[i]/2, 0, 1))^2)^0.5)
#}
#rm(dist01, dist02, e_list)


# Sampling particles of the first population
bvdv_simulation <- function() {
  
  repeat {
    #beta <- runif(1, 0, 1)
    beta <- rbeta(1, 1.1785, 1.5355)
    if (beta >= 0.0000 & beta <= 1.0000) {break}	
  }
  
  repeat {
    #rho <- runif(1, 0, 1)
    rho <- rbeta(1, 1.1386, 7.7907)
    if (rho >= 0.0000 & rho <= 1.0000) {break}	
  }
  
  repeat {
    tau <- round(runif(1, 0, d_test01-1))
    if (tau >= 0 & tau <= d_test01-1) {break}	
  }
  
  .C("bvdv_simulation", as.double(beta), as.double(rho), as.integer(tau), as.integer(n_hfr), 
     as.integer(n_sample), as.double(p_conc), as.integer(d_mate_start), as.integer(d_mate_end), 
     as.integer(d_test01), as.integer(d_test02), as.integer(n_sim_pos01), as.integer(n_sim_pos02), 
     as.double(p_sim01), as.double(p_sim02))
}

dyn.load("bvdv_simulation.dll")
#dyn.unload("bvdv_simulation_sa01.dll")


# ABC-SMC for the first sequence
t_init <- proc.time()
t <- 1; i <- 1

while (i <= max_i) {
  result <- bvdv_simulation()
  
  #if (((result[[13]] - p_obs01)^2 + (result[[14]] - p_obs02)^2) <= epsilon[t]) {
  if (((result[[11]] - n_test01)^2 + (result[[12]] - n_test02)^2) <= epsilon[t]) {
    list_beta <- c(list_beta, result[[1]])
    list_rho <- c(list_rho, result[[2]])
    list_tau <- c(list_tau, result[[3]])
    
    i= i+1
  }
  rm(result)
}

wgt_beta <- rep(1/max_i, max_i)
wgt_rho <- rep(1/max_i, max_i)
wgt_tau <- rep(1/max_i, max_i)


# Sampling particles of the rest population
bvdv_simulation <- function() {
  
  repeat {
    beta_star <- sample(list_beta, size= 1, prob= wgt_beta)
    beta <- runif(1, min= beta_star - (tuner * (max(list_beta) - min(list_beta))), 
                  max= beta_star + (tuner * (max(list_beta) - min(list_beta))))
    if (beta >= 0.0000 & beta <= 1.0000) {break}	
  }
  
  repeat {
    rho_star <- sample(list_rho, size= 1, prob= wgt_rho)
    rho <- runif(1, min= rho_star - (tuner * (max(list_rho) - min(list_rho))), 
                 max= rho_star + (tuner * (max(list_rho) - min(list_rho))))
    if (rho >= 0.0000 & rho <= 1.0000) {break}	
  }
  
  repeat {
    tau_star <- sample(list_tau, size= 1, prob= wgt_tau)
    tau <- runif(1, min= tau_star - (tuner * (max(list_tau) - min(list_tau))), 
                 max= tau_star + (tuner * (max(list_tau) - min(list_tau))))
    if (tau >= 0 & tau <= d_test01-1) {break}	
  }
  
  .C("bvdv_simulation", as.double(beta), as.double(rho), as.integer(tau), as.integer(n_hfr), 
     as.integer(n_sample), as.double(p_conc), as.integer(d_mate_start), as.integer(d_mate_end), 
     as.integer(d_test01), as.integer(d_test02), as.integer(n_sim_pos01), as.integer(n_sim_pos02), 
     as.double(p_sim01), as.double(p_sim02))
}

for (t in 2:length(epsilon)) {
  
  tlist_beta <- NULL
  tlist_rho <- NULL
  tlist_tau <- NULL
  
  tlist_wgt_beta <- NULL
  tlist_wgt_rho <- NULL
  tlist_wgt_tau <- NULL
  
  list_sim01 <- NULL
  list_sim02 <- NULL
  
  i <- 1
  while (i <= max_i) {
    result <- bvdv_simulation()
    
    #if (((result[[13]] - p_obs01)^2 + (result[[14]] - p_obs02)^2) <= epsilon[t]) {
    if (((result[[11]] - n_test01)^2 + (result[[12]] - n_test02)^2) <= epsilon[t]) {
      
      tlist_beta <- c(tlist_beta, result[[1]])
      tlist_rho <- c(tlist_rho, result[[2]])
      tlist_tau <- c(tlist_tau, result[[3]])
      
      #tlist_wgt_beta <- c(tlist_wgt_beta, 1 / sum(wgt_beta))    
      tlist_wgt_beta <- c(tlist_wgt_beta, dbeta(result[[1]], 1.1785, 1.5355) / sum(wgt_beta))    
      
      #tlist_wgt_rho <- c(tlist_wgt_rho, 1 / sum(wgt_rho))    
      tlist_wgt_rho <- c(tlist_wgt_rho, dbeta(result[[2]], 1.1386, 7.7907) / sum(wgt_rho))    
      
      tlist_wgt_tau <- c(tlist_wgt_tau, 1 / sum(wgt_tau))
      
      list_sim01 <- c(list_sim01, result[[11]])
      list_sim02 <- c(list_sim02, result[[12]])
      
      i= i+1
    }
  }
  
  list_beta <- tlist_beta
  list_rho <- tlist_rho
  list_tau <- tlist_tau
  
  wgt_beta <- tlist_wgt_beta
  wgt_rho <- tlist_wgt_rho
  wgt_tau <- tlist_wgt_tau
}

proc.time() - t_init
dyn.unload("bvdv_simulation.dll")
#dyn.unload("bvdv_simulation_sa01.dll")


# Descriptive statistics of the result
quantile(list_beta, probs= c(0.05, 0.25, 0.5, 0.75, 0.95)); hist(list_beta)
quantile(list_rho, probs= c(0.05, 0.25, 0.5, 0.75, 0.95)); hist(list_rho)
quantile(list_tau, probs= c(0.05, 0.25, 0.5, 0.75, 0.95)); hist(list_tau)

quantile(list_sim01, probs= c(0.05, 0.25, 0.5, 0.75, 0.95)); n_test01; hist(list_sim01)
quantile(list_sim02, probs= c(0.05, 0.25, 0.5, 0.75, 0.95)); n_test02; hist(list_sim02)


# Beta:             5%       25%       50%       75%       95% 
## Original: 0.4510471 0.5635389 0.6370501 0.7173260 0.8273098
## 0.01 TI:  0.4310139 0.5455339 0.6297603 0.7130800 0.8286523
## 0.10 TI:  0.4707612 0.5694092 0.6400318 0.7175812 0.8245427
## Uni-Beta: 0.6273355 0.7466245 0.8236772 0.8910183 0.9665029
## Uni-Rho:  0.4668673 0.5762697 0.6581800 0.7392723 0.8495488

# Rho:               5%       25%         50%        75%        95% 
## Original: 0.01556799 0.02158055 0.02662769 0.03135023 0.03638902
## 0.01 TI:  0.01542385 0.02104635 0.02604048 0.03103224 0.03635030
## 0.10 TI:  0.01522382 0.02061453 0.02590789 0.03086326 0.03625318
## Uni-Beta: 0.01521305 0.02089753 0.02619890 0.03102464 0.03598309
## Uni-Rho:  0.01552894 0.02124073 0.02631076 0.03118635 0.03641848

# Tau:        5% 25% 50% 75% 95% 
## Original: 396 481 503 517 533
## 0.01 TI:  381 409 463 496 518
## 0.10 TI:  485 510 523 534 546
## Uni-Beta: 494 511 522 531 543
## Uni-Rho:  399 486 506 519 535

# Export beta, rho, tau for visualisation
write.table(data.frame(list_beta, list_rho, list_tau), paste("Para_", farm, ".csv", sep=""), sep= ",", row.names= F, col.names= F)


# Select the parameter sets with the most probable beta
density(list_beta)$x[which.max(density(list_beta)$y)] # 0.6141476
density(list_rho)$x[which.max(density(list_beta)$y)] # 0.03576595
density(list_tau)$x[which.max(density(list_beta)$y)] # 439
density(list_rho)$x[which.max(density(list_rho)$y)] # 0.02663139
density(list_tau)$x[which.max(density(list_tau)$y)] # 512


# Run a disease simulation using C and Import the result
dat_s <- read.table(paste("Prev_s_", farm, ".csv", sep= ""), header= F, sep= ",")
dat_r <- read.table(paste("Prev_r_", farm, ".csv", sep= ""), header= F, sep= ",")
dat_pi <- read.table(paste("Prev_pi_", farm, ".csv", sep= ""), header= F, sep= ",")

quantile(dat_s[, d_test02], probs= c(0.05, 0.5, 0.95)) # 0.03, 0.10, 0.55

Se <- 0.969; Sp <- 0.978
dat_s <- dat_s * (Se + Sp - 1) - Sp + 1
dat_s <- dat_s * 100
dat_s <- dat_s[, -(d_test02+1)]

dat_r <- dat_r * (Se + Sp - 1) - Sp + 1
dat_r <- dat_r * 100
dat_r <- dat_r[, -(d_test02+1)]

dat_pi <- dat_pi * (Se + Sp - 1) - Sp + 1
dat_pi <- dat_pi * 100
dat_pi <- dat_pi[, -(d_test02+1)]

prev_s_0.05 <- NULL; prev_s_0.50 <- NULL ; prev_s_0.95 <- NULL
prev_r_0.05 <- NULL; prev_r_0.50 <- NULL ; prev_r_0.95 <- NULL
prev_pi_0.05 <- NULL; prev_pi_0.50 <- NULL ; prev_pi_0.95 <- NULL

for (i in 1:d_test02) {
  prev_s_0.05 <- c(prev_s_0.05, quantile(dat_s[,i], probs= 0.05))
  prev_s_0.50 <- c(prev_s_0.50, quantile(dat_s[,i], probs= 0.50))
  prev_s_0.95 <- c(prev_s_0.95, quantile(dat_s[,i], probs= 0.95))
  
  prev_r_0.05 <- c(prev_r_0.05, quantile(dat_r[,i], probs= 0.05))
  prev_r_0.50 <- c(prev_r_0.50, quantile(dat_r[,i], probs= 0.50))
  prev_r_0.95 <- c(prev_r_0.95, quantile(dat_r[,i], probs= 0.95))
  
  prev_pi_0.05 <- c(prev_pi_0.05, quantile(dat_pi[,i], probs= 0.05))
  prev_pi_0.50 <- c(prev_pi_0.50, quantile(dat_pi[,i], probs= 0.50))
  prev_pi_0.95 <- c(prev_pi_0.95, quantile(dat_pi[,i], probs= 0.95))
}

prev <- data.frame("day"= c(1:d_test02), "group"= "S", "prev"= as.numeric(prev_s_0.50),
                   "prev_05"= as.numeric(prev_s_0.05), "prev_95"= as.numeric(prev_s_0.95))
prev <- rbind(prev, data.frame("day"= c(1:d_test02), "group"= "R", "prev"= as.numeric(prev_r_0.50), 
                               "prev_05"= as.numeric(prev_r_0.05), "prev_95"= as.numeric(prev_r_0.95)))
prev <- rbind(prev, data.frame("day"= c(1:d_test02), "group"= "PI", "prev"= as.numeric(prev_pi_0.50), 
                               "prev_05"= as.numeric(prev_pi_0.05), "prev_95"= as.numeric(prev_pi_0.95)))

ap01 <- 100 * (n_test01/n_sample)
ap02 <- 100 * (n_test02/n_sample)

ll01 <- ap01 - 100 * qnorm(1-0.05/2, 0, 1) * ((n_test01/n_sample) * (1 - n_test01/n_sample)/n_sample)^0.5
ul01 <- ap01 + 100 * qnorm(1-0.05/2, 0, 1) * ((n_test01/n_sample) * (1 - n_test01/n_sample)/n_sample)^0.5
ll02 <- ap02 - 100 * qnorm(1-0.05/2, 0, 1) * ((n_test02/n_sample) * (1 - n_test02/n_sample)/n_sample)^0.5
ul02 <- ap02 + 100 * qnorm(1-0.05/2, 0, 1) * ((n_test02/n_sample) * (1 - n_test02/n_sample)/n_sample)^0.5

ll01 <- ifelse(ll01 > 100, 100, ifelse(ll01 < 0, 0, ll01))
ul01 <- ifelse(ul01 > 100, 100, ifelse(ll01 < 0, 0, ul01))
ll02 <- ifelse(ll02 > 100, 100, ifelse(ul02 == 100, 100 * (1 - 3/n_sample), ifelse(ll01 < 0, 0, ll02)))
ul02 <- ifelse(ul02 >= 100, 100, ifelse(ll01 < 0, 0, ul02))

img <- ggplot(data= prev, aes(x= day, y= prev, group= group)) +
  geom_line(aes(color= group), size= 1.2) +
  scale_color_manual(values= c("S"= "blue", "R"= "green", "PI"= "red")) +
  geom_ribbon(data= prev, aes(x= day, ymax= prev_95, ymin= prev_05, fill= group), alpha= 0.3) +
  scale_fill_manual(values= c("S"= "lightblue", "R"= "lightgreen", "PI"= "darksalmon")) +
  geom_segment(aes(x= d_test01, xend= d_test01, y= ll01, yend= ul01)) +
  geom_segment(aes(x= d_test02, xend= d_test02, y= ll02, yend= ul02)) +
  geom_point(aes(x= d_test01, y= ap01), colour= "red", shape= 20, size= 3.5) +
  geom_point(aes(x= d_test02, y= ap02), colour= "red", shape= 20, size= 3.5) +
  labs(x= "Days from weaning", y= "Simulated AP of sero-positive heifers (%)") +
  theme(legend.position= "none",
        panel.background= element_rect(fill= "white", colour= "black"),
        axis.title.x= element_text(size= 15, colour= "black"),
        axis.title.y= element_text(size= 15, colour= "black"),
        axis.text.x = element_text(size= 12.5, colour= "black"),
        axis.text.y = element_text(size= 12.5, colour= "black"),
        strip.text = element_text(size= 15),
        aspect.ratio= 0.5) +
  scale_x_continuous(breaks= c(0, 200, 400, 600), limits= c(-5, d_test02+5)) +
  scale_y_continuous(breaks= c(0, 25, 50, 75, 100), limits= c(-5, 105)) +
  geom_vline(xintercept= c(0, 200, 400, 600), linetype= "solid", colour= "grey84", alpha= 0.5) +
  geom_hline(yintercept= c(0, 25, 50, 75, 100), linetype= "solid", colour= "grey84", alpha= 0.5)
img
ggsave(filename= paste("Prev_", farm, ".tiff", sep= ""), device= "tiff", dpi= 500)
