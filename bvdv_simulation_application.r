# Warning! This code is built for parallel computation. 
# For those who want to use the single core, all the command relevant to "parallel" package should be removed and modified.
# However, running this code using single core would cost ~240 hours depending on the system.

# To run "bvdv_simulation.c" in R, one should compile the .c file into .dll that can be recognised and imported to R.
# Compilation can be done by running "R CMD SHLIB bvdv_simulation.c" at either R terminal or MS-DOS prompt.
# To run "R CMD SHLIB", one may need to install "Rtools35.exe".

rm(list= ls())
library(parallel)
library(TeachingDemos)


detectCores(logical= T)
n_core <- 8 # It depends on the number of core available for the local machine.
cl <- makeCluster(n_core) # Setting up the cores to use for parallel computation.

# Setup seed for RNG: Assign different seeds for different cores.
set.seed(1)
clusterSetRNGStream(cl, iseed= 1)
clusterApply(cl, seq_along(cl), function(i) seed <<- i - 1)

# ABC-SMC for estimating within-herd BVDV transmission parameter
setwd("D:/Temp/Seroconversion")

# Observed value
n_hfr <- c(42, 90, 45, 26, 23, 40, 21, 20, 29) # Number of heifers per farm
init_age <- c(202, 257, 242, 233, 225, 207, 227, 221, 208) # Calibrated age of heifers when they were weaned
n_sample <- c(15, 14, 15, 15, 15, 14, 15, 15, 15) # Number of sampled heifers per farm
f_psm <- c(250, 230, 203, 236, 590, 250, 590, 236, 240) # Simulation day of first planned start date of mating
s_psm <- f_psm + 365 # Simulation day of second planned start date of mating
mate_period <- c(70, 42, 42, 84 ,84, 42, 84, 112, 84) # Mating period
d_test01 <- c(590, 549, 533, 579, 946, 625, 950, 597, 640) # Simulation day of the first test
d_test02 <- c(735, 690, 711, 731, 1097, 752, 1097, 736, 794) # Simulation day of the second test
n_test01 <- c(9, 11, 2, 13, 13, 7, 3, 12, 14) # Number of sero-positive heifers at the first test
n_test02 <- c(6, 1, 6, 2, 2, 4, 9, 3, 1) # Number of sero-positive heifers at the second test

clusterExport(cl, c("n_hfr", "init_age", "n_sample", "f_psm", "s_psm", "mate_period", "d_test01", "d_test02", "n_test01", "n_test02"))


# Setting the variables for ABC-SMC
max_t <- 15 # Number of sequences
max_i <- 2000/n_core # Number of particles
tuner <- 0.68 # Tuning value for adjusting variance to sample in the next sequence

n_sim_pos01 <- rep(0, 9) # Dummy variable to store the simulation result (first test)
n_sim_pos02 <- rep(0, 9) # Dummy variable to store the simulation result (second test)

diff_01_list <- NULL # Dummy list to store the summary statistics
diff_02_list <- NULL # Dummy list to store the summary statistics

clusterExport(cl, c("max_t", "max_i", "tuner", "n_sim_pos01", "n_sim_pos02"))

# Model description for the practise and first sequence
bvdv_simulation <- function(beta, n_hfr, init_age, n_sample, f_psm, s_psm, 
                            mate_period, d_test01, d_test02, n_sim_pos01, n_sim_pos02) {
  
  repeat {rho <- rbeta(1, 1.1386, 7.7907);if (rho >= 0 & rho <= 1) {break}}
  repeat {mu <- runif(1, 0, 1);if (mu >= 0 & mu <= 1) {break}}
  repeat {tau <- round(runif(1, 1, d_test01)); if (tau >= 1 & tau <= d_test01) {break}}
  
  .C("bvdv_simulation", as.double(beta), as.double(rho), as.double(mu), as.integer(tau), as.integer(n_hfr), 
     as.integer(init_age), as.integer(n_sample), as.integer(f_psm), as.integer(s_psm), as.integer(mate_period), 
     as.integer(d_test01), as.integer(d_test02), as.integer(n_sim_pos01), as.integer(n_sim_pos02), as.integer(seed))
}

clusterExport(cl, c('bvdv_simulation'))


t_init <- proc.time()
# Estimating the initial tolerances
result <- 
  clusterEvalQ(cl, {
    
    dyn.load("D:/Temp/Seroconversion/bvdv_simulation.dll")
    i <- 0
    diff_01_list <- NULL; diff_02_list <- NULL
    while (i < max_i) {
      repeat {beta <- rbeta(1, 1.9668, 2.4501); if (beta >= 0 & beta <= 1) {break}}

      for (j in 1:9) {
        assign(paste("result0", j, sep= ""), 
               bvdv_simulation(beta, n_hfr[j], init_age[j], n_sample[j], f_psm[j], s_psm[j], 
                               mate_period[j], d_test01[j], d_test02[j], n_sim_pos01[j], n_sim_pos02[j], seed))
      }
      
      diff01 <- (
        (result01[[13]] - n_test01[1])^2 + (result02[[13]] - n_test01[2])^2 + (result03[[13]] - n_test01[3])^2 + 
          (result04[[13]] - n_test01[4])^2 + (result05[[13]] - n_test01[5])^2 + (result06[[13]] - n_test01[6])^2 +
          (result07[[13]] - n_test01[7])^2 + (result08[[13]] - n_test01[8])^2 + (result09[[13]] - n_test01[9])^2)^0.5
      
      diff02 <- (
        (result01[[14]] - n_test02[1])^2 + (result02[[14]] - n_test02[2])^2 + (result03[[14]] - n_test02[3])^2 + 
          (result04[[14]] - n_test02[4])^2 + (result05[[14]] - n_test02[5])^2 + (result06[[14]] - n_test02[6])^2 +
          (result07[[14]] - n_test02[7])^2 + (result08[[14]] - n_test02[8])^2 + (result09[[14]] - n_test02[9])^2)^0.5
      
      diff_01_list <- c(diff_01_list, diff01)
      diff_02_list <- c(diff_02_list, diff02)
      
      i <- i + 1
    }
    list(diff_01_list, diff_02_list)
  }
  )

for (i in 1:length(result)) {
  diff_01_list <- c(diff_01_list, result[[i]][[1]]) 
  diff_02_list <- c(diff_02_list, result[[i]][[2]])  
}

# Update tolerance
epsilon <- c(as.numeric(quantile(diff_01_list, 0.5)), as.numeric(quantile(diff_02_list, 0.5)))
diff_01_list <- NULL; diff_02_list <- NULL
print(epsilon); clusterExport(cl, "epsilon")


# Setting the dummy lists of weights and sds for sampling in the second sequence of ABC-SMC
beta_list <- rbeta(max_i*n_core, 1.9668, 2.4501)
mu_list <- NULL
rho_list <- NULL
tau_list <- NULL

for (j in 1:9) {
  mu_list[[j]] <- runif(max_i*n_core, 0, 1)
  rho_list[[j]] <- rbeta(max_i*n_core, 1.1386, 7.7907)
  tau_list[[j]] <- round(runif(max_i*n_core, 1, d_test01[j]))
}

beta_wgt <- rep(1/max_i*n_core, max_i*n_core)
rho_wgt<- list(rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), 
               rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), 
               rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core))
mu_wgt<- list(rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), 
              rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), 
              rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core))
tau_wgt<- list(rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), 
               rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), 
               rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core), rep(1/max_i*n_core, max_i*n_core))

# Lists of _tmp is required for updating information between sequences 
beta_list_tmp <- NULL; beta_wgt_tmp <- NULL

rho_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
rho_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
mu_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
mu_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
tau_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
tau_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)

rho_sd <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
mu_sd <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
tau_sd <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)

beta_ess <- NULL
rho_ess <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
mu_ess <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
tau_ess <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)

beta_ess <- c(beta_ess, 1/sum((beta_wgt/sum(beta_wgt))^2)) # Sum of squre-inverse of normalised weights
for (j in 1:9) {
  rho_ess[[j]] <- c(rho_ess[[j]], 1/sum((rho_wgt[[j]]/sum(rho_wgt[[j]]))^2))
  mu_ess[[j]] <- c(mu_ess[[j]], 1/sum((mu_wgt[[j]]/sum(mu_wgt[[j]]))^2))
  tau_ess[[j]] <- c(tau_ess[[j]], 1/sum((tau_wgt[[j]]/sum(tau_wgt[[j]]))^2))
}


# Model description for the rest sequences
bvdv_simulation <- function(beta, rho_list, rho_wgt, rho_sd, mu_list, mu_wgt, mu_sd, tau_list, tau_wgt, tau_sd, 
                            n_hfr, init_age, n_sample, f_psm, s_psm, mate_period, d_test01, d_test02, n_sim_pos01, n_sim_pos02, seed) {
  
  repeat {rho_star <- sample(rho_list, size= 1, prob= rho_wgt)
    rho <- rnorm(1, rho_star, rho_sd); if (rho >= 0 & rho <= 1) {break}}
  repeat {mu_star <- sample(mu_list, size= 1, prob= mu_wgt) 
    mu <- rnorm(1, mu_star, mu_sd); if (mu >= 0 & mu <= 1) {break}}
  repeat {tau_star <- sample(tau_list, size= 1, prob= tau_wgt) 
    tau <- round(rnorm(1, tau_star, tau_sd)); if (tau >= 1 & tau <= d_test01) {break}}
  
  .C("bvdv_simulation", as.double(beta), as.double(rho), as.double(mu), as.integer(tau), as.integer(n_hfr), as.integer(init_age), 
     as.integer(n_sample), as.integer(f_psm), as.integer(s_psm), as.integer(mate_period), 
     as.integer(d_test01), as.integer(d_test02), as.integer(n_sim_pos01), as.integer(n_sim_pos02), as.integer(seed))
}

clusterExport(cl, c('bvdv_simulation'))


# Sampling particles of the rest population
total_count <- 0
for (t in 2:max_t) {
  rm(result)
  clusterExport(cl, c("beta_list", "mu_list", "rho_list", "tau_list", "beta_wgt", "rho_wgt", "mu_wgt", "tau_wgt", "rho_sd", "mu_sd", "tau_sd"))
  
  # Estimate the SD of particles in the previous sequence  
  result <- 
    clusterEvalQ(cl, {
      
      dyn.load("D:/Temp/Seroconversion/bvdv_simulation.dll")
      
      diff_01_list <- NULL; diff_02_list <- NULL
      
      beta_sd <- (tuner * sum(((beta_list - mean(beta_list))^2) * beta_wgt) / sum(beta_wgt))^0.5
      
      for (j in 1:9) {
        rho_sd[[j]] <- (tuner * sum(((rho_list[[j]] - mean(rho_list[[j]]))^2) * rho_wgt[[j]]) / sum(rho_wgt[[j]]))^0.5
        mu_sd[[j]] <- (tuner * sum(((mu_list[[j]] - mean(mu_list[[j]]))^2) * mu_wgt[[j]]) / sum(mu_wgt[[j]]))^0.5
        tau_sd[[j]] <- (tuner * sum(((tau_list[[j]] - mean(tau_list[[j]]))^2) * tau_wgt[[j]]) / sum(tau_wgt[[j]]))^0.5
      }
      
      # Prepare the lists for storing intermediate results 
      beta_list_tmp <- NULL; beta_wgt_tmp <- NULL
      
      rho_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      rho_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      mu_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      mu_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      tau_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      tau_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      
      i <- 0; d <- 1
      while (i < max_i) {
        
        repeat {beta_star <- sample(beta_list, size= 1, prob= beta_wgt)
          beta <- rnorm(1, beta_star, beta_sd); if (beta >= 0 & beta <= 1) {break}}
      
        for (j in 1:9) {
          assign(paste("result0", j, sep= ""), 
                 bvdv_simulation(beta, rho_list[[j]], rho_wgt[[j]], rho_sd[[j]], mu_list[[j]], mu_wgt[[j]], mu_sd[[j]], 
                                 tau_list[[j]], tau_wgt[[j]], tau_sd[[j]], n_hfr[j], init_age[j], n_sample[j], f_psm[j], s_psm[j], 
                                 mate_period[j], d_test01[j], d_test02[j], n_sim_pos01[j], n_sim_pos02[j], seed))
        }
        
        diff01 <- (
          (result01[[13]] - n_test01[1])^2 + (result02[[13]] - n_test01[2])^2 + (result03[[13]] - n_test01[3])^2 + 
            (result04[[13]] - n_test01[4])^2 + (result05[[13]] - n_test01[5])^2 + (result06[[13]] - n_test01[6])^2 +
            (result07[[13]] - n_test01[7])^2 + (result08[[13]] - n_test01[8])^2 + (result09[[13]] - n_test01[9])^2)^0.5
        
        diff02 <- (
          (result01[[14]] - n_test02[1])^2 + (result02[[14]] - n_test02[2])^2 + (result03[[14]] - n_test02[3])^2 + 
            (result04[[14]] - n_test02[4])^2 + (result05[[14]] - n_test02[5])^2 + (result06[[14]] - n_test02[6])^2 +
            (result07[[14]] - n_test02[7])^2 + (result08[[14]] - n_test02[8])^2 + (result09[[14]] - n_test02[9])^2)^0.5
        
        if (diff01 <= epsilon[1] & diff02 <= epsilon[2]) {
          beta_list_tmp <- c(beta_list_tmp, beta)
          beta_wgt_tmp <- c(beta_wgt_tmp, dbeta(beta, 1.9668, 2.4501) / sum(beta_wgt * dnorm(beta, beta_list, beta_sd)))
          
          for (j in 1:9) {
            rho_list_tmp[[j]] <- c(rho_list_tmp[[j]], get(paste("result0", j, sep= ""))[[2]])
            rho_wgt_tmp[[j]] <- c(rho_wgt_tmp[[j]], 
                                  dbeta(get(paste("result0", j, sep= ""))[[2]], 1.1386, 7.7907)/
                                    sum(rho_wgt[[j]] * dnorm(get(paste("result0", j, sep= ""))[[2]], rho_list[[j]], rho_sd[[j]])))
            
            mu_list_tmp[[j]] <- c(mu_list_tmp[[j]], get(paste("result0", j, sep= ""))[[3]])
            mu_wgt_tmp[[j]] <- c(mu_wgt_tmp[[j]], 1/sum(mu_wgt[[j]] * dnorm(get(paste("result0", j, sep= ""))[[3]], mu_list[[j]], mu_sd[[j]])))
            
            tau_list_tmp[[j]] <- c(tau_list_tmp[[j]], get(paste("result0", j, sep= ""))[[4]])
            tau_wgt_tmp[[j]] <- c(tau_wgt_tmp[[j]], 1/sum(tau_wgt[[j]] * dnorm(get(paste("result0", j, sep= ""))[[4]], tau_list[[j]], tau_sd[[j]])))
          }
          
          diff_01_list <- c(diff_01_list, diff01)
          diff_02_list <- c(diff_02_list, diff02)
          
          i<- i + 1
        }
        d <- d + 1
      }
      beta_list <- beta_list_tmp
      beta_wgt <- beta_wgt_tmp
      
      rho_list <- rho_list_tmp; rho_wgt <- rho_wgt_tmp
      mu_list <- mu_list_tmp; mu_wgt <- mu_wgt_tmp
      tau_list <- tau_list_tmp; tau_wgt <- tau_wgt_tmp
      
      list(diff_01_list, diff_02_list, beta_list, beta_wgt, rho_list, rho_wgt, mu_list, mu_wgt, tau_list, tau_wgt, d)
    }
    
    )
  
  diff_01_list <- NULL; diff_02_list <- NULL
  for (i in 1:length(result)) {
    diff_01_list <- c(diff_01_list, result[[i]][[1]]) 
    diff_02_list <- c(diff_02_list, result[[i]][[2]])
    beta_list_tmp <- c(beta_list_tmp, result[[i]][[3]])
    beta_wgt_tmp <- c(beta_wgt_tmp, result[[i]][[4]])
    
    for (j in 1:9) {
      rho_list_tmp[[j]] <- c(rho_list_tmp[[j]], result[[i]][[5]][[j]])
      rho_wgt_tmp[[j]] <- c(rho_wgt_tmp[[j]], result[[i]][[6]][[j]])
      mu_list_tmp[[j]] <- c(mu_list_tmp[[j]], result[[i]][[7]][[j]])
      mu_wgt_tmp[[j]] <- c(mu_wgt_tmp[[j]], result[[i]][[8]][[j]])
      tau_list_tmp[[j]] <- c(tau_list_tmp[[j]], result[[i]][[9]][[j]])
      tau_wgt_tmp[[j]] <- c(tau_wgt_tmp[[j]], result[[i]][[10]][[j]])
      
    }
    
    total_count <- total_count + result[[i]][[11]]
  }
  
  # Effective sample sizes for each parameter
  beta_ess <- c(beta_ess, 1/sum((beta_wgt_tmp/sum(beta_wgt_tmp))^2)) # Sum of squre-inverse of normalised weights
  for (j in 1:9) {
    rho_ess[[j]] <- c(rho_ess[[j]], 1/sum((rho_wgt_tmp[[j]]/sum(rho_wgt_tmp[[j]]))^2))
    mu_ess[[j]] <- c(mu_ess[[j]], 1/sum((mu_wgt_tmp[[j]]/sum(mu_wgt_tmp[[j]]))^2))
    tau_ess[[j]] <- c(tau_ess[[j]], 1/sum((tau_wgt_tmp[[j]]/sum(tau_wgt_tmp[[j]]))^2))
  }
  
  # Update tolerance
  epsilon <- c(as.numeric(quantile(diff_01_list, 0.5)), as.numeric(quantile(diff_02_list, 0.5)))
  
  beta_result <- beta_list
  beta_list <- NULL # Particles of beta
  
  beta_list <- beta_list_tmp
  beta_wgt <- beta_wgt_tmp
  
  rho_list <- rho_list_tmp; rho_wgt <- rho_wgt_tmp
  mu_list <- mu_list_tmp; mu_wgt <- mu_wgt_tmp
  tau_list <- tau_list_tmp; tau_wgt <- tau_wgt_tmp
  
  beta_list_tmp <- NULL; beta_wgt_tmp <- NULL
  
  rho_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  rho_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  mu_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  mu_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  tau_list_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  tau_wgt_tmp <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  
  if (t != max_t) {
    assign(paste("beta", t, sep= ""), beta_list)
    assign(paste("rho", t, sep= ""), rho_list)
    assign(paste("mu", t, sep= ""), mu_list)
    assign(paste("tau", t, sep= ""), tau_list)
  }
  
  clusterExport(cl, "epsilon")
  print(t); print(epsilon); t_mid <- proc.time(); print(t_mid - t_init); Sys.sleep(0.01)
}

stopCluster(cl)
t_end <- proc.time()
t_end - t_init

save.image("bvdv_norm_parll.RData")
