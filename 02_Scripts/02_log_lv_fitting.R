# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Logistic Lotka-Volterra fitting 

library(rstan)
library(readr)
library(dplyr)

# Reference --------------------------------------------------------
# Stan & Lotka-Volterra 
# https://canada1.discourse-cdn.com/flex030/uploads/mc_stan/original/2X/b/b6acad3d9f2cd8e7b925c90809f9fb4e4b005a4f.pdf 
 

# Data.frame --------------------------------------------------------
indv_gr <- read_tsv("01_RawData/individual_strains_growth_curves_filtered.tsv")
indv_gr <- as.data.frame(indv_gr)
View(indv_gr)


# Stan function --------------------------------------------------

loglv <- '
functions {
  vector lvfnc(real t,           // time
               vector z,         // state (OD in this particular case)
               real r,           // r
               real k) {         // k
               
    vector[3] dzdt;
    for (j in 1:3) {
      dzdt[j] = r * z[j] * (1 - z[j] / k);
    }
    return dzdt;
  }
}

data {
  int<lower=1> N;            
  real ts[N];                // time series [without the first time measure]
  vector<lower=0>[3] y0;     // observed initial state 
  real<lower=0> y[N, 3];     // population measures
}

parameters {
  real<lower=0> r;
  real<lower=0> k;
  vector<lower=0>[3] z0;     // initial state (estimation)
  real<lower=0> sigma;
}

transformed parameters {
  // using vectors instead of real[]
  vector[3] z[N] = ode_rk45(lvfnc, z0, 0.0, ts, r, k);
}

model {
  // Priors
  r ~ normal(0.3, 0.05);
  k ~ normal(0.8, 0.05);
  sigma ~ lognormal(0, 0.05); # low error 

  // Prior for the initial state based on the values
  z0 ~ normal(y0, 0.1);
  
  // Likelihood
  for (j in 1:3) {
    y[, j] ~ normal(z[, j], sigma); 
  }
}
'

# Corrections 
# Theta, x_i & x_r are no needed 
# vectors instead of real[] objects 

# Function testing ------------------------------------------------------------
writeLines(loglv, con = "loglv.stan")

# time series 
time <- seq(0, 18, by = 2)

# OD
rep1_CH23 <- subset(indv_gr, Cepa == "CH23" &
                      ((rep == 1 & temp == 30))) %>%
             arrange(temp) %>%
             pull(OD_real)

rep2_CH23 <- subset(indv_gr, Cepa == "CH23" &
                      ((rep == 2 & temp == 30))) %>%
             arrange(temp) %>%
             pull(OD_real)

rep3_CH23 <- subset(indv_gr, Cepa == "CH23" &
                      ((rep == 4 & temp == 30))) %>%
             arrange(temp) %>%
             pull(OD_real)

# OD [ column are the different replicas]
OD_CH23 <- cbind(
  # replica 1
  rep1_CH23,
  # replica 2
  rep2_CH23,
  # replica 3
  rep3_CH23
)

# initial values 
in_val <- c(rep1_CH23[1], rep2_CH23[1], rep3_CH23[1])
N_val <- length(time)

log_CH23_df <- list(
  
    N = 9,
    ts = time[-1],
    y0 = in_val,
    y = OD_CH23[-1,]
  
)

# Specifying initial values for the variables of interest 
# because we want logical growth and capacity values 
CH23_init <- function(){
  list(
    r = 0.3,
    k = 0.8, 
    z0 = as.array(c(0.001, 0.001, 0.001)),
    sigma = 0.1
  )
}

log_CH23fit <- stan(model_code = loglv, 
                    data = log_CH23_df, 
                    save_dso = FALSE, 
                    iter = 2000,  # iterations
                    chains = 4,   # n. chains 
                    init = CH23_init) # the initial values i want stan to start from 
print(log_CH23fit)
plot(log_CH23fit)

# CH23 - 30° results stan function --------------------------------------------
# r = 0.30 
# k = 0.8

# Table interpretation notes --------------------------------------------------
# r & k are the ones that i'm interested in 
# r - fitness 
# k - carrying capacity 
# mean 
# se_mean 
# sd - standard deviation 
# % -> quantiles, confidence intervals ?? these provide an estimate of range where 
# the true value is likely to be 
# n_eff -> number of effective samples / if it is too low the aproximations 
# may lack precision, we want this number to be specially high 
# Rhat - Do the mcmc chains get to the same value? - 1.01 is a good indicator
# if this values goes above from 1.01 or 1.05 then that means the value of interest (r, k or the initial one)
# are not realiable values and we cannot use them in the model 
# Values close to 1.0 indicate that the chains have successfully converged to the same distribution

# Pipeline -------------------------------------------------------------------
indv_gr

stan_ccfunct <- function (df){
  
  spps <- list(unique(df[]))
  
  for (i in 1:N){
    rep1_sp <- subset(df, spp_col == spp[i] &
                        ((reppscol == reppsnum & tempcol == tempnum))) %>%
      arrange(temp) %>%
      pull(OD_real)
    
    
  }
  
  
  
}


indv_gro <- read_tsv("01_RawData/modified_individual_strain_growth_curves_ord_replica - individual_strains_growth_curves_filtered.tsv")
indv_gro <- as.data.frame(indv_gro)

especies <- unique(indv_gro[, 1])
temperatura <- sort(unique(indv_gro[,"temp"]), decreasing = F)
repeticion <- unique(indv_gro[, "ord_replica"])

vector_xrepe <- list() 
p <- 1

for (m in 1:length(especies)) {
  for (o in 1:length(repeticion)) {
    
    df_complete <- subset(indv_gro, 
                      Cepa == especies[m] & 
                        ord_replica == repeticion[o] & 
                        temp %in% c(30, 37, 42))
    
    if (nrow(df_complete) > 0){
      df_filtered <- df_complete %>% 
        arrange(temp) %>%
        pull(OD_real)
      
      df_matrix <- matrix(df_filtered, ncol = 3)
      colnames(df_matrix) <- c("30", "37", "42")
      vector_xrepe[[p]] = df_matrix
      
      names(vector_xrepe)[p] <- paste0(especies[m], "_rep", repeticion[o])
      p <- p + 1
    }
  }
}

