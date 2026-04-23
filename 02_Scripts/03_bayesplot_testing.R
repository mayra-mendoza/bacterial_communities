# Author: Mayra Beatriz Mendoza Velazquez 
# Title: Rstan 

library (bayesplot)
library(ggplot2)
library(hexbin)

# To do list: 
# * bayesplot - pairs 
# * bayesplot - trace 
# * distribution graph 
# r x iter graph 

# Load the data (stan outputs) ---------------------------------------------------
stan_out_allbut <- readRDS("03_Output/rk_by_temp_and_strain_withoutCH29")
CH29_30_10 <- readRDS("03_Output/CH29_T30_10times")
CH29_30_6 <- readRDS("03_Output/CH29_T30_6times")
CH29_37_10 <- readRDS("03_Output/CH29_T37_10times")
CH29_37_6 <- readRDS("03_Output/CH29_T37_6times")
CH29_42_10 <- readRDS("03_Output/CH29_T42_10times")
CH29_42_6 <- readRDS("03_Output/CH29_T42_6times")

# Merge two lists (for posterior manipulation) -------------------------------------
CH29_allist <- list( CH29_T30_10t = CH29_30_10,
                     CH29_T30_6t = CH29_30_6, 
                     CH29_T37_10t = CH29_37_10,
                     CH29_T37_6t = CH29_37_6,
                     CH29_T42_10t = CH29_42_10,
                     CH29_T42_6t = CH29_42_6)

stan_toutputs <- append(stan_out_allbut, CH29_allist)

# TRACE GRAPHS -------------------------------------------------------------------
# Ref: https://mc-stan.org/rstan/reference/stanfit-method-plot.html

# mcmc_trace ----------------------------------------------------------------
pdf("03_Output/Traceplots_CCStrains.pdf", width = 15, height = 7)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]
  
  pst_cp <- as.array(stan_toutputs[[i]])
  
  color_scheme_set("mix-brightblue-gray")
  
  # adding title to each plot 
  f <- mcmc_trace(pst_cp, pars = c("r", "k", "sigma")) +
    labs(title = paste("Traceplot:", strain),
         subtitle = "Sampling (1500)")
  
  print(f)
}

dev.off()


# plot()----------------------------------------------
pdf("03_Output/Traceplots_plotfun.pdf", width = 15, height = 7)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]

  # adding title to each plot 
  x <- plot(stan_toutputs[[i]], plotfun = "trace", pars = c("r", "k", "sigma"), inc_warmup = TRUE) +
    labs(title = paste("Traceplot:", strain),
         subtitle = "Shaded area: Warmup / Unshaded area: Sampling")
  
  print(x)
}

dev.off()


# PAIRS GRAPHS -------------------------------------------------------------------

# pairs func -------------------------------------------------------------------
pdf("03_Output/Pairsplots_pairsfun.pdf", width = 15, height = 7)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]
  
  color_scheme_set("purple")
  # adding title to each plot 
  x <- pairs(stan_toutputs[[i]], pars = c("r", "k", "sigma"), gap = 0, 
             main = paste0("Pairs plot_", strain), pch = 16, cex = 0.5) 
  
  print(x)
}

dev.off()


# mcmc_pairs fun ----------------------------------------------------------------

pdf("03_Output/Pairsplots_mcmcfun.pdf", width = 15, height = 7 , onefile = TRUE)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]
  pst_cp2 <- as.array(stan_toutputs[[i]])
  
  color_scheme_set("brightblue")
  # adding title to each plot 
  x <- mcmc_pairs(pst_cp2, pars = c("r", "k", "sigma"),
                  diag_fun = "hist", off_diag_fun = "hex",
                  condition = pairs_condition(chains = list(1, 2:4)), 
                  grid_args = list(top = paste("Pairs plot:", strain))) 
  
  print(x)
}

dev.off()


# r & k values -----------------------------------------------------------------

lista_r <- lapply(names(stan_toutputs), function(nm) {
  # 
  
  draws <- as_draws_df(stan_toutputs[[nm]]) %>%
    select(r) %>% # select the r value in each spp 
    mutate(ID = nombre) # keep the name 
  return(draws)
})

# Bind all the rows in the object and make two columns one for each variable 
df_boxplot <- bind_rows(lista_r) %>% 
separate(ID, into = c("Strain", "Temperature"), sep = "_T") 

df_boxplot <- df_boxplot %>%
  mutate(

      # generate a new "sampling_times column 
      sampling_times = case_when(
      str_detect(Temperature, "_6t")  ~ 6, # for the CH29 special cases 
      str_detect(Temperature, "_10t") ~ 10,
      TRUE                           ~ 10
    ),
    # to clean the temperature column 
    Temperature = str_replace_all(Temperature, "_6t|_10t", "")
  )

df_boxplot$Temperature <- as.numeric(df_boxplot$Temperature)


