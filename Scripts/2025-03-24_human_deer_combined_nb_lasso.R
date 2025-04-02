

# update and load needed packages
source("./Scripts/package_load.R")


# Read in data and process for analysis -----------------------------------

# read in data list
# dat <- st_read(dsn = "./Data/GIS/hexagons_parksmerged_201801_202001",
#                layer = "hexagons_parksmerged_201801_202001")

dat <- readRDS("./Data/2025-03-24_full_data.rds")

# isolate just the site level covariates
sitecovs <- dat |>
  dplyr::select(pop_total,
                habitat,
                road_dens,
                com_m2,
                res_m2,
                contig,
                trail_dens,
                poly_area) |>
  units::drop_units() |>
  st_drop_geometry()

# do we need to do any transformations
# for(i in 1:ncol(sitecovs)){
#   hist(sitecovs[, i], main = colnames(sitecovs)[i])
# }
# 
# hist(log(sitecovs$com_m2))

# transform the variables needed
sitecovs <- sitecovs |>
  mutate(pop_log = log(pop_total + 1),
         trail_log = log(trail_dens + 1),
         res_log = log(res_m2 + 1),
         com_log = log(com_m2 + 1),
         log_area = log(poly_area))

# check collinearity
# we will use LASSO regression so not as big of a deal, but worth looking at
# no high correlation
#cor_table <- cor(sitecovs)
#write.csv(cor_table, "covs_cor_table.csv")

# scale site covs
sitecovs_scaled <- apply(sitecovs, 2, function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x)
})


# mean surface temp
temp <- dat |>
  dplyr::select(st_mo_1:st_mo_25) |>
  st_drop_geometry()

# scale temp
temp_scaled <- scale(temp)


# Set up model for human activity -----------------------------------------

# isolate human activity
y <- dat |>
  dplyr::select(Jan18_hum:Jan20_hum) |>
  st_drop_geometry()

# number of model coefficients
nbeta <- 7

## Set up model code -------------------------------------------------------

# zero inflated poisson with lasso regression
nb_lasso <- nimbleCode({
  
  ## Priors
  # prior on intercept
  b0 ~ dnorm(0, 100)
  # scaling parameter for double exponential priors
  # smaller mean more mass at zero
  # each beta gets the same lambda
  b.lambda ~ dgamma(0.001, 0.001)
  # create beta parameters for covariates
  for(i in 1:nbeta){
    # lasso regression
    b[i] ~ ddexp(0, b.lambda) 
  }
  # prior for dispersian parameter
  r ~ dunif(0, 50)
  
  
  ## Likelihood
  for(t in 1:T){
    for(i in 1:N){
      y[i,t] ~ dnegbin(p[i,t],r)
      p[i,t] <- r/(r+lambda[i,t]) 
      lambda[i,t] <- exp(log_poi[i] + b0 + b[1]*habitat[i] + b[2]*road_dens[i] + 
                           b[3]*com_building[i] + b[4]*res_building[i] +
                           b[5]*contig[i] + b[6]*temp[i,t] + b[7]*pop_log[i])
    } 
  }
  
})



## Model arguments ---------------------------------------------------------

# constants list
my.constants <- list(habitat = sitecovs_scaled[,"habitat"],
                     road_dens = sitecovs_scaled[,"road_dens"],
                     log_poi = log(dat$POI_count),
                     com_building = sitecovs_scaled[ ,"com_log"],
                     res_building = sitecovs_scaled[ ,"res_log"],
                     contig = sitecovs_scaled[ ,"contig"],
                     #trail_dens = sitecovs_scaled[ ,"trail_dens"],
                     temp = temp_scaled,
                     pop_log = sitecovs_scaled[,"pop_log"],
                     N = nrow(y),
                     T = ncol(y),
                     nbeta = nbeta)

# initial values list

initial.values <- function() list(b0 = rnorm(1), 
                                  b = rnorm(nbeta),
                                  b.lambda = rgamma(1,1,1),
                                  r = dunif(1, 0, 50))

# parameters to monitor
params <- c("b0", "b", "b.lambda", "r")


## Function to run the model in parallel -----------------------------------

run_MCMC_allcode <- function(model.code, data, params, constants, initial_values,
                             niter, nburnin, nthin, seed) {
  
  library(nimble)
  
  myModel <- nimbleModel(code = model.code,
                         data = data,
                         constants = constants,
                         inits = initial_values)
  
  CmyModel <- compileNimble(myModel)
  
  configModel <- configureMCMC(CmyModel, monitors = params, enableWAIC = FALSE)
  
  myMCMC <- buildMCMC(configModel)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, 
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     setSeed = seed,
                     WAIC = FALSE)
  
  return(results)
}


## Fit Model ---------------------------------------------------------------

initial_values <- initial.values()


cl <- makeCluster(4)
nbLasso_human <- parLapply(cl = cl, 
                     X = 1:4, 
                     fun = run_MCMC_allcode, 
                     data = list(y = as.matrix(y)),
                     model.code = nb_lasso,
                     constants = my.constants,
                     initial_values = initial_values,
                     params = params,
                     niter = 35000, 
                     nburnin = 10000, 
                     nthin = 2)

stopCluster(cl)



## Model diagnostics -------------------------------------------------------

str(nbLasso_human)

MCMCtrace(nbLasso_human)
plotAuto(nbLasso_human)
MCMCsummary(nbLasso_human)
MCMCplot(nbLasso_human)


# Deer Model --------------------------------------------------------------

# isolate the deer data
y_Deer <- dat |>
  dplyr::select(Jan18_deer:Jan20_deer) |>
  st_drop_geometry()


## Set up model code -------------------------------------------------------

# zero inflated poisson with lasso regression
nb_lasso <- nimbleCode({
  
  ## Priors
  # prior on intercept
  b0 ~ dnorm(0, 100)
  # scaling parameter for double exponential priors
  # smaller mean more mass at zero
  # each beta gets the same lambda
  b.lambda ~ dgamma(0.001, 0.001)
  # create beta parameters for covariates
  for(i in 1:nbeta){
    # lasso regression
    b[i] ~ ddexp(0, b.lambda) 
  }
  # prior for dispersian parameter
  r ~ dunif(0, 50)
  
  
  ## Likelihood
  for(t in 1:T){
    for(i in 1:N){
      y[i,t] ~ dnegbin(p[i,t],r)
      p[i,t] <- r/(r+lambda[i,t]) 
      lambda[i,t] <- exp(log_area[i] + b0 + b[1]*habitat[i] + b[2]*road_dens[i] + 
                           b[3]*com_building[i] + b[4]*res_building[i] +
                           b[5]*contig[i] + b[6]*temp[i,t] + b[7]*pop_log[i])
    } 
  }
  
})



## Model arguments ---------------------------------------------------------

# constants list
my.constants <- list(habitat = sitecovs_scaled[,"habitat"],
                     road_dens = sitecovs_scaled[,"road_dens"],
                     area = sitecovs_scaled[ ,"poly_area"],
                     log_area = sitecovs[ ,"log_area"],
                     com_building = sitecovs_scaled[ ,"com_log"],
                     res_building = sitecovs_scaled[ ,"res_log"],
                     contig = sitecovs_scaled[ ,"contig"],
                     #trail_dens = sitecovs_scaled[ ,"trail_dens"],
                     temp = temp_scaled,
                     pop_log = sitecovs_scaled[,"pop_log"],
                     N = nrow(y_Deer),
                     T = ncol(y_Deer),
                     nbeta = nbeta)

# initial values list

initial.values <- function() list(b0 = rnorm(1), 
                                  b = rnorm(nbeta),
                                  b.lambda = rgamma(1,1,1),
                                  r = dunif(1, 0, 50))

# parameters to monitor
params <- c("b0", "b", "b.lambda", "r")


## Function to run the model in parallel -----------------------------------

run_MCMC_allcode <- function(model.code, data, params, constants, initial_values,
                             niter, nburnin, nthin, seed) {
  
  library(nimble)
  
  myModel <- nimbleModel(code = model.code,
                         data = data,
                         constants = constants,
                         inits = initial_values)
  
  CmyModel <- compileNimble(myModel)
  
  configModel <- configureMCMC(CmyModel, monitors = params, enableWAIC = FALSE)
  
  myMCMC <- buildMCMC(configModel)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, 
                     niter = niter,
                     nburnin = nburnin,
                     thin = nthin,
                     setSeed = seed,
                     WAIC = FALSE)
  
  return(results)
}


## Fit Model ---------------------------------------------------------------

initial_values <- initial.values()


cl <- makeCluster(4)
nbLasso_deer <- parLapply(cl = cl, 
                     X = 1:4, 
                     fun = run_MCMC_allcode, 
                     data = list(y = as.matrix(y_Deer)),
                     model.code = nb_lasso,
                     constants = my.constants,
                     initial_values = initial_values,
                     params = params,
                     niter = 35000, 
                     nburnin = 10000, 
                     nthin = 2)

stopCluster(cl)


## Model diagnostics -------------------------------------------------------

MCMCtrace(nbLasso_deer)
plotAuto(nbLasso_deer)
MCMCsummary(nbLasso_deer)
MCMCplot(nbLasso_deer)


# Save Results ------------------------------------------------------------

saveRDS(nbLasso_deer, "./Results/2025-03-24_CombinedAnalysis_DeerModelResults.rds")
saveRDS(nbLasso_human, "./Results/2025-03-24_CombinedAnalysis_HumanModelResults.rds")

# Plotting deer human results combined ------------------------------------

# combine chains
deer_mat <- do.call("rbind", nbLasso_deer)

human_mat <- do.call("rbind", nbLasso_human)

# parameter names
param_names <- c('Available habitat', 'Road density', 'Commercial building area', 
                 'Residential building area', 'Connectedness',
                 'Temperature', "Population Density")

# start plotting

{png("./Results/2025-03-24_deer_human_combine_pop.png", width = 10, height = 8, 
     units = "in", res = 300)
  
  par(mar = c(4,2,4,0) + 0.1)
  par(oma = c(0,10,0,0))
  
  plot(1:10, ylim=c(0.75,nbeta+0.75), yaxt="n", xaxt="n", xlim=c(-1,1), 
       xlab="", ylab="",
       pch=16, col='white', bty="n")
  axis(1, at = seq(-1,1,0.5), labels = FALSE, cex.axis = 1.5)
  mtext(as.character(seq(-1,1,0.5)), side = 1, line = 0.5, 
        at = seq(-1,1,0.5), cex = 0.8)
  mtext("Model coefficient", side = 1, line = 1.75, font = 2, cex = 0.8)
  axis(2, at = 1:7, labels = FALSE, outer = FALSE, las = 2)
  for(i in 1:nbeta){
    mtext(param_names[i], side = 2, at = i, line = 0.75, las = 2)
  }
  title(main = "deer = pink, human = yellow")
  abline(v = 0, lwd = 2, lty = 2, col = "grey")
  
  # Create underlines for each parameter
  for(i in 1:nbeta){
    lines(x=c(-1,1),y=c(i,i), lwd = 1, col="lightgrey")
  }
  
  # Plot deer
  for(i in 1:nbeta){
    vioplot(deer_mat[,i], at = i, horizontal = TRUE,
            add = TRUE, side = "right", drawRect = FALSE, wex = 2,
            col = alpha("white", 0.001), border = "black", frame.plot = FALSE,
            h = 0.025)
    vioplot(deer_mat[,i][between(deer_mat[,i],
                                 quantile(deer_mat[,i], probs = 0.025), 
                                 quantile(deer_mat[,i], probs = 0.975))], 
            at = i, horizontal = TRUE, 
            add = TRUE, side = "right", drawRect = FALSE, wex = 2,
            col = alpha("#1f77b4", 0.7), border = "black", frame.plot = FALSE,
            h = 0.025)
  }
  
  # plot human
  for(i in 1:nbeta){
    vioplot(human_mat[,i], at = i, horizontal = TRUE,
            add = TRUE, side = "right", drawRect = FALSE, wex = 2,
            col = alpha("white", 0.001), border = "black", frame.plot = FALSE,
            h = 0.01)
    vioplot(human_mat[,i][between(human_mat[,i],
                                  quantile(human_mat[,i], probs = 0.025), 
                                  quantile(human_mat[,i], probs = 0.975))], 
            at = i, horizontal = TRUE, 
            add = TRUE, side = "right", drawRect = FALSE, wex = 2,
            col = alpha("#ff7f0e", 0.7), border = "black", frame.plot = FALSE,
            h = 0.01)
  }
  
  dev.off()}



