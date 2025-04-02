# function to update and load needed libraries

package_load<-function(packages = NULL, quiet=TRUE, verbose=FALSE, 
                       warn.conflicts=FALSE){
  
  # download required packages if they're not already
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", 
                     quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, 
            warn.conflicts=warn.conflicts)
}

package_load(c('sf', 'terra', 'tidyverse', 'tidycensus', 'nimble', 'MCMCvis', 
               'parallel', 'foreach', 'doParallel', 'scales', 'mapview', 'units', 
               'landscapemetrics', 'runjags','rjags','coda', 'plotMCMC', 
               'vioplot', 'adehabitatHR', 'zoo', 'lubridate'))