#########################################################
# Function to fit model on given data using JAGS
# Author: Julie Sherman
# Email: juliesherman@ucsb.edu
#########################################################
mcmc_norm_run <- function(mydata, maxRake,
                          ni = 10000, nt = 2, 
                          nb = 3000, nc = 3,
                          rake.cluster = TRUE)
{
#Distributions for initial values
  QS_inits <- function() { list(
    "beta_p" = 0,#runif(1,0,5), 
    "alpha0" = sort(rnorm(QS_data$maxRakeScore, 0, 2)), 
    #"sd_p" = runif(1,0,5),
    "phi" = runif(1,0,1), 
    "psi" = runif(1,0,1),
    "sig_within" = runif(1,0,5),
    "bio_intercept" = runif(1,-10,7),
    "sig_site" = runif(1,0,5)
  )
  }
#name of JAGS model  
  model.name <- ifelse(rake.cluster,"QS_model_normal_nobio.jags",
                       "QS_model_normal_noRE.jags")
  
#parameters to return
  QS_params <- c("sd_p","alpha_p","beta_p","phi",
                 "psi","sig_within","sig_site",
                 "bio_intercept","exp.biomass.occquad",
                 "exp.biomass.allquad","prob.nodetect",
                 "tau.pred","random.pred",#"z.pred",
                 "biomass.pred") 
  if(!rake.cluster){
    QS_params<-QS_params[-c(1,13)]
  }

  all<-mydata
#leave one unit out for prediction
  set.seed(89)
  leftout<-sample(1:max(all$site),1)
  rakes<-all$rake[all$site==leftout]

#convert rake scores to matrix 
#each row is a site, each column indicator for rake score 0-5
#will convert to maxRake in the JAGS model
  rakeMatrix <- matrix(0,ncol=6, nrow=nrow(mydata))
  for(i in 1:nrow(mydata)){
    rakeMatrix[i, mydata$rake[i] + 1] <- 1
  }
  mydata$rakeMatrix <- rakeMatrix
  all<-mydata
  
#data and parameters to pass to JAGS
  QS_data <- list(
    y = all$rakeMatrix, 
    rake = all$rake, 
    n = nrow(all),
    maxRakeScore = maxRake, 
    siteIndex = all$site,
    z = all$z,
    occ.idx = which(as.logical(all$z)),
    abs.idx = which(!as.logical(all$z)),
    biomass = all$biomass, #divide by 0.5 to get g/m^2 units 
    numSites = max(all$site),
    numQuads = max(all$quad),
    site.occ = all %>% group_by(site) %>% 
      dplyr::summarize(site.occ = first(site.occ)) %>% 
      .$site.occ,
    site.occ.idx = which(as.logical(all %>% group_by(site) %>% 
                                      dplyr::summarize(site.occ = first(site.occ)) %>% 
                                      .$site.occ)),
 #   siteleftout = leftout,
    rakeleftout = rakes
  )
  set.seed(55442)
out<-as.mcmc.list(autorun.jags(model = model.name,
                   monitor = QS_params,
                   data = QS_data,
                   n.chains = nc,
                   inits = QS_inits,
                   startburnin = nb,
                   startsample = ni,
                   adapt = 3000,
                   thin = nt,
                   method = "parallel",
                   max.time="1d"))
  #call the JAGS model for fitting  
  #out1 <- jags(data = QS_data, 
  #             parameters.to.save = QS_params,  
  #             inits = QS_inits,
  #             model.file = model.name, 
  #             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
  #             working.directory = getwd())
#if hasn't converged, do a few more iterations
  #if(all(c(out1$BUGSoutput$summary["sd_p","Rhat"]<=1.1,
  #        out1$BUGSoutput$summary["beta_p","Rhat"]<=1.1,
  #        out1$BUGSoutput$summary["phi","Rhat"]<=1.1,
  #        out1$BUGSoutput$summary["psi","Rhat"]<=1.1,
  #        out1$BUGSoutput$summary["sig_within","Rhat"]<=1.1,
  #        out1$BUGSoutput$summary["sig_site","Rhat"]<=1.1,
  #        out1$BUGSoutput$summary["bio_intercept","Rhat"]<=1.1,
  #        out1$BUGSoutput$summary["alpha_p[1]","Rhat"]<=1.1,
  #        out1$BUGSoutput$summary["alpha_p[2]","Rhat"]<=1.1))){
  #  out<-out1
  #}else{
  #  out<- autojags(out1,n.thin=nt,n.iter=1000)
  #}
  print("finished model fitting")
#return model output and the site used for prediction  
  return(list(model = out))
}
