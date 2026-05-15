#########################################################
# Chose parameter combos, run the model, and save
# Author: Julie Sherman
# Email: juliesherman@ucsb.edu
#########################################################
#note that with d=100, n.site = 100, n.quad = 5,
#each parameter combination takes about one hour to run on 32 cores
#recommend initial exploration of this script by setting
#d and maxiter to small numbers (e.g. 5 and 5)

#load packages
library(tidyverse)
library(runjags)
library(mcmc)
library(jmvcore)
library(boot)
library(parallel)
library(dplyr)
library(purrr)
library(tibble)
library(tictoc)
library(foreach)
library(doMC)
library(tgp)


#load functions
source("sim_func_normal.R")
source("mcmc_norm_run.R")

#MCMC settings
ni <- 10000   # num interations (kept)
nt <- 10       # thinning number
nb <- 3000    # burn in
nc <- 3       # num mcmc chains
(ni-nb)/n

d<-100 #number of datasets per parameter combo
n.quad<-5
n.site<-100

maxRake<-2

## Using latin hypercube sampling of parameter space 
## Parameter space search sampling alternative starts on line ~160
rect <- rbind(phi = c(0.1,1),               
              psi = c(0.1,1),                  
              mu = c(-5,6),                
              beta = c(0,15),                     
              alpha0 = c(-10,10),  
              alpha1   = c(-5,20),  
              sig.site  = c(0.01,5), 
              sig.within   = c(0.01,5),
              alpha2 = c(-5,20),
              alpha3 = c(-5,25),
              alpha4 = c(0,25))

#define the maximum number of combos you want to run
maxiter<-1
numran<-1

set.seed(444)
combos <- lhs(maxiter, rect)

#with latin  hypercube sampling
while(numran<maxiter+1){
  phi.site = combos[numran,1]
  psi.quad = combos[numran,2]
  mu = combos[numran,3]
  beta = combos[numran,4]
  alpha = sort(combos[numran,c(5,6,9,10,11)])
  sig.site = combos[numran, 7]
  sig.within = combos[numran, 8]
  
  #set file names to hold info on underlying parameters
  #to avoid special characters in file names, replace - with n, . with p
  outname<-paste0(gsub("-","n",
                       gsub(".","p",
                            paste("no_rc_lhs_sim_out_runjags", 
                                  n.site, n.quad, 
                                  round(phi.site,5), 
                                  round(psi.quad,5),
                                  round(sig.within,5),
                                  round(sig.site,5),
                                  round(mu,5), 
                                  round(beta,5),maxRake,
                                  round(alpha[1],5),round(alpha[2],5),
                                  0,#round(re.rake,5), 
                                  ni,nt,nb,nc,d,
                                  sep = "_"),fixed=TRUE),
                       fixed = TRUE),'.RData')
  outnamedata<-paste0(gsub("-","n",
                           gsub(".","p",
                                paste("no_rc_lhs_sim_data_runjags",
                                      n.site, n.quad, 
                                      round(phi.site,5), 
                                      round(psi.quad,5),
                                      round(sig.within,5),
                                      round(sig.site,5),
                                      round(mu,5), 
                                      round(beta,5),maxRake,
                                      round(alpha[1],5),round(alpha[2],5),
                                      0,#round(re.rake,5), 
                                      ni,nt,nb,nc,d,
                                      sep = "_"),fixed=TRUE),
                           fixed = TRUE),'.RData')
  
  
  #See if this parameter set has already been run, 
  #if so, load it. If not, run the model.
  if(isError(tryCatch(readRDS(outname), 
                      error = function(e) e))){
    datalist<-list()
    for (p  in 1:d){
      #simulate the data      
      mydata <- sim_norm_data(n.site, n.quad, phi.site, psi.quad, 
                              mu,sig.site,sig.within, maxRake, 
                              alpha, beta, 
                              rake.cluster = FALSE, re.rake)
      #type is random effects model      
      mydata <- bind_rows(mutate(mydata, type="not RE"))
      #site occupancy is unknown unless a positive rake score is recorded
      mydata$site.occ <- NA
      mydata <- mydata %>% group_by(site, type) %>% 
        mutate(site.occ = ifelse(sum(z) > 0, 1, NA)) %>%   ungroup()
      datalist[[p]]<- mydata
    }
    #save simulated data    
    saveRDS(datalist, file = outnamedata)
    tic()
    #initialize parallel workers    
    registerDoMC(detectCores()-1)
    getDoParWorkers()
    #fit the model for each dataset in parallel    
    out<-foreach(p=1:d, .inorder = FALSE,
                 .packages = c("runjags","dplyr"),
                 .verbose = TRUE,
                 .errorhandling = "pass") %dopar% 
      mcmc_norm_run(datalist[[p]], maxRake = maxRake, 
                    ni = ni, nt = nt, 
                    nb = nb, nc = nc,
                    rake.cluster = FALSE)
    toc()
    #save model ouput
    saveRDS(out,file = outname)
  } else {
    print("Already ran it!")
#    out<-readRDS(outname)
  }
  numran<-numran+1
  print(paste("Starting parameter combo ", numran))
}



##  with the "search" sampling
if(FALSE){
  
  #set bounds of parameter search
  #can adjust to search a subset of parameter space
  phimin<-0.1
  phimax<-1
  psimin<-0.15
  psimax<-1
  mumin<--5
  mumax<-6
  betamin<-0
  betamax<-15
  alpha0min<--10
  alpha0max<-10
  alpha1min<--5
  alpha1max<-20
  sig.sitemin<-0.01
  sig.sitemax<-5
  sig.withinmin<-0.001
  sig.withinmax<-5
  re.rakemin<-0.01
  re.rakemax<-5
  
  while(numran<maxiter){
    #d is the number of datasets per parameter combination
    d<-100
    #randomly resample every 5 to avoid local minima
    if(numran%%5 == 0){
      set.seed(numran+33233)
      #uniformly sample from parameter space  
      phi.site<-runif(1,phimin,phimax)
      psi.quad<-runif(1,psimin,psimax)
      mu <-runif(1,mumin,mumax)
      sig.within <-runif(1,sig.withinmin,sig.withinmax)
      sig.site <-runif(1,sig.sitemin,sig.sitemax)
      beta <-runif(1,betamin,betamax)
      alpha <-sort(c(runif(1,alpha0min,alpha0max),
                     runif(1,alpha1min,alpha1max)))
      re.rake <-runif(1,re.rakemin,re.rakemax)
      #initialize cost value, momentum term, perturbation term    
      co<-rep(100000, length(truth))
      Me<-0
      Ge<-0
    } 
    #set file names to hold info on underlying parameters
    #to avoid special characters in file names, replace - with n, . with p
    outname<-paste0(gsub("-","n",
                         gsub(".","p",
                              paste("rand_sim_out_nobio_ff", 
                                    n.site, n.quad, 
                                    round(phi.site,5), 
                                    round(psi.quad,5),
                                    round(sig.within,5),
                                    round(sig.site,5),
                                    round(mu,5), 
                                    round(beta,5),maxRake,
                                    round(alpha[1],5),round(alpha[2],5),
                                    0,round(re.rake,5), 
                                    ni,nt,nb,nc,d,
                                    sep = "_"),fixed=TRUE),
                         fixed = TRUE),'.RData')
    outnamedata<-paste0(gsub("-","n",
                             gsub(".","p",
                                  paste("sim_data_nobio_ff",
                                        n.site, n.quad, 
                                        round(phi.site,5), 
                                        round(psi.quad,5),
                                        round(sig.within,5),
                                        round(sig.site,5),
                                        round(mu,5), 
                                        round(beta,5),maxRake,
                                        round(alpha[1],5),round(alpha[2],5),
                                        0,round(re.rake,5), 
                                        ni,nt,nb,nc,d,
                                        sep = "_"),fixed=TRUE),
                             fixed = TRUE),'.RData')
    
    
    #See if this parameter set has already been run, 
    #if so, load it. If not, run the model.
    if(isError(tryCatch(readRDS(outname), 
                        error = function(e) e))){
      datalist<-list()
      for (p  in 1:d){
        #simulate the data      
        mydata <- sim_norm_data(n.site, n.quad, phi.site, psi.quad, 
                                mu,sig.site,sig.within, maxRake, 
                                alpha, beta, 
                                RE = TRUE, re.rake)
        #type is random effects model      
        mydata <- bind_rows(mutate(mydata, type="RE"))
        #site occupancy is unknown unless a positive rake score is recorded
        mydata$site.occ <- NA
        mydata <- mydata %>% group_by(site, type) %>% 
          mutate(site.occ = ifelse(sum(z) > 0, 1, NA)) %>%   ungroup()
        datalist[[p]]<- mydata
      }
      #save simulated data    
      saveRDS(datalist, file = outnamedata)
      tic()
      #initialize parallel workers    
      registerDoMC(detectCores()-1)
      getDoParWorkers()
      #fit the model for each dataset in parallel    
      out<-foreach(p=1:d, .inorder = FALSE,
                   .packages = c("runjags","dplyr"),
                   .verbose = TRUE,
                   .errorhandling = "pass") %dopar% 
        mcmc_norm_run(datalist[[p]], maxRake = maxRake, 
                      ni = ni, nt = nt, 
                      nb = nb, nc = nc)
      toc()
      #save model ouput
      saveRDS(out,file = outname)
    } else {
      print("Already ran it!")
      out<-readRDS(outname)
    }
    #assess convergence, error  
    truth<- c(re.rake, beta, phi.site, psi.quad, 
              sig.within, sig.site, mu,alpha)
    mpe<-0
    conv<-rep(0,length(truth))
    for(p in 1:100){
      mod<-out[[p]]$model
      convcurr<-c(mod$BUGSoutput$summary["sd_p","Rhat"]<=1.1,
                  mod$BUGSoutput$summary["beta_p","Rhat"]<=1.1,
                  mod$BUGSoutput$summary["phi","Rhat"]<=1.1,
                  mod$BUGSoutput$summary["psi","Rhat"]<=1.1,
                  mod$BUGSoutput$summary["sig_within","Rhat"]<=1.1,
                  mod$BUGSoutput$summary["sig_site","Rhat"]<=1.1,
                  mod$BUGSoutput$summary["bio_intercept","Rhat"]<=1.1,
                  mod$BUGSoutput$summary["alpha_p[1]","Rhat"]<=1.1,
                  mod$BUGSoutput$summary["alpha_p[2]","Rhat"]<1.1)
      conv<-conv+convcurr
      #if convergent, calculate error    
      if(sum(convcurr)==length(truth)){
        #parameter estimates are mean of posterior distributions      
        est<-c(mod$BUGSoutput$mean$sd_p,mod$BUGSoutput$mean$beta_p,
               mod$BUGSoutput$mean$phi,mod$BUGSoutput$mean$psi,
               mod$BUGSoutput$mean$sig_within,mod$BUGSoutput$mean$sig_site,
               mod$BUGSoutput$mean$bio_intercept,mod$BUGSoutput$mean$alpha_p)
        mpe<-mpe+(est-truth)/abs(truth) 
      }else{
        #not convergent
        #d holds the number of convergent datasets
        d<- d-1
      }
    }
    
    #update previous cost, momentum, perturbation
    Gem1<-Ge
    Mem1<-Me
    com1<-co
    
    #current cost is a combination of error and convergence
    co<-abs(mpe)/d-(conv)/100
    #if cost decreased, momentum term is positive  
    #if cost increased, turn around  
    m<-(2*(co<com1)-1)*0.8
    
    #momentum term in the +/- direction of the last movement  
    Me<-m*(Gem1+Mem1)
    #random perturbation  
    Ge<-c(rnorm(2,0,0.05),rnorm(1,0,0.5),rnorm(2,0,0.1),
          rnorm(2,0,0.5),rnorm(1,0,0.1),rnorm(1,0,0.5))
    #update parameters
    phi.site<-phi.site+Ge[1]+Me[1]
    psi.quad<-psi.quad+Ge[2]+Me[2]
    mu<-mu+Ge[3]+Me[3]
    sig.site<-sig.site+Ge[4]+Me[4]
    sig.within<-sig.within+Ge[5]+Me[5]
    alpha0<-alpha[1]+Ge[6]+Me[6]
    alpha1<-alpha[2]+Ge[7]+Me[7] 
    beta<-beta+Ge[8]+Me[8]
    re.rake<-re.rake+Ge[9]+Me[9]
    alpha<-sort(c(alpha0,alpha1))
    #make sure still within bounds  
    phi.site<-min(phi.site,phimax)
    phi.site<-max(phi.site,phimin)
    psi.quad<-min(psi.quad,psimax)
    psi.quad<-max(psi.quad,psimin)
    sig.site<-min(sig.site,sig.sitemax)
    sig.site<-max(sig.site,sig.sitemin)
    sig.within<-min(sig.within,sig.withinmax)
    sig.within<-max(sig.within,sig.withinmin)
    mu<-min(mu,mumax)
    mu<-max(mu,mumin)
    beta<-min(beta,betamax)
    beta<-max(beta,betamin)
    re.rake<-min(re.rake,re.rakemax)
    re.rake<-max(re.rake,re.rakemin)
    alpha[1]<-min(alpha[1],alpha0max)
    alpha[1]<-max(alpha[1],alpha0min)
    alpha[2]<-min(alpha[2],alpha1max)
    alpha[2]<-max(alpha[2],alpha1min)  
    
    numran<-numran+1
    print(paste("You've ran", numran, "iterations"))
  }
}

