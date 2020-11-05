##############UNITS AS WEEKS###############

#Step function ie. transitions in each time non-MDA timestep # 

stepF_fixed <- function(vals,params,demog,bet){
  
  with(as.list(c(vals,params,demog,bet)),{
    
    #Assign names to variables just so it's easier to track #
    
    IndI<-vals[[1]]
    IndD<-vals[[2]]
    No_Inf<-vals[[3]]
    T_latent<-vals[[4]]
    T_ID<-vals[[5]]
    T_D<-vals[[6]]
    Ind_latent<-vals[[7]]
    Ind_ID_period_base<-vals[[8]]
    Ind_D_period_base<-vals[[9]]
    bact_load<-vals[[10]]
    Age<-vals[[11]]
    
    
    #Selecting individuals to make transistions#
    
    # New I=  Becoming infected, not yet diseased
    Ss <- which(IndI==0) # # susceptible indivs available for infection. Disease positive individuals can also get infected as scaled rate
    #Susceptibles acquiring new infections# This give a lambda for each individual dependent on age and disease status
    lambda_step<-1-exp(-getlambdaStep(params=parameters, Age=Age, demog=demog, bact_load=bact_load, IndD=IndD, bet=bet))
    
    newInf<- Ss[which(runif(Ss) < lambda_step[Ss])] # New infections 
    
    # New ID= Infected people who are infected but not yet clinically diseased, progressing to diseased
    Is <- which(IndI==1 & IndD==0) #Infected but not diseased
    T_latent<-T_latent-1 # with each timestep, subtract 1, those not in latent period go to -1 then should reset to zero, those in latent period should count down with each timestep
    newDis<-which(T_latent==0) #Designated latent period for that individual has expired
    T_latent[T_latent<1]<-0 #Reset all not in latent period to 0 at each timestep
    
    # New D= Clearing infection, becoming disease only, 
    IDs<-which(IndI==1 & IndD==1)
    bact_load[IDs]<-bacterialLoad(No_Inf = No_Inf[IDs],parameters=parameters)
    T_ID<-T_ID-1 # with each timestep, subtract 1, those uninfected should go to -1 then be reset to zero, those infected should count down with each timestep
    newClearInf<-which(T_ID==0) #Designated infectious period for that individual has expired
    T_ID[T_ID<1]<-0 #Reset all uninfected to 0 at each timestep
    
    # New S= Clear disease
    Ds<-which(IndI==0 & IndD==1)
    bact_load[Ds]<-0
    T_D<-T_D-1 #Each timestep, substract 1, those not diseased should go to -1 then be reset to zero, those diseased should count down with each timestep
    newClearDis<-which(T_D==0) #Designated diseased period for that individual has expired
    T_D[T_D<1]<-0 #Reset all not in diseased period to 0
    
    #Tracking infection history
    No_Inf[newInf]<-No_Inf[newInf]+1
    
    ##TRANSITION## become infected##
    IndI[newInf] <- 1 #if they've become infected, become I=1
    ##TRANSITION## become diseased (and infected)
    IndD[newDis] <- 1 # if they've become diseased they become D=1
    ##TRANSITION## Clear infection (remain diseased)
    IndI[newClearInf]<-0 #clear infection they become I=0
    ##TRANSITION## Clear disease (become susceptible again)
    IndD[newClearDis]<-0 #clear disease they become D=0
    
    #Setting duration of infection/disease for individuals#
    T_latent[newInf]<-Ind_latent[newInf] # When individual becomes infected, set their latent period- this is how long they remain in categore I infected but not diseased
    # When individual becomes diseased, set their infected+diseased period- this is how long ID for
    T_ID[newDis]<-ID_period_function(Ind_ID_period_base=Ind_ID_period_base[newDis],No_Inf = No_Inf[newDis],parameters=parameters)
    #When indiv clears infection, their diseased only period is set, this how long D for 
    T_D[newClearInf]<-D_period_function(Ind_D_period_base=Ind_D_period_base[newClearInf],No_Inf = No_Inf[newClearInf],parameters=parameters)
    
    #Update age, all age by 1w at each timestep, and reseting all "reset indivs" age to zero
    #Reset_indivs - Identify individuals who die in this timeset, either reach max age or random death rate
    Age<-Age+1
    reset_indivs<-Reset(Age=Age,demog=demog,params=parameters)
    Age[reset_indivs]<-0
    
    #Resetting new parameters for all new individuals created
    
    IndI[reset_indivs]<-0
    IndD[reset_indivs]<-0
    No_Inf[reset_indivs]<-0 
    T_latent[reset_indivs]<-0 
    T_ID[reset_indivs]<-0
    T_D[reset_indivs]<-0
    
    
    vals <- list(IndI,IndD,No_Inf,T_latent,T_ID,T_D,Ind_latent,Ind_ID_period_base,Ind_D_period_base,bact_load,Age) 
    
    
    return(c(vals))
    
  }
  
  )}

#lambda scaled according to disease status
getlambdaStep = function(params,Age,demog,bact_load,IndD,bet){
  with(as.list(c(params,Age,demog,bact_load,IndD,bet)),{
    totalLoad <- c(0,0,0)
    y_children <- which(Age>=0 & Age<9*52) #Young children
    totalLoad[1] <- sum(bact_load[y_children])/length(y_children)
    o_children <- which(Age>=9*52 & Age<15*52) #Older children
    totalLoad[2] <- sum(bact_load[o_children])/length(o_children)
    adults <- which(Age>=15*52) #Adults
    totalLoad[3] <- sum(bact_load[adults])/length(adults)
    
    prevLambda <- bet * (v_1 * totalLoad + v_2 * totalLoad ^ (phi + 1))
    
    demog_matrix <- t(matrix(rep(c(length(y_children),length(o_children),length(adults))/N,3), nrow=3, ncol=3))
    social_mixing <- (epsilon*diag(3) + (1 - epsilon))*demog_matrix #scales mixing with other groups by 50% less than within same group? not 100% sure here
    
    lambda <- as.vector(social_mixing%*%prevLambda) 
    ## if disease == 1, lambda*.5; if disease == 0, lambda * 1
    return((lambda[findInterval(Age,c(0,9*52,15*52,max_age*52))])*(0.5+0.5*(1-IndD))) 
  })
}

#Function to identify individuals who either die due to background mortality, or who reach max age. 
Reset<-function(Age,demog,params){
  with(as.list(c(Age,demog,params)),{
    reset_indivs<-(which(runif(N)<1-exp(-tau)| Age>max_age))  
    return(reset_indivs)
  })
}


#Function to decide when MDA occurs, return vector of timepoints when MDA carried out. 
#Can be annually or less freq/more freq

Set_t_MDA=function(sim_params){
  with(as.list(sim_params),{
    MDA_t<-c()
    for (i in 1: N_MDA){
      MDA_t[i]<-burnin+((i-1)*52*Freq_MDA) + 52 #Add 52 if want first MDA to be 1 year after burnin
    }
    return(MDA_t)
  })
}

#Decide who is cured during community-based MDA
doMDA<-function(params,Age,MDA_round){
  with(as.list(c(params,Age,MDA_round)),{
    babies<-which(Age<26)
    MDA_round<-MDA_round
    treated_babies<-babies[which(Tx_mat[babies,MDA_round]==1)]
    cured_babies<-treated_babies[which(runif(treated_babies)<(MDA_Eff*0.5))]
    older<-which(Age>26)
    treated_older<-older[which(Tx_mat[older,MDA_round]==1)]
    cured_older<-treated_older[which(runif(treated_older)<MDA_Eff)]
    treated_cured<-c(cured_babies,cured_older)
    return(treated_cured)
  })
}


#Decide which children are cured during children-only MDA#
doMDAkids<-function(params,Age,MDA_round_kids){
  with(as.list(c(params,Age,MDA_round_kids)),{
    kids<-which(Age<10*52 & Age>104)
    treated_kids<-kids[which(Tx_mat_kids[kids,MDA_round_kids]==1)]
    treated_cured<-treated_kids[which(runif(treated_kids) < MDA_Eff)] 
    return(treated_cured)
  })
}


#Create matrix to determine who gets treated at each MDA round, allowing for systematic non-compliance as specified by Dyson
Tx_matrix <- function(params,sim_params){
  with(as.list(c(params,sim_params)),{
    ind_treat <- matrix(0,N,N_MDA)
    for(i in 1: N){
      z <- runif(N)
      ind_treat[,1]<-z < MDA_Cov #Randomly assign first treatment
      for (k in 2:N_MDA){
        ind_treat[i,k]<-rbinom(1,1,((MDA_Cov*(1-rho)+(rho*sum(ind_treat[i,1:k])))/(1+(k-2)*rho))) #Subsequent treatment probs function of previous treatments
        
      }
    }
    return(ind_treat)
  })
}


#This is for timesteps in which community-level MDA occurs#

MDA_timestep <- function(vals,params,demog,MDA_round){
  
  with(as.list(c(vals,params,demog,MDA_times,MDA_round)),{
    
    #Assign names to variables just so it's easier to track #
    
    IndI<-vals[[1]]
    IndD<-vals[[2]]
    No_Inf<-vals[[3]]
    T_latent<-vals[[4]]
    T_ID<-vals[[5]]
    T_D<-vals[[6]]
    Ind_latent<-vals[[7]]
    Ind_ID_period_base<-vals[[8]]
    Ind_D_period_base<-vals[[9]]
    bact_load<-vals[[10]]
    Age<-vals[[11]]
    
    #Id who is treated and cured
    
    treated_cured<-doMDA(params=parameters,Age=Age,MDA_round=MDA_round)
    
    #Set treated/cured indivs infection status and bacterial load to 0
    IndI[treated_cured]<-0 #clear infection they become I=0
    bact_load[treated_cured]<-0 #clear disease they become D=0
    
    
    vals <- list(IndI,IndD,No_Inf,T_latent,T_ID,T_D,Ind_latent,Ind_ID_period_base,Ind_D_period_base,bact_load,Age) 
    
    
    return(c(vals))
    
  }
  
  )} 

#This is for timesteps in which children-only MDA occurs#

MDA_timestep_kids <- function(vals,params,demog,MDA_round_kids){
  
  with(as.list(c(vals,params,demog,MDA_round_kids)),{
    
    #Assign names to variables just so it's easier to track #
    
    IndI<-vals[[1]]
    IndD<-vals[[2]]
    No_Inf<-vals[[3]]
    T_latent<-vals[[4]]
    T_ID<-vals[[5]]
    T_D<-vals[[6]]
    Ind_latent<-vals[[7]]
    Ind_ID_period_base<-vals[[8]]
    Ind_D_period_base<-vals[[9]]
    bact_load<-vals[[10]]
    Age<-vals[[11]]
    
    #Id who is treated and cured
    treated_cured<-doMDAkids(params=parameters, Age=Age,MDA_round_kids = MDA_round_kids )
    
    #Set treated/cured indivs infection status and bacterial load to 0
    IndI[treated_cured]<-0 #clear infection they become I=0
    bact_load[treated_cured]<-0 #clear disease they become D=0
    
    
    vals <- list(IndI,IndD,No_Inf,T_latent,T_ID,T_D,Ind_latent,Ind_ID_period_base,Ind_D_period_base,bact_load,Age) 
    
    
    return(c(vals))
    
  }
  
  
  
  )}




#Function to give duration of active infection 
ID_period_function<-function(Ind_ID_period_base,No_Inf,parameters){
  with(as.list(c(Ind_ID_period_base,No_Inf,parameters)),{
    T_ID<- round((Ind_ID_period_base-min_ID)*exp(-inf_red*(No_Inf-1))+min_ID)
    return(T_ID)
  })
}


#Function to give duration of disease only period
D_period_function<-function(Ind_D_period_base,No_Inf,parameters){
  with(as.list(c(Ind_D_period_base,No_Inf,parameters)),{
    T_D<- round((Ind_D_period_base-min_D)*exp(-dis_red*(No_Inf-1))+min_D)
    return(T_D)
  })
}

#Function to scale bacterial load according to infection history
bacterialLoad = function(No_Inf,parameters){
  with(as.list(c(No_Inf,parameters)),{
    b1 = 1     ## param for bacterial load function from Pinsent 2018
    ep2 = 0.114 ## param for bacterial load function fron Pinsent 2018
    bact_load=(b1*exp((No_Inf-1)*-ep2))
    return(bact_load)}) ## functional form for reduction in bacterial load with repeated infections
}



# Set Initial values: 

with(as.list(c(params,demog)),{
  IndI <- rep(0,N) #Individual's infected status
  IndD <- rep(0,N) #Individual's disease status
  No_Inf<-rep(0,N) #Individual's total number of infections, increases by 1 each time they become newly infected
  T_latent<-rep(0,N) #Duration of latent period (I), i.e. infected but not yet clinically diseased.
  T_ID<-rep(0,N) #Duration of current ID period,set when becomes infected and counts down with each time step
  T_D<-rep(0,N) #Duration individual spends diseased after clearing infection
  Ind_latent<-rep(av_I_duration,N) #Individual's latent period fixed for now
  Ind_ID_period_base<-rpois(N,av_ID_duration) #Individual's baseline ID period (first infection)
  Ind_D_period_base<-rpois(N,av_D_duration) #Individual's baseline diseased period (first infection)
  bact_load<-rep(0,N) #Baseline bacterial load set to zero
  Age<-init_ages(params=parameters,demog=demog)
  
  vals <<- list(IndI,IndD,No_Inf,T_latent,T_ID,T_D,Ind_latent,Ind_ID_period_base,Ind_D_period_base,bact_load,Age) 
  names(vals)<<-c("IndI","IndD","No_Inf","T_latent","T_ID","T_D","Ind_latent","Ind_ID_period_base","Ind_D_period_base","bact_load","Age") 
  
})

}

#Seed infection
Seed_infection<-function(params,vals){
  vals[[1]][1:round(parameters['N']*0.01)] <<- 1 
  Init_infected <- which(vals[[1]]==1)
  vals[[4]][Init_infected]<<-vals[[7]][Init_infected]# set latent period for those infected at start of simulation
  vals[[3]][Init_infected]<<-1 # Set number of infections to 1 for those infected at start of simulation.
}


#Initialise age distribution 
#NOTE- ages are in weeks
init_ages<-function(params,demog){
  with(as.list(c(params,demog)),{
    ages<-1:max_age
    nAges <-length(ages)  
    propAges <- rep(NA, nAges) ## ensure the population is in equilbirum
    propAges[1:(nAges-1)]  <-  exp( - ages[1:(nAges-1)]/mean_age ) - exp( - ages[2:nAges]/mean_age ) 
    propAges[nAges] <- 1 - sum( propAges[1:(nAges-1)] )
    
    return(sample(ages,N,replace = T,prob=propAges))
    
  })
}