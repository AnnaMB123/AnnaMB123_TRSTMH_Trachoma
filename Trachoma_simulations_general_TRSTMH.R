###Running stochastic, individual-based trachoma transmission model###

source("") # File for all model functions

library(iterators)
library(foreach)
library(doParallel)

cores=detectCores()
cl<-makeCluster(cores[1]-1, type="FORK") #
#Parallelise runs#
registerDoParallel(cl)

#Parameters relating to simulation#

sim_params<-c(
  timesim=56*52, #years total duration of simulation (*52 so in weeks) including burn-in
  burnin=40*52, #years burin to be discarded (*52 so in weeks)
  Freq_MDA=1, #1 for annual, 0.5 for biannual, 2 for q 2 years etc
  N_MDA=16, # No. of rounds of MDA to be carried out in simulation
  n_sim=4000 #number of simulations
  
)


#General params relating to transmission of infection
parameters <- c(N=1000, #Population size
                av_I_duration=2,av_ID_duration=200/7,inf_red=0.45,min_ID=11, #Parameters relating to duration of infection period, including ID period
                av_D_duration=300/7,dis_red=0.3,min_D=1, #Parameters relating to duration of disease period
                v_1=1,v_2=2.6,phi=1.4,epsilon=0.5,#Parameters relating to lambda function- calculating force of infection
                #Parameters relating to MDA
                MDA_Cov=0.8,#MDA coverage
                MDA_Eff= 0.85, # Efficacy of treatment
                rho=0 #correlation parameter for systematic non-compliance function
)

#Demography parameters
demog<-c(tau=1/(40*52),max_age=60*52,mean_age=20*52) #tau= death rate in weeks^-1, max_age= maximum age in population, mean_age=mean age in population


#Function to run a single simulation with MDA at time points determined by vector MDA_times
#Output is true prevalence of infection/disease at each timestep in whole population (True_prev_Infection/True_Prev_Disease) or children aged 1-9 (True_Prev_disease_children_1_9/True_Prev_infection_children_1_9)
sim_Ind_MDA<-function(params,vals,timesim,demog,bet){
  with(as.list(c(params,vals,timesim,demog,bet)),{
    
    N_Children_ages_1_9<-array()
    
    N_True_Infected_children_1_9<-array()
    True_Prev_Infection_children_1_9<-array()
    
    N_True_Diseased_children_1_9<-array()
    True_Prev_Disease_children_1_9<-array()
    
    True_Prev_Disease<-array()
    True_prev_infection<-array()
    
    Ages<-array()
    
    N_Inf<-array()
    
    Inf_status<-array()
    
    t<-array()
    
    MDA_round<-array()
    
    MDA_round_kids<-array()
    
    for (i in timesim){
      
      if (i %in% MDA_times){
        
        MDA_round[i]<-which(MDA_times==i)
        
        vals<- MDA_timestep(vals,params,demog,MDA_round=MDA_round[i])
        
        
        t[i]<-i
        
        N_Children_ages_1_9[i]<-length(which(vals[[11]]<10*52 & vals[[11]]>52))
        
        N_True_Infected_children_1_9[i]<-sum(vals[[1]][which(vals[[11]]<10*52 & vals[[11]]>52)])
        True_Prev_Infection_children_1_9[i]<-N_True_Infected_children_1_9[i]/N_Children_ages_1_9[i]
        
        N_True_Diseased_children_1_9[i]<-sum(vals[[2]][which(vals[[11]]<10*52 & vals[[11]]>52)])
        True_Prev_Disease_children_1_9[i]<-N_True_Diseased_children_1_9[i]/N_Children_ages_1_9[i]
        
        True_Prev_Disease[i]<-sum(vals[[2]])/N
        True_prev_infection[i]<-sum(vals[[1]])/N
        
        
        
        
      }
      
      if (i %in% MDA_kids){
        
        MDA_round_kids[i]<-which(MDA_kids==i)
        
        vals<- MDA_timestep_kids(vals,params,demog,MDA_round_kids[i])
        
        t[i]<-i
        
        N_Children_ages_1_9[i]<-length(which(vals[[11]]<10*52 & vals[[11]]>52))
        
        N_True_Infected_children_1_9[i]<-sum(vals[[1]][which(vals[[11]]<10*52 & vals[[11]]>52)])
        True_Prev_Infection_children_1_9[i]<-N_True_Infected_children_1_9[i]/N_Children_ages_1_9[i]
        
        N_True_Diseased_children_1_9[i]<-sum(vals[[2]][which(vals[[11]]<10*52 & vals[[11]]>52)])
        True_Prev_Disease_children_1_9[i]<-N_True_Diseased_children_1_9[i]/N_Children_ages_1_9[i]
        
        True_Prev_Disease[i]<-sum(vals[[2]])/N
        True_prev_infection[i]<-sum(vals[[1]])/N
        
        
        
      }
      
      {
        
        vals<- stepF_fixed(vals,params,demog,bet)
        
        
        
        t[i]<-i-1
        
        N_Children_ages_1_9[i]<-length(which(vals[[11]]<10*52 & vals[[11]]>52))
        
        N_True_Infected_children_1_9[i]<-sum(vals[[1]][which(vals[[11]]<10*52 & vals[[11]]>52)])
        True_Prev_Infection_children_1_9[i]<-N_True_Infected_children_1_9[i]/N_Children_ages_1_9[i]
        
        N_True_Diseased_children_1_9[i]<-sum(vals[[2]][which(vals[[11]]<10*52 & vals[[11]]>52)])
        True_Prev_Disease_children_1_9[i]<-N_True_Diseased_children_1_9[i]/N_Children_ages_1_9[i]
        
        
        
        True_Prev_Disease[i]<-sum(vals[[2]])/N
        True_prev_infection[i]<-sum(vals[[1]])/N
        
        
        
      }
      
    }
    
    
    output<-(list(t,True_Prev_Infection_children_1_9,True_Prev_Disease_children_1_9,True_Prev_Disease,True_prev_infection))
    names(output)<-c("Time","True_Prev_Infection_children_1_9","True_Prev_Disease_children_1_9","True_Prev_Disease","True_prev_infection")
    return(output)
    
  })
}


####SETTING 1 ####

#Longitudinal simulations:
#Set initial conditions and seed infection
set.seed(5)
Set_inits(params = parameters, demog=demog)
set.seed(5)
Seed_infection(params = parameters, vals=vals)

timesim<-sim_params["timesim"]
burnin<-sim_params["burnin"]
nsim<-sim_params["n_sim"]

#Specify which individuals get treated at each MDA, create treatment matrix 
set.seed(5)
Tx_mat<-Tx_matrix(params=parameters,sim_param=sim_params)

#Setting transmission parameter beta
set.seed(5)
bet<-(runif(sim_params["n_sim"],0.165,0.175)) #Corresponds roughly to 40% TF1-9 prevalence baseline (after burnin), filter at baseline to set bounds
hist(bet)

####1A: MDA carried out annually as planned 

#Times for MDA to be carried out
MDA_times<-Set_t_MDA(sim_params)
MDA_kids<-99999999 #Not simulating treatment of kids here

data_store_all_sim1A<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()

data_store_all_sim1A<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  
  
  data_store_all_sim1A[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:timesim,demog=demog,bet=bet[i])
}
end_time<-Sys.time()
end_time-start_time


#### Simulation 1B: 2 rounds, miss 2nd round
#Times for MDA to be carried out
MDA_times<-Set_t_MDA(sim_params)
MDA_times<-MDA_times[!MDA_times %in% MDA_times[3]] 
MDA_kids<-99999999 #Not simulating treatment of kids only here

data_store_all_sim1B<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()

data_store_all_sim1B<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  
  
  data_store_all_sim1B[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:timesim,demog=demog,bet=bet[i])
}
end_time<-Sys.time()
end_time-start_time


#### Simulation 1B 2 year delay: 2 rounds, miss 2nd & 3rd round
#Times for MDA to be carried out
MDA_times<-Set_t_MDA(sim_params)
MDA_times<-MDA_times[!MDA_times %in% MDA_times[3] & !MDA_times %in% MDA_times[4]] 


data_store_all_sim1B_2y<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()

data_store_all_sim1B_2y<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  data_store_all_sim1B_2y[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:timesim,demog=demog,bet=bet[i])
}
end_time<-Sys.time()
end_time-start_time


####Simulation 1C: 2 rounds, miss 2nd round, "make up" round in the year activities are resumed (Mitigation M1)

#Times for MDA to be carried out
MDA_times<-Set_t_MDA(sim_params)
MDA_times<-MDA_times[!MDA_times %in% MDA_times[3]] 
Make_up_MDA<-MDA_times[3]+26
MDA_times<-c(MDA_times, Make_up_MDA)
MDA_times<-sort(MDA_times)
MDA_kids<-99999999 #Not simulating treatment of kids only here


data_store_all_sim1C<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()

data_store_all_sim1C<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  
  data_store_all_sim1C[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:timesim,demog=demog,bet=bet[i])
}
end_time<-Sys.time()
end_time-start_time


####Simulation 1C 2y: 2 rounds, miss 3rd and 4th, "make up" round in the year activities are resumed (Mitigation M1)
#Times for MDA to be carried out
MDA_times<-Set_t_MDA(sim_params)
MDA_times<-MDA_times[!MDA_times %in% MDA_times[3] & !MDA_times %in% MDA_times[4]] 
Make_up_MDA<-MDA_times[3]+26
MDA_times<-c(MDA_times, Make_up_MDA)
MDA_times<-sort(MDA_times)

data_store_all_sim1C_2y<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()

data_store_all_sim1C_2y<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  
  data_store_all_sim1C_2y[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:timesim,demog=demog,bet=bet[i])
}
end_time<-Sys.time()
end_time-start_time


####Simulation 1D 1y: 2 rounds, miss 3rd round, "make up" round in the year activities are resumed treats only kids (Mitigation M2)
#Times for MDA to be carried out
MDA_times<-Set_t_MDA(sim_params)
MDA_times<-MDA_times[!MDA_times %in% MDA_times[3]] 
MDA_kids<-MDA_times[3]+26

#Create an empty matrix for kids only MDA
Tx_mat_kids<-matrix(0,parameters["N"],sim_params["N_MDA"])
##Designate 4th round from original matrix as first kids only round, doesn't matter that the original matrix includes all ages, the do_MDA_kids function will only treat kids
Tx_mat_kids[,1]<-Tx_mat[,4]        

#And then remove the 4th round from original matrix to ensure exact same children aren't treated twice
Tx_mat<-Tx_mat[,-c(4)]

data_store_all_sim1D<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()

data_store_all_sim1D<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  
  data_store_all_sim1D[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:timesim,demog=demog,bet=bet[i])
}
end_time<-Sys.time()
end_time-start_time


####Simulation 1D 2years: 2 rounds, miss 3rd & 4th round, "make up" round in the year activities are resumed treats only kids (Mitigation M2)
#Times for MDA to be carried out
MDA_times<-Set_t_MDA(sim_params)
MDA_times<-MDA_times[!MDA_times %in% MDA_times[3] & !MDA_times %in% MDA_times[4]] 
MDA_kids<-MDA_times[3]+26

#Recreate original treatment matrix
set.seed(5)
Tx_mat<-Tx_matrix(params=parameters,sim_param=sim_params)

#Create an empty matrix for kids only MDA
Tx_mat_kids<-matrix(0,parameters["N"],sim_params["N_MDA"])

##Designate 4th round from original matrix as first kids only round, doesn't matter that the original matrix includes all ages, the do_MDA_kids function will only treat kids
Tx_mat_kids[,1]<-Tx_mat[,4]        

#And then remove the 4th round from original matrix to ensure exact same children aren't treated twice
Tx_mat<-Tx_mat[,-c(4)]

data_store_all_sim1D_2y<-vector(length(nsim),mode='list') #store output from multiple runs including burnin

start_time<-Sys.time()

data_store_all_sim1D_2y<-foreach(i =1:nsim) %dopar% {
  set.seed(i+5)
  
  data_store_all_sim1D_2y[[i]]<-sim_Ind_MDA(params=parameters, vals=vals,timesim=1:timesim,demog=demog,bet=bet[i])
}
end_time<-Sys.time()
end_time-start_time

#Removing burnin

#Create empty matrices length burnin:end of simulation
True_Prev_Infection_children_1_9<-matrix(0,nrow=(timesim-burnin),ncol=nsim)
True_Prev_Disease_children_1_9<-matrix(0,nrow=(timesim-burnin),ncol=nsim)
True_Prev_Disease<-matrix(0,nrow=(timesim-burnin),ncol=nsim)
True_Prev_Infection<-matrix(0,nrow=(timesim-burnin),ncol=nsim)
Time<-seq(0,(timesim-burnin-1),1) #Start time 0 from end of burnin

data_store_all_sim<-data_store_all_sim_x

#Fill with output, 
for(i in 1:nsim){
  True_Prev_Disease[,i]<-data_store_all_sim[[i]]$True_Prev_Disease[(burnin+1):timesim]
  True_Prev_Infection[,i]<-data_store_all_sim[[i]]$True_prev_infection[(burnin+1):timesim]
  
  True_Prev_Infection_children_1_9[,i]<-data_store_all_sim[[i]]$True_Prev_Infection_children_1_9[(burnin+1):timesim]
  True_Prev_Disease_children_1_9[,i]<-data_store_all_sim[[i]]$True_Prev_Disease_children_1_9[(burnin+1):timesim]  
  
}



########### SETTING 2 ###################

# Repeat above with alternative transmission parameter beta, and interruption is in year 2 rather than year 3
set.seed(5)
bet<-(runif(sim_params["n_sim"],0.05,0.12)) #Corresponds to roughly 20% true TF in ages 1-9, filter output depending on category boundaries needed

