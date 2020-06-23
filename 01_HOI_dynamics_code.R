
# R script for HOI models by **Gaurav Baruah** and **Pragya Singh** for the paper 
#Title: "Higher-order interactions and species coexistence." 2020 for Theoretical Ecology.

rm(list=ls())
library(pracma)
library(entropy)
set.seed(1234)
HOI_function_dynamics<-function(time,species,InitialAB=1, 
                                interaction= 'pairwise',
                                alphas, HOI_matrix)
  
{
  
  
  time<-time
  species<-species
  dt<-0.05
  s<- array(dim=c(1,species))
  two_alpha<-array(dim=c(time,species ,species))
  alpha<-array(dim=c(species,species))
  A<-array(dim=c(time,species))
  Horder<-array(dim=c(time,species))
  g<-array(dim=c(time,species))
  H_alpha<-array(dim=c(species,species,species))
  N<-array(dim=c(time,species))

  #inital values
  N[1,]<-InitialAB #initial abundance
  alpha <- alphas
  H_alpha<- HOI_matrix 
    
for (t in 1:(time-1)) {
     
        if (interaction == 'pairwise'){
          
          
       A[t,]<- alpha%*%N[t,] #cumulative pair-wise affect
       N[t+1,]<-N[t,] + (1 - A[t,])*N[t,]*dt 
      
           
       N[t+1,which(N[t+1,]<1e-8)] <- 0
       N[t+1,which(N[t+1,]>1e8)]  <- 0
       
        }
        else if (interaction == 'HOI'){
          
        for(i in 1:species){
          Horder[t,i]<-t(N[t,])%*% H_alpha[i,,]%*%N[t,] #summed HOI effects
          }
          
          A[t,]<- alpha%*%N[t,] #cumulative pair-wise affect
          N[t+1,]<-N[t,] + (1 + Horder[t,] - A[t,])*N[t,i]*dt 
       
          
         #feasible densities: lower max is 1e-8 below which species are assigned 0.
         #and lower max which is 1e12 above which species are assigned 0
         N[t+1,which(N[t+1,]<1e-8)] <- 0
         N[t+1,which(N[t+1,]>1e12)]  <- 0
          
        }

    }

      return(list(alpha= alpha, HOI_alpha = H_alpha ,Density=N))
}
  



#generating the dynamical data that fullfilles Weyl's inequality criteria
#parms : list of parameters
dynamics.data<-function(parms){
mydata <- list()
dat <-HOI_function_dynamics(time = parms$time,species = parms$species,InitialAB = 0.1,
                            interaction = parms$interaction,alphas=parms$alpha, 
                            HOI_matrix=parms$HOI_alpha)
inter.alpha<-dat$alpha
ss<-parms$species
diag(inter.alpha)<-0 # this gives the interspecific competition matrix only with diagonal terms all zero 
eigen_pairwise1<-sort(Re(eigen(-inter.alpha, only.values = T)$values))
intra.compt<-dat$alpha[1,1]

stability_criteria <- -intra.compt + (eigen_pairwise1[ss]) < 0
temp<-list(dat$alpha, dat$HOI_alpha, dat$Density,stability_criteria)
names(temp)<-c("Pairwise.matrix", "HOI.matrix", "Density", "Coexistence.criteria")
mydata<-temp   
  
return(mydata)
}






