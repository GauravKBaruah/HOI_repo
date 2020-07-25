rm(list=ls())
# R script for HOI analysis by **Gaurav Baruah** and **Pragya Singh** for the paper 
#Title: "Higher-order interactions and species coexistence." 2020 for Theoretical Ecology.
# R scripts for evaluating a range of parameter choices of HOI when Weyl's inequality is satisified vs. when it 
# is not satisfied.



#sourcing the source function for the dynamical HOI model
source('~./multispecies_negative_hoi.R', echo=F)
library(gdata)
library(grid)
library(gridExtra)
library(dplyr)




### SYMMETRIC ALPHA MATRICES : 
# Weyl's inequality satisfied: and range of parameter space explored of negative HOIs ----
d<-0.4
parms<-list()
time =3000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model
mean.d<-seq(-0.15, 0.15, 0.05)


#create an empty dataframe
fact.w<- expand.grid(d=0.4, mean.HOI.diff= mean.d,
                     `random_seed`=4327+(1:30)*100) %>%
  as_tibble %>%
  mutate(`Species_richness`=0,
         `Weyls` = "" )

# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for(i in 1: nrow(fact.w)){
  
  set.seed(fact.w$random_seed[i])
  alpha=matrix(NA, nrow = species, ncol = species) 
  #interspecific effects
  temp<-runif(length(upperTriangle(alpha,diag = F)),0.01,0.04)
  upperTriangle(alpha,diag = F,byrow = T)<-temp
  alpha<-t(alpha)
  upperTriangle(alpha,diag = F,byrow = T)<-temp
  
  
  diag(alpha)<- fact.w$d[i] # intraspecfic effects
  for(k in 1:species){
    HOI_a[k,,]<- -matrix(runif(species^2, 0.01,0.09), nrow=species,ncol=species) 
    HOI_a[k,k,]<- -runif(species,0.01+fact.w$mean.HOI.diff[i], 0.09+fact.w$mean.HOI.diff[i]) #-U[0.01 +d, 0.09+d]
  }
  
  parameters<-list(time,species,interaction,alpha,HOI_a)
  names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")
  data.fullfilling.coexistence.criteria<-dynamics.data(parms = parameters)
  
  
  fact.w$Species_richness[i] = length(which(data.fullfilling.coexistence.criteria$Density[3000,]>0))
  fact.w$Weyls[i] = data.fullfilling.coexistence.criteria$Coexistence.criteria
  print(i)
 
  
}

# plotting
g0<-ggplot(fact.w, 
           aes(x=factor(round(mean.HOI.diff,3)),y=Species_richness))+
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_classic()+xlab("d")+ylab("Species richness")+ggtitle("A) Weyl's inequality satisfied")

g0


# Weyl's inequality NOT satisfied and range of parameter space explored of negative HOIs ----

d<-0.01
parms<-list()
time =3000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model
mean.d<-seq(-0.15, 0.15, 0.05)

fact.nw<- expand.grid(d=d, mean.HOI.diff= mean.d,
                      `random_seed`=4327+(1:30)*100) %>%
  as_tibble %>%
  mutate(`Species_richness`=0)

# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))

for(i in 1: nrow(fact.nw)){
  set.seed(fact.nw$random_seed[i])
  
  #interspecific effects
  alpha=matrix(NA, nrow = species, ncol = species)
  temp<-runif(length(upperTriangle(alpha,diag = F)),0.04,0.09)
  upperTriangle(alpha,diag = F,byrow = T)<-temp
  alpha<-t(alpha)
  upperTriangle(alpha,diag = F,byrow = T)<-temp
  #isSymmetric(alpha)
  
  diag(alpha)<- fact.nw$d[i] # intraspecfic effec
  for(k in 1:species){
    HOI_a[k,,]<- -matrix(runif(species^2, 0.01,0.09), nrow=species,ncol=species) 
    HOI_a[k,k,]<- -runif(species,0.01+fact.nw$mean.HOI.diff[i], 0.09+fact.nw$mean.HOI.diff[i])
  }
  
  parameters<-list(time,species,interaction,alpha,HOI_a)
  names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")
  data.fullfilling.coexistence.criteria<-dynamics.data(parms = parameters)
  
  
  fact.nw$Species_richness[i] = length(which(data.fullfilling.coexistence.criteria$Density[3000,]>0))
  fact.nw$Weyls[i] = data.fullfilling.coexistence.criteria$Coexistence.criteria
  
  print(i)
}


# plotting
g1<-ggplot(fact.nw, 
           aes(x=factor(round(mean.HOI.diff,3)),y=Species_richness))+
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_classic()+xlab("d")+ylab("Species richness")+
  ggtitle("B) Weyl's inequality NOT satisfied")

g1

# plotting figure 5 in main-text
plot_grid(g0,g1)


## ASYMMETRIC pairwise ALPHA Matrices:


# Range of parameter space explored of negative HOIs but
#intraspecific pairwise effects > interspecific pairwise effects  ----


d<-0.4
parms<-list()
time =3000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model
diag(alpha)<- d # intraspecfic effects
mean.d<-seq(-0.15, 0.15, 0.05)

fact.w<- expand.grid(d=0.4, mean.HOI.diff= mean.d,
                     `random_seed`=4327+(1:30)*100) %>%
  as_tibble %>%
  mutate(`Species_richness`=0)

# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for(i in 1: nrow(fact.w)){
  set.seed(fact.w$random_seed[i])
  alpha=matrix(runif(species^2, 0.01,0.04), nrow = species, ncol = species) #interspecific
  #interspecific
  diag(alpha)<- fact.w$d[i] # intraspecfic effects
  for(k in 1:species){
    HOI_a[k,,]<- -matrix(runif(species^2, 0.01,0.09), nrow=species,ncol=species) 
    HOI_a[k,k,]<- -runif(species,0+fact.w$mean.HOI.diff[i], 0.09+fact.w$mean.HOI.diff[i])
  }
  
  parameters<-list(time,species,interaction,alpha,HOI_a)
  names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")
  data.fullfilling.coexistence.criteria<-dynamics.data(parms = parameters)
  
  
  fact.w$Species_richness[i] = length(which(data.fullfilling.coexistence.criteria$Density[3000,]>0))
  
  print(i)

  
}

g0<-ggplot(fact.w, 
           aes(x=factor(round(mean.HOI.diff,3)),y=Species_richness))+
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_classic()+xlab("d")+ylab("Species richness")+ggtitle("A) Weyl's inequality satisfied")

g0


# Range of parameter space explored of negative HOIs but
#intraspecific pairwise effects < interspecific pairwise effects  ----


d<-0.01
parms<-list()
time =3000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model
alpha=matrix(runif(species^2, 0.01,0.04), nrow = species, ncol = species) #interspecific
diag(alpha)<- d[1] # intraspecfic effects
mean.d<-seq(-0.15, 0.15, 0.05)
sp.richness<-matrix(NA, nrow=length(d),ncol=length(mean.d))


fact.nw<- expand.grid(d=d, mean.HOI.diff= mean.d,
                      `random_seed`=4327+(1:30)*100) %>%
  as_tibble %>%
  mutate(`Species_richness`=0)

# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for(i in 1: nrow(fact.nw)){
 
   set.seed(fact.nw$random_seed[i])
  alpha=matrix(runif(species^2, 0.01,0.04), nrow = species, ncol = species) #interspecific
 
  diag(alpha)<- fact.nw$d[i] # intraspecfic effects
  for(k in 1:species){
    HOI_a[k,,]<- matrix(runif(species^2, 0.01,0.09), nrow=species,ncol=species) 
    HOI_a[k,k,]<- runif(species,0.01+fact.nw$mean.HOI.diff[i], 0.09+fact.nw$mean.HOI.diff[i])
  }
  
  parameters<-list(time,species,interaction,alpha,HOI_a)
  names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")
  data.fullfilling.coexistence.criteria<-dynamics.data(parms = parameters)
  
  
  fact.nw$Species_richness[i] = length(which(data.fullfilling.coexistence.criteria$Density[3000,]>0))
  
  print(i)

  
}
g1<-ggplot(fact.nw, 
           aes(x=factor(round(mean.HOI.diff,3)),y=Species_richness))+
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_classic()+xlab("d")+ylab("Species richness")+
  ggtitle("B) Weyl's inequality NOT satisfied")

g1

plot_grid(g0,g1)





# Weyl's inequality satisfied and range of parameter space explored of positive HOIs----

d<-0.4
parms<-list()
time =3000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model

mean.d<-seq(-0.0002, 0.0002, 0.0001)

fact.p<- expand.grid(d=d, mean.HOI.diff= mean.d,
                     `random_seed`=4327+(1:30)*100) %>%
  as_tibble %>%
  mutate(`Species_richness`=0,
         `Weyls`= "")


# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for(i in 1: nrow(fact.p)){
  set.seed(fact.p$random_seed[i])
  alpha=matrix(NA, nrow = species, ncol = species) #interspecific
  temp<-runif(length(upperTriangle(alpha,diag = F)),0.04,0.09)
  upperTriangle(alpha,diag = F,byrow = T)<-temp
  alpha<-t(alpha)
  upperTriangle(alpha,diag = F,byrow = T)<-temp
  
  diag(alpha)<- fact.p$d[i] # intraspecfic effec
  for(k in 1:species){
    HOI_a[k,,]<- matrix(runif(species^2, 0.0005,0.0007), nrow=species,ncol=species) 
    HOI_a[k,k,]<- runif(species,0.0005+fact.pw$mean.HOI.diff[i], 0.0007+fact.pw$mean.HOI.diff[i])  }
  
  parameters<-list(time,species,interaction,alpha,HOI_a)
  names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")
  data.fullfilling.coexistence.criteria<-dynamics.data(parms = parameters)
  
  
  fact.p$Species_richness[i] = length(which(data.fullfilling.coexistence.criteria$Density[3000,]>0))
  fact.p$Weyls[i]=data.fullfilling.coexistence.criteria$Coexistence.criteria
  print(i)
  
  
}



a0<-ggplot(fact.p, 
           aes(x=factor(mean.HOI.diff),y=Species_richness))+
  geom_boxplot()+geom_jitter(shape=16, position=position_dodge(1))+ylim(c(0,50))+
  theme_classic()+xlab("d")+ylab("Species richness")+ggtitle("B) Weyl's inequality satisfied")

a0



# Weyl's inequality NOT satisfied and range of parameter space explored of positive HOIs----

d<-0.01
parms<-list()
time =3000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model
mean.d<-seq(-0.0002, 0.0002, 0.0001)



fact.pw<- expand.grid(d=d, mean.HOI.diff= mean.d,
                      `random_seed`=4327+(1:30)*100) %>%
  as_tibble %>%
  mutate(`Species_richness`=0,
         `Weyls` = "")


# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for(i in 1: nrow(fact.pw)){
  set.seed(fact.pw$random_seed[i])
  alpha=matrix(NA, nrow = species, ncol = species) #interspecific
  temp<-runif(length(upperTriangle(alpha,diag = F)),0.04,0.09)
  upperTriangle(alpha,diag = F,byrow = T)<-temp
  alpha<-t(alpha)
  upperTriangle(alpha,diag = F,byrow = T)<-temp
  
  diag(alpha)<- fact.pw$d[i] # intraspecfic effec
  for(k in 1:species){
    HOI_a[k,,]<- matrix(runif(species^2, 0.0005,0.0007), nrow=species,ncol=species) 
    HOI_a[k,k,]<- runif(species,0.0005+fact.pw$mean.HOI.diff[i], 0.0007+fact.pw$mean.HOI.diff[i])
  }
  
  parameters<-list(time,species,interaction,alpha,HOI_a)
  names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")
  data.fullfilling.coexistence.criteria<-dynamics.data(parms = parameters)
  
  
  fact.pw$Species_richness[i] = length(which(data.fullfilling.coexistence.criteria$Density[3000,]>0))
  fact.pw$Weyls[i] =data.fullfilling.coexistence.criteria$Coexistence.criteria
  print(i)
 
  
}


a1<-ggplot(fact.pw, 
           aes(x=factor(mean.HOI.diff),y=Species_richness))+
  geom_boxplot()+geom_jitter(shape=16, position=position_dodge(.5))+ylim(c(0,50))+
  theme_classic()+xlab("d")+ylab("Species richness")+ggtitle("B) Weyl's inequality satisfied")

a1


plot_grid(a0,a1)


