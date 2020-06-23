rm(list=ls())

source('~/Dropbox/HOI_short_paper/HOIcoexistence/JTB or theoretical ecologt/Codes/multispecies_negative_hoi.R')
library(RColorBrewer)
library(gdata)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

### Positive HOIs - Weyls criteria fulfilled ,intrahoi >interhoi ----

parms<-list()
time =5000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model
alpha=matrix(NA, nrow = species, ncol = species) #interspecific
temp<-runif(length(upperTriangle(alpha,diag = F)),0.04,0.09)
upperTriangle(alpha,diag = F,byrow = T)<-temp
alpha<-t(alpha)
upperTriangle(alpha,diag = F,byrow = T)<-temp
isSymmetric(alpha) #interspecific
diag(alpha)<-0.4 # intraspecfic effects


# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for (i in 1:species){
  HOI_a[i,,]<- -matrix(runif(species^2, 0.0003,0.0005), nrow=species,ncol=species) 
  HOI_a[i,i,]<--runif(species,0.0007,0.0009)
}

parameters<-list(time,species,interaction,alpha,HOI_a)
names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")

Weyl_yes_nHOI_1<-dynamics.data(parms = parameters)

par(mfrow=c(2,1))
ts.plot(Weyl_yes_nHOI_1$Density[8:5000,],col=col_vector, ylab="Density",
         main="Weyl's criteria satisfied : Yes")
plot(density(abs(Weyl_yes_nHOI_1$HOI.matrix[1,1,])),
     xlim=c(0.0002,0.00095),ylim=c(0,9000),main=bquote(~ beta[iik] > .(cor) ~ beta[ijk]), ylab="Freq. Distribution")

# here r goes from 1 to 5, where 5 represents 5 species. When the loop runs for 50 species, the figure produced
# becomes very large in size and difficult to work with.
for(r in 1:5){
  for( j in 1:5){
    if(r==j){
      polygon(0,0,main="",xlim=c(0.01,0.1), col='firebrick')
    }
    else if( r!=j){
      polygon(density(abs(Weyl_yes_nHOI_1$HOI.matrix[r,r,])),
              main="",xlim=c(0.0001,0.00095), col='firebrick')
      lines(density(Weyl_yes_nHOI_1$HOI.matrix[r,j,]),
            main="")
      polygon(density(abs(Weyl_yes_nHOI_1$HOI.matrix[r,j,])),
              main="",col = 'salmon2')
    }
  }
}
legend("top", inset= 0.01,legend = c(expression(paste(beta[iik])),
                                     expression(paste(beta[ijk]))), fill=c("firebrick","salmon2"))


### Positive HOIs - Weyls criteria fulfilled ,intrahoi < interhoi ----

parms<-list()
time =5000
species=50
alpha=matrix(NA, nrow = species, ncol = species) #interspecific
temp<-runif(length(upperTriangle(alpha,diag = F)),0.04,0.09)
upperTriangle(alpha,diag = F,byrow = T)<-temp
alpha<-t(alpha)
upperTriangle(alpha,diag = F,byrow = T)<-temp
isSymmetric(alpha) #interspecific
diag(alpha)<-0.4 # intraspecfic effects

# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for (i in 1:species){
  HOI_a[i,,]<- -matrix(runif(species^2, 0.0007,0.0009), nrow=species,ncol=species) 
  HOI_a[i,i,]<--runif(species,0.0003,0.0005)
  
}

parameters<-list(time,species,interaction,alpha,HOI_a)
names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")

Weyl_yes_nHOI_2<-dynamics.data(parms = parameters)

# plotting the snapshots of the simulations

par(mfrow=c(2,1))
ts.plot(Weyl_yes_nHOI_2$Density[8:5000,],col=col_vector, ylab="Density",
        main="Weyl's criteria satisfied : Yes")
plot(density(abs(Weyl_yes_nHOI_2$HOI.matrix[1,1,])),xlim=c(0.0002,0.00095),
     ylim=c(0,8000),main=bquote(~ beta[iik] < .(cor) ~ beta[ijk]), ylab="Freq. Distribution")

for(r in 1:5){
  for( j in 1:5){
    if(r==j){
      polygon(0,0,main="", col='white')
    }
    else if( r!=j){
      polygon(density(abs(Weyl_yes_nHOI_2$HOI.matrix[r,r,])),
              main="", col='firebrick')
      lines(density(abs(Weyl_yes_nHOI_2$HOI.matrix[r,j,])),
            main="")
      polygon(density(abs(Weyl_yes_nHOI_2$HOI.matrix[r,j,])),
              main="",col = 'salmon2')
    }
  }
}
legend("top", inset= 0.01,legend = c(expression(paste(beta[iik])),
                                     expression(paste(beta[ijk]))), fill=c("firebrick","salmon2"))




### Positive HOIs, Weyl's criteria NOT fulfilled but intraHOI> Inter HOI ----
parms<-list()
time =5000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model
alpha=matrix(NA, nrow = species, ncol = species) #interspecific
temp<-runif(length(upperTriangle(alpha,diag = F)),0.01,0.04)
upperTriangle(alpha,diag = F,byrow = T)<-temp
alpha<-t(alpha)
upperTriangle(alpha,diag = F,byrow = T)<-temp
isSymmetric(alpha) #interspecific

diag(alpha)<-0.01 # intraspecfic effects


# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for (i in 1:species){
  HOI_a[i,,]<- -matrix(runif(species^2, 0.0003,0.0005), nrow=species,ncol=species) 
  HOI_a[i,i,]<--runif(species,0.0007,0.0009)
  
}

parameters<-list(time,species,interaction,alpha,HOI_a)
names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")

Weyl_no_nHOI_1<-dynamics.data(parms = parameters)


par(mfrow=c(2,1))
ts.plot(Weyl_no_nHOI_1$Density[1:5000,],col=col_vector,
        ylab="Density", main="Weyl's criteria sastisfied : NO")
plot(density(abs(Weyl_no_nHOI_1$HOI.matrix[1,1,])),
     xlim=c(0.0002,0.00095),ylim=c(0,8000),main=bquote(~ beta[iik] > .(cor) ~ beta[ijk]),
     ylab = "Freq. Distribution")


for(r in 1:5){
  for( j in 1:5){
    if(r==j){
      polygon(0,0,main="",xlim=c(0.0002,0.00095), col='white')
    }
    else if( r!=j){
      polygon(density(abs(Weyl_no_nHOI_1$HOI.matrix[r,r,])),
              main="", col='firebrick')
      lines(density(abs(Weyl_no_nHOI_1$HOI.matrix[r,j,])),
            main="")
      polygon(density(abs(Weyl_no_nHOI_1$HOI.matrix[r,j,])),
              main="",col = 'salmon2')
    }
  }
}
legend("top", inset= 0.01,legend = c(expression(paste(beta[iik])),
                                     expression(paste(beta[ijk]))), fill=c("firebrick","salmon2"))



# Positive HOIs, Weyl's criteria not fulfilled but intraHOI< Inter HOI ----
parms<-list()
time =5000
species=50
interaction= 'HOI' #specificatioin of which type of interaction to model
alpha=matrix(NA, nrow = species, ncol = species) #interspecific
temp<-runif(length(upperTriangle(alpha,diag = F)),0.01,0.04)
upperTriangle(alpha,diag = F,byrow = T)<-temp
alpha<-t(alpha)
upperTriangle(alpha,diag = F,byrow = T)<-temp
isSymmetric(alpha) #interspecific

diag(alpha)<-0.01 # intraspecfic effects


# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for (i in 1:species){
  HOI_a[i,,]<- -matrix(runif(species^2, 0.0007,0.0009), nrow=species,ncol=species) 
  HOI_a[i,i,]<--runif(species,0.0003,0.0005)
  
}

parameters<-list(time,species,interaction,alpha,HOI_a)
names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")

Weyl_no_nHOI_2<-dynamics.data(parms = parameters)


par(mfrow=c(2,1))
ts.plot(Weyl_no_nHOI_2$Density[8:5000,],col=col_vector,
        ylab="Density", main="Weyl's criteria satisfied : NO")
plot(density(abs(Weyl_no_nHOI_2$HOI.matrix[1,1,])),
     xlim=c(0.0001,0.00091),main=bquote(~ beta[iik] < .(cor) ~ beta[ijk]),
     ylab="Freq. Distribution", ylim=c(0,9500))


for(r in 1:5){
  for( j in 1:5){
    if(r==j){
      polygon(0,0,main="",xlim=c(0.01,0.1), col='white')
    }
    else if( r!=j){
      polygon(density(abs(Weyl_no_nHOI_2$HOI.matrix[r,r,])),
              main="", col='firebrick')
      lines(density(abs(Weyl_no_nHOI_2$HOI.matrix[r,j,])),
            main="")
      polygon(density(abs(Weyl_no_nHOI_2$HOI.matrix[r,j,])),
              main="",col = 'salmon2')
    }
  }
}
legend("top", inset= 0.01,legend = c(expression(paste(beta[iik])),
                                     expression(paste(beta[ijk]))), fill=c("firebrick","salmon2"))



# Pairwise interaction, Weyl's criteria fulfilled ----
parms<-list()
time =5000
species=50
interaction= 'pairwise' #specificatioin of which type of interaction to model
alpha=matrix(NA, nrow = species, ncol = species) #interspecific
temp<-runif(length(upperTriangle(alpha,diag = F)),0.01,0.04)
upperTriangle(alpha,diag = F,byrow = T)<-temp
alpha<-t(alpha)
upperTriangle(alpha,diag = F,byrow = T)<-temp
isSymmetric(alpha)

diag(alpha)<-0.4 # intraspecfic effects


# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for (i in 1:species){
  HOI_a[i,,]<- matrix(runif(species^2, 0,0), nrow=species,ncol=species) 
  HOI_a[i,i,]<-runif(species,  0,0)
}

parameters<-list(time,species,interaction,alpha,HOI_a)
names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")

Weyl_yes_pairwise_1<-dynamics.data(parms = parameters)


ts.plot(Weyl_yes_pairwise_1$Density[8:5000,],col=col_vector,
        ylab="Density", main="Weyl's criteria satisfied : YES")


# Pairwise interaction, Weyl's criteria NOT fulfilled ----


parms<-list()
time =5000
species=50
interaction= 'pairwise' #specificatioin of which type of interaction to model
alpha=matrix(NA, nrow = species, ncol = species) #interspecific
temp<-runif(length(upperTriangle(alpha,diag = F)),0.01,0.04)
upperTriangle(alpha,diag = F,byrow = T)<-temp
alpha<-t(alpha)
upperTriangle(alpha,diag = F,byrow = T)<-temp
isSymmetric(alpha) 

diag(alpha)<-0.01 # intraspecfic effects


# hoi terms are picked from random uniform distribution but not used
HOI_a<-array(dim=c(species,species,species))
for (i in 1:species){
  HOI_a[i,,]<- matrix(runif(species^2, 0,0), nrow=species,ncol=species) 
  HOI_a[i,i,]<-runif(species,  0,0)
}

parameters<-list(time,species,interaction,alpha,HOI_a)
names(parameters)<-c("time", "species","interaction", "alpha", "HOI_alpha")

Weyl_no_pairwise_1<-dynamics.data(parms = parameters)


ts.plot(Weyl_no_pairwise_1$Density[8:5000,],col=col_vector,
        ylab="Density", main="Weyl's criteria satisfied : NO")



# plotting timeseries data by Weyl's criteria
par(mfrow=c(3,2))

ts.plot(Weyl_yes_pairwise_1$Density[8:5000,],col=col_vector,
        ylab="Density", main="Weyl's criteria satisfied : YES",xlab="")

ts.plot(Weyl_no_pairwise_1$Density[8:5000,],col=col_vector,
        ylab="Density", main="Weyl's criteria satisfied : NO",xlab="")

ts.plot(Weyl_yes_nHOI_2$Density[8:5000,],col=col_vector,
        ylab="Density",xlab="", main="")

ts.plot(Weyl_no_nHOI_2$Density[8:5000,],col=col_vector,
        ylab="Density",xlab="", main="")

ts.plot(Weyl_yes_nHOI_1$Density[8:5000,],col=col_vector,
        ylab="Density", xlab="Time", main="")


ts.plot(Weyl_no_nHOI_1$Density[8:5000,],col=col_vector,
        ylab="Density", xlab="Time",main="")

