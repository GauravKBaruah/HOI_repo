

# R script for HOI analysis by **Pragya Singh and Gaurav** for the paper: "Higher-order interactions and species coexistence." 2020 for Theoretical Ecology.
# R scripts for reproducing figure 1 in the main-text.

rm(list=ls())
library(entropy)
library(ggplot2)
library(cowplot)

#function for calculating the invasion growth rates
invasion.growth.rate<-function(higher_order_matrices, theta_inv, Neq){
  
  
  new.matr<-0.01*matrix(c(2*sqrt(theta_inv), theta_inv, 0 ,
                          1, 2*sqrt(theta_inv), 0,
                          0, 0, 10), nrow=3, ncol=3, byrow = T) #this is equation 6 in the main-text
  
  #species 1 invasion growth rate for three-way HOIs
  num.sp1 <- new.matr[1,2]*(new.matr[2,2]-higher_order_matrices$Sp2[2,3]*Neq) -
    higher_order_matrices$Sp1[2,2]-higher_order_matrices$Sp1[2,3]*Neq*(new.matr[2,2]-higher_order_matrices$Sp2[2,3]*Neq)
  
  
  denom.sp1 <- (new.matr[2,2]-higher_order_matrices$Sp2[2,3]*Neq)^2
  
  
  #invasion growth rate of species 1 and species 2 for pairwise interaction case 
  rA <- 1- new.matr[1,2]/new.matr[2,2]  
  rD<- 1-new.matr[2,1]/new.matr[1,1]
  
  #species 2 invasion growth rate for three-way HOIs
  num.sp2 <- new.matr[2,1]*(new.matr[1,1]-higher_order_matrices$Sp1[1,3]*Neq) -
    higher_order_matrices$Sp2[1,1]-
    higher_order_matrices$Sp2[1,3]*Neq*(new.matr[1,1]-higher_order_matrices$Sp1[1,3]*Neq)
  denom.sp2 <- (new.matr[1,1]-higher_order_matrices$Sp1[1,3]*Neq)^2
  
  rB <- 1 -  num.sp1/denom.sp1 
  rC <- 1 - num.sp2/denom.sp2
  
  output<-(list(rA,rD,rB, rC))
  names(output) <-c("Invasion_Growth_Rate_pairwise.sp1",
                    "Invasion_Growth_Rate_pairwise.sp2",
                    "Invasion_Growth_Rate_Sp1", "Invasion_sp2")
  
 
  return(output)
}


#function that takes in three parameters and returns invasion growth rates of species 1 and 2 with pairwise and HOI interactions
#hh_alpha1 = HOI matrix for species 1
##hh_alpha2 = HOI matrix for species 2
## alp = a character name for which HOI term is being perturbed from the HOI matrix
different.alpha.invasion.growth.rate<-function(hh_alpha1, hh_alpha2,alp){

species<-3#species

# HOI matrix for species 3 which is always 0
hh_alpha3<-matrix(c(0, 0, 0 ,                                       # a311,312,313
                    0, 0, 0 ,                                       # a321,322,323
                    0, 0, 0) , nrow=3, ncol = 3, byrow = T)         # a331,332,333




higher_order_matrix<- list(hh_alpha1,  hh_alpha2,hh_alpha3)
names(higher_order_matrix) <- c("Sp1", "Sp2", "Sp3")

theta.seq<-seq(0.1,10, 0.01) # this is the range of \theta values for equation 6 in the main text
rr.sp1<-rr.sp2<-rr.pariwise.sp1<-rr.pariwise.sp2<-numeric()

# calculating invasion growth rates with and without HOIs for each theta values
for(i in 1:length(theta.seq))
  
{
  
  rr.sp1[i]<- invasion.growth.rate(higher_order_matrices =higher_order_matrix, theta_inv=theta.seq[i],Neq = 10)$Invasion_Growth_Rate_Sp1
  rr.sp2[i]<- invasion.growth.rate(higher_order_matrices =higher_order_matrix, theta_inv=theta.seq[i],Neq = 10)$Invasion_sp2
  rr.pariwise.sp1[i]<-invasion.growth.rate(higher_order_matrices =higher_order_matrix, theta_inv=theta.seq[i],Neq = 0)$Invasion_Growth_Rate_pairwise.sp1
  rr.pariwise.sp2[i]<-invasion.growth.rate(higher_order_matrices =higher_order_matrix, theta_inv=theta.seq[i],Neq = 0)$Invasion_Growth_Rate_pairwise.sp2
  
}

HOI_3<-data.frame(rbind((cbind(theta.seq,rr.pariwise.sp1,'rr.pariwise.sp1','sp1','pairwise',alp)),
                        (cbind(theta.seq,rr.pariwise.sp2,'rr.pariwise.sp2','sp2','pairwise',alp)),
                        (cbind(theta.seq,rr.sp1,'rr.sp1','sp1','3rd order',alp)),
                        (cbind(theta.seq,rr.sp2,'rr.sp2','sp2','3rd order',alp))))
names(HOI_3)<-c('theta','IGRate','type','species','interaction','hoi.alpha')

output<-HOI_3

return(output)
}

sign=-1



# for the symmetric case which means all the HOI terms have the same value/strength
hh_alpha1a<-sign*matrix(c(0, 0.01, 0.01 ,                                   # a111,112,113
                    0.01, 0.01, 0.01 ,                                # a121,122,123   
                    0, 0, 0) , nrow=3, ncol = 3, byrow = T)        # a131,132,133

hh_alpha2a<-sign*matrix(c(0.01, 0.01, 0.01,                                  # a211,212,213
                    0.01, 0, 0.01,                                     # a221,222,223
                    0, 0,0) , nrow=3, ncol = 3, byrow = T)          # a231,232,233

HOI_3a<-different.alpha.invasion.growth.rate(
  hh_alpha1 = hh_alpha1a,hh_alpha2 = hh_alpha2a,alp = 'symmetric')




da<- -0.01 #by how much each HOI term is being perturbed

#initialising the HOI matrices for species 1 and 2
hh_alpha1a<-1*sign*matrix(c(0, 0, 0 ,                                   # a111,112,113
                          0, 0, 0 ,                                # a121,122,123   
                          0, 0, 0) , nrow=3, ncol = 3, byrow = T)        # a131,132,133

hh_alpha2a<-1*sign*matrix(c(0, 0, 0,                                  # a211,212,213
                          0, 0, 0,                                     # a221,222,223
                          0, 0,0) , nrow=3, ncol = 3, byrow = T)          # a231,232,233

# only a certain HOI term has some strength, rest are zero otherwise
hh_alpha2b<-hh_alpha2a
hh_alpha2b[2,3]<-hh_alpha2b[2,3]+da

hh_alpha1c<-hh_alpha1a
hh_alpha1c[2,2]<-hh_alpha1c[2,2]+da

hh_alpha1d<-hh_alpha1a
hh_alpha1d[2,3]<-hh_alpha1d[2,3]+da

hh_alpha1e<-hh_alpha1a
hh_alpha1e[1,3]<-hh_alpha1e[1,3]+da

hh_alpha2f<-hh_alpha2a
hh_alpha2f[1,1]<-hh_alpha2f[1,1]+da

hh_alpha2g<-hh_alpha2a
hh_alpha2g[1,3]<-hh_alpha2g[1,3]+da


# invasion growth rates in the presence of one HOI term
HOI_3b<-different.alpha.invasion.growth.rate(
  hh_alpha1 = hh_alpha1a,hh_alpha2 = hh_alpha2b,alp = paste('beta_223',hh_alpha2b[2,3]))
HOI_3c<-different.alpha.invasion.growth.rate(
  hh_alpha1 = hh_alpha1c,hh_alpha2 = hh_alpha2a,alp = paste('beta_122',hh_alpha1c[2,2]))
HOI_3d<-different.alpha.invasion.growth.rate(
  hh_alpha1 = hh_alpha1d,hh_alpha2 = hh_alpha2a,alp = paste('beta_123',hh_alpha1d[2,3]))
HOI_3e<-different.alpha.invasion.growth.rate(
  hh_alpha1 = hh_alpha1e,hh_alpha2 = hh_alpha2a,alp = paste('beta_113',hh_alpha1e[1,3]))
HOI_3f<-different.alpha.invasion.growth.rate(
  hh_alpha1 = hh_alpha1a,hh_alpha2 = hh_alpha2f,alp = paste('beta_211',hh_alpha2f[1,1]))
HOI_3g<-different.alpha.invasion.growth.rate(
  hh_alpha1 = hh_alpha1a,hh_alpha2 = hh_alpha2g,alp = paste('beta_213',hh_alpha2g[1,3]))
if(hh_alpha1e[1,3]==hh_alpha2b[2,3])
{HOI_3h<-different.alpha.invasion.growth.rate(
  hh_alpha1 = hh_alpha1e,hh_alpha2 = hh_alpha2b,alp = paste('beta_113_beta_223 ',hh_alpha2b[2,3]))} else{
  HOI_3h<-different.alpha.invasion.growth.rate(
    hh_alpha1 = hh_alpha1e,hh_alpha2 = hh_alpha2b,alp = 
      paste('beta_113_beta_223',hh_alpha1e[1,3],hh_alpha2b[2,3]))
}




#color paletter
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#creating the data frame
HOI_final<-data.frame(rbind(HOI_3a,HOI_3b,HOI_3c,
                            HOI_3c,HOI_3d,HOI_3e,HOI_3f,HOI_3g,HOI_3h))
names(HOI_final)<-c('theta','IGRate','type','species','interaction','hoi.alpha')


HOI_final$IGRate<-as.numeric(levels(HOI_final$IGRate))[HOI_final$IGRate]
HOI_final$theta<-as.numeric(levels(HOI_final$theta))[HOI_final$theta]
HOI_final$type<-as.factor(HOI_final$type)
HOI_final$species<-as.factor(HOI_final$species)
HOI_final$order<-as.factor(HOI_final$interaction)
HOI_final$hoi.alpha<-as.factor(HOI_final$hoi.alpha)

HOI_final_sp1<-subset(HOI_final,HOI_final$species=='sp1')
HOI_final_sp2<-subset(HOI_final,HOI_final$species=='sp2')



p<-ggplot(subset(HOI_final_sp1, hoi.alpha %in% c("symmetric","beta_223 -0.01", "beta_122 -0.01",
                                                         "beta_123 -0.01")), aes(x=(sqrt(theta)), y=IGRate)) +
  geom_line(aes(linetype=interaction),size=1,color="#999999")+
  ggtitle('Species 1')+xlab("Fitness difference")+theme_classic()+theme(legend.position="top")+
  geom_hline(yintercept = 0,linetype='dashed',col='red')+ylab("Invasion growth rate")+labs(color = "Type of interaction")#+
  # scale_linetype_manual(values =c("dotdash","dotdash","solid","solid"))+

mpg3 <- transform(subset(HOI_final_sp1, hoi.alpha %in% c("symmetric","beta_223 -0.01", "beta_122 -0.01",
                                                         "beta_123 -0.01")),
                  hoi.alpha = factor(hoi.alpha, levels=c("symmetric","beta_223 -0.01", "beta_122 -0.01",
                                                         "beta_123 -0.01"), 
                                     labels=c("symmetric", 
                                              "beta[223] == -0.01", 
                                              "beta[122] == -0.01",
                                              "beta[123] == -0.01")))
sp11<-(p %+% mpg3) + facet_wrap(.~ hoi.alpha,ncol=2,nrow=2,
                          scales = 'free', labeller=label_parsed)


sp11

# for species 2

sp22<-ggplot(subset(HOI_final_sp2, hoi.alpha %in% c("symmetric", "beta_113 -0.01",
                                                    "beta_211 -0.01",
                                                    "beta_213 -0.01")),
             aes(x=(1/sqrt(theta)), y=IGRate)) + facet_wrap(.~hoi.alpha,ncol=4,
                                                                     scales = 'free') +
  geom_line(aes(linetype=interaction),size=1,color="#999999")+
  ggtitle('Species 2')+xlab("Fitness difference")+theme_classic()+
  geom_hline(yintercept = 0,linetype='dashed',col='red')+ylab("Invasion growth rate")+
  theme(legend.position=c(0.96, 0.7))+labs(color = "Type of interaction")
 



mpg5 <- transform(subset(HOI_final_sp2, hoi.alpha %in% c("symmetric", "beta_113 -0.01",
                                                         "beta_211 -0.01",
                                                         "beta_213 -0.01")),
                  hoi.alpha = factor(hoi.alpha,  levels=c("symmetric",
                                                          "beta_113 -0.01",
                                                          "beta_211 -0.01",
                                                          "beta_213 -0.01"), 
                                     labels=c("symmetric",  "beta[113] == -0.01",
                                              "beta[211] == -0.01",
                                              "beta[213] == -0.01"
                                     )))
sp222<-(sp22 %+% mpg5) + facet_wrap(.~ hoi.alpha,ncol=2,nrow=2,
                          scales = 'free', labeller=label_parsed)

sp222
