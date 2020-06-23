# R script for HOI analysis by **Pragya Singh and Gaurav** for the paper: "Higher-order interactions and species coexistence." 2020 for Theoretical Ecology.
# R scripts for reproducing analysis.



rm(list=ls())
library(ggplot2)
library(cowplot)

invasion.growth.rate<-function(higher_order_matrices, theta_inv, Neq){
  n3=Neq #equilibrium density for species 3
  n2=5 #equilibrium density for species 2
  n1=5 #equilibrium density for species 1
  hh_alpha1<-matrix(c(0, 0.01, 0.01 ,                                   # a111,112,113
                      0.01, 0.01, 0.01 ,                                 # a121,122,123   
                      0, 0, 0) , nrow=3, ncol = 3, byrow = T)         # a131,132,133
  
  hh_alpha2<-matrix(c(0.01, 0.01, 0.01,                                  # a211,212,213
                      0.01, 0, 0.01,                                    # a221,222,223
                      0, 0,0) , nrow=3, ncol = 3, byrow = T)          # a231,232,233
  
  
  hh_alpha3<-matrix(c(0, 0, 0 ,                                       # a311,312,313
                      0, 0, 0 ,                                       # a321,322,323
                      0, 0, 0) , nrow=3, ncol = 3, byrow = T)         # a331,332,333
  
  #equation 6 matrix where x = 0.01 here.
  new.matr<-0.01*matrix(c(2*sqrt(theta_inv), theta_inv, 0 ,
                          1, 2*sqrt(theta_inv), 0,
                          0, 0, 20), nrow=3, ncol=3, byrow = T)
  alpha<-new.matr
  
  # pairwise
  rA <- 1- new.matr[1,2]/new.matr[2,2]   
  rD<- 1-new.matr[2,1]/new.matr[1,1]
  
  #species 1 with four way HOI
  rB=1-(n2^2*(hh_alpha1[2,2]-n3*(h4_alpha12[2,3]+h4_alpha12[3,2])))-(n2*alpha[1,2])+
    ((n2^3)*h4_alpha12[2,2])+(n2*n3*(hh_alpha1[2,3]+n3*h4_alpha12[3,3]))
  
  #species 2 with four way HOI
  
  rC=1-(n1^2*(hh_alpha2[1,1]-n3*(h4_alpha21[1,3]+h4_alpha21[3,1])))-(n1*alpha[2,1])+
    ((n1^3)*h4_alpha21[1,1])+(n1*n3*(hh_alpha2[1,3]+n3*h4_alpha21[3,3]))
  
  
  output<-(list(rA,rD,rB, rC))
  names(output) <-c("Invasion_Growth_Rate_pairwise.sp1",
                    "Invasion_Growth_Rate_pairwise.sp2",
                    "Invasion_Growth_Rate_Sp1", "Invasion_sp2")
  
  return(output)
}


#takes four way HOI matrix for species 1 and species 2
#alp= four-way HOI character value that is being perturbed
invasion.growth.rate.fourth.order<-function(h4_alpha12,  h4_alpha21,alp){
  
  
  
  
  
  
  higher_order_matrix<- list(h4_alpha12,h4_alpha21)
  names(higher_order_matrix) <- c('h4_alpha12','h4_alpha21')
  
  theta.seq<-seq(0.1,9, 0.01) #sequence of theta values
  rr.sp1<-rr.sp2<-rr.pariwise<-rr.pariwise.sp1<-rr.pariwise.sp2<-numeric()
  for(i in 1:length(theta.seq))
  {
    
    rr.sp1[i]<- invasion.growth.rate(higher_order_matrices =higher_order_matrix, theta_inv=theta.seq[i],Neq = 5)$Invasion_Growth_Rate_Sp1
    rr.sp2[i]<- invasion.growth.rate(higher_order_matrices =higher_order_matrix, theta_inv=theta.seq[i],Neq = 5)$Invasion_sp2
    rr.pariwise.sp1[i]<-invasion.growth.rate(higher_order_matrices =higher_order_matrix, theta_inv=theta.seq[i],Neq = 0)$Invasion_Growth_Rate_pairwise.sp1
    rr.pariwise.sp2[i]<-invasion.growth.rate(higher_order_matrices =higher_order_matrix, theta_inv=theta.seq[i],Neq = 0)$Invasion_Growth_Rate_pairwise.sp2
    
  }
  
  HOI_4<-data.frame(rbind((cbind(theta.seq,rr.pariwise.sp1,'rr.pariwise.sp1','sp1','pairwise',alp)),
                          (cbind(theta.seq,rr.pariwise.sp2,'rr.pariwise.sp2','sp2','pairwise',alp)),
                          (cbind(theta.seq,rr.sp1,'rr.sp1','sp1','4th order',alp)),
                          (cbind(theta.seq,rr.sp2,'rr.sp2','sp2','4th order',alp))))
  names(HOI_4)<-c('theta','IGRate','type','species','interaction','hoi.alpha')
  
  output<-HOI_4
  
  return(output)
}

#positive sign for positive HOI
sign=-1
#symmetric case i.e. all HOIs have the same value
h4_alpha12<-1*sign*matrix(c(0.01, 0.01, 0.01 ,                                # a1211,1212,1213
                            0.01, 0.01, 0.01 ,                                 # a1221,1222,1223   
                            0.01, 0.01, 0.01 ) , nrow=3, ncol = 3, byrow = T)  # a1231,1232,1233


h4_alpha21<-1*sign*matrix(c(0.01, 0.01, 0.01,                                 # a2111,2112,2113
                            0.01, 0.01, 0.01,                                 # a2121,2122,2123
                            0.01, 0.01, 0.01) , nrow=3, ncol = 3, byrow = T)   # a2131,2132,2133



HOI_4a<-invasion.growth.rate.fourth.order(
  h4_alpha12 = h4_alpha12,h4_alpha21 = h4_alpha21,alp = "sym")




# now all the HOI parameters being perturbed one by one
h4_alpha12<-1*sign*matrix(c(0, 0, 0 ,                                # a1211,1212,1213
                            0, 0, 0 ,                                 # a1221,1222,1223   
                            0, 0, 0 ) , nrow=3, ncol = 3, byrow = T)  # a1231,1232,1233


h4_alpha21<-1*sign*matrix(c(0, 0, 0,                                 # a2111,2112,2113
                            0, 0, 0,                                 # a2121,2122,2123
                            0, 0, 0) , nrow=3, ncol = 3, byrow = T)   # a2131,2132,2133


da=-0.01
b.h4_alpha12<-h4_alpha12
b.h4_alpha12[2,2]<-h4_alpha12[2,2]+da #parameter 1222 being perturbed by da increment

c.h4_alpha12<-h4_alpha12
c.h4_alpha12[2,3]<-h4_alpha12[2,3]+da

d.h4_alpha12<-h4_alpha12
d.h4_alpha12[3,2]<-h4_alpha12[3,2]+da

e.h4_alpha12<-h4_alpha12
e.h4_alpha12[3,3]<-h4_alpha12[3,3]+da

f.h4_alpha21<-h4_alpha21
f.h4_alpha21[1,1]<-h4_alpha21[1,1]+da

g.h4_alpha21<-h4_alpha21
g.h4_alpha21[1,3]<-h4_alpha21[1,3]+da

h.h4_alpha21<-h4_alpha21
h.h4_alpha21[3,1]<-h4_alpha21[3,1]+da

i.h4_alpha21<-h4_alpha21
i.h4_alpha21[3,3]<-h4_alpha21[3,3]+da

HOI_4b<-invasion.growth.rate.fourth.order(
  h4_alpha12 = b.h4_alpha12,h4_alpha21 = h4_alpha21,alp = paste('1222=',b.h4_alpha12[2,2]))
HOI_4c<-invasion.growth.rate.fourth.order(
  h4_alpha12 = c.h4_alpha12,h4_alpha21 = h4_alpha21,alp = paste('1223=',c.h4_alpha12[2,3]))
HOI_4d<-invasion.growth.rate.fourth.order(
  h4_alpha12 = d.h4_alpha12,h4_alpha21 = h4_alpha21,alp = paste('1232=',d.h4_alpha12[3,2]))
HOI_4e<-invasion.growth.rate.fourth.order(
  h4_alpha12 = e.h4_alpha12,h4_alpha21 = h4_alpha21,alp = paste('1233=',e.h4_alpha12[3,3]))
HOI_4f<-invasion.growth.rate.fourth.order(
  h4_alpha12 = h4_alpha12,h4_alpha21 = f.h4_alpha21,alp = paste('2111=',f.h4_alpha21[1,1]))
HOI_4g<-invasion.growth.rate.fourth.order(
  h4_alpha12 = h4_alpha12,h4_alpha21 = g.h4_alpha21,alp = paste('2113=',g.h4_alpha21[1,3]))
HOI_4h<-invasion.growth.rate.fourth.order(
  h4_alpha12 = h4_alpha12,h4_alpha21 = h.h4_alpha21,alp = paste('2131=',h.h4_alpha21[3,1]))
HOI_4i<-invasion.growth.rate.fourth.order(
  h4_alpha12 = h4_alpha12,h4_alpha21 = i.h4_alpha21,alp = paste('2133=',i.h4_alpha21[3,3]))

# #################################################PLOTTING



HOI_final4thorder<-data.frame(rbind(HOI_4a,HOI_4b,HOI_4c,HOI_4d,HOI_4e,HOI_4f,HOI_4g,HOI_4h,HOI_4i))
names(HOI_final4thorder)<-c('theta','IGRate','type','species','interaction','hoi.alpha')
HOI_final4thorder$IGRate<-as.numeric(levels(HOI_final4thorder$IGRate))[HOI_final4thorder$IGRate]
HOI_final4thorder$theta<-as.numeric(levels(HOI_final4thorder$theta))[HOI_final4thorder$theta]
HOI_final4thorder$type<-as.factor(HOI_final4thorder$type)
HOI_final4thorder$species<-as.factor(HOI_final4thorder$species)
HOI_final4thorder$interaction<-as.factor(HOI_final4thorder$interaction)


#subseting for species 1 and species 2
HOI_final4thorder_sp1<-subset(HOI_final4thorder,HOI_final4thorder$species=='sp1')
HOI_final4thorder_sp2<-subset(HOI_final4thorder,HOI_final4thorder$species=='sp2')


#png("Figure3 sp2 +alpha n3=.png", width = 12, height = 6, units = 'in', res = 300)
#pdf( file = "4th hoi n3=10 +alpha.pdf")
p1<-ggplot(subset(HOI_final4thorder_sp1, hoi.alpha %in% c("sym", "1222= -0.01","1223= -0.01","1232= -0.01","1233= -0.01")), aes(x=(sqrt(theta)), y=IGRate)) + facet_wrap(.~hoi.alpha,ncol=3) +
  geom_line(aes(linetype=interaction),size=1,color="#999999")+geom_hline(yintercept = 0,linetype='dashed',col='red')+ylab("Invasion growth rate")+
  ggtitle('Species 1')+xlab("Fitness difference")+
  theme(legend.position=c(0.8, 0.2))+labs(color = "Type of interaction")

mpg31 <- transform(subset(HOI_final4thorder_sp1, hoi.alpha %in% c("sym", "1222= -0.01","1223= -0.01","1232= -0.01","1233= -0.01")),
                  hoi.alpha = factor(hoi.alpha, levels=c("sym", "1222= -0.01", "1223= -0.01","1232= -0.01","1233= -0.01"), 
                                     labels=c("symmetric","gamma[1222] == -0.01", "gamma[1223] == -0.01","gamma[1232] == -0.01", "gamma[1233] == -0.01"
                                              )))
(p1 %+% mpg31) + facet_wrap(.~hoi.alpha,ncol=4,
                          scales = 'free', labeller=label_parsed)


p2<-ggplot(subset(HOI_final4thorder_sp2, hoi.alpha %in% c("sym", "2111= -0.01","2113= -0.01","2131= -0.01","2133= -0.01")), aes(x=(1/sqrt(theta)), y=IGRate)) + facet_wrap(.~hoi.alpha,ncol=3) +
  geom_line(aes(linetype=interaction),size=1,color="#999999")+geom_hline(yintercept = 0,linetype='dashed',col='red')+ylab("Invasion growth rate")+
  ggtitle('Species 2')+xlab("Fitness difference")+
  theme(legend.position=c(0.8, 0.2))+labs(color = "Type of interaction")


mpg4 <- transform(subset(HOI_final4thorder_sp2, hoi.alpha %in% c("sym", "2111= -0.01","2113= -0.01","2131= -0.01","2133= -0.01")),
                  hoi.alpha = factor(hoi.alpha, levels=c("sym",
                                                         "2111= -0.01",
                                                         "2113= -0.01",
                                                         "2131= -0.01","2133= -0.01"), 
                                     labels=c("symmetric",  "gamma[2111] == -0.01",
                                              "gamma[2113] == -0.01",
                                              "gamma[2131] == -0.01",
                                              "gamma[2133] == -0.01"
                                     )))
(p2 %+% mpg4) + facet_wrap(.~ hoi.alpha,ncol=4,
                          scales = 'free', labeller=label_parsed)


