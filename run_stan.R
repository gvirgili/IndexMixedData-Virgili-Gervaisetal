library(rstan)


#_____setting up the data___________#

load('transformed_data.RData')
data_dichot<-households_ama[,9:23]
data_cont<-households_ama[,4:8]
eas_asNum<-as.numeric(factor(households_ama$ea_code))


#_____Getting the neighbors_________#

load('number_neighboors_adj.RData')
load('adjacency_matrix.RData')

#_______Setting up Stan____________#


stan_data<-list(nEa=max(unique(eas_asNum)), nObs=nrow(data_cont), 
                nVar=ncol(data_cont)+ncol(data_dichot),
                ea=eas_asNum, x_cont=data_cont, x_dichot=data_dichot, 
                nContinuous = ncol(data_cont),
                W=1*adj, D=diag(num_adj), mu_alpha=0, mu_beta=0,
                W_n=sum(1*adj)/2)

options(mc.cores = 2)
stan_fit<-stan(file='stan_CARtheta_HierarAlpha.stan', data=stan_data,cores = 2, chains=2,
               seed=523, it=10000, refresh=100)

save(stan_fit, file='results.RData')
