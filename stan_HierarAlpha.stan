// Mixed factor analysis
// Using dichotomous and continuous entries
// One factor
// Is equivalent to a mix of IRT probit model and factor analysis.

//The latent factor is thought to vary amongst EAs (enumeration areas), 
//but not between households within EAS.
// y_ijk = a_ik + B_k ' theta_i + e_ijk

//Theta follows a normal distribution.
//Alpha has a hierarchical structure.


//____The continuous outcomes MUST BE CENTERED and come FIRST in the data______

//returns the beta coefficients (discrimination parameters/factor loadings).
//returns the alpha intercept (item difficulty) for each EAs.
//returns the overall alpha intercept (overall item difficulty).
//returns theta, the latent factor for each observation.

data {
  
  int <lower=1> nEa ; //number of EAS
  int <lower=1> nObs ; //number of observations (households)
  int <lower=1> nVar ; //number of variables observed
  int <lower=1> nContinuous ; //number of variables observed

  
  int <lower=1, upper=nEa> ea[nObs]; //EA for observation n
  real x_cont[nObs,nContinuous]; //observed continuous variables
  int <lower=0,upper=1> x_dichot[nObs,nVar-nContinuous];//observed dichotomous variables
  
}

parameters{
  matrix[nEa,nVar-nContinuous] alpha; //difficulty parameters, unique to EAs
  vector <lower=0> [nVar] beta; //discrimination parameter
  vector[nEa] theta; //factor scores for each EA
  vector <lower=0> [nContinuous] sigma_continuous;//scale parameters for normals of continuous variables.
  vector[nVar-nContinuous] alpha_star; //overall mean difficulty

}

model {
  matrix[nObs,nVar-nContinuous] x_star_dichot; // continuous latent form of the dichotomous observations.
  matrix[nObs,nContinuous] x_star_cont;// continuous latent form of the continuous observations.
  matrix[nObs,nVar] Pi; //matrix of probabilities for the observed dichotomous data.

  theta ~ std_normal();
  alpha_star ~ std_normal();
  beta ~ std_normal();
  sigma_continuous~cauchy(0,1);


  for (i in 1:nEa){
    for (j in 1:(nVar-nContinuous)) alpha[i,j]~normal(alpha_star[j],1);
  }

  for (n in 1:nObs){
    for (k in 1:nContinuous){
      x_star_cont[n,k] = beta[k]*theta[ea[n]];
      x_cont[n,k] ~ normal(x_star_cont[n,k], sigma_continuous[k]);
    }
    
    for (k in (nContinuous+1):nVar){
      x_star_dichot[n,k-nContinuous] = alpha[ea[n],k-nContinuous] + beta[k]*theta[ea[n]];
      Pi[n,k] = Phi(x_star_dichot[n,k-nContinuous]);
      x_dichot[n,k-nContinuous] ~ bernoulli(Pi[n,k]);
    }
  }
}




