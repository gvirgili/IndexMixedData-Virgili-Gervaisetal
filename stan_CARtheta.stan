// Mixed factor analysis
// Using dichotomous and continuous entries
// One factor
// Is equivalent to a mix of IRT probit model and factor analysis.

//The latent factor is thought to vary amongst EAs (enumeration areas), 
// but not between households within EAS.
// y_ijk = a_k + B_k ' theta_i + e_ijk

//theta follows an exact CAR distribution with tau set to 1.

//____The continuous outcomes MUST BE CENTERED and come FIRST______

//returns the beta coefficients (discrimination parameters/factor loadings).
//returns the alpha intercept (item difficulty).
//returns theta, the latent factor for each observation.


functions {
  /**Taken from https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}




data {
  
  int <lower=1> nEa ; //number of EAS
  int <lower=1> nObs ; //number of observations (households)
  int <lower=1> nVar ; //number of variables observed
  int <lower=1> nContinuous ; //number of continuous variables observed

  
  int <lower=1, upper=nEa> ea[nObs]; //EA for observation n
  real x_cont[nObs,nContinuous]; //observed continuous variables
  int <lower=0,upper=1> x_dichot[nObs,nVar-nContinuous];//observed dichotomous variables

  matrix <lower = 0, upper = 1> [nEa, nEa] W;//adjencency matrix
  matrix <lower = 0> [nEa, nEa] D;//diagonal matrix with number of neighboors
  int <lower = 0> W_n;// number of adjacent region pairs

}

transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[nEa] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[nEa] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(nEa - 1)) {
      for (j in (i + 1):nEa) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:nEa) D_sparse[i] = sum(W[i]);
  {
    vector[nEa] invsqrtD;  
    for (i in 1:nEa) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}


parameters{
  vector[nVar-nContinuous] alpha; //difficulty parameters
  vector<lower=0> [nVar] beta; //discrimination parameters
  vector[nEa] theta; //factor scores for each EA
  vector <lower=0> [nContinuous] sigma_continuous;//scale parameters for normals of continuous variables.
  real<lower = 0, upper = 1> iota_spatial;//controls spatial dependence.
}


model {
  matrix[nObs,nVar-nContinuous] x_star_dichot; // continuous latent form of the dichotomous observations.
  matrix[nObs,nContinuous] x_star_cont;// continuous latent form of the continuous observations.
  matrix[nObs,nVar] Pi; //matrix of probabilities for the observed dichotomous data.

  theta ~ sparse_car(1, iota_spatial, W_sparse, D_sparse, lambda, nEa, W_n);
  alpha ~ std_normal();
  beta ~ std_normal();
  sigma_continuous~cauchy(0,1);
    //iota has a uniform 0-1 prior

  for (n in 1:nObs){
    
    for (k in 1:nContinuous){
      x_star_cont[n,k] = beta[k]*theta[ea[n]];
      x_cont[n,k] ~ normal(x_star_cont[n,k], sigma_continuous[k]);
    }
    
    for (k in (nContinuous+1):nVar){
      x_star_dichot[n,k-nContinuous] = alpha[k-nContinuous] + beta[k]*theta[ea[n]];
      Pi[n,k-nContinuous] = Phi(x_star_dichot[n,k-nContinuous]);
      x_dichot[n,k-nContinuous] ~ bernoulli(Pi[n,k-nContinuous]);
    }
    
  }

}








