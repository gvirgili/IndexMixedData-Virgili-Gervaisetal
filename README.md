# IndexMixedData-Virgili-Gervaisetal


## Requirements
- [R](https://www.r-project.org) [Free]
- [Stan](https://mc-stan.org/) [Free]

## Contents
- `transformed_data.RData` The Ghana data set which has been previously centered and logged where necessary. 
- `adjacency_matrix.RData` The adjacency matrix of the enumerations areas in the data set. It is only required when running CAR models.
- `number_neighboors_adj.RData` The number of adjacent areas for each area in the data set.  It is only required when running CAR models.
- `run_stan.R` An exemple of code to run the Stan models when using a CAR model.
- `stan_CARtheta_HierarAlpha.stan` The final proposed model.
- `stan_CARtheta.stan` The CAR model with a normal distribution for the difficulty parameters.
- `stan_HierarAlpha.stan` The non-spatial prior model with hierarchical difficulty parameters.
- `stan_simplest.stan` The non-spatial, non-hierarhical model.
