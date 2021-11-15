//
// Author: Michele Nguyen
// Date: 2019-12-27
//
// Joint spatial ordinal regression
//
// Data: Damage state data with coordinates and log(PGA) values.
//
// Note: At the moment, using a frequentist approach and holding non-spatial ordinal parameters as fixed.

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{

  using namespace density;
  using namespace Eigen;

  // ------------------------------------------------------------------------ //
  // Spatial field data
  // ------------------------------------------------------------------------ //

  // Distance matrices:
  DATA_MATRIX(D);

  // ------------------------------------------------------------------------ //
  // Input data
  // ------------------------------------------------------------------------ //

  // Covariates: log(PGA) data
  DATA_VECTOR(x1);
  DATA_VECTOR(x2);

  // Response: Damage state point data
  DATA_MATRIX(damage_ind_1); // Indicator matrix (no_obs x no_state).
  DATA_MATRIX(damage_lower_1); // Indicator matrix for one damage state down (if applicable).
  DATA_VECTOR(damage_state_1);
  DATA_MATRIX(damage_ind_2); // Indicator matrix (no_obs x no_state).
  DATA_MATRIX(damage_lower_2); // Indicator matrix for one damage state down (if applicable).
  DATA_VECTOR(damage_state_2);
  DATA_SCALAR(no_states);

  // ------------------------------------------------------------------------ //
  // Parameters
  // ------------------------------------------------------------------------ //

  // Matern covariance parameters
  PARAMETER(log_phi);
  PARAMETER(log_sigma_2);
  PARAMETER(log_tau_2);
  PARAMETER(log_tau1_2);
  PARAMETER(log_tau2_2);

  DATA_SCALAR(log_kappa);

  // Ordinal model parameters
  PARAMETER(log_slope1);
  PARAMETER(log_slope2);
  PARAMETER_VECTOR(c_factor1); // Cut-off factors (no. of damage states - 1).
  PARAMETER_VECTOR(c_factor2);

  // Convert hyperparameters to natural scale
  // Matern function of geoR package: phi = range; kappa = smoothness; tau_2 = nugget; sigma_2 = partial sill.
  Type phi = exp(log_phi);
  Type kappa = exp(log_kappa);
  Type tau_2 = exp(log_tau_2);
  Type sigma_2 = exp(log_sigma_2);
  Type tau1_2 = exp(log_tau1_2);;
  Type tau2_2 = exp(log_tau2_2);

  Type slope1 = exp(log_slope1);
  Type slope2 = exp(log_slope2);

 // Convert cutoff factors to cutoffs:

 vector<Type> cutoffs1 = c_factor1;
 for(int i=1; i<(no_states-1); i++) {
   for(int j=0; j<i; j++) {
     cutoffs1[i] += c_factor1[j];
   }
 }

 vector<Type> cutoffs2 = c_factor2;
 for(int i=1; i<(no_states-1); i++) {
   for(int j=0; j<i; j++) {
     cutoffs2[i] += c_factor2[j];
   }
 }

  PARAMETER_VECTOR(field);

  // Number of points per category (note at the moment, fixed to be the same - save loops)
  int n_points = x1.size();

  Type nll = Type(0.0);

  // ------------------------------------------------------------------------ //
  // Likelihoods of the random fields.
  // ------------------------------------------------------------------------ //

  matrix<Type> C(D);
  for(int i=0; i<C.rows(); i++) {
    for(int j=0; j<C.cols(); j++){
      C(i,j) = sigma_2*matern(D(i,j), phi, kappa);
      if(i==j){
        C(i,j) += tau_2;
      }
    }
  }

 // Evaluate the negative log density of a mean zero multivariate Gaussian variable with general covariance matrix Sigma.
 nll += density::MVNORM_t<Type>(C)(field);

  // ------------------------------------------------------------------------ //
  // Likelihood from data
  // ------------------------------------------------------------------------ //

  vector<Type> pixel_latent1(n_points);
  vector<Type> reportnll1(n_points);
  vector<Type> pixel_latent2(n_points);
  vector<Type> reportnll2(n_points);
  vector<Type> field_subset1(n_points);
  vector<Type> field_subset2(n_points);

  field_subset1 = field.segment(0, n_points);
  field_subset2 = field.segment(n_points, n_points);
  pixel_latent1 = (x1 + field_subset1.array()) * slope1;
  pixel_latent2 = (x2 + field_subset2.array()) * slope2;

  Type obs_state1 = 0;
  Type upper_prob1 = 0.0;
  Type lower_prob1 = 0.0;
  vector<Type> upper_cut1(n_points);
  vector<Type> lower_cut1(n_points);

  Type obs_state2 = 0;
  Type upper_prob2 = 0.0;
  Type lower_prob2 = 0.0;
  vector<Type> upper_cut2(n_points);
  vector<Type> lower_cut2(n_points);

  upper_cut1 = damage_ind_1*cutoffs1;
  lower_cut1 = damage_lower_1*cutoffs1;
  upper_cut2 = damage_ind_2*cutoffs2;
  lower_cut2 = damage_lower_2*cutoffs2;

  for (int damage_point = 0; damage_point < n_points; damage_point++) {
    obs_state1 = damage_state_1[damage_point];
    obs_state2 = damage_state_2[damage_point];
    if (obs_state1 == no_states){
      upper_prob1 = 1;
    }
    else {
      upper_prob1 = pnorm((upper_cut1[damage_point] - pixel_latent1[damage_point])/sqrt(tau1_2));
    }
    if (obs_state2 == no_states){
      upper_prob2 = 1;
    }
    else {
      upper_prob2 = pnorm((upper_cut2[damage_point] - pixel_latent2[damage_point])/sqrt(tau2_2));
    }
    if (obs_state1 == 1){
     nll -= log(upper_prob1);
     reportnll1[damage_point] = -log(upper_prob1);
    }
    else {
      lower_prob1 = pnorm((lower_cut1[damage_point] - pixel_latent1[damage_point])/sqrt(tau1_2));
      nll -= log(upper_prob1 - lower_prob1);
      reportnll1[damage_point] = -log(upper_prob1 - lower_prob1);
    }
    if (obs_state2 == 1){
      nll -= log(upper_prob2);
      reportnll2[damage_point] = -log(upper_prob2);
    }
    else {
      lower_prob2 = pnorm((lower_cut2[damage_point] - pixel_latent2[damage_point])/sqrt(tau2_2));
      nll -= log(upper_prob2 - lower_prob2);
      reportnll2[damage_point] = -log(upper_prob2 - lower_prob2);
    }
  }

  REPORT(field);
  REPORT(reportnll1);
  REPORT(reportnll2);
  REPORT(nll);

  return nll;
}
