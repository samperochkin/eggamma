#include <Rcpp.h>
using namespace Rcpp;

//* All functions here are constructed so that standard formulas are used
//* whenever nu and xi=1/(nu^2*sigma^2) are outside an epsilon ball around 0 (E0);
//* formulas based on the Taylor series of z=(x/mu)^nu around 0 are used when
//* nu is in E0; and Stirling's approximation of logGamma(xi) when xi is in E0.
//* See inside the first function (dggamma) for some info on the method.

//* Similar formulas of the log-likelihood and gradient based on simple linear 
//* interpolation were texted. However, their use with the 'nlm'
//* function led to errors due to objective/gradient mismatches. 
//* It worked in some instances, but only when specifying a small number 
//* of significant digits. Perhaps using a higher-order (Hermite) polynomial 
//* interpolator would work better. TBD.

// Bernoulli numbers B_even = {B0, B2, B4, ..., B30}
const std::vector<double> B_even = {
  1.0, 1.0 / 6.0, -1.0 / 30.0, 1.0 / 42.0,
  -1.0 / 30.0, 5.0 / 66.0, -691.0 / 2730.0, 7.0 / 6.0,
  -3617.0 / 510.0, 43867.0 / 798.0, -174611.0 / 138.0,
  854513.0 / 2730.0, -236364091.0 / 6.0, 8553103.0 / 870.0,
  -23749461029.0 / 14322.0, 8615841276005.0 / 510.0
};




//------------------------------------------------------------------------------
// plain density (for a single observation) ------------------------------------
//------------------------------------------------------------------------------
// [[Rcpp::export]]
double dggamma_single_cpp(double x, double mu, double sigma, double nu,
                      NumericVector epsilon, int kappa, 
                      bool return_log = true, bool negate_loglik = false) {
  
  // Checked before the function is called
  // if (mu <= 0) stop("mu (in generalized gamma) is not positive");
  // if (sigma <= 0) stop("sigma (in generalized gamma) is not positive");
  
  // Basic quantites
  double sigma2 = sigma * sigma;
  double nu2 = nu * nu;
  double nu_abs = std::abs(nu);
  double l_x_mu = log(x/mu);
  double z = exp(nu * l_x_mu);
  double xi_inv = sigma2 * nu2;
  
  // ll = -log(x) + a + b
  
  // a depends on z, log(z), ...
  // a = w_a*a1 + (1-w_a)*a2, where a1 is standard and a2 its expansion
  
  // b depends on lgamma(xi), log(xi), ...
  // b = w_b*b1 + (1-w_b)*b2, where b1 is standard and b2 uses Stirling's
  
  
  // a part --------------------------------------------------------------------
  double a1 = 0.0; // nu away from zero
  double a2 = 0.0; // nu near zero
  
  if (nu_abs > epsilon[0]) a1 += (log(z) + 1 - z)/xi_inv;
  if (nu_abs < epsilon[1]){
    double log_prod = l_x_mu * l_x_mu / (2 * sigma2);
    a2 -= log_prod; // Part of log-normal log-lik
    for (int k = 3; k <= (nu_abs > 0 ? kappa + 2 : 0); k++) {
      log_prod *= nu * l_x_mu / k;
      a2 -= log_prod;
    }
  } 
  
  double w_a = (nu_abs > epsilon[1]) ? 1 : (nu_abs > epsilon[0] ? (nu_abs-epsilon[0])/(epsilon[1]-epsilon[0]) : 0);
  double a = w_a * a1 + (1-w_a) * a2;
  
  // b part --------------------------------------------------------------------
  double b1 = 0.0; // 1/xi away from zero
  double b2 = 0.0; // 1/xi near zero
  
  if (xi_inv > epsilon[0]) b1 += log(nu_abs) - (log(xi_inv) + 1)/xi_inv - lgamma(1/xi_inv);
  if (xi_inv < epsilon[1]){
    b2 -= log(sigma)+log(2*M_PI)/2; // Part of log-normal log-lik
    
    if(nu_abs > 0){
      double nu_sigma_prod = nu2 * sigma2;
      double nu4_sigma4 = nu_sigma_prod * nu_sigma_prod;
      
      if (kappa >= 2) b2 -= nu_sigma_prod * B_even[1] / 2;
      for (int k = 2; k <= (static_cast<double>(kappa) / 4 + 0.5); k++) {
        nu_sigma_prod *= nu4_sigma4;
        b2 -= nu_sigma_prod * B_even[k] / (2*k * (2*k - 1));
      }
    }
  } 
  
  double w_b = (xi_inv > epsilon[1]) ? 1 : (xi_inv > epsilon[0] ? (xi_inv-epsilon[0])/(epsilon[1]-epsilon[0]) : 0);
  double b = w_b * b1 + (1-w_b) * b2;
  
  
  // Putting it all together ---------------------------------------------------
  double ll = (negate_loglik ? -1 : 1) * (-log(x) + a + b);
  return return_log ? ll : std::exp(ll);
}


//------------------------------------------------------------------------------
// plain density (vectorized, nu and sigma constant) ---------------------------
//------------------------------------------------------------------------------
NumericVector dggamma_vec_cpp(NumericVector x, NumericVector mu, double sigma, double nu,
                          NumericVector epsilon, int kappa, 
                          bool return_log = true, bool negate_loglik = false) {
  
  // Basic quantities
  double sigma2 = sigma * sigma;
  double nu2 = nu * nu;
  double nu_abs = std::abs(nu);
  double xi_inv = sigma2 * nu2;
  double negator = (negate_loglik ? -1 : 1);
  
  // Components that can be computed only once (including the entire b-part)----
  double lgamma_xi_inv = lgamma(1/xi_inv);
  double w_a = (nu_abs > epsilon[1]) ? 1 : (nu_abs > epsilon[0] ? (nu_abs-epsilon[0])/(epsilon[1]-epsilon[0]) : 0);
  
  double b1 = 0.0; double b2 = 0.0;
  if (xi_inv > epsilon[0]){
    b1 = log(nu_abs) - (log(xi_inv) + 1) / xi_inv - lgamma_xi_inv;
  } 
  if (xi_inv < epsilon[1]){
    b2 -= log(sigma)+log(2*M_PI)/2; // Part of log-normal log-lik
    if(nu_abs > 0){
      double nu_sigma_prod = nu2 * sigma2;
      double nu4_sigma4 = nu_sigma_prod * nu_sigma_prod;
      
      if (kappa >= 2) b2 -= nu_sigma_prod * B_even[1] / 2;
      for (int k = 2; k <= (static_cast<double>(kappa) / 4 + 0.5); k++) {
        nu_sigma_prod *= nu4_sigma4;
        b2 -= nu_sigma_prod * B_even[k] / (2*k * (2*k-1));
      }
    }
  }
  double w_b = (xi_inv > epsilon[1]) ? 1 : (xi_inv > epsilon[0] ? (xi_inv-epsilon[0])/(epsilon[1]-epsilon[0]) : 0);
  double b = w_b * b1 + (1-w_b) * b2;
  
  
  
  // Loop ----------------------------------------------------------------------
  NumericVector ll(x.size());
  double a = 0.0; double a1 = 0.0; double a2 = 0.0; 
  double l_x_mu = 0.0;
  double z = 0.0;
  double log_prod = 0.0;
  
  for (int i = 0; i < x.size(); i++) {
    l_x_mu = log(x[i]/mu[i]);
    z = exp(nu * l_x_mu);
    
    // a part ------------------------------------------------------------------
    a1 = 0.0; a2 = 0.0; 
    if (nu_abs > epsilon[0]) a1 += (log(z) + 1 - z) / xi_inv;
    if (nu_abs < epsilon[1]){
      log_prod = l_x_mu * l_x_mu / (2 * sigma2);
      a2 -= log_prod; // Part of log-normal log-lik
      for (int k = 3; k <= (nu_abs > 0 ? kappa + 2 : 0); k++) {
        log_prod *= nu * l_x_mu / k;
        a2 -= log_prod;
      }
    } 
    a = w_a * a1 + (1-w_a) * a2;
    
    ll[i] = negator * (-log(x[i]) + a + b);
    if(!return_log) ll[i] = exp(ll[i]);
  }
  
  return ll;
}


//------------------------------------------------------------------------------
// Cpp wrapper of the density --------------------------------------------------
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector dggamma_cpp(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector nu,
                          NumericVector epsilon, int kappa, 
                          bool return_log = false, bool negate_loglik = false) {
  
  // Different possibilities for the lengths of x, mu, sigma, nu is taken care of in the R wrapper,
  // where some conditions are also checked.
  
  NumericVector ll(x.size());
  
  if(x.size() == sigma.size()){
    // Possibility 1.1 : {x, mu, sigma, nu} were (and are still) all scalars
    // Possibility 1.2 : {x, mu, sigma, nu} were (or were made) vectors of the same length
    for (int i = 0; i < x.size(); i++) 
      ll[i] = dggamma_single_cpp(x[i], mu[i], sigma[i], nu[i], epsilon, kappa, return_log, negate_loglik);
    
  }else if(sigma.size() == 1){
    // Possibility 2 : sigma and nu are scalars, but x and mu are vectors (more efficient)
    ll = dggamma_vec_cpp(x, mu, sigma[0], nu[0], epsilon, kappa, return_log, negate_loglik);
    
  }else{
    // Should not happen
    stop("Error: The varying lengths of the vectors x, mu, sigma, and nu cannot be handled by dggamma_cpp."
           "This should have been caught in R tho...");
  } 
  
  return ll;
}


//------------------------------------------------------------------------------
// negative log-likelihood (w/ gradients and Hessian) --------------------------
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector nllggamma_single_cpp(double x, double mu, double sigma, double nu,
                                   NumericVector epsilon, int kappa) {
  
  // if (mu <= 0) stop("mu (in generalized gamma) is not positive");
  // if (sigma <= 0) stop("sigma (in generalized gamma) is not positive");
  
  
  // setup
  double mu2 = mu * mu;
  double sigma2 = sigma * sigma;
  double nu2 = nu * nu;
  double nu_abs = std::abs(nu);
  double l_x_mu = log(x/mu);
  double z = exp(nu * l_x_mu);
  double xi_inv = sigma2 * nu2;
  
  // components of ll and scores
  double a1 = 0.0; double a1_mu = 0.0; double a1_sigma = 0.0; double a1_nu = 0.0;
  double a2 = 0.0; double a2_mu = 0.0; double a2_sigma = 0.0; double a2_nu = 0.0;
  double b1 = 0.0; double b1_sigma = 0.0; double b1_nu = 0.0;
  double b2 = 0.0; double b2_sigma = 0.0; double b2_nu = 0.0;
  
  // components of negative hessians
  double a1_mm = 0.0; double a2_mm = 0.0;
  double a1_ss = 0.0; double a2_ss = 0.0; double b1_ss = 0.0; double b2_ss = 0.0;
  double a1_nn = 0.0; double a2_nn = 0.0; double b1_nn = 0.0; double b2_nn = 0.0;
  double a1_ms = 0.0; double a2_ms = 0.0;
  double a1_mn = 0.0; double a2_mn = 0.0;
  double a1_sn = 0.0; double a2_sn = 0.0; double b1_sn = 0.0; double b2_sn = 0.0;
  
  
  // a part --------------------------
  
  if (nu_abs > epsilon[0]){
    double xi = 1/xi_inv;
    
    a1 += xi*(log(z) + 1 - z);
    a1_mu += (z-1)/(mu*sigma2*nu);
    a1_sigma += 2*xi*(z-1 - log(z))/sigma;
    a1_nu += 2*xi*((z-1) - (z+1)*log(z)/2)/nu;
    
    a1_mm += (z*(1 + 1/nu) - 1/nu)/(mu2*sigma2);
    a1_ss += 6*xi*(z - log(z) - 1)/sigma2;
    a1_nn += xi*(log(z)*(z*log(z)-4*z-2)+6*z-6)/nu2;
    a1_ms += 2*(z-1)/(mu*sigma2*sigma*nu);
    a1_mn += xi*(z*(1-log(z)) - 1)/mu;
    a1_sn += 2*xi*(z*(2-log(z))-log(z)-2)/(nu*sigma);
  }
  
  if (nu_abs < epsilon[1]){
    
    // ll
    double log_prod = l_x_mu * l_x_mu / (2 * sigma2);
    a2 -= log_prod; // Part of log-normal log-lik
    for (int k = 3; k <= (nu_abs > 0 ? kappa + 2 : 0); k++) {
      log_prod *= nu * l_x_mu / k;
      a2 -= log_prod;
    }
    
    // SCORES -------------------------------------
    // mu
    log_prod = l_x_mu / (mu*sigma2);
    a2_mu += log_prod;
    for (int k = 2; k <= (nu_abs > 0 ? kappa + 1 : 0); k++) {
      log_prod *= nu * l_x_mu / k;
      a2_mu += log_prod;
    }
    
    // sigma
    log_prod = l_x_mu*l_x_mu / (sigma2*sigma);
    a2_sigma += log_prod;
    for (int k = 3; k <= (nu_abs > 0 ? kappa + 2 : 0); k++) {
      log_prod *= nu * l_x_mu / k;
      a2_sigma += log_prod;
    }
    
    // nu
    log_prod = - l_x_mu*l_x_mu*l_x_mu / (6*sigma2);
    a2_nu += log_prod;
    for (int k = 4; k <= (nu_abs > 0 ? kappa + 3 : 0); k++) {
      log_prod *= nu * l_x_mu / k;
      a2_nu += (k-2) * log_prod;
    }
    
    // NEG. HESSIANS ----------------------------------
    // mm
    log_prod = l_x_mu/(mu2*sigma2);
    a2_mm += 1/(mu2*sigma2) + log_prod;
    for (int k = 1; k <= (nu_abs > 0 ? kappa : 0); k++) {
      log_prod *= nu;
      a2_mm += log_prod;
      log_prod *= l_x_mu/(k+1);
      a2_mm += log_prod;
    }
    
    // ss
    log_prod = 3*l_x_mu*l_x_mu/(sigma2*sigma2);
    a2_ss += log_prod;
    for (int k = 3; k <= (nu_abs > 0 ? kappa+2 : 0); k++) {
      log_prod *= nu * l_x_mu / k;
      a2_ss += log_prod;
    }
    
    // nn
    log_prod = l_x_mu*l_x_mu*l_x_mu*l_x_mu/(24*sigma2);
    a2_nn += 2*log_prod;
    for (int k = 5; k <= (nu_abs > 0 ? kappa + 4 : 0); k++) {
      log_prod *= nu*l_x_mu/k;
      a2_nn += (k-2)*(k-3)*log_prod;
    }
    
    // ms
    log_prod = 2*l_x_mu/(mu*sigma2*sigma);
    a2_ms += log_prod;
    for (int k = 2; k <= (nu_abs > 0 ? kappa + 1 : 0); k++) {
      log_prod *= nu*l_x_mu/k;
      a2_ms += log_prod;
    }
    
    // mn
    log_prod = -l_x_mu*l_x_mu/(2*mu*sigma2);
    a2_mn += log_prod;
    for (int k = 3; k <= (nu_abs > 0 ? kappa + 2 : 0); k++) {
      log_prod *= nu*l_x_mu/k;
      a2_mn += log_prod;
    }
    
    // sn
    log_prod = -l_x_mu*l_x_mu*l_x_mu/(3*sigma*sigma2);
    a2_sn += log_prod;
    for (int k = 3; k <= (nu_abs > 0 ? kappa + 3 : 0); k++) {
      log_prod *= nu*l_x_mu/k;
      a2_sn += (k-2)*log_prod;
    }
  }
  
  double w_a = (nu_abs > epsilon[1]) ? 1 : (nu_abs > epsilon[0] ? (nu_abs-epsilon[0])/(epsilon[1]-epsilon[0]) : 0);
  double a = w_a * a1 + (1-w_a) * a2;
  double a_mu = w_a * a1_mu + (1-w_a) * a2_mu;
  double a_sigma = w_a * a1_sigma + (1-w_a) * a2_sigma;
  double a_nu = w_a * a1_nu + (1-w_a) * a2_nu;
  double a_mm = w_a * a1_mm + (1-w_a) * a2_mm;
  double a_ss = w_a * a1_ss + (1-w_a) * a2_ss;
  double a_nn = w_a * a1_nn + (1-w_a) * a2_nn;
  double a_ms = w_a * a1_ms + (1-w_a) * a2_ms;
  double a_mn = w_a * a1_mn + (1-w_a) * a2_mn;
  double a_sn = w_a * a1_sn + (1-w_a) * a2_sn;
  
  
  // b part --------------------------
  
  if (xi_inv > epsilon[0]){
    double xi = 1/xi_inv;
    double tri_xi = R::trigamma(xi);
    double dig_minus_log = R::digamma(xi) - log(xi);
    b1 += log(nu_abs) - xi*(1-log(xi)) - lgamma(xi);
    b1_sigma += 2*xi*dig_minus_log/sigma;
    b1_nu += 2*xi*dig_minus_log/nu + 1/nu;
    
    b1_ss += xi*(4*xi*tri_xi + 6*dig_minus_log-4)/sigma2;
    b1_nn += xi*(4*xi*tri_xi+6*dig_minus_log-4)/nu2 + 1/nu2;
    b1_sn += 2*xi*(2*xi*tri_xi+2*dig_minus_log-2)/(nu*sigma);
  }
  
  if (xi_inv < epsilon[1]){
    
    // kappa >= 0 ----------------------------------------
    b2 -= log(sigma)+log(2*M_PI)/2;
    b2_sigma -= 1/sigma;
    b2_ss -= 1/sigma2;
    b2_nn += sigma2 * B_even[1];
    
    if(nu_abs > 0 && kappa >= 1){
      double nu_sigma_prod = nu2 * sigma2;
      double nu4_sigma4 = nu_sigma_prod * nu_sigma_prod;
      
      // kappa >= 1 ----------------------------------------
      // nu
      nu_sigma_prod = nu * sigma2;
      nu4_sigma4 = nu_sigma_prod * nu_sigma_prod * nu2;
      b2_nu -= nu_sigma_prod * B_even[1];
      for (int k = 2; k <= (static_cast<double>(kappa) / 4 + .75); k++) {
        nu_sigma_prod *= nu4_sigma4;
        b2_nu -= nu_sigma_prod * B_even[k] / k;
      }
      
      // nn
      nu_sigma_prod = sigma2;
      for (int k = 2; k <= (static_cast<double>(kappa) / 4 + 1.0); k++) {
        nu_sigma_prod *= nu4_sigma4;
        b2_nn += nu_sigma_prod * (4*k-3) * B_even[k]/k;
      }
      
      // sn
      nu_sigma_prod = 2*nu*sigma;
      b2_sn += nu_sigma_prod * B_even[1];
      for (int k = 2; k <= (static_cast<double>(kappa) / 4 + 0.75); k++) {
        nu_sigma_prod *= nu4_sigma4;
        b2_sn += nu_sigma_prod * (2*k-1) * B_even[k] / k;
      }
      
      
      if (kappa >= 2){
        // kappa >= 2 ----------------------------------------
        
        // ll
        nu_sigma_prod = nu2 * sigma2;
        b2 -= nu_sigma_prod * B_even[1] / 2;
        for (int k = 2; k <= (static_cast<double>(kappa) / 4 + 0.5); k++) {
          nu_sigma_prod *= nu4_sigma4;
          b2 -= nu_sigma_prod * B_even[k] / (2*k * (2*k - 1));
        }
        
        // sigma
        nu_sigma_prod = nu2 * sigma;
        nu4_sigma4 = nu_sigma_prod * nu_sigma_prod * sigma2;
        b2_sigma -= nu_sigma_prod * B_even[1];
        for (int k = 2; k <= (static_cast<double>(kappa) / 4 + 0.5); k++) {
          nu_sigma_prod *= nu4_sigma4;
          b2_sigma -= nu_sigma_prod * B_even[k] / k;
        }
        
        // ss
        nu_sigma_prod = nu2;
        b2_ss += nu_sigma_prod * B_even[1];
        for (int k = 2; k <= (static_cast<double>(kappa) / 4 + 0.5); k++) {
          nu_sigma_prod *= nu4_sigma4;
          b2_ss += nu_sigma_prod * (4*k-3) * B_even[k] / k;
        }
      }
    }
  }
  
  double w_b = (xi_inv > epsilon[1]) ? 1 : (xi_inv > epsilon[0] ? (xi_inv-epsilon[0])/(epsilon[1]-epsilon[0]) : 0);
  double b = w_b * b1 + (1-w_b) * b2;
  double b_sigma = w_b * b1_sigma + (1-w_b) * b2_sigma;
  double b_nu = w_b * b1_nu + (1-w_b) * b2_nu;
  double b_ss = w_b * b1_ss + (1-w_b) * b2_ss;
  double b_nn = w_b * b1_nn + (1-w_b) * b2_nn;
  double b_sn = w_b * b1_sn + (1-w_b) * b2_sn;
  
  
  
  // Putting it all together
  // negative log-likelihood
  double nll = log(x) - a - b;
  
  // negative scores
  double mu_s = -a_mu;
  double sigma_s = -a_sigma - b_sigma;
  double nu_s = -a_nu - b_nu;
  
  // negative hessians
  double mm = a_mm;
  double ss = a_ss + b_ss;
  double nn = a_nn + b_nn;
  double ms = a_ms;
  double mn = a_mn;
  double sn = a_sn + b_sn;
  
  
  return NumericVector {nll, mu_s, sigma_s, nu_s, mm, ss, nn, ms, mn, sn};
}


//------------------------------------------------------------------------------
// Cpp wrapper of the negative log-likelihood (w/ gradients and Hessians) ------
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector nllggamma_cpp(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector nu,
                            NumericVector weights, NumericVector epsilon, int kappa) {
  // NumericMatrix h1g_prod, NumericMatrix h2g) {
  
  // Different possibilities for the lengths of x, mu, sigma, nu is taken care of in the R wrapper,
  // where some conditions are also checked.
  
  NumericVector nll(10); // nll, three gradients, and six hessians
  for (int i = 0; i < x.size(); i++)
    nll += weights[i]*nllggamma_single_cpp(x[i], mu[i], sigma[i], nu[i], epsilon, kappa);
  
  return nll;
}
