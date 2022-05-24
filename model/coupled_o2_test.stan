//=======================================================================================

// Coupled benthic-pelagic O2 model
// Running the model with a fixed value for the random walk s.d. 
// The P-I equation has been changed to: P = Pmax*(I/Iopt)*exp(1 - I/Iopt)
// k[n] enters as data that is calculated from lake metabolizer package (rather than being calc in transformed data)
// changed priors for c0 and c1
//=======================================================================================




//=======================================================================================


functions {
  
  // Analytical solution to ODE for coupled benthic-pelagic O2 dynamics
  real x_pred_fn(real x1_n, real x2_n, real x_eq_n, real phi1_n, real phi2_n, real k_n, real D_n,
                    real z_n, int l) {
                      
    // Declare variables
    real theta1;
    real theta2;
    real theta3;
    real theta4;
    real theta5;
    real theta6;
    real pred;
    
    // Define aggregatee parameeteres
    theta1 = sqrt(4*D_n^2 + k_n^2);
    theta2 = 2*D_n + k_n;
    theta3 = exp(theta1/(z_n/2));
    theta4 = exp(theta2/(2*(z_n/2)))*exp(theta1/(2*(z_n/2)));
    theta5 = k_n*(x1_n - x_eq_n) - phi1_n - phi2_n;
    theta6 = k_n*D_n*(x2_n - x_eq_n) - D_n*phi1_n - D_n*phi2_n - k_n*phi2_n;
    
    // Pelagic dynamics
    if(l == 1){
     pred = 1/(2*k_n*theta1*theta4)*(theta1*(theta5+theta5*theta3+2*theta4*(k_n*x_eq_n + phi1_n + phi2_n)) 
             + (theta3 - 1)*(2*D_n*k_n*(x2_n - x_eq_n) + k_n^2*(x_eq_n - x1_n) - 2*D_n*(phi1_n + phi2_n) 
                             + k_n*(phi1_n - phi2_n)));
   } 
   
   // Benthic dynamics
   if(l == 2) {
     pred = 1/(2*k_n*D_n*theta1*theta4)*(theta1*(theta6+theta6*theta3
                                                 +2*theta4*(k_n*D_n*x_eq_n + D_n*phi1_n + D_n*phi2_n
                                                            + k_n*phi2_n)) +
              + (theta3 - 1)*(2*D_n^2*k_n*(x1_n - x_eq_n) + D_n*k_n^2*(x2_n - x_eq_n) 
                              - 2*D_n^2*(phi1_n + phi2_n) - D_n*k_n*(phi1_n + phi2_n) - k_n^2*phi2_n));
   }
   
   // Return
   return pred;
  }
}


//=======================================================================================



data {
  
  // Declare variables
  // indices
  int n_obs; // number of observations
  int n_sites; // number of sites
  int n_days; // number of days
  int n_series; // number of time series
  // int map_sites[n_obs]; // mapping of observations to days
  int map_days[n_obs]; // mapping of observations to days
  int days_per_site[n_sites]; // number of days in each site
  int obs_per_series[n_series]; // number of obervsations in each time series
  int obs_per_day[n_days]; // number of observations in each day
  // actual data
  matrix<lower=0>[n_obs, 2] o2_obs; // observed oxygen [mg m^-3]
  real<lower=0> o2_eq[n_obs]; // equilibrium oxygen [mg m^-3] 
  matrix<lower=0>[n_obs, 2] light; // light [umol-photons m^-2 s^-1]
  matrix<lower=0>[n_obs, 2] temp; // temperature [C]
  real<lower=0> wspeed[n_obs]; // used -2.4 when using log(wspeed), wind speeed [m s^-1]
  //real<lower=0> sch_conv[n_obs]; // Schmidt number conversion (REPLACE WITH K GET RID OF K0K1K2)
  real<lower=0> z[1]; // mixing depth [m]
  real<lower=0> temp_ref; // reference temperature [C]
  real<lower=0> k[n_obs]; // gas exchange constant 0
  //real<lower=0> k1; // gas exchange constant 1
  //real<lower=0> k2; // gas exhange constant 2
  real<lower=0> sig_b0; // AM: putting the 2019 s.d. for beta here
  real<lower=0> sig_r; // AM: putting the 2019 s.d. for beta here (can use 2019 s.d. for rho later on)
  
}


//=======================================================================================



transformed data {
  
  // Declare variables
  real kc; // Boltzmann's constant
  real<lower=0> mu; // mean of DO for entire time series
  real<lower=0> tau; // sd of DO for entire time series
  matrix[n_obs, 2] x_obs; // scaled observed DO
  real x_eq[n_obs]; // scaled equilibrium DO
  real<lower=0> eta; // mean of light for entire time series
  matrix<lower=0>[n_obs, 2] lambda; // scaled light 
  //real<lower=0> k[n_obs]; // scaled gas exchange constant (GET RID OF)
  
  kc = 8.50 * 10^-5;
  
  // Scale data
  mu = mean(o2_obs);
  tau = sd(o2_obs);
  eta = mean(light);
  for (n in 1:n_obs){
    for (j in 1:2){
      x_obs[n,j] = (o2_obs[n,j] - mu)/tau;
      x_eq[n] = (o2_eq[n] - mu)/tau;
      lambda[n,j] = light[n,j]/eta;
    } // j
    // k[n] = (1/z[map_sites[n]])*((k0 + k1*wspeed[n]^k2)/100)*sch_conv[n];
    //k[n] = ((k0 + k1*wspeed[n]^k2)/100)*sch_conv[n]; //if using metabolizer just feed it a vector (GET RID OF)
  } // n
  
}


//=======================================================================================


parameters {
  
  // Declare variables
  //real<lower=0> a[2]; // scaled initial slope
  real<lower=0> opt[2]; // scaled optimal light
  real<lower=0> gamma_1[2]; // scaling of gpp with temperature
  real<lower=0> gamma_2[2]; // scaling of er with temperature
  matrix<lower=0>[n_days,2] b0; // scaled max gpp at temp_ref
  matrix<lower=0>[n_days,2] r; // scaled er at temp_ref
  //real<lower=0> sig_b0; // sd of log_b0 random walk (AM commented out)
  //real<lower=0> sig_r; // sd of log_r random walk (AM commented out)
  real<lower=0> sig_proc; // sd of oxygen state process error
  // real<lower=0> D;
  real c0; // basline of diffusion rate
  real c1; // effect  of wind speed on diffusion rate
  // real zeta; // effective depth of gas exchange
  // matrix[n_obs, 2] x; // scaled inferred oxygen state (AM commented out)
  
}


//=======================================================================================


transformed parameters {
  
  // Declare variables
  matrix[n_obs,2] b; // scaled max gpp at high light
  matrix[n_obs,2] chi; // scaled gpp 
  matrix[n_obs,2] kappa; // scaled er 
  matrix[n_obs,2] phi; // nep [g m^-2 h^-1]
  real D[n_obs]; // oxygen exchange rate between pelagic and benthic zones
  matrix[n_obs,2] x_pred; // predicted oxygen [g m^-3]
  
  // Predicted oxygen 
  for (n in 1:n_obs) {
    // loop over benthic and pelagic zones
    for (j in 1:2) {
      // gpp
      b[n,j] = b0[map_days[n],j]*exp(-gamma_1[j] / ((temp[n,j] + 273.15) * kc));
      //chi[n,j] = b[n,j]*tanh((a[j]/b[n,j])*lambda[n,j]); //doing new Iopt eq
      chi[n,j] = b[n,j]*(lambda[n,j]/opt[j])*exp(1-lambda[n,j]/opt[j]); //new equation
      // er
      kappa[n,j] = r[map_days[n],j]*exp(-gamma_2[j] / ((temp[n,j] + 273.15) * kc));
      // nep
      phi[n,j] = chi[n,j] - kappa[n,j];
    } // j 
    // D[n] = exp(c0 + c1*(wspeed[n]-mean(wspeed))/sd(wspeed))/z[map_sites[n]];
    D[n] = exp(c0 + c1*(wspeed[n]-mean(wspeed))/sd(wspeed));
    // calculate predicted oxygen for pelagic and benthic zones (allow exchange between each)
    x_pred[n,1] = x_pred_fn(x_obs[n,1], x_obs[n,2], x_eq[n], phi[n,1], phi[n,2], k[n], D[n], z[1], 1); // AM: changed x to x_obs
    x_pred[n,2] = x_pred_fn(x_obs[n,1], x_obs[n,2], x_eq[n], phi[n,1], phi[n,2], k[n], D[n], z[1], 2); // AM: changed x to x_obs
  } // n
  
}


//=======================================================================================


model {
  
  // Priors
  // zeta ~ normal(0, 1) T[0, ];
  //sig_b0 ~ normal(0, 1) T[0, ]; // AM: removed here and below
  //sig_r ~ normal(0, 1) T[0, ]; // AM: removed here 
  sig_proc ~ gamma(1.5, 1.5 / 0.5); 
  // D ~ normal(0, 1) T[0, ];
  c0 ~ normal(0, 0.5); //was -3, 3
  c1 ~ normal(0, 0.5); //was -3, 3
  for (j in 1:2) { 
    //a[j] ~ normal(0, 1) T[0, ];
    opt[j] ~ exponential(1 / 0.5);
    gamma_1[j] ~ gamma(4, 4 / 0.5); 
    gamma_2[j] ~ gamma(4, 4 / 0.5);    
  }
  
  // Time-varuing rates
  {
      int pos = 1;
      // loop over sites
      for (s in 1:n_sites) {
        for (j in 1:2) {
          // inital value for time varying-parameters
          b0[pos,j] ~ exponential(1 / 0.5);
          r[pos,j] ~ exponential(1 / 0.5);
          // random walk for time varying-parameters
          for (d in (pos+1):(pos+days_per_site[s]-1)){
            b0[d,j] ~ normal(b0[d-1,j], sig_b0) T[0, ];
            r[d,j] ~ normal(r[d-1,j], sig_r) T[0, ];
          } // d
        } // j
        pos = pos + days_per_site[s];
      } // s
    }
    
  // State process
  {
    int pos = 1; 
    // loop over benthic and pelagic zones
      // loop over contiguous time series
      for (t in 1:n_series) {
        for (j in 1:2) {
        // initial value
        //x_obs[pos,j] ~ normal(x_obs[pos,j], 0.01); // AM: assuming the initial value was what was observed, no longer need this (also changed x to x_obs)
        // project oxygen at t based on predicted value at t-1 plus process error
        x_obs[(pos+1):(pos+obs_per_series[t]-1),j] // AM: changed x to x_obs
          ~ normal(x_pred[pos:(pos+obs_per_series[t]-2),j], sig_proc); // AM: changed x to x_obs; serving as likelihood for model now rather than 241
        }
        pos = pos + obs_per_series[t];
    }
  }
  
  // Observation process
  // AM: comment this section out; currently penalizing deviation from x and x_obs
  //for (j in 1:2) {
  //  x_obs[,j] ~ normal(x[,j], 0.01); 
  //}
  
}


//=======================================================================================


generated quantities {
  
  // Declare variables
  //real o2[n_obs, 2]; // inferred oxygen state [mg m^-3] omitted 31mar21
  real o2_pred[n_obs, 2]; // predicted oxygen [mg m^-3]
  real beta0[n_days, 2]; // max gpp at temp_ref [mg m^-2 h^-1]
  //real alpha[2]; // initial slope [mg-O2 s umol-photons-1 m-1 h-1]
  real omega[2]; // optimal light [umol-photons-1 m-2 h-1]
  real rho[n_days, 2]; // er at temp_ref [mg m^-2 h^-1]
  real beta[n_obs, 2]; // max gpp at high light [mg m^-2 h^-1]
  real gpp[n_obs, 2]; // gpp [mg m^-2 h^-1]
  real er[n_obs, 2]; // er [mg m^-2 h^-1]
  real nep[n_obs, 2]; // nep [mg m^-2 h^-1]
  real GPP[n_days, 2]; // total daily flux [g m^-2 d^-1]
  real ER[n_days, 2]; // total daily flux [g m^-2 d^-1]
  real NEP[n_days, 2]; // total daily flux [g m^-2 d^-1]
  // real AIR[n_days]; // total daily flux [g m^-2 d^-1]
  
  // Back-tranformed scaled variables [g m^-2 d^-1]
  for(j in 1:2){
   //alpha[j] = tau*a[j]/eta;
   omega[j] = eta*opt[j];
  }
  for(n in 1:n_obs){
    for(j in 1:2){
      //o2[n,j] = tau*x_obs[n,j] + mu; // AM: changed x to x_obs (omitted on 31mar21)
      o2_pred[n,j] = tau*x_pred[n,j] + mu;
      beta[n,j] = tau*b[n,j];
      gpp[n,j] = tau*chi[n,j];
      er[n,j] = tau*kappa[n,j];
      nep[n,j] = tau*phi[n,j];
    }
  }
  {
    int pos = 1;
    for (d in 1:n_days){
      // AIR[d] = 24*mean(air[pos:(pos + obs_per_day[d] - 1)]);
      for(j in 1:2){
        beta0[d,j] = tau*b0[d,j];
        rho[d,j] = tau*r[d,j];
        GPP[d,j] = 24*mean(gpp[pos:(pos + obs_per_day[d] - 1),j]);
        ER[d,j] = 24*mean(er[pos:(pos + obs_per_day[d] - 1),j]);
        NEP[d,j] = 24*mean(nep[pos:(pos + obs_per_day[d] - 1),j]);
      } // j
      pos = pos + obs_per_day[d]; // advance position counter
  } // d
  }
  
}


//=======================================================================================
