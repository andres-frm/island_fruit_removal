functions{
    
    // function for conducting imputation of NA values
    vector merge_missing(array[] int miss_indx, vector x_obs, vector x_miss) {
      int N = dims(x_obs)[1];
      int N_miss = dims(x_miss)[1];
      vector[N] merge;
      merge = x_obs;
      
      for (i in 1:N_miss) {
        merge[miss_indx[i]] = x_miss[i];
      }
      return merge;
    }
    
    // Gaussian processes function
    matrix GP_quadratic(matrix x, 
                          real eta, 
                          real rho, 
                          real delta) {
                          
                          int N = dims(x)[1];
                          matrix[N, N] K;
                          matrix[N, N] L_K;
                          
                          for (i in 1:(N-1)) {
                            
                            K[i, i] = eta + delta; // small number to the diagonal
                            
                            for (j in (i+1):N) {
                              // quadratic kernell
                              K[i, j] = eta * exp(-rho * square(x[i, j]));
                              // filling lower part of the matrix 
                              K[j, i] = K[i, j];
                            }
                          }
                          
                          K[N, N] = eta + delta; // small number for stability
                          // Cholesky decomposition
                          L_K = cholesky_decompose(K);
                          return L_K;
                          }
      
    }
    
    data {
      int N;
      int N_island;
      int N_country;
      int N_plant;
      int N_grid;
      int N_type_island;
      int N_realm;
      int N_ecoregion;
      int N_biome;
      int N_plant_invasive_rank;
      int N_na_bush;
      // response
      array[N] int total_remotion;
      array[N] int Arthropod;
      array[N] int Bird;
      array[N] int Lizard;
      array[N] int Mammal_non_rodent;
      array[N] int rodent_all; 
      array[N] int dispersion;
      array[N] int predation;
      // propulation effects
      vector[N] lat;
      array[N] int plant_invasive_rank;
      array[N] int island_type;
      vector[N] human_pop;
      vector[N] human_footprint;
      vector[N] distance_mainland;
      vector[N] island_size;
      vector[N] altitude_m;
      vector[N] isolation;
      vector[N] temperature;
      vector[N] native_vegetation;
      vector[N] bush_cover;
      // group level effects
      array[N] int country;
      array[N] int island;
      array[N] int realm;
      array[N] int ecoregion;
      array[N] int biome;
      array[N] int grid;
      array[N] int plant_ID;
      // Na imputation
      array[N_na_bush] int na_bush;
      // matrix of island distances (std)
      matrix[N_island, N_island] dist_island;
    }
    
    parameters {
      
      // imputed bush cover
      // vector[N_na_bush] bush_imputed;
      // real bush_mu;
      // real<lower = 0> bush_sigma;
      
      // intercept
      real alpha;
      
      // population effects
      vector[N_type_island] TI; 
      vector[N_realm] p_realm;
      vector[N_plant_invasive_rank] inv_rank;
      //real beta_lat; // main effect
      // real beta_H_pop;
      real beta_H_foot;
      // real beta_I_mainland;
      //real beta_I_size;
      real beta_I_alt;
      real beta_I_isolation;
      real beta_temp;
      real beta_NV;
      // real beta_bush;
    
      // group level effects
      // GP island
      vector[N_island] z_islands;
      real<lower = 0> eta;
      real<lower = 0> rho;
      
      // country
      vector[N_country] z_country;
      real mu_country;
      real<lower = 0> sigma_country;
      
      // grid
      vector[N_grid] z_grid;
      real mu_grid;
      real<lower = 0> sigma_grid;
      
      // plant
      vector[N_plant] z_plant;
      real mu_plant;
      real<lower = 0> sigma_plant;
      
      // ecoregion
      vector[N_ecoregion] z_ecoR;
      real mu_ecoR;
      real<lower = 0> sigma_ecoR;
      
      // biome
      vector[N_biome] z_biome;
      real mu_biome;
      real<lower = 0> sigma_biome;
      
    }
    
    transformed parameters {
      
      // imputed variable
      // vector[N] bush_merge;
      // bush_merge = merge_missing(na_bush, 
      //                            to_vector(bush_cover), 
      //                            bush_imputed);
      
      // Population effects
      
      // group level effects
      // GP islands
      vector[N_island] p_island;
      matrix[N_island, N_island] L_K_islands;
      L_K_islands = GP_quadratic(dist_island,
                                 eta,
                                 rho,
                                 1e-6);
      p_island = L_K_islands * z_islands;
      
      // country
      vector[N_country] p_country;
      p_country = mu_country + z_country * sigma_country;
      
      // grid
      vector[N_grid] p_grid;
      p_grid = mu_grid + z_grid * sigma_grid;
      
      // plant
      vector[N_plant] p_plant;
      p_plant = mu_plant + z_plant * sigma_plant;
      
      // ecoregion
      vector[N_ecoregion] p_ecoR;
      p_ecoR = mu_ecoR + z_ecoR * sigma_ecoR;
      
      // biome
      vector[N_biome] p_biome;
      p_biome = mu_biome + z_biome * sigma_biome;
      
    }
    
    model {
      
      // imputation
      // bush_mu ~ normal(0, 1);
      // bush_sigma ~ exponential(1);
      // bush_merge ~ normal(bush_mu, bush_sigma);
      
      // intercept and dispersion
      alpha ~ normal(0, 1);
      
      // Population effects
      TI ~ normal(0, 0.5);
      p_realm ~ normal(0, 0.5);
      inv_rank ~ normal(0, 0.25);
      //beta_lat ~ normal(0, 1); 
      // beta_H_pop ~ normal(0, 1);
      beta_H_foot ~ normal(0, 0.25);
      //beta_I_mainland ~ normal(0, 1);
      // beta_I_size ~ normal(0, 1);
      beta_I_alt ~ normal(0, 0.25);
      beta_I_isolation ~ normal(0, 0.25);
      beta_temp ~ normal(0, 0.25);
      beta_NV ~ normal(0, 0.25);
      //beta_bush ~ normal(0, 1);
      
      // GP islands
      eta ~ exponential(3);
      rho ~ exponential(1);
      z_islands ~ normal(0, 1);
      
      // country
      z_country ~ normal(0, 0.5);
      mu_country ~ normal(0, 0.25);
      sigma_country ~ exponential(1.5);
      
      // grid
      z_grid ~ normal(0, 0.5);
      mu_grid ~ normal(0, 0.25);
      sigma_grid ~ exponential(1.5);
      
      // plant
      z_plant ~ normal(0, 0.5);
      mu_plant ~ normal(0, 0.25);
      sigma_plant ~ exponential(1.5);
      
      // Ecoregion
      z_ecoR ~ normal(0, 0.5);
      mu_ecoR ~ normal(0, 0.25);
      sigma_ecoR ~ exponential(1.5);
      
      // Biome
      z_biome ~ normal(0, 0.5);
      mu_biome ~ normal(0, 0.25);
      sigma_biome ~ exponential(1.5);
      
      Bird ~ binomial(1, 
                               inv_logit(
                                 alpha +
                               //beta_lat * lat +
                               // beta_H_pop * human_pop +
                               beta_H_foot * human_footprint +
                               //beta_I_mainland * +
                               // beta_I_size * island_size +
                               beta_I_alt * altitude_m +
                               beta_I_isolation * isolation +
                               beta_temp * temperature +
                               beta_NV * native_vegetation + // main effect
                               //beta_bush * +
                               inv_rank[plant_invasive_rank] +
                               TI[island_type] +
                               p_island[island] +
                               p_country[country] +
                               p_grid[grid] +
                               p_plant[plant_ID] +
                               p_realm[realm] +
                               p_ecoR[ecoregion] +
                               p_biome[biome]));
    }
    
    generated quantities {
      array[N] int ppcheck;
      
      ppcheck = binomial_rng(1, 
                             inv_logit(
                               alpha +
                               //beta_lat * lat +
                               // beta_H_pop * human_pop +
                               beta_H_foot * human_footprint +
                               //beta_I_mainland * +
                               // beta_I_size * island_size +
                               beta_I_alt * altitude_m +
                               beta_I_isolation * isolation +
                               beta_temp * temperature +
                               beta_NV * native_vegetation + // main effect
                               //beta_bush * +
                               inv_rank[plant_invasive_rank] +
                               TI[island_type] +
                               p_island[island] +
                               p_country[country] +
                               p_grid[grid] +
                               p_plant[plant_ID] +
                               p_realm[realm] +
                               p_ecoR[ecoregion] +
                               p_biome[biome]
                               ));
    }