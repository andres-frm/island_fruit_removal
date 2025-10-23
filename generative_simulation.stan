
    functions{
    
      matrix GP_quadratic(matrix x, 
                          real eta, 
                          real rho, 
                          real delta) {
                          
                          int N = dims(x)[1];
                          matrix[N, N] K;
                          matrix[N, N] L_K;
                          
                          for (i in 1:(N-1)) {
                            K[i, i] = eta + delta;
                            for (j in (i+1):N) {
                              K[i, j] = eta * exp(-rho * square(x[i, j]));
                              K[j, i] = K[i, j];
                            }
                          }
                          
                          K[N, N] = eta + delta;
                          L_K = cholesky_decompose(K);
                          return L_K;
                          }
    }
    
    data {
      int N;
      int N_islands;
      int N_country;
      int N_plant;
      int N_grid;
      int N_type_island;
      // response
      array[N] int fruit_removal;
      // propulation effects
      array[N] int type_island;
      vector[N] altitude;
      vector[N] island_isolation;
      vector[N] native_cover;
      vector[N] bush_cover;
      // group level effects
      array[N] int islands_ID;
      array[N] int country_ID;
      array[N] int grid_ID;
      array[N] int plant_ID;
      // matrix of island distances (std)
      matrix[N_islands, N_islands] dist_islands;
    }
    
    parameters {
      
      // intercept
      real alpha;
      
      // population effects
      vector[N_type_island] TI;
      // real beta_alt;
      real beta_iso;
      real beta_NV;
      // real beta_bush;
    
      // group level effects
      // GP island
      vector[N_islands] z_islands;
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
      
    }
    
    transformed parameters {
      
      // Population effects
      
      
      // group level effects
      // GP islands
      vector[N_islands] island;
      matrix[N_islands, N_islands] L_K_islands;
      L_K_islands = GP_quadratic(dist_islands,
                                 eta, 
                                 rho, 
                                 0.001);
      island = L_K_islands * z_islands;
      
      // country
      vector[N_country] country;
      country = mu_country + z_country * sigma_country;
      
      // grid
      vector[N_grid] grid;
      grid = mu_grid + z_grid * sigma_grid;
      
      // plant
      vector[N_plant] plant;
      plant = mu_plant + z_plant * sigma_plant;
    }
    
    model {
      
      // intercept and dispersion
      alpha ~ normal(0, 1);
      
      // Population effects
      TI ~ normal(0, 1);
      //beta_alt ~ normal(0, 1);
      beta_iso ~ normal(0, 1);
      beta_NV ~ normal(0, 1);
      // beta_bush ~ normal(0, 1);
      
      // GP islands
      eta ~ exponential(4);
      rho ~ exponential(1);
      z_islands ~ normal(0, 1);
      
      // country
      z_country ~ normal(0, 1);
      mu_country ~ normal(0, 0.5);
      sigma_country ~ exponential(1);
      
      // grid
      z_grid ~ normal(0, 1);
      mu_grid ~ normal(0, 0.5);
      sigma_grid ~ exponential(1);
      
      // plant
      z_plant ~ normal(0, 1);
      mu_plant ~ normal(0, 0.5);
      sigma_plant ~ exponential(1);
      
      fruit_removal ~ binomial(15, 
                               inv_logit(
                               alpha +
                               // beta_alt * altitude +
                               beta_iso * island_isolation +
                               beta_NV * native_cover +
                               // beta_bush * bush_cover +
                               TI[type_island] +
                               country[country_ID] +
                               island[islands_ID] +
                               grid[grid_ID] +
                               plant[plant_ID]
                               ));
    }
    
    generated quantities {
      array[N] int ppcheck;
      
      ppcheck = binomial_rng(15, 
                             inv_logit(
                               alpha +
                               //beta_alt * altitude +
                               beta_iso * island_isolation +
                               beta_NV * native_cover +
                               // beta_bush * bush_cover +
                               TI[type_island] +
                               country[country_ID] +
                               island[islands_ID] +
                               grid[grid_ID] +
                               plant[plant_ID]
                               ));
    }
    