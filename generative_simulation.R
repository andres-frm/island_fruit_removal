set.seed(23061993)

sapply(c('cmdstanr', 'readxl', 'magrittr', 'dplyr', 'ggplot2', 
         'tidyr', 'tibble', 'forcats', 'geosphere', 
         'MASS', 'rethinking'), library, character.only = T)

source('functions_mod_diagnostics.r')

coords<- read_xlsx('ilhas2025.xlsx', 
                   sheet = 1, col_names = T, 
                   na = 'NA')[, c('lat', 'long', 'island')]

coords <- lapply(split(coords, coords$island), 
                 function(x) x[1, ])

coords <- do.call('rbind', coords)

# shortest distance over the Earth's surface (great-circle distance).
dist_matrix <- distm(coords[, c('long', 'lat')], fun = distHaversine)
dimnames(dist_matrix) <- list(coords$island, 
                              coords$island)

dist_matrix <- dist_matrix/max(dist_matrix)


# ====  indexes and data structure
d <- read_xlsx('ilhas2025.xlsx', sheet = 1, col_names = T, na = 'NA')

sim_data <- 
  d[, c("country", "island", "grid",   
        "realm", "ecoregion", "biome", 
        "island_type")]

sim_data <- 
  lapply(split(sim_data, sim_data$grid), FUN = 
           function(x) {
             x <- x[1:15, ]
             x$plant <- rep(1, 15)
             x
           })

sim_data <- do.call('rbind', sim_data)

for (i in 1:7) sim_data[[i]] <- as.numeric(as.factor(sim_data[[i]]))

sim_data$plant <- 1:nrow(sim_data)
tail(sim_data)

# ===========#

mu_island <- rnorm(nrow(coords), 0, 2) # mu par island

# GP islands 
rho <- 0.25 # Rate of covariance decline between points
curve(exp(-rho * x), from = 0, to = 30, lty = 2)
curve(exp(-rho^2*x^2), add = T)
eta <- 2 # Maximum covariance between consecutive points 
delta <- 0.1 # No variation within sites

quadratic_kernel <- 
  function(m, eta, rho, delta) {
    
    N <- ncol(m)
    mat <- matrix(ncol = N, nrow = N)
    
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        mat[i, j] <- eta * exp(-0.5 * (m[i, j]/rho)^2)
        mat[j, i] <- mat[i, j]
      }
    }
    
    diag(mat) <- eta + delta
    mat
  }

cov_mat <- quadratic_kernel(dist_matrix, eta, rho, delta)
eigen(cov_mat)$values

# spatial correlated residuals
res_cor_res <- 
  mvrnorm(nrow(coords), 
        mu = rep(0, nrow(coords)), 
        Sigma = cov_mat)

dim(res_cor_res)

islands <- 
  res_cor_res + 
  matrix(rep(mu_island, length(mu_island)), 
         nrow = length(mu_island), byrow = T)
# pars 
islands_pars <- islands[1, ] # correlated parameters
plants <- rnorm(nrow(sim_data), 0, 0.01) # pars plant
grid <- rnorm(max(sim_data$grid), 0, 0.015) # pars grid
country <- rnorm(max(sim_data$country), # pars countries
                 0, 0.015)
type_island <- c(1.2, -1, 0) # par type of island

# continuous variables
sim_alt <- 
  do.call('rbind', 
          lapply(1:nrow(sim_data), FUN = 
                   function(x) {
                     
                     i <- sim_data$island_type[x]
                     mu <- c(1.25, -0.5, 0.8)[i]
                     alt <- rnorm(1, mu, 0.25)
                     tibble(altitude = alt, 
                            island_type = i)
                     
                   }))

sim_data$altitude <- sim_alt$altitude

beta_alt <- -0.8

# beta parameters

fruit_removal <- 
  sapply(1:nrow(sim_data), FUN = 
         function(x) {
           p <- sim_data$plant[x]
           g <- sim_data$grid[x]
           i <- sim_data$island[x]
           c <- sim_data$country[x]
           TI <- sim_data$island_type[x]
           alt <- sim_data$altitude[x]
           rbinom(1, 15, # available fruit per bush 
                  inv_logit(alt * beta_alt + 
                            plants[p] + 
                            grid[g] + 
                            islands_pars[i] +
                            country[c] +
                            type_island[TI]))
         }) 
sim_data$fruit_removal <- fruit_removal
plot(density(sim_data$fruit_removal))
plot(sim_data$altitude, 
     sim_data$fruit_removal)

coords$island <- as.factor(coords$island)
levels(coords$island) == colnames(dist_matrix)

dat <- 
  list(N = nrow(sim_data), 
       N_islands = max(sim_data$island),
       N_grid = max(sim_data$grid),
       N_plant = max(sim_data$plant), 
       N_country = max(sim_data$country),
       N_type_island = max(sim_data$island_type),
       altitude = sim_data$altitude,
       fruit_removal = sim_data$fruit_removal,
       type_island = sim_data$island_type,
       country_ID = sim_data$country,
       islands_ID = sim_data$island,
       grid_ID = sim_data$grid, 
       plant_ID = sim_data$plant,
       dist_islands = dist_matrix)


cat(file = 'generative_simulation.stan', 
    '
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
                              K[i, j] = square(eta) * exp(-rho * square(x[i, j]));
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
      real beta_alt;
    
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
      beta_alt ~ normal(0, 1);
      
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
                               beta_alt * altitude +
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
                               beta_alt * altitude +
                               TI[type_island] +
                               country[country_ID] +
                               island[islands_ID] +
                               grid[grid_ID] +
                               plant[plant_ID]
                               ));
    }
    ')


file <- paste0(getwd(), '/generative_simulation.stan')
fit_gen_sim <- cmdstan_model(file, compile = T)

mod_gen_sim <- 
  fit_gen_sim$sample(
    data = dat,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3,
    chains = 4, 
    parallel_chains = 4,
    seed = 06231993
  )

summary_mod_gen <- mod_gen_sim$summary()
mod_diagnostics(mod_gen_sim, summary_mod_gen)

ppcheck <- mod_gen_sim$draws('ppcheck', format = 'matrix')

plot(density(dat$fruit_removal), main = '', 
     xlab = 'Simulated fruit removal')
for (i in 1:200) lines(density(ppcheck[i, ]), lwd = 0.1)
lines(density(dat$fruit_removal), lwd = 2, col = 'red')


posterior_pars <- 
  mod_gen_sim$draws(c('alpha', 
                      'beta_alt',
                      'TI',
                      'country', 
                      'island', 
                      'grid', 
                      'plant'), 
                    format = 'df')


posterior_pars <- 
  lapply(c('alpha', 'beta_alt', 'TI', 'country', 'island', 
         'grid', 'plant'), FUN = 
         function(x) {
           posterior_pars[, grep(x, colnames(posterior_pars))]
         })

names(posterior_pars) <- c('alpha', 'beta_alt', 'TI', 
                           'country', 'island', 
                           'grid', 'plant')

# ==== Parameter recovery =======

# === Random effects 

compare_posterior <- 
  function(post = posterior_pars$island, 
           real = islands_pars, 
           xlab = 'Islands', 
           ylab = 'Posterior mean', 
           main = 'Island parameters') {
    post <- 
      do.call('rbind', 
              lapply(post, FUN = 
                       function(x) {
                         tibble(li = quantile(x, 0), 
                                ls = quantile(x, 1), 
                                mu = mean(x))
                       }))
    post$x <- 1:nrow(post)
    post$real <- real
    
    plot(NULL, xlim = c(0, nrow(post)), ylim = c(-7, 7), 
         ylab = ylab, xlab = xlab, 
         main = main)
    post %$% 
      segments(y0 = li, y1 = ls, x0 = x)
    post %$%
      points(x = x, y = mu)
    post %$%
      points(x = x, y = real, pch = 16, col = 'red')
  }

# countries 

compare_posterior(posterior_pars$country, 
                  country, 
                  xlab = 'Countries', 
                  ylab = 'Posterior mean', 
                  main = 'Country parameters')




# islands 

compare_posterior(posterior_pars$island, 
                  islands_pars, 
                  xlab = 'Islands', 
                  ylab = 'Posterior mean', 
                  main = 'Islands parameters')

# plants 

compare_posterior(posterior_pars$plant, 
                  plants, 
                  xlab = 'Plants', 
                  ylab = 'Posterior mean', 
                  main = 'PLant parameters')

# grids

compare_posterior(posterior_pars$grid, 
                  grid, 
                  xlab = 'Grid', 
                  ylab = 'Posterior mean', 
                  main = 'Grid parameters')


# === Population effects 

compare_posterior(posterior_pars$TI, 
                  type_island, 
                  xlab = 'Type of island', 
                  ylab = 'Posterior mean', 
                  main = 'Type of island parameters')

plot(density(posterior_pars$beta_alt$beta_alt))
abline(v = beta_alt, col = 'red', lwd = 2)
