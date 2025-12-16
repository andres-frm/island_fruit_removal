# ======== Loading packages =======

sapply(c('cmdstanr', 'readxl', 'magrittr', 'dplyr', 'ggplot2', 
         'tidyr', 'tibble', 'forcats', 'rethinking', 
         'cowplot', 'bayesplot', 'patchwork', 'ggh4x'), 
       library, character.only = T)

source('functions_mod_diagnostics.r')

# ======== Data preparation =======

data <- readRDS('data.rds')

codes <- 
  lapply(data$codes, function(x) x[order(x$code), ])

codes$real$island <- c('Afrotropics', 'Australasia', 'Indomalaya', 
                       'Nearctic', 'Neotropics', 'Oceania', 'Palearctic')

d <- data$data
dis_matrix <- data$dist_islands

# Na factors

sapply(d, function(x) mean(is.na(x)))

d[, c('lat', "human_pop", "human_footprint", "distance_mainland", 
      "island_size", "isolation", "altitude_m", "temperature", 
      "native_vegetation")]

dat <- d

for (i in c('lat', "human_pop", "human_footprint", "distance_mainland", 
            "island_size", "isolation", "altitude_m", "temperature", 
            "native_vegetation")) {
  dat[[i]] <- as.vector(scale(dat[[i]], center = T, scale = T))
}

# bush cover 10.5%
dat[which(is.na(dat$bush_cover)), c("realm", "country2", 
                                    "island", "grid")] |> unique()
# plant invasion rank 2.6%
dat$bush_cover[which(!is.na(d$bush_cover))] <- 
  as.vector(scale(dat$bush_cover[which(!is.na(d$bush_cover))]))

indx_na_bush <- which(is.na(d$bush_cover))
dat$bush_cover[indx_na_bush] <- 0

# NA of invasive rank in Cozumel island will be treated as a 

dat[which(is.na(dat$plant_invasive_rank)), c("realm", "country2", 
                                             "island", "grid")] |> unique()

dat$plant_invasive_rank[which(!is.na(dat$plant_invasive_rank))] <- 
  dat$plant_invasive_rank[which(!is.na(dat$plant_invasive_rank))] + 1

unique(dat$plant_invasive_rank)

dat$plant_invasive_rank[which(is.na(dat$plant_invasive_rank))] <- 4

mean(dat$Arthropod > 0)
mean(dat$Rodent > 0)
dat[dat$Unknown == 15, ] |> print(n = 45)
predation <- dat$Arthropod + dat$Rodent
predation[which(dat$Unknown == 15)] <- 15
dat$predation <- predation
mean(dat$Lizard > 0)
mean(dat$Bird > 0)
dat$dispersion <- dat$Lizard + dat$Bird
dat$rodent_all <- dat$Rodent
dat$rodent_all[which(dat$Unknown == 15)] <- 15

dat <- lapply(dat, FUN = function(x) if (is.factor(x)) as.numeric(x) else x)

names(dat)[1] <- 'country'

dat <- dat[-c(grep('season', names(dat)), 
              grep('correct_coord', names(dat)))]

dat$dist_island <- data$dist_islands
dat$N <- length(dat$country)
dat$N_type_island <- max(dat$island_type)
dat$N_country <- max(dat$country)
dat$N_island <- max(dat$island)
dat$N_realm <- max(dat$realm)
dat$N_ecoregion <- max(dat$ecoregion)
dat$N_biome <- max(dat$biome)
dat$N_grid <- max(dat$grid)
dat$N_plant <- max(dat$plant_ID)
dat$N_plant_invasive_rank <- max(dat$plant_invasive_rank)
dat$N_na_bush <- length(indx_na_bush)
dat$na_bush <- indx_na_bush
dat$Lizard <- ifelse(dat$Lizard > 0, 1, 0)
dat$Bird <- ifelse(dat$Bird > 0, 1, 0)

# ======== Custom functions  =======

average_effects <- function(n_levels, 
                            posterior,
                            x_var1,
                            x_var2, 
                            par1, 
                            par2) {
  
  y <- lapply(1:n_levels, FUN =
                function(i) {
                  
                  indx_xvar <- which(dat[[x_var1]] == i)
                  indx_xvar2 <- unique(dat[[x_var2]][indx_xvar])
                  indx_count <- unique(dat$country[indx_xvar])
                  indx_is <- unique(dat$island[indx_xvar])
                  indx_grid <- unique(dat$grid[indx_xvar])
                  indx_plant <- unique(dat$plant_ID[indx_xvar])
                  indx_ecoR <- unique(dat$ecoregion[indx_xvar])
                  indx_biome <- unique(dat$biome[indx_xvar])
                  
                  est <-
                    with(posterior,
                         {
                           inv_logit(alpha[[1]] +
                                       posterior[[par1]][[i]] +
                                       apply(posterior[[par2]][, indx_xvar2],
                                             1, mean) +
                                       apply(p_island[, indx_is], 1, mean) +
                                       apply(p_country[, indx_count], 1, mean) +
                                       apply(p_grid[, indx_grid], 1, mean) +
                                       apply(p_plant[, indx_plant], 1, mean) +
                                       apply(p_ecoR[, indx_ecoR], 1, mean) +
                                       apply(p_biome[, indx_biome], 1, mean))
                         })
                  
                  tibble(y = est,
                         code = i,
                         x = x_var1)
                  
                })
  
  do.call('rbind', y)
  
}

cond_effects <- function(posterior, 
                         x_bar, 
                         slope, 
                         type = 'averaging', 
                         n = 100,
                         causal_effect = F, 
                         intervention1,
                         intervention2
                         ) {
  
  x_seq <- seq(min(x_bar), max(x_bar), length.out = n)
  
  type_is <- apply(posterior$TI, 1, mean)
  island <- apply(posterior$p_island, 1, mean)
  country <- apply(posterior$p_country, 1, mean)
  grid <- apply(posterior$p_grid, 1, mean)
  plant <- apply(posterior$p_plant, 1, mean)
  realm <- apply(posterior$p_realm, 1, mean)
  eco <- apply(posterior$p_ecoR, 1, mean)
  biome <- apply(posterior$p_biome, 1, mean)
  
  if (causal_effect) {
    
    causal_effect <- 
    lapply(c(intervention1, intervention2), FUN = 
             function(x) {
               est <- 
                 with(posterior, 
                      {
                        inv_logit(alpha[[1]] +
                                    beta[[slope]] * x +
                                    type_is +
                                    island +
                                    country +
                                    grid +
                                    plant +
                                    realm +
                                    eco +
                                    biome)
                      })
               
               tibble(x = est)
             })
    
    causal_effect <- do.call('cbind', causal_effect)
    colnames(causal_effect) <- c('intervention_1', 'intervention_2')
    causal_effect$contrast <-
      causal_effect[[2]] - causal_effect[[1]]
    
    colnames(causal_effect)[3] <- paste(intervention1, 'SD', 
                                        '\nvs.', 
                                        intervention2, 'SD')

    as_tibble(causal_effect)
    
  } else {
    
    y <- lapply(seq_along(x_seq), FUN = 
                  function(i) {
                    
                    x <- x_seq[i]
                    
                    est <- 
                      with(posterior, 
                           {
                             inv_logit(alpha[[1]] +
                                         beta[[slope]] * x +
                                         type_is +
                                         island +
                                         country +
                                         grid +
                                         plant +
                                         realm +
                                         eco +
                                         biome)
                           })
                    
                    if (type == 'averaging') {
                      tibble(x = x,
                             y = median(est),
                             li = quantile(est, 0.025),
                             ls = quantile(est, 0.975))
                      
                    } else if (type == 'random') {
                      
                      set.seed(23061993)
                      est <- sample(est, 100, replace = F)
                      pred <- rbinom(length(est), 15, est)
                      tibble(x = x,
                             y = median(pred),
                             li = quantile(pred, 0.025),
                             ls = quantile(pred, 0.975),
                             indx = i,
                             type = 'Predicted')
                      
                    } else {
                      
                      tibble(x = x,
                             y = est,
                             indx = i)
                    }
                    
                  })
    
    do.call('rbind', y)
    
  }
  
  
}



slope_pars <- 
  function(posterior, 
           slope) {
    
    env <- ls(globalenv())

    post1 <- env[grep(posterior, env)]
    
    post <-
      lapply(post1, FUN =
             function(x) {
               p <- get(x)
               beta <- p$beta[[slope]]

               tibble(mu = median(beta),
                      li = quantile(beta, 0.025),
                      ls = quantile(beta, 0.975),
                      `P(beta > 0)` = mean(beta > 0),
                      `P(beta < 0)` = mean(beta < 0),
                      `Fruit consumption` = x)
             })

    post <- do.call('rbind', post)
    post$`Fruit consumption` <- 
      c('Bird frugivory', 'Seed dispersal', 
        'Lizard frugivory', 
        'Seed predation',
        'Total frugivory')
    post
  }




#########
# =============== ***Biogeographical effects*** ===========
########

# =============== *Latitude* ==========

mean(dat$total_remotion)

# =============== Overall frugivory  ======================

file <- paste0(getwd(), '/mod_latitud_total.stan')
fit_latitude_tot <- cmdstan_model(file, compile = T)

mod_latitude_tot <- 
  fit_latitude_tot$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_latitude_tot$save_object('mod_latitude_tot.rds')

mod_latitude_tot <- readRDS('mod_latitude_tot.rds')


mcmc_trace(mod_latitude_tot$draws(c('alpha', 
                                     'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     # 'beta_I_isolation',
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
                                    'TI', 'p_realm', 
                                    'eta', 'rho')))

sum_latitude_tot <- mod_latitude_tot$summary()
mod_diagnostics(mod_latitude_tot, sum_latitude_tot)
ppcheck_latitude_tot <- mod_latitude_tot$draws('ppcheck', format = 'matrix')

plot(density(dat$total_remotion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.1))
for (i in 1:200) lines(density(ppcheck_latitude_tot[i, ]), lwd = 0.1)
lines(density(dat$total_remotion), lwd = 2, col = 'red')


post_latitude_tot <- 
  mod_latitude_tot$draws(c('alpha', 
                       'beta_lat', 
                       # 'beta_H_pop', 
                       # 'beta_H_foot',
                       # 'beta_I_mainland', 
                       # 'beta_I_size',
                       # 'beta_I_alt', 
                       # 'beta_I_isolation',
                       # 'beta_temp', 
                       # 'beta_NV',
                       # 'beta_bush', 
                       # 'inv_rank', 
                       'TI', 'p_island', 
                       'p_country', 'p_grid', 
                       'p_plant', 'p_realm', 
                       'p_ecoR', 'p_biome'), 
                     format = 'df')

post_latitude_tot <- 
  lapply(c('alpha', 
           'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           # 'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_latitude_tot[, grep(x, colnames(post_latitude_tot))]
           })

names(post_latitude_tot) <- c('alpha', 
                          'beta', 
                          # 'beta_H_pop', 
                          # 'beta_H_foot',
                          # 'beta_I_mainland', 
                          # 'beta_I_size',
                          # 'beta_I_alt', 
                          # 'beta_I_isolation',
                          # 'beta_temp', 
                          # 'beta_NV',
                          # 'beta_bush', 
                          # 'inv_rank', 
                          'TI', 'p_island', 
                          'p_country', 'p_grid', 
                          'p_plant', 'p_realm', 
                          'p_ecoR', 'p_biome')

post_latitude_tot

avg_prov <- 
  with(post_latitude_tot, 
     {
       inv_logit(alpha[[1]] +
                   apply(TI, 1, median) +
                   apply(p_realm, 1, median) +
                   apply(p_island, 1, median) +
                   apply(p_country, 1, median) +
                   apply(p_grid, 1, median) +
                   apply(p_plant, 1, median) +
                   apply(p_ecoR, 1, median) +
                   apply(p_biome, 1, median))
     })

plot(density(avg_prov))
quantile(avg_prov, c(0.025, 0.5, 0.975))


# =============== Birds  ======================

file <- paste0(getwd(), '/mod_latitud_bird.stan')
fit_latitude_bird <- cmdstan_model(file, compile = T)

mod_latitude_bird <- 
  fit_latitude_bird$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 4e3,
    thin = 4, 
    seed = 23061993
  )

mod_latitude_bird$save_object('mod_latitude_bird.rds')

mod_latitude_bird <- readRDS('mod_latitude_bird.rds')

mcmc_trace(mod_latitude_bird$draws(c('alpha', 
                                     'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     # 'beta_I_isolation',
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
                                     'TI', 'p_realm', 
                                     'eta', 'rho')))


sum_latitude_bird <- mod_latitude_bird$summary()
mod_diagnostics(mod_latitude_bird, sum_latitude_bird)
ppcheck_latitude_bird <- mod_latitude_bird$draws('ppcheck', format = 'matrix')

plot(density(dat$Bird), main = '', 
     xlab = 'Birds fruits removal', ylim = c(0, 3))
for (i in 1:200) lines(density(ppcheck_latitude_bird[i, ], lwd = 0.1))
lines(density(dat$Bird), lwd = 2, col = 'red')


post_latitude_bird <- 
  mod_latitude_bird$draws(c('alpha', 
                           'beta_lat', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_size',
                           # 'beta_I_alt', 
                           # 'beta_I_isolation',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome'), 
                         format = 'df')

post_latitude_bird <- 
  lapply(c('alpha', 
           'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           # 'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_latitude_bird[, grep(x, colnames(post_latitude_bird))]
           })

names(post_latitude_bird) <- c('alpha', 
                              'beta', 
                              # 'beta_H_pop', 
                              # 'beta_H_foot',
                              # 'beta_I_mainland', 
                              # 'beta_I_size',
                              # 'beta_I_alt', 
                              # 'beta_I_isolation',
                              # 'beta_temp', 
                              # 'beta_NV',
                              # 'beta_bush', 
                              # 'inv_rank', 
                              'TI', 'p_island', 
                              'p_country', 'p_grid', 
                              'p_plant', 'p_realm', 
                              'p_ecoR', 'p_biome')



# =============== Lizards  ======================

file <- paste0(getwd(), '/mod_latitud_lizard.stan')
fit_latitude_lizard <- cmdstan_model(file, compile = T)

mod_latitude_lizard <- 
  fit_latitude_lizard$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

mod_latitude_lizard$save_object('mod_latitude_lizard.rds')

mod_latitude_lizard <- readRDS('mod_latitude_lizard.rds')

mcmc_trace(mod_latitude_lizard$draws(c('alpha', 
                                     'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     # 'beta_I_isolation',
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
                                     'TI', 'p_realm', 
                                     'eta', 'rho')))

sum_latitude_lizard <- mod_latitude_lizard$summary()
mod_diagnostics(mod_latitude_lizard, sum_latitude_lizard)
ppcheck_latitude_lizard <- mod_latitude_lizard$draws('ppcheck', format = 'matrix')

plot(density(dat$Lizard), main = '', 
     xlab = 'Lizard fruits removal', ylim = c(0, 7))
for (i in 1:200) lines(density(ppcheck_latitude_lizard[i, ], lwd = 0.1))
lines(density(dat$Lizard), lwd = 2, col = 'red')


post_latitude_lizard <- 
  mod_latitude_lizard$draws(c('alpha', 
                            'beta_lat', 
                            # 'beta_H_pop', 
                            # 'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_size',
                            # 'beta_I_alt', 
                            # 'beta_I_isolation',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_latitude_lizard <- 
  lapply(c('alpha', 
           'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           # 'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_latitude_lizard[, grep(x, colnames(post_latitude_lizard))]
           })

names(post_latitude_lizard) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_size',
                               # 'beta_I_alt', 
                               # 'beta_I_isolation',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')



# =============== Fruit dispersion  ======================

file <- paste0(getwd(), '/mod_latitud_dispersion.stan')
fit_latitude_disp <- cmdstan_model(file, compile = T)

mod_latitude_disp <- 
  fit_latitude_disp$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_latitude_disp$save_object('mod_latitude_disp.rds')

mod_latitude_disp <- readRDS('mod_latitude_disp.rds')

mcmc_trace(mod_latitude_disp$draws(c('alpha', 
                                     'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     # 'beta_I_isolation',
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
                                     'TI', 'p_realm', 
                                     'eta', 'rho')))

sum_latitude_disp <- mod_latitude_disp$summary()
mod_diagnostics(mod_latitude_disp, sum_latitude_disp)
ppcheck_latitude_disp <- mod_latitude_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_latitude_disp[i, ], lwd = 0.1))
lines(density(dat$dispersion), lwd = 2, col = 'red')

post_latitude_disp <- 
  mod_latitude_disp$draws(c('alpha', 
                           'beta_lat', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_size',
                           # 'beta_I_alt', 
                           # 'beta_I_isolation',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome'), 
                         format = 'df')

post_latitude_disp <- 
  lapply(c('alpha', 
           'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           # 'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_latitude_disp[, grep(x, colnames(post_latitude_disp))]
           })

names(post_latitude_disp) <- c('alpha', 
                              'beta', 
                              # 'beta_H_pop', 
                              # 'beta_H_foot',
                              # 'beta_I_mainland', 
                              # 'beta_I_size',
                              # 'beta_I_alt', 
                              # 'beta_I_isolation',
                              # 'beta_temp', 
                              # 'beta_NV',
                              # 'beta_bush', 
                              # 'inv_rank', 
                              'TI', 'p_island', 
                              'p_country', 'p_grid', 
                              'p_plant', 'p_realm', 
                              'p_ecoR', 'p_biome')



# =============== Fruit predation  ======================

file <- paste0(getwd(), '/mod_latitud_predation.stan')
fit_latitude_pred <- cmdstan_model(file, compile = T)

mod_latitude_pred <- 
  fit_latitude_pred$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

mod_latitude_pred$save_object('mod_latitude_pred.rds')

mod_latitude_pred <- readRDS('mod_latitude_pred.rds')

mcmc_trace(mod_latitude_pred$draws(c('alpha', 
                                     'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     # 'beta_I_isolation',
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
                                     'TI', 'p_realm', 
                                     'eta', 'rho')))

sum_latitude_pred <- mod_latitude_pred$summary()
mod_diagnostics(mod_latitude_pred, sum_latitude_pred)

ppcheck_latitude_pred <- mod_latitude_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_latitude_pred[i, ], lwd = 0.1))
lines(density(dat$predation), lwd = 2, col = 'red')

post_latitude_pred <- 
  mod_latitude_pred$draws(c('alpha', 
                            'beta_lat', 
                            # 'beta_H_pop', 
                            # 'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_size',
                            # 'beta_I_alt', 
                            # 'beta_I_isolation',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_latitude_pred <- 
  lapply(c('alpha', 
           'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           # 'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_latitude_pred[, grep(x, colnames(post_latitude_pred))]
           })

names(post_latitude_pred) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_size',
                               # 'beta_I_alt', 
                               # 'beta_I_isolation',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')



# ======= Plots =====
# 
# ============= Slops ===========

est_latitude_tot <- cond_effects(posterior = post_latitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'random',
                                 n = 100)

plot(dat$lat, dat$total_remotion, cex = 0.1)
est_latitude_tot %$% lines(x, y)
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)

est_latitude_tot <- cond_effects(posterior = post_latitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'averaging',
                                 n = 100)

plot(NULL, xlim = c(-2.2, 1.8), ylim = c(0, 1))
est_latitude_tot %$% lines(x, y)
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)


# error bars

rbind(pivot_longer(post_latitude_tot$beta, 'beta_lat') |> 
        mutate(type = 'Frugivory', 
               effect = 'Latitude'), 
      pivot_longer(post_latitude_disp$beta, 'beta_lat') |> 
        mutate(type = 'Seed dispersion', 
               effect = 'Latitude'), 
      pivot_longer(post_latitude_pred$beta, 'beta_lat') |> 
        mutate(type = 'Seed predation', 
               effect = 'Latitude'), 
      pivot_longer(post_latitude_bird$beta, 'beta_lat') |> 
        mutate(type = 'Birds', 
               effect = 'Latitude'), 
      pivot_longer(post_latitude_lizard$beta, 'beta_lat') |> 
        mutate(type = 'Lizards', 
               effect = 'Latitude')) |> 
  group_by(type) |> 
  transmute(mu = median(value), 
            li = quantile(value, 0.025), 
            ls = quantile(value, 0.975), 
            x = 'Slope', 
            effect = effect) |> 
  unique() |> 
  ggplot(aes(type, mu, ymin = li, ymax = ls)) +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_wrap(~ effect) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs()


# ============= *Type island* ======


TI_tot <- average_effects(n_levels = ncol(post_latitude_tot$TI),
                          posterior = post_latitude_tot, 
                          x_var1 = 'island_type', 
                          x_var2 = 'realm', 
                          par1 = 'TI', 
                          par2 = 'p_realm')

TI_bird <- average_effects(n_levels = ncol(post_latitude_bird$TI),
                          posterior = post_latitude_bird, 
                          x_var1 = 'island_type', 
                          x_var2 = 'realm', 
                          par1 = 'TI', 
                          par2 = 'p_realm')

TI_lizard <- average_effects(n_levels = ncol(post_latitude_lizard$TI),
                          posterior = post_latitude_lizard, 
                          x_var1 = 'island_type', 
                          x_var2 = 'realm', 
                          par1 = 'TI', 
                          par2 = 'p_realm')

TI_disp <- average_effects(n_levels = ncol(post_latitude_disp$TI),
                           posterior = post_latitude_disp, 
                           x_var1 = 'island_type', 
                           x_var2 = 'realm', 
                           par1 = 'TI', 
                           par2 = 'p_realm')

TI_pred <- average_effects(n_levels = ncol(post_latitude_pred$TI),
                           posterior = post_latitude_pred, 
                           x_var1 = 'island_type', 
                           x_var2 = 'realm', 
                           par1 = 'TI', 
                           par2 = 'p_realm')


rbind(full_join(TI_tot, codes$island_type, 'code') |> 
        mutate(type = 'Frugivory'), 
      full_join(TI_bird, codes$island_type, 'code') |> 
        mutate(type = 'Birds'),
      full_join(TI_lizard, codes$island_type, 'code') |> 
        mutate(type = 'Lizards'),
      full_join(TI_disp, codes$island_type, 'code') |> 
        mutate(type = 'Seed dispersal'), 
      full_join(TI_pred, codes$island_type, 'code') |> 
        mutate(type = 'Seed predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> 
  ggplot(aes(island, mu, ymin = li, ymax = ls, color = type)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2)) +
  lims(y = c(0, 0.85)) +
  labs(y = 'P(fruit consumption)', x = 'Type of island')


# ====== *Realm* =====

realm_tot <- average_effects(n_levels = ncol(post_latitude_tot$p_realm),
                          posterior = post_latitude_tot, 
                          x_var1 = 'realm', 
                          x_var2 = 'island_type', 
                          par1 = 'p_realm', 
                          par2 = 'TI')

realm_disp <- average_effects(n_levels = ncol(post_latitude_disp$p_realm),
                           posterior = post_latitude_disp, 
                           x_var1 = 'realm', 
                           x_var2 = 'island_type', 
                           par1 = 'p_realm', 
                           par2 = 'TI')

realm_bird <- average_effects(n_levels = ncol(post_latitude_bird$p_realm),
                              posterior = post_latitude_bird, 
                              x_var1 = 'realm', 
                              x_var2 = 'island_type', 
                              par1 = 'p_realm', 
                              par2 = 'TI')

realm_lizard <- average_effects(n_levels = ncol(post_latitude_lizard$p_realm),
                              posterior = post_latitude_lizard, 
                              x_var1 = 'realm', 
                              x_var2 = 'island_type', 
                              par1 = 'p_realm', 
                              par2 = 'TI')

realm_pred <- average_effects(n_levels = ncol(post_latitude_pred$p_realm),
                           posterior = post_latitude_pred, 
                           x_var1 = 'realm', 
                           x_var2 = 'island_type', 
                           par1 = 'p_realm', 
                           par2 = 'TI')


rbind(full_join(realm_tot, codes$real, 'code') |> 
        mutate(type = 'Frugivory'),
      full_join(realm_bird, codes$real, 'code') |> 
        mutate(type = 'Birds'),
      full_join(realm_lizard, codes$real, 'code') |> 
        mutate(type = 'Lizards'),
      full_join(realm_disp, codes$real, 'code') |> 
        mutate(type = 'Seed dispersal'), 
      full_join(realm_pred, codes$real, 'code') |> 
        mutate(type = 'Seed predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> 
  ggplot(aes(island, mu, ymin = li, ymax = ls)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.3)) +
  geom_point(position = position_dodge(width = 0.3)) +
  labs(y = 'P(fruit consumption)', x = 'Biogeographic realm') +
  facet_wrap(~type, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# ======== ***Regional effects*** ======

# =============== *Isolation* ==========

# =============== Overall frugivory  ======================

file <- paste0(getwd(), '/mod_isolation_total.stan')
fit_isolation_tot <- cmdstan_model(file, compile = T)

mod_isolation_tot <- 
  fit_isolation_tot$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_isolation_tot$save_object('mod_isolation_tot.rds')

mod_isolation_tot <- readRDS('mod_isolation_tot.rds')

mcmc_trace(mod_isolation_tot$draws(c('alpha', 
                                     #'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     'beta_I_isolation'
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
                                     )))

sum_isolation_tot <- mod_isolation_tot$summary()
mod_diagnostics(mod_isolation_tot, sum_isolation_tot)
ppcheck_isolation_tot <- mod_isolation_tot$draws('ppcheck', format = 'matrix')

plot(density(dat$total_remotion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.1))
for (i in 1:200) lines(density(ppcheck_isolation_tot[i, ], lwd = 0.1))
lines(density(dat$total_remotion), lwd = 2, col = 'red')


post_isolation_tot <- 
  mod_isolation_tot$draws(c('alpha', 
                           #'beta_lat', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_size',
                           # 'beta_I_alt', 
                           'beta_I_isolation',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome'), 
                         format = 'df')

post_isolation_tot <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_isolation_tot[, grep(x, colnames(post_isolation_tot))]
           })

names(post_isolation_tot) <- c('alpha', 
                              'beta', 
                              # 'beta_H_pop', 
                              # 'beta_H_foot',
                              # 'beta_I_mainland', 
                              # 'beta_I_size',
                              # 'beta_I_alt', 
                              #'beta_I_isolation',
                              # 'beta_temp', 
                              # 'beta_NV',
                              # 'beta_bush', 
                              # 'inv_rank', 
                              'TI', 'p_island', 
                              'p_country', 'p_grid', 
                              'p_plant', 'p_realm', 
                              'p_ecoR', 'p_biome')



# =============== Fruit dispersion  ======================

file <- paste0(getwd(), '/mod_isolation_dispersion.stan')
fit_isolation_disp <- cmdstan_model(file, compile = T)

mod_isolation_disp <- 
  fit_isolation_disp$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_isolation_disp$save_object('mod_isolation_disp.rds')

mod_isolation_disp <- readRDS('mod_isolation_disp.rds')

mcmc_trace(mod_isolation_disp$draws(c('alpha', 
                                     #'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     'beta_I_isolation'
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
)))

sum_isolation_disp <- mod_isolation_disp$summary()
mod_diagnostics(mod_isolation_disp, sum_isolation_disp)
ppcheck_isolation_disp <- mod_isolation_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_isolation_disp[i, ], lwd = 0.1))
lines(density(dat$dispersion), lwd = 2, col = 'red')

post_isolation_disp <- 
  mod_isolation_disp$draws(c('alpha', 
                            # 'beta_lat', 
                            # 'beta_H_pop', 
                            # 'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_size',
                            # 'beta_I_alt', 
                            'beta_I_isolation',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_isolation_disp <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_isolation_disp[, grep(x, colnames(post_isolation_disp))]
           })

names(post_isolation_disp) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_size',
                               # 'beta_I_alt', 
                               # 'beta_I_isolation',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')


# =============== Birds  ======================

file <- paste0(getwd(), '/mod_isolation_bird.stan')
fit_isolation_bird <- cmdstan_model(file, compile = T)

mod_isolation_bird <- 
  fit_isolation_bird$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_isolation_bird$save_object('mod_isolation_bird.rds')

mod_isolation_bird <- readRDS('mod_isolation_bird.rds')

mcmc_trace(mod_isolation_bird$draws(c('alpha', 
                                     #'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     'beta_I_isolation'
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
)))

sum_isolation_bird <- mod_isolation_bird$summary()
mod_diagnostics(mod_isolation_bird, sum_isolation_bird)
ppcheck_isolation_bird <- mod_isolation_bird$draws('ppcheck', format = 'matrix')

plot(density(dat$Bird), main = '', 
     xlab = 'Bird fruits removal', ylim = c(0, 3))
for (i in 1:200) lines(density(ppcheck_isolation_bird[i, ], lwd = 0.1))
lines(density(dat$Bird), lwd = 2, col = 'red')

post_isolation_bird <- 
  mod_isolation_bird$draws(c('alpha', 
                             # 'beta_lat', 
                             # 'beta_H_pop', 
                             # 'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_size',
                             # 'beta_I_alt', 
                             'beta_I_isolation',
                             # 'beta_temp', 
                             # 'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome'), 
                           format = 'df')

post_isolation_bird <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_isolation_bird[, grep(x, colnames(post_isolation_bird))]
           })

names(post_isolation_bird) <- c('alpha', 
                                'beta', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                # 'beta_I_isolation',
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                'TI', 'p_island', 
                                'p_country', 'p_grid', 
                                'p_plant', 'p_realm', 
                                'p_ecoR', 'p_biome')


# =============== Lizards  ======================

file <- paste0(getwd(), '/mod_isolation_lizard.stan')
fit_isolation_lizard <- cmdstan_model(file, compile = T)

mod_isolation_lizard <- 
  fit_isolation_lizard$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

#mod_isolation_lizard$save_object('mod_isolation_lizard.rds')

mod_isolation_lizard <- readRDS('mod_isolation_lizard.rds')

mcmc_trace(mod_isolation_lizard$draws(c('alpha', 
                                     #'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     'beta_I_isolation'
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
)))

sum_isolation_lizard <- mod_isolation_lizard$summary()
mod_diagnostics(mod_isolation_lizard, sum_isolation_lizard)
ppcheck_isolation_lizard <- mod_isolation_lizard$draws('ppcheck', format = 'matrix')

plot(density(dat$Lizard), main = '', 
     xlab = 'Lizard fruits removal', ylim = c(0, 6))
for (i in 1:200) lines(density(ppcheck_isolation_lizard[i, ], lwd = 0.1))
lines(density(dat$Lizard), lwd = 2, col = 'red')

post_isolation_lizard <- 
  mod_isolation_lizard$draws(c('alpha', 
                             # 'beta_lat', 
                             # 'beta_H_pop', 
                             # 'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_size',
                             # 'beta_I_alt', 
                             'beta_I_isolation',
                             # 'beta_temp', 
                             # 'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome'), 
                           format = 'df')

post_isolation_lizard <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_isolation_lizard[, grep(x, colnames(post_isolation_lizard))]
           })

names(post_isolation_lizard) <- c('alpha', 
                                'beta', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                # 'beta_I_isolation',
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                'TI', 'p_island', 
                                'p_country', 'p_grid', 
                                'p_plant', 'p_realm', 
                                'p_ecoR', 'p_biome')






# =============== Fruit predation  ======================

file <- paste0(getwd(), '/mod_isolation_predation.stan')
fit_isolation_pred <- cmdstan_model(file, compile = T)

mod_isolation_pred <- 
  fit_isolation_pred$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_isolation_pred$save_object('mod_isolation_pred.rds')

mod_isolation_pred <- readRDS('mod_isolation_pred.rds')

mcmc_trace(mod_isolation_pred$draws(c('alpha', 
                                     #'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_size',
                                     # 'beta_I_alt', 
                                     'beta_I_isolation'
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
)))

sum_isolation_pred <- mod_isolation_pred$summary()
mod_diagnostics(mod_isolation_pred, sum_isolation_pred)

ppcheck_isolation_pred <- mod_isolation_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_isolation_pred[i, ], lwd = 0.1))
lines(density(dat$predation), lwd = 2, col = 'red')

post_isolation_pred <- 
  mod_isolation_pred$draws(c('alpha', 
                            # 'beta_lat', 
                            # 'beta_H_pop', 
                            # 'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_size',
                            # 'beta_I_alt', 
                            'beta_I_isolation',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_isolation_pred <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_isolation',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_isolation_pred[, grep(x, colnames(post_isolation_pred))]
           })

names(post_isolation_pred) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_size',
                               # 'beta_I_alt', 
                               # 'beta_I_isolation',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')



# ======= Plots =====
# 
# ============= Slops ===========

est_latitude_tot <- cond_effects(posterior = post_altitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'random',
                                 n = 100)

plot(dat$lat, dat$total_remotion, cex = 0.1)
est_latitude_tot %$% lines(x, y)
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)

est_latitude_tot %$% plot(x, y, type = 'l', ylim = c(0, 1))
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)


est_latitude_tot <- cond_effects(posterior = post_altitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'averaging',
                                 n = 100)

plot(NULL, xlim = c(-2.2, 1.8), ylim = c(0, 1))
est_latitude_tot %$% lines(x, y)
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)


# error bars

rbind(pivot_longer(post_isolation_tot$beta, 'beta_I_isolation') |> 
        mutate(type = 'Frugivory', 
               effect = 'Island isolation'), 
      pivot_longer(post_isolation_bird$beta, 'beta_I_isolation') |> 
        mutate(type = 'Birds', 
               effect = 'Island isolation'), 
      pivot_longer(post_isolation_lizard$beta, 'beta_I_isolation') |> 
        mutate(type = 'Lizards', 
               effect = 'Island isolation'), 
      pivot_longer(post_isolation_disp$beta, 'beta_I_isolation') |> 
        mutate(type = 'Seed dispersion', 
               effect = 'Island isolation'), 
      pivot_longer(post_isolation_pred$beta, 'beta_I_isolation') |> 
        mutate(type = 'Seed predation', 
               effect = 'Island isolation')) |> 
  group_by(type) |> 
  transmute(mu = median(value), 
            li = quantile(value, 0.025), 
            ls = quantile(value, 0.975), 
            x = 'Slope', 
            effect = effect) |> 
  unique() |> 
  ggplot(aes(type, mu, ymin = li, ymax = ls)) +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_wrap(~ effect) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs()





# =============== *Size* ==========

# =============== Overall frugivory  ======================

file <- paste0(getwd(), '/mod_size_total.stan')
fit_size_tot <- cmdstan_model(file, compile = T)

mod_size_tot <- 
  fit_size_tot$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_size_tot$save_object('mod_size_tot.rds')

mod_size_tot <- readRDS('mod_size_tot.rds')

mcmc_trace(mod_size_tot$draws(c('alpha', 
                                #'beta_lat', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                'beta_I_size'
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                )))

sum_size_tot <- mod_size_tot$summary()
mod_diagnostics(mod_size_tot, sum_size_tot)
ppcheck_size_tot <- mod_size_tot$draws('ppcheck', format = 'matrix')

plot(density(dat$total_remotion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.1))
for (i in 1:200) lines(density(ppcheck_size_tot[i, ], lwd = 0.1))
lines(density(dat$total_remotion), lwd = 2, col = 'red')


post_size_tot <- 
  mod_size_tot$draws(c('alpha', 
                            #'beta_lat', 
                            # 'beta_H_pop', 
                            # 'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_size',
                            # 'beta_I_alt', 
                            'beta_I_size',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_size_tot <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_size',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_size_tot[, grep(x, colnames(post_size_tot))]
           })

names(post_size_tot) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_size',
                               # 'beta_I_alt', 
                               #'beta_I_size',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')



# =============== Fruit dispersion  ======================

file <- paste0(getwd(), '/mod_size_dispersion.stan')
fit_size_disp <- cmdstan_model(file, compile = T)

mod_size_disp <- 
  fit_size_disp$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_size_disp$save_object('mod_size_disp.rds')

mod_size_disp <- readRDS('mod_size_disp.rds')

mcmc_trace(mod_size_tot$draws(c('alpha', 
                                #'beta_lat', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                'beta_I_size'
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
)))

sum_size_disp <- mod_size_disp$summary()
mod_diagnostics(mod_size_disp, sum_size_disp)
ppcheck_size_disp <- mod_size_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_size_disp[i, ], lwd = 0.1))
lines(density(dat$dispersion), lwd = 2, col = 'red')

post_size_disp <- 
  mod_size_disp$draws(c('alpha', 
                             # 'beta_lat', 
                             # 'beta_H_pop', 
                             # 'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_size',
                             # 'beta_I_alt', 
                             'beta_I_size',
                             # 'beta_temp', 
                             # 'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome'), 
                           format = 'df')

post_size_disp <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_size',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_size_disp[, grep(x, colnames(post_size_disp))]
           })

names(post_size_disp) <- c('alpha', 
                                'beta', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                # 'beta_I_size',
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                'TI', 'p_island', 
                                'p_country', 'p_grid', 
                                'p_plant', 'p_realm', 
                                'p_ecoR', 'p_biome')


# =============== Birds  ======================

file <- paste0(getwd(), '/mod_size_bird.stan')
fit_size_bird <- cmdstan_model(file, compile = T)

mod_size_bird <- 
  fit_size_bird$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

#mod_size_bird$save_object('mod_size_bird.rds')

mod_size_bird <- readRDS('mod_size_bird.rds')

mcmc_trace(mod_size_bird$draws(c('alpha', 
                                #'beta_lat', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                'beta_I_size'
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
)))

sum_size_bird <- mod_size_bird$summary()
mod_diagnostics(mod_size_bird, sum_size_bird)
ppcheck_size_bird <- mod_size_bird$draws('ppcheck', format = 'matrix')

plot(density(dat$Bird), main = '', 
     xlab = 'Bird fruits removal', ylim = c(0, 3))
for (i in 1:200) lines(density(ppcheck_size_bird[i, ], lwd = 0.1))
lines(density(dat$Bird), lwd = 2, col = 'red')

post_size_bird <- 
  mod_size_bird$draws(c('alpha', 
                        # 'beta_lat', 
                        # 'beta_H_pop', 
                        # 'beta_H_foot',
                        # 'beta_I_mainland', 
                        # 'beta_I_size',
                        # 'beta_I_alt', 
                        'beta_I_size',
                        # 'beta_temp', 
                        # 'beta_NV',
                        # 'beta_bush', 
                        # 'inv_rank', 
                        'TI', 'p_island', 
                        'p_country', 'p_grid', 
                        'p_plant', 'p_realm', 
                        'p_ecoR', 'p_biome'), 
                      format = 'df')

post_size_bird <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_size',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_size_bird[, grep(x, colnames(post_size_bird))]
           })

names(post_size_bird) <- c('alpha', 
                           'beta', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_size',
                           # 'beta_I_alt', 
                           # 'beta_I_size',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome')


# =============== Lizards  ======================

file <- paste0(getwd(), '/mod_size_lizard.stan')
fit_size_lizard <- cmdstan_model(file, compile = T)

mod_size_lizard <- 
  fit_size_lizard$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

#mod_size_lizard$save_object('mod_size_lizard.rds')

mod_size_lizard <- readRDS('mod_size_lizard.rds')

mcmc_trace(mod_size_lizard$draws(c('alpha', 
                                #'beta_lat', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                'beta_I_size'
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
)))

sum_size_lizard <- mod_size_lizard$summary()
mod_diagnostics(mod_size_lizard, sum_size_lizard)
ppcheck_size_lizard <- mod_size_lizard$draws('ppcheck', format = 'matrix')

plot(density(dat$Lizard), main = '', 
     xlab = 'Lizard fruits removal', ylim = c(0, 6))
for (i in 1:200) lines(density(ppcheck_size_lizard[i, ], lwd = 0.1))
lines(density(dat$Lizard), lwd = 2, col = 'red')

post_size_lizard <- 
  mod_size_lizard$draws(c('alpha', 
                        # 'beta_lat', 
                        # 'beta_H_pop', 
                        # 'beta_H_foot',
                        # 'beta_I_mainland', 
                        # 'beta_I_size',
                        # 'beta_I_alt', 
                        'beta_I_size',
                        # 'beta_temp', 
                        # 'beta_NV',
                        # 'beta_bush', 
                        # 'inv_rank', 
                        'TI', 'p_island', 
                        'p_country', 'p_grid', 
                        'p_plant', 'p_realm', 
                        'p_ecoR', 'p_biome'), 
                      format = 'df')

post_size_lizard <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_size',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_size_lizard[, grep(x, colnames(post_size_lizard))]
           })

names(post_size_lizard) <- c('alpha', 
                           'beta', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_size',
                           # 'beta_I_alt', 
                           # 'beta_I_size',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome')






# =============== Fruit predation  ======================

file <- paste0(getwd(), '/mod_size_predation.stan')
fit_size_pred <- cmdstan_model(file, compile = T)

mod_size_pred <- 
  fit_size_pred$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_size_pred$save_object('mod_size_pred.rds')

mod_size_pred <- readRDS('mod_size_pred.rds')

mcmc_trace(mod_size_pred$draws(c('alpha', 
                                #'beta_lat', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                'beta_I_size'
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
)))

sum_size_pred <- mod_size_pred$summary()
mod_diagnostics(mod_size_pred, sum_size_pred)

ppcheck_size_pred <- mod_size_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Predation', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_size_pred[i, ], lwd = 0.1))
lines(density(dat$predation), lwd = 2, col = 'red')

post_size_pred <- 
  mod_size_pred$draws(c('alpha', 
                             # 'beta_lat', 
                             # 'beta_H_pop', 
                             # 'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_size',
                             # 'beta_I_alt', 
                             'beta_I_size',
                             # 'beta_temp', 
                             # 'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome'), 
                           format = 'df')

post_size_pred <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_size',
           # 'beta_I_alt', 
           'beta_I_size',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_size_pred[, grep(x, colnames(post_size_pred))]
           })

names(post_size_pred) <- c('alpha', 
                                'beta', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_size',
                                # 'beta_I_alt', 
                                # 'beta_I_size',
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                'TI', 'p_island', 
                                'p_country', 'p_grid', 
                                'p_plant', 'p_realm', 
                                'p_ecoR', 'p_biome')



# ======= Plots =====
# 
# ============= Slops ===========

est_latitude_tot <- cond_effects(posterior = post_altitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'random',
                                 n = 100)

plot(dat$lat, dat$total_remotion, cex = 0.1)
est_latitude_tot %$% lines(x, y)
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)

est_latitude_tot %$% plot(x, y, type = 'l', ylim = c(0, 1))
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)


est_latitude_tot <- cond_effects(posterior = post_altitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'averaging',
                                 n = 100)

plot(NULL, xlim = c(-2.2, 1.8), ylim = c(0, 1))
est_latitude_tot %$% lines(x, y)
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)


# error bars

rbind(pivot_longer(post_size_tot$beta, 'beta_I_size') |> 
        mutate(type = 'Frugivory', 
               effect = 'Island size'), 
      pivot_longer(post_size_bird$beta, 'beta_I_size') |> 
        mutate(type = 'Birds', 
               effect = 'Island size'), 
      pivot_longer(post_size_lizard$beta, 'beta_I_size') |> 
        mutate(type = 'Lizards', 
               effect = 'Island size'), 
      pivot_longer(post_size_disp$beta, 'beta_I_size') |> 
        mutate(type = 'Seed dispersion', 
               effect = 'Island size'), 
      pivot_longer(post_size_pred$beta, 'beta_I_size') |> 
        mutate(type = 'Seed predation', 
               effect = 'Island size')) |> 
  group_by(type) |> 
  transmute(mu = median(value), 
            li = quantile(value, 0.025), 
            ls = quantile(value, 0.975), 
            x = 'Slope', 
            effect = effect) |> 
  unique() |> 
  ggplot(aes(type, mu, ymin = li, ymax = ls)) +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_wrap(~ effect) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs()







# ============ ***Landscape effects*** =====

# =============== *Altitude* ==========

# =============== Overall frugivory  ======================

file <- paste0(getwd(), '/mod_altitude_total.stan')
fit_altitude_tot <- cmdstan_model(file, compile = T)

mod_altitude_tot <- 
  fit_altitude_tot$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_altitude_tot$save_object('mod_altitude_tot.rds')

mod_altitude_tot <- readRDS('mod_altitude_tot.rds')

sum_altitude_tot <- mod_altitude_tot$summary()
mod_diagnostics(mod_altitude_tot, sum_altitude_tot)
ppcheck_altitude_tot <- mod_altitude_tot$draws('ppcheck', format = 'matrix')

plot(density(dat$total_remotion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.1))
for (i in 1:200) lines(density(ppcheck_altitude_tot[i, ], lwd = 0.1))
lines(density(dat$total_remotion), lwd = 2, col = 'red')


post_altitude_tot <- 
  mod_altitude_tot$draws(c('alpha', 
                       #'beta_lat', 
                       # 'beta_H_pop', 
                       # 'beta_H_foot',
                       # 'beta_I_mainland', 
                       # 'beta_I_altitude',
                       # 'beta_I_alt', 
                       'beta_I_alt',
                       # 'beta_temp', 
                       # 'beta_NV',
                       # 'beta_bush', 
                       # 'inv_rank', 
                       'TI', 'p_island', 
                       'p_country', 'p_grid', 
                       'p_plant', 'p_realm', 
                       'p_ecoR', 'p_biome'), 
                     format = 'df')

post_altitude_tot <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_altitude',
           # 'beta_I_alt', 
           'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_altitude_tot[, grep(x, colnames(post_altitude_tot))]
           })

names(post_altitude_tot) <- c('alpha', 
                          'beta', 
                          # 'beta_H_pop', 
                          # 'beta_H_foot',
                          # 'beta_I_mainland', 
                          # 'beta_I_altitude',
                          # 'beta_I_alt', 
                          #'beta_I_altitude',
                          # 'beta_temp', 
                          # 'beta_NV',
                          # 'beta_bush', 
                          # 'inv_rank', 
                          'TI', 'p_island', 
                          'p_country', 'p_grid', 
                          'p_plant', 'p_realm', 
                          'p_ecoR', 'p_biome')



# =============== Fruit dispersion  ======================

file <- paste0(getwd(), '/mod_altitude_dispersion.stan')
fit_altitude_disp <- cmdstan_model(file, compile = T)

mod_altitude_disp <- 
  fit_altitude_disp$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

mod_altitude_disp$save_object('mod_altitude_disp.rds')

mod_altitude_disp <- readRDS('mod_altitude_disp.rds')

mcmc_trace(mod_altitude_disp$draws(c('alpha', 
                                     # 'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_altitude',
                                     # 'beta_I_alt', 
                                     'beta_I_alt',
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
                                     )))

sum_altitude_disp <- mod_altitude_disp$summary()
mod_diagnostics(mod_altitude_disp, sum_altitude_disp)
ppcheck_altitude_disp <- mod_altitude_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Seed dispersal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_altitude_disp[i, ], lwd = 0.1))
lines(density(dat$dispersion), lwd = 2, col = 'red')

post_altitude_disp <- 
  mod_altitude_disp$draws(c('alpha', 
                        # 'beta_lat', 
                        # 'beta_H_pop', 
                        # 'beta_H_foot',
                        # 'beta_I_mainland', 
                        # 'beta_I_altitude',
                        # 'beta_I_alt', 
                        'beta_I_alt',
                        # 'beta_temp', 
                        # 'beta_NV',
                        # 'beta_bush', 
                        # 'inv_rank', 
                        'TI', 'p_island', 
                        'p_country', 'p_grid', 
                        'p_plant', 'p_realm', 
                        'p_ecoR', 'p_biome'), 
                      format = 'df')

post_altitude_disp <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_altitude',
           # 'beta_I_alt', 
           'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_altitude_disp[, grep(x, colnames(post_altitude_disp))]
           })

names(post_altitude_disp) <- c('alpha', 
                           'beta', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_altitude',
                           # 'beta_I_alt', 
                           # 'beta_I_altitude',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome')


# =============== Birds  ======================

file <- paste0(getwd(), '/mod_altitude_bird.stan')
fit_altitude_bird <- cmdstan_model(file, compile = T)

mod_altitude_bird <- 
  fit_altitude_bird$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

#mod_altitude_bird$save_object('mod_altitude_bird.rds')

mod_altitude_bird <- readRDS('mod_altitude_bird.rds')

mcmc_trace(mod_altitude_bird$draws(c('alpha', 
                                     # 'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_altitude',
                                     # 'beta_I_alt', 
                                     'beta_I_alt',
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
)))

sum_altitude_bird <- mod_altitude_bird$summary()
mod_diagnostics(mod_altitude_bird, sum_altitude_bird)
ppcheck_altitude_bird <- mod_altitude_bird$draws('ppcheck', format = 'matrix')

plot(density(dat$Bird), main = '', 
     xlab = 'Bird seed dispersal', ylim = c(0, 3))
for (i in 1:200) lines(density(ppcheck_altitude_bird[i, ], lwd = 0.1))
lines(density(dat$Bird), lwd = 2, col = 'red')

post_altitude_bird <- 
  mod_altitude_bird$draws(c('alpha', 
                            # 'beta_lat', 
                            # 'beta_H_pop', 
                            # 'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_altitude',
                            # 'beta_I_alt', 
                            'beta_I_alt',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_altitude_bird <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_altitude',
           # 'beta_I_alt', 
           'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_altitude_bird[, grep(x, colnames(post_altitude_bird))]
           })

names(post_altitude_bird) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_altitude',
                               # 'beta_I_alt', 
                               # 'beta_I_altitude',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')


# =============== Lizards  ======================

file <- paste0(getwd(), '/mod_altitude_lizard.stan')
fit_altitude_lizard <- cmdstan_model(file, compile = T)

mod_altitude_lizard <- 
  fit_altitude_lizard$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

#mod_altitude_lizard$save_object('mod_altitude_lizard.rds')

mod_altitude_lizard <- readRDS('mod_altitude_lizard.rds')

mcmc_trace(mod_altitude_lizard$draws(c('alpha', 
                                     # 'beta_lat', 
                                     # 'beta_H_pop', 
                                     # 'beta_H_foot',
                                     # 'beta_I_mainland', 
                                     # 'beta_I_altitude',
                                     # 'beta_I_alt', 
                                     'beta_I_alt',
                                     # 'beta_temp', 
                                     # 'beta_NV',
                                     # 'beta_bush', 
                                     # 'inv_rank', 
)))

sum_altitude_lizard <- mod_altitude_lizard$summary()
mod_diagnostics(mod_altitude_lizard, sum_altitude_lizard)
ppcheck_altitude_lizard <- mod_altitude_lizard$draws('ppcheck', format = 'matrix')

plot(density(dat$Lizard), main = '', 
     xlab = 'Lizard seed dispersal', ylim = c(0, 6))
for (i in 1:200) lines(density(ppcheck_altitude_lizard[i, ], lwd = 0.1))
lines(density(dat$Lizard), lwd = 2, col = 'red')

post_altitude_lizard <- 
  mod_altitude_lizard$draws(c('alpha', 
                            # 'beta_lat', 
                            # 'beta_H_pop', 
                            # 'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_altitude',
                            # 'beta_I_alt', 
                            'beta_I_alt',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_altitude_lizard <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_altitude',
           # 'beta_I_alt', 
           'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_altitude_lizard[, grep(x, colnames(post_altitude_lizard))]
           })

names(post_altitude_lizard) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_altitude',
                               # 'beta_I_alt', 
                               # 'beta_I_altitude',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')





# =============== Fruit predation  ======================

file <- paste0(getwd(), '/mod_altitude_predation.stan')
fit_altitude_pred <- cmdstan_model(file, compile = T)

mod_altitude_pred <- 
  fit_altitude_pred$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_altitude_pred$save_object('mod_altitude_pred.rds')

mod_altitude_pred <- readRDS('mod_altitude_pred.rds')

sum_altitude_pred <- mod_altitude_pred$summary()
mod_diagnostics(mod_altitude_pred, sum_altitude_pred)

ppcheck_altitude_pred <- mod_altitude_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Fruit predation', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_altitude_pred[i, ], lwd = 0.1))
lines(density(dat$predation), lwd = 2, col = 'red')

post_altitude_pred <- 
  mod_altitude_pred$draws(c('alpha', 
                        # 'beta_lat', 
                        # 'beta_H_pop', 
                        # 'beta_H_foot',
                        # 'beta_I_mainland', 
                        # 'beta_I_altitude',
                        # 'beta_I_alt', 
                        'beta_I_alt',
                        # 'beta_temp', 
                        # 'beta_NV',
                        # 'beta_bush', 
                        # 'inv_rank', 
                        'TI', 'p_island', 
                        'p_country', 'p_grid', 
                        'p_plant', 'p_realm', 
                        'p_ecoR', 'p_biome'), 
                      format = 'df')

post_altitude_pred <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_altitude',
           # 'beta_I_alt', 
           'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_altitude_pred[, grep(x, colnames(post_altitude_pred))]
           })

names(post_altitude_pred) <- c('alpha', 
                           'beta', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_altitude',
                           # 'beta_I_alt', 
                           # 'beta_I_altitude',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome')



# ======= Plots =====
# 
# ============= Slops ===========

# error bars

rbind(pivot_longer(post_altitude_tot$beta, 'beta_I_alt') |> 
        mutate(type = 'Frugivory', 
               effect = 'Island altitude'),
      pivot_longer(post_altitude_bird$beta, 'beta_I_alt') |> 
        mutate(type = 'Birds', 
               effect = 'Island altitude'),
      pivot_longer(post_altitude_lizard$beta, 'beta_I_alt') |> 
        mutate(type = 'Lizards', 
               effect = 'Island altitude'),
      pivot_longer(post_altitude_disp$beta, 'beta_I_alt') |> 
        mutate(type = 'Seed dispersion', 
               effect = 'Island altitude'), 
      pivot_longer(post_altitude_pred$beta, 'beta_I_alt') |> 
        mutate(type = 'Seed predation', 
               effect = 'Island altitude')) |> 
  group_by(type) |> 
  transmute(mu = median(value), 
            li = quantile(value, 0.025), 
            ls = quantile(value, 0.975), 
            x = 'Slope', 
            effect = effect) |> 
  unique() |> 
  ggplot(aes(type, mu, ymin = li, ymax = ls)) +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_wrap(~ effect) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs()

# Scatter plots

est_latitude_tot <- cond_effects(posterior = post_altitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'random',
                                 n = 100)

plot(dat$lat, dat$total_remotion, cex = 0.1)
est_latitude_tot %$% lines(x, y)
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)

est_latitude_tot %$% plot(x, y, type = 'l', ylim = c(0, 1))
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)


est_latitude_tot <- cond_effects(posterior = post_altitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'averaging',
                                 n = 100)

plot(NULL, xlim = c(-2.2, 1.8), ylim = c(0, 1))
est_latitude_tot %$% lines(x, y)
est_latitude_tot %$% lines(x, li, lty = 3)
est_latitude_tot %$% lines(x, ls, lty = 3)








# =============== *Humman footprint* ==========

# =============== Overall frugivory  ======================

file <- paste0(getwd(), '/mod_footprint_total.stan')
fit_footprint_tot <- cmdstan_model(file, compile = T)

mod_footprint_tot <- 
  fit_footprint_tot$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_footprint_tot$save_object('mod_footprint_tot.rds')

mod_footprint_tot <- readRDS('mod_footprint_tot.rds')

mcmc_trace(mod_footprint_tot$draws(c('alpha', 
                                        # 'beta_lat', 
                                        # 'beta_H_pop', 
                                        'beta_H_foot',
                                        # 'beta_I_mainland', 
                                        # 'beta_I_footprint',
                                        # 'beta_I_alt', 
                                        # 'beta_I_alt',
                                        # 'beta_temp', 
                                        # 'beta_NV',
                                        # 'beta_bush', 
                                        # 'inv_rank', 
                                        )))

sum_footprint_tot <- mod_footprint_tot$summary()
mod_diagnostics(mod_footprint_tot, sum_footprint_tot)
ppcheck_footprint_tot <- mod_footprint_tot$draws('ppcheck', format = 'matrix')

plot(density(dat$total_remotion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.1))
for (i in 1:200) lines(density(ppcheck_footprint_tot[i, ], lwd = 0.1))
lines(density(dat$total_remotion), lwd = 2, col = 'red')


post_footprint_tot <- 
  mod_footprint_tot$draws(c('alpha', 
                           #'beta_lat', 
                           # 'beta_H_pop', 
                           'beta_H_foot',
                           # 'beta_I_mainland', 
                           #'beta_I_footprint',
                           # 'beta_I_alt', 
                           #'beta_I_alt',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome'), 
                         format = 'df')

post_footprint_tot <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_footprint',
           # 'beta_I_alt', 
           'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_footprint_tot[, grep(x, colnames(post_footprint_tot))]
           })

names(post_footprint_tot) <- c('alpha', 
                              'beta', 
                              # 'beta_H_pop', 
                              # 'beta_H_foot',
                              # 'beta_I_mainland', 
                              # 'beta_I_footprint',
                              # 'beta_I_alt', 
                              #'beta_I_footprint',
                              # 'beta_temp', 
                              # 'beta_NV',
                              # 'beta_bush', 
                              # 'inv_rank', 
                              'TI', 'p_island', 
                              'p_country', 'p_grid', 
                              'p_plant', 'p_realm', 
                              'p_ecoR', 'p_biome')



# =============== Fruit dispersion  ======================

file <- paste0(getwd(), '/mod_footprint_dispersion.stan')
fit_footprint_disp <- cmdstan_model(file, compile = T)

mod_footprint_disp <- 
  fit_footprint_disp$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

#mod_footprint_disp$save_object('mod_footprint_disp.rds')

mod_footprint_disp <- readRDS('mod_footprint_disp.rds')

mcmc_trace(mod_footprint_disp$draws(c('alpha', 
                                        # 'beta_lat', 
                                        # 'beta_H_pop', 
                                        'beta_H_foot',
                                        # 'beta_I_mainland', 
                                        # 'beta_I_footprint',
                                        # 'beta_I_alt', 
                                        # 'beta_I_alt',
                                        # 'beta_temp', 
                                        # 'beta_NV',
                                        # 'beta_bush', 
                                        # 'inv_rank', 
)))

sum_footprint_disp <- mod_footprint_disp$summary()
mod_diagnostics(mod_footprint_disp, sum_footprint_disp)
ppcheck_footprint_disp <- mod_footprint_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Seed dispersal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_footprint_disp[i, ], lwd = 0.1))
lines(density(dat$dispersion), lwd = 2, col = 'red')

post_footprint_disp <- 
  mod_footprint_disp$draws(c('alpha', 
                            # 'beta_lat', 
                            # 'beta_H_pop', 
                            'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_footprint',
                            # 'beta_I_alt', 
                            # 'beta_I_alt',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_footprint_disp <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_footprint',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_footprint_disp[, grep(x, colnames(post_footprint_disp))]
           })

names(post_footprint_disp) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_footprint',
                               # 'beta_I_alt', 
                               # 'beta_I_footprint',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')


# =============== Birds  ======================

file <- paste0(getwd(), '/mod_footprint_bird.stan')
fit_footprint_bird <- cmdstan_model(file, compile = T)

mod_footprint_bird <- 
  fit_footprint_bird$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_footprint_bird$save_object('mod_footprint_bird.rds')

mod_footprint_bird <- readRDS('mod_footprint_bird.rds')

mcmc_trace(mod_footprint_bird$draws(c('alpha', 
                                        # 'beta_lat', 
                                        # 'beta_H_pop', 
                                        'beta_H_foot',
                                        # 'beta_I_mainland', 
                                        # 'beta_I_footprint',
                                        # 'beta_I_alt', 
                                        # 'beta_I_alt',
                                        # 'beta_temp', 
                                        # 'beta_NV',
                                        # 'beta_bush', 
                                        # 'inv_rank', 
)))

sum_footprint_bird <- mod_footprint_bird$summary()
mod_diagnostics(mod_footprint_bird, sum_footprint_bird)
ppcheck_footprint_bird <- mod_footprint_bird$draws('ppcheck', format = 'matrix')

plot(density(dat$Bird), main = '', 
     xlab = 'Bird seed dispersal', ylim = c(0, 3))
for (i in 1:200) lines(density(ppcheck_footprint_bird[i, ], lwd = 0.1))
lines(density(dat$Bird), lwd = 2, col = 'red')

post_footprint_bird <- 
  mod_footprint_bird$draws(c('alpha', 
                             # 'beta_lat', 
                             # 'beta_H_pop', 
                             'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_footprint',
                             # 'beta_I_alt', 
                             # 'beta_I_alt',
                             # 'beta_temp', 
                             # 'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome'), 
                           format = 'df')

post_footprint_bird <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_footprint',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_footprint_bird[, grep(x, colnames(post_footprint_bird))]
           })

names(post_footprint_bird) <- c('alpha', 
                                'beta', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_footprint',
                                # 'beta_I_alt', 
                                # 'beta_I_footprint',
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                'TI', 'p_island', 
                                'p_country', 'p_grid', 
                                'p_plant', 'p_realm', 
                                'p_ecoR', 'p_biome')


# =============== Lizards  ======================

file <- paste0(getwd(), '/mod_footprint_lizard.stan')
fit_footprint_lizard <- cmdstan_model(file, compile = T)

mod_footprint_lizard <- 
  fit_footprint_lizard$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

mod_footprint_lizard$save_object('mod_footprint_lizard.rds')

mod_footprint_lizard <- readRDS('mod_footprint_lizard.rds')

mcmc_trace(mod_footprint_lizard$draws(c('alpha', 
                                        # 'beta_lat', 
                                        # 'beta_H_pop', 
                                        'beta_H_foot',
                                        # 'beta_I_mainland', 
                                        # 'beta_I_footprint',
                                        # 'beta_I_alt', 
                                        # 'beta_I_alt',
                                        # 'beta_temp', 
                                        # 'beta_NV',
                                        # 'beta_bush', 
                                        # 'inv_rank', 
)))

sum_footprint_lizard <- mod_footprint_lizard$summary()
mod_diagnostics(mod_footprint_lizard, sum_footprint_lizard)
ppcheck_footprint_lizard <- mod_footprint_lizard$draws('ppcheck', format = 'matrix')

plot(density(dat$Lizard), main = '', 
     xlab = 'Lizard seed dispersal', ylim = c(0, 6))
for (i in 1:200) lines(density(ppcheck_footprint_lizard[i, ], lwd = 0.1))
lines(density(dat$Lizard), lwd = 2, col = 'red')

post_footprint_lizard <- 
  mod_footprint_lizard$draws(c('alpha', 
                             # 'beta_lat', 
                             # 'beta_H_pop', 
                             'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_footprint',
                             # 'beta_I_alt', 
                             # 'beta_I_alt',
                             # 'beta_temp', 
                             # 'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome'), 
                           format = 'df')

post_footprint_lizard <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_footprint',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_footprint_lizard[, grep(x, colnames(post_footprint_lizard))]
           })

names(post_footprint_lizard) <- c('alpha', 
                                'beta', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_footprint',
                                # 'beta_I_alt', 
                                # 'beta_I_footprint',
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                'TI', 'p_island', 
                                'p_country', 'p_grid', 
                                'p_plant', 'p_realm', 
                                'p_ecoR', 'p_biome')





# =============== Fruit predation  ======================

file <- paste0(getwd(), '/mod_footprint_predation.stan')
fit_footprint_pred <- cmdstan_model(file, compile = T)

mod_footprint_pred <- 
  fit_footprint_pred$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_footprint_pred$save_object('mod_footprint_pred.rds')
mod_footprint_pred <- readRDS('mod_footprint_pred.rds')

mcmc_trace(mod_footprint_pred$draws(c('alpha', 
                                        # 'beta_lat', 
                                        # 'beta_H_pop', 
                                        'beta_H_foot',
                                        # 'beta_I_mainland', 
                                        # 'beta_I_footprint',
                                        # 'beta_I_alt', 
                                        # 'beta_I_alt',
                                        # 'beta_temp', 
                                        # 'beta_NV',
                                        # 'beta_bush', 
                                        # 'inv_rank', 
)))

sum_footprint_pred <- mod_footprint_pred$summary()
mod_diagnostics(mod_footprint_pred, sum_footprint_pred)

ppcheck_footprint_pred <- mod_footprint_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_footprint_pred[i, ], lwd = 0.1))
lines(density(dat$predation), lwd = 2, col = 'red')

post_footprint_pred <- 
  mod_footprint_pred$draws(c('alpha', 
                            # 'beta_lat', 
                            # 'beta_H_pop', 
                            'beta_H_foot',
                            # 'beta_I_mainland', 
                            # 'beta_I_footprint',
                            # 'beta_I_alt', 
                            #'beta_I_alt',
                            # 'beta_temp', 
                            # 'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_footprint_pred <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_footprint',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_footprint_pred[, grep(x, colnames(post_footprint_pred))]
           })

names(post_footprint_pred) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_footprint',
                               # 'beta_I_alt', 
                               # 'beta_I_footprint',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')



# ======= Plots =====
# 
# ============= Slops ===========

# error bars

rbind(pivot_longer(post_footprint_tot$beta, 'beta_H_foot') |> 
        mutate(type = 'Frugivory', 
               effect = 'Human footprint'),
      pivot_longer(post_footprint_bird$beta, 'beta_H_foot') |> 
        mutate(type = 'Birds', 
               effect = 'Human footprint'),
      pivot_longer(post_footprint_lizard$beta, 'beta_H_foot') |> 
        mutate(type = 'Lizard', 
               effect = 'Human footprint'),
      pivot_longer(post_footprint_disp$beta, 'beta_H_foot') |> 
        mutate(type = 'Seed dispersion', 
               effect = 'Human footprint'), 
      pivot_longer(post_footprint_pred$beta, 'beta_H_foot') |> 
        mutate(type = 'Seed predation', 
               effect = 'Human footprint')) |> 
  group_by(type) |> 
  transmute(mu = median(value), 
            li = quantile(value, 0.025), 
            ls = quantile(value, 0.975), 
            x = 'Slope', 
            effect = effect) |> 
  unique() |> 
  ggplot(aes(type, mu, ymin = li, ymax = ls)) +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_wrap(~ effect) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs()

# Scatter plots





# =============== *Native vegetation* ==========

# =============== Overall frugivory  ======================

file <- paste0(getwd(), '/mod_nativeV_total.stan')
fit_nativeV_tot <- cmdstan_model(file, compile = T)

mod_nativeV_tot <- 
  fit_nativeV_tot$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_nativeV_tot$save_object('mod_nativeV_tot.rds')
mod_nativeV_tot <- readRDS('mod_nativeV_tot.rds')

mcmc_trace(mod_nativeV_tot$draws(c('alpha', 
                                      # 'beta_lat', 
                                      # 'beta_H_pop', 
                                      #'beta_H_foot',
                                      # 'beta_I_mainland', 
                                      # 'beta_I_nativeV',
                                      # 'beta_I_alt', 
                                      # 'beta_I_alt',
                                      # 'beta_temp', 
                                      'beta_NV',
                                      # 'beta_bush', 
                                      # 'inv_rank', 
)))

sum_nativeV_tot <- mod_nativeV_tot$summary()
mod_diagnostics(mod_nativeV_tot, sum_nativeV_tot)
ppcheck_nativeV_tot <- mod_nativeV_tot$draws('ppcheck', format = 'matrix')

plot(density(dat$total_remotion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.1))
for (i in 1:200) lines(density(ppcheck_nativeV_tot[i, ], lwd = 0.1))
lines(density(dat$total_remotion), lwd = 2, col = 'red')



post_nativeV_tot <- 
  mod_nativeV_tot$draws(c('alpha', 
                            #'beta_lat', 
                            # 'beta_H_pop', 
                            # 'beta_H_foot',
                            # 'beta_I_mainland', 
                            #'beta_I_nativeV',
                            # 'beta_I_alt', 
                            #'beta_I_alt',
                            # 'beta_temp', 
                            'beta_NV',
                            # 'beta_bush', 
                            # 'inv_rank', 
                            'TI', 'p_island', 
                            'p_country', 'p_grid', 
                            'p_plant', 'p_realm', 
                            'p_ecoR', 'p_biome'), 
                          format = 'df')

post_nativeV_tot <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           #'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_nativeV',
           # 'beta_I_alt', 
           #'beta_I_alt',
           # 'beta_temp', 
           'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_nativeV_tot[, grep(x, colnames(post_nativeV_tot))]
           })

names(post_nativeV_tot) <- c('alpha', 
                               'beta', 
                               # 'beta_H_pop', 
                               # 'beta_H_foot',
                               # 'beta_I_mainland', 
                               # 'beta_I_nativeV',
                               # 'beta_I_alt', 
                               #'beta_I_nativeV',
                               # 'beta_temp', 
                               # 'beta_NV',
                               # 'beta_bush', 
                               # 'inv_rank', 
                               'TI', 'p_island', 
                               'p_country', 'p_grid', 
                               'p_plant', 'p_realm', 
                               'p_ecoR', 'p_biome')



# =============== Fruit dispersion  ======================

file <- paste0(getwd(), '/mod_nativeV_dispersion.stan')
fit_nativeV_disp <- cmdstan_model(file, compile = T)

mod_nativeV_disp <- 
  fit_nativeV_disp$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_nativeV_disp$save_object('mod_nativeV_disp.rds')

mod_nativeV_disp <- readRDS('mod_nativeV_disp.rds')

mcmc_trace(mod_nativeV_disp$draws(c('alpha', 
                                      # 'beta_lat', 
                                      # 'beta_H_pop', 
                                      #'beta_H_foot',
                                      # 'beta_I_mainland', 
                                      # 'beta_I_nativeV',
                                      # 'beta_I_alt', 
                                      # 'beta_I_alt',
                                      # 'beta_temp', 
                                      'beta_NV',
                                      # 'beta_bush', 
                                      # 'inv_rank', 
)))

sum_nativeV_disp <- mod_nativeV_disp$summary()
mod_diagnostics(mod_nativeV_disp, sum_nativeV_disp)
ppcheck_nativeV_disp <- mod_nativeV_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Seed dispersal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_nativeV_disp[i, ], lwd = 0.1))
lines(density(dat$dispersion), lwd = 2, col = 'red')

post_nativeV_disp <- 
  mod_nativeV_disp$draws(c('alpha', 
                             # 'beta_lat', 
                             # 'beta_H_pop', 
                             #'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_nativeV',
                             # 'beta_I_alt', 
                             # 'beta_I_alt',
                             # 'beta_temp', 
                             'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome'), 
                           format = 'df')

post_nativeV_disp <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_nativeV',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_nativeV_disp[, grep(x, colnames(post_nativeV_disp))]
           })

names(post_nativeV_disp) <- c('alpha', 
                                'beta', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_nativeV',
                                # 'beta_I_alt', 
                                # 'beta_I_nativeV',
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                'TI', 'p_island', 
                                'p_country', 'p_grid', 
                                'p_plant', 'p_realm', 
                                'p_ecoR', 'p_biome')

# =============== Birds  ======================

file <- paste0(getwd(), '/mod_nativeV_bird.stan')
fit_nativeV_bird <- cmdstan_model(file, compile = T)

mod_nativeV_bird <- 
  fit_nativeV_bird$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

mod_nativeV_bird$save_object('mod_nativeV_bird.rds')

mod_nativeV_bird <- readRDS('mod_nativeV_bird.rds')

mcmc_trace(mod_nativeV_bird$draws(c('alpha', 
                                      # 'beta_lat', 
                                      # 'beta_H_pop', 
                                      #'beta_H_foot',
                                      # 'beta_I_mainland', 
                                      # 'beta_I_nativeV',
                                      # 'beta_I_alt', 
                                      # 'beta_I_alt',
                                      # 'beta_temp', 
                                      'beta_NV'
                                      # 'beta_bush', 
                                      # 'inv_rank', 
)))

sum_nativeV_bird <- mod_nativeV_bird$summary()
mod_diagnostics(mod_nativeV_bird, sum_nativeV_bird)
ppcheck_nativeV_bird <- mod_nativeV_bird$draws('ppcheck', format = 'matrix')

plot(density(dat$Bird), main = '', 
     xlab = 'Bird seed dispersal', ylim = c(0, 3))
for (i in 1:200) lines(density(ppcheck_nativeV_bird[i, ], lwd = 0.1))
lines(density(dat$Bird), lwd = 2, col = 'red')

post_nativeV_bird <- 
  mod_nativeV_bird$draws(c('alpha', 
                           # 'beta_lat', 
                           # 'beta_H_pop', 
                           #'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_nativeV',
                           # 'beta_I_alt', 
                           # 'beta_I_alt',
                           # 'beta_temp', 
                           'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome'), 
                         format = 'df')

post_nativeV_bird <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_nativeV',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_nativeV_bird[, grep(x, colnames(post_nativeV_bird))]
           })

names(post_nativeV_bird) <- c('alpha', 
                              'beta', 
                              # 'beta_H_pop', 
                              # 'beta_H_foot',
                              # 'beta_I_mainland', 
                              # 'beta_I_nativeV',
                              # 'beta_I_alt', 
                              # 'beta_I_nativeV',
                              # 'beta_temp', 
                              # 'beta_NV',
                              # 'beta_bush', 
                              # 'inv_rank', 
                              'TI', 'p_island', 
                              'p_country', 'p_grid', 
                              'p_plant', 'p_realm', 
                              'p_ecoR', 'p_biome')


# =============== Lizards  ======================

file <- paste0(getwd(), '/mod_nativeV_lizard.stan')
fit_nativeV_lizard <- cmdstan_model(file, compile = T)

mod_nativeV_lizard <- 
  fit_nativeV_lizard$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

mod_nativeV_lizard$save_object('mod_nativeV_lizard.rds')

mod_nativeV_lizard <- readRDS('mod_nativeV_lizard.rds')

mcmc_trace(mod_nativeV_lizard$draws(c('alpha', 
                                      # 'beta_lat', 
                                      # 'beta_H_pop', 
                                      #'beta_H_foot',
                                      # 'beta_I_mainland', 
                                      # 'beta_I_nativeV',
                                      # 'beta_I_alt', 
                                      # 'beta_I_alt',
                                      # 'beta_temp', 
                                      'beta_NV'
                                      # 'beta_bush', 
                                      # 'inv_rank', 
                                      )))

sum_nativeV_lizard <- mod_nativeV_lizard$summary()
mod_diagnostics(mod_nativeV_lizard, sum_nativeV_lizard)
ppcheck_nativeV_lizard <- mod_nativeV_lizard$draws('ppcheck', format = 'matrix')

plot(density(dat$Lizard), main = '', 
     xlab = 'Lizard seed dispersal', ylim = c(0, 6))
for (i in 1:200) lines(density(ppcheck_nativeV_lizard[i, ], lwd = 0.1))
lines(density(dat$Lizard), lwd = 2, col = 'red')

post_nativeV_lizard <- 
  mod_nativeV_lizard$draws(c('alpha', 
                           # 'beta_lat', 
                           # 'beta_H_pop', 
                           #'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_nativeV',
                           # 'beta_I_alt', 
                           # 'beta_I_alt',
                           # 'beta_temp', 
                           'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome'), 
                         format = 'df')

post_nativeV_lizard <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_nativeV',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_nativeV_lizard[, grep(x, colnames(post_nativeV_lizard))]
           })

names(post_nativeV_lizard) <- c('alpha', 
                              'beta', 
                              # 'beta_H_pop', 
                              # 'beta_H_foot',
                              # 'beta_I_mainland', 
                              # 'beta_I_nativeV',
                              # 'beta_I_alt', 
                              # 'beta_I_nativeV',
                              # 'beta_temp', 
                              # 'beta_NV',
                              # 'beta_bush', 
                              # 'inv_rank', 
                              'TI', 'p_island', 
                              'p_country', 'p_grid', 
                              'p_plant', 'p_realm', 
                              'p_ecoR', 'p_biome')




# =============== Fruit predation  ======================

file <- paste0(getwd(), '/mod_nativeV_predation.stan')
fit_nativeV_pred <- cmdstan_model(file, compile = T)

mod_nativeV_pred <- 
  fit_nativeV_pred$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_nativeV_pred$save_object('mod_nativeV_pred.rds')

mod_nativeV_pred <- readRDS('mod_nativeV_pred.rds')

sum_nativeV_pred <- mod_nativeV_pred$summary()
mod_diagnostics(mod_nativeV_pred, sum_nativeV_pred)

mcmc_trace(mod_nativeV_pred$draws(c('alpha', 
                                      # 'beta_lat', 
                                      # 'beta_H_pop', 
                                      #'beta_H_foot',
                                      # 'beta_I_mainland', 
                                      # 'beta_I_nativeV',
                                      # 'beta_I_alt', 
                                      # 'beta_I_alt',
                                      # 'beta_temp', 
                                      'beta_NV'
                                      # 'beta_bush', 
                                      # 'inv_rank', 
)))

ppcheck_nativeV_pred <- mod_nativeV_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_nativeV_pred[i, ], lwd = 0.1))
lines(density(dat$predation), lwd = 2, col = 'red')

post_nativeV_pred <- 
  mod_nativeV_pred$draws(c('alpha', 
                             # 'beta_lat', 
                             # 'beta_H_pop', 
                             # 'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_nativeV',
                             # 'beta_I_alt', 
                             #'beta_I_alt',
                             # 'beta_temp', 
                             'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome'), 
                           format = 'df')

post_nativeV_pred <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_nativeV',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           'beta_NV',
           # 'beta_bush', 
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_nativeV_pred[, grep(x, colnames(post_nativeV_pred))]
           })

names(post_nativeV_pred) <- c('alpha', 
                                'beta', 
                                # 'beta_H_pop', 
                                # 'beta_H_foot',
                                # 'beta_I_mainland', 
                                # 'beta_I_nativeV',
                                # 'beta_I_alt', 
                                # 'beta_I_nativeV',
                                # 'beta_temp', 
                                # 'beta_NV',
                                # 'beta_bush', 
                                # 'inv_rank', 
                                'TI', 'p_island', 
                                'p_country', 'p_grid', 
                                'p_plant', 'p_realm', 
                                'p_ecoR', 'p_biome')



# ======= Plots =====
# 
# ============= Slops ===========

# error bars

rbind(pivot_longer(post_nativeV_tot$beta, 'beta_NV') |> 
        mutate(type = 'Frugivory', 
               effect = 'Native cover'), 
      pivot_longer(post_nativeV_bird$beta, 'beta_NV') |> 
        mutate(type = 'Birds', 
               effect = 'Native cover'), 
      pivot_longer(post_nativeV_lizard$beta, 'beta_NV') |> 
        mutate(type = 'Lizards', 
               effect = 'Native cover'), 
      pivot_longer(post_nativeV_disp$beta, 'beta_NV') |> 
        mutate(type = 'Seed dispersion', 
               effect = 'Native cover'), 
      pivot_longer(post_nativeV_pred$beta, 'beta_NV') |> 
        mutate(type = 'Seed predation', 
               effect = 'Native cover')) |> 
  group_by(type) |> 
  transmute(mu = median(value), 
            li = quantile(value, 0.025), 
            ls = quantile(value, 0.975), 
            x = 'Slope', 
            effect = effect) |> 
  unique() |> 
  ggplot(aes(type, mu, ymin = li, ymax = ls)) +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_wrap(~ effect) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs()

# Scatter plots












# =============== *Bush cover* ==========

# =============== Overall frugivory  ======================

file <- paste0(getwd(), '/mod_bush_total.stan')
fit_bush_tot <- cmdstan_model(file, compile = T)

mod_bush_tot <- 
  fit_bush_tot$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_bush_tot$save_object('mod_bush_tot.rds')

mod_bush_tot <- readRDS('mod_bush_tot.rds')

sum_bush_tot <- mod_bush_tot$summary()
mod_diagnostics(mod_bush_tot, sum_bush_tot)
ppcheck_bush_tot <- mod_bush_tot$draws('ppcheck', format = 'matrix')

plot(density(dat$total_remotion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.1))
for (i in 1:200) lines(density(ppcheck_bush_tot[i, ], lwd = 0.1))
lines(density(dat$total_remotion), lwd = 2, col = 'red')


post_bush_tot <- 
  mod_bush_tot$draws(c('alpha', 
                          #'beta_lat', 
                          # 'beta_H_pop', 
                          # 'beta_H_foot',
                          # 'beta_I_mainland', 
                          #'beta_I_bush',
                          # 'beta_I_alt', 
                          #'beta_I_alt',
                          # 'beta_temp', 
                          # 'beta_NV',
                          'beta_bush',
                          # 'inv_rank', 
                          'TI', 'p_island', 
                          'p_country', 'p_grid', 
                          'p_plant', 'p_realm', 
                          'p_ecoR', 'p_biome'), 
                        format = 'df')

post_bush_tot <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           #'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_bush',
           # 'beta_I_alt', 
           #'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           'beta_bush',
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_bush_tot[, grep(x, colnames(post_bush_tot))]
           })

names(post_bush_tot) <- c('alpha', 
                             'beta', 
                             # 'beta_H_pop', 
                             # 'beta_H_foot',
                             # 'beta_I_mainland', 
                             # 'beta_I_bush',
                             # 'beta_I_alt', 
                             #'beta_I_bush',
                             # 'beta_temp', 
                             # 'beta_NV',
                             # 'beta_bush', 
                             # 'inv_rank', 
                             'TI', 'p_island', 
                             'p_country', 'p_grid', 
                             'p_plant', 'p_realm', 
                             'p_ecoR', 'p_biome')



# =============== Fruit dispersion  ======================

file <- paste0(getwd(), '/mod_bush_dispersion.stan')
fit_bush_disp <- cmdstan_model(file, compile = T)

mod_bush_disp <- 
  fit_bush_disp$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_bush_disp$save_object('mod_bush_dip.rds')

mod_bush_disp <- readRDS('mod_bush_dip.rds')

sum_bush_disp <- mod_bush_disp$summary()
mod_diagnostics(mod_bush_disp, sum_bush_disp)
ppcheck_bush_disp <- mod_bush_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_bush_disp[i, ], lwd = 0.1))
lines(density(dat$dispersion), lwd = 2, col = 'red')

post_bush_disp <- 
  mod_bush_disp$draws(c('alpha', 
                           # 'beta_lat', 
                           # 'beta_H_pop', 
                           #'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_bush',
                           # 'beta_I_alt', 
                           # 'beta_I_alt',
                           # 'beta_temp', 
                           # 'beta_NV',
                           'beta_bush',
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome'), 
                         format = 'df')

post_bush_disp <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_bush',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           'beta_bush',
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_bush_disp[, grep(x, colnames(post_bush_disp))]
           })

names(post_bush_disp) <- c('alpha', 
                              'beta', 
                              # 'beta_H_pop', 
                              # 'beta_H_foot',
                              # 'beta_I_mainland', 
                              # 'beta_I_bush',
                              # 'beta_I_alt', 
                              # 'beta_I_bush',
                              # 'beta_temp', 
                              # 'beta_NV',
                              # 'beta_bush', 
                              # 'inv_rank', 
                              'TI', 'p_island', 
                              'p_country', 'p_grid', 
                              'p_plant', 'p_realm', 
                              'p_ecoR', 'p_biome')


# =============== Birds  ======================

file <- paste0(getwd(), '/mod_bush_bird.stan')
fit_bush_bird <- cmdstan_model(file, compile = T)

mod_bush_bird <- 
  fit_bush_bird$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

mod_bush_bird$save_object('mod_bush_bird.rds')

mod_bush_bird <- readRDS('mod_bush_bird.rds')

sum_bush_bird <- mod_bush_bird$summary()
mod_diagnostics(mod_bush_bird, sum_bush_bird)
ppcheck_bush_bird <- mod_bush_bird$draws('ppcheck', format = 'matrix')

plot(density(dat$Bird), main = '', 
     xlab = 'Bird seed dispersal', ylim = c(0, 3))
for (i in 1:200) lines(density(ppcheck_bush_bird[i, ], lwd = 0.1))
lines(density(dat$Bird), lwd = 2, col = 'red')

post_bush_bird <- 
  mod_bush_bird$draws(c('alpha', 
                        # 'beta_lat', 
                        # 'beta_H_pop', 
                        #'beta_H_foot',
                        # 'beta_I_mainland', 
                        # 'beta_I_bush',
                        # 'beta_I_alt', 
                        # 'beta_I_alt',
                        # 'beta_temp', 
                        # 'beta_NV',
                        'beta_bush',
                        # 'inv_rank', 
                        'TI', 'p_island', 
                        'p_country', 'p_grid', 
                        'p_plant', 'p_realm', 
                        'p_ecoR', 'p_biome'), 
                      format = 'df')

post_bush_bird <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_bush',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           'beta_bush',
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_bush_bird[, grep(x, colnames(post_bush_bird))]
           })

names(post_bush_bird) <- c('alpha', 
                           'beta', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_bush',
                           # 'beta_I_alt', 
                           # 'beta_I_bush',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome')



# =============== Lizards  ======================

file <- paste0(getwd(), '/mod_bush_lizard.stan')
fit_bush_lizard <- cmdstan_model(file, compile = T)

mod_bush_lizard <- 
  fit_bush_lizard$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )


mod_bush_lizard$save_object('mod_bush_lizard.rds')
mod_bush_lizard <- readRDS('mod_bush_lizard.rds')

sum_bush_lizard <- mod_bush_lizard$summary()
mod_diagnostics(mod_bush_lizard, sum_bush_lizard)
ppcheck_bush_lizard <- mod_bush_lizard$draws('ppcheck', format = 'matrix')

plot(density(dat$Lizard), main = '', 
     xlab = 'Lizard seed dispersion', ylim = c(0, 6))
for (i in 1:200) lines(density(ppcheck_bush_lizard[i, ], lwd = 0.1))
lines(density(dat$Lizard), lwd = 2, col = 'red')

post_bush_lizard <- 
  mod_bush_lizard$draws(c('alpha', 
                        # 'beta_lat', 
                        # 'beta_H_pop', 
                        #'beta_H_foot',
                        # 'beta_I_mainland', 
                        # 'beta_I_bush',
                        # 'beta_I_alt', 
                        # 'beta_I_alt',
                        # 'beta_temp', 
                        # 'beta_NV',
                        'beta_bush',
                        # 'inv_rank', 
                        'TI', 'p_island', 
                        'p_country', 'p_grid', 
                        'p_plant', 'p_realm', 
                        'p_ecoR', 'p_biome'), 
                      format = 'df')

post_bush_lizard <- 
  lapply(c('alpha', 
           #'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_bush',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           'beta_bush',
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_bush_lizard[, grep(x, colnames(post_bush_lizard))]
           })

names(post_bush_lizard) <- c('alpha', 
                           'beta', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_bush',
                           # 'beta_I_alt', 
                           # 'beta_I_bush',
                           # 'beta_temp', 
                           # 'beta_NV',
                           # 'beta_bush', 
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome')






# =============== Fruit predation  ======================

file <- paste0(getwd(), '/mod_bush_predation.stan')
fit_bush_pred <- cmdstan_model(file, compile = T)

mod_bush_pred <- 
  fit_bush_pred$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

# mod_bush_pred$save_object('mod_bush_pred.rds')
mod_bush_pred <- readRDS('mod_bush_pred.rds')

sum_bush_pred <- mod_bush_pred$summary()
mod_diagnostics(mod_bush_pred, sum_bush_pred)

ppcheck_bush_pred <- mod_bush_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
for (i in 1:200) lines(density(ppcheck_bush_pred[i, ], lwd = 0.1))
lines(density(dat$predation), lwd = 2, col = 'red')

post_bush_pred <- 
  mod_bush_pred$draws(c('alpha', 
                           # 'beta_lat', 
                           # 'beta_H_pop', 
                           # 'beta_H_foot',
                           # 'beta_I_mainland', 
                           # 'beta_I_bush',
                           # 'beta_I_alt', 
                           #'beta_I_alt',
                           # 'beta_temp', 
                           # 'beta_NV',
                           'beta_bush',
                           # 'inv_rank', 
                           'TI', 'p_island', 
                           'p_country', 'p_grid', 
                           'p_plant', 'p_realm', 
                           'p_ecoR', 'p_biome'), 
                         format = 'df')

post_bush_pred <- 
  lapply(c('alpha', 
           # 'beta_lat', 
           # 'beta_H_pop', 
           # 'beta_H_foot',
           # 'beta_I_mainland', 
           # 'beta_I_bush',
           # 'beta_I_alt', 
           # 'beta_I_alt',
           # 'beta_temp', 
           # 'beta_NV',
           'beta_bush',
           # 'inv_rank', 
           'TI', 'p_island', 
           'p_country', 'p_grid', 
           'p_plant', 'p_realm', 
           'p_ecoR', 'p_biome'), FUN = 
           function(x) {
             post_bush_pred[, grep(x, colnames(post_bush_pred))]
           })

names(post_bush_pred) <- c('alpha', 
                              'beta', 
                              # 'beta_H_pop', 
                              # 'beta_H_foot',
                              # 'beta_I_mainland', 
                              # 'beta_I_bush',
                              # 'beta_I_alt', 
                              # 'beta_I_bush',
                              # 'beta_temp', 
                              # 'beta_NV',
                              # 'beta_bush', 
                              # 'inv_rank', 
                              'TI', 'p_island', 
                              'p_country', 'p_grid', 
                              'p_plant', 'p_realm', 
                              'p_ecoR', 'p_biome')


# ======= Plots =====

rm(list = ls()[grep('ppcheck', ls())])

# error bars (island level)
average_effects
codes$island 

islands_post <- 
  lapply(ls()[grep('^post_latit', ls())], FUN = 
         function(x) {
           
           posterior <- get(x)
           
           est1 <- 
             lapply(codes$island$code, FUN =
                      function(i) {
                        
                        indx_xvar <- which(dat$island == i)
                        indx_realm <- unique(dat$realm[indx_xvar])
                        indx_count <- unique(dat$country[indx_xvar])
                        indx_typeI <- unique(dat$island_type[indx_xvar])
                        indx_grid <- unique(dat$grid[indx_xvar])
                        indx_plant <- unique(dat$plant_ID[indx_xvar])
                        indx_ecoR <- unique(dat$ecoregion[indx_xvar])
                        indx_biome <- unique(dat$biome[indx_xvar])
                        
                        island <- codes$island
                        realm <- codes$real
                        ty_is <- codes$island_type
                        
                        lab_island <- island[island$code == i, ]$island
                        lab_realm <- realm[realm$code == indx_realm, ]$island
                        lab_TI <- ty_is[ty_is$code == indx_typeI, ]$island
                        
                        est <-
                          with(posterior,
                               {
                                 inv_logit(alpha$alpha +
                                             p_island[, i, drop = T] +
                                             TI[, indx_typeI, drop = T] +
                                             p_realm[, indx_realm, drop = T] +
                                             p_country[, indx_count, drop = T] +
                                             apply(p_grid[, indx_grid], 1, median) +
                                             apply(p_plant[, indx_plant], 1, median) +
                                             apply(p_ecoR[, indx_ecoR], 1, median) +
                                             apply(p_biome[, indx_biome], 1, median)
                                           )
                               })
                        
                        print(i)
                        
                        tibble(mu = median(est),
                               li = quantile(est, 0.025),
                               ls = quantile(est, 0.975),
                               island = lab_island,
                               type_island = lab_TI,
                               realm = lab_realm,
                               modelo = gsub('^(post_latitude_)(.*)$',
                                             '\\2', x))

                      })
           
           do.call('rbind', est1)

         })

names(islands_post) <- ls()[grep('^post_latit', ls())]

islands_post <- 
  lapply(seq_along(islands_post), FUN = 
         function(x) {
           d <- islands_post[[x]]
           d$model <- c('Bird frugivory', 
                        'Seed dispersal', 'Lizard frugivory', 
                        'Seed predation', 
                         'Frugivory')[[x]]
           d
         })

names(islands_post) <- ls()[grep('^post_latit', ls())]


plots_islands_realm <- 
  lapply(islands_post, FUN = 
           function(x) {
             ggplot(data = x, 
                    aes(fct_reorder(island, mu, .fun = 'max'), mu, 
                        ymin = li, ymax = ls, 
                        color = realm)) +
               labs(y = 'P(fruit consumption)', x = 'Islands') +
               geom_hline(yintercept = 0.5, linetype = 3) +
               geom_errorbar(width = 0, linewidth = 1, alpha = 0.6) +
               scale_color_manual(values = RColorBrewer::brewer.pal(7, 'Set2')) +
               geom_point(size = 0.5) +
               theme_classic() +
               facet_wrap(~model) +
               theme(axis.text.x = element_blank(), 
                     legend.title = element_blank(), 
                     legend.position = c(0.15, 0.8), 
                     legend.box.background = element_blank(), 
                     legend.key.size = unit(4, 'mm'),
                     text = element_text(family = 'Times New Roman'),
                     strip.background = element_blank(), 
                     axis.line = element_line(linewidth = 0.25), 
                     axis.ticks = element_line(linewidth = 0.25))
           })


plots_islands_realm$post_latitude_tot <- 
  plots_islands_realm$post_latitude_tot +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'none')

plots_islands_realm$post_latitude_disp <- 
  plots_islands_realm$post_latitude_disp +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'none')

plots_islands_realm$post_latitude_pred <- 
  plots_islands_realm$post_latitude_pred +
  labs(y = NULL) +
  theme(legend.position = 'none')

plots_islands_realm$post_latitude_bird <- 
  plots_islands_realm$post_latitude_bird +
  theme(legend.position = 'none')

plots_islands_realm$post_latitude_lizard <- 
  plots_islands_realm$post_latitude_lizard +
  labs(y = NULL) +
  theme(legend.position = 'none')

# plots_islands_realm$post_latitude_tot <- 
#   ggplot(data = islands_post$post_latitude_tot, 
#          aes(fct_reorder(island, mu, .fun = 'max'), mu, 
#              ymin = li, ymax = ls, 
#              color = realm)) +
#   labs(y = 'P(fruit consumption)', x = NULL) +
#   geom_hline(yintercept = 0.5, linetype = 3) +
#   geom_errorbar(width = 0, linewidth = 1.5, alpha = 0.6) +
#   scale_color_manual(values = RColorBrewer::brewer.pal(7, 'Set2')) +
#   geom_point(size = 1.5) +
#   theme_classic() +
#   facet_wrap(~model) +
#   theme(axis.text.x = element_blank(), legend.title = element_blank(), 
#         legend.position = c(0.15, 0.85), 
#         legend.box.background = element_blank(), 
#         legend.key.size = unit(4, 'mm'),
#         text = element_text(family = 'Times New Roman'),
#         strip.background = element_blank(), 
#         axis.line = element_line(linewidth = 0.25), 
#         axis.ticks = element_line(linewidth = 0.25))



layout <- 
  '
  aab
  aac
  de#
'

plots_islands_realm$post_latitude_tot + 
  plots_islands_realm$post_latitude_disp +
  plots_islands_realm$post_latitude_pred +
  plots_islands_realm$post_latitude_bird +
  plots_islands_realm$post_latitude_lizard +
  plot_layout(design = layout, labels) +
  plot_annotation(tag_levels = 'a', 
                  tag_prefix = '(', 
                  tag_suffix = ')')

ggsave('Figure_2.jpeg', width = 18, height = 18, units = 'cm', 
       dpi = 500)


plots_islands_type <- 
  lapply(islands_post, FUN = 
           function(x) {
             ggplot(data = x, 
                    aes(fct_reorder(island, mu, .fun = 'max'), mu, 
                        ymin = li, ymax = ls, 
                        color = type_island)) +
               labs(y = 'P(fruit consumption)', x = 'Islands') +
               geom_hline(yintercept = 0.5, linetype = 3) +
               geom_errorbar(width = 0, linewidth = 1, alpha = 0.6) +
               scale_color_manual(values = 
                                    RColorBrewer::brewer.pal(3, "Accent")) +
               geom_point(size = 0.5) +
               theme_classic() +
               facet_wrap(~model) +
               theme(axis.text.x = element_blank(), 
                     legend.title = element_blank(), 
                     legend.position = c(0.15, 0.8), 
                     legend.box.background = element_blank(), 
                     legend.key.size = unit(4, 'mm'),
                     text = element_text(family = 'Times New Roman'),
                     strip.background = element_blank(), 
                     axis.line = element_line(linewidth = 0.25), 
                     axis.ticks = element_line(linewidth = 0.25))
           })


plots_islands_type$post_latitude_disp <- 
  plots_islands_type$post_latitude_disp +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank())

plots_islands_type$post_latitude_pred <- 
  plots_islands_type$post_latitude_pred +
  labs(y = NULL) +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank())

plots_islands_type$post_latitude_bird <- 
  plots_islands_type$post_latitude_bird +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank())

plots_islands_type$post_latitude_lizard <- 
  plots_islands_type$post_latitude_lizard +
  labs(y = NULL) +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank())

plots_islands_type$post_latitude_tot <- 
  ggplot(data = islands_post$post_latitude_tot, 
         aes(fct_reorder(island, mu, .fun = 'max'), mu, 
             ymin = li, ymax = ls, 
             color = type_island)) +
  labs(y = 'P(fruit consumption)', x = NULL) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  geom_errorbar(width = 0, linewidth = 1.5, alpha = 0.6) +
  scale_color_manual(values = RColorBrewer::brewer.pal(7, 'Set2')) +
  geom_point(size = 1.5) +
  theme_classic() +
  facet_wrap(~model) +
  theme(axis.text.x = element_blank(), legend.title = element_blank(), 
        legend.position = c(0.15, 0.9), 
        legend.box.background = element_blank(), 
        legend.key.size = unit(4, 'mm'),
        text = element_text(family = 'Times New Roman'),
        strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.25), 
        axis.ticks = element_line(linewidth = 0.25))


plots_islands_type$post_latitude_tot + 
  plots_islands_type$post_latitude_disp +
  plots_islands_type$post_latitude_pred +
  plots_islands_type$post_latitude_bird +
  plots_islands_type$post_latitude_lizard +
  plot_layout(design = layout, labels) +
  plot_annotation(tag_levels = 'a', 
                  tag_prefix = '(', 
                  tag_suffix = ')')

ggsave('figure_3.jpeg', width = 18, height = 18, units = 'cm', 
       dpi = 500)


# ============= Slops ===========

# error bars

rbind(pivot_longer(post_bush_tot$beta, 'beta_bush') |> 
        mutate(type = 'Total frugivory', 
               effect = 'Bush cover'),
      pivot_longer(post_bush_bird$beta, 'beta_bush') |> 
        mutate(type = 'Bird frugivory', 
               effect = 'Bush cover'),
      pivot_longer(post_bush_lizard$beta, 'beta_bush') |> 
        mutate(type = 'Lizard frugivory', 
               effect = 'Bush cover'),
      pivot_longer(post_bush_disp$beta, 'beta_bush') |> 
        mutate(type = 'Seed dispersal', 
               effect = 'Bush cover'), 
      pivot_longer(post_bush_pred$beta, 'beta_bush') |> 
        mutate(type = 'Seed predation', 
               effect = 'Bush cover')) |> 
  group_by(type) |> 
  transmute(mu = median(value), 
            li = quantile(value, 0.025), 
            ls = quantile(value, 0.975), 
            x = 'Slope', 
            effect = effect) |> 
  unique() |> 
  ggplot(aes(type, mu, ymin = li, ymax = ls)) +
  geom_point() +
  geom_errorbar(width = 0) +
  facet_wrap(~ effect) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs()

# Scatter plots


mean(post_bush_pred$beta$beta_bush < 0)



# ====== *Final plots and estimations* ====

# categorical effects =================

# Type of island
plot_disp_pred_islands <- 
  rbind(full_join(TI_disp, codes$island_type, 'code') |> 
        mutate(type = 'Seed dispersal'),
        full_join(TI_bird, codes$island_type, 'code') |> 
          mutate(type = "Bird frugivory"),
        full_join(TI_lizard, codes$island_type, 'code') |> 
          mutate(type = "Lizard frugivory"),
      full_join(TI_pred, codes$island_type, 'code') |> 
        mutate(type = 'Seed predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> 
  ggplot(aes(island, mu, ymin = li, ymax = ls, color = type)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.4), 
                linewidth = 1.5, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  scale_color_manual(values = c('#F2C230', '#7A577A', '#8C1F28', '#0897B4')) +
  scale_x_discrete(labels = c('Continental', 'Coralline', 'Volcanic')) +
  labs(y = ' ', x = 'Type of island') +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.background = element_blank(),
    legend.position = c(0.24, 0.8),
    legend.box.background = element_blank(), 
    legend.key.size = unit(3.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )

plot_frugivory_islands <- 
  full_join(TI_tot, codes$island_type, 'code') |> 
  mutate(type = 'Frugivory') |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> 
  ggplot(aes(island, mu, ymin = li, ymax = ls, color = type)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.2), 
                linewidth = 1.5, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.2), size = 2) +
  scale_color_manual(values = c('#D92525')) +
  labs(y = 'P(fruit consumption)', x = 'Type of island') +
  scale_x_discrete(labels = c('Continental', 'Coralline', 'Volcanic')) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.background = element_blank(),
    legend.position = c(0.24, 0.8),
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )


plot_frugivory_realms <- 
full_join(realm_tot, codes$real, 'code') |> 
  mutate(type = 'Total frugivory') |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> 
  ggplot(aes(island, mu, ymin = li, ymax = ls, color = type)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.4), 
                linewidth = 1.5, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  scale_color_manual(values = c('#D92525')) +
  labs(y = 'P(fruit consumption)', x = 'Biogeographic realm') +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


plot_disp_pred_realms <- 
rbind(full_join(realm_disp, codes$real, 'code') |> 
        mutate(type = 'Seed dispersal'),
      full_join(realm_bird, codes$real, 'code') |> 
        mutate(type = 'Bird frugivory'),
      full_join(realm_lizard, codes$real, 'code') |> 
        mutate(type = 'Lizard frugivory'),
      full_join(realm_pred, codes$real, 'code') |> 
        mutate(type = 'Seed predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> 
  ggplot(aes(island, mu, ymin = li, ymax = ls, color = type)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.5), 
                linewidth = 1.5, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  scale_color_manual(values = c('#F2C230', '#7A577A', '#8C1F28', '#0897B4')) +
  labs(y = ' ', x = 'Biogeographic realm') +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


plot_grid(plot_grid(NULL, 
                    plot_frugivory_islands, 
                    plot_disp_pred_islands, 
                    NULL,
                    nrow = 1,  
                    rel_widths = c(0.2, 0.8, 0.8, 0.2), 
                    labels = c('', '(a)', '(b)', ''), 
                    label_size = 9, 
                    label_fontfamily = 'Times New Roman', 
                    label_fontface = 'plain', 
                    label_x = 0.15, 
                    label_y = 0.99), 
          plot_grid(plot_frugivory_realms, 
                    plot_disp_pred_realms, 
                    nrow = 1, 
                    labels = c('(c)', '(d)'), 
                    label_size = 9, 
                    label_fontfamily = 'Times New Roman', 
                    label_fontface = 'plain', 
                    label_x = c(0.12, 0.135), 
                    label_y = 0.99),
          ncol = 1)

ggsave('plot_type_island_realm.jpg', dpi = 1e3, width = 15, height = 10, units = 'cm')


avg_realms_est <- 
  rbind(full_join(realm_tot, codes$real, 'code') |> 
        mutate(type = 'Frugivory'), 
      full_join(realm_disp, codes$real, 'code') |> 
        mutate(type = 'Seed\ndispersal'),
      full_join(realm_bird, codes$real, 'code') |> 
        mutate(type = 'Bird\nfrugivory'),
      full_join(realm_lizard, codes$real, 'code') |> 
        mutate(type = 'Lizard\nfrugivory'),
      full_join(realm_pred, codes$real, 'code') |> 
        mutate(type = 'Seed\npredation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique()

avg_realms_est <- 
  lapply(split(avg_realms_est, avg_realms_est$type), 
       FUN = 
         function(x) {
           x <- x[order(x$mu, decreasing = F), ]
           labs <- x$island
           x$island <- factor(x$island, 
                               levels = labs)
           x
         }) 

avg_realms_est <- do.call('rbind', avg_realms_est)

levels(avg_realms_est$island)



avg_realms_plot <- 
  avg_realms_est |> 
  ggplot(aes(fct_reorder(type, mu, .desc = T), 
             mu, ymin = li, ymax = ls, color = island)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.8), 
                linewidth = 1.5, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.8), size = 2) +
  scale_color_manual(values = c("#A6D854", "#66C2A5", "#FC8D62", 
                                "#E5C494", "#E78AC3", "#FFD92F", "#8DA0CB")) +
  labs(y = 'P(fruit consumption)', x = '') +
  theme_classic() +
  lims(y = c(0, 1)) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.position = c(0.9, 0.85), 
    legend.background = element_blank(),
    legend.box.background = element_blank(), 
    legend.key.size = unit(4, 'mm'), 
    text = element_text(family = 'Times New Roman')
  )



layout <- 
  '
  aab
  aac
  dfg
'

avg_realms_plot +
  plots_islands_realm$post_latitude_tot +
  plots_islands_realm$post_latitude_disp +
  plots_islands_realm$post_latitude_pred + 
  labs(y = 'P(fruit consumption)', x = 'Islands') +
  plots_islands_realm$post_latitude_bird +
  labs(y = '') +
  plots_islands_realm$post_latitude_lizard +
  plot_layout(design = layout, 
              widths = c(1, 1, 1, 1, 1, 0.5)) +
  plot_annotation(tag_levels = 'a', 
                  tag_prefix = '(', 
                  tag_suffix = ')')

ggsave('figure_2.jpeg', width = 20, height = 18, units = 'cm', 
       dpi = 500)



# estimations

# Island coefficients

island_tables <- do.call('rbind', islands_post)

island_tables1 <- split(island_tables, list(island_tables$type_island, 
                                            island_tables$model))

lapply(island_tables1, FUN = 
         function(x) {
           x[order(x$mu, decreasing = T), ]
         })

island_tables2 <- split(island_tables, 
                        list(island_tables$realm, 
                             island_tables$model))

island_tables2 <- 
  lapply(island_tables2, FUN = 
         function(x) {
           x[order(x$mu, decreasing = T), ]
         })

island_tables2$`Neotropical.Total frugivory` |> print(n = 100)

island_tables2$Australasia$island
# probability among types of islands
rbind(full_join(TI_tot, codes$island_type, 'code') |> 
        mutate(type = 'Total frugivory'),
      full_join(TI_disp, codes$island_type, 'code') |> 
        mutate(type = 'Seed dispersal'),
      full_join(TI_bird, codes$island_type, 'code') |> 
        mutate(type = "Bird frugivory"),
      full_join(TI_lizard, codes$island_type, 'code') |> 
        mutate(type = "Lizard frugivory"),
      full_join(TI_pred, codes$island_type, 'code') |> 
        mutate(type = 'Seed predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() 

lapply(1, FUN = 
         function(x) {
           
           d <- list(TI_tot, TI_disp, TI_pred, TI_bird, TI_lizard)
           d <- lapply(d, FUN = 
                         function(z) {
                           d <- lapply(1:3, FUN = 
                                         function(k) {
                                           z[z$code == k, 'y']
                                         })
                           d <- do.call('cbind', d) 
                           colnames(d) <- codes$island_type$island
                           d
                         })
          
           names(d) <- c('Total frugivory', 
                         'Seed dispersal', 'Seed predation', 
                         'Bird frugivory', 'Lizard frugivory')
           
           d <- 
             lapply(d, FUN = 
                      function(jj) {
                        
                        prob <-
                          lapply(1:(3-1), FUN =
                                   function(i) {
                                     df <-
                                       lapply((i+1):3, FUN =
                                                function(j) {
                                                  p <- mean(jj[[i]] > jj[[j]])
                                                  dif <- jj[[i]] - jj[[j]]
                                                  contrast <-
                                                    paste0('P(', 
                                                          codes$island_type$island[i],
                                                          ' > ',
                                                          codes$island_type$island[j], 
                                                          ')')
                                                  
                                                  tibble(Comparison = contrast,
                                                         probability = p,
                                                         contrast = mean(dif),
                                                         li = quantile(dif, 0.025),
                                                         ls = quantile(dif, 0.975))
                                                  
                                                })
                                     do.call('rbind', df)
                                   })
                        
                        do.call('rbind', prob)
                        
                      })
           
           for (i in seq_along(d)) {
             d[[i]]$type <- names(d)[i]
           }
          
           d <- do.call('rbind', d)
           
           
           d
         })[[1]] |> print(n = 29) 


# frugivory among realms

# average probabilities 
rbind(full_join(realm_tot, codes$real, 'code') |> 
        mutate(type = 'Total frugivory'),
      full_join(realm_disp, codes$real, 'code') |> 
        mutate(type = 'Seed dispersal'),
      full_join(realm_bird, codes$real, 'code') |> 
        mutate(type = 'Bird frugivory'),
      full_join(realm_lizard, codes$real, 'code') |> 
        mutate(type = 'Lizard frugivory'),
      full_join(realm_pred, codes$real, 'code') |> 
        mutate(type = 'Seed predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> print(n = 40)


# contrast among realms 

lapply(1, FUN = 
         function(x) {
           
           d <- list(realm_tot, realm_disp, realm_pred, 
                     realm_bird, realm_lizard)
           d <- lapply(d, FUN = 
                         function(z) {
                           d <- lapply(codes$real$code, FUN = 
                                         function(k) {
                                           z[z$code == k, 'y']
                                         })
                           d <- do.call('cbind', d) 
                           colnames(d) <- codes$real$island
                           d
                         })
           
           names(d) <- c('Total frugivory', 
                         'Seed dispersal', 'Seed predation', 
                         'Bird frugivory', 'Lizard frugivory')
           
           lapply(d, FUN = 
                    function(z) {
                      
                      m <- 
                      sapply(codes$real$code, FUN = 
                               function(i) {
                                 sapply(codes$real$code, FUN = 
                                          function(j) {
                                            mean(z[[i]] > z[[j]])
                                          })
                               })
                      
                      dimnames(m) <- list(codes$real$island, 
                                          codes$real$island)
                      m[upper.tri(m)] <- NA
                      m
                    })
         })

cont_realms_islands <- 
  lapply(1, FUN = 
         function(x) {
           
           d <- list(realm_tot, realm_disp, realm_pred,
                     realm_bird, realm_lizard)
           d <- lapply(d, FUN = 
                         function(z) {
                           d <- lapply(codes$real$code, FUN = 
                                         function(k) {
                                           z[z$code == k, 'y']
                                         })
                           d <- do.call('cbind', d) 
                           colnames(d) <- codes$real$island
                           d
                         })
           
           names(d) <- c('Total frugivory', 'Seed dispersal',
                         'Seed predation', 
                         'Bird frugivory', 'Lizard frugivory')
           
           d <-
             lapply(d, FUN =
                      function(jj) {

                        prob <-
                          lapply(1:(length(codes$real$code)-1), FUN =
                                   function(i) {
                                     df <-
                                       lapply((i+1):length(codes$real$code), 
                                              FUN =
                                                function(j) {
                                                  p <- mean(jj[[i]] > jj[[j]])
                                                  dif <- jj[[i]] - jj[[j]]
                                                  contrast <-
                                                    paste('P(', 
                                                          codes$real$island[i],
                                                          ' > ',
                                                          codes$real$island[j],
                                                          ')',
                                                          sep = '')

                                                  tibble(Comparison = contrast,
                                                         probability = p,
                                                         contrast = mean(dif),
                                                         li = quantile(dif, 0.025),
                                                         ls = quantile(dif, 0.975))

                                                })
                                     do.call('rbind', df)
                                   })

                        do.call('rbind', prob)

                      })

           for (i in seq_along(d)) {
             d[[i]]$type <- names(d)[i]
           }

           d <- do.call('rbind', d)

           d[d$probability != 0, ]
           d
         })[[1]] 


lapply(split(cont_realms_islands, 
             cont_realms_islands$type), FUN = 
         function(x) {
           lapply(codes$real$island[c(5)], FUN = 
                    function(j) {
                      print(j)
                      x <- x[grep(j, x$Comparison), ]
                      x[order(x$contrast, decreasing = T), ]
                    })
         })$`Bird frugivory`



cont_realms_islands
# Continuous variabes ====================

# latitude (mu and CI)

all_betas <- 
  rbind(slope_pars('post_latitude', 'beta_lat') |> 
        mutate(predictor = 'Latitude'), 
      slope_pars('post_isolation', 'beta_I_isolation') |> 
        mutate(predictor = 'Island isolation'), 
      slope_pars('post_size', 'beta_I_size') |> 
        mutate(predictor = 'Island area'), 
      slope_pars('post_altitude', 'beta_I_alt') |> 
        mutate(predictor = 'Island altitude'), 
      slope_pars('post_footprint', 'beta_H_foot') |> 
        mutate(predictor = 'Human footprint'), 
      slope_pars('post_nativeV', 'beta_NV') |> 
        mutate(predictor = 'Native vegetation'), 
      slope_pars('post_bush', 'beta_bush') |> 
        mutate(predictor = 'Bush cover'))

all_betas |> print(n = 35)

# latitude

est_latitude_tot <- cond_effects(posterior = post_latitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'averaging',
                                 n = 100)

causal_effect_latitude <-  # causal effect of latitude
  lapply(1:4, FUN = 
         function(x) {
           
           #interventions in SD
           i <- c(-1, 0, 1, quantile(dat$lat, 0.1))[x] 
           j <- c(0, 1, 2, quantile(dat$lat, 0.9))[x]
           
           d <-
             cond_effects(posterior = post_latitude_tot,
                          x_bar = dat$lat,
                          slope = 'beta_lat',
                          type = 'averaging',
                          causal_effect = T,
                          intervention1 = i,
                          intervention2 = j)
           
           d[, 3]
           
         })

causal_effect_latitude <- as_tibble(do.call('cbind', causal_effect_latitude))

colnames(causal_effect_latitude)[4] <- paste('10th vs. 90th', ' percentile', 
                                             sep = '\n')

causal_effect_latitude <- 
  pivot_longer(causal_effect_latitude, colnames(causal_effect_latitude), 
               names_to = 'intervention', values_to = 'probability')

causal_effect_latitude$frugivory_type <- 'Total frugivory'
causal_effect_latitude$factor <- 'Latitude'
#"#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" 
# lizard  predation. total      bird.    'disp

beta_latitude <- 
  post_latitude_tot$beta |> 
  mutate(frugivory_type = 'Total frugivory') |> 
  rename(x = beta_lat) |> 
  ggplot(aes(x)) +
  geom_density(fill = "#4DAF4A", alpha = 0.45, 
               linewidth = 0) +
  labs(x = expression(beta['latitude']), 
       y = 'Density') +
  geom_vline(xintercept = 0, linetype = 3) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35), 
        axis.title.x = element_text(size = 15), 
        text = element_text(family = 'Times New Roman'))

plot_causal_effect_latitude <- 
  causal_effect_latitude |> 
  group_by(intervention) |> 
  transmute(mu = median(probability), 
            li = quantile(probability, 0.025), 
            ls = quantile(probability, 0.975), 
            var = 'Latitude') |> 
  unique() |> 
  ggplot(aes(intervention, mu, ymin = li, ymax = ls)) +
  geom_errorbar(width = 0, linewidth = 2, 
                color = "#E41A1C", alpha = 0.5) +
  geom_point(size = 1.5, color = "#E41A1C") +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(~var) +
  labs(x = 'z-score difference', 
       y = expression(Delta ~'P(Fruit consumption)')) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35),
        text = element_text(family = 'Times New Roman'))
  

# isolation

names(post_isolation_tot)[2] <- 'beta'

est_isolation_tot <- cond_effects(posterior = post_isolation_tot,
                                  x_bar = dat$isolation,
                                  slope = 'beta_I_isolation',
                                  type = 'averaging',
                                  n = 100)


est_isolation_pred <- cond_effects(posterior = post_isolation_pred,
                                   x_bar = dat$isolation,
                                   slope = 'beta_I_isolation',
                                   type = 'averaging',
                                   n = 100)

est_isolation_lizard <- cond_effects(posterior = post_isolation_lizard,
                                   x_bar = dat$isolation,
                                   slope = 'beta_I_isolation',
                                   type = 'averaging',
                                   n = 100)

causal_effect_isolation <- 
  lapply(1:3, FUN = 
         function(x) {
           
           model <- list(post_isolation_tot, 
                         post_isolation_pred, 
                         post_isolation_lizard)[[x]]
           model_lab <- c('Total frugivory', 
                          'Seed predation', 
                          'Lizard frugivory')[x]
           
           d <-  # causal effect of latitude
             lapply(1:4, FUN = 
                      function(x) {
                        
                        #interventions in SD
                        i <- c(-1, 0, 1, quantile(dat$isolation, 0.1))[x] 
                        j <- c(0, 1, 2, quantile(dat$isolation, 0.9))[x]
                        
                        d <-
                          cond_effects(posterior = model,
                                       x_bar = dat$isolation,
                                       slope = 'beta_I_isolation',
                                       type = 'averaging',
                                       causal_effect = T,
                                       intervention1 = i,
                                       intervention2 = j)
                        
                        d[, 3]
                        
                      })
           
           d <- as_tibble(do.call('cbind', d))
           
           colnames(d)[4] <- paste('10th vs. 90th', ' percentile', 
                                   sep = '\n')
           
           d <- 
             pivot_longer(d, colnames(d), 
                          names_to = 'intervention', values_to = 'probability')
           
           d$frugivory_type <- model_lab
           d$factor <- 'Island isolation'
           d
         })

causal_effect_isolation <- do.call('rbind', causal_effect_isolation)

tribble(~code_color, ~group, 
        "#E41A1C", 'lizard',
        "#377EB8", 'seed predation',
        "#4DAF4A", 'total frugivory',
        "#984EA3", 'bird frugivory',
        "#FF7F00", 'seed dispersal')

beta_isolation <- 
rbind(post_isolation_tot$beta |> 
        mutate(frugivory_type = 'Total frugivory') |> 
        rename(x = beta_I_isolation), 
      post_isolation_pred$beta |> 
        mutate(frugivory_type = 'Seed predation') |> 
        rename(x = beta_I_isolation), 
      post_isolation_lizard$beta |> 
        mutate(frugivory_type = 'Lizard frugivory') |> 
        rename(x = beta_I_isolation)) |> 
  ggplot(aes(x, fill = frugivory_type)) +
  geom_density(linewidth = 0, alpha = 0.4) +
  labs(x = expression(beta['island isolation']), 
       y = 'Density') +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35), 
        axis.title.x = element_text(size = 15), 
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none')

plot_causal_isolation <- 
causal_effect_isolation |> 
  group_by(intervention, frugivory_type) |> 
  transmute(mu = median(probability), 
            li = quantile(probability, 0.025), 
            ls = quantile(probability, 0.975), 
            var = 'Island isolation') |> 
  unique() |> 
  ggplot(aes(intervention, mu, ymin = li, ymax = ls, 
             color = frugivory_type)) +
  geom_errorbar(width = 0, linewidth = 2, alpha = 0.5,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 1.5, 
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(~var) +
  labs(x = 'z-score difference', 
       y = expression(Delta ~'P(Fruit consumption)')) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35),
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none')


# size of island

est_size_pred <- cond_effects(posterior = post_size_pred,
                              x_bar = dat$island_size,
                              slope = 'beta_I_size',
                              type = 'averaging',
                              n = 100)

est_size_lizard <- cond_effects(posterior = post_size_lizard,
                              x_bar = dat$island_size,
                              slope = 'beta_I_size',
                              type = 'averaging',
                              n = 100)


causal_effect_size <- 
  lapply(1:2, FUN = 
           function(x) {
             
             model <- list(post_size_pred, 
                           post_size_lizard)[[x]]
             model_lab <- c('Seed predation', 
                            'Lizard frugivory')[x]
             
             d <-  # causal effect of latitude
               lapply(1:4, FUN = 
                        function(x) {
                          
                          #interventions in SD
                          i <- c(-1, 0, 1, quantile(dat$island_size, 0.1))[x] 
                          j <- c(0, 1, 2, quantile(dat$island_size, 0.9))[x]
                          
                          d <-
                            cond_effects(posterior = model,
                                         x_bar = dat$island_size,
                                         slope = 'beta_I_size',
                                         type = 'averaging',
                                         causal_effect = T,
                                         intervention1 = i,
                                         intervention2 = j)
                          
                          d[, 3]
                          
                        })
             
             d <- as_tibble(do.call('cbind', d))
             
             colnames(d)[4] <- paste('10th vs. 90th', ' percentile', 
                                     sep = '\n')
             
             d <- 
               pivot_longer(d, colnames(d), 
                            names_to = 'intervention', values_to = 'probability')
             
             d$frugivory_type <- model_lab
             d$factor <- 'Island size'
             d
           })

causal_effect_size <- do.call('rbind', causal_effect_size)

tribble(~code_color, ~group, 
        "#E41A1C", 'lizard',
        "#377EB8", 'seed predation',
        "#4DAF4A", 'total frugivory',
        "#984EA3", 'bird frugivory',
        "#FF7F00", 'seed dispersal')

beta_size <- 
  rbind(post_size_pred$beta |> 
          mutate(frugivory_type = 'Seed predation') |> 
          rename(x = beta_I_size), 
        post_size_lizard$beta |> 
          mutate(frugivory_type = 'Lizard frugivory') |> 
          rename(x = beta_I_size)) |> 
  ggplot(aes(x, fill = frugivory_type)) +
  geom_density(linewidth = 0, alpha = 0.4) +
  labs(x = expression(beta['island size']), 
       y = 'Density') +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_fill_manual(values = c(c("#E41A1C", "#377EB8"))) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35), 
        axis.title.x = element_text(size = 15), 
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none'
        )

plot_causal_size <- 
causal_effect_size |> 
  group_by(intervention, frugivory_type) |> 
  transmute(mu = median(probability), 
            li = quantile(probability, 0.025), 
            ls = quantile(probability, 0.975), 
            var = 'Island size') |> 
  unique() |> 
  ggplot(aes(intervention, mu, ymin = li, ymax = ls, 
             color = frugivory_type)) +
  geom_errorbar(width = 0, linewidth = 2, alpha = 0.5,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 1.5, 
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(~var) +
  labs(x = 'z-score difference', 
       y = expression(Delta ~'P(Fruit consumption)')) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35),
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none'
        )

beta_latitude + plot_causal_effect_latitude +
  beta_isolation + plot_causal_isolation + 
   beta_size + plot_causal_size +
  plot_layout(design = 
                '
                ab
                cd
                fg
              ')

# landscape ==

# altitude

slope_pars('post_altitude', 'beta_I_alt')


est_alt_pred <- cond_effects(posterior = post_altitude_pred,
                             x_bar = dat$altitude_m,
                             slope = 'beta_I_alt',
                             type = 'averaging',
                             n = 100)


est_alt_disp <- cond_effects(posterior = post_altitude_disp,
                             x_bar = dat$altitude_m,
                             slope = 'beta_I_alt',
                             type = 'averaging',
                             n = 100)

est_alt_bird <- cond_effects(posterior = post_altitude_bird,
                             x_bar = dat$altitude_m,
                             slope = 'beta_I_alt',
                             type = 'averaging',
                             n = 100)

est_alt_lizard <- cond_effects(posterior = post_altitude_lizard,
                             x_bar = dat$altitude_m,
                             slope = 'beta_I_alt',
                             type = 'averaging',
                             n = 100)


causal_effect_alt <- 
  lapply(1:4, FUN = 
           function(k) {
             
             model <- list(post_altitude_pred, 
                           post_altitude_disp,
                           post_altitude_bird,
                           post_altitude_lizard)[[k]]
             
             model_lab <- c('Seed predation', 
                            'Seed dispersal',
                            'Bird frugivory',
                            'Lizard frugivory')[k]
             
             d <-  # causal effect of latitude
               lapply(1:4, FUN = 
                        function(x) {
                          
                          #interventions in SD
                          i <- c(-1, 0, 1, quantile(dat$altitude_m, 0.1))[x] 
                          j <- c(0, 1, 2, quantile(dat$altitude_m, 0.9))[x]
                          
                          d <-
                            cond_effects(posterior = model,
                                         x_bar = dat$altitude_m,
                                         slope = 'beta_I_alt',
                                         type = 'averaging',
                                         causal_effect = T,
                                         intervention1 = i,
                                         intervention2 = j)
                          
                          d[, 3]
                          
                        })
             
             d <- as_tibble(do.call('cbind', d))
             
             colnames(d)[4] <- paste('10th vs. 90th', 'percentile', 
                                     sep = '\n')
             
             d <- 
               pivot_longer(d, colnames(d), 
                            names_to = 'intervention', values_to = 'probability')
             
             d$frugivory_type <- model_lab
             d$factor <- 'Altitude'
             d
           })

causal_effect_alt <- do.call('rbind', causal_effect_alt)

tribble(~code_color, ~group, 
        "#E41A1C", 'lizard',
        "#377EB8", 'seed predation',
        "#4DAF4A", 'total frugivory',
        "#984EA3", 'bird frugivory',
        "#FF7F00", 'seed dispersal')

beta_altitude <- 
  rbind(post_altitude_pred$beta |> 
          mutate(frugivory_type = 'Seed predation') |> 
          rename(x = beta_I_alt),
        post_altitude_disp$beta |> 
          mutate(frugivory_type = 'Seed dispersal') |> 
          rename(x = beta_I_alt),
        post_altitude_bird$beta |> 
          mutate(frugivory_type = 'Bird frugivory') |> 
          rename(x = beta_I_alt), 
        post_altitude_lizard$beta |> 
          mutate(frugivory_type = 'Lizard frugivory') |> 
          rename(x = beta_I_alt)) |> 
  ggplot(aes(x, fill = frugivory_type)) +
  geom_density(linewidth = 0, alpha = 0.4) +
  labs(x = expression(beta['island altitude']), 
       y = 'Density') +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_fill_manual(values = c("#984EA3", "#E41A1C", '#FF7F00', '#377EB8')) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35), 
        axis.title = element_text(size = 11), 
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none'
  )

plot_causal_alt <- 
causal_effect_alt |> 
  group_by(intervention, frugivory_type) |> 
  transmute(mu = median(probability), 
            li = quantile(probability, 0.025), 
            ls = quantile(probability, 0.975), 
            type = ifelse(grepl('^Min', intervention), 
                          'b', 'a'), 
            var = 'Island altitude') |> 
  unique() |> 
  ggplot(aes(intervention, mu, ymin = li, ymax = ls, 
             color = frugivory_type)) +
  geom_errorbar(width = 0, linewidth = 2, alpha = 0.5,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 1.5, 
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#984EA3", "#E41A1C", '#FF7F00', '#377EB8')) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(~var) +
  # facet_wrap(~type, nrow = 1, scales = 'free') +
  # force_panelsizes(
  #   cols = unit(c(3, 2), "in"),
  #   respect = T
  # ) +
  labs(x = 'z-score difference', 
       y = expression(Delta ~'P(Fruit consumption)')) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        #strip.text = element_blank(),
        axis.line = element_line(linewidth = 0.35),
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none',
        axis.title = element_text(size = 10)
        # axis.title.y = element_text(margin = margin(r = 1)),
        # # Adjust plot margins if needed
        # plot.margin = margin(10, 10, 10, 40)
  )

beta_latitude + plot_causal_effect_latitude +
  beta_isolation + plot_causal_isolation + 
  beta_size + plot_causal_size +
  beta_altitude + plot_causal_alt +
  plot_layout(design = 
                '
                ab
                cd
                fg
                hi
              ')




# human footprint

est_foot_bird <- cond_effects(posterior = post_footprint_bird,
                               x_bar = dat$human_footprint,
                               slope = 'beta_H_foot',
                               type = 'averaging',
                               n = 100)

causal_effect_foot <- 
  lapply(1, FUN = 
           function(k) {
             
             model <- list(post_footprint_bird)[[k]]
             
             model_lab <- c('Bird frugivory')[k]
             
             d <-  # causal effect of latitude
               lapply(1:4, FUN = 
                        function(x) {
                          
                          #interventions in SD
                          i <- c(-1, 0, 1, quantile(dat$human_footprint, 0.1))[x] 
                          j <- c(0, 1, 2, quantile(dat$human_footprint, 0.9))[x]
                          
                          d <-
                            cond_effects(posterior = model,
                                         x_bar = dat$human_footprint,
                                         slope = 'beta_H_foot',
                                         type = 'averaging',
                                         causal_effect = T,
                                         intervention1 = i,
                                         intervention2 = j)
                          
                          d[, 3]
                          
                        })
             
             d <- as_tibble(do.call('cbind', d))
             
             colnames(d)[4] <- paste('10th vs. 90th', 'percentile', 
                                     sep = '\n')
             
             d <- 
               pivot_longer(d, colnames(d), 
                            names_to = 'intervention', values_to = 'probability')
             
             d$frugivory_type <- model_lab
             d$factor <- 'Human footprint'
             d
           })

causal_effect_foot <- do.call('rbind', causal_effect_foot)


tribble(~code_color, ~group, 
        "#E41A1C", 'lizard',
        "#377EB8", 'seed predation',
        "#4DAF4A", 'total frugivory',
        "#984EA3", 'bird frugivory',
        "#FF7F00", 'seed dispersal')

beta_foot <- 
  rbind(post_footprint_bird$beta |> 
          mutate(frugivory_type = 'Bird frugivory') |> 
          rename(x = beta_H_foot)) |> 
  ggplot(aes(x, fill = frugivory_type)) +
  geom_density(linewidth = 0, alpha = 0.4) +
  labs(x = expression(beta['Human footprint']), 
       y = 'Density') +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_fill_manual(values = c("#984EA3")) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35), 
        axis.title = element_text(size = 11), 
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none'
  )

plot_causal_foot <- 
causal_effect_foot |> 
  group_by(intervention, frugivory_type) |> 
  transmute(mu = median(probability), 
            li = quantile(probability, 0.025), 
            ls = quantile(probability, 0.975), 
            type = ifelse(grepl('^Min', intervention), 
                          'b', 'a'),
            var = 'Human footprint') |> 
  unique() |> 
  ggplot(aes(intervention, mu, ymin = li, ymax = ls, 
             color = frugivory_type)) +
  geom_errorbar(width = 0, linewidth = 2, alpha = 0.5,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 1.5, 
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#984EA3")) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(~var) +
  # facet_wrap(~type, nrow = 1, scales = 'free') +
  # force_panelsizes(
  #   cols = unit(c(3, 2), "in"),
  #   respect = T
  # ) +
  labs(x = 'z-score difference', 
       y = expression(Delta ~'P(Fruit consumption)')) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        #strip.text = element_blank(),
        axis.line = element_line(linewidth = 0.35),
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none',
        axis.title = element_text(size = 10)
        # axis.title.y = element_text(margin = margin(r = 1)),
        # # Adjust plot margins if needed
        # plot.margin = margin(10, 10, 10, 40)
  )

beta_latitude + plot_causal_effect_latitude +
  beta_isolation + plot_causal_isolation + 
  beta_size + plot_causal_size +
  beta_altitude + plot_causal_alt +
  beta_foot + plot_causal_foot
  plot_layout(design = 
                '
                abcd
                fghi
                jk##
              ')


# native vegetation
# ****no effect

# bush cover

slope_pars('post_bush', 'beta_bush')

bush_merge <- mod_bush_pred$draws('bush_merge', format = 'matrix')

bush_merge <- apply(bush_merge, 2, median)

est_bush_pred <- cond_effects(posterior = post_bush_pred,
                             x_bar = bush_merge,
                             slope = 'beta_bush',
                             type = 'averaging',
                             n = 100)

est_bush_bird <- cond_effects(posterior = post_bush_bird,
                              x_bar = bush_merge,
                              slope = 'beta_bush',
                              type = 'averaging',
                              n = 100)


causal_effect_bush <- 
  lapply(1:2, FUN = 
           function(k) {
             
             model <- list(post_bush_pred, 
                           post_bush_bird)[[k]]
             
             model_lab <- c('Seed predation', 
                            'Bird frugivory')[k]
             
             d <-  # causal effect of latitude
               lapply(1:4, FUN = 
                        function(x) {
                          
                          #interventions in SD
                          i <- c(-1, 0, 1, quantile(bush_merge, 0.1))[x] 
                          j <- c(0, 1, 2, quantile(bush_merge, 0.9))[x]
                          
                          d <-
                            cond_effects(posterior = model,
                                         x_bar = bush_merge,
                                         slope = 'beta_bush',
                                         type = 'averaging',
                                         causal_effect = T,
                                         intervention1 = i,
                                         intervention2 = j)
                          
                          d[, 3]
                          
                        })
             
             d <- as_tibble(do.call('cbind', d))
             
             colnames(d)[4] <- paste('10th vs. 90th', 'percentile', 
                                     sep = '\n')
             
             d <- 
               pivot_longer(d, colnames(d), 
                            names_to = 'intervention', values_to = 'probability')
             
             d$frugivory_type <- model_lab
             d$factor <- 'Bush cover'
             d
           })

causal_effect_bush <- do.call('rbind', causal_effect_bush)


tribble(~code_color, ~group, 
        "#E41A1C", 'lizard',
        "#377EB8", 'seed predation',
        "#4DAF4A", 'total frugivory',
        "#984EA3", 'bird frugivory',
        "#FF7F00", 'seed dispersal')

beta_bush <- 
  rbind(post_bush_bird$beta |> 
          mutate(frugivory_type = 'Bird frugivory') |> 
          rename(x = beta_bush), 
        post_bush_pred$beta |> 
          mutate(frugivory_type = 'Seed predation') |> 
          rename(x = beta_bush)) |> 
  ggplot(aes(x, fill = frugivory_type)) +
  geom_density(linewidth = 0, alpha = 0.4) +
  labs(x = expression(beta['Bush cover']), 
       y = 'Density') +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_fill_manual(values = c("#984EA3", "#377EB8")) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        axis.line = element_line(linewidth = 0.35), 
        axis.title = element_text(size = 11), 
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none'
  )

plot_causal_bush <- 
causal_effect_bush |> 
  group_by(intervention, frugivory_type) |> 
  transmute(mu = median(probability), 
            li = quantile(probability, 0.025), 
            ls = quantile(probability, 0.975), 
            type = ifelse(grepl('^Min', intervention), 
                          'b', 'a'), 
            var = 'Bush cover') |> 
  unique() |> 
  ggplot(aes(intervention, mu, ymin = li, ymax = ls, 
             color = frugivory_type)) +
  geom_errorbar(width = 0, linewidth = 2, alpha = 0.5,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 1.5, 
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#984EA3", "#377EB8")) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(~var) +
  # facet_wrap(~type, nrow = 1, scales = 'free') +
  # force_panelsizes(
  #   cols = unit(c(3, 2), "in"),
  #   respect = T
  # ) +
  labs(x = 'z-score difference', 
       y = expression(Delta ~'P(Fruit consumption)')) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        #strip.text = element_blank(),
        axis.line = element_line(linewidth = 0.35),
        text = element_text(family = 'Times New Roman'), 
        legend.position = 'none',
        axis.title = element_text(size = 10)
        # axis.title.y = element_text(margin = margin(r = 1)),
        # # Adjust plot margins if needed
        # plot.margin = margin(10, 10, 10, 40)
  )

beta_latitude + plot_causal_effect_latitude +
  beta_isolation + plot_causal_isolation + 
  beta_size + plot_causal_size +
  beta_altitude + plot_causal_alt +
  beta_foot + plot_causal_foot +
  beta_bush + plot_causal_bush
plot_layout(design = 
              '
                abcd
                fghi
                jklm
              ')



plot_scatter_isolation <- 
  rbind(est_isolation_tot |> 
        mutate(var = 'Island isolation', 
               type = 'Total frugivory'), 
      est_isolation_pred |> 
        mutate(var = 'Island isolation', 
               type = 'Seed predation'),
      est_isolation_lizard |> 
        mutate(var = 'Island isolation', 
               type = 'Lizard frugivory')) |> 
  mutate(x_reverted = mean(d$isolation) + x * sd(d$isolation)) |> 
  ggplot(aes(x_reverted, y, ymin = li, ymax = ls)) +
  geom_ribbon(aes(fill = type), alpha = 0.5) +
  geom_line(aes(color = type)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, 'Set1')[1:3]) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, 'Set1')[1:3]) +
  labs(y = 'P(fruit consumption)', 
       x = 'Island isolation (km)') +
  #facet_wrap(~var, scales = 'free', nrow = 1) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )

plot_scatter_altitude <- 
  rbind(est_alt_pred |> 
        mutate(var = 'Island altitude', 
               type = 'Seed predation'), 
      est_alt_bird |> 
        mutate(var = 'Island altitude', 
               type = 'Bird frugivory'), 
      est_alt_lizard |> 
        mutate(var = 'Island altitude', 
               type = 'Lizard frugivory'), 
      est_alt_disp |> 
        mutate(var = 'Island altitude', 
               type = 'Seed dispersion')) |> 
  mutate(x_reverted = mean(d$altitude_m) + x * sd(d$altitude_m)) |> 
  ggplot(aes(x_reverted, y, ymin = li, ymax = ls)) +
  geom_ribbon(aes(fill = type), alpha = 0.5) +
  geom_line(aes(color = type)) +
  scale_color_manual(values = c("#984EA3", "#E41A1C", "#FF7F00", "#377EB8")) +
  scale_fill_manual(values = c("#984EA3", "#E41A1C", "#FF7F00", "#377EB8")) +
  labs(y = 'P(fruit consumption)', 
       x = 'Altitude (m)') +
  #facet_wrap(~var, scales = 'free', nrow = 1) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )

plot_scatter_bush <- 
  rbind(est_bush_pred |> 
        mutate(var = 'Bush cover', 
               type = 'Seed predation'), 
      est_bush_bird |> 
        mutate(var = 'Bush cover', 
               type = 'Bird frugivory')) |> 
  mutate(x_revert = mean(d$bush_cover, na.rm = T) + x * 
           sd(d$bush_cover, na.rm = T))|> 
  ggplot(aes(x_revert, y, ymin = li, ymax = ls)) +
  geom_ribbon(aes(fill = type), alpha = 0.5) +
  geom_line(aes(color = type)) +
  scale_color_manual(values = c("#984EA3", "#377EB8")) +
  scale_fill_manual(values = c("#984EA3", "#377EB8")) +
  labs(y = 'P(fruit consumption)', 
       x = 'Bush cover (%)') +
  #facet_wrap(~var, scales = 'free', nrow = 1) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )
  

plot_scatter_latitude <- 
  rbind(est_latitude_tot |> 
        mutate(var = 'Latitude', 
               type = 'Total frugivory')) |> 
  mutate(x_revert = mean(d$lat) + x * sd(d$lat))|> 
  ggplot(aes(x_revert, y, ymin = li, ymax = ls)) +
  geom_ribbon(alpha = 0.5, fill = "#4DAF4A") +
  geom_line(color = "#4DAF4A") +
  #scale_color_manual(values = c('#D92525', '#7A577A')) +
  #scale_fill_manual(values = c('#D92525', '#7A577A')) +
  labs(y = 'P(fruit consumption)', 
       x = 'Latitude') +
  #facet_wrap(~var, scales = 'free', nrow = 1) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )


plot_scatter_size <- 
  rbind(est_size_lizard |> 
        mutate(var = 'Island size', 
               type = 'Lizard frugivory'), 
      est_size_pred |> 
        mutate(var = 'Island size', 
               type = 'Seed predation')) |> 
  mutate(x_revert = (mean(d$island_size) + x * sd(d$island_size))) |> 
  ggplot(aes(x_revert, y, ymin = li, ymax = ls)) +
  geom_ribbon(aes(fill = type), alpha = 0.5) +
  geom_line(aes(color = type)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  labs(y = 'P(fruit consumption)', 
       x = 'Island size') +
  labs(y = 'P(fruit consumption)', 
       x = expression('Island size (km)'^2)) +
  #facet_wrap(~var, scales = 'free', nrow = 1) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )

plot_scatter_humanF <- 
  rbind(est_foot_bird |> 
        mutate(var = 'Human footprint', 
               type = 'Bird frugivory')) |> 
    mutate(x_revert = mean(d$human_footprint) + x * sd(d$human_footprint)) |> 
  ggplot(aes(x_revert, y, ymin = li, ymax = ls)) +
  geom_ribbon(fill = '#984EA3', alpha = 0.5) +
  geom_line(color = '#984EA3') +
  labs(y = 'P(fruit consumption)', 
       x = 'Human footprint') +
  #facet_wrap(~var, scales = 'free', nrow = 1) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.position = 'none',
    legend.box.background = element_blank(), 
    legend.key.size = unit(2.5, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )


all_betas <- 
  rbind(slope_pars('post_latitude', 'beta_lat') |> 
          mutate(predictor = 'Latitude'), 
        slope_pars('post_isolation', 'beta_I_isolation') |> 
          mutate(predictor = 'Island  \nisolation'), 
        slope_pars('post_size', 'beta_I_size') |> 
          mutate(predictor = 'Island\narea  '), 
        slope_pars('post_altitude', 'beta_I_alt') |> 
          mutate(predictor = 'Island \naltitude'), 
        slope_pars('post_footprint', 'beta_H_foot') |> 
          mutate(predictor = 'Human \nfootprint'), 
        slope_pars('post_nativeV', 'beta_NV') |> 
          mutate(predictor = 'Native   \nvegetation'), 
        slope_pars('post_bush', 'beta_bush') |> 
          mutate(predictor = 'Bush\ncover'))

all_betas$significant <- 
  all_betas$`P(beta > 0)` >= 0.8 | all_betas$`P(beta < 0)` >= 0.8

all_betas$`Fruit consumption`[grep('^Total', all_betas$`Fruit consumption`)] <- 
  'Frugivory'

all_betas[all_betas$significant == T, ]

unique(all_betas[all_betas$significant == T, ]$predictor)

split(all_betas[all_betas$significant == T, ], 
      all_betas[all_betas$significant == T, ]$`Fruit consumption`)
#"#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" 
# lizard  predation. total      bird.    'disp
plot_beta_continuous <- 
  all_betas |> 
  ggplot(aes(predictor, mu, ymin = li, ymax = ls, 
             color = `Fruit consumption`)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.6), 
                linewidth = 1.5, 
                alpha = ifelse(all_betas$significant, 1, 0.2)) +
  geom_point(position = position_dodge(width = 0.6), size = 2, 
             alpha = ifelse(all_betas$significant, 1, 0.3)) +
  scale_color_manual(values = c("#984EA3", "#E41A1C", "#FF7F00", 
                                "#377EB8", "#4DAF4A")) +
  # annotate("segment",
  #          x = c(1.13, 
  #                4.13, 
  #                4, 
  #                5-0.13, 
  #                5+0.13, 
  #                3+0.13, 
  #                6-0.13), 
  #          xend = c(1.13,
  #                   4.13, 
  #                   4, 
  #                   5-0.13, 
  #                   5+0.13, 
  #                   3+0.13, 
  #                   6-0.13), # y
  #          y = c(-0.8, 
  #                -0.35, 
  #                -1.1, 
  #                -0.65, 
  #                -0.83, 
  #                -0.87, 
  #                -0.73), 
  #          yend = c(-0.5, 
  #                   -0.05, 
  #                   -0.8, 
  #                   -0.35, 
  #                   -0.53, 
  #                   -0.57, 
  #                   -0.43), # x
  #          arrow = arrow(type = "closed", 
  #                        length = unit(1.5, "mm"),
  #                        angle = 30),
  #          color = c('#7A577A', 
  #                    '#7A577A', 
  #                    '#F2C230', 
  #                    '#D92525', 
  #                    '#7A577A', 
  #                    '#7A577A',
  #                    '#D92525'),
  #          linewidth = 0.5) +
  labs(y = expression(beta['slope']), x = 'Predictors') +
  theme_classic() +
  coord_flip() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    legend.position = c(0.2, 0.92),
    legend.box.background = element_blank(), 
    legend.key.size = unit(3, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 12)
  )

plot_beta_continuous
  
ggsave('figure_3.jpg', width = 10, height = 12.5, units = 'cm', dpi = 500)

layout2 <- 
  '
  aaa
  aaa
  bcd
  fgh
'

plot_beta_continuous +
  plot_scatter_isolation + 
  plot_scatter_altitude + labs(y = NULL) + 
  plot_scatter_bush + labs(y = NULL) +
  plot_scatter_latitude + 
  plot_scatter_size + labs(y = NULL) + 
  plot_scatter_humanF + labs(y = NULL) +
  plot_layout(design = layout2) +
  plot_annotation(tag_levels = 'a', 
                  tag_prefix = '(', 
                  tag_suffix = ')')
  

ggsave('effects_plot.jpg', width = 15, height = 22, units = 'cm', dpi = 1e3)

p_layout <- 
  '
  bbcc
  ddee
  ffgg
'

plot_causal_effect_latitude + labs(x = "") + 
  plot_causal_isolation + labs(y = "", x = "") + 
  plot_causal_size + labs(x = "") + 
  plot_causal_alt + labs(y = "", x = '') + 
  plot_causal_foot + plot_causal_bush + labs(y = "") + 
  plot_layout(design = p_layout) +
  plot_annotation(tag_levels = 'a', 
                  tag_prefix = '(', 
                  tag_suffix = ')') &
  theme(plot.tag.position = c(0.1, 1))



tribble(~code_color, ~group, 
        "#E41A1C", 'lizard',
        "#377EB8", 'seed predation',
        "#4DAF4A", 'total frugivory',
        "#984EA3", 'bird frugivory',
        "#FF7F00", 'seed dispersal')

c("#984EA3", "#E41A1C", "#FF7F00", 
  "#377EB8", "#4DAF4A")

causal_df <- 
  rbind(causal_effect_latitude |> 
        group_by(intervention) |> 
        transmute(frugivory_type = 'Total frugivory',
                  mu = median(probability), 
                  li = quantile(probability, 0.025), 
                  ls = quantile(probability, 0.975), 
                  var = 'Latitude') |> 
        unique(),
      causal_effect_isolation |> 
        group_by(intervention, frugivory_type) |> 
        transmute(mu = median(probability), 
                  li = quantile(probability, 0.025), 
                  ls = quantile(probability, 0.975), 
                  var = 'Island isolation') |> 
        unique(), 
      causal_effect_size |> 
        group_by(intervention, frugivory_type) |> 
        transmute(mu = median(probability), 
                  li = quantile(probability, 0.025), 
                  ls = quantile(probability, 0.975), 
                  var = 'Island size') |> 
        unique(), 
      causal_effect_alt |> 
        group_by(intervention, frugivory_type) |> 
        transmute(mu = median(probability), 
                  li = quantile(probability, 0.025), 
                  ls = quantile(probability, 0.975),
                  var = 'Island altitude') |> 
        unique(), 
      causal_effect_foot |> 
        group_by(intervention, frugivory_type) |> 
        transmute(mu = median(probability), 
                  li = quantile(probability, 0.025), 
                  ls = quantile(probability, 0.975),
                  var = 'Human footprint') |> 
        unique(), 
      causal_effect_bush |> 
        group_by(intervention, frugivory_type) |> 
        transmute(mu = median(probability), 
                  li = quantile(probability, 0.025), 
                  ls = quantile(probability, 0.975), 
                  var = 'Bush cover') |> 
        unique())

causal_df$frugivory_type[grep('^Total', causal_df$frugivory_type)] <- 
  'Frugivory'

causal_df$intervention[grep("10th vs. 90th\npercentile", 
                            causal_df$intervention)] <- 
  '10th vs. 90th\n percentile'

unique(causal_df$intervention)

ggplot(causal_df, aes(intervention, mu, ymin = li, ymax = ls, 
           color = frugivory_type)) +
  geom_errorbar(width = 0, linewidth = 2, alpha = 0.5,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 1.5, 
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#984EA3", "#E41A1C", "#FF7F00", 
                                "#377EB8", "#4DAF4A")) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(~var, scales = 'free_y') +
  # facet_wrap(~type, nrow = 1, scales = 'free') +
  # force_panelsizes(
  #   cols = unit(c(3, 2), "in"),
  #   respect = T
  # ) +
  labs(x = '\nz-score difference of predictor variable', 
       y = expression(Delta ~'P(Fruit consumption)')) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 9),
        axis.line = element_line(linewidth = 0.35),
        text = element_text(family = 'Times New Roman', size = 10), 
        legend.position = 'top',
        legend.title = element_blank(),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 9.5)
        # axis.title.y = element_text(margin = margin(r = 1)),
        # # Adjust plot margins if needed
        # plot.margin = margin(10, 10, 10, 40)
  )


ggsave('figure_4.jpg', width = 20, height = 12, units = 'cm', dpi = 500)


all_betas[all_betas$significant == T, ]

unique(all_betas[all_betas$significant == T, ]$predictor)

split(all_betas[all_betas$significant == T, ], 
      all_betas[all_betas$significant == T, ]$`Fruit consumption`)


causal_effects_estimates <- 
  lapply(ls()[grep('^causal', ls())][-c(1, 8)], FUN = 
         function(x) {
           
           d <- get(x)
           
           d |>
             group_by(intervention, frugivory_type) |>
             transmute(factor = factor,
                       mean_effect = median(probability),
                       li = quantile(probability, 0.025),
                       ls = quantile(probability, 0.975),
                       sd_effect = sd(probability),
                       `p(causal effect > 0)` = mean(probability > 0),
                       `p(causal effect < 0)` = mean(probability < 0)) |>
             unique()
         })

causal_effects_estimates <- do.call('rbind', causal_effects_estimates)


zzz <- 
  lapply(ls()[grep('^causal', ls())][-c(1, 8)], FUN = 
         function(x) {
           
           d <- get(x)
           
           d <- d[-grep('^10', d$intervention), ]
           
           d |>
             group_by(frugivory_type, factor) |>
             transmute(mu = median(probability), 
                       sd = sd(probability),
                       li = quantile(probability, 0.025),
                       ls = quantile(probability, 0.975),
                       `p > 0` = mean(probability > 0),
                       `p < 0` = mean(probability < 0)) |>
             unique()
           
         })

zzz <- do.call('rbind', zzz)

zzz[order(zzz$frugivory_type), ]

sessionInfo()
