# ======== Loading packages =======

sapply(c('cmdstanr', 'readxl', 'magrittr', 'dplyr', 'ggplot2', 
         'tidyr', 'tibble', 'forcats', 'rethinking', 
         'cowplot', 'bayesplot'), 
       library, character.only = T)

source('functions_mod_diagnostics.r')

# ======== Data preparation =======

data <- readRDS('data.rds')

codes <- 
  lapply(data$codes, function(x) x[order(x$code), ])

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
                         n = 100) {
  
  x_seq <- seq(min(x_bar), max(x_bar), length.out = n)
  
  y <- lapply(seq_along(x_seq), FUN = 
                function(i) {
                  
                  x <- x_seq[i]
                  
                  est <- 
                    with(posterior, 
                         {
                           inv_logit(alpha[[1]] +
                                       beta[[slope]] * x +
                                       apply(TI, 1, mean) +
                                       apply(p_island, 1, mean) +
                                       apply(p_country, 1, mean) +
                                       apply(p_grid, 1, mean) +
                                       apply(p_plant, 1, mean) +
                                       apply(p_realm, 1, mean) +
                                       apply(p_ecoR, 1, mean) +
                                       apply(p_biome, 1, mean))
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
                           y = median(est),
                           indx = i)
                  }
                  
                })
  
  do.call('rbind', y)
  
}



slope_pars <- 
  function(posterior, 
           slope) {
    
    env <- ls(globalenv())

    post <- env[grep(posterior, env)]
    
    post <-
      lapply(post, FUN =
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
      c('Seed dispersion', 'Seed predation', 'Frugivory')
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
for (i in 1:200) lines(density(ppcheck_latitude_tot[i, ], lwd = 0.1))
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

realm_pred <- average_effects(n_levels = ncol(post_latitude_pred$p_realm),
                           posterior = post_latitude_pred, 
                           x_var1 = 'realm', 
                           x_var2 = 'island_type', 
                           par1 = 'p_realm', 
                           par2 = 'TI')


rbind(full_join(realm_tot, codes$real, 'code') |> 
        mutate(type = 'Frugivory'), 
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
  facet_wrap(~type, scales = 'fixed') +
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

mod_isolation_tot$save_object('mod_isolation_tot.rds')

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

mod_isolation_disp$save_object('mod_isolation_disp.rds')

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

mod_isolation_pred$save_object('mod_isolation_pred.rds')

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

mod_size_tot$save_object('mod_size_tot.rds')

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

mod_size_disp$save_object('mod_size_disp.rds')

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

mod_size_pred$save_object('mod_size_pred.rds')

sum_size_pred <- mod_size_pred$summary()
mod_diagnostics(mod_size_pred, sum_size_pred)

ppcheck_size_pred <- mod_size_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
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

mod_altitude_tot$save_object('mod_altitude_tot.rds')
mod_altitude_disp$save_object('mod_altitude_dip.rds')
mod_altitude_pred$save_object('mod_altitude_pred.rds')

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

sum_altitude_disp <- mod_altitude_disp$summary()
mod_diagnostics(mod_altitude_disp, sum_altitude_disp)
ppcheck_altitude_disp <- mod_altitude_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
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

sum_altitude_pred <- mod_altitude_pred$summary()
mod_diagnostics(mod_altitude_pred, sum_altitude_pred)

ppcheck_altitude_pred <- mod_altitude_pred$draws('ppcheck', format = 'matrix')

plot(density(dat$predation), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
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

mod_footprint_tot$save_object('mod_footprint_tot.rds')
mod_footprint_disp$save_object('mod_footprint_dip.rds')
mod_footprint_pred$save_object('mod_footprint_pred.rds')

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

sum_footprint_disp <- mod_footprint_disp$summary()
mod_diagnostics(mod_footprint_disp, sum_footprint_disp)
ppcheck_footprint_disp <- mod_footprint_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
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

mod_nativeV_tot$save_object('mod_nativeV_tot.rds')
mod_nativeV_disp$save_object('mod_nativeV_dip.rds')
mod_nativeV_pred$save_object('mod_nativeV_pred.rds')

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

sum_nativeV_disp <- mod_nativeV_disp$summary()
mod_diagnostics(mod_nativeV_disp, sum_nativeV_disp)
ppcheck_nativeV_disp <- mod_nativeV_disp$draws('ppcheck', format = 'matrix')

plot(density(dat$dispersion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.4))
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

sum_nativeV_pred <- mod_nativeV_pred$summary()
mod_diagnostics(mod_nativeV_pred, sum_nativeV_pred)

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

mod_bush_tot$save_object('mod_bush_tot.rds')
mod_bush_disp$save_object('mod_bush_dip.rds')
mod_bush_pred$save_object('mod_bush_pred.rds')

# ======= Plots =====
# 
# ============= Slops ===========

# error bars

rbind(pivot_longer(post_bush_tot$beta, 'beta_bush') |> 
        mutate(type = 'Frugivory', 
               effect = 'Bush cover'), 
      pivot_longer(post_bush_disp$beta, 'beta_bush') |> 
        mutate(type = 'Seed dispersion', 
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
        mutate(type = 'Dispersion'), 
      full_join(TI_pred, codes$island_type, 'code') |> 
        mutate(type = 'Predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> 
  ggplot(aes(island, mu, ymin = li, ymax = ls, color = type)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.4), 
                linewidth = 1.5, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  scale_color_manual(values = c('#F2C230', '#7A577A')) +
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
    legend.key.size = unit(2.5, 'mm'), 
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
                    label_x = 0.16, 
                    label_y = 0.99), 
          plot_grid(plot_frugivory_realms, 
                    plot_disp_pred_realms, 
                    nrow = 1, 
                    labels = c('(c)', '(d)'), 
                    label_size = 9, 
                    label_fontfamily = 'Times New Roman', 
                    label_fontface = 'plain', 
                    label_x = c(0.13, 0.145), 
                    label_y = 0.99),
          ncol = 1)

ggsave('plot_type_island_realm.jpg', dpi = 1e3, width = 15, height = 10, units = 'cm')


# estimations

# probability among types of islands
rbind(full_join(TI_tot, codes$island_type, 'code') |> 
        mutate(type = 'Frugivory'), 
      full_join(TI_disp, codes$island_type, 'code') |> 
        mutate(type = 'Dispersion'), 
      full_join(TI_pred, codes$island_type, 'code') |> 
        mutate(type = 'Predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() 

lapply(1, FUN = 
         function(x) {
           
           d <- list(TI_tot, TI_disp, TI_pred)
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
          
           names(d) <- c('Frugivory', 'Dispersion', 'Predation')
           
           d <- 
             lapply(d, FUN = 
                      function(jj) {
                        
                        prob <-
                          lapply(1:3, FUN =
                                   function(i) {
                                     df <-
                                       lapply(1:3, FUN =
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
           
           d[d$probability != 0, ]
         })


# frugivory among realms

# average probabilities 
rbind(full_join(realm_tot, codes$real, 'code') |> 
        mutate(type = 'Frugivory'), 
      full_join(realm_disp, codes$real, 'code') |> 
        mutate(type = 'Seed dispersal'), 
      full_join(realm_pred, codes$real, 'code') |> 
        mutate(type = 'Seed predation')) |> 
  group_by(type, island) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975)) |> 
  unique() |> print(n = 21)


# contrast among realms 

lapply(1, FUN = 
         function(x) {
           
           d <- list(realm_tot, realm_disp, realm_pred)
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
           
           names(d) <- c('Frugivory', 'Dispersion', 'Predation')
           
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

lapply(1, FUN = 
         function(x) {
           
           d <- list(realm_tot, realm_disp, realm_pred)
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
           
           names(d) <- c('Frugivory', 'Dispersion', 'Predation')
           
           d <-
             lapply(d, FUN =
                      function(jj) {

                        prob <-
                          lapply(codes$real$code, FUN =
                                   function(i) {
                                     df <-
                                       lapply(codes$real$code, FUN =
                                                function(j) {
                                                  p <- mean(jj[[i]] > jj[[j]])
                                                  dif <- jj[[i]] - jj[[j]]
                                                  contrast <-
                                                    paste(codes$real$island[i],
                                                          codes$real$island[j],
                                                          sep = ' > ')

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
           
         })[[1]] |> print(n = 200) 

plot_frugivory_realms <- 
  full_join(realm_tot, codes$real, 'code') |> 
        mutate(type = 'Frugivory') |> 
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
  scale_color_manual(values = c('#F2C230', '#7A577A')) +
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



# Continuous variabes ====================

# latitude (mu and CI)

rbind(slope_pars('post_latitude', 'beta_lat') |> 
        mutate(predictor = 'Latitude'), 
      slope_pars('post_isolation', 'beta_I_isolation') |> 
        mutate(predictor = 'Island isolation'), 
      slope_pars('post_size', 'beta_I_size') |> 
        mutate(predictor = 'Island size'), 
      slope_pars('post_altitude', 'beta_I_alt') |> 
        mutate(predictor = 'Island altitude'), 
      slope_pars('post_footprint', 'beta_H_foot') |> 
        mutate(predictor = 'Human footprint'), 
      slope_pars('post_nativeV', 'beta_NV') |> 
        mutate(predictor = 'Native vegetation'), 
      slope_pars('post_bush', 'beta_bush') |> 
        mutate(predictor = 'Bush cover')) |> 
  print(n = 21)


est_latitude_tot <- cond_effects(posterior = post_latitude_tot,
                                 x_bar = dat$lat,
                                 slope = 'beta_lat',
                                 type = 'averaging',
                                 n = 100)

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

# size of island

est_size_pred <- cond_effects(posterior = post_size_pred,
                              x_bar = dat$island_size,
                              slope = 'beta_I_size',
                              type = 'averaging',
                              n = 100)

# landscape ============

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

# human footprint
# *****no effects 




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

plot_scatter_plot1 <- 
  rbind(est_isolation_tot |> 
        mutate(var = 'Island isolation', 
               type = 'Frugivory'), 
      est_isolation_pred |> 
        mutate(var = 'Island isolation', 
               type = 'Predation'),
      est_alt_pred |> 
        mutate(var = 'Island altitude', 
               type = 'Predation'), 
      est_alt_disp |> 
        mutate(var = 'Island altitude', 
               type = 'Dispersion'), 
      est_bush_pred |> 
        mutate(var = 'Bush cover', 
               type = 'Predation')) |> 
  ggplot(aes(x, y, ymin = li, ymax = ls)) +
  geom_ribbon(aes(fill = type), alpha = 0.5) +
  geom_line(aes(color = type)) +
  scale_color_manual(values = c('#F2C230', '#D92525', '#7A577A')) +
  scale_fill_manual(values = c('#F2C230', '#D92525', '#7A577A')) +
  labs(y = 'P(fruit consumption)', 
       x = 'z-scores') +
  facet_wrap(~var, scales = 'free', nrow = 1) +
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

plot_scatter_plot2 <- 
  rbind(est_latitude_tot |> 
        mutate(var = 'Latitude', 
               type = 'Frugivory'), 
      est_size_pred |> 
        mutate(var = 'Island size', 
               type = 'Predation')) |> 
  ggplot(aes(x, y, ymin = li, ymax = ls)) +
  geom_ribbon(aes(fill = type), alpha = 0.5) +
  geom_line(aes(color = type)) +
  scale_color_manual(values = c('#D92525', '#7A577A')) +
  scale_fill_manual(values = c('#D92525', '#7A577A')) +
  labs(y = 'P(fruit consumption)', 
       x = 'z-scores') +
  facet_wrap(~var, scales = 'free', nrow = 1) +
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

plot_beta_continuous <- 
  rbind(slope_pars('post_latitude', 'beta_lat') |> 
          mutate(predictor = 'Latitude'), 
        slope_pars('post_isolation', 'beta_I_isolation') |> 
          mutate(predictor = 'Island  \nisolation'), 
        slope_pars('post_size', 'beta_I_size') |> 
          mutate(predictor = 'Island\nsize  '), 
        slope_pars('post_altitude', 'beta_I_alt') |> 
          mutate(predictor = 'Island \naltitude'), 
        slope_pars('post_footprint', 'beta_H_foot') |> 
          mutate(predictor = 'Human \nfootprint'), 
        slope_pars('post_nativeV', 'beta_NV') |> 
          mutate(predictor = 'Native   \nvegetation'), 
        slope_pars('post_bush', 'beta_bush') |> 
          mutate(predictor = 'Bush\ncover')) |> 
  ggplot(aes(predictor, mu, ymin = li, ymax = ls, 
             color = `Fruit consumption`)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.4), 
                linewidth = 1.5, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  scale_color_manual(values = c('#D92525', '#F2C230', '#7A577A')) +
  annotate("segment",
           x = c(1.13, 
                 4.13, 
                 4, 
                 5-0.13, 
                 5+0.13, 
                 3+0.13, 
                 6-0.13), 
           xend = c(1.13,
                    4.13, 
                    4, 
                    5-0.13, 
                    5+0.13, 
                    3+0.13, 
                    6-0.13), # y
           y = c(-0.8, 
                 -0.35, 
                 -1.1, 
                 -0.65, 
                 -0.83, 
                 -0.87, 
                 -0.73), 
           yend = c(-0.5, 
                    -0.05, 
                    -0.8, 
                    -0.35, 
                    -0.53, 
                    -0.57, 
                    -0.43), # x
           arrow = arrow(type = "closed", 
                         length = unit(1.5, "mm"),
                         angle = 30),
           color = c('#7A577A', 
                     '#7A577A', 
                     '#F2C230', 
                     '#D92525', 
                     '#7A577A', 
                     '#7A577A',
                     '#D92525'),
           linewidth = 0.5) +
  labs(y = expression(beta), x = 'Predictors') +
  theme_classic() +
  coord_flip() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.background = element_blank(),
    legend.position = c(0.85, 0.1),
    legend.box.background = element_blank(), 
    legend.key.size = unit(4, 'mm'), 
    text = element_text(family = 'Times New Roman', size = 9)
  )


plot_grid(plot_beta_continuous, 
          plot_grid(plot_scatter_plot1, 
                    plot_grid(NULL, 
                              plot_scatter_plot2, 
                              NULL, 
                              nrow = 1, 
                              rel_widths = c(0.2, 0.8, 0.2)), 
                    nrow = 2), 
          ncol = 1, 
          rel_heights = c(0.6, 0.4), 
          labels = c('(a)', '(b)'), 
          label_fontface = 'plain', 
          label_fontfamily = 'Times New Roman', 
          label_y = c(1, 1.05))

ggsave('effects_plot.jpg', width = 10, height = 20, units = 'cm', dpi = 1e3)


sessionInfo()
