# ======== Loading packages =======

sapply(c('cmdstanr', 'readxl', 'magrittr', 'dplyr', 'ggplot2', 
         'tidyr', 'tibble', 'forcats', 'rethinking', 
         'cowplot'), 
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

rbind(pivot_longer(post_latitude_tot$beta, 'beta_lat') |> 
        mutate(type = 'Frugivory', 
               effect = 'Latitude'), 
      pivot_longer(post_latitude_disp$beta, 'beta_lat') |> 
        mutate(type = 'Seed dispersion', 
               effect = 'Latitude'), 
      pivot_longer(post_latitude_pred$beta, 'beta_lat') |> 
        mutate(type = 'Seed predation', 
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
  lims(y = c(0, 0.75)) +
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








