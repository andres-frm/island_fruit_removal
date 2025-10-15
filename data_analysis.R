sapply(c('cmdstanr', 'readxl', 'magrittr', 'dplyr', 'ggplot2', 
         'tidyr', 'tibble', 'forcats', 'rethinking', 
         'cowplot'), 
       library, character.only = T)

source('functions_mod_diagnostics.r')

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

# =============== Macroecological processes ===========

# =============== Effects of latitude ==========

file <- paste0(getwd(), '/mod_latitud_total.stan')
fit_latitude_tot <- cmdstan_model(file, compile = T)

# =============== all frugivores  ======================

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


post_altitude_tot <- 
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

post_altitude_tot <- 
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
             post_altitude_tot[, grep(x, colnames(post_altitude_tot))]
           })

names(post_altitude_tot) <- c('alpha', 
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


# ======== functions for extracting effects =======

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


# est_latitude_tot <- cond_effects(posterior = post_altitude, 
#                              x_bar = dat$lat, 
#                              slope = 'beta_lat', 
#                              type = 'random', 
#                              n = 100)
# 
# plot(dat$lat, dat$total_remotion, cex = 0.1)
# est_latitude_tot %$% lines(x, y)
# est_latitude_tot %$% lines(x, li, lty = 3)
# est_latitude_tot %$% lines(x, ls, lty = 3)
# 
# est_latitude_tot %$% plot(x, y, type = 'l', ylim = c(0, 1))
# est_latitude_tot %$% lines(x, li, lty = 3)
# est_latitude_tot %$% lines(x, ls, lty = 3)
# 
# plot(NULL, xlim = c(-3, 3), ylim = c(0, 0.8))
# for(i in seq_along(post_altitude$p_realm)) {
#   lines(density(post_altitude$p_realm[[i]]), col = i)
#}

# ======= Plots =====
# 
# ============= Slops ===========

rbind(pivot_longer(post_altitude_tot$beta, 'beta_lat') |> 
        mutate(type = 'Frugivory', 
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
  geom_hline(yintercept = 0, linetype = 3)


# ============= type island ======


TI_tot <- average_effects(n_levels = ncol(post_altitude_tot$TI),
                          posterior = post_altitude_tot, 
                          x_var1 = 'island_type', 
                          x_var2 = 'realm', 
                          par1 = 'TI', 
                          par2 = 'p_realm')


TI_tot <- full_join(TI_tot, codes$island_type, 'code')


rbind(TI_tot |> 
        mutate(type = 'Frugivory')) |> 
  group_by(type, island, x) |> 
  transmute(mu = median(y), 
            li = quantile(y, 0.025), 
            ls = quantile(y, 0.975),
            type = type,
            island = island
            ) |> 
  unique() |> 
  ggplot(aes(island, mu, ymin = li, ymax = ls)) +
  geom_errorbar(width = 0) +
  geom_point() +
  facet_wrap(~type) +
  lims(y = c(0, 1)) +
  labs(y = 'P(fruit consumption)', x = 'Type of island')
