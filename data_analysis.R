sapply(c('cmdstanr', 'readxl', 'magrittr', 'dplyr', 'ggplot2', 
         'tidyr', 'tibble', 'forcats', 'rethinking', 
         'cowplot'), 
       library, character.only = T)

source('functions_mod_diagnostics.r')

data <- readRDS('data.rds')

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

# fruit predation
dat$rodent_all <- ifelse(dat$Unknown == 15, 1, 0)
dat$rodent_all <- ifelse((dat$rodent_all + dat$Rodent) > 0, 1, 0)

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

file <- paste0(getwd(), '/mod_latitude.stan')
fit_latitude <- cmdstan_model(file, compile = T)

# =============== all frugivores 

mod_latitude <- 
  fit_latitude$sample(
    data = dat, 
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3, 
    seed = 23061993
  )

sum_latitude <- mod_latitude$summary()
mod_diagnostics(mod_latitude, sum_latitude)
ppcheck_latitude <- mod_latitude$draws('ppcheck', format = 'matrix')

plot(density(dat$total_remotion), main = '', 
     xlab = 'Total fruits removal', ylim = c(0, 0.1))
for (i in 1:200) lines(density(ppcheck_latitude[i, ], lwd = 0.1))
lines(density(dat$total_remotion), lwd = 2, col = 'red')


post_altitude <- 
  mod_latitude$draws(c('alpha', 
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

post_altitude <- 
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
             post_altitude[, grep(x, colnames(post_altitude))]
           })

names(post_altitude) <- c('alpha', 
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


est_latitude <- cond_effects(posterior = post_altitude, 
                             x_bar = dat$lat, 
                             slope = 'beta_lat', 
                             type = 'random', 
                             n = 100)

plot(dat$lat, dat$total_remotion, cex = 0.1)
est_latitude %$% lines(x, y)
est_latitude %$% lines(x, li, lty = 3)
est_latitude %$% lines(x, ls, lty = 3)


est_latitude %$% plot(x, y, type = 'l', ylim = c(0, 1))
est_latitude %$% lines(x, li, lty = 3)
est_latitude %$% lines(x, ls, lty = 3)





