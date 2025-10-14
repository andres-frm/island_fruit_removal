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

# =============== all_frugivory 







