sapply(c('readxl', 'magrittr', 'dplyr', 'ggplot2', 
         'tidyr', 'tibble', 'forcats', 'lubridate', 
         'maps', 'mapdata', 'geosphere'), 
       library, character.only = T)

source('functions_mod_diagnostics.r')

# ============ Loading data ====

d <- read_xlsx('ilhas2025.xlsx', sheet = 1, col_names = T, na = 'NA')

str(d)

summary(d)

apply(d, 2, function(x) mean(is.na(x)))

d$country2 <- d$country

d$country2[which(d$country2 == 'US Virgin Islands')] <- 'USA'

# =============== Map ===============

coords_islands <- unique(d[, c("island", "long", 'lat', "realm", 
                               "ecoregion", "biome", "island_type", 
                               "country2")])


coords_islands <- lapply(split(coords_islands, 
                               coords_islands$island), function(x) x[1, ])

coords_islands <- do.call('rbind', coords_islands)

lapply(unique(coords_islands$realm), FUN = 
         function(x) {
           paste(x, '=', mean(coords_islands$realm == x) * 100)
         })

lapply(unique(coords_islands$country2), FUN = 
         function(x) {
           paste(x, '=', mean(coords_islands$country2 == x) * 100)
         })

coords_islands$realm <- as.factor(coords_islands$realm)
coords_islands$real_code <- as.numeric(coords_islands$realm)
levels(coords_islands$realm)

#jpeg('map.jpeg', width = 20, height = 20, units = 'cm', res = 500)
par(mar = c(1, 1, 1, 1))
map("world", col="gray90", fill=TRUE, bg="white", border="gray50")
coords_islands %$% points(long, lat, col = coords_islands$real_code)
legend(x = -50, y = 190, legend = levels(coords_islands$realm), 
       lty = 1, col = 1:length(levels(coords_islands$realm)), 
       cex = 0.8)
#dev.off()

# ======== distance among islands =====

coords_islands

# shortest distance over the Earth's surface (great-circle distance).
islands_dist <- distm(coords_islands[, c("long", 'lat')], fun = distHaversine)
dimnames(islands_dist) <- list(coords_islands$island, 
                               coords_islands$island)

islands_dist[1:5, 1:5]

islands_dist <- islands_dist/max(islands_dist)

# ======== Distance among grids ======== 

coords_grids <- unique(d[, c("country2", "grid", "long", 'lat', "realm", 
                             "ecoregion", "biome", "island_type")])


coords_grids <- lapply(split(coords_grids, 
                             coords_grids$grid), function(x) x[1, ])

coords_grids <- do.call('rbind', coords_grids)

grid_dist <- distm(coords_grids[, c("long", 'lat')], fun = distHaversine)
dimnames(grid_dist) <- list(coords_grids$grid, 
                               coords_grids$grid)

grid_dist <- grid_dist/max(grid_dist)

grid_dist[1:5, 1:5]

# ========= Data cleaning ===========

d$fruit_fate[which(d$fruit_fate == 'eaten')] <- 'Eaten'

unique(d[, c("bush_ID", "fruit_ID")])

unique(d$bush_ID)
d[, c("fruit_fate", "fate_bin")]

colnames(d)

d$plant_ID <- d %$% paste0(grid, '--', bush_ID)

d <- 
  d[, c("country2", "island", "correct_coord", 'lat',
      'long', "realm", "ecoregion", "biome", "island_type", 
      "grid", "plant_ID", "fruit_ID", "fate_bin", 
      "frugivore", "season", "plant_invasive_rank", "human_pop", 
      "human_footprint", "distance_mainland", "island_size", "isolation", 
      "native_vegetation")]

plot(density(d$distance_mainland))

str(d)

d[] <- lapply(d, function(x) if (is.character(x)) as.factor(x) else x)

d <- split(d, d$plant_ID)

# islands where plant ID has to be corrected
d[which(unlist(lapply(d, nrow), use.names = F) != 15)]

d$`FN01--2` <- 
  rbind(d$`FN01--2`, d$`FN01--4`[1, ])

d$`FN01--2`$plant_ID[15] <- 'FN01--2'
d$`FN01--2`
(d$`FN01--4` <- d$`FN01--4`[-1, ])


d$`FN02--2` <- 
  rbind(d$`FN02--2`, d$`FN02--4`[1, ])
d$`FN02--2`$plant_ID[15] <- 'FN02--2'
d$`FN02--2`
(d$`FN02--4` <- d$`FN02--4`[-1, ])

d[which(unlist(lapply(d, nrow), use.names = F) != 15)]


# ===== Data: total and per-group fruit removal =========

d_total <- 
  lapply(d, FUN = 
         function(x) {
           remotion <- sum(x$fate_bin) # total remoremovaltion
           
           v <- 
             sapply(levels(x$frugivore), FUN = 
                    function(i) {
                      # removal per group
                      sum(x$frugivore == i)
                    })
           v1 <- as_tibble(matrix(v, nrow = 1))
           colnames(v1) <- names(v)
           
           # removing unnecessary grouping factors
           indx <- c(grep('fate_bin', colnames(x)), 
                     grep('frugivore', colnames(x)), 
                     grep('fruit_ID', colnames(x)))
           x$lat <- x$lat[1]
           x$long <- x$long[1]
           x <- unique(x[, -c(indx)])
           x$total <- 15 # total available fruits
           # adding total and per-group fruit removal
           x$total_remotion <- remotion 
           as_tibble(cbind(x, v1))
         })

d_total[unlist(lapply(d_total, function(x) nrow(x) > 1), use.names = F)]

d_total <- do.call('rbind', d_total)

unique(d_total$plant_invasive_rank)

apply(d_total, 2, FUN = function(x) mean(is.na(x)))

colnames(d_total)

# fixing island label 
d_total$island_type[which(d_total$island == 'Sicily' &
                            d_total$island_type == 'Coraline')] <- 'Continental'

codes1 <- d_total[, c("country2", "island", "grid", "plant_ID", 
                     "realm", "ecoregion", "biome", "island_type")]

codes <- codes1

codes[] <- lapply(codes1, function(x) as.numeric(as.factor(x)))

codes_labes <- 
  list(country = 
       unique(tibble(country = codes1$country2, 
                     code = codes$country2)), 
     island = 
       unique(tibble(island = codes1$island, 
                     code = codes$island)), 
     grid = unique(tibble(island = codes1$grid, 
                          code = codes$grid)), 
     real = unique(tibble(island = codes1$realm, 
                          code = codes$realm)), 
     ecoregion = unique(tibble(island = codes1$ecoregion, 
                               code = codes$ecoregion)), 
     biome = unique(tibble(island = codes1$biome, 
                           code = codes$biome)), 
     island_type = unique(tibble(island = codes1$island_type, 
                                 code = codes$island_type)))

# ====== Data structure for generative simulation ======

saveRDS(list(data_structure = codes,
             labels = codes1,
             dist_mat = islands_dist), 'data_structure_simulation.rds')



# ======= Data for models ========

saveRDS(list(codes = codes_labes,
             dist_islands = islands_dist, 
             data = d_total), 'data.rds')



