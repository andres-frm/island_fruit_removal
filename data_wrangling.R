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

# ===== Data 1: total removal =========

d_total <- d

d_total <- split(d_total, d_total$plant_ID)

d_total[which(unlist(lapply(d_total, nrow), use.names = F) != 15)]

d_total$`FN01--2` <- 
  rbind(d_total$`FN01--2`, d_total$`FN01--4`[1, ])

d_total$`FN01--2`$plant_ID[15] <- 'FN01--2'
d_total$`FN01--2`
(d_total$`FN01--4` <- d_total$`FN01--4`[-1, ])


d_total$`FN01--2` <- 
  rbind(d_total$`FN01--2`, d_total$`FN01--4`[1, ])

d_total$`FN01--2`$plant_ID[15] <- 'FN01--2'
d_total$`FN01--2` 


d_total$`FN02--2` <- 
  rbind(d_total$`FN02--2`, d_total$`FN02--4`[1, ])
d_total$`FN02--2`$plant_ID[15] <- 'FN02--2'
d_total$`FN02--2`
(d_total$`FN02--4` <- d_total$`FN02--4`[-1, ])

d_total[which(unlist(lapply(d_total, nrow), use.names = F) != 15)]

length(d_total)

