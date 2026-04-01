impute_cov_P1 <- function(cov.P1, GridBAUs.Points, data_popolazione_smoothed_nodes) {
  # 1. Trova gli indici con valori mancanti
  idx_na <- cov.P1 %>% is.na() %>% which()

  # Se non ci sono NA, restituisco subito
  if (length(idx_na) == 0) {
    return(cov.P1^0.5)
  }

  # 2. Seleziona le coordinate corrispondenti
  coords_na <- as.data.frame(coordinates(GridBAUs.Points))[idx_na, ]

  # 3. Crea sf POINT
  points_sf <- st_as_sf(as.data.frame(coords_na),
                        coords = c("x", "y"),
                        crs = st_crs(GridBAUs.Points)) %>%
    st_transform(st_crs(data_popolazione_smoothed_nodes))

  # 4. Trova il punto più vicino in data_popolazione_smoothed_nodes
  nearest_idx <- st_nearest_feature(points_sf, data_popolazione_smoothed_nodes)

  # 5. Estrai i valori corrispondenti
  val_nearest <- data_popolazione_smoothed_nodes$density[nearest_idx]

  # 6. Riempie i valori mancanti
  cov.P1[idx_na] <- val_nearest

  # 7. Applica trasformazione
  cov.P1.trans <- cov.P1^0.5

  return(cov.P1.trans)
}


impute_cov_Elev <- function(cov.Elev, GridBAUs.Points, GridBAUs1) {
  library(terra)
  spdf.cov.Elev <- SpatialPixelsDataFrame(points = GridBAUs.Points, data = data.frame(value = cov.Elev))
  r <- rast(spdf.cov.Elev["value"])
  filled <- focal(r, w = matrix(1, 5, 5), fun = mean, na.policy = "only", fillvalue = NA)
  r_raster <- raster(filled)
  cov.Elev <- Interp_OnBaus_Raster(GridBAUs1, r_raster)
  cov.Elev
}


