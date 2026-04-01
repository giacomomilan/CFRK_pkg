# CFRK Utilities: Helper functions for CFRK spatial prediction
library(magrittr)
library(compositions)
library(lubridate)
library(FRK)
library(sf)
library(ggplot2)

# Function to create a bounding box from coordinates
# Parameters:
# - crs: Coordinate reference system
# - xmin, xmax, ymin, ymax: Bounding box limits
GetBBox <- function(crs, xmin=0, xmax=1, ymin=0, ymax=1) {
  st_as_sfc(st_bbox(
    c(xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax), crs = crs
    )
  )
}

# Function to create a bounding box from an sf geometry
# Parameters:
# - sf_geometry: sf object with geometry
GetBBox_FromGeometry <- function(sf_geometry) {
  library(sf)
  limits = unlist(st_bbox(sf_geometry))
  sfbbox = GetBBox(st_crs(sf_geometry),
      xmin=limits[[1]], xmax=limits[[3]],
      ymin=limits[[2]], ymax=limits[[4]])
  sfbbox
}

# Function to create a spatial grid from a bounding box
# Parameters:
# - this_bbox: Bounding box sf object
# - xnum, ynum: Number of grid cells in x and y directions
GetSpatialGrid_FromBBox <- function(this_bbox, xnum=100, ynum=100) {
  xmin = min(st_coordinates(this_bbox)[,1])
  xmax = max(st_coordinates(this_bbox)[,1])
  ymin = min(st_coordinates(this_bbox)[,2])
  ymax = max(st_coordinates(this_bbox)[,2])

  xstep = (xmax - xmin) / xnum
  ystep = (ymax - ymin) / ynum

  x = GridTopology(c(xmin+0.5*xstep,ymin+0.5*ystep), c(xstep,ystep), c(xnum-1,ynum-1))
  xlarge = GridTopology(c(xmin-xstep,ymin-ystep), c(xstep,ystep), c(xnum+2,ynum+2))
  y = SpatialGrid(grid = xlarge, proj4string = CRS(paste("EPSG:", st_crs(this_bbox)$epsg)))
  coordnames(y) <- c('x', 'y')
  y
}

# Algebraic methods

# Function to interpolate a vector on a regular grid
# Parameters:
# - nodes.coords: Coordinates of nodes
# - nodes.v: Values at nodes
# - eval.coords: Coordinates for evaluation
# - extrap: Extrapolation flag
InterpThisVector <- function(nodes.coords, nodes.v, eval.coords, extrap = FALSE){
  library(fields)
  library(akima)
  # Linear interpolation on regular grid
  interp_res <- interp(x = nodes.coords[,1], y = nodes.coords[,2], z = nodes.v,
                      linear = TRUE, extrap = extrap)
  # Interpolation at arbitrary points
  valori_interp <- interp.surface(interp_res, eval.coords)
  valori_interp
}

# PCA Transform

# Function to transform variables to PCA scores
# Parameters:
# - dataframe: Data frame with variables
# - colnames: Column names for PCA
VarsToScores <- function(dataframe, colnames) {
  pca_output <- prcomp(as.data.frame(dataframe)[,colnames])
  list(scores=pca_output$x, pca_output=pca_output)
}

# Function to transform PCA scores back to variables
# Parameters:
# - dataframe_scores: Data frame with PCA scores
# - pca_result: PCA result object
# - colnames: Original column names
ScoresToVars <- function(dataframe_scores, pca_result, colnames) {
  original_data <- as.data.frame(dataframe_scores %*% t(pca_result$rotation) +
                                   matrix(rep(pca_result$center,
                                              nrow(dataframe_scores)), nrow = nrow(dataframe_scores), byrow = TRUE)) %>% set_colnames(colnames)
  original_data
}

# ILR Transform

# Function to transform compositional data to ILR coordinates
# Parameters:
# - df_cx: Compositional data
# - cnames: Component names
# - df_stations: Station data
# - domain: Spatial domain
TransformILR <- function(df_cx, cnames, df_stations, domain) {
  colnames <- sapply(1:(length(cnames) - 1), FUN=function(i){paste0("Y", as.character(i))})
  df_comp = ilr(df_cx[,cnames]) %>% set_colnames(colnames)

  df1 <- df_comp %>% as.data.frame %>% tibble::rownames_to_column("id")
  df2 <- df_stations %>% tibble::rownames_to_column("id")
  df_comp <- left_join(df1, df2, by = "id")

  df_coord <- st_as_sf(df_comp)
  st_crs(df_coord) <- st_crs(df_stations)
  df_coord <- st_transform(df_coord, st_crs(domain))

  # Add inversion columns (normalized by sum)
  df_comp <- as.data.frame(setNames(ilrInv(ilr(df_cx[,cnames])), paste0(cnames, "_perc"))) %>% cbind(df_comp)
  #coordinates(df_comp) <- c('Longitude','Latitude') # mimic SpatialPointsDataFrame
  df_comp_out = SpatialPointsDataFrame(data = as.data.frame(df_comp),
    coords = st_coordinates(df_coord), proj4string= CRS(paste("EPSG:", st_crs(domain)$epsg)))
  coordnames(df_comp_out) <- c('x', 'y')
  df_comp_out
}

# Function to revert ILR transformation
# Parameters:
# - dataframe: ILR data
# - x, y: Coordinates
RevertILR <- function(dataframe, x, y) {
 as.data.frame(matrix(as.numeric(ilrInv(dataframe)), ncol=dim(dataframe)[[2]]+1)) %>%
    set_colnames(paste("c",1:(dim(dataframe)[[2]]+1), sep='')) %>%
    cbind(x=x, y=y)
}

# Logit Transform

# Logit function
logit <- function(p) {
  log(p / (1 - p))
}

# Inverse logit function
inverse_logit <- function(x) {
  1 / (1 + exp(-x))
}

# ILR for two components
ilr_two_comp <- function(x1, x2) {
  sqrt(1/2) * log(sqrt(x1*x2)/x2)
}

# Function to transform data using logit
# Parameters:
# - df_cx: Data frame
# - cnames: Component names
# - df_stations: Station data
# - domain: Spatial domain
# - data_p: Probability data
TransformLogit <- function(df_cx, cnames, df_stations, domain, data_p) {
  #colnames <- sapply(1:(length(cnames) - 1), FUN=function(i){paste0("Y", as.character(i))})
  #df_comp = ilr(df_cx[,cnames]) %>% set_colnames(colnames)
  df_comp <- cbind(data_p, logit=log(data_p$p / (1 - data_p$p)))

  df1 <- df_comp %>% as.data.frame %>% tibble::rownames_to_column("id")
  df2 <- df_stations %>% tibble::rownames_to_column("id")
  df_comp <- left_join(df1, df2, by = "id")

  df_coord <- st_as_sf(df_comp)
  st_crs(df_coord) <- st_crs(df_stations)
  df_coord <- st_transform(df_coord, st_crs(domain))

  # Add inversion columns (normalized by sum)
  #df_comp <- data.frame(logit=log(p / (1 - p))) %>% cbind(df_comp)
  #coordinates(df_comp) <- c('Longitude','Latitude') # mimic SpatialPointsDataFrame
  df_comp_out = SpatialPointsDataFrame(data = as.data.frame(df_comp),
    coords = st_coordinates(df_coord), proj4string= CRS(paste("EPSG:", st_crs(domain)$epsg)))
  coordnames(df_comp_out) <- c('x', 'y')
  df_comp_out
}

# Function to revert logit transformation
# Parameters:
# - logit: Logit values
# - x, y: Coordinates
RevertLogit <- function(logit, x, y) {
 as.data.frame(cbind(1-inverse_logit(logit), inverse_logit(logit))) %>%
    set_colnames(c("c1", "c2")) %>%
    cbind(x=x, y=y)
}

# FRK objects

# Function to get a spatial grid within a domain
# Parameters:
# - domain.geometry: Domain geometry
# - grid.xnum, grid.ynum: Grid dimensions
GetGrid <- function(domain.geometry, grid.xnum=80, grid.ynum = 80) {
  library(sp)
  sfbbox = GetBBox_FromGeometry(domain.geometry)
  grid = GetSpatialGrid_FromBBox(sfbbox, xnum=grid.xnum, ynum=grid.ynum)
  grid.df = data.frame(grid)
  grid = SpatialPoints(
    grid.df[unlist(st_contains(domain.geometry, st_as_sf(x=grid.df, coords=c('x','y'), crs = st_crs(domain.geometry)))),],
    proj4string = CRS(paste("EPSG:", st_crs(sfbbox)$epsg)))
  grid
}

# Function to create BAU grid
# Parameters:
# - dataframe: Data for BAU creation
# - cellsize: Cell size for BAUs
GetGridBAUs <- function(dataframe, cellsize) {
  set.seed(1)
  # 1. Set of BAUs
  GridBAUs1 <- auto_BAUs(manifold = plane(), # 2D plane
                          cellsize = c(cellsize, cellsize), # BAU cellsize
                          type = "grid", # grid (not hex)
                          data = dataframe, # data around which to create BAUs
                          convex = -0.05, # border buffer factor
                          nonconvex_hull = FALSE) # convex hull

  # 2. BAU-specific variation
  # Variation at BAU level for heteroscedasticity:
  # For the ith BAU, we also need to attribute the element vi that describes
  # the heteroscedasticity of the fine-scale variation for that BAU
  # (can be typically altitude)
  GridBAUs1$fs <- 1 # fine-scale variation at BAU level
  GridBAUs1
}

# Function to subset intersecting grid with polygons
# Parameters:
# - GridBAUs1: BAU grid
# - polygons_sf_cover: Covering polygons
Subset_Intersecating_Grid <- function(GridBAUs1, polygons_sf_cover){
  # Assuming spdf is a SpatialPixelsDataFrame
  pixels_sf <- st_as_sf(as(GridBAUs1, "SpatialPointsDataFrame"))
  polygons_sf <- st_make_valid(polygons_sf_cover)  # if necessary
  # Calculate intersections
  hits <- st_intersects(pixels_sf, polygons_sf, sparse = TRUE)
  # hits is a list of vectors: for each pixel, the indices of polygons it intersects
  intersecting_pixels <- lengths(hits) > 0
  GridBAUs1[intersecting_pixels,]
}

# Function to get FRK output folder
# Parameters:
# - figures_folder: Base folder
# - DestName.Pollutant: Pollutant name
# - Global_Use_Cov_Population: Use population covariate
# - Global_Use_Cov_Elevation: Use elevation covariate
# - frk.args: FRK arguments
# - flag_use_pca: PCA flag
# - flag_obs_fs: Observation fs flag
# - baus.cellsize: BAU cell size
# - DestName.IndexDrop: Index drop
# - foldername: Folder name (optional)
GetFRKFolderOutput <- function(
  figures_folder, DestName.Pollutant,
  Global_Use_Cov_Population, Global_Use_Cov_Elevation,
  frk.args, flag_use_pca, flag_obs_fs, baus.cellsize, DestName.IndexDrop, foldername = "") {
  conf = paste0(figures_folder, "/FRK_Pollutant_", DestName.Pollutant,
    "_Global_Use_Cov_Population_", Global_Use_Cov_Population,  "_Global_Use_Cov_Elevation_", Global_Use_Cov_Elevation,
    "_prediction_levels", as.character(frk.args["basis_nres"]) , "_", frk.args["basis_type"],
    "_PCA_", as.character(flag_use_pca), "_obsFs_", as.character(flag_obs_fs), "_BausCellsize_", baus.cellsize,
    "_IndexDrop_", as.character(DestName.IndexDrop))
  if (nchar(foldername) == 0) {
    foldname <- format(Sys.time(), "%Y%m%d_%H%M%S")
    foldname <- "output"
  }
  folder_output = paste(figures_folder, foldname, sep='/')
  if(!dir.exists(folder_output)){dir.create(folder_output)}
  # Write to file "output.txt"
  cat(conf, file = paste(folder_output, 'configuration.txt', sep = '/'))
  folder_output
}

# Function to get basis functions for FRK
# Parameters:
# - dataframe: Data for basis creation
# - frk.args: FRK arguments
GetBasis <- function(dataframe, frk.args) {
  if ("basis_type" %in% names(frk.args)) {
    type = frk.args[["basis_type"]]
  } else {type = "Gaussian"}
  if ("type" %in% names(frk.args)) {
    type = frk.args[["type"]]
  }
  if ("basis_nres" %in% names(frk.args)) {
    nres = frk.args[["basis_nres"]]
  } else {nres = 3}
  if ("nres" %in% names(frk.args)) {
    nres = frk.args[["nres"]]
  }
  if ("regular" %in% names(frk.args)) {
    regular = frk.args[["regular"]]
  } else {regular = 0}
  cat(">> AUTO BASIS: using Type =", type,'\n')
  cat(">> AUTO BASIS: using Nres =", nres,'\n')
  cat(">> AUTO BASIS: using Regular =", regular,'\n')

  # 3. Basis functions for covariance modeling
  G <- auto_basis(manifold = plane(), # 2D plane
                  data = dataframe, # meuse data
                  nres = nres, # number of resolutions
                  type = type, # type of basis function
                  regular = regular) # place regularly in domain
  G
}

# Function to interpolate raster values on BAUs
# Parameters:
# - GridBAUs: BAU grid
# - raster: Raster object
Interp_OnBaus_Raster <- function(GridBAUs, raster) {
  library(raster)
  grid_r <- st_transform(st_as_sf(as.data.frame(coordinates(GridBAUs)),
                                  coords = c('x', 'y'), crs = st_crs(GridBAUs)),
                         crs = st_crs(raster))
  # Covariate conversion for method
  cov = raster::extract(raster, grid_r)
  cov
}

# Function to interpolate population density on BAUs
# Parameters:
# - data_air: Air data
# - GridBAUs1: BAU grid
# - data_popolazione_smoothed_nodes: Population data
Interp_OnBaus_DensityPopOnNodes <- function(data_air, GridBAUs1, data_popolazione_smoothed_nodes) {
  cov.P1 <- InterpThisVector(st_coordinates(data_popolazione_smoothed_nodes),
       data_popolazione_smoothed_nodes$density, as.data.frame(GridBAUs1)[,c('x', 'y')])
  cov.P1
}

# Function to interpolate elevation on BAUs
# Parameters:
# - data_air: Air data
# - GridBAUs1: BAU grid
# - raster.elevazione: Elevation raster
Interp_OnBaus_Elevation <- function(data_air, GridBAUs1, raster.elevazione) {
  grid_r <- st_transform(st_as_sf(as.data.frame(coordinates(GridBAUs1)),
    coords = c('x', 'y'), crs = st_crs(data_air$grid)),
               crs = st_crs(raster.elevazione))
  # Covariate conversion for method
  cov.Elev = raster::extract(raster.elevazione, grid_r)
  cov.Elev
}

# Core FRK functions

# Function to fit and predict with formula
# Parameters:
# - formula: Model formula
# - dataframe: Data
# - GridBAUs: BAU grid
# - BasisG: Basis functions
# - pars.fit.no_iterations: Number of EM iterations
# - tol: Tolerance
FitAndPredict_withFormula <- function(formula, dataframe, GridBAUs, BasisG, pars.fit.no_iterations, tol=0.01) {
    library(FRK) # for carrying out FRK
    boot_data <- as.data.frame(dataframe)
    row.names(boot_data) <- 1:nrow(boot_data)
    coordinates(boot_data) <- ~ x + y
    proj4string(boot_data) <- proj4string(dataframe)
    S <- SRE(f = as.formula(formula), # formula
                data = list(boot_data), # list of datasets
                #method="EM",
                BAUs = GridBAUs, # BAUs
                basis = BasisG, # basis functions
                est_error = TRUE, # estimation measurement error
                average_in_BAU = FALSE) # do not average data over BAUs

    # EM method
    S <- SRE.fit(S, # SRE model
                    n_EM = pars.fit.no_iterations, # max. no. of EM iterations
                    tol = tol, # tolerance at which EM is assumed to have converged
                    print_lik=TRUE) # print log-likelihood at each iteration

    BAUs_prediction <- predict(S, obs_fs = flag_obs_fs)
    BAUs_prediction_df <- as(BAUs_prediction,"data.frame")
    list(BAUs_prediction_df=BAUs_prediction_df, S=S)
}

# Function to fit and predict with formula and K type
# Parameters:
# - formula: Model formula
# - dataframe: Data
# - GridBAUs: BAU grid
# - BasisG: Basis functions
# - pars.fit.no_iterations: Number of EM iterations
# - pars.fit.Ktype: K type
FitAndPredict_withFormula_K <- function(formula, dataframe, GridBAUs, BasisG, pars.fit.no_iterations, pars.fit.Ktype) {
    cat("K_type ", pars.fit.Ktype, "\n")
    library(FRK) # for carrying out FRK
    S <- SRE(f = as.formula(formula), # formula
                data = list(dataframe), # list of datasets
                #method="EM",
                BAUs = GridBAUs, # BAUs
                basis = BasisG, # basis functions
                est_error = TRUE, # estimation measurement error
                average_in_BAU = FALSE) # do not average data over BAUs

    # EM method
    S <- SRE.fit(S, # SRE model
                    n_EM = pars.fit.no_iterations, # max. no. of EM iterations
                    tol = 0.01, # tolerance at which EM is assumed to have converged
                    print_lik=TRUE) # print log-likelihood at each iteration

    BAUs_prediction <- predict(S, obs_fs = flag_obs_fs)
    BAUs_prediction_df <- as(BAUs_prediction,"data.frame")
    list(BAUs_prediction_df=BAUs_prediction_df, S=S)
}

# FRK steps functions

# Function to prepare data for kriging
# Parameters:
# - data.comp: Compositional data
# - colnames.euclid: Euclidean column names
# - flag_use_pca: PCA flag
PrepareData_For_Kriging <- function(data.comp, colnames.euclid, flag_use_pca) {
  if (flag_use_pca)
  {
    cat('YES pca\n')
    pca_result <- VarsToScores(as.data.frame(data.comp)[,colnames.euclid], colnames.euclid)
    pca_output <- pca_result$pca_output

    data.krig = SpatialPointsDataFrame(
      coordinates(data.comp),
      as.data.frame(pca_result$scores) %>%
        set_colnames(colnames.euclid),
      proj4string = CRS(paste("EPSG:",st_crs(data.comp)$epsg)))
  } else {
    cat('NO pca\n')
    data.for.spdf = as.data.frame(as.data.frame(data.comp)[, colnames.euclid])
    colnames(data.for.spdf) <- colnames.euclid
    data.krig = SpatialPointsDataFrame(
      coordinates(data.comp),
      data.for.spdf,
      proj4string = CRS(paste("EPSG:", st_crs(data.comp)$epsg)), match.ID = FALSE)
    pca_output = {}
  }
  list(data.krig = data.krig, pca_output = pca_output)
}

# Function to get average values over municipalities
# Parameters:
# - dataframe.prediction: Prediction data
# - comuni: Municipal boundaries
# - field: Field to average
# - data_air_dataset: Air dataset
GetAvg <- function(dataframe.prediction, comuni, field, data_air_dataset) {
  library(sp)
  library(sf)
  library(raster)
  spdf = SpatialPixelsDataFrame(coordinates(dataframe.prediction)[, c('x', 'y')], data = as.data.frame(dataframe.prediction[[field]]))
  raster_obj <- raster(spdf)
  crs(raster_obj) <- CRS(paste("EPSG:", st_crs(data_air_dataset)$epsg))
  library(exactextractr)
  exact_extract(raster_obj, comuni, fun="mean")
}

# Function to perform FRK prediction
# Parameters:
# - domain: Spatial domain
# - GridBAUs: BAU grid
# - frk.args: FRK arguments
# - data.comp: Compositional data
# - data.krig: Kriging data
# - colnames.euclid: Euclidean column names
# - comuni: Municipal boundaries
# - Global_Use_Cov_Population: Use population covariate
# - Global_Use_Cov_Elevation: Use elevation covariate
# - figures_folder: Figures folder
# - no_iterations: Number of iterations
# - fitted.locations: Fitted locations
# - tol: Tolerance
DoFRKPrediction <- function(domain, GridBAUs, frk.args, data.comp, data.krig,
                            colnames.euclid, comuni,
                            Global_Use_Cov_Population, Global_Use_Cov_Elevation,
                            figures_folder,
                            no_iterations, fitted.locations = NULL, tol=0.01, flag_no_plot = TRUE) {
  colnames_pred = paste0(colnames.euclid, '.pred')

  if (is.null(GridBAUs)) {
    prediction.grid = data_air$grid
    GridBAUs1 <- GetGridBAUs(prediction.grid, cellsize=baus.cellsize)

    if (is.null(cov.P1)) {
      # Add population to GridBAUs
      layer_popolazione.Spat = as(st_transform(layer_popolazione, st_crs(GridBAUs1)), "Spatial")
      cov.P1 = over(SpatialPoints(GridBAUs1, proj4string = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")), layer_popolazione.Spat)$Popolazione.P1
    }

    grid_r <- st_transform(st_as_sf(as.data.frame(coordinates(GridBAUs1)),
      coords = c('x', 'y'), crs = st_crs(prediction.grid)),
                crs = st_crs(raster.elevazione))
    # Covariate conversion for method
    cov.Elev = raster::extract(raster.elevazione,
      grid_r)

    cat("cov.P1", length(cov.P1), "\n")
    cat("cov.Elev", length(cov.Elev), "\n")
    discard_tiles = is.na(cov.P1) | is.na(cov.Elev)
    cat("dim(GridBAUs1)", dim(GridBAUs1[!discard_tiles,]), "\n")

    GridBAUs = GridBAUs1[!discard_tiles,]

    GridBAUs$Elevazione = cov.Elev[!discard_tiles]
    GridBAUs$Popolazione.P1 = cov.P1[!discard_tiles]
  }

  if (!flag_no_plot) {
    coogrid = as.data.frame(coordinates(GridBAUs1))[ ,c('x', 'y')]
    ggplot() + geom_sf(data=domain, aes()) + geom_point(data=coogrid, aes(x, y))
    ggsave(paste(folder_output, "baus.png", sep="/"))
  }

  domain.prediction.grid <- GetGrid(domain) # TODO construct GetGridBAUs without domain.prediction.gridv
  G <- GetBasis(domain.prediction.grid, frk.args)
  cat("done basis G\n")

  # Show basis functions
  if (!flag_no_plot) {
    g <- ggplot() + geom_sf(data=domain, aes())
    show_basis(G, g = g) + # coord_fixed() + # fix aspect ratio
      xlab("Easting (m)") + # x-label
      ylab("Northing (m)") # y-label
    ggsave(paste(folder_output, "basis_functions.png", sep="/"))
  }

  # Output info
  cat("INFO: BASE\n")
  cat("Number of Basis Functions = ", attributes(G)$n, "\n\n")
  cat("Number of observations = ", dim(data.comp)[1], "\n\n")
  cat('\n')

  # ============================================================================
  # 1   FIT AND PREDICTION
  # 1.1 Prediction on BAUs

  n <- nrow(coordinates(GridBAUs))
  prediction <- SpatialPixelsDataFrame(
    coordinates(GridBAUs),
    data = data.frame(id = seq_len(n)),
    proj4string = CRS(paste0("EPSG:", st_crs(data.comp)$epsg))
  )

  Slist = list()
  # tutte le predizioni
  for (i in (1:length(colnames.euclid)))
  {
    formula = paste0(colnames.euclid[i], " ~ 1")
    if(Global_Use_Cov_Population)
    {
      formula = paste0(colnames.euclid[i], " ~ Popolazione.P1")
    }
    if((Global_Use_Cov_Population) & (Global_Use_Cov_Elevation))
    {
      formula = paste0(colnames.euclid[i], " ~ Popolazione.P1 + Elevazione")
    }
    cat("using formula in frk: ", formula, "\n")
    cat("using tolerance in frk: ", tol, "\n")

    if(TRUE) {
      x11(width = 8, height = 8)
      png(paste0(figures_folder, '/convergence_loglik_index_', i, '.png'))
    }

    if (exists("Flag_K_type_exponential") && Flag_K_type_exponential) {
      result <- FitAndPredict_withFormula_K(formula, data.krig, GridBAUs, G, pars.fit.no_iterations = no_iterations,
        pars.fit.Ktype = "block-exponential", tol = tol)
    } else {
      result <- FitAndPredict_withFormula(formula, data.krig, GridBAUs, G, pars.fit.no_iterations = no_iterations,
        tol = tol)
    }
    Slist = append(Slist, list(result))
    result.fit.pred = result
    dev.off()
    prediction[[paste0(colnames.euclid[i],".pred")]] = result$BAUs_prediction_df$mu
    prediction[[paste0(colnames.euclid[i],".sd")]] = result$BAUs_prediction_df$sd
  }

  geometries = st_transform(comuni$geometry, st_crs(GridBAUs1))

  comuniresultall.meaneucl = comuni

  # Aggregate predictions to municipal level
  for (i in (1:length(colnames_pred)))
  {
    datafill = GetAvg(prediction, comuni, colnames_pred[i], data.comp)
    comuniresultall.meaneucl[colnames_pred[i]] = datafill
  }

  # Prediction on observed locations (fitted)
  Extract <- function(prediction, fitted.locations, field) {
    library(raster)
    raster::extract(raster(SpatialPixelsDataFrame(
      coordinates(prediction),
      data = as.data.frame(prediction[[field]]))), fitted.locations )
  }


  if (is.null(fitted.locations)) {
    fitted.locations <- SpatialPoints(as.data.frame(data.comp)[c("x", "y")], proj4string = CRS(paste("EPSG:", st_crs(data.comp)$epsg)))
  }


  fitted.values = sapply(colnames_pred, FUN=function(field){Extract(prediction, fitted.locations, field)})
  fitted.values = fitted.values %>% matrix(nrow = length(fitted.locations)) %>% data.frame %>% set_colnames(colnames_pred)
  fitted.values = SpatialPointsDataFrame(coordinates(fitted.locations), data = fitted.values, proj4string = CRS(paste("EPSG:", st_crs(data.comp)$epsg)), match.ID = FALSE)

  list(
    predictionbaus.values = prediction,
    fitted.values = fitted.values,
    comuni.values = comuniresultall.meaneucl,
    Slist=Slist)
}

## Plot

PlotSymbols = function(dataset, colname, Do.Log = T, Use.Pop = T, title="", labels=c(), layer_popolazione=NULL) {
  plotsym_vals_shape = c(1, 1, 2, 5, 5)
  plotsym_vals_size = c(3, 1, 1, 1, 3)

  dataset$interval = as.integer(findInterval(dataset[[colname]], c(-Inf,unlist(quantile(dataset[[colname]], c(0.2,0.4,0.6, 0.8))), Inf)))
  dataset$size = plotsym_vals_size[dataset$interval]
  dataset$shape = as.factor(plotsym_vals_shape[dataset$interval])
  if (Use.Pop) {
    g <- ggplot() + geom_sf(data=layer_popolazione$geom, aes(fill=1e6 * layer_popolazione$Popolazione.P1 / as.numeric(st_area(layer_popolazione$geom))), color = "white", linewidth = 0.001, alpha=0.7) +
      scale_fill_distiller(palette = "Blues", direction=-1)
  } else {
    g <- ggplot() + geom_sf(data=regioni$geometry[3], aes())
  }
  g <- g + geom_point(data=dataset, aes(x, y), shape=dataset$shape, size=dataset$size, color="black", stroke=0.5) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust=0.5)) +
    guides(fill=guide_legend(title=""))
  if (length(labels) > 0) {
    g <- g + geom_text(aes(x=dataset$x, y=dataset$y, label = labels), vjust = -1)
  }
  g
}


CFRKPlotPrediction <- function(data.comp, domain, prediction.conc, comp.fitted.df, comp.comuni.conc,
    colnames.aitch, colnames.aitch.pred, DoPLOT) {
  # numbers
  prediction.conc.n = prediction.conc[, paste0("c", as.character(1:length(colnames.aitch)))] %>%
    set_colnames(colnames.aitch.pred)
  prediction.conc.n$x = prediction.conc$x
  prediction.conc.n$y = prediction.conc$y

  if (!is.null(comp.fitted.df)){
    comp.fitted.df.n = comp.fitted.df[, paste0("c", as.character(1:length(colnames.aitch)))] %>%
      set_colnames(colnames.aitch.pred)
    comp.fitted.df.n$x = comp.fitted.df$x
    comp.fitted.df.n$y = comp.fitted.df$y
  }

  if (DoPLOT) # PLOT su  BAUS
  {
    cat("Printing maps on BAUS\n")
    library(scales)  # Per la funzione label_number()


    x11(width = 8, height = 7)
    for (i in (1:length(colnames.aitch.pred)))
    {
      datafill = prediction.conc.n[[colnames.aitch.pred[i]]]
      datatitle = paste0("Prediction of interval I", as.character(i))
      ggpred =  ggplot() + geom_sf(data=domain$geom, aes()) +
        geom_tile(data = prediction.conc.n, aes(x,y,fill=get(colnames.aitch.pred[i]))) + ggtitle(datatitle) +
        guides(fill=guide_legend(title="")) +
        theme(plot.margin = unit(c(0.8, 0.5, 0.5, 1), "cm"),
              plot.title = element_text(size = 20, hjust=0.5),
              legend.text = element_text(size = 12),   # Dimensione delle etichette della legenda
              axis.title.x = element_blank(),  # Nascondere l'etichetta dell'asse X
              axis.title.y = element_blank(),   # Nascondere l'etichetta dell'asse Y
              #legend.title = element_text(size = 14),   # Dimensione del titolo della legenda
              #axis.title.x = element_text(size = 15),  # Dimensione del titolo dell'asse x
              #axis.title.y = element_text(size = 15),  # Dimensione del titolo dell'asse y
              axis.text.x = element_text(size = 12),   # Dimensione dei valori sull'asse x
              axis.text.y = element_text(size = 12)    # Dimensione dei valori sull'asse y
        ) +
        scale_fill_viridis_c(option = "C", direction = 1,
                             breaks = seq(min(datafill), max(datafill),
                            length.out = 10), labels = label_number(0.0001)) # legend name
      ggpred
      ggsave(paste(folder_output, paste0("prediction_cofrk_bausnopoints_component_", as.character(i),
          '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
          '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png"), sep="/"))
      obs_name = paste0("R", as.character(i), "_perc")
      ggpred = ggpred + geom_point(data = data.frame(data.comp), aes(x,y, fill=get(obs_name)), colour="black", pch=21, size=4)
      ggpred + geom_text(aes(x = data.comp$x, y = data.comp$y, label = data.comp$stationid), vjust = -1, size=2)
      ggsave(paste(folder_output, paste0("prediction_cofrk_baus_component_", as.character(i),
          '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
          '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png"), sep="/"))
    }
    graphics.off()
  }

  if (TRUE) # PLOT su base COMUNALE
  {
    comuniresult = comuni
    for (i in (1:length(colnames.aitch.pred)))
    {
      datafill = as.numeric(as.data.frame(comp.comuni.conc)[,colnames.aitch.pred[i]])
      comuniresult[colnames.aitch.pred[i]] <- datafill
      comuniresult$datafill = datafill
      if (DoPLOT) {
        x11(width = 8, height = 7)
        breaks = seq(min(datafill, na.rm = T), max(datafill, na.rm = T), length.out = 10)
        limits = c(min(datafill, na.rm = T), max(datafill, na.rm = T))
        if (Use_Global_Colors)
        {
          colors.pred.breaks = Colors_Breaks
          breaks = seq(colors.pred.breaks[[i]][1], colors.pred.breaks[[i]][2], length.out = 10)
          limits = c(colors.pred.breaks[[i]][1], colors.pred.breaks[[i]][2])
        }
        datatitle = paste0("CFRK prediction of relative frequencies \n in interval I", as.character(i),
                          " = ",  data_air$levels[i], " μg/m3 on municipal base")
        gcomuni = ggplot() + geom_sf(data = comp.comuni.conc$geometry, aes(fill=datafill), linewidth = 0.1) + ggtitle(datatitle) +
          guides(fill=guide_legend(title="Frequency")) +
          theme(plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm"),
                plot.title = element_text(size = 20, hjust=0.5),
                legend.text = element_text(size = 14),   # Dimensione delle etichette della legenda
                legend.title = element_text(size = 14),   # Dimensione del titolo della legenda
                axis.title.x = element_text(size = 15),  # Dimensione del titolo dell'asse x
                axis.title.y = element_text(size = 15),  # Dimensione del titolo dell'asse y
                axis.text.x = element_text(size = 12),   # Dimensione dei valori sull'asse x
                axis.text.y = element_text(size = 12)    # Dimensione dei valori sull'asse y
          ) + scale_fill_viridis_c(option = "C", direction = 1,
                                  breaks = breaks,
                                  limits = limits,
                                  labels = label_number(0.0001)) +# legend name
          labs(fill = "", x = "Longitude", y = "Latitude")

        gcomuni

        ggsave(paste(folder_output, paste0("prediction_cofrk_comunale_component_", as.character(i), ".png"), sep="/"))
      }
      # st_write(comuniresult, paste0(folder_output,'/prediction_cofrk_comunale_component_', i,
      #                               '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
      #                               '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, '.gpkg'), append = FALSE)
    }
    graphics.off()
  }
  comuniresult
}


CFRKSave <- function(comuniresult, exporting_colnames, folder_output,
    Global_Use_Cov_Population, Global_Use_Cov_Elevation, Flag_Export_only_PRO_COM = TRUE, suffix="") {

  filepath <- paste0(folder_output,'/prediction_cofrk_comunale_component_all', suffix)
  comuni_export <- comuniresult

  if (Flag_Export_only_PRO_COM) {
    comuni_export <- as.data.frame(comuniresult)[,c("PRO_COM", exporting_colnames)]
    cat(colnames(comuni_export),"\n\n")
    cat("Export", paste0(filepath, '.csv'), "\n")
    write.csv(comuni_export, paste0(filepath, '.csv'))
  } else {
    cat(colnames(comuni_export),"\n\n")
    cat("Export", paste0(filepath, '.gpkg'), "\n")
    st_write(comuniresult, paste0(filepath, '.gpkg'), append = FALSE)
  }
}

