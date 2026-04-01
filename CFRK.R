if (exists("cfrk_package_dir")) {
  source(paste(cfrk_package_dir, 'CFRK_utils.R',  sep='/'))
} else {
  source(paste(CFRK_PACKAGE_DIR, 'CFRK_utils.R',  sep='/'))
}

# CFRK method implementation
no_iterations = 20

# Main CFRK function for spatial prediction using Fixed Rank Kriging
# Parameters:
# - data.comp: Compositional data for kriging
# - flag_use_pca: Boolean to use PCA transformation
# - flag_bau_comunali: Flag for BAU at municipal level
# - flag_obs_fs: Flag for observation fine-scale variation
# - baus.cellsize: Cell size for BAUs
# - frk.args: Arguments for FRK basis functions
# - domain: Spatial domain for prediction
# - comuni: Municipal boundaries
# - cov.P1: Population covariate (optional)
# - GridBAUs: BAU grid (optional)
# - DestName.IndexDrop: Index for dropping components
# - DestName.Pollutant: Pollutant name for output
# - fitted.locations: Locations for fitted values (optional)
# - comp.levels: Compositional levels
# - colnames.aitch: Column names for aitch components
# - colnames.euclid: Column names for Euclidean components
# - figures_folder: Folder for output figures
# - Global_Use_Cov_Population: Use population covariate
# - Global_Use_Cov_Elevation: Use elevation covariate
# - Use_Global_Colors: Use global color scales
# - Colors_Breaks: Color breaks for plots
# - DoPLOT: Boolean to generate plots
# - flag_no_plot: Flag to suppress plots
CFRK <- function(data.comp, flag_use_pca, flag_bau_comunali, flag_obs_fs, baus.cellsize, frk.args,
  domain, comuni, cov.P1 = NULL, GridBAUs = NULL,
  DestName.IndexDrop = 0, DestName.Pollutant = "", fitted.locations = NULL,
  comp.levels = c(), colnames.aitch = c(), colnames.euclid = c(), figures_folder = "./figures/",
  Global_Use_Cov_Population = TRUE, Global_Use_Cov_Elevation = FALSE,
  Use_Global_Colors = TRUE, Colors_Breaks = NULL, DoPLOT = TRUE, flag_no_plot = FALSE)
{
  # Prepare column names for predictions
  colnames_pred = paste0(colnames.euclid, '.pred')
  colnames.aitch.pred = paste0(colnames.aitch, '.pred')

  # Get output folder for FRK results
  folder_output = GetFRKFolderOutput (figures_folder,
   DestName.Pollutant, Global_Use_Cov_Population, Global_Use_Cov_Elevation,
   frk.args, flag_use_pca, flag_obs_fs, baus.cellsize, 0)

  cat("\nDoFixedRankCokriging: Using PCA = ", flag_use_pca, "\n")
  cat("Folder output = ", folder_output, "\n")

  # Load required libraries
  library(magrittr)
  library(FRK) # for carrying out FRK

  # 1. PCA transformation (ensure 'data_air' is defined beforehand)
  data_for_kriging <- PrepareData_For_Kriging(
    data.comp, colnames.euclid, flag_use_pca)

  data.krig <- data_for_kriging[["data.krig"]]
  pca_output <- data_for_kriging[["pca_output"]]

  # ============================================================================
  # Correlation check on kriging data
  cor(as.data.frame(data.krig)[colnames.euclid])
  # ============================================================================

  # Set fitted locations if not provided
  if (is.null(fitted.locations)) {
    fitted.locations <- SpatialPoints(
      setNames(as.data.frame(data.comp)[c("x", "y")], c("X", "Y")),
      proj4string = CRS(paste("EPSG:", st_crs(data.comp)$epsg)))
  }

  # Perform FRK prediction
  resultFRK <- DoFRKPrediction(domain, GridBAUs, frk.args, data.comp,
    data_for_kriging[['data.krig']], colnames.euclid, comuni,
    Global_Use_Cov_Population, Global_Use_Cov_Elevation,
    figures_folder=figures_folder, no_iterations,
    fitted.locations = fitted.locations)


  predictionbaus.values = resultFRK$predictionbaus.values
  fitted.values = resultFRK$fitted.values
  comuni.values = resultFRK$comuni.values

  # Invert PCA transformation if PCA was used
  if(flag_use_pca)
  {
    # Transform BAU predictions back to original basis
    originalbasis.data <- ScoresToVars(
      as.matrix(as.data.frame(predictionbaus.values)[,colnames_pred]),
      pca_output, colnames=colnames.euclid) %>%
      set_colnames(colnames_pred)

    # Transform fitted values back to original basis
    originalbasis.fitted.values <- ScoresToVars(
      as.matrix(as.data.frame(fitted.values)[,colnames_pred]),
      pca_output, colnames=colnames.euclid) %>%
      set_colnames(colnames_pred)

    # Transform municipal values back to original basis
    originalbasis.municipal.values <- ScoresToVars(
      as.matrix(as.data.frame(comuni.values)[,colnames_pred]),
      pca_output, colnames=colnames.euclid) %>%
      set_colnames(colnames_pred)

    predictionbaus.values <- SpatialPixelsDataFrame(
      coordinates(predictionbaus.values), data = originalbasis.data)
    fitted.values <- SpatialPointsDataFrame(
      coordinates(fitted.locations), data = originalbasis.fitted.values)
  } else {
    originalbasis.municipal.values = comuni.values
  }

  # Invert ILR transformation
  # 1/3 Invert ILR for BAU predictions
  cat("Inverting ILR for prediction BAUs...\n")
  data.to.revert.prediction =
    as.data.frame(predictionbaus.values)[,colnames_pred] %>% as.data.frame
  colnames(data.to.revert.prediction) <- colnames_pred
  prediction.conc = RevertILR(data.to.revert.prediction,
    predictionbaus.values$x, predictionbaus.values$y)
  # Select only pixels within the domain (e.g., Lombardy region)
  prediction.conc <- prediction.conc[
    unlist(st_contains(data_air$domain,
    st_as_sf(prediction.conc, coords=c("x","y"),
    crs=st_crs(data.comp)$epsg)$geometry)),]

  # 2/3 Invert ILR for fitted values
  cat("Inverting ILR for fitted locations...\n")
  data.to.revert.fitted =
    as.data.frame(fitted.values)[colnames_pred] %>% as.data.frame
  colnames(data.to.revert.fitted) <- colnames_pred
  comp.fitted.df = RevertILR(data.to.revert.fitted,
    as.data.frame(fitted.values)$X,  as.data.frame(fitted.values)$Y)
  comp.fitted.df = cbind(comp.fitted.df, as.data.frame(fitted.values)[colnames_pred])

  # 3/3 Invert ILR for municipal values
  cat("Inverting ILR for municipal values...\n")
  data.to.revert.comuni =
    as.data.frame(originalbasis.municipal.values)[,colnames_pred] %>% as.data.frame
  colnames(data.to.revert.comuni) <- colnames_pred
  comp.comuni.conc = comuni
  comp.comuni.conc[,colnames.aitch.pred] = ilrInv(data.to.revert.comuni)

  # Return results as a list
  list(
    comp.prediction.values = prediction.conc,
    comp.fitted.values = comp.fitted.df,
    comp.comuni.values = comp.comuni.conc,
    Slist = resultFRK$Slist)
}

# Function to extract standard deviation from CFRK results
# Parameters:
# - result: CFRK result object
# - data_air_dataset: Air quality dataset
# - comuni: Municipal boundaries
# - index: Component index (default 1)
CFRK_GetSD <- function(result, data_air_dataset, comuni, index = 1) {
  prediction = SpatialPixelsDataFrame(
    as.matrix(as.data.frame(result$Slist[[index]]$BAUs_prediction_df)[,c('x','y')]),
    data = data.frame(Y.pred.sd = result$Slist[[index]]$BAUs_prediction_df$sd),
    proj4string = CRS(paste("EPSG:", st_crs(data_air_dataset)$epsg)))

  GetAvg(prediction, comuni, "Y.pred.sd", data_air_dataset)
}

# Function to compute prediction error intervals from CFRK results
# Parameters:
# - result: CFRK result object
# - data_air_dataset: Air quality dataset
# - comuni: Municipal boundaries
# - index: Component index (default 1)
CFRK_GetErrorInterval <- function (result, data_air_dataset, comuni, index = 1) {
  prediction = SpatialPixelsDataFrame(
    as.matrix(as.data.frame(result$Slist[[index]]$BAUs_prediction_df)[,c('x','y')]),
    data = data.frame(
      Y.pred.upper = result$Slist[[index]]$BAUs_prediction_df$mu + result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.lower = result$Slist[[index]]$BAUs_prediction_df$mu - result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.upper.95 = result$Slist[[index]]$BAUs_prediction_df$mu + qnorm(0.975) * result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.lower.95 = result$Slist[[index]]$BAUs_prediction_df$mu - qnorm(0.975) * result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.upper.99 = result$Slist[[index]]$BAUs_prediction_df$mu + qnorm(0.995) * result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.lower.99 = result$Slist[[index]]$BAUs_prediction_df$mu - qnorm(0.995) * result$Slist[[index]]$BAUs_prediction_df$sd
      ),
    proj4string = CRS(paste("EPSG:", st_crs(data_air_dataset)$epsg)))

  comp.comuni.conc = comuni
  # Compute upper bound
  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.upper", data_air_dataset)
  comp.comuni.conc[,"Y.pred.upper"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  # Compute lower bound
  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.lower", data_air_dataset)
  comp.comuni.conc[,"Y.pred.lower"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  # Compute 95% upper bound
  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.upper.95", data_air_dataset)
  comp.comuni.conc[,"Y.pred.upper.95"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  # Compute 95% lower bound
  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.lower.95", data_air_dataset)
  comp.comuni.conc[,"Y.pred.lower.95"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  # Compute 99% upper bound
  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.upper.99", data_air_dataset)
  comp.comuni.conc[,"Y.pred.upper.99"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  # Compute 99% lower bound
  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.lower.99", data_air_dataset)
  comp.comuni.conc[,"Y.pred.lower.99"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  list(prediction.baus = prediction, prediction.comuni = comp.comuni.conc)
}

