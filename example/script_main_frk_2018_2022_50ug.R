rm(list = ls())
## -----------------------------------------------------------------------------
## Configure FLAGS and MACRO for program
flag_no_plot = TRUE
flag_plot_covariata = FALSE

# current_dir can be passed via command-line argument or environment var
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && nzchar(args[1])) {
  current_dir <- args[1]
} else {
  current_dir <- Sys.getenv("CURRENT_DIR", unset = "")
}

if (!nzchar(current_dir)) {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    stop("current_dir non impostato. Passalo come argomento Rscript o come variabile d'ambiente CURRENT_DIR.")
  }
}

source(paste0(current_dir, "/configure.conf"))
setwd(frk_working_dir)

cat("output_main_dir =", output_main_dir, "\n")

library(sf)
library(sp)
library(FRK)
library(ggplot2)

source(paste(cfrk_package_dir, 'CFRK_utils.R',  sep='/'))
source(paste(cfrk_package_dir, 'CFRK.R', sep='/'))

# ------------------------------------------------------------------------------
# config output path

dir_figures = paste0(output_main_dir, '/figures-quantiles-2018-2022-50ug-frk')
if(!dir.exists(dir_figures)) {dir.create(dir_figures, recursive = T)}


## Load data "data_air"
data_air <- readRDS(paste0(data_dir, 'data_air_2026.rds'))

#saveRDS(data_air, file=paste0(output_main_dir, "data_air.rds"))

colnames.aitch = paste0(data_air$colnames.aitch.abs, "_perc")
colnames.euclid = sapply(1:(length(data_air$levels) - 1),
  FUN=function(i){paste0("Y", as.character(i))})

# ------------------------------------------------------------------------------
# Initialize PCA
flag_use_pca = FALSE
# Initialize FRK
flag_bau_comunali = FALSE
flag_obs_fs = TRUE

# BAUs
baus.cellsize = 1600 # Dimensione delle BAU (sono le unità di massima risoluzione del modello)
# Basis
frk.args = list()
# opts_FRK$set("verbose", 1L)
frk.args["basis_type"] = "Gaussian"
# frk.args["basis_type"] = "Matern32"
frk.args["basis_nres"] = 3
Use_Global_Colors = FALSE

# ==============================================================================
### Covariate
comuni <- readRDS(paste0(data_dir, 'comuni.rds'))
library(terra)
data_elevazione <- rast(paste0(data_dir, 'data_elevazione.tif'))
data_popolazione_smoothed_nodes <- readRDS(paste0(data_dir, 'data_popolazione_smoothed_nodes.rds'))


# ==============================================================================
### Config

Global_Use_Cov_Population = TRUE
Global_Use_Cov_Elevation = TRUE
Use_Global_Colors = FALSE

# Flag_K_type_exponential = T
no_iterations = 200

# ------------------------------------------------------------------------

# Baus
GridBAUs0 <- GetGridBAUs(data_air$grid, cellsize=baus.cellsize)
GridBAUs1 <- Subset_Intersecating_Grid(GridBAUs0, comuni)
GridBAUs.Points <- SpatialPoints(coordinates(GridBAUs1), proj4string = CRS(proj4string(GridBAUs1)))


# popolazione
cov.P1 <- Interp_OnBaus_DensityPopOnNodes(data_air, GridBAUs1, data_popolazione_smoothed_nodes)
#cov.P1 <- Read_OnBaus_DensityPopOnNodes()
source(paste0(current_dir, "/utils_cov.R"))
cov.P1.trans <- impute_cov_P1(cov.P1, GridBAUs.Points, data_popolazione_smoothed_nodes)
#cov.P1.trans <- cov.P1^0.5
#spdf.cov.P1 <- SpatialPixelsDataFrame(points = GridBAUs.Points, data = data.frame(value = cov.P1.trans))



# elevazione
library(raster)
raster.elevazione <- raster(data_elevazione)
cov.Elev <- Interp_OnBaus_Elevation(data_air, GridBAUs1, raster.elevazione)
# cov.Elev.trans <- cov.Elev
# spdf.cov.Elev <- SpatialPixelsDataFrame(points = GridBAUs.Points, data = data.frame(value = cov.Elev.trans))
cov.Elev.trans <- impute_cov_Elev(cov.Elev, GridBAUs.Points, GridBAUs1)


DoBaus_WithCov <- function(GridBAUs1, cov.P1 = NULL, cov.Elev = NULL) {
  # cov.P1 <- Interp_OnBaus_DensityPopOnNodes(data_air, GridBAUs1, data_popolazione_smoothed_nodes)
  # if (is.null(cov.P1)) {
  #   cov.P1 <- array(Read_OnBaus_DensityPopOnNodes())[[1]]
  # }
  # if (is.null(cov.Elev)) {
  #   cov.Elev <- Interp_OnBaus_Elevation(data_air, GridBAUs1, raster.elevazione)
  # }

  cat("cov.P1", length(cov.P1), "\n")
  cat("cov.Elev", length(cov.Elev), "\n")
  discard_tiles = is.na(cov.P1) | is.na(cov.Elev)
  cat("dim(GridBAUs1)", dim(GridBAUs1[!discard_tiles,]), "\n")

  GridBAUs = GridBAUs1[!discard_tiles,]
  GridBAUs$Elevazione = cov.Elev[!discard_tiles]
  GridBAUs$Popolazione.P1 = cov.P1[!discard_tiles]
  GridBAUs
}
GridBAUs <- DoBaus_WithCov(GridBAUs1, cov.P1 = cov.P1.trans, cov.Elev = cov.Elev.trans) # calcolo BAUs con covariate trasformate

library(sp)

# Plot dell'intero SpatialPixelsDataFrame
# spplot(GridBAUs, "Elevazione")  # sostituisci "variable_name" con la colonna da visualizzare
# spplot(GridBAUs, "Popolazione.P1")  # sostituisci "variable_name" con la colonna da visualizzare
# saveRDS(GridBAUs, file=paste0(output_main_dir, "/GridBAUs.rds"))

# ----------------------------------------------------------------------------------
# frk.args.temp <- frk.args
# frk.args.temp$basis_nres <- 3
# domain.prediction.grid <- GetGrid(data_air$domain) # TODO construct GetGridBAUs without domain.prediction.gridv
# G.temp <- GetBasis(domain.prediction.grid, frk.args.temp)
#cat("done basis G\n")
#cat("Numero Funzioni Base = ", attributes()$n, "\n\n")
nrow(coordinates(GridBAUs))
# >> AUTO BASIS: using Nres =  4
# >> Numero Funzioni Base =  2889
# >> AUTO BASIS: using Nres =  3
# >> Numero Funzioni Base =  502
# cellcentre.offset cellsize cells.dim
# x          209247.4     1600       487
# y         4778943.9     1600       319

# ----------------------------------------------------------------------------------
pullutant_name = "PM10"
result <- CFRK (data.comp = data_air$dataset, flag_use_pca, flag_bau_comunali,
       flag_obs_fs, baus.cellsize, frk.args,
       domain = data_air$domain, comuni=comuni, cov.P1 = cov.P1, GridBAUs = GridBAUs,
       DestName.Pollutant = pullutant_name,
       comp.levels = c(), colnames.aitch = colnames.aitch, colnames.euclid = colnames.euclid,
       figures_folder = dir_figures,
       Global_Use_Cov_Population = Global_Use_Cov_Population, Global_Use_Cov_Elevation = Global_Use_Cov_Elevation,
       Use_Global_Colors = Use_Global_Colors, Colors_Breaks = Colors_Breaks,
       DoPLOT = TRUE, flag_no_plot = flag_no_plot)


for (name in paste0(colnames.aitch, ".pred")) {
  if(sum(is.na(result$comp.comuni.values[[name]])) > 0) {
    stop("Field ", name, " has ", sum(is.na(result$comp.comuni.values[[name]])),
    " municipalities with na \n",
    "❌ Errore: condizione non rispettata, esecuzione interrotta.")
  }
}


folder_output = GetFRKFolderOutput (dir_figures,
   "PM10", Global_Use_Cov_Population, Global_Use_Cov_Elevation,
   frk.args, flag_use_pca, flag_obs_fs, baus.cellsize, 0)

cat("K matrix type ", result$Slist[[1]]$S@K_type, "\n")

colnames.aitch.pred <- paste0(colnames.aitch, '.pred')

# plot
comuniresult <- CFRKPlotPrediction (
  data.comp = data_air$dataset, domain = data_air$domain,
  result$comp.prediction.values, result$comp.fitted.values, result$comp.comuni.values,
  colnames.aitch, colnames.aitch.pred, DoPLOT = TRUE)

res.confintervals <- CFRK_GetErrorInterval (result, data_air$dataset, comuni, index = 1)

confintervals.prediction <- res.confintervals$prediction.baus


upper <- GetAvg(confintervals.prediction, comuni, "Y.pred.upper", data_air$dataset)
ilrInv(upper[1])
as.data.frame(res.confintervals$prediction.comuni)[which(comuni$COMUNE == "Imperia"), "Y.pred.upper"]

comuniresult$sd <- CFRK_GetSD(result, data_air$dataset,comuni, index = 1)


comuniresult.2 <- cbind(comuniresult,
  as.data.frame(res.confintervals$prediction.comuni)[,
  c("Y.pred.upper","Y.pred.lower","Y.pred.upper.95","Y.pred.lower.95","Y.pred.upper.99","Y.pred.lower.99")])

# saveRDS(result, file=paste0(OUTPUT_MAIN_DIR, "/result_cfrk.rds"))

# save
CFRKSave (comuniresult.2, exporting_colnames =
  c(colnames.aitch.pred,"Y.pred.upper","Y.pred.lower","Y.pred.upper.95","Y.pred.lower.95","Y.pred.upper.99","Y.pred.lower.99"), folder_output,
  Global_Use_Cov_Population, Global_Use_Cov_Elevation, Flag_Export_only_PRO_COM = TRUE)



