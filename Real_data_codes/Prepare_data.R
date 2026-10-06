## Builds Real_data_codes/EPA_data/monitors.csv from the public source files. The
## file is provided, so this script is needed only to rebuild it.
##
## Sources (download into Real_data_codes/EPA_data/):
##   EPA Air Quality System, https://aqs.epa.gov/aqsweb/airdata/download_files.html:
##     annual_conc_by_monitor_2022.csv and aqs_sites.csv
##   WorldPop population density 2020, 1 km, unconstrained, UN-adjusted, United States,
##     https://hub.worldpop.org/: usa_pd_2020_1km.tif
##   TerraClimate 10 m wind speed (ws), monthly 2022,
##     https://www.climatologylab.org/terraclimate.html: TerraClimate_ws_2022.nc
library(readr); library(dplyr); library(sf); library(terra)
src <- function(f) file.path("Real_data_codes/EPA_data", f)

## PM2.5 (parameter 88101), daily samples, at least 250 observations, exceptional events
## kept, Alaska, Hawaii and Puerto Rico dropped, one monitor per site (lowest POC)
m <- read_csv(src("annual_conc_by_monitor_2022.csv"), show_col_types = FALSE) |>
  filter(`Parameter Code` == 88101,
         `Sample Duration` %in% c("24 HOUR", "24-HR BLK AVG"),
         `Observation Count` >= 250,
         `Event Type` %in% c("Events Included", "No Events"),
         !`State Code` %in% c("02", "15", "72")) |>
  arrange(POC) |>
  distinct(`State Code`, `County Code`, `Site Num`, .keep_all = TRUE) |>
  transmute(station_id = paste(`State Code`, `County Code`, `Site Num`, sep = "_"),
            State.Code = as.integer(`State Code`), County.Code = as.integer(`County Code`),
            Site.Number = as.integer(`Site Num`),
            lat = Latitude, lon = Longitude, pm25 = `Arithmetic Mean`)
sites <- read.csv(src("aqs_sites.csv")) |>
  transmute(State.Code = as.integer(State.Code), County.Code = as.integer(County.Code),
            Site.Number = as.integer(Site.Number), Elevation, Land.Use, Location.Setting)
m <- left_join(m, sites, by = c("State.Code", "County.Code", "Site.Number"))

## Albers equal-area coordinates in km, shifted to the origin and divided by the largest
## coordinate (4547 km), so that the monitors lie in the unit square
xy <- st_as_sf(m, coords = c("lon", "lat"), crs = 4326) |> st_transform(crs = 5070) |>
  st_coordinates()
shifted <- sweep(xy, 2, apply(xy, 2, min), "-")
m$V11 <- shifted[, 1] / max(shifted); m$V12 <- shifted[, 2] / max(shifted)

## land use and location setting as indicators; the other levels, and missing values,
## are the reference
m <- m |> mutate(
  commercial  = as.integer(Land.Use %in% "COMMERCIAL"),
  residential = as.integer(Land.Use %in% "RESIDENTIAL"),
  industrial  = as.integer(Land.Use %in% "INDUSTRIAL"),
  desert      = as.integer(Land.Use %in% "DESERT"),
  forest      = as.integer(Land.Use %in% "FOREST"),
  mobile      = as.integer(Land.Use %in% "MOBILE"),
  suburban    = as.integer(Location.Setting %in% "SUBURBAN"),
  rural       = as.integer(Location.Setting %in% "RURAL"))

## population density (log1p, since some monitors are in zero-population cells) and the
## 2022 annual mean wind speed, extracted at the monitor coordinates
pts <- vect(st_as_sf(m, coords = c("lon", "lat"), crs = 4326))
m$pop_density <- terra::extract(rast(src("usa_pd_2020_1km.tif")), pts)[, 2]
m$log_pop <- log1p(m$pop_density)
m$wind_ms <- terra::extract(mean(rast(src("TerraClimate_ws_2022.nc"))), pts)[, 2]

## the continuous covariates are standardized
for (v in c("Elevation", "log_pop", "wind_ms")) m[[v]] <- as.numeric(scale(m[[v]]))
write.csv(m, src("monitors.csv"), row.names = FALSE)
