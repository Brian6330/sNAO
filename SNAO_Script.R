####### Read EKF400 #######
library(ncdf4) # Read NC Files
library(maps) # Plot Land Masks

# Set working directory
setwd("C:/Drive/_Klima 3/Exercises/R Code")
setwd("/Users/louiszwyssig/Desktop/Uni/7. Semester/Climatology III/R_work_directory")

# Open dataset
ekf400 <- nc_open("EKF400v2.0_ensmean.nc")
lon <- ekf400$dim[[1]]$vals    # Dimension longitude of ekf400
lat <- ekf400$dim[[2]]$vals    # Dimension latitude	of	ekf400

####### Define Reference Period #######
start_year <- 1950 # first possible year is 1603
end_year <- 2000 # last possible year is 2002
start <- start_year - 1602 # year 1603 is start = 1
end <- end_year - 1602 # year 2002 is end = 400

####### Define Climate Variables #######
# Calculate annual mean for sea-level pressure (slp)
slp_annualmean <- array(NA, dim = c(192, 96, 400))
for (i in 1603:2002) {
	slp_annualmean[, , i - 1602] <-
		apply(ncvar_get(
			ekf400, varid = "slp",
			start = c(1, 1, (i - 1602) * 12 + 1), # i_1: (1603-1602)*12+1 = 13
			count = c(-1, -1, 12)
		), c(1, 2), mean)
}

# Calculate high Summer (Jul. / Aug.) mean for slp
slp_highSum_mean <- array(NA, c(192,96,400))
for (i in 1603:2002) {
  slp_highSum_mean[, , i - 1602] <-
		apply(ncvar_get(
			ekf400, varid = "slp",
			start = c(1, 1, (i - 1602) * 12 + 7), # (1603-1602)*12+1 = 13
			count = c(-1, -1, 2)
		), c(1, 2), mean)
}

# Calculate high Summer (Jul. / Aug.) mean for temperature
temp_highSum_mean <- array(NA, dim = c(192, 96, 400))
for (i in 1603:2002) {
	temp_highSum_mean[, , i - 1602] <-
		apply(ncvar_get(
			ekf400, varid = "t2m",
			start = c(1, 1, (i - 1602) * 12 + 7), # 7 = July
			count = c(-1, -1, 2) # 2 = July & August
		), c(1, 2), mean)
}

# Calculate high Summer (Jul. / Aug.) mean for precipitation
prec_highSum_mean <- array(NA, dim = c(192, 96, 400))
for (i in 1603:2002) {
	prec_highSum_mean[, , i - 1602] <-
		apply(ncvar_get(
			ekf400, varid = "prec",
			start = c(1, 1, (i - 1602) * 12 + 7), # 7 = July
			count = c(-1, -1, 2) # 2 = July & August
		), c(1, 2), mean)
}

# Calculate high Summer (Jul. / Aug.) mean for geopotential height at 500hPa
h500_highSum_mean <- array(NA, dim = c(192, 96, 400))
for (i in 1603:2002) {
	h500_highSum_mean[, , i - 1602] <-
		apply(ncvar_get(
			ekf400, varid = "hgt500",
			start = c(1, 1, (i - 1602) * 12 + 7), # 7 = July
			count = c(-1, -1, 2) # 2 = July & August
		), c(1, 2), mean)
}

# Extract global timeseries (1603-2003) for slp from ekf400
slp_globaltimeseries <- ncvar_get(
 	ekf400, varid = "slp", start = c(1, 1, 13), count = c(-1, -1, -1)
)

####### Principal Component #######
# Coordinates over region 25.8-70.8N, 70.8W-50.8E.
sNAO_lon_min <- 289.2 #E = 360-70.8W
sNAO_lon_max <- 50.8  #E (0-50.8)
sNAO_lat_min <- 25.8  #N
sNAO_lat_max <- 70.8  #N

# Select gridcell of sNAO:
pc_lon_index <- ((0 < lon) & (lon < sNAO_lon_max)) | ((sNAO_lon_min < lon) & (lon < 360))
pc_lat_index <- (sNAO_lat_min < lat) & (sNAO_lat_max > lat)
lon_adjusted <- lon[pc_lon_index]
lat_adjusted <- lat[pc_lat_index]

# Calculate anomaly for the SNAO relevant lon and lat
slp_highSummer_filtered <- slp_highSum_mean[pc_lon_index, pc_lat_index,start:end]
slp_highSummer_filtered_anom <-
  slp_highSummer_filtered - mean(slp_highSummer_filtered)

# Principal component analysis: standardize, weighting (sqrt(cos(lat)))
y_scaled <- apply(slp_highSummer_filtered_anom, c(1, 2), scale)
slp_highSummer_sd <- apply(slp_highSummer_filtered_anom, c(1, 2), sd)
y_weighted <- apply(y_scaled, c(1, 2), "*", sqrt(cos(pi * lat_adjusted / 180)))
y_perm <- aperm(y_weighted, c(2, 3, 1))
y_mat <- matrix(y_perm, nrow = length(start:end), ncol = length(lon_adjusted) * length(lat_adjusted))
pc <- prcomp(y_mat, retx = T, scale. = F)

# Plot time series
plottime <- pc$x[, 1] # 1 = first pc
plot(start_year:end_year, plottime, type = "l")

# PC back-weighted and -standardized
plotpc <- t(apply(
		matrix(pc$rotation[, 1], nrow = dim(lon_adjusted), ncol = dim(lat_adjusted)),
		1, "/", sqrt(cos(pi * lat[pc_lat_index] / 180))
	)) * slp_highSummer_sd

x_pc <- c(lon_adjusted[lon_adjusted > 180] - 360,
					lon_adjusted[lon_adjusted <= 180])
y_pc <- rev(lat_adjusted)
z_pc <- rbind(plotpc[lon_adjusted > 180, rev(1:length(lat_adjusted))],
						plotpc[lon_adjusted <= 180, rev(1:length(lat_adjusted))])

mycol <- c("blue2", "cornflowerblue", "cadetblue2", "azure2", "bisque1", "burlywood1", "brown1", "brown3")
minmax <- max(abs(c(min(z_pc), max(z_pc))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_pc, y_pc, z_pc,
							 levels = mylevs, col = mycol,
							 plot.axes = {
							 	map("world", interior = F, add = T)
							 }, main = "Principal Component Analysis SLP & sNAO 1603-2002")

# Extract coordinates of centers of action for Regression
max_index <- which(plotpc == max(plotpc), arr.ind = TRUE)
min_index <- which(plotpc == min(plotpc), arr.ind = TRUE)
max_pos <- c(lon_adjusted[max_index[1,1]], lat_adjusted[max_index[1,2]])
min_pos <- c(lon_adjusted[min_index[1,1]], lat_adjusted[min_index[1,2]])

####### Extract time-series for SLP (Minimum Center of Action) #######
# Use coordinates from principal component analysis
c_minCoA_lon <- min_pos[1]
c_minCoA_lat <- min_pos[2]

dlon <- lon[2] - lon[1] # Longitudinal grid space
dlat <- lat[1] - lat[2] # Latitudinal	 grid space

# Select gridcell:
minCoA_lon_index <- (lon < (c_minCoA_lon + (dlon / 2))) & (lon > (c_minCoA_lon - (dlon / 2)))
minCoA_lat_index <- (lat < (c_minCoA_lat + (dlat / 2))) & (lat > (c_minCoA_lat - (dlat / 2)))

# Select time series
slp_series_minCoA <- slp_globaltimeseries[minCoA_lon_index, minCoA_lat_index,]

# Calculate anomalies for every month for slp
slp_anom_minCoA <- array(NA, dim = 4800)
for (i in 1:4800) {
	slp_anom_minCoA[i] <-
		slp_series_minCoA[i] - mean(slp_annualmean[minCoA_lon_index, minCoA_lat_index, start:end])
}

# Calculate standardized anomalies for every month for slp
stdev_minCoA <- array(sd(slp_series_minCoA), dim = 4800)
slp_anom_stand_minCoA <- slp_anom_minCoA / stdev_minCoA

####### Extract time-series for SLP (Maximum Center of Action) #######
# Use coordinates from principal component analysis
c_maxCoA_lon <- max_pos[1]
c_maxCoA_lat <- max_pos[2]

# Select gridcell:
maxCoA_lon_index <- (lon < (c_maxCoA_lon + (dlon / 2))) & (lon > (c_maxCoA_lon - (dlon / 2)))
maxCoA_lat_index <- (lat < (c_maxCoA_lat + (dlat / 2))) & (lat > (c_maxCoA_lat - (dlat / 2)))

# Select time series
slp_series_maxCoA <- slp_globaltimeseries[maxCoA_lon_index, maxCoA_lat_index,]

# Calculate anomalies for every month for slp
slp_anom_maxCoA <- array(NA, dim = 4800)
for (i in 1:4800) {
	slp_anom_maxCoA[i] <-
		slp_series_maxCoA[i] - mean(slp_annualmean[maxCoA_lon_index, maxCoA_lat_index, start:end])
}

# Calculate standardized anomalies for every month for slp
stdev_maxCoA <- array(sd(slp_series_maxCoA), dim = 4800)
slp_anom_stand_maxCoA <- slp_anom_maxCoA / stdev_maxCoA

####### NAO & sNAO Index #######
# NAO index with SLPasc: MaxCenterOfAction - MinCenterOfAction
NAO_index_asc <- slp_anom_stand_maxCoA - slp_anom_stand_minCoA

# Plot the last years
time_offset <- 12 # 1 year = 12 months
years <- 20
plot(NAO_index_asc[(length(NAO_index_asc) - years * time_offset):length(NAO_index_asc)],
		 type = "l", main = "NAO Index", xlab = "Months", ylab = ""
)

# Aggregate by Boreal Winter and Boreal Summer
NAO_index_asc_win <- array(NA, dim = 400)
start_month <- 12 # December
years <- 400 # 1603-2002
for (i in seq(from = start_month, to = years * time_offset, by = time_offset)) {
	winter <- array(c(NAO_index_asc[i], NAO_index_asc[i + 1], NAO_index_asc[i + 2]), dim = 3)
	NAO_index_asc_win[ceiling(i/12)] <- mean(winter, na.rm = TRUE) # na.rm required due to last January = NA
}

NAO_index_asc_sum <- array(NA, dim = 400)
start_month <- 6 # June
years <- 400 # 1603-2002
for (i in seq(from = start_month, to = years * time_offset, by = time_offset)) {
	summer <- array(c(NAO_index_asc[i], NAO_index_asc[i + 1], NAO_index_asc[i + 2]), dim = 3)
	NAO_index_asc_sum[ceiling(i/12)] <- mean(summer)
}

sNAO_index_asc <- array(NA, dim = 400)
start_month <- 7 # July of the second year
years <- 400 # 1603-2002
for (i in seq(from = start_month, to = years * time_offset, by = time_offset)) {
	summer <- array(c(NAO_index_asc[i], NAO_index_asc[i + 1]), dim = 2)
	sNAO_index_asc[ceiling(i/12)] <- mean(summer)
}

# Plot the last years for winter NAO, summer NAO and sNAO
years <- 20
plot(NAO_index_asc_win[(length(NAO_index_asc_win) - years):length(NAO_index_asc_win)],
		 type = "l", main = "NAO Index Winter",
		 xlab = "Year", ylab = ""
)
plot(NAO_index_asc_sum[(length(NAO_index_asc_sum) - years):length(NAO_index_asc_sum)],
		 type = "l", main = "NAO Index Summer",
		 xlab = "Year", ylab = ""
)
plot(sNAO_index_asc[(length(sNAO_index_asc) - years):length(sNAO_index_asc)],
		 type = "l", main = "sNAO Index",
		 xlab = "Year", ylab = ""
)

####### Regression #######
# Read-in predictors
x_SNAO <- as.matrix(cbind(1, sNAO_index_asc[start:end]))

# Regression Calculation
xm_SNAO <- solve(t(x_SNAO) %*% x_SNAO) %*% t(x_SNAO)
b_temp <- apply(temp_highSum_mean[,,start:end], c(1, 2), "%*%", t(xm_SNAO))
b_prec <- apply(prec_highSum_mean[,,start:end], c(1, 2), "%*%", t(xm_SNAO))
b_h500 <- apply(h500_highSum_mean[,,start:end], c(1, 2), "%*%", t(xm_SNAO))
b_slp <- apply(slp_highSum_mean[,,start:end], c(1, 2), "%*%", t(xm_SNAO))

# General Data for Plots
lon_limit <- c(52,282)
lat_limit <- c(22,82)
lon_range <- ((0 < lon) & (lon < lon_limit[1])) | ((lon_limit[2] < lon) & (lon < 360))
lat_range <- (lat_limit[1] < lat) & (lat_limit[2] > lat)
lon_shortened <- lon[lon_range]
lat_shortened <- lat[lat_range]
x_lon <- c(lon_shortened[lon_shortened > 180] - 360,
					 lon_shortened[lon_shortened <= 180])
y_lat <- rev(lat_shortened)
mycol <- c("blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3")

# Plot for Temperature and SNAO
field_temp <- b_temp[2, lon_range, lat_range]
minmax <- max(abs(c(min(field_temp), max(field_temp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
z_temp <- rbind(field_temp[lon_shortened > 180, rev(1:length(lat_shortened))],
								field_temp[lon_shortened <= 180, rev(1:length(lat_shortened))])
filled.contour(x_lon, y_lat, z_temp,
							 levels = mylevs, col = mycol,
							 plot.axes = {map("world", interior = F, add = T)}
							 , main = "sNAO Regression Temperature from 1603-2002",)

# Plot for Precipitation and SNAO
field_prec <- b_prec[2, lon_range, lat_range]
minmax <- max(abs(c(min(field_prec), max(field_prec))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
z_prec <- rbind(field_prec[lon_shortened > 180, rev(1:length(lat_shortened))],
								field_prec[lon_shortened <= 180, rev(1:length(lat_shortened))])
filled.contour(x_lon, y_lat, z_prec,
							 levels = mylevs, col = rev(mycol),
							 plot.axes = {map("world", interior = F, add = T)}
							 , main = "sNAO Regression Precipitation from 1603-2002",)

# Plot for Geopotential Height and SNAO
field_h500 <- b_h500[2, lon_range, lat_range]
minmax <- max(abs(c(min(field_h500), max(field_h500))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
z_h500 <- rbind(field_h500[lon_shortened > 180, rev(1:length(lat_shortened))],
								field_h500[lon_shortened <= 180, rev(1:length(lat_shortened))])
filled.contour(x_lon, y_lat, z_h500,
							 levels = mylevs, col = mycol,
							 plot.axes = {map("world", interior = F, add = T)}
							 , main = "sNAO Regression Geopotential Height from 1603-2002",)

# Plot for Sea-level Pressure and SNAO
field_slp <- b_slp[2, lon_range, lat_range]
minmax <- max(abs(c(min(field_slp), max(field_slp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
z_slp <- rbind(field_slp[lon_shortened > 180, rev(1:length(lat_shortened))],
               field_slp[lon_shortened <= 180, rev(1:length(lat_shortened))])
filled.contour(x_lon, y_lat, z_slp,
               levels = mylevs, col = mycol,
               plot.axes = {map("world", interior = F, add = T)}
               , main = "sNAO Regression Sea-level Pressure from 1603-2002",)

####### Correlation #######
# Calculate Correlation Fields
SNAO_matrix <- as.matrix(cbind(sNAO_index_asc[start:end]))
cor_field_sNAO_temp <- apply(temp_highSum_mean[,,start:end], c(1, 2), cor, SNAO_matrix)
cor_field_sNAO_prec <- apply(prec_highSum_mean[,,start:end], c(1, 2), cor, SNAO_matrix)
cor_field_sNAO_h500 <- apply(h500_highSum_mean[,,start:end], c(1, 2), cor, SNAO_matrix)
cor_field_sNAO_slp <- apply(slp_highSummer_mean[,,start:end], c(1, 2), cor, SNAO_matrix)

# General Data for Plots
lon_limit <- c(52,282)
lat_limit <- c(22,82)
lon_range <- ((0 < lon) & (lon < lon_limit[1])) | ((lon_limit[2] < lon) & (lon < 360))
lat_range <- (lat_limit[1] < lat) & (lat_limit[2] > lat)
lon_shortened <- lon[lon_range]
lat_shortened <- lat[lat_range]
x_lon <- c(lon_shortened[lon_shortened > 180] - 360,
           lon_shortened[lon_shortened <= 180])
y_lat <- rev(lat_shortened)
mycol <- c("blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3")

# Plot for Temperature and SNAO
cor_field_temp <- cor_field_sNAO_temp[lon_range, lat_range]
minmax <- max(abs(c(min(cor_field_temp), max(cor_field_temp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
z_cor_temp <- rbind(cor_field_temp[lon_shortened > 180, rev(1:length(lat_shortened))],
                    cor_field_temp[lon_shortened <= 180, rev(1:length(lat_shortened))])
filled.contour(x_lon, y_lat, z_cor_temp,
               levels = mylevs, col = mycol,
               plot.axes = {map("world", interior = F, add = T)}
               , main = "Correlation sNAO and Temperature 1603-2002")

# Plot for Precipitation and SNAO
cor_field_prec <- cor_field_sNAO_prec[lon_range, lat_range]
minmax <- max(abs(c(min(cor_field_prec), max(cor_field_prec))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
z_cor_prec <- rbind(cor_field_prec[lon_shortened > 180, rev(1:length(lat_shortened))],
                    cor_field_prec[lon_shortened <= 180, rev(1:length(lat_shortened))])
filled.contour(x_lon, y_lat, z_cor_prec,
               levels = mylevs, col = rev(mycol),
               plot.axes = {map("world", interior = F, add = T)}
               , main = "Correlation sNAO and Precipitation 1603-2002")

# Plot for Geopotential Height and SNAO
cor_field_h500 <- cor_field_sNAO_h500[lon_range, lat_range]
minmax <- max(abs(c(min(cor_field_h500), max(cor_field_h500))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
z_cor_h500 <- rbind(cor_field_h500[lon_shortened > 180, rev(1:length(lat_shortened))],
                    cor_field_h500[lon_shortened <= 180, rev(1:length(lat_shortened))])
filled.contour(x_lon, y_lat, z_cor_h500,
               levels = mylevs, col = mycol,
               plot.axes = {map("world", interior = F, add = T)}
               , main = "Correlation sNAO and Geopotential Height 1603-2002")

# Plot for Sea-level Pressure and SNAO
cor_field_slp <- cor_field_sNAO_slp[lon_range, lat_range]
minmax <- max(abs(c(min(cor_field_slp), max(cor_field_slp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
z_cor_slp <- rbind(cor_field_slp[lon_shortened > 180, rev(1:length(lat_shortened))],
                   cor_field_slp[lon_shortened <= 180, rev(1:length(lat_shortened))])
filled.contour(x_lon, y_lat, z_cor_slp,
               levels = mylevs, col = mycol,
               plot.axes = {map("world", interior = F, add = T)}
               , main = "Correlation sNAO and Sea-level Pressure 1603-2002")

####### Case Study Timeseries #######
# Coordinates of London: 51.5N, 0E
c_case_lon <- 0
c_case_lat <- 51.5
dlon <- lon[2]-lon[1] #longitudinal grid space
dlat <- lat[1]-lat[2] #latitudinal grid space

# Index to select gridcell
case_lon_index <- (lon < (c_case_lon + (dlon / 2))) & (lon > (c_case_lon - (dlon / 2)))
case_lat_index <- (lat < (c_case_lat + (dlat / 2))) & (lat > (c_case_lat - (dlat / 2)))

# Temperature timeseries for case study
case_timeseries_temp <- temp_highSum_mean[case_lon_index,case_lat_index,] - mean(temp_highSum_mean[case_lon_index,case_lat_index,])
case_timeseries_temp_stand <- case_timeseries_temp/sd(case_timeseries_temp)

# Precipitation timeseries for case study
case_timeseries_prec <- prec_highSum_mean[case_lon_index,case_lat_index,] - mean(prec_highSum_mean[case_lon_index,case_lat_index,])
case_timeseries_prec_stand <- case_timeseries_prec/sd(case_timeseries_prec)

# 500hPa-height timeseries for case study
case_timeseries_h500 <- h500_highSum_mean[case_lon_index,case_lat_index,] - mean(h500_highSum_mean[case_lon_index,case_lat_index,])
case_timeseries_h500_stand <- case_timeseries_h500/sd(case_timeseries_h500)

# Sea-level Pressure timeseries for case study
case_timeseries_slp <- slp_highSum_mean[case_lon_index,case_lat_index,] - mean(slp_highSum_mean[case_lon_index,case_lat_index,])
case_timeseries_slp_stand <- case_timeseries_slp/sd(case_timeseries_slp)

# See section "Define Reference Period" (line 14) to set reference period

# Define vector for corresponding period to be plotted
sNAO_index_asc_plot <- sNAO_index_asc[start:end]
case_temp_plot <- case_timeseries_temp_stand[start:end]
case_prec_plot <- case_timeseries_prec_stand[start:end]
case_h500_plot <- case_timeseries_h500_stand[start:end]
case_slp_plot <- case_timeseries_slp_stand[start:end]

# Plot sNAO, temperature & precipitation timeseries
plot(sNAO_index_asc_plot,
     type = "l", col = "black",
     main = "sNAO Index, Precipitation & Temperature", xlab = , ylab = "",
     ylim = c(min(c(sNAO_index_asc_plot, case_temp_plot, case_prec_plot)),
              max(c(sNAO_index_asc_plot, case_temp_plot, case_prec_plot)))
)
lines(case_temp_plot, col = "darkred")
lines(case_prec_plot*-1, col = "darkblue") # Precipitation and sNAO index correlate negatively
      
# Correlation coefficients of the timeseries and sNAO index for the reference period
cor(sNAO_index_asc_plot, case_temp_plot)
cor(sNAO_index_asc_plot, case_prec_plot)
cor(sNAO_index_asc_plot, case_h500_plot)
cor(sNAO_index_asc_plot, case_slp_plot)

####### Case Study Pattern #######
# Year with maximum index
max <- which(sNAO_index_asc == max(sNAO_index_asc), arr.ind = TRUE)
# Year with minimum index
min <- which(sNAO_index_asc == min(sNAO_index_asc), arr.ind = TRUE)
# Year with "normal" index (average index)
normal <- which(sNAO_index_asc == min(abs(sNAO_index_asc)), arr.ind = TRUE)

# Fields for Temperature for years with max, min & normal index
temp_max_year_field <- temp_highSum_mean[lon_range, lat_range, max]
temp_min_year_field <- temp_highSum_mean[lon_range, lat_range, min]
temp_normal_year_field <-   temp_highSum_mean[lon_range, lat_range, normal]

# Fields for Precipitation for years with max, min & normal index
prec_max_year_field <- prec_highSum_mean[lon_range, lat_range, max]
prec_min_year_field <- prec_highSum_mean[lon_range, lat_range, min]
prec_normal_year_field <- prec_highSum_mean[lon_range, lat_range, normal]

# Fields for 500hPa-height for years with max, min & normal index
h500_max_year_field <- h500_highSum_mean[lon_range, lat_range, max]
h500_min_year_field <- h500_highSum_mean[lon_range, lat_range, min]
h500_normal_year_field <- h500_highSum_mean[lon_range, lat_range, normal]

# Fields for Sea-level Pressure for years with max, min & normal index
slp_max_year_field <- slp_highSum_mean[lon_range,lat_range,max]
slp_min_year_field <- slp_highSum_mean[lon_range,lat_range,min]
slp_normal_year_field <- slp_highSum_mean[lon_range,lat_range,normal]

# Fields for Temperature Anomaly for years with max, min & normal index
case_max_temp_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_max_temp_field[i,j] <- temp_max_year_field[i,j] - mean(temp_highSum_mean[i,j,])}}
case_min_temp_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_min_temp_field[i,j] <- temp_min_year_field[i,j] - mean(temp_highSum_mean[i,j,])}}
case_normal_temp_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_normal_temp_field[i,j] <- temp_normal_year_field[i,j] - mean(temp_highSum_mean[i,j,])}}

# Fields for Precipitation Anomaly for years with max, min & normal index
case_max_prec_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_max_prec_field[i,j] <- prec_max_year_field[i,j] - mean(prec_highSum_mean[i,j,])}}
case_min_prec_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_min_prec_field[i,j] <- prec_min_year_field[i,j] - mean(prec_highSum_mean[i,j,])}}
case_normal_prec_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_normal_prec_field[i,j] <- prec_normal_year_field[i,j] - mean(prec_highSum_mean[i,j,])}}

# Fields for 500hPa-height Anomaly for years with max, min & normal index
case_max_h500_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_max_h500_field[i,j] <- h500_max_year_field[i,j] - mean(h500_highSum_mean[i,j,])}}
case_min_h500_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_min_h500_field[i,j] <- h500_min_year_field[i,j] - mean(h500_highSum_mean[i,j,])}}
case_normal_h500_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_normal_h500_field[i,j] <- h500_normal_year_field[i,j] - mean(h500_highSum_mean[i,j,])}}

# Fields for Precipitation Anomaly for years with max, min & normal index
case_max_slp_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_max_slp_field[i,j] <- slp_max_year_field[i,j] - mean(slp_highSum_mean[i,j,])}}
case_min_slp_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_min_slp_field[i,j] <- slp_min_year_field[i,j] - mean(slp_highSum_mean[i,j,])}}
case_normal_slp_field <- array(NA, dim = c(68,32))
for (i in 1:68) {for (j in 1:32) {
  case_normal_slp_field[i,j] <- slp_normal_year_field[i,j] - mean(slp_highSum_mean[i,j,])}}

# Plot for Temperature in year with max index
z_case_max_temp <- rbind(case_max_temp_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                    case_max_temp_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_max_temp), max(z_case_max_temp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_max_temp,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study max index temperature")

# Plot for Temperature in year with min index
z_case_min_temp <- rbind(case_min_temp_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                         case_min_temp_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_min_temp), max(z_case_min_temp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_min_temp,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study min index temperature")

# Plot for Temperature in year with normal index
z_case_normal_temp <- rbind(case_normal_temp_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                       case_normal_temp_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_normal_temp), max(z_case_normal_temp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_normal_temp,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study normal index temperature")

# Plot for Precipitation in year with max index
z_case_max_prec <- rbind(case_max_prec_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                         case_max_prec_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_max_prec), max(z_case_max_prec))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_max_prec,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study max index precipitation")

# Plot for Precipitation in year with min index
z_case_min_prec <- rbind(case_min_prec_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                         case_min_prec_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_min_prec), max(z_case_min_prec))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_min_prec,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study min index precipitation")

# Plot for Precipitation in year with normal index
z_case_normal_prec <- rbind(case_normal_prec_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                            case_normal_prec_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_normal_prec), max(z_case_normal_prec))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_normal_prec,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study normal index precipitation")

# Plot for 500hPa-height in year with max index
z_case_max_h500 <- rbind(case_max_h500_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                         case_max_h500_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_max_h500), max(z_case_max_h500))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_max_h500,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study max index 500hPa-height")

# Plot for 500hPa-heigt in year with min index
z_case_min_h500 <- rbind(case_min_h500_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                         case_min_h500_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_min_h500), max(z_case_min_h500))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_min_h500,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study min index 500hPa-height")

# Plot for 500hPa-heigt in year with normal index
z_case_normal_h500 <- rbind(case_normal_h500_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                            case_normal_h500_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_normal_h500), max(z_case_normal_h500))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_normal_h500,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study normal index 500hPa-height")

# Plot for Sea-level Pressure in year with max index
z_case_max_slp <- rbind(case_max_slp_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                        case_max_slp_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_max_slp), max(z_case_max_slp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_max_slp,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study max index Sea-level Pressure")

# Plot for Sea-level Pressure in year with min index
z_case_min_slp <- rbind(case_min_slp_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                        case_min_slp_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_min_slp), max(z_case_min_slp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_min_slp,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study min index Sea-level Pressure")

# Plot for Sea-level Pressure in year with normal index
z_case_normal_slp <- rbind(case_normal_slp_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                           case_normal_slp_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_normal_slp), max(z_case_normal_slp))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_normal_slp,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study normal index Sea-level Pressure")

####### Plots #######
# Pattern of the sNAO
mycol <- c("blue2", "cornflowerblue", "cadetblue2", "azure2", "bisque1", "burlywood1", "brown1", "brown3")
minmax <- ceiling(max(abs(c(min(z_pc), max(z_pc)))))
mylevs <-  0.25* minmax * (c(1:9) - 5)
par(bg = NA, lwd = 1, bty = "o" )
filled.contour(x_pc, y_pc, z_pc,
               levels = mylevs, col = mycol,
               plot.axes = {
                 axis(1, at = seq(-60, 40, 20), labels = c("60W", "40W", "20W", "0W", "20E", "40E"), font = 2, cex.axis = 1.1);
                 axis(2, at = seq(30, 60, 10), labels = c("30N", "40N", "50N", "60N"), font = 2, cex.axis = 1.1);
                 map("world", interior = F, add = T)},
                 key.axes = axis(4, at = seq(-minmax, minmax, 4), font = 2, cex.axis = 1.5)
)

# Regression sNAO against Temperature
field_temp <- b_temp[2, lon_range, lat_range]
mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <-  seq(-1.25,1.25,0.25)
z_temp <- rbind(field_temp[lon_shortened > 180, rev(1:length(lat_shortened))],
                field_temp[lon_shortened <= 180, rev(1:length(lat_shortened))])
par(bg = NA, lwd = 1, bty = "o")
filled.contour(x_lon, y_lat, z_temp,
               levels = mylevs, col = mycol,
               plot.axes = {
                 axis(1, at = seq(-60, 40, 20), labels = c("60W", "40W", "20W", "0W", "20E", "40E"), font = 2, cex.axis = 1.1);
                 axis(2, at = seq(30, 80, 10), labels = c("30N", "", "50N", "", "70N", ""), font = 2, cex.axis = 1.1);
                 map("world", interior = F, add = T)},
                 key.axes = axis(4, at = c(-1, -0.5, 0, 0.5, 1), font = 2, cex.axis = 1.5)
)

# Regression sNAO against Precipitation
field_prec <- b_prec[2, lon_range, lat_range]
mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <-  seq(-30,30,6)
z_temp <- rbind(field_prec[lon_shortened > 180, rev(1:length(lat_shortened))],
                field_prec[lon_shortened <= 180, rev(1:length(lat_shortened))])
par(bg = NA, lwd = 1, bty = "o")
filled.contour(x_lon, y_lat, z_prec,
               levels = mylevs, col = rev(mycol),
               plot.axes = {
                 axis(3, at = seq(-60, 40, 20), labels = c("", "", "", "", "", ""), font = 2, cex.axis = 1.1);
                 axis(2, at = seq(30, 80, 10), labels = c("30N", "", "50N", "", "70N", ""), font = 2, cex.axis = 1.1);
                 map("world", interior = F, add = T)},
                 key.axes = axis(4, at = c(-24, -12, 0, 12, 24), font = 2, cex.axis = 1.5)
)

# Case study pattern precipitation year with max index
z_case_max_prec <- rbind(case_max_prec_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                         case_max_prec_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <-  c(-250,seq(-200, 200, 50),250)
par(bg = NA, lwd = 1, bty = "o")
filled.contour(x_lon, y_lat, z_case_max_prec,
               levels = mylevs, col = rev(mycol),
               plot.axes = {
                 axis(1, at = seq(-60, 40, 20), labels = c("60W", "40W", "20W", "0W", "20E", "40E"), font = 2, cex.axis = 1.1);
                 axis(2, at = seq(30, 80, 10), labels = c("30N", "", "50N", "", "70N", ""), font = 2, cex.axis = 1.1);
                 map("world", interior = F, add = T)},
                 key.axes = axis(4, at = seq(-200, 200, 100),
                               labels = c("<-200", "-100", "0", "100", ">200")
                               , font = 2, cex.axis = 1.4)
)

# Case study pattern precipitation year with min index
z_case_min_prec <- rbind(case_min_prec_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                         case_min_prec_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
mycol <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
mylevs <-  c(-250,seq(-200, 200, 50),250)
par(bg = NA, lwd = 1, bty = "o")
filled.contour(x_lon, y_lat, z_case_min_prec,
               levels = mylevs, col = rev(mycol),
               plot.axes = {
                 axis(3, at = seq(-60, 40, 20), labels = c("", "", "", "", "", ""), font = 2, cex.axis = 1.1);
                 axis(2, at = seq(30, 80, 10), labels = c("30N", "", "50N", "", "70N", ""), font = 2, cex.axis = 1.1);
                 map("world", interior = F, add = T)},
                 key.axes = axis(4, at = seq(-200, 200, 100),
                               labels = c("<-200", "-100", "0", "100", ">200")
                               , font = 2, cex.axis = 1.4)
)

# Plot for Precipitation in year with min index
z_case_min_prec <- rbind(case_min_prec_field[lon_shortened > 180, rev(1:length(lat_shortened))],
                         case_min_prec_field[lon_shortened <= 180, rev(1:length(lat_shortened))])
minmax <- max(abs(c(min(z_case_min_prec), max(z_case_min_prec))))
mylevs <-  0.25 * minmax * (c(1:9) - 5)
filled.contour(x_lon, y_lat, z_case_min_prec,
               levels = mylevs, col = mycol,
               plot.axes = {
                 map("world", interior = F, add = T)
               }, main = "Case study min index precipitation")

# Timeseries sNAO, temperature & precipitation
par(bg = NA, bty = "L")
plot(sNAO_index_asc_plot,
     type = "l", lwd = 4, lty = 1, col = "black",
     main = "", xlab = "", ylab = "",
     ylim = c(min(c(sNAO_index_asc_plot, case_temp_plot, case_prec_plot*-1)), max(c(sNAO_index_asc_plot, case_temp_plot, case_prec_plot*-1))),
     xaxt = "none", yaxt = "none")
legend(x = 6.8, y = 4.1, legend = c("sNAO index", "Precipitation", "Temperature"),
       cex = 1.4, col = c("black", "blue2", "brown3"), lty = 1, bty = "n", lwd = 8, y.intersp = 0.45,
       x.intersp = 0.3, seg.len = 0.5)
axis(1, at = seq(start-start, end-start, 10),
     labels = c("1950", "1960", "1970", "1980", "1990", "2000"),
     font = 2, cex.axis = 1.4)
axis(2, font = 2, cex.axis = 1.4 )
mtext(side=1, line=2.5, "Year", font=2,cex=1.4)
mtext(side=2, line=2.5, "Standardized", font=2,cex=1.4)
lines(case_prec_plot*-1, col = "blue2", lwd = 2.5) # Precipitation and sNAO index correlate negatively
lines(case_temp_plot, col = "brown3", lwd = 2.5 )

# Separate negative and positive loadings to plot them differently
z_pc_title_pos <- array(NA, dim = c(64,24))
for (i in 1:64) {for (j in 1:24) {
  z_pc_title_pos[i,j] <- 
    (if (z_pc[i, j] >= 0) {
    z_pc[i, j]
  } else {
    0
  })
}}

z_pc_title_neg <- array(NA, dim = c(64,24))
for (i in 1:64) {for (j in 1:24) {
  z_pc_title_neg[i,j] <- (if (z_pc[i,j] < 0) {z_pc[i,j]}
                          else {NA})
}}

# Pattern of the sNAO (Title)
minmax <- ceiling(max(abs(c(min(z_pc), max(z_pc)))))
mylevs <-  0.25* minmax * (c(1:9) - 5)
par(bg = NA, fg = NA, bty = "n")
filled.contour(x_pc, y_pc, z_pc,
               levels = mylevs, col = "transparent",
               plot.axes = {
                 axis(1, at = 180); axis(2, at = FALSE);
                 c(map("world", interior = F, add = T, lwd = 1, col = "white"),
                   contour(x_pc, y_pc, z_pc_title_pos,
                           levels = mylevs,
                           #labels = c("0.05", "0.1", "0.15", "0.2"),
                           col = "white", drawlabels = F,
                           axes = F, frame.plot = F,      #?
                           add = TRUE, lwd = 4, lty = 6),
                   contour(x_pc, y_pc, z_pc_title_neg,
                           levels = mylevs,
                           #labels = c("0.05", "0.1", "0.15", "0.2"),
                           col = "white", drawlabels = F,
                           axes = F, frame.plot = F,      #?
                           add = TRUE, lwd = 4, lty = 1))
               }, key.axes = FALSE
)

# Timeseries of the sNAO index (Title)
time_offset = 12 # One year = 12 months
years <- 6
par(bg=NA,fg=NA)
plot(NAO_index_asc[(length(NAO_index_asc) - years * time_offset):length(NAO_index_asc)],
     type = "h", main = "", xlab = "", ylab = "", col="white",fg = "transparent",
     bg = "transparent", col.axis = "transparent", col.lab = "transparent", lwd=4)
dev.copy(png,"myplot.png")
dev.off()
