######################
# Data Preprocessing #
######################
data(GEF14)
GEF14$doy <- unclass(as.POSIXlt(GEF14$date))$yday + 1
GEF14$month <- unclass(as.POSIXlt(GEF14$date))$mon + 1
GEF14 <- GEF14[-c(1 : 47, 60528),]

d <- 24
nobs <- nrow(GEF14)
GEF14_load <- GEF14_load24 <- GEF14_temp <- GEF14_temp95 <- matrix(0, nobs/d, d)
GEF14_doy <- rep(0, nobs/d)

for(i in unique(GEF14$tod)){
  GEF14_load[, i + 1] <- GEF14$load[GEF14$tod == i]
  GEF14_load24[, i + 1] <- GEF14$load24[GEF14$tod == i]
  GEF14_temp[, i + 1] <- GEF14$temp[GEF14$tod == i]
  GEF14_temp95[, i + 1] <- GEF14$temp95[GEF14$tod == i]
}
GEF14_data <- as.data.frame(cbind(GEF14_load, GEF14_load24, GEF14_temp, GEF14_temp95))

colnames(GEF14_data)[1 : d] <- paste0("load_h", 0 : 23)
colnames(GEF14_data)[(d + 1) : (2 * d)] <- paste0("load24_h", 0 : 23)
colnames(GEF14_data)[(2 * d + 1) : (3 * d)] <- paste0("temp_h", 0 : 23)
colnames(GEF14_data)[(3 * d + 1) : (4 * d)] <- paste0("temp95_h", 0 : 23)

GEF14_data$doy <- GEF14$doy[GEF14$tod == 0]
GEF14_data$dow <- GEF14$dow[GEF14$tod == 0]
GEF14_data$month <- GEF14$month[GEF14$tod == 0]
GEF14_data$year <- GEF14$year[GEF14$tod == 0]
GEF14_data$progtime <- 1 : nrow(GEF14_data)

##############################################################################################
