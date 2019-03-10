#load libraries
library(tidyverse)
library(stringr)

#read data (Estonian traffic violations data)
df_src <- read_csv('liiklusjarelevalve_1_LEST.csv')

#set 1st converstion function
rad2deg <- function(rad) {(rad * 180) / (pi)}

#set 2nd conversion function
lest2geo <- function(x,y) {
    a <- 6378137.00000000000
    F <- 1/298.257222100883
    ESQ <- (F + F - F * F)
    B0 <- ((57.00000000000 + 31.0000000000 / 60.000000000000 + 3.19414800000 / 3600.00000000000) / rad2deg(1))
    L0 <- (24.00000000000 / rad2deg(1))
    FN <- 6375000.00000000000
    FE <- 500000.00000000000
    B2 <- ((59.00000000000 + 20.00000000000 / 60.00000000000) / rad2deg(1))
    B1 <- (58.00000000000 / rad2deg(1))
    xx <- (x - FN)
    yy <- (y - FE)
    t0 <- sqrt((1.00000000000 - sin(B0)) / (1.00000000000 + sin(B0)) * `^`(((1.00000000000 + sqrt(ESQ) * sin(B0)) / (1.00000000000 - sqrt(ESQ) * sin(B0))), sqrt(ESQ)))
    t1 <- sqrt((1.00000000000 - sin(B1)) / (1.00000000000 + sin(B1)) * `^`(((1.00000000000 + sqrt(ESQ) * sin(B1)) / (1.00000000000 - sqrt(ESQ) * sin(B1))), sqrt(ESQ)))
    t2 <- sqrt((1.00000000000 - sin(B2)) / (1.00000000000 + sin(B2)) * `^`(((1.00000000000 + sqrt(ESQ) * sin(B2)) / (1.00000000000 - sqrt(ESQ) * sin(B2))), sqrt(ESQ)))
    m1 <- (cos(B1) / `^`((1.00000000000 - ESQ * sin(B1) * sin(B1)), 0.50000000000))
    m2 <- (cos(B2) / `^`((1.00000000000 - ESQ * sin(B2) * sin(B2)), 0.50000000000))
    n1 <- ((log(m1) - log(m2)) / (log(t1) - log(t2)))
    FF <- (m1 / (n1 * `^`(t1, n1)))
    p0 <- (a * FF * `^`(t0, n1))
    p <- `^`((yy * yy + (p0 - xx) * (p0 - xx)), 0.50000000000)
    t <- `^`((p / (a * FF)), (1.00000000000 / n1))
    FII <- atan(yy / (p0 - xx))
    LON <- (FII / n1 + L0)
    u <- ((pi / 2.00000000000) - (2.00000000000 * atan(t)))
    LAT <- (u + (ESQ / 2.00000000000 + (5.00000000000 * `^`(ESQ, 2) / 24.00000000000) + (`^`(ESQ, 3) / 12.00000000000) +
     (13.00000000000 * `^`(ESQ, 4) / 360.00000000000)) * sin(2.00000000000 * u) +
     ((7.00000000000 * `^`(ESQ, 2) / 48.00000000000) + (29.00000000000 * `^`(ESQ, 3) / 240.00000000000) +
     (811.00000000000 * `^`(ESQ, 4) / 11520.00000000000)) * sin(4.00000000000 * u) +
     ((7.00000000000 * `^`(ESQ, 3) / 120.00000000000) + (81.00000000000 * `^`(ESQ, 4) / 1120.00000000000)) * sin(6.00000000000 * u) +
     (4279.00000000000 * `^`(ESQ, 4) / 161280.00000000000) * sin(8.00000000000 * u))
    LAT <- rad2deg(LAT);
    LON <- rad2deg(LON);
    #LON = LON * rad2dedegrees(1);
    paste(LAT,LON, sep = ',')
}  

# use functions to convert LEST coordinates to Lat/Lon
# create a copy of the original datasource
df_fin <- df_src
# string replacement on LEST coordinates so R can coerce the chr into num
df_fin$Lest_X <- str_replace(string = df_fin$Lest_X, pattern = '-', replacement = '.') 
df_fin$Lest_Y <- str_replace(string = df_fin$Lest_Y, pattern = '-', replacement = '.')
# coerce chr into num
df_fin$Lest_X <- as.numeric(df_fin$Lest_X)
df_fin$Lest_Y <- as.numeric(df_fin$Lest_Y)
# create a new column and populate it with the results from lest2geo
#df_fin <- mutate(.data = df_fin, LatLon = lest2geo(df_fin$Lest_X,df_fin$Lest_Y))
df_fin$LatLon <- lest2geo(df_fin$Lest_X, df_fin$Lest_Y)
# split new LatLon column into two columns 
df_fin <- separate(data = df_fin, col = LatLon, into = c('Latitude','Longitude'), sep =',', remove = TRUE)
# as the output of lest2geo is str, convert new split columns to num
df_fin$Latitude <- as.numeric(df_fin$Latitude)
df_fin$Longitude <- as.numeric(df_fin$Longitude)
# clean up the output columns
df_fin <- select(.data = df_fin, JuhtumId,Latitude,Longitude)
str(df_fin)
head(df_fin)
