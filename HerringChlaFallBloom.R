################################################################################
########LINKING HERRING RECRUITMENT TO PHYTOPLANKTON BLOOMS#####################
################################################################################

#Project with Rafa and Brian

#Required packages
#install.packages("R.matlab")
library(R.matlab)
library(dplyr)
library(ggplot2)
library(reshape2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#Read satellite .mat file into R
mat_data <- readMat("OC_timeseries_NorthSea_Giacomo.mat")

# view list of objects in .mat file
names(mat_data)

#Rename objects in .mat file
CHLW <- mat_data$CHLW
CHLWsmall <- mat_data$CHLW.small
CHLWavg <- mat_data$CHLWavg
CHLWavgsmall <- mat_data$CHLWavg.small
DATE <- mat_data$DATEW
LAT <- mat_data$LAT
LONG <- mat_data$LON
LONGsmall <- mat_data$LON.small

#Create data frame with avg CHLW and date
chldate <- cbind(mat_data$CHLWavg,mat_data$DATEW)
chldate <- data.frame(Chl=chldate[,1], Year=chldate[,2], Day = chldate[,4],
                      Month = chldate[,3])
chldatena <- na.omit(chldate)

#Chlorophyll a conc. label µg L-1
chlaconc <- expression(paste("Chlorophyll a ( ",mu,"g ",L^-1,")"))

#Load herring recruitment data
herr <- read.csv("herring.csv")
her <- herr[,c("Year","Recruitment","SSB")]
her$rec1 <- c(her$Recruitment[2:76], NA)
her$recssb <- her$rec1/her$SSB
her <- na.omit(her)
vecrec <- her$recssb[52:75]

#save(her, file = "her.RData")
#write.csv(her, "NSherringRecSSB.csv", row.names=FALSE)


################ Analysis 1: Max Chla of the year #############################

chldatena <- subset(chldatena, Year!=2022 & Year!=1997)

#Create data frame with chl max of fall bloom
hchl <- chldatena %>% group_by(Year) %>% summarize(Chlmax = max(Chl))
hchl$recssb1 <- vecrec

#Linear model
mchlh <- lm(recssb1 ~ Chlmax, hchl)
summary(mchlh)
par(mfrow=c(2,2))
plot(mchlh)

#Time series plot
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(hchl$Year, hchl$Chlmax, type='l', lwd=3, col = "green", xlab = "Year", ylab = "")
axis(2, col = "green", col.ticks = "green", col.axis = "green", col.lab="green")
mtext("Maximum annual Chl concentration [/]", side = 2, line = 3, col = "green") 
par(new = TRUE) 
plot(hchl$Year, hchl$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext("ln(Recruitment/SSB [thousands/tonnes])", side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "ln(Rec/SSB)"), col = c("green", "blue"), lty = 1, cex = 0.8)
#800 vs 450 aspect ratio

#Plot model results
newchl <- data.frame("Chlmax"=seq(min(hchl$Chlmax),max(hchl$Chlmax),length=1000))
confrec <- predict(mchlh,newchl, interval="confidence")
confrec

plot(hchl$Chlmax,hchl$recssb1, ylim=c(6,45), xlab = "Maximum Chl concentration [/]", ylab = "ln(Recruitment/SSB [thousands/tonnes])")
text(hchl$Chlmax,hchl$recssb1, labels = hchl$Year, pos = 3, cex=0.5)
matlines(newchl,confrec[,1:3],lty=c(1,2,2), col=c(1,3,3),lwd=c(1.5,1.5,1.5))
#650 450

### #### #### Add week number
hchl$nweek_max_year <- rep(0,nrow(hchl))

for (i in 1:length(yphen)) {
  nwmhchl <- which.max(chldatena$Chl[chldatena$Year==yphen[i]])
  hchl$nweek_max_year[hchl$Year==yphen[i]] <- chldatena$nweek[chldatena$Year==yphen[i]][nwmhchl]
}

hchl

#no relationship with year or with herring rec
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(hchl$Year, hchl$nweek_max_year, type='l', lwd=2, col = "grey", xlab = "Year", ylab = "")
axis(2, col = "grey", col.ticks = "grey", col.axis = "grey", col.lab="grey")
mtext("Week number at which yearly max is observed", side = 2, line = 3, col = "grey") 
par(new = TRUE) 
plot(hchl$Year, hchl$recssb1, type='l', lwd=2, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext("ln(Recruitment/SSB [thousands/tonnes])", side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "ln(Rec/SSB)"), col = c("grey", "blue"), lty = 1, lwd=2, cex = 0.8, bty="n")


################ Analysis 2: Max Chla of the fall #############################

#Calculate start of bloom for every year
#Calculate year median chl
chldatena$medians <- rep(0, length(chldatena$Chl))
mediansy <- chldatena %>% group_by(Year) %>% summarize(median = median(Chl))
str(chldatena)
chldatena$YearF <- as.factor(chldatena$Year)

chldatena$start <- rep("No", length(chldatena$Chl))

#Put medians in
yrs <- levels(chldatena$YearF)
#yrs <- yrs[-1]
for (i in 1:length(yrs)){
  chldatena <- chldatena %>%
    mutate( medians = if_else(YearF == yrs[i], mediansy$median[mediansy$Year==yrs[i]], medians))
  
}

#chl above yearly median?
for (i in 1:length(yrs)){
for (j in 1:nrow(chldatena)){
      if (chldatena$YearF[j]==yrs[i] & chldatena$Chl[j] >= max(chldatena$medians[chldatena$YearF==yrs[i]]*1.05) )
        chldatena$start[j] <- "Yes"
    }
}

#Subset for fall
chldatenafall <- subset(chldatena, Month>=8)


#Create data frame with chl max of fall bloom
hchlfall <- chldatenafall %>% group_by(Year) %>% summarize(Chlmax = max(Chl))
hchlfall <- subset(hchlfall, Year!=2022)
hchlfall$recssb1 <- vecrec

#Plot of time series
chlalab <- expression(paste("Maximum fall bloom Chl-a ( ",mu," g ",L^-1,")"))
rssblab <- expression("ln(Recruitment/SSB [thousands/tonnes])")
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(hchlfall$Year, hchlfall$Chlmax, type='l', lwd=3, col = "green", xlab = "Year", ylab = "", main="Entire North Sea")
axis(2, col = "green", col.ticks = "green", col.axis = "green", col.lab="green")
mtext(chlalab, side = 2, line = 3, col = "green") 
par(new = TRUE) 
plot(hchlfall$Year, hchlfall$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext(rssblab, side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "ln(R/SSB)"), col = c("green", "blue"), lty = 1, cex = 0.8, bty="n")
#800 vs 450 aspect ratio
#450 420

#Linear model
mchlhfall <- lm(recssb1 ~ Chlmax, hchlfall)
summary(mchlhfall)
par(mfrow=c(2,2))
plot(mchlhfall)

#######################  Analysis 3: Mean chl of the year  #####################

meanchl <- chldatena %>% group_by(Year) %>% summarize(Chlmean = mean(Chl))
meanchl$recssb1 <- vecrec

#Plot
par(mfrow=c(1,1))
plot(meanchl$Chlmean, meanchl$recssb1)

#Plot of time series
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(meanchl$Year, meanchl$Chlmean, type='l', lwd=3, col = "green", xlab = "Year", ylab = "")
axis(2, col = "green", col.ticks = "green", col.axis = "green", col.lab="green")
mtext("Mean annual Chl concentration [/]", side = 2, line = 3, col = "green") 
par(new = TRUE) 
plot(meanchl$Year, meanchl$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext("ln(Recruitment/SSB [thousands/tonnes])", side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "ln(Rec/SSB)"), col = c("green", "blue"), lty = 1, cex = 0.8)
#800 vs 450 aspect ratio


#######################  Analysis 4: Mean chl of the fall bloom  #####################

chldatenafall

hchlmfall <- subset(chldatenafall, start=="Yes")
hchlmfall <- hchlmfall %>% group_by(Year) %>% summarize(Chlmean = mean(Chl))
hchlmfall$recssb1 <- vecrec

#Plot of time series
meanchlalab <- expression(paste("Mean fall bloom Chl-a ( ",mu," g ",L^-1,")"))
meanrssblab <- expression("ln(Recruitment/SSB [thousands/tonnes])")
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(hchlmfall$Year, hchlmfall$Chlmean, type='l', lwd=3, col = "green", xlab = "Year", ylab = "", main="Entire North Sea")
axis(2, col = "green4", col.ticks = "green", col.axis = "green", col.lab="green")
mtext(meanchlalab, side = 2, line = 3, col = "green") 
par(new = TRUE) 
plot(hchlmfall$Year, hchlmfall$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext(meanrssblab, side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "ln(R/SSB)"), col = c("green", "blue"), lty = 1, cex = 0.8, bty="n")
#800 vs 450 aspect ratio



##################################For the small box#############################

#Create chl data frame
CHLWavgsmall
dates <- mat_data$DATEW
df <- data.frame(year = dates[,1], month = dates[,2], day = dates[,3])
df$date_var <- as.Date(paste(df$year, df$month, df$day, sep = "-"), format = "%Y-%m-%d")
df

chldatesmall <- data.frame(Chl=CHLWavgsmall, date=df$date_var, Year=df$year, month=df$month)
chldatenasmall <- na.omit(chldatesmall)

chldatenasmall <- subset(chldatenasmall, Year!=2022 & Year!=1997)

#Create chl max of fall bloom data frame
hchlsmall <- chldatenasmall %>% group_by(Year) %>% summarize(Chlmax = max(Chl))
hchlsmall$recssb1 <- vecrec

hchlsmall <- subset(hchlsmall, Year!=1997)

#Statistical model
mchlhsmall <- lm(recssb1 ~ Chlmax, hchlsmall)
summary(mchlhsmall)

plot(hchlsmall$Chlmax,hchlsmall$recssb1, ylim=c(6,45))

text(hchlsmall$Chlmax,hchlsmall$recssb1, labels = hchlsmall$Year, pos = 3, cex=0.5)

######## Fall maximum 

chldatenafallsmall <- subset(chldatenasmall, month>=8)

hchlsmallfall <- chldatenafallsmall %>% group_by(Year) %>% summarize(Chlmax = max(Chl))
hchlsmallfall$recssb1 <- vecrec

meanrssblab <- expression("ln(Recruitment/SSB [thousands/tonnes])")
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(hchlsmallfall$Year, hchlsmallfall$Chlmax, type='l', lwd=3, col = "green3", xlab = "Year", ylab = "", main="Western North Sea")
axis(2, col = "green3", col.ticks = "green3", col.axis = "green3", col.lab="green3")
mtext(chlalab, side = 2, line = 3, col = "green3") 
par(new = TRUE) 
plot(hchlsmallfall$Year, hchlsmallfall$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext(meanrssblab, side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "ln(R/SSB)"), col = c("green3", "blue"), lty = 1, cex = 0.8, bty="n")


####Annual mean

meanchlsmall <- chldatenasmall %>% group_by(Year) %>% summarize(Chlmean = mean(Chl))
meanchlsmall <- subset(meanchlsmall, Year!=1997 & Year!=2022)
meanchlsmall$recssb1 <- vecrec

#Plot
par(mfrow=c(1,1))
plot(meanchlsmall$Chlmean, meanchlsmall$recssb1)

#Plot of time series
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(meanchlsmall$Year, meanchlsmall$Chlmean, type='l', lwd=3, col = "green", xlab = "Year", ylab = "")
axis(2, col = "green", col.ticks = "green", col.axis = "green", col.lab="green")
mtext("Mean annual Chl concentration [/]", side = 2, line = 3, col = "green") 
par(new = TRUE) 
plot(meanchlsmall$Year, meanchlsmall$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext("ln(Recruitment/SSB [thousands/tonnes])", side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "ln(Rec/SSB)"), col = c("green", "blue"), lty = 1, cex = 0.8)
#800 vs 450 aspect ratio



####fall mean

#Calculate start of bloom for every year
#Calculate year median chl
chldatenasmall$medians <- rep(0, length(chldatenasmall$Chl))
mediansysmall <- chldatenasmall %>% group_by(Year) %>% summarize(median = median(Chl))
chldatenasmall$YearF <- as.factor(chldatenasmall$Year)

chldatenasmall$start <- rep("No", length(chldatenasmall$Chl))

#Put medians in
yrs <- levels(chldatenasmall$YearF)
#yrs <- yrs[-1]
for (i in 1:length(yrs)){
  chldatenasmall <- chldatenasmall %>%
    mutate( medians = if_else(YearF == yrs[i], mediansysmall$median[mediansysmall$Year==yrs[i]], medians))
  
}

#chl above yearly median?
for (i in 1:length(yrs)){
  for (j in 1:nrow(chldatenasmall)){
    if (chldatenasmall$YearF[j]==yrs[i] & chldatenasmall$Chl[j] >= max(chldatenasmall$medians[chldatenasmall$YearF==yrs[i]]*1.05) )
      chldatenasmall$start[j] <- "Yes"
  }
}

#Subset for fall
chldatenafallsmall <- subset(chldatenasmall, month>=8)

hchlmfallsmall <- subset(chldatenafallsmall, start=="Yes")
hchlmfallsmall <- hchlmfallsmall %>% group_by(Year) %>% summarize(Chlmean = mean(Chl))
hchlmfallsmall$recssb1 <- vecrec


#Plot
par(mfrow=c(1,1))
plot(hchlmfallsmall$Chlmean, hchlmfallsmall$recssb1)

#Plot of time series ######### mean fall small box
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(hchlmfallsmall$Year, hchlmfallsmall$Chlmean, type='l', lwd=3, col = "green3", xlab = "Year", ylab = "", main="Western North Sea")
axis(2, col = "green3", col.ticks = "green3", col.axis = "green3", col.lab="green3")
mtext(meanchlalab, side = 2, line = 3, col = "green3") 
par(new = TRUE) 
plot(hchlmfallsmall$Year, hchlmfallsmall$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext(meanrssblab, side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "ln(R/SSB)"), col = c("green3", "blue"), lty = 1, cex = 0.8, bty="n")
#800 vs 450 aspect ratio







####################################Creating spatial map#########################

#Average matrix large cell
dim(CHLW[1, ,])[1]

avg_mat <- matrix(0, nrow = nrow(CHLW[1, ,]), ncol=ncol(CHLW[1, ,]))

for (i in 1:nrow(CHLW[1, ,])) {
  for (j in 1:ncol(CHLW[1, ,])) {
    avg_mat[i, j] <- mean(CHLW[1:1319, i, j], na.rm = TRUE)
  }
}


#Add latitudes and longitudes
colnames(avg_mat) <- LONG
rownames(avg_mat) <- LAT

# Reshape matrix into data frame
avgchlmap <- melt(avg_mat, varnames = c("latitude", "longitude"), value.name = "chl")

# View resulting data frame
ggplot() +
  geom_tile(data = avgchlmap, aes(x = longitude, y = latitude, fill = sqrt(chl))) +
  scale_fill_gradient(low = "white", high = "green") +
  theme_void()

#####Map with countries

# Define the North Sea box
northsea_bbox <- st_bbox(c(xmin = min(avgchlmap$longitude), ymin = min(avgchlmap$latitude), xmax = max(avgchlmap$longitude), ymax = max(avgchlmap$latitude)),crs = st_crs(4326))
northsea_bbox <- as.array(northsea_bbox)

# Get the countries around the North Sea from the `rnaturalearth` package
northsea_countries <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid() %>%
  st_intersection(st_as_sfc(northsea_bbox))

# Get the landmasses from the `rnaturalearthdata` package
landmasses <- ne_download(type = "land", category = "physical", returnclass = "sf")

###Good plot
ggplot() +
  geom_tile(data = avgchlmap, aes(x = longitude, y = latitude, fill = log10(chl))) +
  scale_fill_gradient(low = "white", high = "green") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries, color = "black") +
  coord_sf(xlim = c(northsea_bbox$xmin, northsea_bbox$xmax), ylim = c(northsea_bbox$ymin, northsea_bbox$ymax), expand = FALSE) +
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")
#670 x 420

log10chlalab <- expression(paste(log[10],"(Chl-a (",mu,"g ",L^-1,"))"))
###Good plot2
ggplot() +
  geom_tile(data = avgchlmap, aes(x = longitude, y = latitude, fill = log10(chl))) +
  scale_fill_gradientn(colors = topo.colors(50)[5:45], name = log10chlalab) +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countries, color = "black") +
  coord_sf(xlim = c(northsea_bbox$xmin, northsea_bbox$xmax), ylim = c(northsea_bbox$ymin, northsea_bbox$ymax), expand = FALSE) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "grey"))+
  labs(x = "Longitude (°)", y = "Latitude (°)")
#750 x 420

#####
# create the breaks- and label vectors
xcoor <- seq(-3,northsea_bbox$xmax,5)
ycoor <- seq(northsea_bbox$ymin,northsea_bbox$ymax,5)
ewlbls <- unlist(lapply(xcoor, function(x) ifelse(x < 0, paste(x, "°E"), ifelse(x > 0, paste(x, "°W"),x))))
nslbls <- unlist(lapply(ycoor, function(x) ifelse(x < 0, paste(x, "°S"), ifelse(x > 0, paste(x, "°N"),x))))

scale_x_continuous(breaks = xcoor, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = ycoor, labels = nslbls, expand = c(0, 0)) +
  theme(axis.text = element_text(size=12))




#####

#Average matrix for small cell
dim(CHLWsmall[1, ,])[1]

avg_mat_small <- matrix(0, nrow = nrow(CHLWsmall[1, ,]), ncol=ncol(CHLWsmall[1, ,]))

for (i in 1:nrow(CHLWsmall[1, ,])) {
  for (j in 1:ncol(CHLWsmall[1, ,])) {
    avg_mat_small[i, j] <- mean(CHLWsmall[1:1319, i, j], na.rm = TRUE)
  }
}

#Add latitudes and longitudes
colnames(avg_mat_small) <- LONGsmall
rownames(avg_mat_small) <- LAT

# Reshape matrix into data frame
avgchlmapsmall <- melt(avg_mat_small, varnames = c("latitude", "longitude"), value.name = "chl")

# View resulting data frame
ggplot() +
  geom_tile(data = avgchlmapsmall, aes(x = longitude, y = latitude, fill = log10(chl))) +
  scale_fill_gradient(low = "white", high = "green")+
  theme_minimal()


#####Map with countries

# Define the North Sea box
northsea_bboxs <- st_bbox(c(xmin = min(avgchlmapsmall$longitude), ymin = min(avgchlmapsmall$latitude), xmax = max(avgchlmapsmall$longitude), ymax = max(avgchlmapsmall$latitude)),crs = st_crs(4326))
northsea_bboxs <- as.array(northsea_bboxs)

# Get the countries around the North Sea from the `rnaturalearth` package
northsea_countriess <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid() %>%
  st_intersection(st_as_sfc(northsea_bboxs))

# Get the landmasses from the `rnaturalearthdata` package
landmasses <- ne_download(type = "land", category = "physical", returnclass = "sf")

###Good plot
ggplot() +
  geom_tile(data = avgchlmapsmall, aes(x = longitude, y = latitude, fill = log10(chl))) +
  scale_fill_gradient(low = "white", high = "green") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countriess, color = "black") +
  coord_sf(xlim = c(northsea_bboxs$xmin, northsea_bboxs$xmax), ylim = c(northsea_bboxs$ymin, northsea_bboxs$ymax), expand = FALSE) +
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")


#Just right cell

#Create average matrix
avg_mat_large <- avg_mat[,145:288]

#Add latitudes and longitudes
colnames(avg_mat_large) <- LONG[145:288]
rownames(avg_mat_large) <- LAT

# Reshape matrix into data frame
avgchlmapslarge <- melt(avg_mat_large, varnames = c("latitude", "longitude"), value.name = "chl")


# Define the North Sea box
northsea_bboxl <- st_bbox(c(xmin = min(avgchlmapslarge$longitude), ymin = min(avgchlmapslarge$latitude), xmax = max(avgchlmapslarge$longitude), ymax = max(avgchlmapslarge$latitude)),crs = st_crs(4326))
northsea_bboxl <- as.array(northsea_bboxl)

# Get the countries around the North Sea from the `rnaturalearth` package
northsea_countriesl <- ne_countries(scale = "medium", returnclass = "sf") %>% st_make_valid() %>%
  st_intersection(st_as_sfc(northsea_bboxl))

# Get the landmasses from the `rnaturalearthdata` package
landmasses <- ne_download(type = "land", category = "physical", returnclass = "sf")

###Good plot
ggplot() +
  geom_tile(data = avgchlmapslarge, aes(x = longitude, y = latitude, fill = log10(chl))) +
  scale_fill_gradient(low = "white", high = "green") +
  geom_sf(data = landmasses, fill = "grey90") +
  geom_sf(data = northsea_countriesl, color = "black") +
  coord_sf(xlim = c(northsea_bboxl$xmin, northsea_bboxl$xmax), ylim = c(northsea_bboxl$ymin, northsea_bboxl$ymax), expand = FALSE) +
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")












CHLWavgsmall <- mat_data$CHLWavg.small
#small box
chldatesmall <- cbind(mat_data$CHLWavg.small,mat_data$DATEW)
chldatesmall <- data.frame(Chl=chldatesmall[,1], Year=chldatesmall[,2], Day = chldatesmall[,4], Month = chldatesmall[,3])
chldatenasmall <- na.omit(chldatesmall)

#calculate start of bloom for every year
mediansysmall <- chldatenasmall %>% group_by(Year) %>% summarize(median = median(Chl))
str(chldatenasmall)
chldatenasmall$YearF <- as.factor(chldatenasmall$Year)

chldatenasmall$medians <- rep(0, length(chldatenasmall$Chl))
chldatenasmall$start <- rep("No", length(chldatenasmall$Chl))

#Check and remove years 1997 and 2022 as they may not have enought data
#remove year 1997, only afeter mid november
chldatenasmall <- subset(chldatenasmall, YearF!="1997")
chldatena2022small <- subset(chldatenasmall, YearF=="2022")
#year 2022 is ok

#Put medians in
for (i in 1:length(yrs)){
  chldatenasmall <- chldatenasmall %>%
    mutate( medians = if_else(YearF == yrs[i], mediansysmall$median[mediansysmall$Year==yrs[i]], medians))
  
}

yrs1 <- levels(chldatenasmall$YearF)
#yrs1 <- yrs1[-1] #remove 1997

#label bloom
#chldatena2021 <- subset(chldatena, YearF=="2021")
#chldatena2021 <- chldatena2021 %>%
#  mutate( start = if_else(Chl>medians*1.05, "Yes", start))


for (i in 1:length(yrs1)){
  for (j in 1:nrow(chldatenasmall)){
    if (chldatenasmall$YearF[j]==yrs1[i] & chldatenasmall$Chl[j] >= max(chldatenasmall$medians[chldatenasmall$YearF==yrs1[i]]*1.05) )
      chldatenasmall$start[j] <- "Yes"
  }
}

#subset for 
chldatenafallsmall <- subset(chldatenasmall, Month>=8)


chldatenafallsmall$recssb1 <- rep(0, length(chldatenafallsmall$Chl))

chldatenafallsmall$YearF <- as.factor(chldatenafallsmall$Year)
her$YearF <- as.factor(her$Year)
yrs <- levels(chldatenafallsmall$YearF)

for (l in 1:length(yrs)){
  chldatenafallsmall <- chldatenafallsmall %>%
    mutate(recssb1 = if_else(Year == yrs[l], her$recssb[her$Year==yrs[l]], recssb1))
  
}

hchlfallsmall <- chldatenafallsmall %>% group_by(Year) %>% summarize(Chlmax = max(Chl), recssb1=mean(recssb1))
hchlfallsmall <- na.omit(hchlfallsmall)
hchlfallsmall_dat <- data.frame(Year=hchlfallsmall$Year, Chlmax=hchlfallsmall$Chlmax, recssb1=hchlfallsmall$recssb1)
hchlfallsmall_dat <- subset(hchlfallsmall_dat, Year!="2022")

mchlhfallsmall <- lm(recssb1 ~ Chlmax, hchlfallsmall_dat)
summary(mchlhfallsmall)

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(hchlfallsmall_dat$Year, hchlfallsmall_dat$Chlmax, type='l', lwd=3, col = "green", xlab = "Year", ylab = "")
axis(2, col = "green", col.ticks = "green", col.axis = "green", col.lab="green")
mtext("Maximum Chl concentration [/]", side = 2, line = 3, col = "green") 
par(new = TRUE) 
plot(hchlfallsmall_dat$Year, hchlfallsmall_dat$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext("Recruitment/SSB [thousands/tonnes]", side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "Rec/SSB"), col = c("green", "blue"), lty = 1, cex = 0.8)
#800 vs 450 aspect ratio

newchlfallsmall <- data.frame("Chlmax"=seq(min(hchlfallsmall_dat$Chlmax),max(hchlfallsmall_dat$Chlmax),length=1000))
confrecfallsmall <- predict(mchlhfallsmall,newchlfallsmall, interval="confidence")
confrecfallsmall

plot(hchlfallsmall_dat$Chlmax,hchlfallsmall_dat$recssb1, ylim=c(6,50), xlab = "Maximum Chl concentration [/]", ylab = "Recruitment/SSB [thousands/tonnes]")
text(hchlfallsmall_dat$Chlmax,hchlfallsmall_dat$recssb1, labels = hchlfallsmall_dat$Year, pos = 3, cex=0.5)
matlines(newchlfallsmall,confrecfallsmall[,1:3],lty=c(1,2,2), col=c(1,3,3),lwd=c(1.5,1.5,1.5))


#try with mean
hchlfallsmall <- subset(chldatenafallsmall, start=="Yes")
hchlfallsmall <- chldatenafallsmall %>% group_by(Year) %>% summarize(Chlmean = mean(Chl), recssb1=mean(recssb1))
hchlfallsmall <- na.omit(hchlfallsmall)
hchlfallsmall_dat <- data.frame(Year=hchlfallsmall$Year, Chlmean=hchlfallsmall$Chlmean, recssb1=hchlfallsmall$recssb1)
hchlfallsmall_dat <- subset(hchlfallsmall_dat, Year!="2022")

mchlhfallsmall <- lm(recssb1 ~ Chlmean, hchlfallsmall_dat)
summary(mchlhfallsmall)

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(hchlfallsmall_dat$Year, hchlfallsmall_dat$Chlmean, type='l', lwd=3, col = "green", xlab = "Year", ylab = "")
axis(2, col = "green", col.ticks = "green", col.axis = "green", col.lab="green")
mtext("Mean Chl concentration [/]", side = 2, line = 3, col = "green") 
par(new = TRUE) 
plot(hchlfallsmall_dat$Year, hchlfallsmall_dat$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext("Recruitment/SSB [thousands/tonnes]", side = 4, line = 4, col = "blue") 
legend("top", legend = c("Chl", "Rec/SSB"), col = c("green", "blue"), lty = 1, cex = 0.8)
#800 vs 450 aspect ratio

newchlfallsmall <- data.frame("Chlmean"=seq(min(hchlfallsmall_dat$Chlmean),max(hchlfallsmall_dat$Chlmean),length=1000))
confrecfallsmall <- predict(mchlhfallsmall,newchlfallsmall, interval="confidence")
confrecfallsmall

plot(hchlfallsmall_dat$Chlmean,hchlfallsmall_dat$recssb1, ylim=c(6,50), xlab = "Mean Chl concentration [/]", ylab = "Recruitment/SSB [thousands/tonnes]")
text(hchlfallsmall_dat$Chlmean,hchlfallsmall_dat$recssb1, labels = hchlfallsmall_dat$Year, pos = 3, cex=0.5)
matlines(newchlfallsmall,confrecfallsmall[,1:3],lty=c(1,2,2), col=c(1,3,3),lwd=c(1.5,1.5,1.5))









############### Bloom Phenology ################################

#Look at average year

#Use data
chldatena

#Create nweek row in data
pehn <- subset (chldatena, Year!=1997)
yphen <- levels(as.factor(pehn$Year))
pehn_names <- paste0("pheny_", 1:length(yphen))

for (i in 1:length(yphen)){
  yp <- subset (chldatena, Year==yphen[i])
  assign(pehn_names[i], yp)
}

#pehn_frames <- ls(pattern = "^pheny_\\d+$")

lengthsyphens <- data.frame("Year"=yphen, "nweeks"=rep(0,length(yphen)))

for (i in 1:length(yphen)) {
  n_obs <- nrow(get(pehn_names[i]))
  lengthsyphens$nweeks[i] <- n_obs
}

#Create nweek vector to attach
vec <- numeric()

for (i in lengthsyphens$nweeks) {
  vec <- c(vec, 1:i)
}

nrow(chldatena)
length(vec)

chldatena$nweek <- vec
chldatena

##Create plots

#Average for each week acros years
avweekchl <- chldatena %>% group_by(nweek) %>% summarize(mean_week_chl = mean(Chl))

nweeklab <- expression(paste("Week number"))

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(avweekchl$nweek, avweekchl$mean_week_chl, type='l', lwd=3, col = "green", xlab = nweeklab, ylab = chlaconc)
#820 x 420

#Create plots for every year separately
par(mfrow=c(5,5))
par(mar = c(2, 2, 2, 2))
plot(chldatena$nweek[chldatena$Year==yphen[1]], chldatena$Chl[chldatena$Year==yphen[1]], type='l', lwd=3, ylim=c(0.8,5.5), col = "green", xlab = "", ylab = "", main=yphen[1])

for (i in 2:24){
  plot(chldatena$nweek[chldatena$Year==yphen[i]], chldatena$Chl[chldatena$Year==yphen[i]], type='l', lwd=3, ylim=c(0.8,5.5), col = "green", xlab = "", ylab = "", main=yphen[i])
  
}
#900 x 600

#1125 x 750

#colfunc1 <- colorRampPalette(c("lightgrey", "darkred"))
colfunc <- topo.colors(24)
colfunc <- c(rep(1,5), rep(3,14), rep(2,5))

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(chldatena$nweek[chldatena$Year==yphen[1]], chldatena$Chl[chldatena$Year==yphen[1]], type='l', lwd=1, ylim=c(0.8,4), col = colfunc[1], xlab = "", ylab = "")

for (i in 2:24){
  points(chldatena$nweek[chldatena$Year==yphen[i]], chldatena$Chl[chldatena$Year==yphen[i]], type='l', lwd=1, ylim=c(0.8,4), col = colfunc[i], xlab = "", ylab = "")
  
}


#Plot with monthly points
monthsy <- chldatena %>% group_by(Year, Month) %>% summarize(mean_Chl_month = mean(Chl), month=mean(Month), sd_Chl_month = sd(Chl),
                                                             se_Chl_month = sd_Chl_month / sqrt(n()))
ggplot(monthsy, aes(x = month, y = mean_Chl_month)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_Chl_month - se_Chl_month,
                    ymax = mean_Chl_month + se_Chl_month),
                width = 0.2) +
  labs(x = "Month", y = "Chl") +
  theme_classic()

colfunc1 <- colorRampPalette(c("white", "darkgreen"))
colfunc <- colfunc1(28)[2:26]
colfunc <- topo.colors(24)

monthsy1 <- monthsy
monthsy1$Year <- as.factor(monthsy1$Year)

ggplot(monthsy, aes(x = month, y = mean_Chl_month, color = factor(Year), name="Year")) +
  geom_point(aes(size = se_Chl_month)) +
  geom_line()+
  scale_size_continuous(name = "Standard Error") +
  #geom_errorbar(aes(ymin = mean_Chl_month - se_Chl_month, ymax = mean_Chl_month + se_Chl_month), width = 0.2) +
  labs(x = "Month", y = expression(paste("Chl-a (", mu, "g ",L^-1,")"))) +
  scale_color_manual(values = colfunc, name="Year") +
  theme_minimal()


#Only fall
ggplot(monthsy, aes(x = month, y = mean_Chl_month, color = factor(Year))) +
  geom_point(aes(size = se_Chl_month)) +
  geom_line()+
  scale_size_continuous(name = "Standard Error") +
  #geom_errorbar(aes(ymin = mean_Chl_month - se_Chl_month, ymax = mean_Chl_month + se_Chl_month), width = 0.2) +
  labs(x = "Month", y = expression(paste("Chl-a (", mu, "g ",L^-1,")"))) +
  scale_color_manual(values = colfunc, name="Year") +
  scale_x_continuous(breaks = 1:12) +
  theme_minimal()
#850 x 420


#Model for indiviual month trends
impm <- chldatena %>% group_by(Month, Year) %>% summarize(mean_Chl_month = mean(Chl), month=mean(Month), sd_Chl_month = sd(Chl),
                                                             se_Chl_month = sd_Chl_month / sqrt(n()))
impmonths <- subset(impm, month>=8)

myonth <- lm(mean_Chl_month ~ Year + as.factor(month), impmonths)
summary(myonth)

ggplot(impmonths, aes(x = Year, y = mean_Chl_month, color = factor(Month))) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_Chl_month - se_Chl_month, ymax = mean_Chl_month + se_Chl_month), width = 0.2)+
  geom_line()+
  labs(x = "Year", y = expression(paste("Chl-a (", mu, "g ",L^-1,")"))) +
  scale_color_manual(values = rainbow(5), name="Year") +
  theme_minimal()

#Improved <- but uses wrong model
ggplot(impmonths, aes(x = Year, y = mean_Chl_month, color = factor(Month))) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_Chl_month - se_Chl_month, ymax = mean_Chl_month + se_Chl_month), width = 0.2)+
  geom_line()+
  geom_smooth(aes(group = month), method = "lm", se = TRUE) +
  labs(x = "Year", y = expression(paste("Chl-a (", mu, "g ",L^-1,")"))) +
  scale_color_manual(values = rainbow(5), name="Month") +
  theme_minimal()
#850 x 420


ggplot(impmonths, aes(x = Year, y = mean_Chl_month, group_by(factor(Month)),color = factor(Month))) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_Chl_month - se_Chl_month, ymax = mean_Chl_month + se_Chl_month), width = 0.2)+
  geom_line()+
  geom_smooth(method = "lm", formula='y ~ x') +
  labs(x = "Year", y = expression(paste("Chl-a (", mu, "g ",L^-1,")"))) +
  scale_color_manual(values = rainbow(5), name="Month") +
  theme_minimal()



#Remove 2022 from chldatena since we won't use 2022 when we have recssbs anyway
### MISTAKE ### pehn1 <- subset (chldatena, Month>=8)
pehn1 <- subset (chldatena, Month>=8)
nrow(pehn1)

#Create timing of max series for fall bloom
pehn1$nweek_max <- rep(0,nrow(pehn1))

for (i in 1:length(yphen)) {
  nwm <- which.max(pehn1$Chl[pehn1$Year==yphen[i]])
  pehn1$nweek_max[pehn1$Year==yphen[i]] <- pehn1$nweek[pehn1$Year==yphen[i]][nwm]
}

pehn1

#Plot time series of nweek max chl for fall
maxnweekyears <- pehn1 %>% group_by(Year) %>% summarize(nweek_maxchlfall = mean(nweek_max))

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(maxnweekyears$Year, maxnweekyears$nweek_maxchlfall, type='l', lwd=3, col = "darkred", xlab = "Year", ylab = "Week number of fall chlorophyll maximum")
#820 x 420

#With model preds
newyearnweekmax <- data.frame("Year"=seq(min(maxnweekyears$Year),max(maxnweekyears$Year),length=1000))
confnmax <- predict(m_nweekmax_year,newyearnweekmax, interval="confidence")
confnmax
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(maxnweekyears$Year, maxnweekyears$nweek_maxchlfall, type='l', lwd=3, col = "darkred", xlab = "Year", ylab = "Week number of fall chlorophyll maximum")
matlines(newyearnweekmax,confnmax[,1:3],lty=c(1,2,2), col=c(1,3,3),lwd=c(1.5,1.5,1.5))
#820 x 420

#Withn ggplot
# plot using ggplot
ggplot(maxnweekyears, aes(x = Year, y = nweek_maxchlfall)) +
  geom_line(color = "darkred", linewidth = 2) +
  geom_smooth(method = "lm", formula='y ~ x') +
  xlab("Year") +
  ylab("Week number of fall chlorophyll maximum") +
  theme_minimal()




#Check relationship with rec ssb
maxnweekyears$recssb1 <- rep(0, nrow(maxnweekyears))

for (l in 1:length(yphen)){
  maxnweekyears <- maxnweekyears %>%
    mutate(recssb1 = if_else(Year == yphen[l], her$recssb[her$Year==yphen[l]], recssb1))
  
}

m_nweek_maxchlfall <- lm(recssb1~nweek_maxchlfall, maxnweekyears)
summary(m_nweek_maxchlfall)

#Changes in phenology
m_nweekmax_year <- lm(nweek_maxchlfall~Year, maxnweekyears)
summary(m_nweekmax_year)

#Check without 2001!!!!!
maxnweekyears2001 <- subset(maxnweekyears, Year=!2001)
m_nweekmax_year2001 <- lm(nweek_maxchlfall~Year, maxnweekyears2001)
summary(m_nweekmax_year2001)

plot(maxnweekyears$nweek_maxchlfall, maxnweekyears$recssb1)
#no relationship

#Create anomalies, based on long term average




#Create derivateives as change between
#Question 1: i+1-i or i-(i-1)?

pehn2 <- chldatena

pehn2$delta_chl <- rep(0, nrow(pehn2))

for (i in 1:(nrow(pehn2)-1)){
  pehn2$delta_chl[i] <- pehn2$Chl[i+1]-pehn2$Chl[i]
}

pehn2$delta_chl[nrow(pehn2)] <- y2022chldatena$Chl[1]-pehn2$Chl[nrow(pehn2)]

#Fill in last value using 2022
chldate1 <- cbind(mat_data$CHLWavg,mat_data$DATEW)
chldate1 <- data.frame(Chl=chldate1[,1], Year=chldate1[,2], Day = chldate1[,4],
                      Month = chldate1[,3])
y2022chldatena <- subset(chldate1, Year==2022)


#Plot time series of derivative
#can do that with if statements!!

##Without labs
par(mfrow=c(5,5))
par(mar = c(2, 2, 2, 2))
plot(pehn2$nweek[pehn2$Year==yphen[1]], pehn2$delta_chl[pehn2$Year==yphen[1]], type='l', lwd=2, ylim=c(-2,2), col = "darkorange", xlab = "", ylab = "", main=yphen[1])

for (i in 2:24){
  plot(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$delta_chl[pehn2$Year==yphen[i]], type='l', lwd=2, ylim=c(-2,2), col = "darkorange", xlab = "", ylab = "", main=yphen[i])
  
}

#for entire add
plot(changeweekyears$nweek, changeweekyears$mean_delta_chl, type='l', lwd=2, col = "darkorange", xlab = "", ylab = "", main="Mean 1998-2021")

####plot rafa
colfunc <- topo.colors(24)
colfunc <- c(rep(1,5), rep(3,14), rep(2,5))

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(pehn2$nweek[chldatena$Year==yphen[1]], pehn2$Chl_GAM[pehn2$Year==yphen[1]], type='l', lwd=1, ylim=c(1,3), col = colfunc[1], xlab = "", ylab = "")

for (i in 2:24){
  points(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$Chl_GAM[pehn2$Year==yphen[i]], type='l', lwd=1, ylim=c(1,3), col = colfunc[i], xlab = "", ylab = "")
  
}



#1125 x 750

#Make plot of average change across years
changeweekyears <- pehn2 %>% group_by(nweek) %>% summarize(mean_delta_chl = mean(delta_chl))

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(changeweekyears$nweek, changeweekyears$mean_delta_chl, type='l', lwd=3, col = "darkorange", xlab = "Year", ylab = "Mean change in chlorophyll concentration")
#820 x 420


#Find start of the
falldeltas <- subset(pehn2, Month>=8)

maxchangeyears <- falldeltas %>% group_by(Year) %>% summarize(max_delta_chl = max(delta_chl))

#Changes in phenology
m_maxchange_year <- lm(max_delta_chl~Year, maxchangeyears)
summary(m_maxchange_year)

#Check relationship with rec ssb
maxchangeyears$recssb1 <- rep(0, nrow(maxchangeyears))

for (l in 1:length(yphen)){
  maxchangeyears <- maxchangeyears %>%
    mutate(recssb1 = if_else(Year == yphen[l], her$recssb[her$Year==yphen[l]], recssb1))
  
}

#Again no relationship with herring rec
m_maxchange_rec <- lm(recssb1~max_delta_chl, maxchangeyears)
summary(m_maxchange_rec)
plot(maxchangeyears$max_delta_chl, maxchangeyears$recssb1)
#no relationship

#Time series plots
maxchangelab <- expression(paste("Maximum change of fall Chl-a ( ",mu," g ",L^-1,week^-1,")"))
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(maxchangeyears$Year, maxchangeyears$max_delta_chl, type='l', lwd=3, col = "darkorange", xlab = "Year", ylab = "", main="Entire North Sea")
axis(2, col = "darkorange", col.ticks = "darkorange", col.axis = "darkorange", col.lab="darkorange")
mtext(maxchangelab, side = 2, line = 3, col = "darkorange") 
par(new = TRUE) 
plot(maxchangeyears$Year, maxchangeyears$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext(meanrssblab, side = 4, line = 4, col = "blue") 
legend("top", legend = c("Change in Chl", "ln(R/SSB)"), col = c("darkorange", "blue"), lty = 1, cex = 0.8, bty="n")
#800 vs 450 aspect ratio


####### Week number of start of fall bloom ####

#Create timing of max series for fall bloom
falldeltas$nweek_maxchange <- rep(0,nrow(falldeltas))

for (i in 1:length(yphen)) {
  nwm <- which.max(falldeltas$delta_chl[falldeltas$Year==yphen[i]])
  falldeltas$nweek_maxchange[falldeltas$Year==yphen[i]] <- falldeltas$nweek[falldeltas$Year==yphen[i]][nwm]
}

falldeltas

#Plot time series of nweek max change chl
maxchangenweekyears <- falldeltas %>% group_by(Year) %>% summarize(nweek_maxchange1 = mean(nweek_maxchange))

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(maxchangenweekyears$Year, maxchangenweekyears$nweek_maxchange1, type='l', lwd=3, col = "darkgoldenrod4", xlab = "Year", ylab = "Week number of max chl change (fall)")
#820 x 420

#Check relationship with rec ssb
maxchangenweekyears$recssb1 <- rep(0, nrow(maxchangenweekyears))

for (l in 1:length(yphen)){
  maxchangenweekyears <- maxchangenweekyears %>%
    mutate(recssb1 = if_else(Year == yphen[l], her$recssb[her$Year==yphen[l]], recssb1))
  
}

#No relationship with recssb
m_nweekmaxchange_rec <- lm(recssb1~nweek_maxchange1, maxchangenweekyears)
summary(m_nweekmaxchange_rec)

#Check at temporal changes in phenology
m_nweekmaxchange_year <- lm(nweek_maxchange1~Year, maxchangenweekyears)
summary(m_nweekmaxchange_year)
#No change in phenology








####Try fitting a spline with a GAM to get derivatives

#Take one year as example
testy <- subset(pehn2, Year==2009)

mtesty <- gam(Chl ~ s(nweek), data=testy)
summary(mtesty)
pred_testy <- data.frame(nweek=testy$nweek)
pred_testy$Chl <- predict(mtesty, newdata = pred_testy)

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(testy$nweek, testy$Chl, type='l', lwd=3, col = "darkgreen", xlab = nweeklab, ylab = chlaconc, main="2009")
points(pred_testy$nweek, pred_testy$Chl, type='l', lwd=3, lty=3, col = "darkgreen", xlab = "", ylab = "")
legend("top", legend = c("Observation", "GAM spline"), col = c("darkgreen", "darkgreen"), lty = c(3,1), cex = 0.8)

#created modelled Chl column
vets <- numeric()

for (i in 1:length(yphen)) {
  yearm <- subset(pehn2, Year==yphen[i])
  myearm <- gam(Chl ~ s(nweek), data=yearm)
  pred_yearm <- data.frame(nweek=yearm$nweek)
  pred_yearm$Chl_GAM <- predict(myearm, newdata = pred_yearm)
  vets <- c(vets, pred_yearm$Chl_GAM)
}

pehn2$Chl_GAM <- vets

#Recreate plots of all years
par(mfrow=c(5,5))
par(mar = c(2, 2, 2, 2))
plot(pehn2$nweek[pehn2$Year==yphen[1]], pehn2$Chl[pehn2$Year==yphen[1]], type='l', lwd=2, ylim=c(0.8,4), col = "green", xlab = "", ylab = "", main=yphen[1])
points(pehn2$nweek[pehn2$Year==yphen[1]], pehn2$Chl_GAM[pehn2$Year==yphen[1]], type='l', lwd=2, lty=3, col = "darkgreen", xlab = "", ylab = "")
points(pehn2$nweek[pehn2$Year==yphen[1]], pehn2$mediansGAM[pehn2$Year==yphen[1]], type='l', lwd=1, lty=3, col="grey10")
points(pehn2$nweek[pehn2$Year==yphen[1]], pehn2$medians[pehn2$Year==yphen[1]], type='l', lwd=1, col="grey65")
legend("topright", legend = c("GAM spline", "Observation"), col = c("darkgreen", "green"), lty = c(3,1), lwd=2, cex = 0.7, bty = "n")

plot(pehn2$nweek[pehn2$Year==yphen[2]], pehn2$Chl[pehn2$Year==yphen[2]], type='l', lwd=2, ylim=c(0.8,4), col = "green", xlab = "", ylab = "", main=yphen[1])
points(pehn2$nweek[pehn2$Year==yphen[2]], pehn2$Chl_GAM[pehn2$Year==yphen[2]], type='l', lwd=2, lty=3, col = "darkgreen", xlab = "", ylab = "")
points(pehn2$nweek[pehn2$Year==yphen[2]], pehn2$mediansGAM[pehn2$Year==yphen[2]], type='l', lwd=1, lty=3, col="grey10")
points(pehn2$nweek[pehn2$Year==yphen[2]], pehn2$medians[pehn2$Year==yphen[2]], type='l', lwd=1, col="grey65")
legend("topleft", legend = c("5% above median", "5% above median (GAM)"), col = c("grey65","grey10"), lty = c(1,3), lwd=2, cex = 0.7, bty = "n")


for (i in 3:24){
  plot(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$Chl[pehn2$Year==yphen[i]], type='l', lwd=2, ylim=c(0.8,4), col = "green", xlab = "", ylab = "", main=yphen[i])
  points(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$Chl_GAM[pehn2$Year==yphen[i]], type='l', lwd=2, lty=3, col = "darkgreen", xlab = "", ylab = "")
  points(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$mediansGAM[pehn2$Year==yphen[i]], type='l', lwd=1, lty=3, col="grey10")
  points(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$medians[pehn2$Year==yphen[i]], type='l', lwd=1, col="grey65")
  
}

#Average for each week acros years
#avweekchl <- chldatena %>% group_by(nweek) %>% summarize(mean_week_chl = mean(Chl))
m_avweekchl <- gam(mean_week_chl ~ s(nweek), data=avweekchl)
pred_avweekchl <- data.frame(nweek=avweekchl$nweek)
pred_avweekchl$mean_week_chl_GAM <- predict(m_avweekchl, newdata = pred_avweekchl)
avweekchl$mean_week_chl_GAM <- pred_avweekchl$mean_week_chl_GAM


plot(avweekchl$nweek, avweekchl$mean_week_chl, type='l', lwd=2, col = "green", xlab = "", ylab = "", main="Mean 1998-2021")
points(avweekchl$nweek, avweekchl$mean_week_chl_GAM, type='l', lwd=2, lty=3, col = "darkgreen", xlab = "", ylab = "")
#1125 x 750

#points(fer2019GAM$nweek, median2019GAM, type='l', lwd=2)

#Calculate the change in chl based on GAM splines
#Without doing that for the last week of each year
pehn3 <- pehn2

pehn3$delta_chl_GAM <- rep(NA, nrow(pehn2))

for (j in 1:length(yphen)) {
  yearxd <- subset(pehn2, Year==yphen[j])
  
  for (i in 1:(nrow(yearxd)-1)){
    pehn3$delta_chl_GAM[pehn3$Year==yphen[j]][i]<- yearxd$Chl_GAM[yearxd$Year==yphen[j]][i+1]-yearxd$Chl_GAM[yearxd$Year==yphen[j]][i]
  }
  
}

pehn2$delta_chl_GAM <- pehn3$delta_chl_GAM

#Plots of these changes

###############combine plots
par(mfrow=c(5,5))
par(mar = c(2, 2, 2, 2))
plot(pehn2$nweek[pehn2$Year==yphen[1]], pehn2$delta_chl_GAM[pehn2$Year==yphen[1]], type='l', lty=3, lwd=2, ylim=c(-0.4,0.4), col = "darkorange3", xlab = "", ylab = "", main=yphen[1])
plot(pehn2$nweek[pehn2$Year==yphen[1]], pehn2$delta_chl[pehn2$Year==yphen[1]], type='l', lty=1, lwd=2, ylim=c(-0.4,0.4), col = "orange", xlab = "", ylab = "", main=yphen[1])

for (i in 2:24){
  plot(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$delta_chl_GAM[pehn2$Year==yphen[i]], type='l', lwd=2, lty=3, ylim=c(-0.4,0.4), col = "darkorange3", xlab = "", ylab = "", main=yphen[i])
  points(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$delta_chl[pehn2$Year==yphen[i]], type='l', lty=1, lwd=2, ylim=c(-0.4,0.4), col = "orange", xlab = "", ylab = "", main=yphen[i])
}
changeGAM <- pehn2 %>% group_by(nweek) %>% summarize(mean_delta_chl_GAM = mean(delta_chl_GAM))

plot(changeGAM$nweek, changeGAM$mean_delta_chl_GAM, type='l', lwd=2,lty=3, col = "darkorange3", xlab = "", ylab = "", main="Mean 1998-2021")

#Take fall max change = start of bloom and check against recssb
##Without labs
par(mfrow=c(5,5))
par(mar = c(2, 2, 2, 2))
plot(pehn2$nweek[pehn2$Year==yphen[1]], pehn2$delta_chl_GAM[pehn2$Year==yphen[1]], type='l', lwd=2, ylim=c(-0.4,0.4), col = "darkorange", xlab = "", ylab = "", main=yphen[1])

for (i in 2:24){
  plot(pehn2$nweek[pehn2$Year==yphen[i]], pehn2$delta_chl_GAM[pehn2$Year==yphen[i]], type='l', lwd=2, ylim=c(-0.4,0.4), col = "darkorange", xlab = "", ylab = "", main=yphen[i])
  
}

#average change across years
changeGAM <- pehn2 %>% group_by(nweek) %>% summarize(mean_delta_chl_GAM = mean(delta_chl_GAM))

plot(changeGAM$nweek, changeGAM$mean_delta_chl_GAM, type='l', lwd=2, col = "darkorange", xlab = "", ylab = "", main="Mean 1998-2021")

#timing of bloom: max change
#Find start of the
fallGAMdeltas <- subset(pehn2, Month>=8)

maxGAMdelta <- fallGAMdeltas %>% group_by(Year) %>% summarize(max_GAM_delta = max(delta_chl_GAM, na.rm = TRUE))

#Changes in phenology
m_maxGAMdelta <- lm(max_GAM_delta~Year, maxGAMdelta)
summary(m_maxGAMdelta)

#Changes in phenology without 2002
maxGAMdelta2002 <- subset(maxGAMdelta, Year=!2002)
m_maxGAMdelta2002 <- lm(max_GAM_delta~Year, maxGAMdelta2002)
summary(m_maxGAMdelta2002)

#Check relationship with rec ssb
maxGAMdelta$recssb1 <- rep(0, nrow(maxGAMdelta))

for (l in 1:length(yphen)){
  maxGAMdelta <- maxGAMdelta %>%
    mutate(recssb1 = if_else(Year == yphen[l], her$recssb[her$Year==yphen[l]], recssb1))
  
}

#REC SSB
m_maxGAMdelta_rec <- lm(recssb1 ~ max_GAM_delta, maxGAMdelta)
summary(m_maxGAMdelta_rec)


####### Week number of start of fall bloom for GAM spline values####

#Create timing of max series for fall bloom
fallGAMdeltas$nweek_maxchange_GAM <- rep(0,nrow(fallGAMdeltas))

for (i in 1:length(yphen)) {
  nwm1 <- which.max(fallGAMdeltas$delta_chl_GAM[fallGAMdeltas$Year==yphen[i]])
  fallGAMdeltas$nweek_maxchange_GAM[fallGAMdeltas$Year==yphen[i]] <- fallGAMdeltas$nweek[fallGAMdeltas$Year==yphen[i]][nwm1]
}

fallGAMdeltas

#Plot time series of nweek max change chl
nweekGAMmax <- fallGAMdeltas %>% group_by(Year) %>% summarize(nweek_maxchange_GAM1 = mean(nweek_maxchange_GAM))

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(nweekGAMmax$Year, nweekGAMmax$nweek_maxchange_GAM1, type='l', lwd=3, col = "darkgoldenrod2", xlab = "Year", ylab = "Week number of max chl change in the fall (GAM)")
#820 x 420

#With model preds
new_GAMyearnweekm <- data.frame("Year"=seq(min(nweekGAMmax$Year),max(nweekGAMmax$Year),length=1000))
confnmaxGAM <- predict(m_nmaxGAMdelta,new_GAMyearnweekm, interval="confidence")
confnmaxGAM
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(nweekGAMmax$Year, nweekGAMmax$nweek_maxchange_GAM1, type='l', lwd=3, col = "purple3", xlab = "Year", ylab = "Week number of fall max Chla change (GAM)")
matlines(new_GAMyearnweekm,confnmaxGAM[,1:3],lty=c(1,2,2), col=c(1,3,3),lwd=c(1.5,1.5,1.5))

##Ggplot
ggplot(nweekGAMmax, aes(x = Year, y = nweek_maxchange_GAM1)) +
  geom_line(color = "purple3", linewidth = 2) +
  geom_smooth(method = "lm", formula='y ~ x') +
  xlab("Year") +
  ylab("Week number of fall max Chla change (GAM)") +
  theme_minimal()




#Check relationship with rec ssb
nweekGAMmax$recssb1 <- rep(0, nrow(nweekGAMmax))

for (l in 1:length(yphen)){
  nweekGAMmax <- nweekGAMmax %>%
    mutate(recssb1 = if_else(Year == yphen[l], her$recssb[her$Year==yphen[l]], recssb1))
  
}

#Changes in phenology -> week significance, p value is 0.03, but we use Bonferroni correction
m_nmaxGAMdelta <- lm(nweek_maxchange_GAM1~Year, nweekGAMmax)
summary(m_nmaxGAMdelta)

#REC SSB -> no relationship
m_nmaxGAMdelta_rec <- lm(recssb1 ~ nweek_maxchange_GAM1, nweekGAMmax)
summary(m_nmaxGAMdelta_rec)



#Calculate start of bloom like ferreira et al, i.e. > 5%>median

pehn2

ferreira <- data.frame("Year"=1998:2021, "week_start"=rep(0,nrow(nweekGAMmax)))

for (i in 1:length(yphen)){
  fersta <- subset(pehn2, Year==yphen[i])
  fersta <- subset(fersta, Month>=9)
  row_start <- which(fersta$Chl > 1.05*fersta$medians[fersta$Year==yphen[i]][1])[1]
  ferreira$week_start[ferreira$Year==yphen[i]] <- fersta$nweek[row_start]
  
}

#Test
fer2019 <- subset(pehn2, Year==2019)

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(fer2019$nweek, fer2019$Chl, type='l', col='green', lwd=2, main="2019", xlab=nweeklab, ylab=chlaconc)
points(fer2019$nweek, 1.05*fer2019$medians, type='l', lwd=2)
legend("top", legend = c("5% above annual median", "Observed Chla"), col = c("black", "green"), lty = c(1,1), lwd=2, cex = 0.7, bty = "n")


plot(ferreira$Year,ferreira$week_start, type='l', lwd=2, col="red4", xlab="Year", ylab="Start of the fall bloom (week number)")

###Time series plot

plot(ferreira$Year,ferreira$week_start, type='l', lwd=2, col="red4", xlab="Year", ylab="Start of the fall bloom (week number)")
#Time series plots
nstart <- expression("Start of the fall bloom (week number)")
par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(ferreira$Year, ferreira$week_start, type='l', lwd=3, col = "red4", xlab = "Year", ylab = "", main="Entire North Sea")
axis(2, col = "red4", col.ticks = "red4", col.axis = "red4", col.lab="red4")
mtext(nstart, side = 2, line = 3, col = "red4") 
par(new = TRUE) 
plot(maxchangeyears$Year, maxchangeyears$recssb1, type='l', lwd=3, col = "blue", xlab = "", ylab = "", xaxt="n", yaxt="n")
axis(side=4, col = "blue", col.ticks = "blue", col.axis = "blue", col.lab="blue")
mtext(meanrssblab, side = 4, line = 4, col = "blue") 
legend("top", legend = c("Start of bloom", "ln(R/SSB)"), col = c("red4", "blue"), lty = 1, cex = 0.8, bty="n")
#800 vs 450 aspect ratio



mferr <- lm(week_start ~ Year, ferreira)
summary(mferr)

#Add herring rec
ferreira$recssb1 <- rep(0, nrow(ferreira))

for (l in 1:length(yphen)){
  ferreira <- ferreira %>%
    mutate(recssb1 = if_else(Year == yphen[l], her$recssb[her$Year==yphen[l]], recssb1))
  
}

mstartfer <- lm(recssb1 ~ week_start, ferreira)
summary(mstartfer)

###############ferreira start GAM
pehn2

texty1998 <- subset(pehn2, Year==1998)
median(texty1998$Chl_GAM)

pehn2$mediansGAM <- rep(0, length(pehn2$Chl_GAM))
mediansyGAM <- pehn2 %>% group_by(Year) %>% summarize(median = median(Chl_GAM))

yrsGAM <- levels(pehn2$YearF)

for (i in 1:length(yrsGAM)){
  pehn2 <- pehn2 %>%
    mutate( mediansGAM = if_else(YearF == yrsGAM[i], mediansyGAM$median[mediansyGAM$Year==yrsGAM[i]], mediansGAM))
  
}


ferreiraGAM <- data.frame("Year"=1998:2021, "week_start_GAM"=rep(0,nrow(nweekGAMmax)))

for (i in 1:length(yphen)){
  ferstaGAM <- subset(pehn2, Year==yphen[i])
  ferstaGAM <- subset(ferstaGAM, Month>=9)
  row_startGAM <- which(ferstaGAM$Chl_GAM > 1.05*ferstaGAM$mediansGAM[ferstaGAM$Year==yphen[i]][1])[1]
  ferreiraGAM$week_start_GAM[ferreiraGAM$Year==yphen[i]] <- ferstaGAM$nweek[row_start]
  
}

plot(ferreiraGAM$Year,ferreiraGAM$week_start_GAM, type='l', lwd=2, col="orange3", xlab="Year", ylab="Start of the fall bloom (week number, GAM)", ylim=c(46,49))


#test 2019 GAM
#Test
fer2019GAM <- subset(pehn2, Year==2019)
xxx <- rep(1, length(fer2019GAM$Chl_GAM))
median2019GAM <- 1.05*median(fer2019GAM$Chl_GAM)*xxx

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(fer2019GAM$nweek, fer2019GAM$Chl_GAM, type='l', col='darkgreen', lty=3, lwd=2, main="2019", xlab=nweeklab, ylab=chlaconc)
points(fer2019GAM$nweek, median2019GAM, type='l', lwd=2)
legend("top", legend = c("5% above annual median", "Observed Chla"), col = c("black", "darkgreen"), lty = c(1,3), lwd=2, cex = 0.7, bty = "n")

fer2004GAM <- subset(pehn2, Year==2004)
xxx <- rep(1, length(fer2004GAM$Chl_GAM))
median2004GAM <- 1.05*median(fer2004GAM$Chl_GAM)*xxx

par(mfrow=c(1,1))
par(mar = c(4, 4, 4, 4) + 1.5) 
plot(fer2004GAM$nweek, fer2004GAM$Chl_GAM, type='l', col='darkgreen', lty=3, lwd=2, main="2004", xlab=nweeklab, ylab=chlaconc)
points(fer2004GAM$nweek, median2004GAM, type='l', lwd=2)
legend("topright", legend = c("5% above annual median", "Observed Chla"), col = c("black", "darkgreen"), lty = c(1,3), lwd=2, cex = 0.7, bty = "n")

