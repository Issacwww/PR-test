library(sp)
library(sf)
library(ggmap)
library(rgdal)
library(rgeos)
library(dplyr)
library(tidyr)
library(tmap)
library(raster)
library(jpeg)
library(ggplot2)
library(ggspatial)
library(tools)
library(OIdata)
library(RColorBrewer)
library(classInt)
library(maps)
library(mapproj)
library(GISTools)
library(farver)

library(rnaturalearth)
library(rnaturalearthdata)
library(usmap)
library(mapdata)

library(spatstat)

#part A
setwd("~/Courses/SDA/Midterm/PACoals")
# setwd("~/Courses/SDA/Midterm/PALocs")


input <- readOGR(dsn= getwd(), layer = "PACoals")
#input <- readOGR(dsn= getwd(), layer = "PALocs")

InputData <- spTransform(input,CRS("+proj=lcc +lat_1=40.88333333333333 +lat_2=41.95 +lat_0=40.16666666666666 +lon_0=-77.75 +x_0=600000 +y_0=0
+ellps=GRS80 +units=us-ft +no_defs "))
boundary = extent(InputData)
Xmin = boundary@xmin
Xmax = boundary@xmax
Ymin = boundary@ymin
Ymax = boundary@ymax

Xlen = (Xmax - Xmin) 
Ylen = (Ymax - Ymin)
area = Xlen * Ylen
 


quardrants_size = 6000
quardrants_num =1500
#while(quardrants_size > 0){
  quad_area = area / (quardrants_size)
  dim_x = sqrt(quad_area)
  dim_y = dim_x
  print(dim_x)
  plot(InputData,main="Random Quadrat Sampling of PACoals")
  #plot(InputData,main="Random Quadrat Sampling of PALocs")
  
  quads <- matrix(1:quardrants_num, quardrants_num, 1)
  i = 1
  while(i <= quardrants_num){
    centX <- runif(1, Xmin, Xmax)
    centY <- runif(1, Ymin, Ymax)
    xleft = centX - dim_x / 2
    ybutton = centY - dim_y / 2
    xright = centX + dim_x / 2
    ytop = centY + dim_y / 2
    k <- nrow(InputData@data[InputData@coords[,1]< xright &
                               InputData@coords[,1] > xleft &
                               InputData@coords[,2] < ytop &
                               InputData@coords[,2] > ybutton,])
    quads[i] = k
    rect(xleft, ybutton, xright, ytop, border = "red")
    i = i + 1
  }
  
  
  
  mu = length(InputData) / quardrants_num
  #calculating the variance
  #calculating the variance
  #Calculating the Variance mean ratio for regular sampling method
  
  K=as.numeric(names(table(c(quads)))) #Number of events K
  X=c(table(c(quads))) #Number of quadrats X
  table=data.frame(K,X) #Create data frame for table with K and X 
  rownames(table)=NULL#Setting the first row name to NULL for easy data collection
  table$K_mu=table$K-mu #K-mu
  table$K_mu2=table$K_mu^2 #(K-mu)^2
  table$XK_mu2=table$X*table$K_mu2 #X(K-mu)^2
  total=sum(table$XK_mu2) #total 
  variance =total/(quardrants_num - 1)
  VMR=variance/mu 
  print(VMR)
  
  
  quardrants_size = quardrants_size - 100
#}


write.csv(table,paste("~/Courses/SDA/Midterm/PACoals.csv"))
#write.csv(table,paste("~/Courses/SDA/Midterm/PALocs.csv"))

######## part b

dis <- nndist(InputData@coords, k=1)
mx <- max(dis) + 100

x <- seq(from = 0, to = mx, by = 500)
y <- c()
for (i in x) {
  nr <- (sum(dis < i)/nrow(InputData))
  y <- c(y, nr)
}

plot(x, y, type="l",
     xlim=c(0,25000), #(xlim=c(0,600) for F function)
     ylim=c(0,1),
     lwd=1,
     xlab="distance",
     ylab="G",
     main="G")

Gdata <- as(InputData, "ppp")
plot(Gest(Gdata))


random_points <- spsample(InputData,500,"random")
plot(InputData)
plot(random_points,add=TRUE,col="red")

count = 0

x2 <- seq(from = 0, to = mx, by = 500)
y2 <- c()

for(i in x2){
  count = 0
  for (j in seq(1,500)){
    point = random_points[j,]
    IDs = which(spDists(point, InputData) <= i)
    if(length(IDs) > 0)
      count = count + 1
  }
  y2 <- c(y2, count/500)
}

plot(x2, y2, type="l",
     xlim=c(0,30000), #(xlim=c(0,600) for F function)
     ylim=c(0,1),
     lwd=1,
     xlab="distance",
     ylab="F",
     main="F")

plot(Fest(Gdata))


