library(adehabitatLT)
library(adehabitatHR)
library(tidyverse)
library(move)
library(circular)
library(sp)
library(maptools)
library(raster)
library(rgdal)
library(sf)
library(stars)
library(terra)
library(dplyr)
library(Rmpfr)
library(leaflet)
library(lubridate)  


deer_data1 <- read.csv("C:/R/R-4.2.2/DEER/T4/15D/T4-15D2.csv",header=TRUE,sep=",")
# Combine longitude, latitude and time data into a data.frame
deer_move1 <- move(x=deer_data1$location.long, y=deer_data1$location.lat,
                  time=as.POSIXct(deer_data1$t, format="%Y/%m/%d %H:%M", 
                  ),
                  data=deer_data1, proj=CRS("+proj=longlat +ellps=WGS84"),
                  animal=deer_data1$id, sensor=deer_data1$sensor)


deer_move1 <- spTransform(deer_move1, CRS("+proj=utm +zone=52 +datum=WGS84"))
deer_dbbmm1 <- brownian.bridge.dyn(object=deer_move1, location.error=20, window.size=15,
                                margin=5, dimSize=100,time.step=60)

deer_data2 <- read.csv("C:/R/R-4.2.2/DEER/T4/15D/T4-15D2.csv",header=TRUE,sep=",")
# Combine longitude, latitude and time data into a data.frame
deer_move2 <- move(x=deer_data2$location.long, y=deer_data2$location.lat,
                   time=as.POSIXct(deer_data2$t, format="%Y/%m/%d %H:%M", 
                   ),
                   data=deer_data2, proj=CRS("+proj=longlat +ellps=WGS84"),
                   animal=deer_data2$id, sensor=deer_data2$sensor)


deer_move2 <- spTransform(deer_move2, CRS("+proj=utm +zone=52 +datum=WGS84"))
deer_dbbmm2 <- brownian.bridge.dyn(object=deer_move2, location.error=20, window.size=15,
                                   margin=5, dimSize=100,time.step=60)


deer_data3 <- read.csv2("T1data_4new.csv",
                        header = TRUE,
                        sep = ",")
deer_move3 <- move(x = as.numeric(deer_data3$location_long),
                   y = as.numeric(deer_data3$location_lat),
                   time=as.POSIXct(deer_data3$t, format="%Y/%m/%d %H:%M"), 
                   data=deer_data3, 
                   proj=CRS("+proj=longlat +ellps=WGS84"),
                   animal=deer_data3$id,
                   sensor=deer_data3$sensor)
table(deer_data3$id)  # 查看每个季节有多少条记录

deer_move3 <- spTransform(deer_move3,
                          CRS("+proj=utm +zone=52 +datum=WGS84"))
#ext <- as.numeric(c(min(deer_move3$location_long), 
#        max(deer_move3$location_long), 
#        min(deer_move3$location_lat), 
#        max(deer_move3$location_lat)))
#str(ext)
deer_dbbmm3 <- brownian.bridge.dyn(object=deer_move3, 
                                   location.error=20, 
                                   window.size=31,
                                   margin=9, 
                                   dimSize=100,
                                   time.step=60,
                                   ext = 0.3
                                  )

deer_ltraj <- as(deer_move3,"ltraj")
xccc#分离stack
split(deer_dbbmm3)

class(deer_dbbmm1)
summary(deer_dbbmm3)
str(deer_move3)


show(deer_dbbmm3)
summary(deer_dbbmm3)
plot(deer_move3)
plot(deer_dbbmm3)
plot(deer_dbbmm3@layers[[3]])


raster(deer_dbbmm3)



 contour(deer_dbbmm3, levels=0.95, add=TRUE)
 contour(deer_dbbmm3, levels=0.5, add=TRUE)
 contour(deer_dbbmm3[[3]], levels=0.95, add=TRUE)
contour(deer_dbbmm1, levels=c(.5), add=TRUE)
contour(deer_dbbmm2, levels=c(.5), add=TRUE)

#leaflet包，互动地图，坐标得是经纬度
#导入tif、shp等地图数据
#需要经纬度坐标，不能用：deer_move3 <- spTransform(deer_move3, CRS("+proj=utm +zone=52 +datum=WGS84"))
map_raster <- raster("E:/ARCgisdata/TPG2023/ASTGTMV003_N48E130_dem.tif")
bounds <- as.vector(bbox(extent(cont1_shape)))
base_map <- leaflet() %>% 
  fitBounds(bounds[1], bounds[2], bounds[3], bounds[4]) %>% addTiles("https://webst01.is.autonavi.com/appmaptile?style=6&x={x}&y={y}&z={z}")
map2 <- base_map %>% 
  addPolylines(data =  as(cont1_shape,'SpatialLines'), color ="grey") %>% 
  addCircles(data =cont1_shape,fillOpacity=0.3, opacity = 0.5, color="blue")
map3 <- map2 %>% 
  addLegend(position= "topright", colors=c("grey","blue"), labels=c("lines","points") ,opacity = 0.7, title = "deer_move3")
map4 <- map3 %>% 
  ddScaleBar(position="bottomleft",options=scaleBarOptions(maxWidth = 100, metric = TRUE, imperial = F, updateWhenIdle = TRUE))
map4

## to optimize the calculation, the cells outside of the 99.99% UD contour 
# are removed by setting them to zero.
#values(deer_dbbmm3)[values(getVolumeUD(deer_dbbmm3))>.999999]<-0
## transform each layer to a probability surface (i.e. sum of their values is 1)
#stk<-(deer_dbbmm3/cellStats(deer_dbbmm3,sum))
#str(stk)


# emd ---------------------------------------------------------------------
#查看少了哪些X7D
levels(deer_dbbmm3@DBMvar@trackIdUnUsedRecords)

# 创建一个空的向量来存储结果
# 创建一个空的向量来存储结果

emd1_values <- numeric(0)
#结果存储到向量中
for (i in 1:(length(deer_dbbmm3@layers) - 1)) {                         
  # emd的循环编号仅代表stack的数量，如中间缺失几个，这些缺失值的序号将继承给接下来的值
  emd1_result <- emd(deer_dbbmm3[[i]], deer_dbbmm3[[i+1]])  # 假设stk是你的数据列表
  emd1_values <- c(emd1_values, emd1_result)
}
# 循环运行emd函数并将
# 将结果存储到数据框中
emd1_values <- data.frame(emd = emd1_values)
# 将结果保存为CSV文件
write.csv(emd1_values, file = "T1/emd_STACK.csv", row.names = FALSE)




emd(stk, threshold=10000)
emd(deer_dbbmm3[[1]],deer_dbbmm3[[2]])
emd(deer_dbbmm1,deer_dbbmm2)
emd(deer_dbbmm3)


# 位移中心到补饲点距离 --------------------------------------------------------------

# 按组名计算经纬度坐标的均值
str(deer_data3)
#center <- aggregate(cbind(as.numeric(location_long),
#                          as.numeric(location_lat)) ~ X7D, deer_data3, mean)
center <- aggregate(list(location_long = as.numeric(deer_data3$location_long), 
                         location_lat = as.numeric(deer_data3$location_lat)), 
                    by = list(X7D = deer_data3$X7D), 
                    FUN = mean)



# 将中心坐标的经纬度保留到小数点后6位
center$location_long <- round(center$location_long, 6)
center$location_lat <- round(center$location_lat, 6)

# 给定的坐标点
given_point <- c(130.7285, 48.17981)

# 将经纬度坐标转换为WGS84坐标系
center_wgs84 <- SpatialPointsDataFrame(coords = center[, c("location_long", "location_lat")],
                                       data = center,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))
# 计算中心坐标到给定坐标点的距离
center$distance <- distGeo(center_wgs84, given_point)

# 对组名进行排序
center <- center %>% 
  arrange(as.numeric(sub(".*-([0-9]+)", "\\1", X7D)))

write.csv(center, "T1/T1_dist.csv")

# 用amt计算质心 ----------------------------------------------------------------


library(amt)
deer_track <- make_track(deer_data3, 
           .x = location_long,
           .y = location_lat,
           .t = t)
centroid(deer_move3)

# 导出shp -------------------------------------------------------------------


#画轨迹线
par(mar = c(5, 5, 1, 1))
plot(deer_move3, xlab="Longitude", ylab="Latitude", type="l", pch=16, lwd=0.5)
points(deer_move3, pch=20, cex=0.5)
#画轮廓线
cont1 <- raster2contour(deer_dbbmm3, level=c(0.95))
cont2 <- raster2contour(deer_dbbmm3@layers[[1]], level=c(0.95))
cont3 <- raster2contour(deer_dbbmm1, level=c(0.5))
cont4 <- raster2contour(deer_dbbmm2, level=c(0.5))
plot(cont1,col = c("black"))
plot(cont2,col = c("black"))
plot(cont3,col = c("black"))
plot(cont4,col = c("black"))
#导出轮廓线到arcgis
cont1_shape <- as(cont1, "SpatialLinesDataFrame")  # 将轮廓线数据转换为SpatialLinesDataFrame对象
cont2_shape <- as(cont2, "SpatialLinesDataFrame")
cont3_shape <- as(cont3, "SpatialLinesDataFrame")
cont4_shape <- as(cont4, "SpatialLinesDataFrame")
# 将SpatialLinesDataFrame对象转换为sf对象
cont1_sf <- st_as_sf(cont1_shape)
cont2_sf <- st_as_sf(cont2_shape)
cont3_sf <- st_as_sf(cont3_shape)
cont4_sf <- st_as_sf(cont4_shape)
# 导出为线shp文件
st_write(cont1_sf, "E:/ARCgisdata/TPG2023/T4/T4_10.shp")
st_write(cont2_sf, "E:/ARCgisdata/TPG2023/T4/T4dBBMMall.shp")
st_write(cont3_sf, "E:/ARCgisdata/TPG2023/contour_linesbao150.shp")
st_write(cont4_sf, "E:/ARCgisdata/TPG2023/contour_linesbao250.shp")

writeRaster(deer_dbbmm3@layers[[1]], 
            filename = "E:/ARCgisdata/TPG2023/T4/T4dbbmm1.shp", 
            format = "GTiff") 
# 将RasterLayer转换为SpatialPolygonsDataFrame
deer_polygons <- rasterToPolygons(deer_dbbmm3@layers[[1]], 
                                  dissolve = TRUE)

# 导出为shapefile
writeOGR(deer_polygons,
         "E:/ARCgisdata/TPG2023/T4/", 
         "T4dbbmm1", 
         driver = "ESRI Shapefile")



# 计算面积 和gis算出来不一样--------------------------------------------------------------------
deer_data3 <- read.csv2("T4data_25.csv",
                        header = TRUE,
                        sep = ",")
deer_move3 <- move(x = as.numeric(deer_data3$location.long),
                   y = as.numeric(deer_data3$location.lat),
                   time=as.POSIXct(deer_data3$t, format="%Y/%m/%d %H:%M"), 
                   data=deer_data3, 
                   proj=CRS("+proj=longlat +ellps=WGS84"),
                   animal=deer_data3$id,
                   sensor=deer_data3$sensor)
table(deer_data3$id)  # 查看每个季节有多少条记录

deer_move3 <- spTransform(deer_move3,
                          CRS("+proj=utm +zone=52 +datum=WGS84"))
deer_dbbmm3 <- brownian.bridge.dyn(object=deer_move3, 
                                   location.error=20, 
                                   window.size=31,
                                   margin=9, 
                                   dimSize=100,
                                   time.step=60,
                                   ext = 0.5
)
cont1 <- raster2contour(deer_dbbmm3, level=c(0.95))
cont1_sf <- st_as_sf(cont1)
plot(deer_dbbmm3)
contour(deer_dbbmm3, levels=0.95, add=TRUE)
# 将RasterLayer转换为SpatialPolygonsDataFrame
deer_polygons <- rasterToPolygons(deer_dbbmm3, 
                                  dissolve = TRUE)

# 导出为栅格shapefile
writeOGR(deer_polygons,
         "E:/ARCgisdata/TPG2023/T4/", 
         "EMD_13", 
         driver = "ESRI Shapefile")
# 导出为等值线shp文件
st_write(cont1_sf, "E:/ARCgisdata/TPG2023/T4/EMD_25.shp")

# 画EMD重叠图 -----------------------------------------------------------------
library(ggspatial)
library(ggsn)

deer_data <- read.csv2("T4/T4data_36.csv",
                        header = TRUE,
                        sep = ",")
deer_move <- move(x = as.numeric(deer_data$location_long),
                   y = as.numeric(deer_data$location_lat),
                   time=as.POSIXct(deer_data$t, format="%Y/%m/%d %H:%M"), 
                   data=deer_data, 
                   proj=CRS("+proj=longlat +ellps=WGS84"),
                   animal=deer_data$id,
                   sensor=deer_data$sensor)
table(deer_data$id)  # 查看每个季节有多少条记录

deer_move <- spTransform(deer_move,
                          CRS("+proj=utm +zone=52 +datum=WGS84"))
deer_dbbmm <- brownian.bridge.dyn(object=deer_move, 
                                   location.error=20, 
                                   window.size=31,
                                   margin=9, 
                                   dimSize=100,
                                   time.step=60,
                                   ext = 0.3
)
cont <- raster2contour(deer_dbbmm, level=c(0.95))
cont_sf <- st_as_sf(cont)
deer_data3 <- read.csv2("T4/T4data_37.csv",
                        header = TRUE,
                        sep = ",")
deer_move3 <- move(x = as.numeric(deer_data3$location_long),
                   y = as.numeric(deer_data3$location_lat),
                   time=as.POSIXct(deer_data3$t, format="%Y/%m/%d %H:%M"), 
                   data=deer_data3, 
                   proj=CRS("+proj=longlat +ellps=WGS84"),
                   animal=deer_data3$id,
                   sensor=deer_data3$sensor)
table(deer_data3$id)  # 查看每个季节有多少条记录

deer_move3 <- spTransform(deer_move3,
                          CRS("+proj=utm +zone=52 +datum=WGS84"))
deer_dbbmm3 <- brownian.bridge.dyn(object=deer_move3, 
                                   location.error=20, 
                                   window.size=31,
                                   margin=9, 
                                   dimSize=100,
                                   time.step=60,
                                   ext = 0.5
)
cont1 <- raster2contour(deer_dbbmm3, level=c(0.95))
cont1_sf <- st_as_sf(cont1)

polygons <- st_cast(cont_sf, "POLYGON")
polygons1 <- st_cast(cont1_sf, "POLYGON")

# 创建ggplot对象并绘制多边形
# 用dbbmmstack画图尝试
cont2 <- raster2contour(deer_dbbmm3, level=c(0.95))
cont2_sf <- st_as_sf(cont2)
polygons <- st_cast(cont2_sf, "POLYGON")
plot(polygons)
# 检查CRS
st_crs(polygons)
st_crs(polygons) <- 4326
plot(polygons$geometry[[2]])
plot(deer_dbbmm3@layers[[2]])

ggplot() +
  geom_sf(data = polygons, aes(fill = "period 68"), alpha = 0.6) +
  geom_sf(data = polygons1, aes(fill = "period 69"), alpha = 0.6) +
  scale_fill_manual(values = c("period 68" = "blue", "period 69" = "red"),
                    name = "seven-day period") +
  ggtitle("d) EMD = 4101") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 90, vjust = 0.5),
        axis.ticks.x = element_blank()) +
  theme(legend.position = "bottom") +
  labs(x = "Longitude", y = "Latitude") +
  annotation_north_arrow(location = "br",
                         which_north = "true",
                         pad_x = unit(1.9, "in"),
                         pad_y = unit(2.4, "in"), 
                         style = north_arrow_fancy_orienteering,
                         scale = 0.5) +
  annotation_scale(location = "br", 
                   width_hint = 0.5)

plot(all_deer_dbbmm[[61]])
plot(deer_dbbmm)

move::emd(deer_dbbmm,deer_dbbmm3)
#大圆
#emd(deer_dbbmm, deer_dbbmm3, gc = TRUE)


# emd代码（循环版 ---------------------------------------------------------------

# 获取当前工作目录中所有匹配"T4data_*.csv"模式的文件
files <- files <- list.files("T4", pattern = "T4data_\\d+\\.csv", full.names = TRUE)

# 初始化列表和向量来存储输出和EMD结果
all_cont_sf <- list()
all_id_counts <- list()
all_deer_dbbmm <- list()
emd_results <- numeric(length(files) - 1)  # 假设有80个文件，因此有79个EMD结果

previous_dbbmm <- NULL  # 用于存储前一个deer_dbbmm对象

# 循环处理每个文件
for (i in seq_along(files)) {
  file_name <- files[i]
  deer_data <- read.csv2(file_name, header = TRUE, sep = ",")
  
  deer_move <- move(x = as.numeric(deer_data$location_long),
                    y = as.numeric(deer_data$location_lat),
                    time = as.POSIXct(deer_data$t, format="%Y/%m/%d %H:%M"), 
                    data = deer_data, 
                    proj = CRS("+proj=longlat +ellps=WGS84"),
                    animal = deer_data$id,
                    sensor = deer_data$sensor)
  
  deer_move <- spTransform(deer_move, CRS("+proj=utm +zone=52 +datum=WGS84"))
  
  deer_dbbmm <- brownian.bridge.dyn(object = deer_move, 
                                    location.error = 20, 
                                    window.size = 17,
                                    margin = 5, 
                                    dimSize = 100,
                                    time.step = 60,
                                    ext = 0.5)
  
  # 生成等高线并转换为sf对象
  cont <- raster2contour(deer_dbbmm, level = c(0.95))
  cont_sf <- st_as_sf(cont)
  
  # 保存每个文件的结果
  all_cont_sf[[i]] <- cont_sf
  all_deer_dbbmm[[i]] <- deer_dbbmm
  all_id_counts[[i]] <- table(deer_data$id)
  

}

#library(emdist)

# 初始化一个向量来存储EMD结果
emd_results <- numeric(length(all_deer_dbbmm) - 1)

# 计算连续的EMD
for (i in 1:(length(all_deer_dbbmm) - 1)) {
  # 由于deer_dbbmm对象是根据BBMM计算的，我们需要提取适合计算EMD的数据
  # 假设我们直接可以计算两个deer_dbbmm对象之间的EMD（这可能需要根据实际情况调整）
  emd_results[i] <- move::emd(all_deer_dbbmm[[i]], all_deer_dbbmm[[i + 1]])
}
# 保存ID计数和EMD结果等

write.csv(emd_results, "T4/emd_results.csv")

plot(all_deer_dbbmm[[37]])


# 转面法2 --------------------------------------------------------------------

# 检查几何类型
geom_type <- st_geometry_type(cont1_sf)
print(geom_type)

# 如果cont1_sf确实包含MULTILINESTRING，下面的示例展示了如何尝试转换
# 请注意，以下代码是假定所有MULTILINESTRING都是闭合且可以构成有效多边形的示例
# 实际情况可能需要根据你的数据进行调整

# 假设cont1_sf只包含一个MULTILINESTRING且它是闭合的
if (all(geom_type == "MULTILINESTRING")) {
  # 尝试将MULTILINESTRING转换为POLYGON
  polygons <- st_cast(cont1_sf, "POLYGON")
  
  # 计算每个转换后的POLYGON的面积
  areas <- st_area(polygons)
  
  print(areas)
} else {
  print("cont1_sf contains geometry types other than MULTILINESTRING or conversion might not be straightforward.")
}

# 凸包法 ---------------------------------------------------------------------

#因未cont是多段线，不是面，多所以算不出面积
# 线转面
# 转换MULTILINESTRING到MULTIPOINT
cont1_points <- st_cast(cont1_sf, "MULTIPOINT")

# 创建凸包
convex_hull <- st_convex_hull(cont1_points)

# 检查凸包的几何类型
print(st_geometry_type(convex_hull))

# 计算凸包的面积
area_convex_hull <- st_area(convex_hull)
print(paste("Area of convex hull:", area_convex_hull, "square units"))
sum(area_convex_hull/1000000)


# 利用分布 --------------------------------------------------------------------
deer_data3 <- read.csv2("T4data_11.csv",
                        header = TRUE,
                        sep = ",")
deer_move3 <- move(x = as.numeric(deer_data3$location.long),
                   y = as.numeric(deer_data3$location.lat),
                   time=as.POSIXct(deer_data3$t, format="%Y/%m/%d %H:%M"), 
                   data=deer_data3, 
                   proj=CRS("+proj=longlat +ellps=WGS84"),
                   animal=deer_data3$id,
                   sensor=deer_data3$sensor)
table(deer_data3$id)  # 查看每个季节有多少条记录

deer_move3 <- spTransform(deer_move3,
                          CRS("+proj=utm +zone=52 +datum=WGS84"))
deer_dbbmm3 <- brownian.bridge.dyn(object=deer_move3, 
                                   location.error=20, 
                                   window.size=31,
                                   margin=9, 
                                   dimSize=100,
                                   time.step=60,
                                   ext = 0.5
)
deer_ud <- getVolumeUD(deer_dbbmm3)
deer_ud <- getVolumeUD(deer_dbbmm)

par(mfrow=c(1,1))
plot(deer_ud, main="UD")
plot(deer_dbbmm3)
## also a contour can be added
plot(deer_ud, main="UD and contour lines")
contour(deer_ud, 
        levels=c(0.5, 0.95), 
        add=TRUE, 
        lwd=c(0.5, 0.5), 
        lty=c(2,1))
par(mfrow=c(1,1))

## mantaining the lower probabilities
ud95 <- deer_ud
ud95[ud95>.95] <- NA
plot(ud95, main="UD95")
contour(ud95, 
        levels=c(0.5, 0.94), 
        add=TRUE, 
        lwd=c(0.5, 0.5), 
        lty=c(2,1))
## or extracting the area with a given probability, where cells that belong to the given probability will get the value 1 while the others get 0
ud95 <- deer_ud<=.95
plot(ud95, main="UD95")

ud50 <- deer_ud<=.5
plot(ud50, main="UD50")
contour(ud95, 
        levels=c(0.5, 0.95), 
        add=TRUE, 
        lwd=c(0.5, 0.5), 
        lty=c(2,1))
par(mfrow=c(1,1))
# 将RasterLayer转换为SpatialPolygonsDataFrame
ud_polygons <- rasterToPolygons(ud95, 
                                  dissolve = TRUE)

# 导出为shapefile
writeOGR(ud_polygons,
         "E:/ARCgisdata/TPG2023/T4/", 
         "EMD_12", 
         driver = "ESRI Shapefile")
str(ud95)
deer <- deer_data3[1:100]
summary(timeLag(deer,"mins")) 

#计算布朗运动方差
deer_dbmvar <- brownian.motion.variance.dyn(object = deer_move3,
                                            location.error = 20,
                                            window.size=31,
                                            margin=9)

## the intended GPS fix rate of leroy was 15min, so we will ignore for example all segments that have a larger time lag than 5hours. The 'dBMvariance' object resulting from the function above, contains the slot '@interest' in which those segments marked as FALSE won't be included in the calculation of the dBBMM. Therefore we set all segments with time lag larger than 300mins to false
deer_dbmvar@interest[timeLag(deer_move3,"hours")>2] <- FALSE

## then we use the 'dBMvariance' object to calculate the dBBMM
dbb.corrected <- brownian.bridge.dyn(deer_dbmvar, raster=100, ext=.45,location.error=20)
#计算布朗运动方差
dbbvar_1 <- getMotionVariance(deer_dbbmm3)
dbbvar_2 <- getMotionVariance(dbb.corrected)
log_dbbvar_1 <- log(dbbvar_1)
plot(dbbvar_1)
boxplot(log_dbbvar_1)
## now the UD makes more sense
ud.corrected <- getVolumeUD(dbb.corrected)

par(mfrow=c(1,1))
plot(ud.corrected, main="UD")
contour(ud.corrected, levels=c(0.5, 0.95), add=TRUE, lwd=c(0.5, 0.5), lty=c(2,1))

plot(ud.corrected, main="UD with locations")
points(deer_move3, col="red", cex=.1, pch=20)
contour(ud.corrected, levels=c(0.5, 0.95), add=TRUE, lwd=c(0.5, 0.5), lty=c(2,1))

# 插值 ----------------------------------------------------------------------
# 利用分布 --------------------------------------------------------------------
deer_data3 <- read.csv2("T1data_4new.csv",
                        header = TRUE,
                        sep = ",")
deer_move3 <- move(x = as.numeric(deer_data3$location_long),
                   y = as.numeric(deer_data3$location_lat),
                   time=as.POSIXct(deer_data3$t, format="%Y/%m/%d %H:%M"), 
                   data=deer_data3, 
                   proj=CRS("+proj=longlat +ellps=WGS84"),
                   animal=deer_data3$id,
                   sensor=deer_data3$sensor)
## 按照固定时间间隔插值
interp1hour <- interpolateTime(deer_move3, 
                               time=as.difftime(1, units="hours"), 
                               spaceMethod='rhumbline')
plot(deer_move3, 
     col="red",
     pch=20, 
     main="By time interval")
points(interp1hour)
lines(deer_move3, col="red")
legend("bottomleft",
       c("True locations", "Interpolated locations"), 
       col=c("red", "black"), 
       pch=c(20,1))
summary(timeLag(interp1hour, "hours"))
write.csv(interp1hour@coords, "T1data_4inter1hour.csv")
#自定义插值
ts <- as.POSIXct(c("2023-02-12 00:00", "2023-02-12 01:00", "2023-02-12 02:00",
                   "2023-02-13 00:00", "2023-02-13 01:00", "2009-02-13 02:00"),
                 format="%Y-%m-%d %H:%M", tz="Asia/Shanghai")

interpCusom <- interpolateTime(deer_move3, time=ts, spaceMethod='euclidean')

plot(deer_move3, col="red",pch=20, main="By custom timestamps")
points(interpCusom)
lines(deer_move3, col="red")
legend("bottomleft", c("True locations", "Interpolated locations"), col=c("red", "black"), pch=c(20,1))



