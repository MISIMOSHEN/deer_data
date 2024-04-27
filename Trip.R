library(trip)
library(tidyverse)
library(move)
library(adehabitatHR)
library(adehabitatLT)
library(adehabitatMA)
library(sp)
library(sf)
library(terra)
library(stars)

# 导入数据，转化格式 ---------------------------------------------------------------
deer_data <- read.csv2("T1data_4new.csv",
                       header = TRUE,
                       sep = ",")
# 将 "location_long" 和 "location_lat" 列转换为数值类型
deer_data$location_long <- as.numeric(deer_data$location_long)
deer_data$location_lat <- as.numeric(deer_data$location_lat)

#转为POSIXct
deer_data$t <- as.POSIXct(strptime(as.character(x = deer_data$t),
                                   format = "%Y/%m/%d %H:%M",
                                   tz = "Asia/Shanghai"))
#转为POSIXlt
#deer_data$t <- as.POSIXlt(deer_data$t)

#转为UTM
#将deer_data转换为sf对象并设置EPSG 4326（WGS84）
deer_sf <- st_as_sf(deer_data, 
                    coords = c("location_long", "location_lat"), 
                    crs = 4326)

# 转换坐标为UTM Zone 52
deer_sf_utm <- st_transform(deer_sf,
                            crs = 32652)  # 32652是UTM Zone 52的EPSG代码

# 如果需要，您可以将UTM坐标重新添加到deer_data数据框
deer_data$utm_x <- st_coordinates(deer_sf_utm)[, 1]
deer_data$utm_y <- st_coordinates(deer_sf_utm)[, 2]
deer_data$utm_x <- deer_data$utm_x/1000#m变成km
deer_data$utm_y <- deer_data$utm_y/1000

str(deer_data)
#转为ltraj
coordinates <- deer_data[, c("location_long", "location_lat")]
coordinates_utm <- deer_data[, c("utm_x", "utm_y")]
deer_ltraj_utm <- as.ltraj(xy = coordinates_utm,
                          date = deer_data$t, 
                          id = deer_data$id,
                          proj4string = CRS("+proj=utm +zone=52 +ellps=WGS84"))
deer_ltraj <- as.ltraj(xy = coordinates,
                      date = deer_data$t, 
                      id = deer_data$id,
                      proj4string = CRS("+proj=longlat +ellps=WGS84")
                      )

deer_ltraj
deer_ltraj_utm
plotltr(deer_ltraj, "dist")
plotltr(deer_ltraj_utm, "dist")
head(deer_ltraj)
summary(deer_ltraj[[1]])
summary(deer_ltraj_utm[[1]])#与moveHMM比较可以发现UTM版是正确的
trajdyn(deer_ltraj_utm)#互动
#转为ltraj.trip
#deer_ltraj <- as.ltraj.trip(xy = deer_data[,c("location.long","location.lat")])

#转为trip
deer_trip <- as.trip(deer_ltraj)
deer_trip_utm <- as.trip(deer_ltraj_utm)
trip <- trip(deer_data,TORnames = "DeerTrip")

summary(deer_trip)
summary(deer_trip_utm)
#对比一下dist
summary(deer_trip[[4]])
summary(deer_trip_utm[[4]])
plot(deer_trip[[4]])
plot(deer_trip_utm[[4]])
#对比一下rel.angle
summary(deer_trip[[8]])
summary(deer_trip_utm[[8]])
plot(deer_trip[[8]])
plot(deer_trip_utm[[8]])


# 过滤，启动！ ------------------------------------------------------------------
deer_filtered <- sda(x = deer_trip,
                     smax = 4, #超过此速度被删
                     ang = c(1, 75), #在此范围内的角度被删
                     distlim = c(0, 0.2),#不在此范围内的步长被删
                     pre = NULL)
#想运行下面一步就不能把value转为dataframe
#deer_filtered <- as.data.frame(deer_filtered)
summary(deer_filtered)

write.csv(x = deer_filtered,
          file = "T26/T26_filtered.csv")
# deer_filtered 数据框（其中包含逻辑值 TRUE 和 FALSE）作为逻辑条件来子集 deer_data 数据框。
#这将返回一个新的数据框 deer_data_filtered，
#其中只包含 deer_filtered 中为 TRUE 的行。
#deer_data <- read.csv2("T23/T23data.csv",
#                       header = TRUE,
 #                      sep = ",")
deer_data_filtered <- deer_data[deer_filtered, ]

str(deer_data_filtered)
deer_ltraj_filtered <- as.ltraj(xy = deer_data_filtered[, c("location_long", "location_lat")],
                       date = deer_data_filtered$t, 
                       id = deer_data_filtered$id,
                       proj4string = CRS("+proj=longlat +ellps=WGS84"))

deer_trip_filtered <- as.trip(deer_ltraj_filtered)
#比较清洗后的数据
plot(deer_trip_filtered)
plot(deer_trip)

plot(deer_trip)
plot(deer_trip_filtered)
#导出csv
write.csv(x = deer_data_filtered,
          file = "T26/T26filtered_data.csv")


# LEAFLET 地图 --------------------------------------------------------------
# 在线交互地图地图
require(move)
library(leaflet)
library(sp)
library(lubridate)
library(shiny)
#清洗完之后的数据记得先在excel里调整时间格式！
deer_data4 <- read.csv("T32/T32data.csv",
                       header=TRUE,
                       sep=",")

deer_move4 <- move(x = deer_data4$location_long, 
                   y = deer_data4$location_lat,
                   time = as.POSIXct(deer_data4$t, format="%Y/%m/%d %H:%M"), 
                   data = deer_data4,
                   proj = CRS("+proj=longlat +ellps=WGS84"),
                   animal = deer_data4$id, 
                   sensor = deer_data4$sensor)


#需要经纬度坐标，不能用：deer_move3 <- spTransform(deer_move3, CRS("+proj=utm +zone=52 +datum=WGS84"))
bounds <- as.vector(bbox(extent(deer_move4)))
base_map <- leaflet() %>%
  fitBounds(bounds[1], bounds[2], bounds[3], bounds[4]) %>%
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery)
map2 <- base_map %>% 
  addPolylines(data = as(deer_move4, 
                         "SpatialLines"), 
               color = "grey") %>%
  addCircles(data = deer_move4, 
             fillOpacity = 0.3,
             opacity = 0.5, 
             color = "blue")
map3 <- map2 %>% addLegend(position = "topright",
                           colors = c("grey", "blue"),
                           labels = c("lines", "points"), 
                           opacity = 0.7,
                           title = "DEERT32")
map4 <- map3 %>% addScaleBar(position = "topleft",
                             options = scaleBarOptions(maxWidth = 100, 
                                                       metric = TRUE, 
                                                       imperial = F, 
                                                       updateWhenIdle = TRUE))
map4
base_map

# 筛选速度 --------------------------------------------------------------------


deer_speedfilter <- speedfilter(x = deer_trip,
                                max.speed = 4,
                                test = FALSE)


deer_speedfilter <- as.data.frame(deer_speedfilter)

write.csv(x = deer_filtered,
          file = "deer_speedfilter.csv"
)

vignette("trip")

