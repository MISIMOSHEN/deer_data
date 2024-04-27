# Load package
library(ctmm)
library(sp)
library(adehabitatLT)
library(terra)
library(sf)
library(stars)
library(move)
library(raster)
library(ggplot2)
# Load example buffalo data
deer_data <- read.csv2("T4data_9.csv",
  header = TRUE,
  sep = ","
)
# 将 "location.long" 和 "location.lat" 列转换为数值类型
deer_data$location.long <- as.numeric(deer_data$location.long)
deer_data$location.lat <- as.numeric(deer_data$location.lat)
# Combine longitude, latitude and time data into a data.frame
deer_move <- move(
  x = as.numeric(deer_data$location.long),
  y = as.numeric(deer_data$location.lat),
  time = as.POSIXct(deer_data$t, format = "%Y/%m/%d %H:%M"),
  data = deer_data,
  proj = CRS("+proj=longlat +ellps=WGS84"),
  animal = deer_data$id,
  sensor = deer_data$sensor
)

deer_data$t <- as.POSIXct(strptime(as.character(deer_data$t), "%Y/%m/%d %H:%M", tz = "Asia/Shanghai"))
deer_telemetry <- as.telemetry(deer_move, 
                               timeformat = "POSIXct",
                               tz = "Asia/Shanghai")

attr(deer_data$t, "tz")


data(buffalo)
Pepper <- buffalo$Peppe

data("buffalo")
# Extract data for buffalo 1, Cilla
cilla <- buffalo[[1]]
# Plot the positions

# ctmm，启动！ ----------------------------------------------------------------


plot(deer_telemetry,color = deer_telemetry$season)
# 使用ggplot2为散点图设置颜色
ggplot(data = deer_data, aes(x = location.long, y = location.lat, color = season)) +
  geom_point(size = 0.5) +
  theme_minimal()  # 可选：添加一个简单的主题

# Calculate variogram
vg.deer <- variogram(deer_telemetry,
                     dt = 3600)
# Plot up to 50% of the maximum lag in the data
level <- c(0.5, 0.95)
plot(vg.deer, level = level)
# Zoom in on the shortest lags
plot(vg.deer, fraction = 0.005, level = level)

# Use the sliders provided by variogram.fit to specifystarting values.
# The default choices are usually acceptable.
variogram.fit(vg.deer, fraction = 1)
# using the initial parameter values obtained from
#variogram.fit(vg.deer)

fitted.mods <- ctmm.select(deer_telemetry,
  CTMM = GUESS,
  verbose = TRUE
)

summary(fitted.mods)

# Extract the fitted anisotropic version of IID, OU,and OUF.
iid <- fitted.mods[[5]]
ou <- fitted.mods[[3]]
ouf <- fitted.mods[[2]]
plot(vg.deer, CTMM = iid, col.CTMM = "#1b9e77")
plot(vg.deer,
  CTMM = iid, col.CTMM = "#1b9e77",
  fraction = 0.005
)
plot(vg.deer, CTMM = ou, col.CTMM = "#1b9e77")
plot(vg.deer,
  CTMM = ou, col.CTMM = "#1b9e77",
  fraction = 0.005
)
plot(vg.deer, CTMM = ouf, col.CTMM = "#1b9e77")
plot(vg.deer,
  CTMM = ouf, col.CTMM = "#1b9e77",
  fraction = 0.005
)

summary(iid)
summary(ou)
summary(ouf)
# Conventional KDE estimate
kde.deer <- akde(deer_telemetry, CTMM = iid, dt = 3600)
# Autocorrelated KDE estimate
akde.deer <- akde(deer_telemetry, CTMM = ouf, weights = T, dt = 3600)
summary(kde.deer)
summary(akde.deer)
plot(deer_telemetry, UD = kde.deer)
title(expression("IID KDE"))
plot(kde.deer)
plot(deer_telemetry, UD = akde.deer)
title(expression("OUF AKDE"))
plot(akde.deer)

summary(akde.deer)


#导出tif，暂时不需要
#writeRaster(akde.deer, filename = "E:/ARCgisdata/TPG2023/T4/T4ctmm10.shp",format = "GTiff")
#转为sp
spdf_deer <- SpatialPolygonsDataFrame.UD(akde.deer, UD.level = 0.95)
# 指定输出路径和文件名
output_path <- "E:/ARCgisdata/TPG2023/T4/T4ctmm9"

# 使用writeOGR函数将SpatialPolygonsDataFrame对象写入为Shapefile
writeOGR(spdf_deer, 
         dsn = dirname(output_path), 
         layer = basename(output_path), 
         driver = "ESRI Shapefile")





# 画图尝试 --------------------------------------------------------------------
# 加载必要的包
library(leaflet)
library(sp)

bounds <- as.vector(bbox(extent(deer_move)))
base_map <- leaflet() %>%
  fitBounds(bounds[1], bounds[2], bounds[3], bounds[4]) %>%
  addTiles() %>% 
  addProviderTiles(providers$Esri.WorldImagery)
map2 <- base_map %>% 
  addPolylines(data = as(deer_move, 
                         "SpatialLines"), 
               color = "grey") %>%
  addCircles(data = deer_move, 
             fillOpacity = 0.3,
             opacity = 0.5, 
             color = "blue") %>%
  addPolygons(data = akde.deer,
              fillOpacity = 0.3,
              opacity = 0.5, 
              color = "red")
map3 <- map2 %>% addLegend(position = "topright",
                           colors = c("grey", "blue"),
                           labels = c("lines", "points"), 
                           opacity = 0.7,
                           title = "DEER")
map4 <- map3 %>% addScaleBar(position = "bottomleft",
                             options = scaleBarOptions(maxWidth = 100, 
                                                       metric = TRUE, 
                                                       imperial = F, 
                                                       updateWhenIdle = TRUE))
map4




# _______________________________________________________
data(buffalo)


Pepper <- buffalo$Pepper
M.IID <- ctmm.fit(Pepper) # no autocorrelation timescales
GUESS <- ctmm.guess(Pepper, interactive = FALSE) # automated model guess
M.OUF <- ctmm.fit(Pepper, GUESS) # in general, use ctmm.select instead
KDE <- akde(Pepper, M.IID) # KDE
AKDE <- akde(Pepper, M.OUF) # AKDE
wAKDE <- akde(Pepper, M.OUF, weights = TRUE) # weighted AKDE
# calculate one extent for all UDs
EXT <- extent(list(KDE, AKDE, wAKDE), level = 0.95)

plot(Pepper, UD = KDE, xlim = EXT$x, ylim = EXT$y)
title(expression("IID KDE"["C"]))
plot(Pepper, UD = AKDE, xlim = EXT$x, ylim = EXT$y)
title(expression("OUF AKDE"["C"]))
plot(Pepper, UD = wAKDE, xlim = EXT$x, ylim = EXT$y)
title(expression("weighted OUF AKDE"["C"]))
