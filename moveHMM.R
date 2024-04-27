
library(moveHMM)
deer_data <-read.csv("T23/T23data.csv",
           header = TRUE,
           sep = ",")

deer_Geo <- deer_data
#wgs84转UTM
xy <-
  data.frame(x = deer_Geo$location_long, y = deer_Geo$location_lat)
coordinates(xy) <- c("x", "y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
deer_utm <-
  suppressWarnings(spTransform(xy, CRS("+proj=utm +zone=52 +ellps=WGS84")))


deer_utm <- as.data.frame(deer_utm)

deer_Geo$location_long <- deer_utm$x
deer_Geo$location_lat <- deer_utm$y
colnames(deer_Geo)[which(names(deer_Geo) == "location_long")] <- "x"
colnames(deer_Geo)[which(names(deer_Geo) == "location_lat")] <- "y"
deer_utm <- deer_Geo
deer_Geo <- deer_data
str(deer_utm)
str(deer_Geo)

deer_utm$x <- deer_utm$x / 1000
deer_utm$y <- deer_utm$y / 1000
head(deer_utm)
#utm预处理
processeddeer_utm <- prepData(deer_utm, type = "UTM")
#经纬度预处理
processeddeer_Geo <- prepData(deer_Geo,
                              type = "LL",
                              coordNames = c("location_long", "location_lat"))
summary(processeddeer_utm$step)
summary(processeddeer_Geo$step)
head(processeddeer_utm$step)
head(processeddeer_Geo$step)


head(processeddeer)
hist(processeddeer$step, xlab = "step length")
write.csv(x = processeddeer, file = "C:/R/R-4.2.2/step.csv")
hist(processeddeer$angle,
     breaks = seq(-pi, pi, length = 15),
     xlab = "turning angle")
summary(processeddeer)
plot(processeddeer)
write.csv(processeddeer, file = "C:/R/R-4.2.2/WDS2.csv")


## initial parameters for gamma and von Mises distributions
mu0 <- c(0.05,0.4,0.8) # step mean (two parameters: one for each state)
sigma0 <- c(0.05,0.4,0.8) # step SD
zeromass0 <- c(0.05,0.05,0.05) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)
angleMean0 <- c(pi,pi,0) # angle mean
kappa0 <- c(1,1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)
## call to fitting function
m0 <- fitHMM(data=processeddeer,nbStates=3,stepPar0=stepPar0,
            anglePar0=anglePar0)

mu1 <- c(0.05,0.28,0.6) # step mean (two parameters: one for each state)
sigma1 <- c(0.01,0.28,0.6) # step SD

stepPar1 <- c(mu1,sigma1)
angleMean1 <- c(pi,pi,0) # angle mean
kappa1 <- c(1.514845,1.514845,1.514845) # angle concentration
anglePar1 <- c(angleMean1,kappa1)
## call to fitting function
m1 <- fitHMM(data=processeddeer,nbStates=3,stepPar0=stepPar1,
             anglePar0=anglePar1)

mu2 <- c(0.05,0.2,0.5) # step mean (two parameters: one for each state)
sigma2 <- c(0.05,0.2,0.5) # step SD
zeromass2 <- c(0.01,0.01,0.01) # step zero-mass
stepPar2 <- c(mu2,sigma2,zeromass2)
angleMean2 <- c(pi,pi,0) # angle mean
kappa2 <- c(0.05,0.2,0.5) # angle concentration
anglePar2 <- c(angleMean2,kappa2)
## call to fitting function
m2 <- fitHMM(data=processeddeer,nbStates=3,stepPar0=stepPar2,
             anglePar0=anglePar2)

mu3 <- c(0.3,0.7) # step mean (two parameters: one for each state)
sigma3 <- c(0.3,0.6) # step SD
zeromass3 <- c(0.01,0.01) # step zero-mass
stepPar3 <- c(mu3,sigma3,zeromass3)
angleMean3 <- c(pi,0) # angle mean
kappa3 <- c(3,1) # angle concentration
anglePar3 <- c(angleMean3,kappa3)
## call to fitting function
m3 <- fitHMM(data=processeddeer,nbStates=2,stepPar0=stepPar3,
             anglePar0=anglePar3)


AIC(m0,m1,m2,m3)
AIC(m1)
# Estimated step length parameters
stepMean <- m1$mle$stepPar["mean",]
stepSD <- m1$mle$stepPar["sd",]
# Estimated turning angle parameters
angleMean <- m1$mle$anglePar["mean",]
angleCon <- m1$mle$anglePar["concentration",]
stepShape <- stepMean^2/stepSD^2
stepRate <- stepMean/stepSD^2
# Most likely state sequence
states <- viterbi(m1)

m1
CI(m1)
plot(m1, plotCI=TRUE)

 
states <- viterbi(m1)

states[1:25]
sp <- stateProbs(m1)
head(sp)
plot(m1)
plotStates(m1,animals="deer")
plotStationary(m1, plotCI=TRUE)
tail(states)
# compute the pseudo-residuals
pr <- pseudoRes(m1)
# time series, qq-plots, and ACF of the pseudo-residuals
plotPR(m1)
utm <- read.csv("C:/R/R-4.2.2/utm.csv",header=TRUE,sep=",")
summary(utm)

library(sp)
library(maptools)
library(rgdal)
library(rgdal)
# Extract the state variable
states <- statePath(m1)
# Add the state variable as an attribute to the processeddeer data frame
processeddeer$state <- state
# Convert processeddeer to a SpatialPointsDataFrame
coordinates(processeddeer) <- c("x","y")
proj4string(processeddeer) <- CRS("+proj=utm +zone=10 +datum=WGS84")
# Write the processeddeer SpatialPointsDataFrame to a Shapefile file
writeOGR(obj = processeddeer, 
         dsn = "C:/R/R-4.2.2/deer_states.shp", 
         layer = "deer_states", 
         driver = "ESRI Shapefile", 
         overwrite_layer = TRUE)

angles <- ifelse(processeddeer$angle < 0, pi + processeddeer$angle, processeddeer$angle)
library(circular)
angles_circ <- as.circular(na.omit(angles), units = "radians", template = "none",modulo="asis")
fit <- mle.vonmises(angles_circ)
fit$kappa


library(circular)
data_circ <- circular(runif(100, -pi, pi))
fit <- mle.vonmises(data_circ)
fit$kappa
