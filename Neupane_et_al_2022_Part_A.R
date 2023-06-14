
#######     PAPER REFERENCE #####################################################################
# Neupane, N., Peruzzi, M., Arab, A., Mayor, S.J., Withey, J.C., Ries, L. and Finley, A.O., 2022.
# A novel model to accurately predict continental-scale timing of forest green-up. 
# International Journal of Applied Earth Observation and Geoinformation, p.102747.
##################################################################################################
rm(list = ls())

library(sp)
library(geoR)
library(spNNGP)
library(fastDummies)

#k<-2001
for (k in 2002:2017)
{
        ls()
        rm(list = setdiff(ls(),"k"))
        ls()
        # file1 <- paste("~/Google Drive File Stream/My Drive/PredictGreenBayesian/DataUpto2017/CommonLocations/","EastOF100W_CommonLocations_",k,".csv", sep = "")
        file1 <- paste("~/Google Drive File Stream/My Drive/SpatTempGreenPred/Data/EastOF100W_CommonLocations_",k,".csv", sep = "")
        mydata <- read.csv(file1)
names(mydata)
head(mydata)
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

coords.ll <- mydata[,1:2]
coordinates(coords.ll) <- c("Lon", "Lat")
proj4string(coords.ll) <- "+proj=longlat +ellps=GRS80"
coords.aea <- coordinates(spTransform(coords.ll, CRS(aea)))
#plot(coords.aea)

##check for duplicates
any(duplicated(coords.aea))

set.seed(1)
n <- nrow(mydata)
n
## X is matrix of regressors
#FD <- fastDummies::dummy_cols(mydata$DominantCover)
#dim(FD)
#head(FD)
#X <- as.matrix(cbind(1, mydata$GDDMinus5, FD[,2:4]))
X <- as.matrix(cbind(1, mydata$GDDMinus5))

y <- mydata[,4]

sigma.sq <- 40
tau.sq <- 10

set.seed(1)
# About 20% sampling
ho <- sample(1:n, 22000)
y.ho <- y[ho]
x.ho <- X[ho,,drop=FALSE]
coords.ho <- coords.aea[ho,]

y <- y[-ho]
x <- X[-ho,,drop=FALSE]
coords <- coords.aea[-ho,]

# number of MCMC iterations, used in spNNGP

n.samples <- 10000

starting <- list("phi"=3/1000, "sigma.sq"=sigma.sq, "tau.sq"=tau.sq)

tuning <- list("phi"=0.1, "sigma.sq"=0.05, "tau.sq"=0.05)

priors <- list("phi.Unif"=c(0.0001/1, 1/0.1), "sigma.sq.IG"=c(2, sigma.sq), "tau.sq.IG"=c(2, tau.sq))

#priors <- list("phi.Unif"=c(3/2500, 3/25), "sigma.sq.IG"=c(2, sigma.sq), "tau.sq.IG"=c(2, tau.sq))

cov.model <- "exponential"
         
n.report <- 1000
#n.report <- 100

m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=4, n.report=n.report)

summary(m.r)

#plot(m.r$p.theta.samples)

burn.in <- floor(0.1*n.samples)

## newest version on cran uses predict in place of spPredict
#p.r <- spPredict(m.r, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=4, verbose=FALSE)
p.r <- spPredict(m.r, X.0 = X, coords.0 = coords.aea, n.omp.threads=4, verbose=FALSE)
mydata$BaysResd <- p.r$residuals
sum(mydata$BaysResd)
head(mydata)
#plot(apply(p.r$p.y.0, 1, mean), y, pch=19, col="blue")

summary(m.r)
Post_Mean = apply(p.r$p.y.0, 1, mean)
Post_SD = apply(p.r$p.y.0, 1, sd)
mydata$Post_Mean <- round(Post_Mean, digits = 0)
mydata$Post_SD <- Post_SD

##Fit linear Model just using GDD as predictor
model1 <- lm(mydata$MODIS~mydata$GDDMinus5)
mydata$residual1 <- model1$residuals
mydata$fittedVal <- model1$fitted.values
#mydata$residual1 <- round(mydata$residual1, digits = 0)
summary(model1)
sum(mydata$residual1)

#Fit linear Model using GDD and dominant cover types, as factors (Categorical) as predictor
#model2 <- lm(mydata$MODIS~mydata$GDDMinus5+as.factor(mydata$DominantCover))
#mydata$residual2 <- model2$residuals
#summary(model2)
#sum(mydata$residual2)

#Difference between the two linear model residuals
#mydata$ResidualM2minusM1 <- mydata$residual2 - mydata$residual1

#Difference between Post mean and Observed MODIS date of onset 
mydata$PostMeanMinusMODIS <- Post_Mean - mydata$MODIS

#DecidForest <- subset(mydata, DominantCover == "DeciduousForest")
#EvrgnForest <- subset(mydata, DominantCover == "EvergreenForest")
#MixedForest <- subset(mydata, DominantCover == "MixedForest")

#nrow(DecidForest)
#nrow(EvrgnForest)
#nrow(MixedForest)

#######################################################################################################
#######                         The END, below this line are some plots                     ###########
#######################################################################################################
head(mydata)
mydata
names(mydata)
#k<-2001
myWD <-paste0("/Volumes/GoogleDrive/My Drive/PredictGreenBayesian/15Jun2020/Figures/",k,"/")
setwd(myWD)
#k<-2001
#setwd("~/Google Drive File Stream/My Drive/PredictGreenBayesian/15Jun2020/Figures")
getwd()

#Scatterplot (Compare BayesianModelPredicted VS MODIS observed Date of onset of greenness !)
BayesianPredictedOnset <- mydata$Post_Mean
MODISobservedOnset <- mydata$MODIS

all = c(BayesianPredictedOnset, MODISobservedOnset)
range = c(min(all), max(all))
mytitle <- paste0("BayesianFit_",k);mytitle

jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
        plot(BayesianPredictedOnset, MODISobservedOnset, xlim=range, ylim=range, pch=1, col="blue",cex = .3, main = mytitle, asp = 1)
        abline(fit <- lm(BayesianPredictedOnset ~ MODISobservedOnset), col='red')
        legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
dev.off()

######################
#Scatterplot (Compare Model 1 predicted NO LULC Onset VS MODIS observed Date of onset of greenness !)
mydata$M1fitted <- model1$fitted.values
MODISobservedOnset <- mydata$MODIS
M1_Predicted <- mydata$M1fitted

mytitle <- paste0("M1.GDD_only_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
all = c(M1_Predicted, MODISobservedOnset)
range = c(min(all), max(all))
plot(M1_Predicted, MODISobservedOnset, xlim=range, ylim=range, pch=1, col="blue",cex = .3, main = mytitle)
abline(fit <- lm(M1_Predicted ~ MODISobservedOnset), col='red')
legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
dev.off()
######################
#Scatterplot (Compare Model 2 predicted with LULC Onset VS MODIS observed Date of onset of greenness !)
#mydata$M2fitted <- model2$fitted.values
#MODISobservedOnset <- mydata$MODIS
#M2_Predicted <- mydata$M2fitted

#all = c(M2_Predicted, MODISobservedOnset)
#range = c(min(all), max(all))
#plot(M2_Predicted, MODISobservedOnset, xlim=range, ylim=range, pch=1, col="blue",cex = .3)
#abline(fit <- lm(M2_Predicted ~ MODISobservedOnset), col='red')
#legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))

##Making some maps in 2D space, raster layers

#df.month.tx <- subset(df, month == 4 | month == 5 & State == "OH" ) 

###############Scatter plots below are just GDDminus5 Vs MODIS Observed Scatter plots #######

#mytitle <- paste0("GDDonlyDeciduous_",k);mytitle
#jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
#        all = c(DecidForest$GDDMinus5, DecidForest$MODIS)
#        range = c(min(all), max(all))
#        plot(DecidForest$GDDMinus5, DecidForest$MODIS, xlim=range, ylim=range, ylab = "MODIS", xlab = "GDDminus5", main = mytitle, pch=1, col="blue",cex = .3)
#        abline(fit <- lm(DecidForest$MODIS ~ DecidForest$GDDMinus5), col='red')
#        legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
#dev.off()

#mytitle <- paste0("GDDonlyEvergreen_",k);mytitle
#jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
#        all2 = c(EvrgnForest$GDDMinus5, EvrgnForest$MODIS)
#        range = c(min(all2), max(all2))
#        plot(EvrgnForest$GDDMinus5, EvrgnForest$MODIS, xlim=range, ylim=range, ylab = "MODIS", xlab = "GDDminus5", main = mytitle, pch=1, col="blue",cex = .3)
#        abline(fit <- lm(EvrgnForest$MODIS ~ EvrgnForest$GDDMinus5), col='red')
#        legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
#dev.off()

#mytitle <- paste0("GDDonlyMixed_",k);mytitle
#jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
#        all3 = c(MixedForest$GDDMinus5, MixedForest$MODIS)
#        range = c(min(all3), max(all3)) 
#        plot(MixedForest$GDDMinus5, MixedForest$MODIS, xlim=range, ylim=range, ylab = "MODIS", xlab = "GDDminus5", main = mytitle, pch=1, col="blue",cex = .3)
#        abline(fit <- lm(MixedForest$MODIS ~ MixedForest$GDDMinus5), col='red')
#        legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
#dev.off()

##################
#plot(EvrgnForest$GDDMinus5, EvrgnForest$MODIS, ylab = "MODIS", xlab = "GDDminus5", pch=1, col="blue",cex = .3)
#abline(fit <- lm(EvrgnForest$MODIS ~ EvrgnForest$GDDMinus5), col='red')
#legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
###############Scatter plots below are Bayesian Model Predicted Mean Vs MODIS Observed Scatter plots #######
#mytitle <- paste0("Bays.Deciduous_",k);mytitle
#jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
#        all4 = c(DecidForest$Post_Mean, DecidForest$MODIS)
#        range = c(min(all4), max(all4)) 
#        plot(DecidForest$Post_Mean, DecidForest$MODIS, xlim=range, ylim=range, ylab = "MODIS", xlab = "Post_Mean", main = mytitle, pch=1, col="blue",cex = .3)
#        abline(fit <- lm(DecidForest$MODIS ~ DecidForest$Post_Mean), col='red')
#        legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
#dev.off()

#mytitle <- paste0("Bays.Evergreen_",k);mytitle
#jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
#        all5 = c(EvrgnForest$Post_Mean, EvrgnForest$MODIS)
#        range = c(min(all5), max(all5)) 
#        plot(EvrgnForest$Post_Mean, EvrgnForest$MODIS, xlim=range, ylim=range, ylab = "MODIS", xlab = "Post_Mean", main = mytitle, pch=1, col="blue",cex = .3)
#        abline(fit <- lm(EvrgnForest$MODIS ~ EvrgnForest$Post_Mean), col='red')
#        legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
#dev.off()

#mytitle <- paste0("Bays.Mixed_",k);mytitle
#jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
#        all6 = c(MixedForest$Post_Mean, MixedForest$MODIS)
#        range = c(min(all6), max(all6))
#        plot(MixedForest$Post_Mean, MixedForest$MODIS, xlim=range, ylim=range, ylab = "MODIS", xlab = "Post_Mean", main = mytitle, pch=1, col="blue",cex = .3)
#        abline(fit <- lm(MixedForest$MODIS ~ MixedForest$Post_Mean), col='red')
#        legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=2)))
#dev.off()
##########################################


library(rgdal)
library(rgeos)
library(maptools)
library(gridExtra) 
library(lattice)
library(maps)
library(fields)
library(raster)

head(mydata)
mydata$Lon <- as.numeric(mydata$Lon)
mydata$Lat <- as.numeric(mydata$Lat)
mydata.sp <- mydata
coordinates(mydata.sp) <- ~Lon+Lat
head(mydata)
mydata.sp

library(RColorBrewer)
#Lets plot some 2D raster layers 
#myraster <- raster(, xmn=-99.91305542, xmx=-67.25476074, ymn=25.15625, ymx=52.59024811, resolution=1/4) 
#myraster
#rm(r.count)
#Model_POST_SD <- rasterize(mydata[, c('Lon', 'Lat')], myraster, mydata[, 'Post_SD'], fun=mean)
#mytitle <- paste0("PosteriorSD_",k);mytitle
#jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
#r.range <- c(2,8)
#plot(Model_POST_SD, 
     #breaks = c(0 to 300 in increment of 5)
#     breaks = seq(2, 8, by=.5),  #define start & end, and interval
#     col = rainbow(12), #col = terrain.colors (13),#col = topo.colors(60),#col=rainbow(60),
#     main=mytitle, 
     #axis.args=list(at=seq(r.range[1], r.range[2], 30),labels=seq(r.range[1], r.range[2], 30),cex.axis=1))
#     axis.args=list(at=seq(r.range[1], r.range[2], 1),labels=seq(r.range[1], r.range[2], 1),cex.axis=1))
#map("state", add=T)
#dev.off()

myraster <- raster(, xmn=-99.91306, xmx=-67.16306, ymn=25.09025, ymx=52.59025, resolution=1/4) 
myraster
rm(r.count)
GDDpred <- rasterize(mydata[, c('Lon', 'Lat')], myraster, mydata[, 'GDDMinus5'], fun=mean)
r.range <- c(59,269)
GDDpred[GDDpred > 269] <- 269
GDDpred[GDDpred < 59] <- 59
mytitle <- paste0("GDDminus5_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
plot(GDDpred, 
     #breaks = c(0 to 300 in increment of 5)
     breaks = seq(59, 269, by=15),  #define start & end, and interval
     col = rainbow(14), #col = terrain.colors (13),#col = topo.colors(60),#col=rainbow(60),
     main=mytitle, 
     #axis.args=list(at=seq(r.range[1], r.range[2], 30),labels=seq(r.range[1], r.range[2], 30),cex.axis=1))
     axis.args=list(at=seq(r.range[1], r.range[2], 30),labels=seq(r.range[1], r.range[2], 30),cex.axis=1))
map("state", add=T)
dev.off()

myraster <- raster(, xmn=-99.91306, xmx=-67.16306, ymn=25.09025, ymx=52.59025, resolution=1/4) 
myraster
rm(r.count)
postmean <- rasterize(mydata[, c('Lon', 'Lat')], myraster, mydata[, 'Post_Mean'], fun=mean)
r.range <- c(59,269)
postmean[postmean > 269] <- 269
postmean[postmean < 59] <- 59
mytitle <- paste0("PostMean_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
plot(postmean, 
     #breaks = c(0 to 300 in increment of 5)
     breaks = seq(59, 269, by=15),  #define start & end, and interval
     col = rainbow(14), #col = terrain.colors (13),#col = topo.colors(60),#col=rainbow(60),
     main=mytitle, 
     #axis.args=list(at=seq(r.range[1], r.range[2], 30),labels=seq(r.range[1], r.range[2], 30),cex.axis=1))
     axis.args=list(at=seq(r.range[1], r.range[2], 30),labels=seq(r.range[1], r.range[2], 30),cex.axis=1))
map("state", add=T)
dev.off()

myraster <- raster(, xmn=-99.91306, xmx=-67.16306, ymn=25.09025, ymx=52.59025, resolution=1/4) 
myraster
rm(r.count)
modisOBS <- rasterize(mydata[, c('Lon', 'Lat')], myraster, mydata[, 'MODIS'], fun=mean)
r.range <- c(59,269)
modisOBS[modisOBS > 269] <- 269
modisOBS[modisOBS < 59] <- 59
mytitle <- paste0("MODIS_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
plot(modisOBS, 
     #breaks = c(0 to 300 in increment of 5)
     breaks = seq(59, 269, by=15),  #define start & end, and interval
     col = rainbow(14), #col = terrain.colors (13),#col = topo.colors(60),#col=rainbow(60),
     main=mytitle, 
     #axis.args=list(at=seq(r.range[1], r.range[2], 30),labels=seq(r.range[1], r.range[2], 30),cex.axis=1))
     axis.args=list(at=seq(r.range[1], r.range[2], 30),labels=seq(r.range[1], r.range[2], 30),cex.axis=1))
map("state", add=T)
dev.off()

library(RColorBrewer)
myraster <- raster(, xmn=-99.91306, xmx=-67.16306, ymn=25.09025, ymx=52.59025, resolution=1/4) 
myraster
ResidualOne <- rasterize(mydata[, c('Lon', 'Lat')], myraster, mydata[, 'residual1'], fun=mean)
ResidualOne[ResidualOne > 30] <- 30
ResidualOne[ResidualOne < -30] <- -30
#r.range <- c(minValue(ResidualOne), maxValue(ResidualOne))
r.range <-c(-30,30)
#cuts = c(minValue(ResidualOne), -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, maxValue(ResidualOne)) #breaks = seq(-15, 15, by=3)
cuts = c(-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30) #breaks = seq(-15, 15, by=3)
#pal <- colorRamp(c("red", "blue"))
pal <- colorRampPalette(c("cyan","purple","blue","red","pink","tan"))
#pal <- colorRampPalette(c("cyan","yellow","blue","red","green","orange","magenta","pink","turquoise","tan","purple","maroon"))
pal(12)
mytitle <- paste0("ResidualM1_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
plot(ResidualOne,
     breaks = cuts,
     col = pal(12),
     main=mytitle, 
     axis.args=list(at=seq(r.range[1], r.range[2], 5),labels=seq(r.range[1], r.range[2], 5),cex.axis=1,legend.shrink=.5))
map("state", add=T)
sum(mydata$residual1)
dev.off()

library(RColorBrewer)
myraster <- raster(, xmn=-99.91306, xmx=-67.16306, ymn=25.09025, ymx=52.59025, resolution=1/4) 
myraster
mydata$PostMeanMinusMODIS <- round(mydata$PostMeanMinusMODIS, digits = 0)
max(mydata$PostMeanMinusMODIS)
min(mydata$PostMeanMinusMODIS)
PostMinusMODIS <- rasterize(mydata[, c('Lon', 'Lat')], myraster, mydata[, 'PostMeanMinusMODIS'], fun=mean)
PostMinusMODIS[PostMinusMODIS > 14] <- 14
PostMinusMODIS[PostMinusMODIS < -14] <- -14
#r.range <- c(minValue(PostMinusMODIS), maxValue(PostMinusMODIS))
r.range <-c(-14,14)
#cuts = c(minValue(PostMinusMODIS), -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, maxValue(PostMinusMODIS)) #breaks = seq(-15, 15, by=3)
cuts = c(-14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14) #breaks = seq(-15, 15, by=3)
#pal <- colorRamp(c("red", "blue"))
pal <- colorRampPalette(c("cyan","purple","blue","red","pink","tan"))
#pal <- colorRampPalette(c("cyan","yellow","blue","red","green","orange","magenta","pink","turquoise","tan","purple","maroon"))
pal(14)
mytitle <- paste0("PostMeanMinMODS_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
plot(PostMinusMODIS,
     breaks = cuts,
     col = pal(14),
     main=mytitle, 
     axis.args=list(at=seq(r.range[1], r.range[2], 2),labels=seq(r.range[1], r.range[2], 2),cex.axis=1,legend.shrink=.5))
map("state", add=T)
sum(mydata$residual1)
dev.off()

mytitle <- paste0("HistPostMeanMinusMODIS_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
hist(mydata$PostMeanMinusMODIS, 
     main = mytitle,
     xlab = "Difference in Days",
     border = "red",
     col = "lightblue",
     xlim = c(-110, 80),
     las = 1,
     breaks = 40)
dev.off()


mytitle <- paste0("HistResidualM1_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
hist(mydata$residual1, 
     main = mytitle,
     xlab = "Difference in Days",
     border = "red",
     col = "lightblue",
     xlim = c(-110, 80),
     las = 1,
     breaks = 40)
dev.off()


mytitle <- paste0("HistPostSD_",k);mytitle
jpeg(paste0(mytitle,".jpg"),width =4, height =4, unit = "in",res = 300)
hist(mydata$Post_SD, 
     main = mytitle,
     xlab = "Posterior SD",
     border = "red",
     col = "lightblue",
     xlim = c(0, 8),
     las = 1,
     breaks = 16)
dev.off()

}


