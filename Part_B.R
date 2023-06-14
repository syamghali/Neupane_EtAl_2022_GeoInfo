rm(list=ls())
gc()
library(colorout)
library(tidyverse)
library(magrittr)
#library(meshgp)

covertype <- read.csv("data/covertypes_processed.csv")

#analysis_name <- "svi_2021"
analysis_name <- "svi_2021_forest"

load("temp/Greenness_svi_2021_processed.RData" %>% glue::glue())
analysis_name <- "svi_2021_forest"
load("out/Greenness_{analysis_name}.RData" %>% glue::glue())

results <- outdf %>% 
  rename(Lon=Var1, Lat=Var2, year=Var3,
         ModelPrediction = y_hat,
         ModelPrediction_low = y_low95,
         ModelPrediction_high = y_hi95,
         SVI = w1_hat) 
         
## check results out of sample
validation <- df_validate %>% left_join(results)
validation %>% with(abs(MODIS-ModelPrediction)) %>% mean() 
# 5.44 days average absolute error in days
validation %>% with((MODIS-ModelPrediction)^2) %>% mean() %>% sqrt()
# 10.05 days root mean square error in days
covg95 <- validation %>% with((ModelPrediction_low < MODIS) * (MODIS < ModelPrediction_high)) %>% mean() %>% multiply_by(100)
validation %>% with(abs(ModelPrediction_high-ModelPrediction_low)) %>% mean()
cat("Nominal 95% coverage? We got {covg95}%" %>% glue::glue(), "\n")


# rebuild original data
results %<>% right_join(dfbig) %>%
  rename(MODIS_observed = MODIS) %>%
  left_join(covertype)

data2017 <- results %>% filter(year == 2017) 
data2017 %>%
  with(abs(ModelPrediction - MODIS_original)) %>% mean() 
data2017 %>%
  with((ModelPrediction - MODIS_original)^2) %>% mean() %>% sqrt()
data2017 %>% with(abs(ModelPrediction_low - ModelPrediction_high)) %>% summary()
data2017 %>% with((ModelPrediction_low < MODIS_original) * (MODIS_original < ModelPrediction_high)) %>% mean()


y <- 2018

dfbig <- read.csv(glue::glue("data/GDD/AccGDD_{y}.csv")) %>%
  mutate(year=y) %>% dplyr::rename(MODIS_original = MODIS)

data2018 <- results %>% filter(year == 2018) %>% 
  dplyr::select(-MODIS_observed, -MODIS_original) %>% 
  left_join(dfbig) %>%
  filter(complete.cases(MODIS_original))
data2018 %>%
  with(abs(ModelPrediction - MODIS_original)) %>% mean() 
data2018 %>%
  with((ModelPrediction - MODIS_original)^2) %>% mean() %>% sqrt()
data2018 %>% with(abs(ModelPrediction_low - ModelPrediction_high)) %>% summary()
data2018 %>% with((ModelPrediction_low < MODIS_original) * (MODIS_original < ModelPrediction_high)) %>% mean()

plottest <- ggplot(data2018 %>% 
                    mutate(anyforest = dc_DeciduousForest + dc_EvergreenForest + dc_MixedForest) %>%
                    mutate(Mdiff = MODIS_original-ModelPrediction) %>%
                    mutate(Mdiff = ifelse(Mdiff > 10, 10, ifelse(Mdiff < -10, -10, Mdiff))) %>% 
                    filter(anyforest>0), 
                    aes(Lon, Lat, fill=Mdiff)) +
geom_raster() +
scale_fill_viridis_c() +
theme_minimal() +
labs(fill="MODIS - Prediction")
theme(legend.position="bottom") 
ggsave(plot=plottest, filename= "plottest2018.png")




# residuals
model_residuals <- results %>% with(MODIS_observed-ModelPrediction) 

save(list=c("results","lambda_mcmc", "beta_mcmc", "theta_mcmc", "tausq_mcmc", "X"), file="export/MeshGP_Export_{analysis_name}_{forest}_{hasgdd}.RData" %>% glue::glue())



# map predictions for next year
results %>% 
  filter(year==2004) %>%
  ggplot(aes(Lon, Lat, fill=w_hat)) + 
  geom_raster() + 
  scale_fill_viridis_c() + 
  theme_minimal()# + 
  #theme(legend.position="none")

# how does the prediction for 2018 compare to 2017 ?
results17 <- results %>% filter(year==2017) %>% dplyr::select(Lon, Lat, ModelPrediction)
results18 <- results %>% filter(year==2018) %>% dplyr::select(Lon, Lat, ModelPrediction)
results_nextyear <- results17 %>% rename(ModelPrediction17 = ModelPrediction) %>% 
  left_join(results18 %>% rename(ModelPrediction18 = ModelPrediction))

plotsave <- results_nextyear %>% ggplot(aes(Lon, Lat, fill=ModelPrediction18-ModelPrediction17)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme_minimal() 
ggsave(plotsave, file="plot/temp.png")


