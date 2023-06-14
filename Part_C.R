rm(list=ls())
library(magrittr)
library(spmeshed)
library(tidyverse)


process <- F
analysis_name <- "svi_2021" %>% glue::glue()

if(process){
  ## LOAD MODIS with GDD
  dflist <- list()
  #if(grepl("newGDD", analysis_name)){
    yearlist <- 2001:2019
  #} else {
    #yearlist <- 2002:2019
  #}
  
  for(yn in 1:length(yearlist)){
    y <- yearlist[yn]
    cat(y, "\n")
    df <- read.csv(glue::glue("data/GDD/AccGDD_{y}.csv")) %>% 
      #read.csv(glue::glue("data/PredJustFromPrecedingYrEastUS_{y}.csv")) %>% 
      mutate(year=y)
    if("dd2200" %in% colnames(df)){
      df %<>% rename(gdd2200 = dd2200)
    }
    dflist[[yn]] <- df %>% 
      filter(MODIS > 0, Lon > -100)#, Lat < 48, Lon < -75)
  }
  
  dfbig <- bind_rows(dflist) 
  
  #if(grepl("newGDD", analysis_name)){
    gdd_names <- colnames(dfbig)[colnames(dfbig) %>% grepl("gdd", .)][-1:-3]
    drop_gdd <- colnames(dfbig)[colnames(dfbig) %>% grepl("gdd", .)][1:3]
    keepcols <- colnames(dfbig) %>% setdiff(drop_gdd)
    dfbig %<>% dplyr::select(!!!syms(keepcols))
  #} else {
  #  gdd_names <- "GDD"
  #}
  
  
  
  covertype <- read.csv("data/covertypes_processed.csv", stringsAsFactor=F) %>% unique()
  
  if(grepl("_forest", analysis_name)){
    covertype %<>% filter(DominantCover %in% c("EvergreenForest", "DeciduousForest", "MixedForest"))
  }
  
  covertype %<>% 
    mutate(Value=1) %>% spread(DominantCover, Value, fill=0)
  colnames(covertype)[-1:-2] <- paste0("dc_", colnames(covertype)[-1:-2])
  
  if(T){
    # keep only locations with GDD for each year
    coords_common <- dflist[[1]] %>% 
      filter(complete.cases(!!!syms(gdd_names))) %>% 
      dplyr::select(Lon, Lat) %>% unique()
    
    for(i in 2:length(dflist)){
      coords_here <- dflist[[i]] %>% 
        filter(complete.cases(!!!syms(gdd_names))) %>% 
        dplyr::select(Lon, Lat) %>% unique()
      
      coords_common %<>% inner_join(coords_here, by=c("Lon"="Lon", "Lat"="Lat"))
      cat(nrow(coords_common), "\n")
    }
    
    coords_common_by_year <- list()
    for(i in 1:length(dflist)){
      year <- dflist[[i]]$year[1]
      coords_common_by_year[[i]] <- coords_common %>% mutate(year=year)
    }
    coords_common <- bind_rows(coords_common_by_year)
    
    
    dfbig <- coords_common %>% left_join(dfbig)
  }
  
  dfbig %<>% mutate(MODIS_original = MODIS)
  dfbig[dfbig$year %in% c(2017, 2018, 2019), "MODIS"] <- NA
  
  
  covertype_big <- dfbig %>% dplyr::select(Lon, Lat, year) %>% left_join(covertype)
  covertype_big[is.na(covertype_big)] <- 0
  
  dfbig %<>% left_join(covertype_big)
  
  which_available <- which(!is.na(dfbig$MODIS))
  
  set.seed(20200923)
  test_set_size <- 50000
  test_set <- which_available %>% sample(test_set_size, replace=F)
  
  df_validate <- dfbig[test_set,] %>% arrange(year, Lon, Lat)
  dftarget <- dfbig
  dftarget[test_set, "MODIS"] <- NA # will be predicted
  
  dftarget %<>% rename(Var1=Lon, Var2=Lat, Var3=year) %>% 
    arrange(Var1, Var2, Var3)
  
  # preprocessing: remove larger areas of NAs (lakes, ocean...)
  n_years <- dftarget$Var3 %>% unique() %>% length()
  
  ####################
  
  # plot
  (plottest <- ggplot(dftarget %>% filter(Var3==2013), aes(x=Var1, y=Var2, fill=MODIS)) +
    geom_raster() + scale_fill_viridis_c() +
    theme_minimal() +
    theme(legend.position="none"))
  #ggsave(plot=plottest, filename= "temp.png", width=9, height=6)
  
  ####################
  save.image(file="temp/Greenness_{analysis_name}_processed.RData" %>% glue::glue())
} else {
  load("temp/Greenness_{analysis_name}_processed.RData" %>% glue::glue())
}


uLon <- dftarget$Var1 %>% unique() %>% sort()
uLat <- dftarget$Var2 %>% unique() %>% sort()
subLon <- uLon[seq(1, length(uLon), 1)]
subLat <- uLat[seq(1, length(uLat), 1)]

dftarget_sub <- dftarget %>% filter(Var1 %in% subLon, Var2 %in% subLat)
dftarget_sub %<>% arrange(Var1, Var2, Var3)

coords <- dftarget_sub %>% dplyr::select(Var1, Var2, Var3) %>% as.matrix()
y <- dftarget_sub %>% pull(MODIS)
ymean <- mean(y, na.rm=T)

time_trend <- dftarget_sub %>% dplyr::select(Var3) %>% as.matrix()

Xgdd <- dftarget_sub %>% 
  dplyr::select(contains("gdd")) #, contains("forest"))
X <- time_trend %>% matrix(ncol=1) #cbind(time_trend, Xgdd) #XH)
colnames(X)[1] <- c("time_trend") 

coords %>% apply(2, function(x) unique(x) %>% length())

axis_partition <- c(72, 60, 19)
nrow(coords)/prod(axis_partition)

mcmc_keep <- 1000
mcmc_burn <- 10000
mcmc_thin <- 10

forest <- ifelse(sum(grepl("Forest", colnames(X)))>0, "forest", "noforest")
hasgdd <- ifelse(sum(grepl("gdd", colnames(X)))>0, "gdd", "nogdd")

set.seed(1)
mesh_total_time <- system.time({
  meshout1 <- meshed(matrix(y, ncol=1), #matrix(y-ymean, ncol=1), 
                    family=c("gaussian"), X, coords, k = 1,
                    axis_partition=axis_partition,
                    n_samples = mcmc_keep, n_burn = mcmc_burn, n_thin = mcmc_thin, 
                    n_threads = 30,
                    prior = list(btmlim= .01, toplim=1e3, phi=c(.01, 1e6), nu=c(.5, .5)),
                    settings = list(adapting=T, forced_grid=F, cache=T, saving=T, ps=T),
                    print_every=1,
                    debug=list(sample_beta=T, sample_tausq=T, 
                               sample_theta=T, sample_w=T, sample_lambda=T,
                               verbose=F, debug=F)
  )})

beta_mcmc <- meshout1$beta_mcmc
tausq_mcmc <- meshout1$tausq_mcmc
theta_mcmc <- meshout1$theta_mcmc
lambda_mcmc <- meshout1$lambda_mcmc

yhat <- meshout1$yhat_mcmc %>% meshgp::list_mean() %>% add(ymean)
ylow95 <- meshout1$yhat_mcmc %>% meshgp::list_qtile(0.025) %>% add(ymean)
yhi95 <- meshout1$yhat_mcmc %>% meshgp::list_qtile(0.975) %>% add(ymean)
ylow90 <- meshout1$yhat_mcmc %>% meshgp::list_qtile(0.05) %>% add(ymean)
yhi90 <- meshout1$yhat_mcmc %>% meshgp::list_qtile(0.95) %>% add(ymean)
ylow75 <- meshout1$yhat_mcmc %>% meshgp::list_qtile(.25/2) %>% add(ymean)
yhi75 <- meshout1$yhat_mcmc %>% meshgp::list_qtile(1-.25/2) %>% add(ymean)

what <- meshout1$w_mcmc %>% meshgp::list_mean()
wlow95 <- meshout1$w_mcmc %>% meshgp::list_qtile(0.025)
whi95 <- meshout1$w_mcmc %>% meshgp::list_qtile(0.975)
wlow90 <- meshout1$w_mcmc %>% meshgp::list_qtile(0.05)
whi90 <- meshout1$w_mcmc %>% meshgp::list_qtile(.95)
wlow75 <- meshout1$w_mcmc %>% meshgp::list_qtile(.25/2)
whi75 <- meshout1$w_mcmc %>% meshgp::list_qtile(1 - .25/2)

q <- 1 #ncol(Zs)
outdf <- meshout1 %>% with(cbind(coordsdata %>% dplyr::select(-forced_grid), 
                                 ylow95, ylow90, ylow75, yhat, yhi75, yhi90, yhi95, 
                                 wlow95, wlow90, wlow75, what, whi75, whi90, whi95)) %>% as.data.frame()
colnames(outdf)[-1:-3] <- c("y_low95", "y_low90", "y_low75", "y_hat", "y_hi75", "y_hi90", "y_hi95",
                            paste0("w", 1:q, "_low95"), paste0("w", 1:q, "_low90"), paste0("w", 1:q, "_low75"), 
                            paste0("w", 1:q, "_hat"), paste0("w", 1:q, "_hi75"), paste0("w", 1:q, "_hi90"), paste0("w", 1:q, "_hi95")
)

results <- dftarget_sub %>% left_join(outdf)

temp_plot <- results %>% filter(Var3==2019) %>% 
  ggplot(aes(Var1, Var2, fill=abs(y_hat))) + geom_raster()
ggsave(temp_plot, filename="temp1.png")

#others <- list(H = "H", wavelet_code='wavethresh::GenW(n=16, filter.number=10, family="DaubLeAsymm", bc="periodic")')
others <- list(description="X includes time, {hasgdd}, {forest}" %>% glue::glue())

save(file="out/Greenness_{analysis_name}_{hasgdd}_{forest}_beta.RData" %>% glue::glue(), 
     list=c("outdf", "mesh_total_time", "others", "forest", "hasgdd", "analysis_name",
            "X", "ymean", "axis_partition", "beta_mcmc", "tausq_mcmc", "theta_mcmc", "lambda_mcmc"))