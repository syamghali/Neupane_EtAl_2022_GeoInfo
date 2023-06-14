download_noaa_files_s3 <- function(siteID, date, cycle, local_directory){
  
  Sys.setenv("AWS_DEFAULT_REGION" = "data",
             "AWS_S3_ENDPOINT" = "ecoforecast.org")
  
  object <- aws.s3::get_bucket("drivers", prefix=paste0("noaa/NOAAGEFS_1hr/",siteID,"/",date,"/",cycle))
  
  for(i in 1:length(object)){
    aws.s3::save_object(object[[i]], bucket = "drivers", file = file.path(local_directory, object[[i]]$Key))
  }
}

download_noaa_files_s3(siteID = "BART", date = "2021-05-18", cycle = "00", local_directory <- "/Volumes/GoogleDrive/My Drive/VaughnData/PhenoPrediction/")

#BART UKFS HARV DELA SCBI CLBJ STEI GRSM 
