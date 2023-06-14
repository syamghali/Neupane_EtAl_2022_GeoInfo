# PhenoPhunctions

# Load required libraries
require(tidyverse)

# Set constant GDD value thresholds per site for phase transitions
Phase2Trans <- c("BART"=380, "SCBI"=860, "HARV"=800, "GRSM"=1010, "STEI"=415,
                 "CLBJ"=1100, "DELA"=1275, "UKFS"=1415)
Phase3Trans <- c("BART"=156, "SCBI"=146, "HARV"=150, "GRSM"=112, "STEI"=161,
                 "CLBJ"=101, "DELA"=93, "UKFS"=127)

# Set constant value for GDD already accumulated in 2021
GDD.cum.2021 <- c("BART"=28.3, "CLBJ"=273, "DELA"=372,
                  "GRSM"=226, "HARV"=55, "SCBI"=108,
                  "STEI"=15, "UKFS"=180)

