library(tidyverse)
library(devtools)

if ("rcldf" %in% installed.packages()) {
  print("rcldf is installed")
} else{
  print("installing rcldf......")
  install_github("SimonGreenhill/rcldf", dependencies = TRUE)
}

library(rcldf)

# this script fetches the asjp data from zenodo and saves it as an rdf

asjp = rcldf::cldf(mdpath = "https://zenodo.org/records/7079637/files/lexibank/asjp-v20.zip")
write_rds(asjp, "data/asjp_database.rds")


