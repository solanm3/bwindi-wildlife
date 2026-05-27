# main.R — run this file to execute the full analysis

library(tidyverse)
library(lubridate)
library(iNEXT)
library(vegan)
library(overlap)
library(effsize)
library(glue)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(performance)
library(cowplot)
library(ggrepel)
library(ncdf4)
library(ggbeeswarm)


setwd()

# Base path for analysis output
DATA_DIR <- " "

source("bwindi_theme.R") # Theme for plots
source("load_data.R")  
source("H1_richness.R")
source("H2_composition.R")
source("H3_temporal.R")

sessionInfo()
