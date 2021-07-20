library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(chron)
library(lemon)
library(GGally)
library(RColorBrewer)

setwd("C:/Jessica/UofT Y4/Research/Coding/Lab Data")

get_dynamic_sd_data <- function(old_filename, filename) {
  old_data <- read.csv(file = old_filename)
  
  sds <- tapply(old_data$ODBA, old_data$Time, sd)
  data <- data_frame(ODBA = sds, Time = dimnames(sds)[[1]])
  data$Time <- chron(times = data$Time)
  data$X_dynamic <- tapply(old_data$X_dynamic, old_data$Time, sd)
  data$Y_dynamic <- tapply(old_data$Y_dynamic, old_data$Time, sd)
  data$Z_dynamic <- tapply(old_data$Z_dynamic, old_data$Time, sd)
  data$Behavior <- tapply(old_data$Behavior, old_data$Time, unique)
  data$Prey <- tapply(old_data$Prey, old_data$Time, unique)
  
  write.csv(data, file = filename, row.names = FALSE)
}

get_static_sd_data <- function(old_filename, filename) {
  old_data <- read.csv(file = old_filename)
  
  sds <- tapply(old_data$X_static, old_data$Time, sd)
  data <- data_frame(X_static = sds, Time = dimnames(sds)[[1]])
  data$Time <- chron(times = data$Time)
  data$Y_static <- tapply(old_data$Y_static, old_data$Time, sd)
  data$Z_static <- tapply(old_data$Z_static, old_data$Time, sd)
  data$Behavior <- tapply(old_data$Behavior, old_data$Time, unique)
  data$Prey <- tapply(old_data$Prey, old_data$Time, unique)
  
  write.csv(data, file = filename, row.names = FALSE)
}

get_all_sd_data <- function(names) {
  n <- length(names)
  for (i in 1:n){
    name <- names[i]
    
    old_dynamic_filename <- paste("Custom", name, "dynamic.csv", sep = "_")
    dynamic_filename <- paste("SD", name, "dynamic.csv", sep = "_")
    get_dynamic_sd_data(old_dynamic_filename, dynamic_filename)
    
    old_static_filename <- paste("Custom", name, "static.csv", sep = "_")
    static_filename <- paste("SD", name, "static.csv", sep = "_")
    get_static_sd_data(old_static_filename, static_filename)
  }
}

get_dynamic_sd_plots <- function(names){
  n <- length(names)
  for (i in 1:n){
    name <- names[i]
    filename <- paste("SD", name, "dynamic.csv", sep = "_")
    data <- read.csv(filename)
    labelled_data <- data %>% filter(!is.na(Behavior))
    dynamic_data <- data %>% select(X_dynamic, Y_dynamic, Z_dynamic)
    
    line_plot(data, paste(name, "sd", "line_plot", sep = "_"))
    hist_plot(data, paste(name, "sd", "histogram", sep = "_"))
    acf_plot(data, paste(name, "sd", "acf", sep = "_"))
    pacf_plot(data, paste(name, "sd", "pacf", sep = "_"))
    behavior_plot(data, paste(name, "sd", "plot_behavior", sep = "_"))
    behavior_plot(labelled_data, paste(name, "sd", "plot_behavior_filtered", sep = "_"))
    filtered_hist(labelled_data, paste(name, "sd", "histogram_filtered", sep = "_"))
    behavior_hist(labelled_data, paste(name, "sd", "histogram_behavior", sep = "_"))
    pairs_plot(dynamic_data, paste(name, "sd", "correlation_dynamic.png", sep = "_"))
    behavior_pairs_plot_dynamic(data, paste(name, "sd", "dynamic", sep = "_"))
  }
}

get_static_sd_plots <- function(names){
  n <- length(names)
  for (i in 1:n){
    name <- names[i]
    filename <- paste("SD", name, "static.csv", sep = "_")
    data <- read.csv(filename)
    labelled_data <- data %>% filter(!is.na(Behavior))
    static_data <- data %>% select(X_static, Y_static, Z_static)
    
    line_plot_static(data, paste(name, "sd", "line_plot", sep = "_"))
    hist_plot_static(data, paste(name, "sd", "histogram", sep = "_"))
    acf_plot_static(data, paste(name, "sd", "acf", sep = "_"))
    pacf_plot_static(data, paste(name, "sd", "pacf", sep = "_"))
    behavior_plot_static(data, paste(name, "sd", "plot_behavior", sep = "_"))
    behavior_plot_static(labelled_data, paste(name, "sd", "plot_behavior_filtered", sep = "_"))
    filtered_hist_static(labelled_data, paste(name, "sd", "histogram_filtered", sep = "_"))
    behavior_hist_static(labelled_data, paste(name, "sd", "histogram_behavior", sep = "_"))
    pairs_plot(static_data, paste(name, "sd", "corxrelation_static.png", sep = "_"))
    behavior_pairs_plot_static(data, paste(name, "sd", "static", sep = "_"))
  }
}

names <- c("BigDaddy_3Apr17",
           "BigDaddy_20Mar17",
           "BigGuy_15Feb18",
           "Eliza_7Sept17",
           "Eliza_20Sept17",
           "Lady_27Mar17")
get_all_sd_data(names)
get_dynamic_sd_plots(names)
get_static_sd_plots(names)

