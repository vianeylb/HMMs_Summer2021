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

setwd("C:/Jessica/UofT Y4/Research/Coding")
setwd("C:/Jessica/UofT Y4/Research/Coding/Field Data")

line_plot_static <- function(data, filename) {
  plotx <- ggplot(data, aes(x = Time, y = X_static)) +
    geom_line(colour = "dark red") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ploty <- ggplot(data, aes(x = Time, y = Y_static)) +
    geom_line(colour = "navy") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  plotz <- ggplot(data, aes(x = Time, y = Z_static)) +
    geom_line(colour = "dark green") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  plot <- grid.arrange(plotx, ploty, plotz, nrow = 3)
  ggsave(paste(filename, "png", sep = "."), plot)
}

hist_plot_static <- function(data, filename) {
  plotx <- ggplot(data, aes(X_static)) +
    geom_histogram(colour = "dark red", fill = "salmon") +
    theme_minimal()
  ploty <- ggplot(data, aes(Y_static)) +
    geom_histogram(colour = "navy", fill = "light blue") +
    theme_minimal()
  plotz <- ggplot(data, aes(Z_static)) +
    geom_histogram(colour = "dark green", fill = "light green") +
    theme_minimal()
  plot <- grid.arrange(plotx, ploty, plotz, nrow = 3)
  ggsave(paste(filename, "png", sep = "."), plot)
}

acf_plot_static <- function(data, filename) {
  plotx <- ggacf(data$X_static) +
    theme_minimal()
  ploty <- ggacf(data$Y_static) +
    theme_minimal()
  plotz <- ggacf(data$Z_static) +
    theme_minimal()
  plot <- grid.arrange(plotx, ploty, plotz, nrow = 3)
  ggsave(paste(filename, "png", sep = "."), plot)
}

pacf_plot_static <- function(data, filename) {
  plotx <- ggpacf(data$X_static) +
    theme_minimal()
  ploty <- ggpacf(data$Y_static) +
    theme_minimal()
  plotz <- ggpacf(data$Z_static) +
    theme_minimal()
  plot <- grid.arrange(plotx, ploty, plotz, nrow = 3)
  ggsave(paste(filename, "png", sep = "."), plot)
}

get_plots <- function(indicies){
  n <- length(indicies)
  for (i in 1:n){
    index <- indicies[i]
    filename <- paste("Track", index, "_Static_25Hz.csv", sep = "")
    data <- read.csv(filename)
    static_data <- data %>% select(X_static, Y_static, Z_static)
    line_plot_static(data, "Track03_static_line_plot")
    hist_plot_static(data, "Track03_static_histogram")
    acf_plot_static(data, "Track03_static_acf")
    pacf_plot_static(data, "Track03_static_pacf")
    pairs_plot(static_data, "Track03_static_correlation_static.png")
    
    line_plot(data, paste("Track", index, "_dynamic_line_plot", sep = ""))
    hist_plot(data, paste("Track", index, "_dynamic_histogram", sep = "")) 
    acf_plot(data, paste("Track", index, "_dynamic_acf", sep = ""))
    pacf_plot(data, paste("Track", index, "_dynamic_pacf", sep = ""))
    pairs_plot(dynamic_data, paste("Track", index, "_dynamic_correlation.png", sep = "")) 
  }
}

indicies <- c("01", "03", "05", "06", "07", "12")

