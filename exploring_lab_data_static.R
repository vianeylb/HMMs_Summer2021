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

source("univariate_normal_hmm_functions.R")

setwd("C:/Jessica/UofT Y4/Research/Coding/Lab Data")

get_custom_data_static <- function(data_file, timestamps_file, filename) {
  data <- read.csv(file = data_file)
  timestamps <- read.csv(file = timestamps_file)
  
  # Alter and rename columns to be consistent
  if ("Date.and.Time" %in% colnames(data)) {
    data <- data %>%
      separate(col = Date.and.Time, into = c("Date", "Time"), sep = " ")
  }
  if ("Prey.Type.Notes" %in% colnames(timestamps)) {
    colnames(timestamps)[which(names(timestamps) ==
                                 "Prey.Type.Notes")] <- "Prey"
  }
  if ("Prey.Type" %in% colnames(timestamps)) {
    colnames(timestamps)[which(names(timestamps) ==
                                 "Prey.Type")] <- "Prey"
  }
  if ("Time.End" %in% colnames(timestamps)) {
    colnames(timestamps)[which(names(timestamps) ==
                                 "Time.End")] <- "Time.Stop"
  }
  
  # Add Behavior and Prey column to main data set, based on data in timestamps
  data$Time <- chron(times = data$Time)
  timestamps$Time.Start <- chron(times = timestamps$Time.Start)
  timestamps$Time.Stop <- chron(times = timestamps$Time.Stop)
  data$Behavior <- ifelse(timestamps$Time.Start <= data$Time &
                            timestamps$Time.Stop >= data$Time,
                          timestamps$Behavior, NA
  )
  data$Prey <- ifelse(timestamps$Time.Start <= data$Time &
                        timestamps$Time.Stop >= data$Time,
                      timestamps$Prey, NA
  )
  # Write to CSV
  write.csv(data, file = filename, row.names = FALSE)
}

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
  ggsave(paste(filename, "static.png", sep = "_"), plot)
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
  ggsave(paste(filename, "static.png", sep = "_"), plot)
}

acf_plot_static <- function(data, filename) {
  plotx <- ggacf(data$X_static) +
    theme_minimal()
  ploty <- ggacf(data$Y_static) +
    theme_minimal()
  plotz <- ggacf(data$Z_static) +
    theme_minimal()
  plot <- grid.arrange(plotx, ploty, plotz, nrow = 3)
  ggsave(paste(filename, "static.png", sep = "_"), plot)
}

pacf_plot_static <- function(data, filename) {
  plotx <- ggpacf(data$X_static) +
    theme_minimal()
  ploty <- ggpacf(data$Y_static) +
    theme_minimal()
  plotz <- ggpacf(data$Z_static) +
    theme_minimal()
  plot <- grid.arrange(plotx, ploty, plotz, nrow = 3)
  ggsave(paste(filename, "static.png", sep = "_"), plot)
}

behavior_plot_static <- function(data, filename) {
  plotx <- ggplot(data, aes(x = Time, y = X_static, colour = Behavior)) +
    geom_point(size = 0.5) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = 7)
    ) +
    scale_color_brewer(palette = "Set1")
  ploty <- ggplot(data, aes(x = Time, y = Y_static, colour = Behavior)) +
    geom_point(size = 0.5) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = 7)
    ) +
    scale_color_brewer(palette = "Set1")
  plotz <- ggplot(data, aes(x = Time, y = Z_static, colour = Behavior)) +
    geom_point(size = 0.5) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = 7)
    ) +
    scale_color_brewer(palette = "Set1")
  plot <- grid_arrange_shared_legend(plotx, ploty, plotz,
                                     nrow = 3, ncol = 1,
                                     position = "right")
  ggsave(paste(filename, "static.png", sep = "_"), plot)
}

filtered_hist_static <- function(data, filename) {
  plotx <- ggplot(data,
                  aes(x = X_static, colour = Behavior, fill = Behavior)) +
    geom_histogram() +
    facet_wrap(~Behavior, ncol = 2) +
    theme_minimal() +
    scale_fill_brewer(palette = "Pastel1") +
    scale_color_brewer(palette = "Set1")
  ggsave(paste(filename, "X_static.png", sep = "_"), plotx)
  ploty <- ggplot(data,
                  aes(x = Y_static, colour = Behavior, fill = Behavior)) +
    geom_histogram() +
    facet_wrap(~Behavior, ncol = 2) +
    theme_minimal() +
    scale_fill_brewer(palette = "Pastel1") +
    scale_color_brewer(palette = "Set1")
  ggsave(paste(filename, "Y_static.png", sep = "_"), ploty)
  plotz <- ggplot(data,
                  aes(x = Z_static, colour = Behavior, fill = Behavior)) +
    geom_histogram() +
    facet_wrap(~Behavior, ncol = 2) +
    theme_minimal() +
    scale_fill_brewer(palette = "Pastel1") +
    scale_color_brewer(palette = "Set1")
  ggsave(paste(filename, "Z_static.png", sep = "_"), plotz)
  
  behaviors <- unique(data$Behavior)
  for (i in seq_len(length(behaviors))) {
    behavior <- behaviors[i]
    subdata <- data %>% dplyr::filter(Behavior == behavior)
    
    plotx <- ggplot(subdata, aes(x = X_static)) +
      geom_histogram(colour = "dark red", fill = "salmon") +
      theme_minimal()
    ploty <- ggplot(subdata, aes(x = Y_static)) +
      geom_histogram(colour = "navy", fill = "light blue") +
      theme_minimal()
    plotz <- ggplot(subdata, aes(x = Z_static)) +
      geom_histogram(colour = "dark green", fill = "light green") +
      theme_minimal()
    plot <- grid.arrange(plotx, ploty, plotz, nrow = 3)
    ggsave(paste(filename, behavior, "static.png", sep = "_"), plot)
  }
}

behavior_hist_static <- function(data, filename) {
  # Create new column indicating each subinterval of behavior
  n <- length(data$Behavior)
  indicies <- c(1, which(data$Behavior != lag(data$Behavior)), n)
  m <- length(indicies)
  foo <- numeric(n)
  for (i in 1:(m - 1)) {
    foo[indicies[i]:indicies[i + 1]] <-
      rep(i, indicies[i + 1] - indicies[i] + 1)
  }
  data$BehaviorIndex <- as.character(foo)
  
  behaviors <- unique(data$Behavior)
  for (i in seq_len(length(behaviors))) {
    behavior <- behaviors[i]
    subdata <- data %>% filter(Behavior == behavior)
    
    plotx <- ggplot(data = subdata,
                    aes(x = X_static,
                        colour = BehaviorIndex,
                        fill = BehaviorIndex)) +
      geom_histogram() +
      facet_wrap(~BehaviorIndex, ncol = 2) +
      theme_minimal() +
      scale_fill_brewer(palette = "Pastel1") +
      scale_color_brewer(palette = "Set1")
    ggsave(paste(filename, behavior, "X_static.png", sep = "_"), plotx)
    
    ploty <- ggplot(data = subdata,
                    aes(x = Y_static,
                        colour = BehaviorIndex,
                        fill = BehaviorIndex)) +
      geom_histogram() +
      facet_wrap(~BehaviorIndex, ncol = 2) +
      theme_minimal() +
      scale_fill_brewer(palette = "Pastel1") +
      scale_color_brewer(palette = "Set1")
    ggsave(paste(filename, behavior, "Y_static.png", sep = "_"), ploty)
    
    plotz <- ggplot(data = subdata,
                    aes(x = Z_static,
                        colour = BehaviorIndex,
                        fill = BehaviorIndex)) +
      geom_histogram() +
      facet_wrap(~BehaviorIndex, ncol = 2) +
      theme_minimal() +
      scale_fill_brewer(palette = "Pastel1") +
      scale_color_brewer(palette = "Set1")
    ggsave(paste(filename, behavior, "Z_static.png", sep = "_"), plotz)
  }
}

pairs_plot <- function(data, filename) {
  plot <- ggpairs(data,
                  lower = list(continuous = wrap("smooth", alpha = 0.2, size = 0.1)),
                  diag = list(continuous = "bar")
  ) +
    theme_minimal()
  ggsave(filename, plot)
}

pseudo_residual_plot <- function(data, filename) {
  # Index plot of pseudo-residuals
  plot_index <- ggplot(data) +
    geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
    theme_minimal()
  # Histogram of pseudo-residuals
  plot_hist <- ggplot(data, aes(npsr)) +
    geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
    stat_function(fun = dnorm, colour = "red") +
    theme_minimal()
  # QQ plot of pseudo-residuals
  plot_qq <- ggplot(data, aes(sample = npsr)) +
    stat_qq() +
    stat_qq_line() +
    theme_minimal()
  # ACF of pseudo-residuals
  plot_acf <- ggacf(data$npsr) +
    theme_minimal()
  plot <- grid.arrange(plot_index, plot_hist, plot_qq, plot_acf,
                       nrow = 2, ncol = 2)
  ggsave(paste(filename, "png", sep = "."), plot)
}

data_file <- "Lady_27Mar17_Static_25Hz.csv"
timestamps_file <- "TimestampedData_Lady_27Mar17.csv"
filename <- "Custom_Lady_27Mar17_static.csv"
get_custom_data(data_file, timestamps_file, filename)
data <- read.csv(filename)
labelled_data <- data %>% filter(!is.na(Behavior))
static_data <- data %>% select(X_static, Y_static, Z_static)

line_plot_static(data, "Lady_27Mar17_line_plot")
hist_plot_static(data, "Lady_27Mar17_histogram")
acf_plot_static(data, "Lady_27Mar17_acf")
pacf_plot_static(data, "Lady_27Mar17_pacf")
behavior_plot_static(data, "Lady_27Mar17_plot_behavior")
behavior_plot_static(labelled_data, "Lady_27Mar17_plot_behavior_filtered")
filtered_hist_static(labelled_data, "Lady_27Mar17_histogram_filtered")
behavior_hist_static(labelled_data, "Lady_27Mar17_histogram_behavior")
pairs_plot(static_data, "Lady_27Mar17_correlation_static.png")