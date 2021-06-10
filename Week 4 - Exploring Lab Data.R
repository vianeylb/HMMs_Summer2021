library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(chron)
library(lemon)

get_custom_data <- function(data_file, timestamps_file, filename){
  data <- read.csv(file = data_file)
  timestamps <- read.csv(file = timestamps_file)
  
  #Alter and rename columns to be consistent
  if ("Date.and.Time" %in% colnames(data)) {
    data <- data %>% separate(col = Date.and.Time, into = c("Date", "Time"), sep=" ")
  }
  if ("Prey.Type.Notes" %in% colnames(timestamps)) {
    colnames(timestamps)[which(names(timestamps) == "Prey.Type.Notes")] <- "Prey"
  }
  if ("Prey.Type" %in% colnames(timestamps)) {
    colnames(timestamps)[which(names(timestamps) == "Prey.Type")] <- "Prey"
  }
  if ("Time.End" %in% colnames(timestamps)) {
    colnames(timestamps)[which(names(timestamps) == "Time.End")] <- "Time.Stop"
  }
  
  #Add Behavior and Prey column to main data set, based on data in timestamps 
  data$Time <- chron(times=data$Time)
  timestamps$Time.Start <- chron(times=timestamps$Time.Start)
  timestamps$Time.Stop <- chron(times=timestamps$Time.Stop)
  data$Behavior <- ifelse(timestamps$Time.Start <= data$Time & timestamps$Time.Stop >= data$Time,
                          timestamps$Behavior, NA)
  data$Prey <- ifelse(timestamps$Time.Start <= data$Time & timestamps$Time.Stop >= data$Time,
                      timestamps$Prey, NA)
  #Write to CSV
  write.csv(data, file=filename, row.names = FALSE)
}

line_plot <- function(data, filename){
  plotx <- ggplot(data, aes(x=Time, y=X_dynamic))+
    geom_line(colour="navy")+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ploty <- ggplot(data, aes(x=Time, y=Y_dynamic))+
    geom_line(colour="dark red")+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  plotz <- ggplot(data, aes(x=Time, y=Z_dynamic))+
    geom_line(colour="dark green")+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  plot <- grid.arrange(plotx, ploty, plotz, nrow=3)
  ggsave(filename, plot)
}

hist_plot <- function(data, filename){
  plotx <- ggplot(data, aes(X_dynamic))+
    geom_histogram(colour="navy", fill="light blue")+
    theme_minimal()
  ploty <- ggplot(data, aes(Y_dynamic))+
    geom_histogram(colour="dark red", fill="salmon")+
    theme_minimal()
  plotz <- ggplot(data, aes(Z_dynamic))+
    geom_histogram(colour="dark green", fill="light green")+
    theme_minimal()
  plot <- grid.arrange(plotx, ploty, plotz, nrow=3)
  ggsave(filename, plot)
}

acf_plot <- function(data, filename){
  plotx <-ggacf(data$X_dynamic)+
    theme_minimal()
  ploty <- ggacf(data$Y_dynamic)+
    theme_minimal()
  plotz <- ggacf(data$Z_dynamic)+
    theme_minimal()
  plot <- grid.arrange(plotx, ploty, plotz, nrow=3)
  ggsave(filename, plot)
}

filtered_line_plot <- function(data, filename){
  plotx <- ggplot(data, aes(x=Time, y=X_dynamic, colour = Behavior))+
    geom_point(size=0.5)+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size = 7))
  ploty <- ggplot(data, aes(x=Time, y=Y_dynamic, colour = Behavior))+
    geom_point(size=0.5)+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size = 7))
  plotz <- ggplot(data, aes(x=Time, y=Z_dynamic, colour = Behavior))+
    geom_point(size=0.5)+
    theme_minimal()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size = 7))
  plot <- grid_arrange_shared_legend(plotx, ploty, plotz, nrow=3, ncol=1, position = "right")
  ggsave(filename, plot)
}

data_file <- "DynamicRaw_BigDaddy_3Apr17.csv"
timestamps_file <- "Timestamps_BigDaddy_20170403.csv"
filename <- "Custom_Lady_27Mar17.csv"
get_custom_data(data_file, timestamps_file, filename)
data <- read.csv(filename)
filtered_data <- data %>% filter(!is.na(Behavior))

line_plot(data, "BigDaddy_3Apr17_line_plot.png")
hist_plot(data, "BigDaddy_3Apr17_histogram.png")
acf_plot(data, "BigDaddy_3Apr17_acf.png")
filtered_line_plot(filtered_data, "Lady_27Mar17_line_plot_filtered.png")
