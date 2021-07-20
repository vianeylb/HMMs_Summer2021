line_plot_log_odba <- function(data, filename) {
  plot <- ggplot(data, aes(x = Time, y = log(ODBA))) +
    geom_line(colour = "dark red") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ggsave(paste(filename, "log_ODBA.png", sep = "_"), plot)
}

hist_plot_log_odba <- function(data, filename) {
  plot <- ggplot(data, aes(log(ODBA))) +
    geom_histogram(colour = "dark red", fill = "salmon") +
    theme_minimal()
  ggsave(paste(filename, "log_ODBA.png", sep = "_"), plot)
}

acf_plot_log_odba<- function(data, filename) {
  plot <- ggacf(log(data$ODBA)) +
    theme_minimal()
  ggsave(paste(filename, "log_ODBA.png", sep = "_"), plot)
}

pacf_plot_log_odba<- function(data, filename) {
  plot <- ggpacf(log(data$ODBA)) +
    theme_minimal()
  ggsave(paste(filename, "log_ODBA.png", sep = "_"), plot)
}

behavior_plot_log_odba <- function(data, filename) {
  plot <- ggplot(data, aes(x = Time, y = log(ODBA), colour = Behavior)) +
    geom_line() +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = 7)
    ) +
    scale_color_brewer(palette = "Set1")
  ggsave(paste(filename, "log_ODBA.png", sep = "_"), plot)
}

filtered_hist_log_odba<- function(data, filename) {
  plot <- ggplot(data, aes(x = log(ODBA), colour = Behavior, fill = Behavior)) +
    geom_histogram() +
    facet_wrap(~Behavior, ncol = 2) +
    theme_minimal() +
    scale_fill_brewer(palette = "Pastel1") +
    scale_color_brewer(palette = "Set1")
  ggsave(paste(filename, "log_ODBA.png", sep = "_"), plot)
  
  behaviors <- unique(data$Behavior)
  for (i in seq_len(length(behaviors))) {
    behavior <- behaviors[i]
    subdata <- data %>% dplyr::filter(Behavior == behavior)
    
    plot <- ggplot(subdata, aes(x = log(ODBA))) +
      geom_histogram(colour = "dark red", fill = "salmon") +
      theme_minimal()
    ggsave(paste(filename, behavior, "log_ODBA.png", sep = "_"), plot)
  }
}

behavior_hist_log_odba <- function(data, filename) {
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
    
    plot <- ggplot(data = subdata,
                        aes(x = log(ODBA),
                            colour = BehaviorIndex,
                            fill = BehaviorIndex)) +
      geom_histogram() +
      facet_wrap(~BehaviorIndex, ncol = 2) +
      theme_minimal() +
      scale_fill_brewer(palette = "Pastel1") +
      scale_color_brewer(palette = "Set1")
    ggsave(paste(filename, behavior, "log_ODBA.png", sep = "_"), plot)
  }
}

get_plots_log_odba <- function(names){
  n <- length(names)
  for (i in 1:n){
    name <- names[i]
    filename <- paste("Custom", name, "dynamic.csv", sep = "_")
    data <- read.csv(filename)
    labelled_data <- data %>% filter(!is.na(Behavior))
    
    line_plot_log_odba(data, paste(name, "line_plot", sep = "_"))
    hist_plot_log_odba(data, paste(name, "histogram", sep = "_"))
    acf_plot_log_odba(data, paste(name, "acf", sep = "_"))
    pacf_plot_log_odba(data, paste(name, "pacf", sep = "_"))
    behavior_plot_log_odba(data, paste(name, "plot_behavior", sep = "_"))
    behavior_plot_log_odba(labelled_data, paste(name, "plot_behavior_filtered", sep = "_"))
    filtered_hist_log_odba(labelled_data, paste(name, "histogram_filtered", sep = "_"))
    behavior_hist_log_odba(labelled_data, paste(name, "histogram_behavior", sep = "_"))
  }
}

names <- c("BigDaddy_3Apr17",
           "BigDaddy_20Mar17",
           "BigGuy_15Feb18",
           "Eliza_7Sept17",
           "Eliza_20Sept17",
           "Lady_27Mar17")
get_plots_log_odba(names)
