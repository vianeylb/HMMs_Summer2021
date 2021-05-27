library(ggplot2)
library(forecast)
library(gridExtra)

#timeseries of observations (works for all types of distributions)
timeseries <- function(hmmdata,
                          title = 'Time Series of Simulated HMM',
                          xlabel = 'Time',
                          ylabel = 'Observation') {
  ggplot(hmmdata, aes(x = Time, y = Observation)) +
    theme_light() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(aes(color = State))+
    geom_line(colour='grey', alpha=0.8, lwd=0.4) +
    labs(x = xlabel, y = ylabel)
}


timeseriesfit <- function(hmmdata,
                       title = 'Time Series of Simulated HMM',
                       xlabel = 'Time',
                       ylabel = 'Observation') {
  ggplot(hmmdata, aes(x = Time, y = Observation)) +
    theme_light() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(aes(color = GuessState))+
    geom_line(colour='grey', alpha=0.8, lwd=0.4) +
    labs(x = xlabel, y = ylabel)
}

#example
#timeseries(pois3sdatasm)


#function to plot time series of observations with true states and fit states (poisson and normal)
#hmmdata has columns: "Time", "Observation", "State" (true state if known), "GuessState"
#use fit and true lambda for poisson and fit and true mu for normal in place of fitparam truparam
timeseries_states <- function(m, hmmdata, fitparam, truparam = NULL){
  p <- ggplot(hmmdata, aes(x = Time, y = Observation)) +
    ggtitle('Timeseries of Observations and Underlying States') +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = 'Time', y = 'Observation') +
    geom_line()
  if (!is.null(truparam)){
    p <- p + geom_hline(yintercept = truparam, color = 2:(m+1), alpha = 0.7, lwd = 0.4) + 
      geom_point(aes(x = Time, y = truparam[State], colour = State)) +
      labs(colour = "State") +
      geom_hline(yintercept = fitparam, linetype = 2, alpha = 0.7, color = 2:(m+1), lwd = 0.4) + 
      geom_point(aes(Time, y = fitparam[GuessState], colour = State)) +
      labs(colour = "State")}
  else {
    p <- p + geom_hline(yintercept = fitparam, linetype = 2, alpha = 0.7, color = 2:(m+1), lwd = 0.4) + 
      geom_point(aes(Time, y = fitparam[GuessState], colour = GuessState)) +
      labs(colour = "State")}
  return(p)
}

#example
#timeseries_states(3, pois3sdatasm, modpois3s$lambda, pois3s$lambda)
#timeseries_states(3, normdatahead, modnorm3s$mu, norm3s$mu)
  

pois_hist_dist <- function(m, hmmdata, mod,
                              title = 'Histogram of HMM with State Distributions',
                              xlabel = 'Observation', 
                              ylabel = 'Frequency', 
                              numbreaks=25){
  observ = hmmdata$Observation
  state = hmmdata$State
  h <- hist(observ, main=title, xlab=xlabel, ylab=ylabel, breaks=numbreaks)
  xfit <- seq(min(observ), max(observ))
  colours = c('blue', 'red', 'green', 'darkgreen', 'gold', 'deeppink', 'purple', 'brown')
  name = c('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8')
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dpois(xfit, mod$lambda[i])
    yfit <- yfit * sum(state==i) * diff(h$mids[1:2])
    lines(xfit, yfit, lwd = 2, col=colours[i])
    marginal <- marginal + yfit
  }
  legend("topright", name[1:m], col=colours[1:m], lwd=2)
  lines(xfit, marginal, lwd = 2, col="black")
}




#plot histogram of observations with overlayed fit distributions poisson
pois_hist_fitdist <- function(m, hmmdata, mod,
                           title = 'Histogram of HMM with Fitted State Distributions',
                           xlabel = 'Observation', 
                           ylabel = 'Frequency', 
                           numbreaks=25){
  observ = hmmdata$Observation
  fitstate = hmmdata$GuessState
  h <- hist(observ, main=title, xlab=xlabel, ylab=ylabel, breaks=numbreaks)
  xfit <- seq(min(observ), max(observ))
  colours = c('blue', 'red', 'green', 'darkgreen', 'gold', 'deeppink', 'purple', 'brown')
  name = c('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8')
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dpois(xfit, mod$lambda[i])
    yfit <- yfit * sum(fitstate==i) * diff(h$mids[1:2])
    lines(xfit, yfit, lwd = 2, col=colours[i])
    marginal <- marginal + yfit
  }
  legend("topright", name[1:m], col=colours[1:m], lwd=2)
  lines(xfit, marginal, lwd = 2, col="black")
}
  

#example
#pois_hist_dist(3, pois3sdata, modpois3s)
#pois_hist_dist(2, pois3sdata, modpois2s)


#plot histogram of observations with overlayed fit distributions normal
norm_hist_dist <- function(m, hmmdata, mod,
                           title = 'Histogram of HMM with Fitted State Distributions',
                           xlabel = 'Observation', 
                           ylabel = 'Frequency', 
                           numbreaks=25){
  observ = hmmdata$Observation
  state = hmmdata$State
  h <- hist(observ, main=title, xlab=xlabel, ylab=ylabel, breaks=numbreaks)
  xfit <- seq(min(observ), max(observ), length.out = 100)
  colours = c('blue', 'red', 'green', 'darkgreen', 'gold', 'deeppink', 'purple', 'brown')
  name = c('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8')
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dnorm(xfit, mod$mu[i], mod$sigma[i])
    yfit <- yfit * sum(state==i) * diff(h$mids[1:2])
    lines(xfit, yfit, lwd = 2, col=colours[i])
    marginal <- marginal + yfit
  }
  legend("topright", name[1:m], col=colours[1:m], lwd=2)
  lines(xfit, marginal, lwd = 2, col="black")
}


#plot histogram of observations with overlayed fit distributions normal
norm_hist_fitdist <- function(m, hmmdata, mod,
                           title = 'Histogram of HMM with Fitted State Distributions',
                           xlabel = 'Observation', 
                           ylabel = 'Frequency', 
                           numbreaks=25){
  observ = hmmdata$Observation
  fitstate = hmmdata$GuessState
  h <- hist(observ, main=title, xlab=xlabel, ylab=ylabel, breaks=numbreaks)
  xfit <- seq(min(observ), max(observ), length.out = 100)
  colours = c('blue', 'red', 'green', 'darkgreen', 'gold', 'deeppink', 'purple', 'brown')
  name = c('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8')
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dnorm(xfit, mod$mu[i], mod$sigma[i])
    yfit <- yfit * sum(fitstate==i) * diff(h$mids[1:2])
    lines(xfit, yfit, lwd = 2, col=colours[i])
    marginal <- marginal + yfit
  }
  legend("topright", name[1:m], col=colours[1:m], lwd=2)
  lines(xfit, marginal, lwd = 2, col="black")
}


#example
#norm_hist_dist(3, normdata, modnorm3s)
#norm_hist_dist(2, normdata, modnorm2s)


#plot histogram of observations with overlayed fit distributions gamma
gam_hist_dist <- function(m, hmmdata, mod,
                           title = 'Histogram of HMM with Fitted State Distributions',
                           xlabel = 'Observation', 
                           ylabel = 'Frequency', 
                           numbreaks=40){
  observ = hmmdata$Observation
  state = hmmdata$State
  h <- hist(observ, main=title, xlab=xlabel, ylab=ylabel, breaks=numbreaks)
  xfit <- seq(min(observ), max(observ), length.out = 100)
  colours = c('blue', 'red', 'green', 'darkgreen', 'gold', 'deeppink', 'purple', 'brown')
  name = c('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8')
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dgamma(xfit, shape=mod$alpha[i], scale=mod$theta[i])
    yfit <- yfit * sum(state==i) * diff(h$mids[1:2])
    lines(xfit, yfit, lwd = 2, col=colours[i])
    marginal <- marginal + yfit
  }
  legend("topright", name[1:m], col=colours[1:m], lwd=2)
  lines(xfit, marginal, lwd = 2, col="black")
}


#plot histogram of observations with overlayed fit distributions gamma
gam_hist_fitdist <- function(m, hmmdata, mod,
                          title = 'Histogram of HMM with Fitted State Distributions',
                          xlabel = 'Observation', 
                          ylabel = 'Frequency', 
                          numbreaks=40){
  observ = hmmdata$Observation
  fitstate = hmmdata$GuessState
  h <- hist(observ, main=title, xlab=xlabel, ylab=ylabel, breaks=numbreaks)
  xfit <- seq(min(observ), max(observ), length.out = 100)
  colours = c('blue', 'red', 'green', 'darkgreen', 'gold', 'deeppink', 'purple', 'brown')
  name = c('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8')
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dgamma(xfit, shape=mod$alpha[i], scale=mod$theta[i])
    yfit <- yfit * sum(fitstate==i) * diff(h$mids[1:2])
    lines(xfit, yfit, lwd = 2, col=colours[i])
    marginal <- marginal + yfit
  }
  legend("topright", name[1:m], col=colours[1:m], lwd=2)
  lines(xfit, marginal, lwd = 2, col="black")
}


#plot normal pseudo residuals
pr.plot.discr <- function(data, labs=TRUE){
  p <- ggplot(data)+
    #geom_point(aes(x=Time, y=PR_lo), size=0.5, colour="blue")+
    #geom_point(aes(x=Time, y=PR_hi), size=0.5, colour="red")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_segment(aes(x = Time, xend = Time, y = PR_lo, yend = PR_hi), alpha=0.6)+
    geom_hline(yintercept=0, linetype=2, alpha = 0.9)+
    ylab("Normal Pseudo Residuals")
  if (labs==FALSE){
    p <- p + theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank())
  }
  return(p)
}

pr.plot.cont <- function(data, labs=TRUE){
  p <- ggplot(data)+
    geom_point(aes(x=Time, y=PR), size=0.04, alpha=0.6)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_hline(yintercept=0, linetype=2, alpha = 0.9)+
    ylab("Normal Pseudo Residuals")
  if (labs==FALSE){
    p <- p + theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank())
  }
  return(p)
}


#plot normal histogram of residuals
pr.hist <- function(data, labs=TRUE){
  p <- ggplot(data, aes(PR))+
    geom_histogram(aes(y=..density..), binwidth=0.5, colour="white", fill="grey")+
    stat_function(fun=dnorm, colour="black")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab('mid-pseudo-residuals')+
    xlim(c(-4,4))
  if (labs==FALSE){
    p <- p + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank())
  }
  return(p)
}


#quantile-quantile plot
pr.qq <- function(data, labs=TRUE){
  p <- ggplot(data, aes(sample=PR, alpha=0.5))+
    stat_qq(size=0.1)+
    stat_qq_line()+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none")
  if (labs==FALSE){
    p <- p + theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   legend.position = "none",
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank())
  }
  return(p)
}


#acf plot mid-pseudo-residuals
pr.acf <- function(pseudo, labs=TRUE){
  p <- ggAcf(pseudo)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_blank())
  if (labs==FALSE){
    p <- p + theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   plot.title = element_blank(),
                   axis.title.y = element_blank())
  }
  return(p)
}



  
