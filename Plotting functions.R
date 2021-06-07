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
timeseries_states <- function(m, hmmdata, fitparam, truparam = NULL, CIup = NULL, CIlow=NULL){
  p <- ggplot() +
    theme_light() +
    ggtitle('Timeseries of Observations and Underlying States') +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = 'Time', y = 'Observation') +
    geom_line(data=hmmdata, aes(x = Time, y = Observation)) 
  
  if (!is.null(CIup)){
    for (i in 1:m){
      df <- data.frame('x'= hmmdata$Time, 'ymin' = rep(CIlow[i], length(hmmdata$Time)),
                       'ymax' = rep(CIup[i], length(hmmdata$Time)))
      p <- p + geom_ribbon(data=df, aes(x = x, ymin= ymin, ymax=ymax), fill=(i+1), alpha=0.3)
    }
  }
  if (!is.null(truparam)){
    p <- p + geom_hline(yintercept = truparam, color = 2:(m+1), alpha = 0.7, lwd = 0.4) +
      geom_point(data=hmmdata, aes(x = Time, y = truparam[State], colour = State)) +
      labs(colour = "State") +
      geom_hline(yintercept = fitparam, linetype = 2, alpha = 0.7, color = 2:(m+1), lwd = 0.4) +
      geom_point(data=hmmdata, aes(Time, y = fitparam[GuessState], colour = State)) +
      labs(colour = "State")}
  else {
    p <- p + geom_hline(yintercept = fitparam, linetype = 2, alpha = 0.7, color = 2:(m+1), lwd = 0.4) +
      geom_point(data=hmmdata, aes(Time, y = fitparam[GuessState], colour = GuessState)) +
      labs(colour = "State")}
  return(p)
}

#example
#timeseries_states(3, pois3sdatasm, modpois3s$lambda, pois3s$lambda)
#timeseries_states(3, normdatahead, modnorm3s$mu, norm3s$mu)


########################## POISSON ########################

pois_hist_dist <- function(m, hmmdata, mod, width=1){
  observ = hmmdata$Observation
  state = hmmdata$State
  
  h <- ggplot() + 
    geom_histogram(data=hmmdata, 
                   aes(x=Observation), 
                   binwidth = width,
                   colour="cornsilk4",
                   fill="white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  xfit <- seq(min(observ), max(observ))
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dpois(xfit, mod$lambda[i])
    yfit <- yfit * sum(state==i) * width
    df <- data.frame('xfit' = xfit, 'yfit' = yfit, col = as.factor(rep(i, length(xfit))))
    h <- h + geom_line(data=df,aes(xfit, yfit, colour=col), lwd=0.7)
    marginal <- marginal + yfit
  }
  h <- h + labs(color = "State")
  df <- data.frame('xfit' = xfit, 'yfit' = marginal)
  h <- h + geom_line(data=df,aes(xfit, yfit), col="black", lwd=0.7)
  
  return(h)
}


pois_hist_dist_CI <- function(m, hmmdata, mod, CI, width=1){
  observ = hmmdata$Observation
  fitstate = hmmdata$GuessState
  
  h <- ggplot() + 
    geom_histogram(data=hmmdata, 
                   aes(x=Observation), 
                   binwidth = width,
                   colour="cornsilk4",
                   fill="white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  xfit <- seq(min(observ), max(observ))
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dpois(xfit, mod$lambda[i])
    yfit <- yfit * sum(fitstate==i) * width
    df <- data.frame('xfit' = xfit, 'yfit' = yfit, col = as.factor(rep(i, length(xfit))))
    h <- h + geom_line(data=df,aes(xfit, yfit, colour=col), lwd=0.7)
    marginal <- marginal + yfit
  }
  h <- h + labs(color = "State")
  df <- data.frame('xfit' = xfit, 'yfit' = marginal)
  h <- h + geom_line(data=df,aes(xfit, yfit), col="black", lwd=0.7)
  
  for (i in 1:m){
    up <- dpois(xfit, CI$lambda.upper.conf[i])
    up <- up * sum(fitstate==i) * width
    down <- dpois(xfit, CI$lambda.lower.conf[i])
    down <- down * sum(fitstate==i) * width
    df <- data.frame('x' = xfit, 'upper' = up, 'lower' = down)
    h <- h + geom_ribbon(data=df, aes(x = xfit, ymin= lower, ymax=upper), fill=(i+1), alpha=0.4)
  }
  return(h)
}


############################ NORMAL ######################

#plot histogram of observations with overlayed fit distributions normal
norm_hist_dist <- function(m, hmmdata, mod, width=1){
  observ = hmmdata$Observation
  state = hmmdata$State
  
  h <- ggplot() + 
    geom_histogram(data=hmmdata, 
                   aes(x=Observation), 
                   binwidth = width,
                   colour="cornsilk4",
                   fill="white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  xfit <- seq(min(observ), max(observ))
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dnorm(xfit, mod$mu[i], mod$sigma[i])
    yfit <- yfit * sum(state==i) * width
    df <- data.frame('xfit' = xfit, 'yfit' = yfit, col = as.factor(rep(i, length(xfit))))
    h <- h + geom_line(data=df,aes(xfit, yfit, colour=col), lwd=0.7)
    marginal <- marginal + yfit
  }
  h <- h + labs(color = "State")
  df <- data.frame('xfit' = xfit, 'yfit' = marginal)
  h <- h + geom_line(data=df,aes(xfit, yfit), col="black", lwd=0.7)
  
  return(h)
}


norm_hist_dist_CI <- function(m, hmmdata, mod, CI, width=1){
  observ = hmmdata$Observation
  fitstate = hmmdata$GuessState
  
  h <- ggplot() + 
    geom_histogram(data=hmmdata, 
                   aes(x=Observation), 
                   binwidth = width,
                   colour="cornsilk4",
                   fill="white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  xfit <- seq(min(observ), max(observ))
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dnorm(xfit, mod$mu[i], mod$sigma[i])
    yfit <- yfit * sum(fitstate==i) * width
    df <- data.frame('xfit' = xfit, 'yfit' = yfit, col = as.factor(rep(i, length(xfit))))
    h <- h + geom_line(data=df,aes(xfit, yfit, colour=col), lwd=0.7)
    marginal <- marginal + yfit
  }
  h <- h + labs(color = "State")
  df <- data.frame('xfit' = xfit, 'yfit' = marginal)
  h <- h + geom_line(data=df,aes(xfit, yfit), col="black", lwd=0.7)
  
  for (k in 1:m){
    upper <- CI$upper[k,] * sum(fitstate==k)
    lower <- CI$lower[k,] * sum(fitstate==k)
    df <- data.frame('x' = CI$range, 'upper' = upper, 'lower' = lower)
    h <- h + geom_ribbon(data=df, aes(x = x, ymin=lower, ymax=upper), fill=(k+1), alpha=0.4)
  }
  return(h)
}

###################### GAMMA ####################

#plot histogram of observations with overlayed fit distributions gamma
gam_hist_dist <- function(m, hmmdata, mod, width=1){
  observ = hmmdata$Observation
  state = hmmdata$State
  
  h <- ggplot() + 
    geom_histogram(data=hmmdata, 
                   aes(x=Observation), 
                   binwidth = width,
                   colour="cornsilk4",
                   fill="white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  xfit <- seq(min(observ), max(observ))
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dgamma(xfit, shape=mod$alpha[i], scale=mod$theta[i])
    yfit <- yfit * sum(state==i) * width
    df <- data.frame('xfit' = xfit, 'yfit' = yfit, col = as.factor(rep(i, length(xfit))))
    h <- h + geom_line(data=df,aes(xfit, yfit, colour=col), lwd=0.7)
    marginal <- marginal + yfit
  }
  h <- h + labs(color = "State")
  df <- data.frame('xfit' = xfit, 'yfit' = marginal)
  h <- h + geom_line(data=df,aes(xfit, yfit), col="black", lwd=0.7)
  
  return(h)
}


gam_hist_dist_CI <- function(m, hmmdata, mod, CI, width=1){
  observ = hmmdata$Observation
  fitstate = hmmdata$GuessState
  
  h <- ggplot() + 
    geom_histogram(data=hmmdata, 
                   aes(x=Observation), 
                   binwidth = width,
                   colour="cornsilk4",
                   fill="white") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  xfit <- seq(min(observ), max(observ))
  marginal <- numeric(length(xfit))
  
  for (i in 1:m) {
    yfit <- dgamma(xfit, shape=mod$alpha[i], scale=mod$theta[i])
    yfit <- yfit * sum(fitstate==i) * width
    df <- data.frame('xfit' = xfit, 'yfit' = yfit, col = as.factor(rep(i, length(xfit))))
    h <- h + geom_line(data=df,aes(xfit, yfit, colour=col), lwd=0.7)
    marginal <- marginal + yfit
  }
  h <- h + labs(color = "State")
  df <- data.frame('xfit' = xfit, 'yfit' = marginal)
  h <- h + geom_line(data=df,aes(xfit, yfit), col="black", lwd=0.7)
  
  for (k in 1:m){
    upper <- CI$upper[k,] * sum(fitstate==k)
    lower <- CI$lower[k,] * sum(fitstate==k)
    df <- data.frame('x' = CI$range, 'upper' = upper, 'lower' = lower)
    h <- h + geom_ribbon(data=df, aes(x = x, ymin=lower, ymax=upper), fill=(k+1), alpha=0.4)
  }
  return(h)
}

######################### PSEUDO RESIDUALS #####################

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
