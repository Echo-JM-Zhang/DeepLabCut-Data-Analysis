library(ggplot2)
library(ggpubr)
library(car)
library(grid)
library(sp)
library(ptinpoly)
library(imputeTS)

readDLCdata <- function(file, fps = 24.99, start.time = NULL, end.time = NULL){
  out <- list()
  
  raw.data <- read.csv(file, skip = 2, header = TRUE)
  body.part <- as.character(read.csv(file, nrows = 1, header = TRUE))
  body.part <- data.frame(body.part)
  body.part <- body.part[-1,]
  n.body <- length(body.part) %/% 3
  
  n.frames <- dim(raw.data)[1]
  
  frames <- raw.data$coords
  x <- raw.data[, 2]
  y <- raw.data[, 3]
  likelihood <- raw.data[, 4]
  bodypart <- rep(body.part[3], n.frames)
  
  for (i in 2:n.body){
    x <- c(x, raw.data[, 3*i-1])
    y <- c(y, raw.data[, 3*i])
    likelihood <- c(likelihood, raw.data[, 3*i+1])
    bodypart <- c(bodypart, rep(body.part[i*3], n.frames))
    frames <- c(frames, raw.data$coords)
  }
  
  out$DLCdata <- data.frame(frames, bodypart, x, y, likelihood)
  out$totalframe <- n.frames
  out$fps <- fps
  out$videolength <- round(n.frames / fps)
  
  if (is.null(start.time)){
        start.time <- 0
      }
    
      if (is.null(end.time)){
        end.time <- out$videolength
      }
  
  out$start.time <- start.time
  out$end.time <- end.time
  
  return(out)
}

#cutVideo <- function(mydata, start.time = NULL, end.time = NULL){
#  out <- list()
  
#  if (is.null(start.time)){
#    start.time <- 0
#  }
  
#  if (is.null(end.time)){
#    end.time <- mydata$end.time
#  }
#  
#  start.frame <- start.time * mydata$fps
#  end.frame <- end.time * mydata$fps
#  new.dlc <- mydata$DLCdata[which(mydata$DLCdata$frames >=start.frame & mydata$DLCdata$frames <= end.frame),]
#  
#  out$DLCdata <- new.dlc
#  out$totalframe <-mydata$totalframe
#  out$fps <- mydata$fps
#  out$videolength <- end.time - start.time
#  out$start.time <- start.time
#  out$end.time <- end.time
#  
#  return(out)
#}

drawPointPlot <- function(mydata, trackingPoint = "median", pcutoff = 0){
  point <- mydata$DLCdata[which(mydata$DLCdata$bodypart == trackingPoint),]
  
  point <- point[which(point$likelihood >= pcutoff),]
  
  
  p <- ggplot(point, aes(x, y))+geom_point(aes(color = likelihood))+ggtitle(trackingPoint)
  return(p)
}

# If the point is not in the identified zone, return 0; others return 1
zone.check <- function(mydata, x, y){
  # out <- pip2d(Vertices = cbind(x,y), Queries = cbind(mydata$x, mydata$y))
  out <- point.in.polygon(point.x = mydata$x, point.y = mydata$y, pol.x = x, pol.y = y)
  i <- which(out != 0)
  for (j in 1:length(i)){
    out[i[j]] <- 1
  }
  #i <- which(out == -1)
  #for (j in 1:length(i)){
  #  out[i[j]] <- 0
  #}
  return(out)
}

integratevector <- function(x){
  if(length(x) < 2){
    stop("can  not integrate a vector of length < 2")
  }
  append(0, x[2:length(x)] - x[1:(length(x)-1)])
}

dataCleaning <- function(mydata, zone.x, zone.y, likelihood.cutoff = 0.8, maxdelta = NULL){
  process <- which(mydata$DLCdata$likelihood < likelihood.cutoff)
  
  zone <- list()
  inzone <- rep(0, dim(mydata$DLCdata)[1])
  for (i in 1:length(zone.x)){
    zone <- zone.check(mydata$DLCdata, zone.x[[i]], zone.y[[i]])
    inzone <- inzone + zone.check(mydata$DLCdata, zone.x[[i]], zone.y[[i]])
  }
  not.in.zone <- which(inzone == 0)
  
  distance <- sqrt(integratevector(mydata$DLCdata$x)^2+integratevector(mydata$DLCdata$y)^2)
  jump <- which(distance >=maxdelta)
  
  process <- c(process, not.in.zone, jump)
  mydata$DLCdata$x[process] <- NA
  mydata$DLCdata$y[process] <- NA
  
  mydata$DLCdata <- na_interpolation(mydata$DLCdata)
  
  median.x <- c()
  median.y <- c()
  median.likelihood <- c()
  for (i in 0:max(mydata$DLCdata$frames)){
    position <- which(mydata$DLCdata$frames == i)
    median.x <- c(median.x, median(mydata$DLCdata$x[position]))
    median.y <- c(median.y, median(mydata$DLCdata$y[position]))
    median.likelihood <- c(median.likelihood, median(mydata$DLCdata$likelihood[position]))
  }
  i <-0:max(mydata$DLCdata$frames)
  median.track <- data.frame(frames = i, bodypart = rep("median", max(mydata$DLCdata$frames)+1), x = median.x, y = median.y, likelihood = median.likelihood)
  mydata$DLCdata <- rbind(mydata$DLCdata, median.track)
  
  return(mydata)
}

#approach-avoidance analysis

#set timepoint for actual time in arm
#time.interval is the time (in seconds) required to identify that rat has entered the maze (door closed)
entry.exit.time <- function(mydata, trackingpart = "median", hub.x, hub.y, time.interval = 150){
  trackingdata <- mydata$DLCdata[which(mydata$DLCdata$bodypart == trackingpart),]
  n <- dim(trackingdata)[1]
  trackingdata$inhub <- zone.check(trackingdata, hub.x, hub.y)
  t <- zone.check(mydata$DLCdata, hub.x, hub.y)
  for (i in 1:length(t)){
    if (t[i] == 1){
      t[i] <- "hub"
    }else{
      t[i] <- "arm"
    }
  }
  mydata$DLCdata$zone <- t
  
  # entering time
  find <- FALSE
  t <- mydata$start.time
  i <- 0
  while (! find){
    i <- i+1
    if (trackingdata$inhub[i] == 0){
      find <- TRUE
      enter <- i
      for (j in (i+1):(i+1+time.interval*mydata$fps)){
        if (trackingdata$inhub[j] == 1){
          find <- FALSE
        }
      }
    }
  }
  
  # exiting time
  find <- FALSE
  t <- mydata$end.time
  i <- n
  while(! find){
    i <- i-1
    if (trackingdata$inhub[i] == 0){
      find <- TRUE
      exit <- i
      for (j in (i-1-time.interval*mydata$fps):(i-1)){
        if (trackingdata$inhub[j] == 1){
          find <- FALSE
        }
      }
    }
  }
  mydata$enter <- enter
  mydata$exit <- exit
  return(mydata)
}
  

#latency to enter

AAconflict.latency <- function(mydata, trackingpart = "median", type = "enter"/"exit"){
  trackingdata <- mydata$DLCdata[which(mydata$DLCdata$bodypart == trackingpart),]
  if (type == "enter"){
    t <- mydata$start.time
    k <- round(t*mydata$fps)+1
    temp <- trackingdata$zone[k:mydata$enter]
    in.hub <- length(which(temp == "hub"))
    time <- in.hub/mydata$fps
  }else if (type == "exit"){
    t <- mydata$end.time
    k <- round(t*mydata$fps)-1
    temp <- trackingdata$zone[k:mydata$exit]
    in.arm <- length(which(temp != "hub"))
    time <- in.arm / mydata$fps
    
  }
  
  return(round(time, 2))
}

#time spent in approach/avoidance

time.approach.avoidance <- function(mydata, trackingpart = "median", arm.x, arm.y){
  x.av <- c(arm.x[1], mean(arm.x[1:2]), mean(arm.x[3:4]), arm.x[4])
  y.av <- c(arm.y[1], mean(arm.y[1:2]), mean(arm.y[3:4]), arm.y[4])
  
  avoidance <- zone.check(mydata$DLCdata, x.av, y.av)
  for (i in 1:length(avoidance)){
    if (mydata$DLCdata$zone[i] != "hub"){
      if (avoidance[i] == 1){
        mydata$DLCdata$zone[i] <- "av"
      }else{
        mydata$DLCdata$zone[i] <- "ap"
      }
    }
  }
  
  
  maze.start <- mydata$enter
  maze.end <- trunc(mydata$end.time*mydata$fps)+1
  process <- mydata$DLCdata[which(mydata$DLCdata$bodypart == trackingpart),]
  
  process <- process[maze.start:maze.end,]
  av.frames <- length(which(process$zone == "av"))
  ap.frames <- length(which(process$zone == "ap"))
  mydata$av.time <- round(av.frames/mydata$fps, 2)
  mydata$ap.time <- round(ap.frames/mydata$fps, 2)
  
  return(mydata)
}

plotZoneVisit <- function(mydata, bodypart = "median"){
  fps <- mydata$fps
  process <- mydata$DLCdata[which(mydata$DLCdata$bodypart == bodypart),]
  process$seconds <- process$frames/fps
  plot <- ggplot(data = process, aes(seconds, zone, color = zone)) + 
    geom_point(size = 4, shape = 124) + 
    ggtitle(bodypart) + theme_bw()
  
  return(plot)
}

generateHeatMap <- function(mydata){
  xbreaks <- seq(floor(min(mydata$DLCdata$x)), ceiling(max(mydata$DLCdata$x)), by = 0.1)
  ybreaks <- seq(floor(min(mydata$DLCdata$y)), ceiling(max(mydata$DLCdata$y)), by = 0.1)
  mydata$DLCdata$latbin <- xbreaks[cut(mydata$DLCdata$x, breaks = xbreaks, labels=F)]
  mydata$DLCdata$longbin <- ybreaks[cut(mydata$DLCdata$y, breaks = ybreaks, labels=F)]
  
  p <- ggplot(data = mydata$DLCdata, aes(latbin,longbin)) + 
    stat_density_2d(aes(fill=..density..), geom = "raster", contour = FALSE)+
    scale_fill_gradient(name = "Time Density", low = "blue", high = "yellow")
  return(p)
}

AddZones <- function(plot, x, y){
  for(i in 1:length(x)){
    process <- data.frame(x=c(x[[i]], x[[i]][1]), y=c(y[[i]], y[[i]][1]))
    p <- p + geom_path(data=process,aes(x,y))
  }
  return(p)
}




setwd("D:/JingminZhang/rutsuko's lab/behavior scoring/batch1M")

x.arm <- c(306.2622, 764.9192, 751.3294, 291.5399)
y.arm <- c(156.9536, 219.2233, 321.1461, 265.6543)
x.hub <- c(163.5689625, 60.51270239, 76.36751163, 198.6760401, 306.2622, 291.5399)
y.hub <- c(339.2658723, 259.9918261, 128.6234067, 75.39654704, 156.9536, 265.6543)
x.maze <- list(x.arm, x.hub)
y.maze <- list(y.arm, y.hub)
x.middle <- c(mean(c(306.2622, 764.9192)), mean(c(751.3294, 291.5399)))
y.middle <- c(mean(c(156.9536, 219.2233)), mean(c(321.1461, 265.6543)))
num.zone <- 2


filename <- c("RAT31.csv", "RAT32.csv", "RAT33.csv", "RAT34.csv", "RAT35.csv", "RAT36.csv", "RAT37.csv", "RAT39.csv")
analyzedDLC <- read.csv("DLCanalysis.csv")
k <- 1
latency.enter <- c()
latency.exit <- c()
av.time <- c()
ap.time <- c()

for (file in filename){
  starttime <- analyzedDLC$Start.time[k]
  endtime <- analyzedDLC$End.time[k]
  k <- k+1
  data <- readDLCdata(file, fps = 24.99, start.time = starttime, end.time = endtime)
  #drawPointPlot(data, "tail", 0.5)
  
  cleaned.data <- dataCleaning(data, x.maze, y.maze, likelihood.cutoff = 0.3)
  
  #drawPointPlot(cleaned.data, "median", 0.3)
  median <- entry.exit.time(cleaned.data, trackingpart = "median", hub.x = x.hub, hub.y = y.hub)
  latency_to_enter <- AAconflict.latency(median, trackingpart = "median", type = "enter")
  latency_to_exit <- AAconflict.latency(median, trackingpart = "median", type = "exit")
  analyzed <- time.approach.avoidance(median, trackingpart = "median", arm.x = x.arm, arm.y = y.arm)
  # write.csv(analyzed, "D:/JingminZhang/rutsuko's lab/behavior scoring/batch1M/test1.csv")
  latency.enter <- c(latency.enter, latency_to_enter)
  latency.exit <- c(latency.exit, latency_to_exit)
  av.time <- c(av.time, analyzed$av.time)
  ap.time <- c(ap.time, analyzed$ap.time)
  
  p1 <- plotZoneVisit(analyzed)
  ggsave(filename = paste(substr(file, 1, 5), "zoneVisit.jpeg", sep = "_"), plot = p1, device = "jpeg", width = 4, height = 2.5)
  
  p <- generateHeatMap(analyzed)
  
  p <- AddZones(p, x.maze, y.maze)+geom_line(data = data.frame(x=x.middle, y=y.middle), aes(x, y))
  ggsave(filename = paste(substr(file, 1, 5), "heatmap.jpeg", sep = "_"), plot = p, device = "jpeg", width = 6, height = 2.5)
}
analyzedDLC$latency.enter <- latency.enter
analyzedDLC$latency.exit <- latency.exit
analyzedDLC$ap.time <- ap.time
analyzedDLC$av.time <- av.time

write.csv(analyzedDLC, file = "analyzedDLC.csv")

    
