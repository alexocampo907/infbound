# library(ggplot2)

# The Basics --------------------------------------------------------------

expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(p) log(p/(1-p))

myhist <- function(vector,binwidth=NULL){
  plot.data <- data.frame(values=vector)
  if(is.null(binwidth)) binwidth <- (max(vector,na.rm = T) - min(vector,na.rm = T))/30
  ggplot(plot.data,aes(values)) + geom_histogram(binwidth=binwidth,aes(fill = ..count..) ) + theme_classic()
}

mytable <- function(factor,digits=2,useNA="always",order=T){
  
  if(!order){x <- table(factor,useNA = useNA)}
  else{
    x <- table(factor,useNA = useNA)[order(table(factor),decreasing = T)]  
  }
  
  return(data.frame(Variable = names(x),Frequency = as.numeric(x),Percent = paste0(100*round(as.numeric(x)/sum(as.numeric(x),na.rm = T),digits),"%")))
}

# MDP Pertinent -----------------------------------------------------------

binaryPermutations  <- function(n, vals = 0:1) do.call(expand.grid,rep(list(vals),n))

# bound.inf -------------------------------------------------------

# formula = y ~ x1 + x2 + x3;data=data;family = "binomial";parameter = 2


bound.inf <- function(formula,parameter = 2,data,family = "binomial"){
  
  bet <- glm(formula,family = family,data = data)$coef
  
  ind.cc <- which(complete.cases(data))
  ind.mis <- which(!complete.cases(data))
  
  # combos!
  udata <- as.matrix(distinct(data[ind.mis,]))
  datmat <- as.matrix(data)
  
  Z <- numeric(nrow(data))
  for(i in 1:nrow(udata)) Z[apply(datmat,1,function(x) identical(udata[i,],x))] <- i+1
  Z[complete.cases(data)] <- 1
  
  max.data <- data
  min.data <- data
  
  for(z in 2:length(unique(Z))){
    
    tmp <- binaryPermutations(sum(is.na(udata[z-1,])))
    deltas <- numeric(nrow(tmp))
    
    for(i in 1:nrow(tmp)){
      new.obs <- as.numeric(udata[z-1,])
      new.obs[is.na(new.obs)] <- as.numeric(tmp[i,])
      deltas[i] <- inf.logistic(new.obs,dat = data[ind.cc,],bet = bet)[parameter]
    }
    
    min.data[Z==z,is.na(min.data[Z==z,][1,])] <- tmp[which.min(deltas),]
    max.data[Z==z,is.na(max.data[Z==z,][1,])] <- tmp[which.max(deltas),]
  }
  
  # Store Results
  result <- c(glm(y ~ x1 + x2 + x3,family = "binomial",data=min.data)$coef[parameter],
              glm(y ~ x1 + x2 + x3,family = "binomial",data=max.data)$coef[parameter])
  
  names(result) <- c("Lower Bound","Upper Bound")
  
  return(list(bounds = result,min.data = min.data,max.data = max.data,bet = bet))
}


# EIF ---------------------------------------------------------------------

# formula = y ~ x1 + x2 + x3;data=data;family = "binomial";parameter = 2

bound.eif <- function(formula,parameter = 2,data,family = "binomial"){

  bet <- glm(formula,family = family,data = data)$coef
  
  ind.cc <- which(complete.cases(data))
  ind.mis <- which(!complete.cases(data))
  
  max.data <- data
  min.data <- data
  
  # inf code
  # combos!
  udata <- as.matrix(distinct(data[ind.mis,]))
  datmat <- as.matrix(data)
  
  Z <- numeric(nrow(data))
  for(i in 1:nrow(udata)) Z[apply(datmat,1,function(x) identical(udata[i,],x))] <- i+1
  Z[complete.cases(data)] <- 1
  
  max.data <- data
  min.data <- data
  
  for(z in 2:length(unique(Z))){
    
    tmp <- binaryPermutations(sum(is.na(udata[z-1,])))
    deltas <- numeric(nrow(tmp))
    
    for(i in 1:nrow(tmp)){
      imp.data <- data
      imp.data[Z==z,is.na(min.data[Z==z,][1,])] <- tmp[i,]
      deltas[i] <- glm(y ~ x1 + x2 + x3,family = "binomial",data=imp.data)$coef[parameter]  - bet[parameter]
    }
    
    min.data[Z==z,is.na(min.data[Z==z,][1,])] <- tmp[which.min(deltas),]
    max.data[Z==z,is.na(max.data[Z==z,][1,])] <- tmp[which.max(deltas),]
  }
  
  ## Store Results
  result <- c(glm(y ~ x1 + x2 + x3,family = "binomial",data=min.data)$coef[2],
                         glm(y ~ x1 + x2 + x3,family = "binomial",data=max.data)$coef[2])
  
  names(result) <- c("Lower Bound","Upper Bound")
  
  return(list(bounds = result,min.data = min.data,max.data = max.data,bet = bet))
}

# Score -------------------------------------------------------------------

g.logit = function(xx){exp(xx)/(1+exp(xx))}
dg.logit = function(xx){g.logit(xx)*(1-g.logit(xx))}

U.fun = function(bet,dat){
  yy = dat[,1]; xx.vec = cbind(1,dat[,-1])
  c(t(c(yy - g.logit(xx.vec%*%bet)))%*%xx.vec)/length(yy)
}

A.fun = function(bet,dat){
  yy = dat[,1]; xx.vec = cbind(1,dat[,-1])
  -t(c(dg.logit(xx.vec%*%bet))*xx.vec)%*%xx.vec/length(yy)
}


# Score and  ------------------------------------------------------

eta.sum = function(bet,dat){
  yy = dat[,1]; xx.vec = as.matrix(cbind(1,dat[,-1]))
  c(t(c(yy - g.logit(xx.vec%*%bet)))%*%xx.vec)
}

deta.sum = function(bet,dat){
  yy = dat[,1]; xx.vec = as.matrix(cbind(1,dat[,-1]))
  -t(c(dg.logit(xx.vec%*%bet))*xx.vec)%*%xx.vec
}

eta.new <- function(new.obs,bet){
  yy <- new.obs[1]; xx.vec <- t(as.matrix(c(1,new.obs[-1])))
  c(t( c(yy - g.logit(xx.vec%*%bet)) )%*%xx.vec)
}

deta.new <- function(new.obs,bet){
  yy <- new.obs[1]; xx.vec <- t(as.matrix(c(1,new.obs[-1])))
  -t(c(dg.logit(xx.vec%*%bet))*xx.vec)%*%xx.vec
}


# Influence Functions -----------------------------------------------------

inf.logistic <- function(new.obs,dat,bet){
  (eta.new(new.obs,bet)) %*% -solve(deta.sum(bet,dat))
}


# Plot Functions ----------------------------------------------------------

# theta.cc
# res.thetas
library(ggthemes)
library(gridExtra)
sim.coverage.plot <- function(res.thetas,theta.cc,truth = 1){
  
  plot.list = list()
    
    result <- data.frame(Estimate = theta.cc,LB = res.thetas[,1],UB = res.thetas[,2])
    
    result$informative <- factor(ifelse(result$LB * result$UB >0,"informative","non-informative"))
    
    ggplot(result, aes(x=1:nrow(result), y=Estimate, colour=informative)) + 
      geom_errorbar(aes(ymin=LB, ymax=UB), width=.1,alpha=.5) +
      geom_point() + theme_minimal() + scale_color_calc() + geom_hline(yintercept = 0) + geom_hline(yintercept = truth,lty = 2) + xlab(" ") +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
    
}

# sim.coverage.plot(res.thetas,theta.cc)

# Create another function but now it has the bounds and the CI around the bounds as well. Once I figure out the variance estimation thing
