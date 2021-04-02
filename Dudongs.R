
mu <- function(x,alpha,beta,gamma){
  return (alpha-beta*gamma^x)
}

dudongs <- function(N,x,Y,prop.sd=c(0.3,0.2,0.2,0.2),nchain=10^4){
  
  init <- c(1,1,1,0.9)
  tau0 <- 10^-6
  chain <- matrix(NA,nchain+1,4)
  colnames(chain) <- c('alpha','beta','tau','gamma')
  chain[1,] <- init
  acc_rates <- rep(0,4)
  
  
  for (iter in 1:nchain){
    
    current <- chain[iter,]
    
  ##Mise a jour de alpha
    
    prop <- current
    prop[1] <- rlnorm(1,log(current[1]),prop.sd[1])
    top <- -(prop[1]^2*tau0*0.5 +0.5*sum((Y-mu(x,prop[1],prop[2],
                                            prop[4]))^2)*prop[3])

    bottom <- -(current[1]^2*tau0*0.5 +0.5*sum((Y-mu(x,current[1],current[2],
                                               current[4]))^2)*current[3])
    
    prop.ratio <- prop[1] / current[1]
    acc_prob <- exp(top-bottom)*prop.ratio
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[1] <- acc_rates[1] + 1
    }
    
    #mise a jour de beta
    prop <- current
    prop[2] <- rlnorm(1,log(current[2]),prop.sd[2])
    top <- -(prop[2]^2*tau0*0.5 +0.5*sum((Y-mu(x,prop[1],prop[2],
                                            prop[4]))^2)*prop[3])
    
    bottom <- -(current[2]^2*tau0*0.5 +0.5*sum((Y-mu(x,current[1],current[2],
                                                  current[4]))^2)*current[3])
    
    prop.ratio <- prop[2] / current[2]
    acc_prob <- exp(top-bottom)*prop.ratio
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[2] <- acc_rates[2] + 1
    }
    
    ##mise a jour de gamma
    prop <- current
    prop[4] <- rlnorm(1,log(current[4]),prop.sd[4])
    top <- -(0.5*sum((Y-mu(x,prop[1],prop[2],
                             prop[4]))^2)*prop[3])
    
    bottom <- -(0.5*sum((Y-mu(x,current[1],current[2],
                                current[4]))^2)*current[3])
    
    prop.ratio <- prop[4] / current[4]
    acc_prob <- exp(top-bottom)*prop.ratio
    if (runif(1) < acc_prob){
      current <- prop
      acc_rates[4] <- acc_rates[4] + 1
    }
    
    ##mise a jour de tau
    updated_shape <- 10^(-3) + N/2
    updated_rate <- 10^(-3)+ 0.5*sum((Y-mu(x,current[1],current[2],
                                        current[4]))^2)
    tauest <- rgamma(1,shape=updated_shape,rate=updated_rate)
    current[3] <- tauest
   
    ## 
    chain[iter+1,] <- current
    
  }
  return(list(chain = chain, acc_rates = acc_rates / nchain))
}

##Application

"N" <- 27  
"x" <-
  c(1, 1.5, 1.5, 1.5, 2.5, 4, 5, 5, 7, 8, 8.5, 9, 9.5, 9.5, 10, 
    12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29, 31.5)
"Y" <-
  c(1.8, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 
    2.26, 2.4, 2.39, 2.41, 2.5, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 
    2.47, 2.64, 2.56, 2.7, 2.72, 2.57)

out <- dudongs(N,x,Y,prop.sd=c(0.005,0.05,0.1,0.015))
out$chain <- out$chain[-(1:2000),]
plot(out$chain[,4], type = "l", main = "")
print(out$acc_rates)
print(colMeans(out$chain))
