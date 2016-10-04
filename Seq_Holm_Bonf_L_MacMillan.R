#Holm-Method Multi-hypothesis Testing
#Author: Laurel MacMillan
#Created: 11-1-2015

library(Rlab)

#Calculating lambda value for a vector of 1-d hypothesis
lambda_Bern <- function(thetaN, theta1, sigmaSq, x.mat, i, d){
  lambda.v <- vector() #Initialize vector to return
  for(j in 1:d){
    k <- sum(x.mat[,j])
    numObs <- length(x.mat[,j])
    #This is the formula for bernoulli distributed random variables
    lambda.v[j] <- k*log(theta1[j]/thetaN[j])+(numObs-k)*log((1-theta1[j])/(1-thetaN[j]))
  }
  return(lambda.v) #Returns a vector of lambda values
}
lambda_Norm <- function(thetaN, theta1, sigmaSq, x.mat, i, d){
  lambda.v <- vector() #Initialize vector to return
  for(j in 1:d){
    #This is the formula for normally distributed random variables
    lambda.v[j] <- (-1/(2*sigmaSq[d]))*(i*(theta1[j]^2-thetaN[j]^2)+2*(thetaN[j]-theta1[j])*sum(x.mat[,j]))
  }
  return(lambda.v) #Returns a vector of lambda values
}

#Creating a dx2 matrix of upper and lower bounds, with Holm stepdown procedure
create_Holm_Bounds <- function(alpha, beta, d){
  a <- vector()
  b <- vector()
  for(k in 1:d){
    a[k] <- log(1/(alpha/(d-k+1)))
    b[k] <- log(beta/k)
  }
  return(cbind(a, b)) #Returns an ordered Holm procedure set of boundaries
}

#Creating a dx2 matrix of upper and lower bounds with Bonferroni procedure
create_Bonf_Bounds <- function(alpha, beta, d){
  a <- rep(log(1/(alpha/d)), d)
  b <- rep(log((beta/d)), d)
  return(cbind(a, b)) #Returns a set of similar row values for each column (A, B) 
}

#Input a vector of lambda values (sorted) and check against the upper and lower bounds
check_Bounds <- function(lambda.v, bound.v, d){
  result <- vector()
  for(l in 1:d) {
    if(lambda.v[l,1] >= bound.v[l,1] || lambda.v[l,1] <= bound.v[l,2]){
      result[l] <- TRUE
    }
    else{
      result[l] <- FALSE
    }
  }
  return(result) #Returns a vector of results (1 if outside of bound, 0 if inside)
}

#Create a vector of decisions (whether the hypothesis is true or not based on true mu)
create_Decision <- function(theta0, thetaA, mu, d){
  decision <- vector()
  for(k in 1:d){
    decision[k] <- (mu[k]==theta0[k])
  }
  return(decision) #Returns a vector of 1 if the decision should be to accept null, 0 if rejecting null
}

#Create a vector of size d that checks whether the lambda values exceeded the correct boundaries
check_Decision <- function(lambda.m, ab.bounds, decision, d){
  err1 <- 0
  err2 <- 0
  for(x in 1:d){
    if(lambda.m[x,5]==1 && lambda.m[x,1]>0){ #If the decision is accepting null, checks lambda value
      err1 <- 1 #If lambda > 0, incorrect rejection of null, familywise type 1 error counter gets a 1
    }
    else if(lambda.m[x,5]==0 && lambda.m[x,1]<0){ #If the decision is rejecting null, checks lambda value
      err2 <- 1 #If lambda < 0, incorrect acceptance of null, familywise type 2 error counter gets a 1
    }
  }
  type1 <- rep(err1, d) #creates a vector for correct output dimension
  type2 <- rep(err2, d) #creates a vector for correct output dimension
  
  return(cbind(type1, type2)) #Returns a matrix of type 1 and type 2 error indicators
}

sim_study <- function(d, theta0, thetaA, mu, sigma, alpha, beta, type, num){
  
  n <- 2000 #Set number of possible iterations to simulate a sequential data set
  numNorm <- num[1]
  numBern <- num[2]
  for(t in 1:numNorm){
    norm.mat<-mapply(function(x,y){rnorm(x,y,n=n)},x=mu[1:numNorm],y=sigma[1:numNorm])
  }
  for(v in 1:numBern){
    bern.mat<-mapply(function(x,y){rbern(n=n,x)},x=mu[(numNorm+1):(numNorm+numBern)]) 
  }
  
  #Create two sets of bounds
  holm.bounds <- create_Holm_Bounds(alpha, beta, d)  
  bonf.bounds <- create_Bonf_Bounds(alpha, beta, d)    
  
  error.dec <- create_Decision(theta0, thetaA, mu, d)
  #For each additional data point until decision is reached...
  for(i in 1:n){
    lambda.vecH1 <- lambda_Norm(theta0[1:numNorm], thetaA[1:numNorm], sigma[1:numNorm], head(norm.mat, i) ,i,numNorm)
    lambda.vecH2 <- lambda_Bern(theta0[(numNorm+1):(numNorm+numBern)], thetaA[(numNorm+1):(numNorm+numBern)], sigma, head(bern.mat, i) ,i,numBern) #calculate lambda value with added Xi
    lambda.vec.H <- c(lambda.vecH1, lambda.vecH2)
    lambda.mH <- cbind(lambda.vec.H, theta0, thetaA, mu, error.dec) #combine lambda with hypothesis information
    
    #Sort matrix by descending lambda values (not necessary for Bonferroni method but happens regardless)
    lambda.mH <- lambda.mH[order(lambda.mH[,1],decreasing=TRUE),]
    rH <- check_Bounds(lambda.mH, holm.bounds, d)
    if(all(rH)){ #When this happens, the sequential data collection will stop
      resH <- check_Decision(lambda.mH, holm.bounds, error.dec, d) #Did tests conclude correct decisions?
      stopH <- i
      break
    }
  }
  
  for(i in 1:n){
    lambda.vecB1 <- lambda_Norm(theta0[1:numNorm], thetaA[1:numNorm], sigma[1:numNorm], head(norm.mat, i) ,i,numNorm)
    lambda.vecB2 <- lambda_Bern(theta0[(numNorm+1):(numNorm+numBern)], thetaA[(numNorm+1):(numNorm+numBern)], sigma, head(bern.mat, i) ,i,numBern) #calculate lambda value with added Xi
    lambda.vec.B <- c(lambda.vecB1, lambda.vecB2)
    lambda.mB <- cbind(lambda.vec.B, theta0, thetaA, mu, error.dec) #combine lambda with hypothesis information
    #Sort matrix by descending lambda values (not necessary for Bonferroni method but happens regardless)
    lambda.mB <- lambda.mB[order(lambda.mB[,1],decreasing=TRUE),]   
    rB <- check_Bounds(lambda.mB, bonf.bounds, d)
    if(all(rB)){ #When this happens, the sequential data collection will stop
      resB <- check_Decision(lambda.mB, bonf.bounds, error.dec, d) #Did tests conclude correct decisions?
      stopB <- i
      break
    }
  }
  
  return(rbind(cbind(stopH, lambda.mH, holm.bounds, resH), cbind(stopB, lambda.mB, bonf.bounds, resB))) #Returns all of the testing information, including decision
}

#Input all of the same information as the sim_study function above
monte_study <- function(dTr, th0, thA, muTr, ssTr, aTr, bTr, typeTest, numTests){
  iter <- 60000 #Number of simulations
  fwH.1 <- vector()
  fwH.2 <- vector()
  fwB.1 <- vector()
  fwB.2 <- vector()
  avgH <- vector()
  avgB <- vector()
  comp <- dTr+1
  #Run a Bonferroni or Holm method simulated sequential test
  for(m in 1: iter){
    result <- sim_study(dTr, th0, thA, muTr, ssTr, aTr, bTr, typeTest, numTests)
    avgH[m] <- result[1,1] #Vector of stopping times
    avgB[m] <- result[comp,1]
    fwH.1[m] <- result[1,9] #Vector of family-wise type I error rate indicators
    fwH.2[m] <- result[1,10] #Vector of family-wise type II error rate indicators
    fwB.1[m] <- result[comp,9]
    fwB.2[m] <- result[comp,10]
  }
  fwH.error1 <- round(sum(fwH.1)/iter, 4) #Calculate type I error rate
  fwH.error2 <- round(sum(fwH.2)/iter, 4) #Calculate type II error rate
  fwB.error1 <- round(sum(fwB.1)/iter, 4) #Calculate type I error rate
  fwB.error2 <- round(sum(fwB.2)/iter, 4) #Calculate type II error rate
  return(rbind(c(mean(avgB), fwB.error1, fwB.error2),
               c(mean(avgH), fwH.error1, fwH.error2))) #Returns average stopping time and family-wise error rates
}

#Parameters in monte carlo function
# D-number of hypothesis tests
# {Null1, Null2, ... , NullD} Vector of null mu/p
# {Alt1, Alt2, ... , AltD} Vector of alternative mu/p
# {mu1, mu2, ... , muD} Vector of true mu/p values
# {sigma1, sigma2, ... , sigmaD} Vector of variance values
# alpha
# beta
# Vector of normal or bernoulli trial
# Number of how many normal trials and how many bernoulli trials
# I.E. c(2, 1) if testing 2 normal and 1 bernoulli hypotheses

monte_study(2,c(0,0.1),c(0.25,0.2),c(0,0.1),rep(1,2),.05,.10, c("n", "b"), c(1,1))
monte_study(2,c(0,0.1),c(0.25,0.2),c(0.25,0.1),rep(1,2),.05,.10, c("n", "b"), c(1,1))
monte_study(2,c(0,0.1),c(0.25,0.2),c(0.25,0.2),rep(1,2),.05,.10, c("n", "b"), c(1,1))

monte_study(3,c(0,0,0.1),c(0.25, 0.25, 0.2),c(0, 0, 0.1),rep(1,3),.05,.10, c("n", "n", "b"), c(2,1))
monte_study(3,c(0,0,0.1),c(0.25, 0.25, 0.2),c(0.25, 0, 0.1),rep(1,3),.05,.10, c("n", "n", "b"), c(2,1))
monte_study(3,c(0,0,0.1),c(0.25, 0.25, 0.2),c(0.25, 0.25, 0.1),rep(1,3),.05,.10, c("n", "n", "b"), c(2,1))
monte_study(3,c(0,0,0.1),c(0.25, 0.25, 0.2),c(0.25, 0.25, 0.2),rep(1,3),.05,.10, c("n", "n", "b"), c(2,1))
monte_study(3,c(0,0,0.1),c(0.25, 0.25, 0.2),c(0, 0, 0.2),rep(1,3),.05,.10, c("n", "n", "b"), c(2,1))

monte_study(4,c(0,0,0.1,0.1),c(0.25, 0.25, 0.2, 0.2),c(0, 0, 0.1, 0.1),rep(1,4),.05,.10, c("n", "n", "b", "b"), c(2,2))
monte_study(4,c(0,0,0.1,0.1),c(0.25, 0.25, 0.2, 0.2),c(0.25, 0, 0.1, 0.1),rep(1,4),.05,.10, c("n", "n", "b", "b"), c(2,2))
monte_study(4,c(0,0,0.1,0.1),c(0.25, 0.25, 0.2, 0.2),c(0.25, 0.25, 0.1, 0.1),rep(1,4),.05,.10, c("n", "n", "b", "b"), c(2,2))
monte_study(4,c(0,0,0.1,0.1),c(0.25, 0.25, 0.2, 0.2),c(0, 0, 0.2, 0.2),rep(1,4),.05,.10, c("n", "n", "b", "b"), c(2,2))
monte_study(4,c(0,0,0.1,0.1),c(0.25, 0.25, 0.2, 0.2),c(0.25, 0, 0.2, 0.2),rep(1,4),.05,.10, c("n", "n", "b", "b"), c(2,2))
monte_study(4,c(0,0,0.1,0.1),c(0.25, 0.25, 0.2, 0.2),c(0.25, 0.25, 0.2, 0.2),rep(1,4),.05,.10, c("n", "n", "b", "b"), c(2,2))

monte_study(2,c(0.1,0.1),c(0.2,0.2),c(0.2,0.1),rep(1,2),.05,.10)
monte_study(2,c(0.1,0.1),c(0.2,0.2),c(0.2,0.2),rep(1,2),.05,.10)

monte_study(3, c(0.1,0.1,0.1), c(0.2,0.2,0.2), c(0.1,0.1,0.1), c(1,1,1), 0.05, 0.10)
monte_study(3, c(0.1,0.1,0.1), c(0.2,0.2,0.2), c(0.2,0.1,0.1), c(1,1,1), 0.05, 0.10)
monte_study(3, c(0.1,0.1,0.1), c(0.2,0.2,0.2), c(0.2,0.2,0.1), c(1,1,1), 0.05, 0.10)
monte_study(3, c(0.1,0.1,0.1), c(0.2,0.2,0.2), c(0.2,0.2,0.2), c(1,1,1), 0.05, 0.10)

monte_study(4, c(0.1,0.1,0.1,0.1), c(0.2,0.2,0.2,0.2), c(0.1,0.1,0.1,0.1), c(1,1,1,1), 0.05, 0.10)
monte_study(4, c(0.1,0.1,0.1,0.1), c(0.2,0.2,0.2,0.2), c(0.2,0.1,0.1,0.1), c(1,1,1,1), 0.05, 0.10)
monte_study(4, c(0.1,0.1,0.1,0.1), c(0.2,0.2,0.2,0.2), c(0.2,0.2,0.1,0.1), c(1,1,1,1), 0.05, 0.10)
monte_study(4, c(0.1,0.1,0.1,0.1), c(0.2,0.2,0.2,0.2), c(0.2,0.2,0.2,0.1), c(1,1,1,1), 0.05, 0.10)
monte_study(4, c(0.1,0.1,0.1,0.1), c(0.2,0.2,0.2,0.2), c(0.2,0.2,0.2,0.2), c(1,1,1,1), 0.05, 0.10)
