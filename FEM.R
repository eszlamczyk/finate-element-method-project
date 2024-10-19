calcXi <- function(h,i){
  #creates point xi
  return(h*i)
}


calcEi <- function(x,h,i){
  #calculates value for given x and i with spacing of h
  xi <- calcXi(h,i)
  xiMin <- calcXi(h,i-1)
  xiMax <- calcXi(h,i+1)
  if(x > xiMax || x < xiMin){
    return (0)
  }
  else{
    if (x < xi){
      return ((x-xiMin)/h)
    } else{
     return ((xiMax-x)/h)
    }
  }
}

calcEiPrime <- function(x,h,i){
  #the same that calcEiPrime but calculates the derivative instead
  xi <- calcXi(h,i)
  xiMin <- calcXi(h,i-1)
  xiMax <- calcXi(h,i+1)
  if(x > xiMax || x < xiMin){
    return (0)
  }
  else{
    if (x < xi){
      return (1/h)
    } else{
      return (-1/h)
    }
  }
}

calcF <- function(x,G){
  #calculates the G*4pi*rho function
  if(x < 1) return (0)
  else{
    if(x > 2) 
      return (0)
    else
      return (G*4*3.14)
  }
}

integrateL <- function(a,b,h,i){
  #integrates ei dx by given bounds and i index
  negativeScalar <- (b-a)/2
  positiveScalar <- (a+b)/2
  return(negativeScalar*
      (
      calcEi((negativeScalar/(sqrt(3)) + positiveScalar),h,i) 
      +
      calcEi((negativeScalar/(-1*sqrt(3)) + positiveScalar),h,i)
      )
    )
}

integrateB <- function(a,b,h,i,j){
  #integrates ei'*ej' dx for given i and j using Gauss-Laplace
  #integration method
  negativeScalar <- (b-a)/2
  positiveScalar <- (a+b)/2
  return(
    negativeScalar*
    (
      (
      calcEiPrime((negativeScalar/(sqrt(3)) + positiveScalar),h,i)*
      calcEiPrime((negativeScalar/(sqrt(3)) + positiveScalar),h,j)
      )  
      +
      (  
      calcEiPrime((negativeScalar/(-1*sqrt(3)) + positiveScalar),h,i)*
      calcEiPrime((negativeScalar/(-1*sqrt(3)) + positiveScalar),h,j)
      )
    )
  )
}

calcB <- function(h,i,j,a,b){
  #creates proper bounds for given i and j
  if (i == j){
    lower <- max(a,calcXi(h,i-1))
    upper <- min(b,calcXi(h,i+1))
  } else {
    lower <- min(calcXi(h,i),calcXi(h,j))
    upper <- max(calcXi(h,i),calcXi(h,j))
  }
  return(integrateB(lower,upper,h,i,j))
}


calcL <- function(G,h,i){
  #creates bounds for given i
  xi <- calcXi(h,i)
  xiMin <- max(1,calcXi(h,i-1))
  xiMax <- min(2,calcXi(h,i+1))
  f = calcF(xi,G)
  return(4*pi*G*integrateL(xiMin,xiMax,h,i))
}

createLVector <- function(G,a,b,n){
  #creates left side vector of Li elements,
  #distribution of values looks correct, can't determine if 
  #they are proper values tho
  l <- b-a
  h <- l/n
  L <- c()
  for(i in 1:(n-1)){
    L <- c(L,calcL(G,h,i))
  }
  return(L)
}

createBMatrix <- function(a,b,n){
  #creates matrix of Bij elements
  l <- b-a
  h <- l/n
  B <- matrix(nrow = n-1, ncol = n-1)
  for(i in 1:(n-1)){
    for(j in 1:(n-1)){
      if(i == 1 && j == 1) {
        B[i,j] <- 1
      }else if(i == 1) {
        B[i,j] <- 0
      }else if(j == 1) {
        B[i,j] <- 0
      }else {
        B[i,j] <- (-1)*calcB(h,i,j,a,b)
      }
    }
  }
  return(B)
}

calcUroof <- function(a,b,u1,u2,x){
  #function used to fix Dirichlet boundary condition
  A <- (u2-u1)/b-a
  B <- u1 - A*a
  return(A*x+B)
}

#main method
differentialSolver <- function(n,a,b,u1,u2,G){
  #useful values
  l <- b-a
  h <- l/n
  x <- c()
  #create column of xi points for later function plot
  for (i in 1:(n-1)){
    x <- c(x,calcXi(h,i))
  }
  #get L
  L <- createLVector(G,a,b,n)
  #get B
  B <- createBMatrix(a,b,n)
  #solve Bw = L to get weights
  weight <- solve(B,L)
  #u = w + u_roof
  y <- c()
  for (i in 1:(n-1)){
    y <- c(y,weight[i]+calcUroof(a,b,u1,u2,x[i]))
  }
  
  plot(x,y, type="l",col="blue",ylab="Φ(x)", )#xlim=c(-0.2,3.2), ylim=c(3.8,5.2))
  #line(x,calcUroof(a,b,u1,u2,x), type="l", col="red")
  #line(x,weight, type="l", col="green")
  #legend("bottomleft",legend = c("fi", "u z daszkiem", "w"), col = c("blue","red","green"), lty=1,cex=0.8)
  
}
