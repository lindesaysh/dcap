#' This function converts parameters of N(\mu, \sigma^2) to parameters of lognormal (\exp{\mu, \sigma^2})
#'
#' @param thetan list object containting the mean ({\code mu}) and the variance ({\code sigma2}) of a normal distribution.
#'
#' @return
#' list object containing {\code muln} which is the mean of the log normal and {\code sigma2ln} which is the variance of the log normal.
#'
#' @author Lindesay Scott-Hayward
#'
fromNorm2logNorm<- function(thetan){

  mu = thetan$mu
  sigma2 = thetan$sigma2
  mux = exp(mu+sigma2/2)
  sigma2x = (exp(sigma2)-1)*exp(2*mu+sigma2)

  return(list(muln=mux, sigma2ln=sigma2x))

}



#' This function evaluates a Gaussian kernel.
#'
#' @param x single value or vector of point to evaluate kernel at.
#' @param X single value or vector of second point to evaluate kernel at.
#' @param h single value or vector which relates to the effective radius of the kernel.
#'
#' @details
#' $$\frac{1}{\sqrt(2*\pi)}\exp(\frac{-\frac{x-X}{h^2}}{2}$$
#'
#' @return
#' number or vector of evaluated points.
#'
#' @examples
#' x<-seq(1,10, length=100)
#' X<-5
#' h<-2
#'
#' kern<-evaluateGaussianKernel(x, X, h)
#' plot(x, kern, type='l')
#'
#' @author Lindesay Scott-Hayward
#'
#'
evaluateGaussianKernel<- function(x, X, h){

  return((1/(sqrt(2*pi)))*exp((-((x-X)/h)^2)/2))

}


#' Autocorrelation function for ordered data.
#'
#' This function is used to find the autocorrelation of different rows of data. These are then divided by three to find sd (sigma) bandwidth of kernel to which independence occurs
#'
#' @param y vector of densities to find autocorrelation.
#'
#' @author Lindesay Scott-Hayward
#'
#'
evaluateAutocorrelation<-function(y){

  lag<- seq(1:length(y)-1)
  out<- length(lag)
  #flag
  pos<- TRUE
  #counter
  k<- 1
  while(pos){
    top<- vector(length=(length(y)-lag[k]))
    bot<- vector(length=length(y))
    for(i in 1:(length(y)-lag[k])){
      top[i]<- (y[i]-mean(y))*(y[(i+lag[k])]-mean(y))
    }
    for(a in 1:length(y)){
      bot[a]<- (y[a]-mean(y))^2
    }
    out[k]<- sum(top)/sum(bot)
    if(out[k]<0){
      pos=FALSE
    }
    k=k+1
  }

  return(c(k-2, out[k-2]))

}


