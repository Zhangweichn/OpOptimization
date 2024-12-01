\name{Optimization}
\alias{Optimization}
\title{Optimization}
\usage{
result<-optimization(S0,N,alpha,num_w,num_p,C,K,rf,T,Mu,Sigm)
}
\description{
Prints 'An Optimization Method of Buy-write Strategy Based on Linear Programming'
}
\arguments{
  \item{S0}{An array recording the initial prices of different assets(eg:S0 <- c(20804.11,5969.34,2406.67))}
  \item{N}{Number of Brownian motion predictions}
  \item{ST}{The end-of-period price of each asset is predicted N times by Brownian motion, with num_w rows and N columns}
  \item{alpha}{alpha}
  \item{num_w}{Number of assets}
  \item{num_p}{The number of options on each asset, an array of length num_w}
  \item{C}{The current market price of the option, with num_w rows and sum_p columns}
  \item{K}{Strike price for each asset, with num_w rows and sum_p columns}
  \item{rf}{Risk-free rate of return}
  \item{T}{Duration of asset price change research, using days/365}
  \item{Mu}{Each asset Mu in the Brownian motion, an array of length num_w}
  \item{Sigm}{The volatility of each asset in the Brownian motion, an array of length num_w}
}
\examples{
}

