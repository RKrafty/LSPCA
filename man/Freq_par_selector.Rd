\name{Freq_par_selector}
\alias{Freq_par_selector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Frequency parameter selector
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Selects the localization (frequency) parameter
}
\usage{
Freq_par_selector(X, f_xx1 = NULL, fxx_evec_s, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
nxp data matrix
}
  \item{f_xx1}{
%%     ~~Describe \code{f_xx1} here~~
estimated spectral density matrix. If set to NULL, it will estimated it by the multitapper estimator using 10 sine tapper functions.
}
  \item{fxx_evec_s}{
%%     ~~Describe \code{fxx_evec_s} here~~
The sparse eigenvector series estimate
}
  \item{d}{
%%     ~~Describe \code{d} here~~
Dinesion of the principal subspace
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{freq_par_AIC}{frequency parameter selected by AIC}
\item{freq_par_AICc}{frequency parameter selected by AICc}
\item{freq_par_BIC}{frequency parameter selected by BIC}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Jamshid Namdari
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(dplR)
library(mvtnorm)
library(dplR)
library(waveslim)
library(astsa)

# Data
p <- 64
n <- 1024


phi_1 <- 1
phi_2 <- -.9


Xt10 <- arima.sim(list(order=c(2,0,0), ar=c(phi_1,phi_2)), n = n)*3
Xt1 <- pass.filt(Xt10, W=c(0.1, 0.2), type="pass", method = "Butterworth")+rnorm(n)
Xt2 <- Xt1 + rnorm(n)
Xt3 <- 1.1*Xt1 + rnorm(n)
Xt4 <- 1.2*Xt1 + rnorm(n)
Xt5 <- 1.15*Xt1 + rnorm(n)

X_omega_2 <- cbind(Xt1,Xt2,Xt3,Xt4,Xt5)
X_wn <- rmvnorm(n, sigma = diag(p-(ncol(X_omega_2))))
X1 <- cbind(X_omega_2, X_wn)

## Multitapper estimate of the spectral denisty matrix
U <- sine.taper(n,20)

X_tp <- apply(U, MARGIN = 2, function(u) u*X1, simplify = FALSE)
F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

len_freq <- n/2
F_tp1 <- array(0, c(p, p, len_freq))
for (ell in 1:len_freq) {
  for(j in 1:length(F_tp_list)){
    F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
  }
  F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
}
plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral density")
for(a in 2:p){
  lines(Re(F_tp1[a,a,])*n, col= a)
}
f_xx1 <- F_tp1*n
rm(U)
rm(X_tp)
#rm(X)
rm(F_tp_list)
gc()


Ex1 <- LSPCA(n,p, f_xx1, Pi <- matrix(0, nrow = 2*p, ncol=2*p), d=1, lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
                                               control = list(fan_maxinc = 10, verbose = 0), eta=500, s=5, n_iter = 20, nu=0)


eta_hat <- Freq_par_selector(X1,f_xx1,Ex1[[1]], d=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
