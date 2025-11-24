\name{selector}
\alias{selector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Frequency Slector
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It selects the eta frequency components with the largest leading eigenvalue.
}
\usage{
selector(fxx_grad_evec, f_xx, len_freq, eta, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fxx_grad_evec}{
%%     ~~Describe \code{fxx_grad_evec} here~~
The estimated principal component loadings over the fundamental frequencies.
}
  \item{f_xx}{
%%     ~~Describe \code{f_xx} here~~
An estimate of the spectral density matrix.
}
  \item{len_freq}{
%%     ~~Describe \code{len_freq} here~~
len_freq = n/2 is the number of fundamental frequency components.
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
Number of frequency compoennts to be selected.
}
  \item{p}{
%%     ~~Describe \code{p} here~~
Dimension
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
\item{LSDPCA}{The localized and sparse eigenvector series}
\item{best_sol}{0-1 valued vector representing the selected frequencies}
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


# Population Spectral density matrix of the Underlying AR(2) process
p <- 10
n <- 1024
omega <- seq(0,.5, length.out = n/2)
f_xx <- array(0,dim=c(p,p,length(omega)))
f_x_omega <- 9/(2.81 - 3.8*cos(2*pi*omega) +(1.8)*cos(4*pi*omega))
f_e <- 1
for(ell in 1:length(omega)){
  f_xx[,,ell] <- diag(p)
  f_xx[1:2,1:2,ell] <- matrix(c(f_x_omega[ell],f_x_omega[ell], f_x_omega[ell], f_x_omega[ell]+f_e), nrow=2, byrow = TRUE)
}

# First eigenvector
f_evec11 <- matrix(0, nrow=p, ncol = length(omega))
for(ell in 1:length(omega)){
  f_evec_ell <- eigen(f_xx[,,ell])$vectors[,1]
  f_evec11[,ell] <- f_evec_ell
  gc()
}


# Data
p <- 10
n <- 1024


phi_1 <- 1
phi_2 <- -.9


Xt1 <- arima.sim(list(order=c(2,0,0), ar=c(phi_1,phi_2)), n = n)*3
Xt2 <- Xt1 + rnorm(n)
X_omega_2 <- cbind(Xt1,Xt2)
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
rm(F_tp_list)
gc()


#LSPCA
Ex1 <- LSPCA(n,p, f_xx1, Pi <- matrix(0, nrow = 2*p, ncol=2*p), d=1, lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
                                               control = list(fan_maxinc = 10, verbose = 0), eta=500, s=5, n_iter = 20, nu=0)


## Frequency selection
# eta = 500
LSPCA500 <- selector(Ex1[[1]], f_xx1, n/2,500,p)

plot(abs(LSPCA500[[1]][1,]), main = "1st coordinate - Mod 1st sample evec over freq - eta = 500", pch=18, ylab="", ylim = c(0,1), type="l")
points(Mod(f_evec11[1,]), pch=18, ylab="", col = "green", type="l")
legend("topright", inset=.002, legend=c("Estimated", "Population"),
       col=c("black", "green"), lty=c(1,1))


# eta = 100
LSPCA100 <- selector(Ex1[[1]], f_xx1, n/2,100,p)

plot(abs(LSPCA100[[1]][1,]), main = "1st coordinate - Mod 1st sample evec over freq - eta = 100", pch=18, ylab="", ylim = c(0,1), type="l")
points(Mod(f_evec11[1,]), pch=18, ylab="", col = "green", type="l")

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
