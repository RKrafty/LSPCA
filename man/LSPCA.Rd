\name{LSPCA}
\alias{LSPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
    Localized Sparse Principal Component Analysis
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  LSPCA(n, p, f_xx, Pi, d, lambda, lr = 0.02, maxiter = 60, control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    %%     ~~Describe \code{n} here~~
      Sample size
  }
  \item{p}{
    %%     ~~Describe \code{p} here~~
      Dimension of the time series
  }
  \item{f_xx}{
    %%     ~~Describe \code{f_xx} here~~
      Estimated spectral density matrix
  }
  \item{d}{
    %%     ~~Describe \code{d} here~~
      Dimension of the principal subspace
  }
  \item{lambda}{
    %%     ~~Describe \code{lambda} here~~
      The sparsity parameter of the ADMM step for finding the initial estimate of the principal subspace. Suggested velue is \eqn{0.5\times\sqrt{(\log(p) / n)}}}.
  \item{lr}{
    %%     ~~Describe \code{lr} here~~
  }
  \item{maxiter}{
    %%     ~~Describe \code{maxiter} here~~
      Maximum number of iterations of the ADMM step
  }
  \item{control}{
    %%     ~~Describe \code{control} here~~
      Additional parameters to control the optimization process of the ADMM step
  }
  \item{eta}{
    %%     ~~Describe \code{eta} here~~
      Frequency parameter
  }
  \item{s}{
    %%     ~~Describe \code{s} here~~
      Sparsity level
  }
  \item{n_iter}{
    %%     ~~Describe \code{n_iter} here~~
      Number of iterations of the SOAP step
  }
  \item{nu}{
    %%     ~~Describe \code{nu} here~~
      Smoothing parameter
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
  \item{Evec_series_k}{The k-th eigenvector series, k=1,2,...}
  \item{LSDPCA}{The localized and sparse eigenvector series}
  \item{best_sol}{0-1 valued vector representing the selected frequencies}
  \item{f_deflated}{Deflated spectral denisty matrix (only returned when d = 1)}
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

  ##########################################################################################
  ## Example 1 : This example reproduces Bottom left panel of Figure 2 in the manuscript
  ##########################################################################################

set.seed(1000*3+12)

####################
## Data Generation
####################
p <- 64
n <- 1024
phi_1 <- 1.5
phi_2 <- -.75
a1 <- 1/20
a2 <- -1/1.15
omega <- seq(0,.5, length.out = n/2)
f_xx <- array(0,dim=c(p,p,length(omega)))
c1 <- 3
f_x_omega0 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
phi_1 <- 1.5 +.05
f20 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
phi_1 <- 1.5 -.05
f30 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
phi_1 <- 1.5 +.15
f40 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
phi_1 <- 1.5 -.15
f50 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec

f_x_omega <- f_x_omega0
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f_x_omega[ell] <- 0
}
f2 <- f20
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f2[ell] <- 0
}
f3 <- f30
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f3[ell] <- 0
}
f4 <- f40
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f4[ell] <- 0
}
f5 <- f50
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f5[ell] <- 0
}

f_e <- 1/4
for(ell in 1:length(omega)){
  f_xx[,,ell] <- diag(f_e, p)
  f_xx[1:5,1:5,ell] <- matrix(c(f_x_omega[ell],(2*1.1)*f_x_omega[ell],1.2*f_x_omega[ell],1.25*f_x_omega[ell], 3*.75*f_x_omega[ell],
                                2*1.1*f_x_omega[ell], 4*1.1*1.1*f_x_omega[ell]+f_e+1*f2[ell], 2*1.1*1.2*f_x_omega[ell],2*1.1*1.25*f_x_omega[ell],2*1.1*3*.75*f_x_omega[ell],
                                1.2*f_x_omega[ell], 2*1.1*1.2*f_x_omega[ell], 1.2*1.2*f_x_omega[ell]+f_e+f3[ell], 1.2*1.25*f_x_omega[ell],1.2*3*.75*f_x_omega[ell],
                                1.25*f_x_omega[ell],2*1.1*1.25*f_x_omega[ell],1.2*1.25*f_x_omega[ell],1.25*1.25*f_x_omega[ell]+f_e+1*f4[ell],1.25*3*.75*f_x_omega[ell],
                                3*.75*f_x_omega[ell],2*1.1*3*.75*f_x_omega[ell],1.2*3*.75*f_x_omega[ell],1.25*3*.75*f_x_omega[ell],9*.75*.75*f_x_omega[ell]+f_e+9*f5[ell]), nrow=5, byrow = TRUE)
}



f_evec11 <- matrix(0, nrow=p, ncol = length(omega))
for(ell in 1:length(omega)){
  f_evec_ell <- eigen(f_xx[,,ell])$vectors[,1]
  f_evec11[,ell] <- f_evec_ell
  gc()
}


len_freq <- n/2
Coord10 <- matrix(0,nrow = 5, ncol = len_freq)
Coord20 <- matrix(0,nrow = 5, ncol = len_freq)
Coord30 <- matrix(0,nrow = 5, ncol = len_freq)
Coord40 <- matrix(0,nrow = 5, ncol = len_freq)
Coord50 <- matrix(0,nrow = 5, ncol = len_freq)

Coord1 <- matrix(0,nrow = 5, ncol = len_freq)
Coord2 <- matrix(0,nrow = 5, ncol = len_freq)
Coord3 <- matrix(0,nrow = 5, ncol = len_freq)
Coord4 <- matrix(0,nrow = 5, ncol = len_freq)
Coord5 <- matrix(0,nrow = 5, ncol = len_freq)


phi_1 <- 1.5
phi_2 <- -.75
a1 <- 1/20
a2 <- -1/1.15


Xt10 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
Xt1 <- pass.filt(Xt10, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp1 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

phi_1 <- 1.5 +.05
yt20 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
yt2 <- pass.filt(yt20, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp2 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

phi_1 <- 1.5 -.05
yt30 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
yt3 <- pass.filt(yt30, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp3 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

phi_1 <- 1.5 +.15
yt40 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
yt4 <- pass.filt(yt40, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp4 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

phi_1 <- 1.5 -.15
yt50 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
yt5 <- pass.filt(yt50, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp5 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

Xt2 <- 2*1.1*Xt1 + 1*yt2 + rnorm(n,0,1/2)
Xt3 <- 1.2*Xt1 + yt3 + rnorm(n,0,1/2)
Xt4 <- 1.25*Xt1 + 1*yt4 + rnorm(n,0,1/2)
Xt5 <- 3*0.75*Xt1 + 3*yt5 + rnorm(n,0,1/2)

X_omega_2 <- cbind(Xt1,Xt2,Xt3,Xt4,Xt5)
X_wn <- rmvnorm(n, sigma = diag(p-(ncol(X_omega_2))))
X_original <- cbind(X_omega_2, X_wn)
#shuffled_index <- sample.int(p)
shuffled_index <- 1:p
X1 <- X_original[,shuffled_index]


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
f_xx1 <- F_tp1*n
rm(U)
rm(X_tp)
rm(F_tp_list)
gc()

plot(omega, f_xx1[1,1,1:512], type="l", ylim = c(0,30), ylab="", xlab="Frequency")
for (k in 2:64) {
  points(omega, f_xx1[k,k,1:512], type="l", col=k)
}
abline(v=c(0.05,0.25))


##################################################
## Estimate of 1-dimensional principal subspaces
##################################################
nu_v <- c(0,.2,.4,.6,.8,1)
for(k in 1){
  nu <- nu_v[k]

  LSDPCA_ADMM_SOAP_Ex1 <- LSPCA(n,p, f_xx1, lambda = 0.5 * sqrt(log(p) / n), d=1, lr = 0.02, maxiter = 60,
                                                 control = list(fan_maxinc = 10, verbose = 0), eta=200, s=5, n_iter = 20, nu=nu_v[k])


  gc()
}


f_MT_evec11 <- matrix(0, nrow=p, ncol = length(omega))
for(ell in 1:length(omega)){
  f_MT_evec_ell <- eigen(f_xx1[,,ell])$vectors[,1]
  f_MT_evec11[,ell] <- f_MT_evec_ell
  gc()
}


###########################
## Plots
###########################

library(latex2exp)


xlab <- "Hz"
ylab <- "Coordinate"
legend_title = ""
asp = .5 #0.2
bar_height = 10/2
font_size = 20/2



Localized_Est <- selector(LSDPCA_ADMM_SOAP_Ex1[[1]],f_xx1,n/2,(2/5)*512,p)
evecs <- t(Mod(Localized_Est[[1]]))
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors, limits = c(0, 1)) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = TeX("$\\theta=0$")) +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), labels=function(x) x*64/n)



  ###################
  ## Extra Examples
  ###################

  ##########################################################
  ## Extra Example 1: with 1-dimensional principal subspace
  ##########################################################

  # Data
  p <- 10
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
  rm(F_tp_list)
  gc()

  #LSPCA
  Ex1 <- LSPCA(n,p, f_xx1, d=1, lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
               control = list(fan_maxinc = 10, verbose = 0), eta=500, s=5, n_iter = 20, nu=0)


  # Plot of the first coordinate of the first PC loadings
  plot(Mod(Ex1[[1]][1,]), type="l")

  ##########################################################
  ## Extra Example 2: with 2-dimensional principal subspace
  ##########################################################
  set.seed(1234)
  # Generating data

  p <- 64
  n <- 1024

  phi_1 <- 1
  phi_2 <- -.9


  Xt1 <- arima.sim(list(order=c(2,0,0), ar=c(phi_1,phi_2)), n = n)*3
  Xt2 <- 5*Xt1 + rnorm(n)
  Xt3 <- 3*Xt1 + rnorm(n)

  Xt4 <- arima.sim(list(order=c(2,0,0), ar=c(1.5,-.75)), n = n)*3
  Xt5 <- 2*Xt4 + rnorm(n)


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
  plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral densities for the first ralization")
  for(a in 2:p){
    lines(Re(F_tp1[a,a,])*n, col= a)
  }
  f_xx1 <- F_tp1*n
  rm(U)
  rm(X_tp)
  rm(F_tp_list)
  gc()


  # LSPCA
  # To observe the effect of smoothing, compare results for nu = 0 with nu = 0.5.
  Ex2 <- LSPCA(n,p, f_xx1, d=2, lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
               control = list(fan_maxinc = 10, verbose = 0), eta=500, s=5, n_iter = 20, nu=0.5)


  ## True eigenvectors
  omega <- seq(0,.5, length.out = n/2)
  f_xx <- array(0,dim=c(p,p,length(omega)))
  f_x_omega_1 <- 9/(2.81 - 3.8*cos(2*pi*omega) +(1.8)*cos(4*pi*omega))
  f_x_omega_2 <- 9/(3.8125 - 5.25*cos(2*pi*omega) +(1.5)*cos(4*pi*omega))
  f_e <- 1
  for(ell in 1:length(omega)){
    f_xx[,,ell] <- diag(f_e, p)
    f_xx[1:5,1:5,ell] <- matrix(c(f_x_omega_1[ell],5*f_x_omega_1[ell], 3*f_x_omega_1[ell],0,0,
                                  5*f_x_omega_1[ell], 25*f_x_omega_1[ell]+f_e, (15)*f_x_omega_1[ell],0,0,
                                  3*f_x_omega_1[ell], (15)*f_x_omega_1[ell], 9*f_x_omega_1[ell]+f_e,0,0,
                                  0,0,0, f_x_omega_2[ell], 2*f_x_omega_2[ell],
                                  0,0,0,2*f_x_omega_2[ell],4*f_x_omega_2[ell]+f_e), nrow=5, byrow = TRUE)


  }


  f_evec11 <- matrix(0, nrow=p, ncol = length(omega))
  for(ell in 1:length(omega)){
    f_evec_ell <- eigen(f_xx[,,ell])$vectors[,1]
    f_evec11[,ell] <- f_evec_ell
    gc()
  }


  f_evec12 <- matrix(0, nrow=p, ncol = length(omega))
  for(ell in 1:length(omega)){
    f_evec_ell <- eigen(f_xx[,,ell])$vectors[,2]
    f_evec12[,ell] <- f_evec_ell
    gc()
  }


  ## Comparing the true v.s. estimated PC loadings
  plot(abs(Ex2[[1]][1,]), main = "1st coordinate - Mod 1st sample evec over freq", pch=18, ylab="", ylim = c(0,1), type="l")
  points(Mod(f_evec11[1,]), pch=18, ylab="", col = "green", type="l")
  legend("topright", inset=.002, legend=c("Estimated", "Population"),
         col=c("black", "green"), lty=c(1,1))


  plot(abs(Ex2[[1]][2,]), main = "2nd coordinate - Mod 1st sample evec over freq", pch=18, ylab="", ylim = c(0,1), type="l")
  points(Mod(f_evec11[2,]), pch=18, ylab="", col = "green", type="l")
  legend("topleft", inset=.002, legend=c("Estimated", "Population"),
         col=c("black", "green"), lty=c(1,1))



  plot(abs(Ex2[[2]][1,]), main = "1st coordinate - Mod 2nd sample evec over freq", pch=18, ylab="", ylim = c(0,1), type="l")
  points(Mod(f_evec12[1,]), pch=18, ylab="", col = "green", type="l")
  legend("topright", inset=.002, legend=c("Estimated", "Population"),
         col=c("black", "green"), lty=c(1,1))


  plot(abs(Ex2[[2]][2,]), main = "2nd coordinate - Mod 2nd sample evec over freq", pch=18, ylab="", ylim = c(0,1), type="l")
  points(Mod(f_evec12[2,]), pch=18, ylab="", col = "green", type="l")
  legend("topright", inset=.002, legend=c("Estimated", "Population"),
         col=c("black", "green"), lty=c(1,1))


  ##########################################################
  ## Extra Example 3: with 1-dimensional principal subspace
  ##########################################################

  set.seed(100*3+15)

  # Generating data
  p <- 64
  n <- 1024
  phi_1 <- 1.5
  phi_2 <- -.75
  a1 <- 1/20
  a2 <- -1/1.15
  omega <- seq(0,.5, length.out = n/2)
  f_xx <- array(0,dim=c(p,p,length(omega)))
  c1 <- 3
  f_x_omega0 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
  phi_1 <- 1.5 +.05
  f20 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
  phi_1 <- 1.5 -.05
  f30 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
  phi_1 <- 1.5 +.15
  f40 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
  phi_1 <- 1.5 -.15
  f50 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec

  f_x_omega <- f_x_omega0
  for(ell in 1:length(omega)) {
    o_ell <- omega[ell]
    if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f_x_omega[ell] <- 0
  }
  f2 <- f20
  for(ell in 1:length(omega)) {
    o_ell <- omega[ell]
    if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f2[ell] <- 0
  }
  f3 <- f30
  for(ell in 1:length(omega)) {
    o_ell <- omega[ell]
    if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f3[ell] <- 0
  }
  f4 <- f40
  for(ell in 1:length(omega)) {
    o_ell <- omega[ell]
    if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f4[ell] <- 0
  }
  f5 <- f50
  for(ell in 1:length(omega)) {
    o_ell <- omega[ell]
    if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f5[ell] <- 0
  }

  f_e <- 1/4
  for(ell in 1:length(omega)){
    f_xx[,,ell] <- diag(f_e, p)
    f_xx[1:5,1:5,ell] <- matrix(c(f_x_omega[ell],(2*1.1)*f_x_omega[ell],1.2*f_x_omega[ell],1.25*f_x_omega[ell], 3*.75*f_x_omega[ell],
                                  2*1.1*f_x_omega[ell], 4*1.1*1.1*f_x_omega[ell]+f_e+1*f2[ell], 2*1.1*1.2*f_x_omega[ell],2*1.1*1.25*f_x_omega[ell],2*1.1*3*.75*f_x_omega[ell],
                                  1.2*f_x_omega[ell], 2*1.1*1.2*f_x_omega[ell], 1.2*1.2*f_x_omega[ell]+f_e+f3[ell], 1.2*1.25*f_x_omega[ell],1.2*3*.75*f_x_omega[ell],
                                  1.25*f_x_omega[ell],2*1.1*1.25*f_x_omega[ell],1.2*1.25*f_x_omega[ell],1.25*1.25*f_x_omega[ell]+f_e+1*f4[ell],1.25*3*.75*f_x_omega[ell],
                                  3*.75*f_x_omega[ell],2*1.1*3*.75*f_x_omega[ell],1.2*3*.75*f_x_omega[ell],1.25*3*.75*f_x_omega[ell],9*.75*.75*f_x_omega[ell]+f_e+9*f5[ell]), nrow=5, byrow = TRUE)
  }




  len_freq <- n/2
  Coord10 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord20 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord30 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord40 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord50 <- matrix(0,nrow = 5, ncol = len_freq)

  Coord1 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord2 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord3 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord4 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord5 <- matrix(0,nrow = 5, ncol = len_freq)


  phi_1 <- 1.5
  phi_2 <- -.75
  a1 <- 1/20
  a2 <- -1/1.15


  Xt10 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
  Xt1 <- pass.filt(Xt10, W=c(0.05, 0.25), type="pass", method = "Butterworth")
  temp1 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

  phi_1 <- 1.5 +.05
  yt20 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
  yt2 <- pass.filt(yt20, W=c(0.05, 0.25), type="pass", method = "Butterworth")
  temp2 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

  phi_1 <- 1.5 -.05
  yt30 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
  yt3 <- pass.filt(yt30, W=c(0.05, 0.25), type="pass", method = "Butterworth")
  temp3 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

  phi_1 <- 1.5 +.15
  yt40 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
  yt4 <- pass.filt(yt40, W=c(0.05, 0.25), type="pass", method = "Butterworth")
  temp4 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

  phi_1 <- 1.5 -.15
  yt50 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
  yt5 <- pass.filt(yt50, W=c(0.05, 0.25), type="pass", method = "Butterworth")
  temp5 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

  Xt2 <- 2*1.1*Xt1 + 1*yt2 + rnorm(n,0,1/2)
  Xt3 <- 1.2*Xt1 + yt3 + rnorm(n,0,1/2)
  Xt4 <- 1.25*Xt1 + 1*yt4 + rnorm(n,0,1/2)
  Xt5 <- 3*0.75*Xt1 + 3*yt5 + rnorm(n,0,1/2)

  X_omega_2 <- cbind(Xt1,Xt2,Xt3,Xt4,Xt5)
  X_wn <- rmvnorm(n, sigma = diag(p-(ncol(X_omega_2))))
  X_original <- cbind(X_omega_2, X_wn)
  shuffled_index <- 1:p
  X1 <- X_original[,shuffled_index]


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
  f_xx1 <- F_tp1*n
  rm(U)
  rm(X_tp)
  rm(F_tp_list)
  gc()


  len_freq <- n/2
  Coord10 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord20 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord30 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord40 <- matrix(0,nrow = 5, ncol = len_freq)
  Coord50 <- matrix(0,nrow = 5, ncol = len_freq)


  nu_v <- c(0,.2,.4,.6,.8,1)
  for(k in 1:6){
    nu <- nu_v[k]

    LSDPCA_ADMM_SOAP_Ex1 <- LSDPC_ADMM_SOAP_PD_fun(n,p, f_xx1, lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
                                                   control = list(fan_maxinc = 10, verbose = 0), eta=200, s=5, n_iter = 20, nu=nu_v[k])
    Coord10[k,] <- abs(LSDPCA_ADMM_SOAP_Ex1[[1]][1,])
    Coord20[k,] <- abs(LSDPCA_ADMM_SOAP_Ex1[[1]][2,])
    Coord30[k,] <- abs(LSDPCA_ADMM_SOAP_Ex1[[1]][3,])
    Coord40[k,] <- abs(LSDPCA_ADMM_SOAP_Ex1[[1]][4,])
    Coord50[k,] <- abs(LSDPCA_ADMM_SOAP_Ex1[[1]][5,])

    gc()
  }



  f_evec11 <- matrix(0, nrow=p, ncol = length(omega))
  for(ell in 1:length(omega)){
    f_evec_ell <- eigen(f_xx[,,ell])$vectors[,1]
    f_evec11[,ell] <- f_evec_ell
    gc()
  }

  f_MT_evec11 <- matrix(0, nrow=p, ncol = length(omega))
  for(ell in 1:length(omega)){
    f_MT_evec_ell <- eigen(f_xx1[,,ell])$vectors[,1]
    f_MT_evec11[,ell] <- f_MT_evec_ell
    gc()
  }

  #################################################################################################################
  ## Plots
  #################################################################################################################
  library(latex2exp)

  ## Plot of the first coordinate of the 1st eigenvector across frequency
  col_map <- c("blue", "brown", "green", "red", "pink", "purple")
  plot(Coord10[1,], main = "", pch=18, ylab="", ylim = c(0,1), type="l", col="blue")
  points(Coord10[4,], pch=18, ylab="", ylim = c(0,1), type="l", col="red")
  points(Mod(f_evec11[1,]), pch=18, ylab="", col = "black", type="l")
  points(Mod(f_MT_evec11[1,]), pch=18, ylab="", col = "black", type="l", lty=2)
  legend("topright", legend = c(TeX("$\\theta=0$"), TeX("$\\theta=0.6$"), "Population"), lty = c(1,1,1), col=c("blue", "red", "black"))



  ## Plot of the second coordinate of the 1st eigenvector across frequency
  col_map <- c("blue", "brown", "green", "red", "pink", "purple")
  plot(Coord20[1,], main = "", pch=18, ylab="", ylim = c(0,1), type="l", col="blue")
  points(Coord20[4,], pch=18, ylab="", ylim = c(0,1), type="l", col="red")
  points(Mod(f_evec11[2,]), pch=18, ylab="", col = "black", type="l")
  points(Mod(f_MT_evec11[2,]), pch=18, ylab="", col = "black", type="l", lty=2)
  legend("topright", legend = c(TeX("$\\theta=0$"), TeX("$\\theta=0.6$"), "Population"), lty = c(1,1,1), col=c("blue", "red", "black"))



  ## Plot of the third coordinate of the 1st eigenvector across frequency
  col_map <- c("blue", "brown", "green", "red", "pink", "purple")
  plot(Coord30[1,], main = "", pch=18, ylab="", ylim = c(0,1), type="l", col="blue")
  points(Coord30[4,], pch=18, ylab="", ylim = c(0,1), type="l", col="red")
  points(Mod(f_evec11[3,]), pch=18, ylab="", col = "black", type="l")
  points(Mod(f_MT_evec11[3,]), pch=18, ylab="", col = "black", type="l", lty=2)
  legend("topright", legend = c(TeX("$\\theta=0$"), TeX("$\\theta=0.6$"), "Population"), lty = c(1,1,1), col=c("blue", "red", "black"))


  ## Plot of the fourth coordinate of the 1st eigenvector across frequency
  col_map <- c("blue", "brown", "green", "red", "pink", "purple")
  plot(Coord40[1,], main = "", pch=18, ylab="", ylim = c(0,1), type="l", col="blue")
  points(Coord40[4,], pch=18, ylab="", ylim = c(0,1), type="l", col="red")
  points(Mod(f_evec11[4,]), pch=18, ylab="", col = "black", type="l")
  points(Mod(f_MT_evec11[4,]), pch=18, ylab="", col = "black", type="l", lty=2)
  legend("topright", legend = c(TeX("$\\theta=0$"), TeX("$\\theta=0.6$"), "Population"), lty = c(1,1,1), col=c("blue", "red", "black"))


  ## Plot of the fifth coordinate of the 1st eigenvector across frequency
  col_map <- c("blue", "brown", "green", "red", "pink", "purple")
  plot(Coord50[1,], main = "", pch=18, ylab="", ylim = c(0,1), type="l", col="blue")
  points(Coord50[4,], pch=18, ylab="", ylim = c(0,1), type="l", col="red")
  points(Mod(f_evec11[5,]), pch=18, ylab="", col = "black", type="l")
  points(Mod(f_MT_evec11[5,]), pch=18, ylab="", col = "black", type="l", lty=2)
  legend("topright", legend = c(TeX("$\\theta=0$"), TeX("$\\theta=0.6$"), "Population"), lty = c(1,1,1), col=c("blue", "red", "black"))



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
