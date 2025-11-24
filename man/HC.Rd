\name{HC}
\alias{HC}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
EEG recording of a Healthy Control subject
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The dataframe contains 64 chanell EEG recording of a Healthy Control subject.
}
\usage{data("HC")}
\format{
  The format is:
 num [1:2048, 1:64] -3.85 -1.32 1.43 8.39 11.23 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:64] "V1" "V2" "V3" "V4" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

library(astsa)
library(waveslim)
library(ggplot2)


## Data

X <- HC
p <- ncol(X)
n <- nrow(X)


## Multitaper estimate of the spectral density matrix
U <- sine.taper(n,10)
X_tp <- apply(U, MARGIN = 2, function(u) u*X, simplify = FALSE)
F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

len_freq <- n/2
F_tp1 <- array(0, c(p, p, len_freq))
for (ell in 1:len_freq) {
  for(j in 1:length(F_tp_list)){
    F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
  }
  F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
}
plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral densities for the first ralization", ylim=c(1,2000))
for(a in 2:p){
  lines(Re(F_tp1[a,a,])*n, col= a)
}
f_xx1 <- F_tp1*n
rm(U)
rm(X_tp)
rm(F_tp_list)
gc()

## Localized and sparse PCA
LSDPCA_ADMM_SOAP_Ex2 <- LSPCA(n,p, f_xx1, lambda = 0.5 * sqrt(log(p) / n), d=2, lr = 0.02, maxiter = 60,
                                               control = list(fan_maxinc = 10, verbose = 0), eta=52, s=8, n_iter = 20, nu=0.6)



## Plot

xlab <- "Hz"
ylab <- "Coordinate"
legend_title = ""
asp = .5 #0.2
bar_height = 10/2
font_size = 20/2



Localized_Est <- selector(LSDPCA_ADMM_SOAP_Ex2[[1]],f_xx1,n/2,52,p)
evecs <- t(Re(Localized_Est[[1]]))
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
                                                                                                                            colors = colors) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = "Real part of the first PC loadings - Healthy Control") +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), labels=function(x) x*64/n)



}
\keyword{datasets}
