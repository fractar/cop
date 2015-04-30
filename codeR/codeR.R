
# packages
library("gumbel")
library("Kendall")
library("gplots")
library("CDVine")
library("copula")
library("fCopulae")
library("copBasic")

# chargement des données
data(windEchirolles)
data(windStMartin)
n <- min(NROW(windStMartin), NROW(windEchirolles))
id2keep <- !is.na(windStMartin$WIND.HIGH[1:n]) & !is.na(windEchirolles$WIND.HIGH[1:n])
x <- windStMartin$WIND.HIGH[1:n][id2keep]/3.6
y <- windEchirolles$WIND.HIGH[1:n][id2keep]/3.6
n <- length(x)

plot(x,y)


# Analyse des corrélations

# coefficient de pearson (corrélation linéaire)
cor(x,y)
# rho de spearman (corrélation non-linéaire)
cor(x,y,method=c("spearman"))
# taux de kendall 
Kendall(x,y) 



# rank / n
Rx = rank(x,ties.method = c("random"))
Ry = rank(y,ties.method = c("random"))

# rank-rank plot
plot(Rx/n,Ry/n)

# histogramme 3D

h2d=hist2d(Rx/n,Ry/n,show=FALSE, same.scale=TRUE, nbins=c(10,10))
persp(h2d$x,h2d$y,h2d$counts,col="lightblue",theta=50)
#hist3D(h2d$x,h2d$y,h2d$counts)

my_hist3d <- function(x, y, freq=FALSE, nclass="auto") {
  n<-length(x)
  if (nclass == "auto") { nclass<-ceiling(sqrt(nclass.Sturges(x))) }
  breaks.x <- seq(min(x),max(x),length=(nclass+1))
  breaks.y <- seq(min(y),max(y),length=(nclass+1))
  h <- NULL
  for (i in 1:nclass) 
    for (j in 1:nclass) 
      h <- c(h, sum(x <= breaks.x[j+1] & x >= breaks.x[j] & y <= breaks.y[i+1] & y >= breaks.y[i] ) )
  if (freq) h <- h / n
  xx <- as.factor(round(mean(breaks.x[1:2])+(0:(nclass-1))*diff(breaks.x[1:2]), 1))
  yy <- as.factor(round(mean(breaks.y[1:2])+(0:(nclass-1))*diff(breaks.y[1:2]), 1))
  res <- cbind(expand.grid(xx,yy), h)
  colnames(res) <- c(deparse(substitute(x)),deparse(substitute(y)),'Frequency')
  formu <- as.formula(paste("Frequency ~ ", paste(colnames(res)[1:2], collapse= "+")))
  rotate.cloud(formu, res, panel.3d.cloud=panel.3dbars, col.facet='lightblue', 
       xbase=1, ybase=1, scales=list(arrows=FALSE, col=1), 
        par.settings = list(axis.line = list(col = "transparent")),theta=40)
}

StMartin <- rx
Echirolles <- ry
my_hist3d(StMartin, Echirolles, nclass=10)


# Recherche de la copule 
BiCopSelect(Rx/n,Ry/n,familyset=c(1:10), selectioncrit ="AIC",indeptest=FALSE,level=0.05)

# copule empirique - Deheuvels (1979)
Cop_Deheuvels<-function(u,v) {
  return(sum( (Rx/n <= u) & (Ry/n <= v) ) / n)
}
# on "vectorise" la fonction Cop_Deheuvels
Cop_Deheuvels_vect <- Vectorize(Cop_Deheuvels, SIMPLIFY = TRUE)

u<-(1:100)/100
v<-u
fdr<-outer(u,v,Cop_Deheuvels_vect)
contour(u,v,fdr,col="blue")
persp(u, v, fdr, theta = 30, phi = 15, expand = 0.7, col = "blue", main = "fdr de la copule empirique (Deheuvels)") 


# copule empirique  bivariée Nelsen, 2006
uv <- data.frame(Rx/n,Ry/n)
u<-(1:100)/100
v<-u

fdrr<-function(u,v){
  return(EMPIRcop(u,v,para=uv))
}

fdr<-outer(u,v,fdrr)
contour(u,v,fdr,col="blue")
persp(u, v, fdr, theta = 30, phi = 15, expand = 0.7, col = "blue", main = "fdr de la copule empirique (Deheuvels)") 
























