
# packages
library("gumbel")
library("Kendall")
library("gplots")
library("CDVine")
library("copula")
library("fCopulae")
library("copBasic")
library("fractal")
library("akima")
library("pspearman")
library("kSamples")
library("TeachingDemos")
library("tcltk2")


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

# rank-rank plot (dépendogramme)
plot(Rx/n,Ry/n)

# histogramme 3D

h2d=gplots::hist2d(Rx/n,Ry/n,show=FALSE, same.scale=TRUE, nbins=c(10,10))
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

StMartin <- Rx/n
Echirolles <- Ry/n
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


# densité de la copule empirique 

#Noyau gaussien
nn<-100
zz<-density2d(x=Rx/n,y=Ry/n,n=nn,limits=c(0,1,0,1))
persp(zz$x,zz$y,z=zz$z, theta = 30, phi = 15, expand = 0.7, 
	col = "blue", main = "densité noyau gaussien de la copule",
	xlab="StMartin", ylab="Echirolles", zlab="") 

#Noyau gaussien
zz<-kde2d(x=Rx/n,y=Ry/n)
persp(zz$x,zz$y,z=zz$z, theta = 30, phi = 15, expand = 0.7, 
	col = "blue", main = "densité noyau gaussien de la copule",
	xlab="StMartin", ylab="Echirolles", zlab="") 

# noyau gaussien
d=2
s = sqrt(d/12)
K_gaussien<-function(x,y) {
  return((2*pi)^(-d/2) / (s^d) * exp(- (x^2+y^2)/(2*s^2) ))
}

# densité de la copule
h=0.2
density_estim<-function(u,v) {
 return(1/(n*h^2) * sum(K_gaussien( (u-Rx/n)/h , (v-Ry/n)/h ) ) )
}

# on "vectorise" la fonction density_estim
density_estim_vect <- Vectorize(density_estim, SIMPLIFY = TRUE)

u<-(1:100)/100
v<-u
fdr<-outer(u,v,density_estim_vect)
contour(u,v,fdr,col="blue")
persp(u, v, fdr, theta = 30, phi = 15, expand = 0.7, col = "blue", main = "densité de la copule") 

#Noyau Epanechnikov
mat<-cbind(Rx/n,Ry/n)
res<-KDE(mat)
plot(res,style="perspective")

#Khi-Plot

Khiplot<-function(x,y,n) {

Hi<-rep(NA,n)
Fi<-rep(NA,n)
Gi<-rep(NA,n)
Xi<-rep(NA,n)
Li<-rep(NA,n)

for (i in 1:n){
Hi[i]<-sum(x[i]>=x[-i] & y[i]>=y[-i])/(n-1)
Fi[i]<-sum(x[i]>=x[-i])/(n-1)
Gi[i]<-sum(y[i]>=y[-i])/(n-1)	
Xi[i]<-(Hi[i]-Fi[i]*Gi[i])/(sqrt(Fi[i]*(1-Fi[i])*Gi[i]*(1-Gi[i])))
Li[i]<-4*sign((Fi[i]-0.5)*(Gi[i]-0.5))*max((Fi[i]-0.5)^2,(Gi[i]-0.5)^2)
}
ind=which(abs(Li) < 4*(1/(n-1)-0.5)^2)	
plot(Li[ind],Xi[ind],main="Khi plot empirique",col="blue")
# abline(h=4*(1/(length(rendFO)-1)-0.5)^2)
# abline(h=-4*(1/(length(rendFO)-1)-0.5)^2)
#dev.off()
}
Khiplot(x,y,n)

#K-Plot

BiCopKPlot(Rx/n,Ry/n)


#test classique d'indépendance
statT = cor(x,y)*sqrt(n-2)/sqrt(1-cor(x,y)^2)
abs(statT) > qt(p=1-0.05/2,df=n-2)

statN = sqrt(n-3)*0.5*log( (1+cor(x,y))/(1-cor(x,y)) ) 
abs(statN) > qnorm(p=1-0.05)

# Test de Spearman
cor.test(x,y,method="spearman")

#Test de Kendall
cor.test(x,y,method="kendall")

#Test de Van Der Waerden
qn.test(x,y,test="vdW")


fitdistr(x,"normal")
ks.test(x,"pnorm",mean=-0.01599282,sd=0.95044161)


####### Copule de Clayton #######


#Estimation du tau de Kendall par la méthode des moments
tauKendall<-Kendall(x,y)

# Estimation semi parametrique (CML) 
paramClayton<-2*tauKendall$tau/(1-tauKendall$tau)
myClaytonCopula<-archmCopula(family="clayton",param=paramClayton)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myClaytonCopula)
thetaN<-coef(CML)

fitClaytonCopula<-claytonCopula(thetaN,dim=2)

#Khi-plot de la copule de Clayton estimé

number = 1000
randomEstime<-rCopula(number,fitClaytonCopula)
XClayton<-randomEstime[,1]
YClayton<-randomEstime[,2]
Khiplot(XClayton,YClayton,number)

#K-plot de la copule de Clayton estimé
KplotXRank<-rank(XClayton,ties.method="random")
KplotYRank<-rank(YClayton,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number))

#Bootstrap paramétrique

Ui = Rx/(n+1)
Vi = Ry/(n+1)

# copule empirique  bivariée Nelsen, 2006
uv <- data.frame(Ui,Vi)
u<-(1:100)/100
v<-u

fdrEmpirique<-function(u,v){
  return(EMPIRcop(u,v,para=uv))
}

UV<-cbind(Ui,Vi)
N=100
randomX = NULL
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitClaytonCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,fitClaytonCopula)
    rankX<-rank(randomXY[,1],ties.method="random")
    rankY<-rank(randomXY[,2],ties.method="random")
    Ui = rankX/(n+1)
    Vi = rankY/(n+1)
    uv <- data.frame(Ui,Vi)    
    fdrEmpirique<-function(u,v){
	return(EMPIRcop(u,v,para=uv))
    }
    tauKendall=Kendall(randomXY[,1],randomXY[,2])

    # Estimation semi parametrique (CML) 
    paramClayton<-2*tauKendall$tau/(1-tauKendall$tau)
    myClaytonCopula<-archmCopula(family="clayton",param=paramClayton)
    CML<-fitCopula(data=cbind(Ui,Vi),copula=myClaytonCopula)
    thetaN<-coef(CML)
    fitClaytonCopula<-claytonCopula(thetaN,dim=2)
    
    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitClaytonCopula))^2))
}

alpha = 0.05
L = sort(Dnk)[floor((1-alpha)*N)]
#Règle de décision
Dn > L
#p-value
sum(Dnk > Dn)/length(Dnk)


###### Copule gaussienne ######


#Bootstrap paramétrique
myNormalCopula<-ellipCopula(family="normal",param=cor(x,y))
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myNormalCopula)
thetaN<-coef(CML)

fitNormalCopula<-normalCopula(thetaN,dim=2)


Ui = Rx/(n+1)
Vi = Ry/(n+1)

# copule empirique  bivariée Nelsen, 2006
uv <- data.frame(Ui,Vi)
u<-(1:100)/100
v<-u

fdrEmpirique<-function(u,v){
  return(EMPIRcop(u,v,para=uv))
}

UV<-cbind(Ui,Vi)
N=100
randomX = NULL
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitNormalCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,fitNormalCopula)
    rankX<-rank(randomXY[,1],ties.method="random")
    rankY<-rank(randomXY[,2],ties.method="random")
    Ui = rankX/(n+1)
    Vi = rankY/(n+1)
    uv <- data.frame(Ui,Vi)    
    fdrEmpirique<-function(u,v){
	return(EMPIRcop(u,v,para=uv))
    }

    # Estimation semi parametrique (CML) 
    myNormalCopula<-ellipCopula(family="normal",param=cor(randomXY[,1],randomXY[,2]))
    CML<-fitCopula(data=cbind(Ui,Vi),copula=myNormalCopula)
    thetaN<-coef(CML)
    fitNormalCopula<-normalCopula(thetaN,dim=2)
    
    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitNormalCopula))^2))
}

alpha = 0.05
L = sort(Dnk)[floor((1-alpha)*N)]
#Règle de décision
Dn > L
#p-value
sum(Dnk > Dn)/length(Dnk)