
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
library("latticeExtra")
library(fitdistrplus)
library("tcltk2")

# chargement des données
data(windEchirolles)
data(windStMartin)
n <- min(NROW(windStMartin), NROW(windEchirolles))
id2keep <- !is.na(windStMartin$WIND.HIGH[1:n]) & !is.na(windEchirolles$WIND.HIGH[1:n])
x <- windStMartin$WIND.HIGH[1:n][id2keep]/3.6
y <- windEchirolles$WIND.HIGH[1:n][id2keep]/3.6
n <- length(x)

# rangs 
Rx = rank(x,ties.method = c("random"))
Ry = rank(y,ties.method = c("random"))

#postscript("scatter_et_rankrank.ps")
par(mfrow=c(1,2))
# scatter plot (diagramme de dispersion)
plot(x,y,xlab="Saint-Martin-en-Haut",ylab="Echirolles",main="Diagramme de dispersion",col="blue")
# rank-rank plot (dépendogramme)
plot(Rx/n,Ry/n,xlab="Saint-Martin-en-Haut",ylab="Echirolles",main="Graphique des rangs",col="blue")
#dev.off()


# Mesures de dépendance

# coefficient de pearson (corrélation linéaire)
cor(x,y)
# rho de spearman (corrélation non-linéaire)
cor(x,y,method=c("spearman"))
# taux de kendall 
Kendall(x,y) 

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
#postscript("copule_empirique.ps")
persp(u, v, fdr, theta = 30, phi = 15, expand = 0.7, col = "blue", main = "fdr de la copule empirique (Deheuvels)") 
#dev.off()

# copule empirique  bivariée Nelsen (2006)
uv <- data.frame(Rx/n,Ry/n)
u<-(1:100)/100
v<-u

fdrr<-function(u,v){
  return(EMPIRcop(u,v,para=uv))
}

fdr<-outer(u,v,fdrr)
contour(u,v,fdr,col="blue")
persp(u, v, fdr, theta = 30, phi = 15, expand = 0.7, col = "blue", main = "fdr de la copule empirique (Nelsen)") 


# étude de la dépendance de queue
u<-(1:(n-1))/n
lambda_L = fdrr(u,u) /u
lambda_U = (1-2*u+fdrr(u,u)) / (1-u)
#postscript("dependance_queue_empir.ps")
par(mfrow=c(1,2))
plot(u,lambda_L,type="l",col="blue",xlab="u",ylab="lambda_L",main="Estimation des dépendances de queue à gauche")
abline(h=0.22)
plot(u,lambda_U,type="l",col="blue",xlab="u",ylab="lambda_U",main="Estimation des dépendances de queue à droite")
abline(h=0.31)
#dev.off()

# densité de la copule empirique 

#Noyau gaussien
nn<-100
zz<-density2d(x=Rx/n,y=Ry/n,n=nn,limits=c(0,1,0,1))
#postscript("densite_empir_gaussien.ps")
persp(zz$x,zz$y,z=zz$z, theta = 30, phi = 15, expand = 0.7, 
	col = "blue", main = "densité empirique de la copule (noyau gaussien)",
	xlab="St-Martin", ylab="Echirolles", zlab="c(u,v)") 
#dev.off()

#Noyau gaussien
zz<-kde2d(x=Rx/n,y=Ry/n)
persp(zz$x,zz$y,z=zz$z, theta = 30, phi = 15, expand = 0.7, 
	col = "blue", main = "densité empirique de la copule (noyau gaussien)",
	xlab="St-Martin", ylab="Echirolles", zlab="c(u,v)") 

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
persp(u, v, fdr, theta = 30, phi = 15, expand = 0.7, col = "blue", main = "densité empirique de la copule (noyau gaussien)",xlab="St-Martin", ylab="Echirolles", zlab="c(u,v)") 

#Noyau Epanechnikov
mat<-cbind(Rx/n,Ry/n)
res<-KDE(mat)
#postscript("densite_empir_Epanechnikov.ps")
par(mfrow=c(1,1))
plot(res,style="perspective",col="blue", theta = 30, phi = 15, main = "densité empirique de la copule (noyau d'Epanechnikov)",xlab="St-Martin", ylab="Echirolles", zlab="c(u,v)")
#dev.off()



#Khi-Plot

Khiplot<-function(x,y,n,xtitle,ytitle,title,ylim) {
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
plot(Li[ind],Xi[ind],main=title,col="blue",xlab=xtitle,ylab=ytitle,ylim=ylim)
abline(h=1.78/sqrt(n))
abline(h=-1.78/sqrt(n))
}
#postscript("chi_plot_empir.ps")
Khiplot(x,y,n,"Li","Xi","Khi-plot empirique",ylim=c(-0.5,1))
#dev.off()

#K-Plot
#postscript("K_plot_empir.ps")
BiCopKPlot(Rx/n,Ry/n,col="blue",main="K-plot empirique")
#dev.off()

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
  colnames(res) <- c(deparse(substitute(x)),deparse(substitute(y)),'Frequence')
  formu <- as.formula(paste("Frequence ~ ", paste(colnames(res)[1:2], collapse= "+")))
  cloud(formu, res, panel.3d.cloud=panel.3dbars, col.facet='lightblue', 
       xbase=1, ybase=1, scales=list(arrows=FALSE, col=1), 
        par.settings = list(axis.line = list(col = "transparent")),theta=40,screen=list(x=-60,y=0,z=0))
}

StMartin <- Rx/n
Echirolles <- Ry/n
my_hist3d(StMartin, Echirolles, nclass=10)

################################### Copule de Clayton #########################

##############    1) METHODE DES MOMENTS  #######################

#Estimation par la méthode des moments (inversion du tau de Kendall)

tauKendall<-Kendall(x,y)
paramClayton<-2*tauKendall$tau/(1-tauKendall$tau)
myClaytonCopula<-claytonCopula(param=paramClayton)


##############    2) METHODE CML  #######################

#Estimation du tau de Kendall par la méthode des moments
tauKendall<-Kendall(x,y)
paramClayton<-2*tauKendall$tau/(1-tauKendall$tau)

# Estimation semi parametrique (CML) 

myClaytonCopula<-claytonCopula(param=paramClayton)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myClaytonCopula,start=paramClayton)
thetaN<-coef(CML)

fitClaytonCopula<-claytonCopula(thetaN,dim=2)

##############    3) METHODE MV  #######################

dgumbel  <-  function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel  <-  function(q,a,b) exp(-exp((a-q)/b))
qgumbel  <-  function(p,a,b) a-b*log(-log(p))

#Choix des lois marginales
xpar_gamma <- mledist(x, "gamma")
ypar_gamma <- mledist(y, "gamma")
xpar_exp <- mledist(x, "exp")
ypar_exp <- mledist(y, "exp")
xpar_norm <- mledist(x, "norm")
ypar_norm <- mledist(y, "norm")
xpar_weibull <- mledist(x, "weibull")
ypar_weibull <- mledist(y, "weibull")
xpar_gumbel <- mledist(x, "gumbel",start=list(a=10,b=5))
ypar_gumbel <- mledist(y, "gumbel",start=list(a=10,b=5))

#Function to determine the parameters of the distribution
distr.choose <- function(x,xpar_distr,pdistr,leg="Distribution",wind)
{
  x <-sort(x)
  x.emp <- ecdf(x)
  plot(x.emp, col="grey", main=wind)
  nb_param = length(xpar_distr$estimate)
  if (nb_param == 1)
  {
    param1 <- xpar_distr$estimate[1] 
    lines(x,pdistr(x, param1), col="blue4")
  } else if (nb_param == 2)
  {
    param1 <- xpar_distr$estimate[1] 
    param2 <- xpar_distr$estimate[2]
    lines(x,pdistr(x, param1, param2), col="blue4")
  }
  legend("topleft",legend=c("Empirical distr",leg), col=c("grey","blue4"),lty=c(1,1)) 
}

#Gamma distribution
distr.choose(x,xpar_gamma,pgamma,"Gamma distr","St Martin Wind Distribution")
distr.choose(y,ypar_gamma,pgamma,"Gamma distr","Echirolles Wind Distribution")
#Normal distribution
distr.choose(x,xpar_norm,pnorm,"Normal distr","St Martin Wind Distribution")
distr.choose(y,ypar_norm,pnorm,"Normal distr","Echirolles Wind Distribution")
#Exponential distribution
distr.choose(x,xpar_exp,pexp,"Exp distr","St Martin Wind Distribution")
distr.choose(y,ypar_exp,pexp,"Exp distr","Echirolles Wind Distribution")
#Weibull distribution
distr.choose(x,xpar_weibull,pweibull,"Weibull distr","St Martin Wind Distribution")
distr.choose(y,ypar_weibull,pweibull,"Weibull distr","Echirolles Wind Distribution")
#Gumbel distribution
distr.choose(x,xpar_gumbel,pgumbel,"Gumbel distr","St Martin Wind Distribution")
distr.choose(y,ypar_gumbel,pgumbel,"Gumbel distr","Echirolles Wind Distribution")

nx <- length(unique(x))
ny <- length(unique(y))
qqplot(sort(unique(x)),qnorm(seq(1:nn)/nn,xpar_norm$estimate[1],xpar_norm$estimate[2]), type="l",main="QQ-PLOT")
qqplot(sort(unique(y)),qnorm(seq(1:nn)/nn,ypar_norm$estimate[1],ypar_norm$estimate[2]), type="l")

myClaytonMvd<-mvdc(copula=archmCopula(family="clayton",param=0.5),margins=c("gumbel","gamma"),paramMargins=list(list(a=xpar_gumbel$estimate[1],b=xpar_gumbel$estimate[2]),list(shape=ypar_gamma$estimate[1],scale=ypar_gamma$estimate[2])))
tauKendall<-Kendall(x,y)
paramClayton<-2*tauKendall$tau/(1-tauKendall$tau)
start<-c(xpar_gumbel$estimate[1],xpar_gumbel$estimate[2],ypar_gamma$estimate[1],ypar_gamma$estimate[2],paramClayton)
fitTest<-fitMvdc(cbind(x,y),myClaytonMvd,start=start,optim.control=list(trace=TRUE,maxit=2000))

########################## TEST GRAPHIQUE ET STATISTIQUE #######################

# Empirique
Khiplot(x,y,n,"Li","Xi","Khi-plot empirique",ylim=c(-0.2,0.8))

fitClaytonCopula<-claytonCopula(0.779,dim=2)
testClaytonMvd<-mvdc(copula=fitClaytonCopula,margins=c("gumbel","gamma"),paramMargins=list(list(a=5.024,b=1.840),list(shape=5.003,scale=1.319)))
randomEstime<-rMvdc(n,testClaytonMvd)

#Khi-plot de la copule de Clayton estimé

XClayton<-randomEstime[,1]
YClayton<-randomEstime[,2]
Khiplot(XClayton,YClayton,n,"Li","Xi","Khi-plot de la copule de Clayton estimée",c(-0.2,0.8))

#K-plot de la copule de Clayton estimé
KplotXRank<-rank(XClayton,ties.method="random")
KplotYRank<-rank(YClayton,ties.method="random")
BiCopKPlot(KplotXRank/n,KplotYRank/n,main="K-plot de la copule de Clayton estimée")

#Dépendogramme

plot(KplotXRank/n,KplotYRank/n,xlab="Saint-Martin-en-Haut",ylab="Echirolles",main="Graphique des rangs pour la copule de Clayton estimée")

#### METHODE 1: FONCTION D'UN PACKAGE DE R ####

gofCopula(fitClaytonCopula,cbind(Rx/(n+1),Ry/(n+1)))

#### METHODE 2: CODE MAISON ####

#Bootstrap paramétrique par CML

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
    CML<-fitCopula(data=cbind(Ui,Vi),copula=myClaytonCopula,start=paramClayton)
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

################################### Copule de Gumbel #########################

##############    1) METHODE DES MOMENTS  #######################

#Estimation par la méthode des moments (inversion du tau de Kendall)
tauKendall<-Kendall(x,y)
paramGumbel<-1/(1-tauKendall$tau) 
myGumbelCopula<-gumbelCopula(param=paramGumbel)

##############    2) METHODE CML  #######################

#Estimation du tau de Kendall par la méthode des moments
tauKendall<-Kendall(x,y)
paramGumbel<-1/(1-tauKendall$tau)

myGumbelCopula<-gumbelCopula(param=paramGumbel)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myGumbelCopula)
thetaN<-coef(CML)

fitGumbelCopula<-gumbelCopula(thetaN,dim=2)

##############    3) METHODE MV  #######################

myGumbelMvd<-mvdc(copula=archmCopula(family="gumbel",param=1.5),margins=c("gumbel","gamma"),paramMargins=list(list(a=xpar_gumbel$estimate[1],b=xpar_gumbel$estimate[2]),list(shape=ypar_gamma$estimate[1],scale=ypar_gamma$estimate[2])))
tauKendall<-Kendall(x,y)
paramGumbel<-1/(1-tauKendall$tau)
start<-c(xpar_gumbel$estimate[1],xpar_gumbel$estimate[2],ypar_gamma$estimate[1],ypar_gamma$estimate[2],paramGumbel)
fitTest<-fitMvdc(cbind(x,y),myGumbelMvd,start=start,optim.control=list(trace=TRUE,maxit=2000))

########################## TEST GRAPHIQUE ET STATISTIQUE #######################

fitGumbelCopula<-gumbelCopula(1.491,dim=2)
testGumbelMvd<-mvdc(copula=fitGumbelCopula,margins=c("gumbel","gamma"),paramMargins=list(list(a=5.021,b=1.818),list(shape=5.007,scale=1.319)))
randomEstime<-rMvdc(n,testGumbelMvd)

#Khi-plot de la copule de Gumbel estimé

XGumbel<-randomEstime[,1]
YGumbel<-randomEstime[,2]
Khiplot(XGumbel,YGumbel,n,"Li","Xi","Khi-plot de la copule de Gumbel estimée",ylim=c(-0.2,0.8))

#K-plot de la copule de Gumbel estimé
KplotXRank<-rank(XGumbel,ties.method="random")
KplotYRank<-rank(YGumbel,ties.method="random")
BiCopKPlot(KplotXRank/n,KplotYRank/n,main="K-plot de la copule de Gumbel estimée")

#Dépendogramme

plot(KplotXRank/n,KplotYRank/n,xlab="Saint-Martin-en-Haut",ylab="Echirolles",main="Graphique des rangs pour la copule de Gumbel estimée")

#### METHODE 1: FONCTION D'UN PACKAGE DE R ####

gofCopula(fitGumbelCopula,cbind(Rx/(n+1),Ry/(n+1)))

#### METHODE 2: CODE MAISON ####


#Bootstrap paramétrique par la méthode CML

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
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitGumbelCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,fitGumbelCopula)
    rankX<-rank(randomXY[,1],ties.method="random")
    rankY<-rank(randomXY[,2],ties.method="random")
    Ui = rankX/(n+1)
    Vi = rankY/(n+1)
    uv <- data.frame(Ui,Vi)    
    fdrEmpirique<-function(u,v){
	return(EMPIRcop(u,v,para=uv))
    }

    # Estimation semi parametrique (CML) 
    tauKendall<-Kendall(randomXY[,1],randomXY[,2])
    paramGumbel<-1/(1-tauKendall$tau)
    myGumbelCopula<-gumbelCopula(param=paramGumbel)
    CML<-fitCopula(data=cbind(Ui,Vi),copula=myGumbelCopula,start=paramGumbel)
    thetaN<-coef(CML)
    fitGumbelCopula<-gumbelCopula(thetaN,dim=2)
    
    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitGumbelCopula))^2))
}

alpha = 0.05
L = sort(Dnk)[floor((1-alpha)*N)]
#Règle de décision
Dn > L
#p-value
sum(Dnk > Dn)/length(Dnk)


################################### Copule de Franck #########################

##############    1) METHODE DES MOMENTS  #######################

#Estimation du tau de Kendall par la méthode des moments

myFrankCopula<-frankCopula(param=1.5)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myFrankCopula,method="itau")
thetaN<-coef(CML)

fitFrankCopula<-frankCopula(thetaN,dim=2)

##############    2) METHODE CML  #######################


#Estimation du tau de Kendall par la méthode des moments

myFrankCopula<-frankCopula(param=1.5)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myFrankCopula)
thetaN<-coef(CML)

fitFrankCopula<-frankCopula(thetaN,dim=2)

##############    3) METHODE MV  #######################

myFrankMvd<-mvdc(copula=archmCopula(family="frank",param=1.5),margins=c("gumbel","gamma"),paramMargins=list(list(a=xpar_gumbel$estimate[1],b=xpar_gumbel$estimate[2]),list(shape=ypar_gamma$estimate[1],scale=ypar_gamma$estimate[2])))
myFrankCopula<-frankCopula(param=1.5)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myFrankCopula,method="itau")
thetaN<-coef(CML)
start<-c(xpar_gumbel$estimate[1],xpar_gumbel$estimate[2],ypar_gamma$estimate[1],ypar_gamma$estimate[2],thetaN)
fitTest<-fitMvdc(cbind(x,y),myFrankMvd,start=start,optim.control=list(trace=TRUE,maxit=2000))

########################## TEST GRAPHIQUE ET STATISTIQUE #######################

fitFrankCopula<-frankCopula(3.228,dim=2)
testFrankMvd<-mvdc(copula=fitFrankCopula,margins=c("gumbel","gamma"),paramMargins=list(list(a=5.012,b=1.815),list(shape=4.935,scale=1.339)))
randomEstime<-rMvdc(n,testFrankMvd)

#Khi-plot de la copule de Frank estimé

XFrank<-randomEstime[,1]
YFrank<-randomEstime[,2]
Khiplot(XFrank,YFrank,n,"Li","Xi","Khi-plot de la copule de Frank estimée",ylim=c(-0.2,0.8))

#K-plot de la copule de Frank estimé
KplotXRank<-rank(XFrank,ties.method="random")
KplotYRank<-rank(YFrank,ties.method="random")
BiCopKPlot(KplotXRank/n,KplotYRank/n,main="K-plot de la copule de Frank estimée")

#Dépendogramme

plot(KplotXRank/n,KplotYRank/n,xlab="Saint-Martin-en-Haut",ylab="Echirolles",main="Graphique des rangs pour la copule de Frank estimée")

#### METHODE 1: FONCTION D'UN PACKAGE DE R ####

gofCopula(fitFrankCopula,cbind(Rx/(n+1),Ry/(n+1)))

#### METHODE 2: CODE MAISON ####


#Bootstrap paramétrique par la méthode CML

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
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitFrankCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,fitFrankCopula)
    rankX<-rank(randomXY[,1],ties.method="random")
    rankY<-rank(randomXY[,2],ties.method="random")
    Ui = rankX/(n+1)
    Vi = rankY/(n+1)
    uv <- data.frame(Ui,Vi)    
    fdrEmpirique<-function(u,v){
	return(EMPIRcop(u,v,para=uv))
    }

    # Estimation semi parametrique (CML) 
    myFrankCopula<-frankCopula(param=1.5)
    CML<-fitCopula(data=cbind(Ui,Vi),copula=myFrankCopula)
    thetaN<-coef(CML)
    fitFrankCopula<-frankCopula(thetaN,dim=2)
    
    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitFrankCopula))^2))
}

alpha = 0.05
L = sort(Dnk)[floor((1-alpha)*N)]
#Règle de décision
Dn > L
#p-value
sum(Dnk > Dn)/length(Dnk)

################################### Copule gaussienne #########################

##############    1) METHODE DES MOMENTS  #######################

#Estimation par la méthode des moments (inversion du tau de Kendall)
#PS: Inversion rho de Spearman: sin(pi*cor(x,y,method="spearman")/6)*2

tauKendall<-Kendall(x,y)
paramNormal<-sin(pi*tauKendall$tau/2) 
myNormalCopula<-normalCopula(param=paramNormal)


##############    2) METHODE CML  #######################

myNormalCopula<-ellipCopula(family="normal",param=cor(x,y))
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myNormalCopula)
thetaN<-coef(CML)

fitNormalCopula<-normalCopula(thetaN,dim=2)

##############    3) METHODE MV  #######################

myNormalMvd<-mvdc(copula=ellipCopula(family="normal",param=0.5),margins=c("gumbel","gamma"),paramMargins=list(list(a=xpar_gumbel$estimate[1],b=xpar_gumbel$estimate[2]),list(shape=ypar_gamma$estimate[1],scale=ypar_gamma$estimate[2])))
tauKendall<-Kendall(x,y)
paramNormal<-sin(pi*tauKendall$tau/2) 
start<-c(xpar_gumbel$estimate[1],xpar_gumbel$estimate[2],ypar_gamma$estimate[1],ypar_gamma$estimate[2],paramNormal)
fitTest<-fitMvdc(cbind(x,y),myNormalMvd,start=start,optim.control=list(trace=TRUE,maxit=2000))

########################## TEST GRAPHIQUE ET STATISTIQUE #######################

fitNormalCopula<-normalCopula(0.518,dim=2)
testNormalMvd<-mvdc(copula=fitNormalCopula,margins=c("gumbel","gamma"),paramMargins=list(list(a=5.020,b=1.814),list(shape=4.984,scale=1.323)))
randomEstime<-rMvdc(n,testNormalMvd)

#Khi-plot de la copule gaussienne estimé

XNormal<-randomEstime[,1]
YNormal<-randomEstime[,2]
Khiplot(XNormal,YNormal,n,"Li","Xi","Khi-plot de la copule normale estimée",ylim=c(-0.2,0.8))

#K-plot de la copule gaussienne estimé
KplotXRank<-rank(XNormal,ties.method="random")
KplotYRank<-rank(YNormal,ties.method="random")
BiCopKPlot(KplotXRank/n,KplotYRank/n,main="K-plot de la copule normale estimée")

#Dépendogramme

plot(KplotXRank/n,KplotYRank/n,xlab="Saint-Martin-en-Haut",ylab="Echirolles",main="Graphique des rangs pour la copule normale estimée")

#### METHODE 1: FONCTION D'UN PACKAGE DE R ####

gofCopula(fitNormalCopula,cbind(Rx/(n+1),Ry/(n+1)))

#### METHODE 2: CODE MAISON ####

#Bootstrap paramétrique par la méthode CML

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
    tauKendall<-Kendall(randomXY[,1],randomXY[,2])
    paramNormal<-sin(pi*tauKendall$tau/2) 
    myNormalCopula<-ellipCopula(family="normal",param=paramNormal)
    CML<-fitCopula(data=cbind(Ui,Vi),copula=myNormalCopula,start=paramNormal)
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

################################### Copule de Student #########################

##############    1) METHODE DES MOMENTS  #######################

#Estimation du premier paramètre par la méthode des moments (inversion du tau de Kendall)
#Estimation du second paramètre par maximum de vraisemblance
tauKendall<-Kendall(x,y)
paramStudent<-sin(pi*tauKendall$tau/2) 
myStudentCopula<-tCopula(param=paramStudent,df.fixed=FALSE,dim=2)
params<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myStudentCopula,method="ml")
df<-coef(params)[2]
myStudentCopula<-tCopula(param=paramStudent,df=floor(df),dim=2)


##############    2) METHODE CML  #######################

#Estimation du premier paramètre par la méthode des moments (inversion du tau de Kendall)
#Estimation du second paramètre par maximum de vraisemblance

tauKendall<-Kendall(x,y)
paramStudent<-sin(pi*tauKendall$tau/2) 
myStudentCopula<-tCopula(param=paramStudent,df.fixed=FALSE,dim=2)
params<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myStudentCopula,method="ml")
df<-coef(params)[2]
myStudentCopula<-tCopula(param=paramStudent,df=df,dim=2)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myStudentCopula)

thetaN<-coef(CML)[1]
df<-coef(CML)[2]

fitStudentCopula<-tCopula(param=thetaN,dim=2,df=floor(df))

##############    3) METHODE MV  #######################

myStudentMvd<-mvdc(copula=ellipCopula(family="t",param=0.5),margins=c("gumbel","gamma"),paramMargins=list(list(a=xpar_gumbel$estimate[1],b=xpar_gumbel$estimate[2]),list(shape=ypar_gamma$estimate[1],scale=ypar_gamma$estimate[2])))
tauKendall<-Kendall(x,y)
paramStudent<-sin(pi*tauKendall$tau/2) 
myStudentCopula<-tCopula(param=paramStudent,df.fixed=FALSE,dim=2)
params<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myStudentCopula,method="ml")
df<-coef(params)[2]
start<-c(xpar_gumbel$estimate[1],xpar_gumbel$estimate[2],ypar_gamma$estimate[1],ypar_gamma$estimate[2],paramStudent,df)
fitTest<-fitMvdc(cbind(x,y),myStudentMvd,start=start,optim.control=list(trace=TRUE,maxit=2000))

########################## TEST GRAPHIQUE ET STATISTIQUE #######################

fitStudentCopula<-tCopula(0.513,dim=2)
testStudentMvd<-mvdc(copula=fitStudentCopula,margins=c("gumbel","gamma"),paramMargins=list(list(a=5.017,b=1.816),list(shape=4.997,scale=1.320)))
randomEstime<-rMvdc(n,testStudentMvd)

#Khi-plot de la copule de Student estimé

XStudent<-randomEstime[,1]
YStudent<-randomEstime[,2]
Khiplot(XStudent,YStudent,n,"Li","Xi","Khi-plot de la copule de Student estimée",ylim=c(-0.2,0.6))

#K-plot de la copule de Student estimé
KplotXRank<-rank(XStudent,ties.method="random")
KplotYRank<-rank(YStudent,ties.method="random")
BiCopKPlot(KplotXRank/n,KplotYRank/n,main="K-plot de la copule de Student estimée")

#Dépendogramme

plot(KplotXRank/n,KplotYRank/n,xlab="Saint-Martin-en-Haut",ylab="Echirolles",main="Graphique des rangs pour la copule de Student estimée")

######## EMV pour les paramètres lois gumbel et gamma #########

EMVgumbelb=sqrt(6)*sd(x)/pi
EMVgumbela=mean(x)-EMVgumbelb*(-digamma(1))
EMVgammashape= (mean(y)^2) /(mean(y^2)-mean(y)^2)
EMVgammascale= (mean(y^2) - mean(y)^2)/mean(y) 

##############    4) METHODE IFM  #######################

myMvd <- mvdc(copula = ellipCopula(family = "normal",param = 0.5),margins = c("gamma", "gamma"), 
              paramMargins = list(list(shape = xpar_gamma[1],scale =xpar_gamma[2] ), list(shape = ypar_gamma[1], scale = ypar_gamma[2])))
#n <- 200
#dat <- rMvdc(myMvd, n)
a.0 <- sin(cor(x, y, method = "kendall") * pi/2)
#Data with the gamma function for margin
udat <- cbind(pgamma(x,xpar_gamma$estimate[1],xpar_gamma$estimate[2]),pgamma(y,ypar_gamma$estimate[1],ypar_gamma$estimate[2]))

gumbel.cop <- gumbelCopula(5, dim=2)
clayton.cop <- claytonCopula(3,dim=2)
frank.cop <- frankCopula(3,dim=2)
normal.cop <- normalCopula(param=0.5,dim=2)
fit.ml.gumbel <- fitCopula(gumbel.cop, udat, method="ml")
fit.ml.clayton <- fitCopula(clayton.cop, udat, method="ml")
fit.ml.frank <- fitCopula(frank.cop, udat, method="ml")
fit.ml.normal <- fitCopula(normal.cop, udat, method="ml")

#Conclusion: 1) Normal, 2)Gumbel 3)Franck 4)Clayton

#Data with the normal distribution for margin
udat <- cbind(pnorm(x,xpar_norm$estimate[1],xpar_norm$estimate[2]),pnorm(y,ypar_norm$estimate[1],ypar_norm$estimate[2]))

gumbel.cop <- gumbelCopula(5, dim=2)
clayton.cop <- claytonCopula(3,dim=2)
frank.cop <- frankCopula(3,dim=2)
normal.cop <- normalCopula(0.5,dim=2)
fit.ml.gumbel <- fitCopula(gumbel.cop, udat, method="ml")
fit.ml.clayton <- fitCopula(clayton.cop, udat, method="ml")
fit.ml.frank <- fitCopula(frank.cop, udat, method="ml")
fit.ml.normal <- fitCopula(normal.cop, udat, method="ml")

#Conclusion: 1) Normal (101.76) , 2)Clayton (95.5) 3) Gumbel (80.5) 4)Franck (78.8)

#Data with the gumbel distribution for margin
udat <- cbind(pgumbel(x,xpar_gumbel$estimate[1],xpar_gumbel$estimate[2]),pgumbel(y,ypar_gumbel$estimate[1],ypar_gumbel$estimate[2]))


gumbel.cop <- gumbelCopula(5, dim=2)
clayton.cop <- claytonCopula(3,dim=2)
frank.cop <- frankCopula(3,dim=2)
normal.cop <- normalCopula(0.5,dim=2)
fit.ml.gumbel <- fitCopula(gumbel.cop, udat, method="ml")
fit.ml.clayton <- fitCopula(clayton.cop, udat, method="ml")
fit.ml.frank <- fitCopula(frank.cop, udat, method="ml")
fit.ml.normal <- fitCopula(normal.cop, udat, method="ml")

#Conclusion: 1) Normal ( 97.00561) , 2)umbel (94.24512)  3)Franck (77.87459) G4)Clayton (70.17998)

#Data with the weibull distribution for margin
udat <- cbind(pweibull(x,xpar_weibull$estimate[1],xpar_weibull$estimate[2]),pweibull(y,ypar_weibull$estimate[1],ypar_weibull$estimate[2]))

gumbel.cop <- gumbelCopula(5, dim=2)
clayton.cop <- claytonCopula(3,dim=2)
frank.cop <- frankCopula(3,dim=2)
normal.cop <- normalCopula(0.5,dim=2)
fit.ml.gumbel <- fitCopula(gumbel.cop, udat, method="ml")
fit.ml.clayton <- fitCopula(clayton.cop, udat, method="ml")
fit.ml.frank <- fitCopula(frank.cop, udat, method="ml")
fit.ml.normal <- fitCopula(normal.cop, udat, method="ml")

#Conclusion: 1) Normal (102.5) ,2)Clayton (89.6025)  3)Gumbel (82.23003)  4)Franck (80.8784) 

xpar_gumbel <- mledist(x, "gumbel",start=list(a=10,b=5))
ypar_gamma <- mledist(y, "gamma")

udat <- cbind(pgumbel(x,xpar_gumbel$estimate[1],xpar_gumbel$estimate[2]),pgamma(y,ypar_gamma$estimate[1],ypar_gamma$estimate[2]))

gumbel.cop <- gumbelCopula(5, dim=2)
clayton.cop <- claytonCopula(3,dim=2)
frank.cop <- frankCopula(3,dim=2)
normal.cop <- normalCopula(0.5,dim=2)
student.cop<- tCopula(0.5,dim=2)
fit.ml.gumbel <- fitCopula(gumbel.cop, udat, method="ml")
fit.ml.clayton <- fitCopula(clayton.cop, udat, method="ml")
fit.ml.frank <- fitCopula(frank.cop, udat, method="ml")
fit.ml.normal <- fitCopula(normal.cop, udat, method="ml")
fit.ml.student <- fitCopula(student.cop, udat, method="ml")

xpar_gumbel
ypar_gamma

summary(fit.ml.gumbel)
summary(fit.ml.clayton)
summary(fit.ml.franck)
summary(fit.ml.normal)
summary(fit.ml.student)
