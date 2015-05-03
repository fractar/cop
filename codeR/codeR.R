
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

Khiplot<-function(x,y,n,xtitle,ytitle,title) {

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
plot(Li[ind],Xi[ind],main=title,col="blue",xlab=xtitle,ylab=ytitle)
}
#postscript("chi_plot_empir.ps")
Khiplot(x,y,n,"Li","Xi","Khi-plot empirique")
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

#Khi-plot de la copule de Clayton estimé

number = 1000
randomEstime<-rCopula(number,myClaytonCopula)
XClayton<-randomEstime[,1]
YClayton<-randomEstime[,2]
Khiplot(XClayton,YClayton,number)

#K-plot de la copule de Clayton estimé
KplotXRank<-rank(XClayton,ties.method="random")
KplotYRank<-rank(YClayton,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number))

#Bootstrap paramétrique par méthode des moments

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
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, myClaytonCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,myClaytonCopula)
    rankX<-rank(randomXY[,1],ties.method="random")
    rankY<-rank(randomXY[,2],ties.method="random")
    Ui = rankX/(n+1)
    Vi = rankY/(n+1)
    uv <- data.frame(Ui,Vi)    
    fdrEmpirique<-function(u,v){
	return(EMPIRcop(u,v,para=uv))
    }
    tauKendall=Kendall(randomXY[,1],randomXY[,2])

    # Estimation par méthode des moments
    paramClayton<-2*tauKendall$tau/(1-tauKendall$tau)
    myClaytonCopula<-claytonCopula(param=paramClayton)
    
    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, myClaytonCopula))^2))
}

alpha = 0.05
L = sort(Dnk)[floor((1-alpha)*N)]
#Règle de décision
Dn > L
#p-value
sum(Dnk > Dn)/length(Dnk)


##############    2) METHODE CML  #######################

#Estimation du tau de Kendall par la méthode des moments
tauKendall<-Kendall(x,y)
paramClayton<-2*tauKendall$tau/(1-tauKendall$tau)

# Estimation semi parametrique (CML) 

myClaytonCopula<-claytonCopula(param=paramClayton)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myClaytonCopula,start=paramClayton)
thetaN<-coef(CML)

fitClaytonCopula<-claytonCopula(thetaN,dim=2)


#Khi-plot de la copule de Clayton estimé

number = 1000
randomEstime<-rCopula(number,fitClaytonCopula)
XClayton<-randomEstime[,1]
YClayton<-randomEstime[,2]
Khiplot(XClayton,YClayton,number,"Li","Xi","Khi-plot de la copule de Clayton estimée (méthode CML)")

#K-plot de la copule de Clayton estimé
KplotXRank<-rank(XClayton,ties.method="random")
KplotYRank<-rank(YClayton,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number),main="K-plot de la copule de Clayton estimée (méthode CML)")

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

################################### Copule gaussienne #########################

##############    1) METHODE DES MOMENTS  #######################

#Estimation par la méthode des moments (inversion du tau de Kendall)
#PS: Inversion rho de Spearman: sin(pi*cor(x,y,method="spearman")/6)*2

tauKendall<-Kendall(x,y)
paramNormal<-sin(pi*tauKendall$tau/2) 
myNormalCopula<-normalCopula(param=paramNormal)

#Khi-plot de la copule gaussienne estimé

number = 1000
randomEstime<-rCopula(number,myNormalCopula)
XNormal<-randomEstime[,1]
YNormal<-randomEstime[,2]
Khiplot(XNormal,YNormal,number)

#K-plot de la copule gaussienne estimé
KplotXRank<-rank(XNormal,ties.method="random")
KplotYRank<-rank(YNormal,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number))


#Bootstrap paramétrique par méthode des moments

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
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, myNormalCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,myNormalCopula)
    rankX<-rank(randomXY[,1],ties.method="random")
    rankY<-rank(randomXY[,2],ties.method="random")
    Ui = rankX/(n+1)
    Vi = rankY/(n+1)
    uv <- data.frame(Ui,Vi)    
    fdrEmpirique<-function(u,v){
	return(EMPIRcop(u,v,para=uv))
    }

    # Estimation par la méthode des moments
    tauKendall<-Kendall(randomXY[,1],randomXY[,2])
    paramNormal<-sin(pi*tauKendall$tau/2) 
    myNormalCopula<-ellipCopula(family="normal",param=paramNormal)
 
    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, myNormalCopula))^2))
}

alpha = 0.05
L = sort(Dnk)[floor((1-alpha)*N)]
#Règle de décision
Dn > L
#p-value
sum(Dnk > Dn)/length(Dnk)

##############    2) METHODE CML  #######################

myNormalCopula<-ellipCopula(family="normal",param=cor(x,y))
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myNormalCopula)
thetaN<-coef(CML)

fitNormalCopula<-normalCopula(thetaN,dim=2)

#Khi-plot de la copule gaussienne estimé

number = 1000
randomEstime<-rCopula(number,fitNormalCopula)
XNormal<-randomEstime[,1]
YNormal<-randomEstime[,2]
Khiplot(XNormal,YNormal,number,"Li","Xi","Khi-plot de la copule normale estimée (méthode CML)")

#K-plot de la copule gaussienne estimé
KplotXRank<-rank(XNormal,ties.method="random")
KplotYRank<-rank(YNormal,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number),main="K-plot de la copule normale estimée (méthode CML)")

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

#Khi-plot de la copule de Student estimé

number = 1000
randomEstime<-rCopula(number,myStudentCopula)
XStudent<-randomEstime[,1]
YStudent<-randomEstime[,2]
Khiplot(XStudent,YStudent,number)

#K-plot de la copule de Student estimé
KplotXRank<-rank(XStudent,ties.method="random")
KplotYRank<-rank(YStudent,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number))


#Bootstrap paramétrique par méthode des moments

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
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, myStudentCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,myStudentCopula)
    rankX<-rank(randomXY[,1],ties.method="random")
    rankY<-rank(randomXY[,2],ties.method="random")
    Ui = rankX/(n+1)
    Vi = rankY/(n+1)
    uv <- data.frame(Ui,Vi)    
    fdrEmpirique<-function(u,v){
	return(EMPIRcop(u,v,para=uv))
    }

    # Estimation par la méthode des moments

    tauKendall<-Kendall(randomXY[,1],randomXY[,2])
    paramStudent<-sin(pi*tauKendall$tau/2) 
    myStudentCopula<-tCopula(param=paramStudent,df.fixed=FALSE,dim=2)
    params<-fitCopula(data=cbind(Ui,Vi),copula=myStudentCopula,method="ml")
    df<-coef(params)[2]
    myStudentCopula<-tCopula(param=paramStudent,df=floor(df),dim=2)
 
    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, myStudentCopula))^2))
}

alpha = 0.05
L = sort(Dnk)[floor((1-alpha)*N)]
#Règle de décision
Dn > L
#p-value
sum(Dnk > Dn)/length(Dnk)


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

#Khi-plot de la copule de Student estimé

number = 1000
randomEstime<-rCopula(number,fitStudentCopula)
XStudent<-randomEstime[,1]
YStudent<-randomEstime[,2]
Khiplot(XStudent,YStudent,number,"Li","Xi","Khi-plot de la copule de Student estimée (méthode CML)")

#K-plot de la copule de Student estimé
KplotXRank<-rank(XStudent,ties.method="random")
KplotYRank<-rank(YStudent,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number),main="K-plot de la copule de Student estimée (méthode CML)")

#### METHODE 1: FONCTION D'UN PACKAGE DE R ####

gofCopula(fitStudentCopula,cbind(Rx/(n+1),Ry/(n+1)))

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
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitStudentCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,fitStudentCopula)
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
    paramStudent<-sin(pi*tauKendall$tau/2) 
    myStudentCopula<-tCopula(param=paramStudent,df.fixed=FALSE,dim=2)
    params<-fitCopula(data=cbind(Ui,Vi),copula=myStudentCopula,method="ml")
    df<-coef(params)[2]
    myStudentCopula<-tCopula(param=paramStudent,df=df,dim=2)
    CML<-fitCopula(data=cbind(Ui,Vi),copula=myStudentCopula,start=c(paramStudent,df))
    thetaN<-coef(CML)[1]
    df<-coef(CML)[2]
    fitStudentCopula<-tCopula(param=thetaN,dim=2,df=floor(df))
    
    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, fitStudentCopula))^2))
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

#Khi-plot de la copule de Gumbel estimé

number = 1000
randomEstime<-rCopula(number,myGumbelCopula)
XGumbel<-randomEstime[,1]
YGumbel<-randomEstime[,2]
Khiplot(XGumbel,YGumbel,number)

#K-plot de la copule de Gumbel estimé
KplotXRank<-rank(XGumbel,ties.method="random")
KplotYRank<-rank(YGumbel,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number))


#Bootstrap paramétrique par méthode des moments

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
Dn = sum((fdrEmpirique(Ui,Vi) - pCopula(UV, myGumbelCopula))^2)
Dnk = NULL
for (k in 1:N) {
    print(k)
    randomXY <-rCopula(n,myGumbelCopula)
    rankX<-rank(randomXY[,1],ties.method="random")
    rankY<-rank(randomXY[,2],ties.method="random")
    Ui = rankX/(n+1)
    Vi = rankY/(n+1)
    uv <- data.frame(Ui,Vi)    
    fdrEmpirique<-function(u,v){
	return(EMPIRcop(u,v,para=uv))
    }

    # Estimation par la méthode des moments
    tauKendall<-Kendall(randomXY[,1],randomXY[,2])
    paramGumbel<-1/(1-tauKendall$tau) 
    myGumbelCopula<-gumbelCopula(param=paramGumbel)

    UV<-cbind(Ui,Vi)
    Dnk = c(Dnk,sum((fdrEmpirique(Ui,Vi) - pCopula(UV, myGumbelCopula))^2))
}

alpha = 0.05
L = sort(Dnk)[floor((1-alpha)*N)]
#Règle de décision
Dn > L
#p-value
sum(Dnk > Dn)/length(Dnk)

##############    2) METHODE CML  #######################

#Estimation du tau de Kendall par la méthode des moments
tauKendall<-Kendall(x,y)
paramGumbel<-1/(1-tauKendall$tau)

myGumbelCopula<-gumbelCopula(param=paramGumbel)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myGumbelCopula)
thetaN<-coef(CML)

fitGumbelCopula<-gumbelCopula(thetaN,dim=2)

#Khi-plot de la copule de Gumbel estimé

number = 1000
randomEstime<-rCopula(number,fitGumbelCopula)
XGumbel<-randomEstime[,1]
YGumbel<-randomEstime[,2]
Khiplot(XGumbel,YGumbel,number,"Li","Xi","Khi-plot de la copule de Gumbel estimée (méthode CML)")

#K-plot de la copule de Gumbel estimé
KplotXRank<-rank(XGumbel,ties.method="random")
KplotYRank<-rank(YGumbel,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number),main="K-plot de la copule de Gumbel estimée (méthode CML)")

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

#Khi-plot de la copule de Frank estimé

number = 1000
randomEstime<-rCopula(number,fitFrankCopula)
XFrank<-randomEstime[,1]
YFrank<-randomEstime[,2]
Khiplot(XFrank,YFrank,number)

#K-plot de la copule de Frank estimé
KplotXRank<-rank(XFrank,ties.method="random")
KplotYRank<-rank(YFrank,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number))


#Bootstrap paramétrique par méthode des moments

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

    # Estimation par la méthode des moments
    myFrankCopula<-frankCopula(param=1.5)
    CML<-fitCopula(data=cbind(Ui,Vi),copula=myFrankCopula,method="itau")
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


##############    2) METHODE CML  #######################


#Estimation du tau de Kendall par la méthode des moments

myFrankCopula<-frankCopula(param=1.5)
CML<-fitCopula(data=cbind(Rx/(n+1),Ry/(n+1)),copula=myFrankCopula)
thetaN<-coef(CML)

fitFrankCopula<-frankCopula(thetaN,dim=2)

#Khi-plot de la copule de Frank estimé

number = 1000
randomEstime<-rCopula(number,fitFrankCopula)
XFrank<-randomEstime[,1]
YFrank<-randomEstime[,2]
Khiplot(XFrank,YFrank,number,"Li","Xi","Khi-plot de la copule de Frank estimée (méthode CML)")

#K-plot de la copule de Frank estimé
KplotXRank<-rank(XFrank,ties.method="random")
KplotYRank<-rank(YFrank,ties.method="random")
BiCopKPlot(KplotXRank/(number),KplotYRank/(number),main="K-plot de la copule de Frank estimée (méthode CML)")

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

