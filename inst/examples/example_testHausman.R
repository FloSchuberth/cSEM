# It is still based on composites!!!!!!!!!!!!1



# Indicator names
Indnames <- c(paste('x',1:12,sep=''),paste('y',1:6,sep=''))

loadXI1  = c(0.8,0.7,0.75)
loadXI2  = c(0.75,0.7,0.75)
loadXI3  = c(0.8,0.8,0.8)
loadXI4  = c(0.7,0.7,0.7)
loadETA1 = c(0.6,0.6,0.6)
loadETA2 = c(0.85,0.75,0.75)



# all loadings in one matrix
LOAD=Matrix::bdiag(loadXI1,loadXI2,loadXI3,loadXI4,loadETA1,loadETA2)
LOAD=as.matrix(LOAD)

# path coefficients
B<-matrix(data=c(0, 0.25, 0.5,0),ncol=2, nrow=2, byrow=TRUE)
Gamma<-matrix(data=c(-0.3, 0.5, 0  ,0,0,0,0.5,0.25),ncol=4, nrow=2,byrow=T)
I=diag(ncol(B))
IB=solve(I-B)

# VCV between the exogenous composites. Exogenous refers here to the structural model
PHI=matrix(c(
  1, 0.5, 0.5, 0.5, 
  0.5,1,0.5,0.5,
  0.5,0.5,1,0.5,
  0.5,0.5,0.5,1),ncol=4)

# VCV of the structural residuals
# PSI= matrix(c(0.5189,-0.0295,
#               -0.0295,0.1054),ncol=2)

# PSI= matrix(c(0.5189,-0.5,
#               -0.5,0.1054),ncol=2)

# Correlation matrix of the error terms of the structural model
CORREndo=sqrt(seq(0,0.6,by = 0.1))

Sigma=lapply(CORREndo, function(CorrEndo){
  
  ENDOCORR=matrix(c(1,CorrEndo,CorrEndo,1),ncol=2)
  PSI=(I-B)%*%ENDOCORR%*%t(I-B)-Gamma%*%PHI %*%t(Gamma)
  
  # PSI= matrix(c(0.5189,Psi12,
  #               Psi12,0.1054),ncol=2)
  # 
  # #Correlation between the composites
  # ENDOCORR1=IB%*%Gamma%*%PHI%*%t(Gamma)%*%t(IB)+IB%*%PSI%*%t(IB)
  # diag(ENDOCORR)=1
  EXOENDOCORR=PHI%*%t(Gamma)%*%t(IB)
  
  COMPCORR=cbind(PHI,EXOENDOCORR)
  COMPCORR=rbind(COMPCORR,cbind(t(EXOENDOCORR),ENDOCORR))
  
  # Sigma without the (co-)variances of the error terms
  Sigmawithout=LOAD%*%COMPCORR%*%t(LOAD)
  
  # Matrix with block correlations
  # IndCORR=Matrix::bdiag(replicate(6,Blockcorr,simplify = F))
  # IndCORR=as.matrix(IndCORR)
  # 
  # Sigmawithout[IndCORR!=0]=IndCORR[IndCORR!=0]
  diag(Sigmawithout)=1
  Sigmacomp=Sigmawithout
  dimnames(Sigmacomp)=list(Indnames,Indnames)
  # ensure that matrix is symmetric
  Sigmacomp[lower.tri(Sigmacomp)]=t(Sigmacomp)[lower.tri(Sigmacomp)]
  Sigmacomp
})

# check Sigma matrices
lapply(Sigma,matrixcalc::is.positive.definite)

library(cSEM)
## model specifictaion ----
model='
ETA1 ~ ETA2 + XI1 + XI2
ETA2 ~ETA1 + XI3 +XI4

XI1 =~ x1 + x2+x3
XI2 =~ x4+x5+x6
XI3 =~ x7 + x8 + x9
XI4 =~ x10 + x11 + x12
ETA1 =~ y1 + y2 + y3
ETA2 =~ y4+y5+y6

'

model1='
ETA1 ~ ETA2 + XI1 + XI2
ETA2 ~ XI3 +XI4

XI1 =~ x1 + x2+x3
XI2 =~ x4+x5+x6
XI3 =~ x7 + x8 + x9
XI4 =~ x10 + x11 + x12
ETA1 =~ y1 + y2 + y3
ETA2 =~ y4+y5+y6

'

library(matrixpls);library(cSEM)
# Comparison matrixpls vs csem
# OLS
dataFis=lapply(Sigma,function(x){MASS::mvrnorm(n=300,mu = rep(0,18),Sigma = x,empirical = T)})
OLScsem=csem(.data = dataFis[[1]],.model = model)
# OLSmat=matrixpls(cor(dataFis[[1]]), model=model,weightFun = weightFun.pls,
# parametersInner = estimator.ols)

summarize(OLScsem)

# calculate difference
OLScsem$Estimates$Path_estimates[dimnames(attr(OLSmat,'inner'))[[1]],dimnames(attr(OLSmat,'inner'))[[2]]]-
  attr(OLSmat,'inner')

# 2SLS
TwoSLScsem=csem(.data = dataFis[[1]],.model = model,.approach_paths= "2SLS",
                .instruments = list(ETA1=c('XI1','XI2','XI3','XI4'), ETA2=c('XI1','XI2','XI3','XI4')))

summarize(TwoSLScsem)


testh <- testHausman(TwoSLScsem, .R = 200)

TwoSLScsem=csem(.data = dataFis[[1]],.model = model1,.approach_paths= "2SLS",
                .instruments = list(ETA1=c('XI1','XI2','XI3','XI4')))

TwoSLScsem$Estimates$Path_estimates
TwoSLScsem$Estimates$Loading_estimates

