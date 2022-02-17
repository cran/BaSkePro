


BaSkePro <- function (x)
{
  Man <- c(2,1.07,0,0,2,2,2,2,2,0,0)
  Atl <- c(1,0.49,0,1,1,1,1,1,0,0,0)
  Axi <-c(1,0.62,0,1,1,1,1,1,0,0,0)
  Cer <- c(5,0.45,0,5,5,5,5,5,0,0,0)
  Tho <- c( 13,0.53,0,0,13,13,13,0,0,0,0)
  Lum <- c( 6,0.51,6,6,6,6,6,6,6,0,0)
  Rib <- c(26,0.96,26,26,26,26,26,26,26,0,0)
  Sac <- c(1,0.4,0,0,1,1,1,0,0,0,0)
  Sca <- c(2,1.04,0,0,0,2,2,2,0,0,0)
  Hum <- c(2,1.12,0,0,0,2,2,2,2,2,0)
  Rad <- c(2,1.09,0,0,0,2,2,2,2,2,0)
  Pel <- c(2,1.02,0,0,2,2,2,0,0,0,0)
  Mca <- c(2,1.1,0,0,0,2,2,2,2,0,0)
  Fem <- c(2,1.15,0,0,0,0,2,2,2,2,2)
  Tib <- c(2,1.13,0,0,0,0,2,2,2,2,2)
  Mta <- c(2,1.1,0,0,0,0,2,2,2,2,0)
  config <- base::as.data.frame(base::rbind(Man,Atl,Axi,Cer,Tho, Lum,Rib,Sac,Sca,Hum,Rad,Mca,Pel,Fem,Tib,Mta))



  L <- function(theta,config,x)
  {
    value<-c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
    factor <- stats::dnorm(value,theta[1],0.25)
    factor <- factor/base::sum(factor)
    Comb <- factor[1]*config[,3]+factor[2]*config[,4]+factor[3]*config[,5]+factor[4]*config[,6]+factor[5]*config[,7]+factor[6]*config[,8]+factor[7]*config[,9]+factor[8]*config[,10]+factor[9]*config[,11]


    NME <- Comb*base::exp(-theta[2]*0.53/config[,2])
    MAU <- NME/config[,1]

    valid <- base::which(x>=0)
    PMAUsim <- MAU[valid]/sum(MAU[valid])


    PMAUobs <- x[valid,1]/sum(x[valid,1])

    x2 <- ((PMAUsim-PMAUobs)^2)/PMAUsim
    p <- 1/base::sum(x2)*base::length(base::which(x2<0.01))/base::length(valid)
    return (p)
  }


  Lbound<-c(-1, 0)
  Ubound<-c(1, 10)


  theta<-Lbound+ (Ubound-Lbound)*base::matrix(stats::runif(20000),10000,2)
  C<-stats::cov(theta)
  N<-100000
  theta<-base::matrix(0,N,2)
  p<-base::array(0,N)

  thetaOld <- Lbound +(Ubound-Lbound)*stats::runif(2)
  naccept <- 0
  pOld <- L(thetaOld,config,x)
  base::cat('Sampling with MCMC\n')


  cont<-1
  for (t in seq(1,N))
  {if (t>N*cont/20)
  {base::cat('Progress: ',cont*5,'%\n')
    cont <- cont+1}
    if (t<N/3 & base::floor(t/10000)==t/10000){
      C <- stats::cov(theta[(t-9999):t,])}
    out<-1
    while (out>0)
    {thetaNew<- MASS::mvrnorm(1,thetaOld,C)
    out<-0
    for (k in seq(1,2))
    {if (thetaNew[k]>Ubound[k]|thetaNew[k]<Lbound[k])
    {out=out+1}
    }
    }
    thetaNew<- Lbound + (Ubound-Lbound)*stats::runif(2)
    pNew<-L(thetaNew,config,x)
    r<-base::min(1,pNew/pOld)
    if (pOld==0|pNew==0)
    {r=0}
    u <- stats::runif(1)
    if (u<r)
    {theta[t,] = thetaNew
    p[t]=pNew
    naccept <- naccept + 1
    pOld = pNew
    thetaOld<-thetaNew}
    else
    {theta[t,]<- thetaOld
    p[t]<-pOld
    }
  }


  base::cat('Progress: 100%\n')
  base::cat()

  a <- base::round(naccept/N,3)
  b <- base::round(stats::median(theta[,1]),2)
  c <- c(base::round(stats::quantile(theta[,1],probs=0.025),2),base::round(stats::quantile(theta[,1],probs=0.975),2))
  d <- base::round(stats::median(theta[,2]),2)
  e <- c(base::round(stats::quantile(theta[,2],probs=0.025),2),base::round(stats::quantile(theta[,2],probs=0.975),2))

  Table_res <- base::data.frame(a,b,c,d,e)
  base::colnames(Table_res)<- c("Ratio of aceptance","Parameter alpha median","Parameter alpha 95%CI", "Parameter beta median","Parameter beta 95%CI")



  graphics::par(mfrow = c(2,1), mar = c(3,4,2,1))

  on.exit(graphics::par(graphics::par(mfrow = c(2,1), mar = c(3,4,2,1))), add = TRUE)


  graphics::hist(theta[,1], breaks=40 , xaxt = "n",probability=TRUE , ylab="Density", col=grDevices::rgb(1,0,0,0.5) , main = expression(paste(alpha," parameter",sep="" )))
  graphics::axis(side=1, at=seq(-1,1, 0.25))
  graphics::abline(v = stats::median(theta[,1]),col = "red",lwd = 4)
  graphics::mtext("Transport strategy", 1, line=2)

   graphics::hist(theta[,2], breaks=50  , xaxt = "n", probability =TRUE,ylab = "Density", col=grDevices::rgb(0,0,1,0.5) , main = expression(paste(beta," parameter",sep="" )))
  graphics::axis(side=1, at=seq(0,10))
  graphics::abline(v = stats::median(theta[,2]),col = "blue",lwd = 4)
  graphics::mtext("Attrition degrees", 1,  line=2)


   return(Table_res)
}

