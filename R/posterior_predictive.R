BetaBinom <- Vectorize(function(rp){
  log.val <- lchoose(np, rp) + lbeta(rp+a+r,b+n-r+np-rp) - lbeta(a+r,b+n-r)
  return(exp(log.val))
})

#####################
# Small sample
#####################

# Example 1
n <- 8; r <- 2; a <- 0; b <- 0; np <- 40
plot(0:10,BetaBinom(0:10),type="b",xlab="r*",ylab="P(R=r*|Data)", main = "Posterior predictive: a=0, b=0",cex.axis= 1.5,cex.lab=1.5,lwd=4, ylim=c(0,0.3))
AUC00=trapz(1:5,BetaBinom(1:5))# 5>=1; 0.3930899
AUC00
AUC00=trapz(1:4,BetaBinom(1:4)) #0.3156377
AUC00
AUC00=trapz(1:3,BetaBinom(1:3)) # 0.2254646
AUC00
AUC00=trapz(1:2,BetaBinom(1:2)) # 0.1208757
AUC00

n <- 8; r <- 1; a <- 0; b <- 1; np <- 40
#plot(0:10,BetaBinom(0:10),type="b",xlab="r*",ylab="P(R=r*|Data)", main = "Posterior predictive: a=0, b=1",cex.axis= 1.5,cex.lab=1.5,lwd=4, ylim=c(0,0.3))
AUC01=trapz(1:5,BetaBinom(1:5))# 5>=1; 0.3930899
AUC01
AUC01=trapz(1:4,BetaBinom(1:4))
AUC01
AUC01=trapz(1:3,BetaBinom(1:3))
AUC01
AUC01=trapz(1:2,BetaBinom(1:2))
AUC01
n <- 8; r <- 1; a <- 0; b <- 2; np <- 40
#plot(0:10,BetaBinom(0:10),type="b",xlab="r*",ylab="\'b' P(R=r*|Data)", main = "Posterior predictive: a=0, b=2",cex.axis= 1.5,cex.lab=1.5,lwd=4, ylim=c(0,0.3))
AUC02=trapz(1:5,BetaBinom(1:5))# 5>=1; 0.3930899
AUC02
AUC02=trapz(1:4,BetaBinom(1:4))
AUC02
AUC02=trapz(1:3,BetaBinom(1:3))
AUC02
AUC02=trapz(1:2,BetaBinom(1:2))
AUC02
n <- 8; r <- 1; a <- 0; b <- 3; np <- 40
#plot(0:10,BetaBinom(0:10),type="b",xlab="r*",ylab="P(R=r*|Data)", main = "Posterior predictive: a=0, b=3",cex.axis= 1.5,cex.lab=1.5,lwd=4, ylim=c(0,0.3))
AUC03=trapz(1:5,BetaBinom(1:5))# 5>=1; 0.3930899
AUC03
AUC03=trapz(1:4,BetaBinom(1:4))
AUC03
AUC03=trapz(1:3,BetaBinom(1:3))
AUC03
AUC03=trapz(1:2,BetaBinom(1:2))
AUC03
n <- 8; r <- 1; a <- 0; b <- 4; np <- 40
plot(0:10,BetaBinom(0:10),type="b",xlab="r*",ylab="P(R=r*|Data)", main = "Posterior predictive: a=0, b=4",cex.axis= 1.5,cex.lab=1.5,lwd=4, ylim=c(0,0.3))
AUC04=trapz(1:5,BetaBinom(1:5))# 5>=1; 0.3930899
AUC04
AUC04=trapz(1:4,BetaBinom(1:4))
AUC04
AUC04=trapz(1:3,BetaBinom(1:3))
AUC04
AUC04=trapz(1:2,BetaBinom(1:2))
AUC04


trapz(BetaBinom(1:39), 1:40)
x<-0:39
y<-BetaBinom(1:39)
id<-order(x)
sum(diff(x[id])*rollmean(y[id],2))
