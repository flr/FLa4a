library(FLa4a)
data(rfLen)
# Make an empty cor matrix
cm <- diag(c(1,1,1))
# k and linf are negatively correlated while t0 is independent
cm[1,2] <- cm[2,1] <- -0.5
# scale cor to var using CV=0.2
cv <- 0.2
p <- c(linf=60, k=0.09, t0=-0.01)
vc <- matrix(1, ncol=3, nrow=3)
l <- vc
l[1,] <- l[,1] <- p[1]*cv
k <- vc
k[,2] <- k[2,] <- p[2]*cv
t <- vc
t[3,] <- t[,3] <- p[3]*cv
mm <- t*k*l
diag(mm) <- diag(mm)^2
mm <- mm*cm
# check that we have the intended correlation
all.equal(cm, cov2cor(mm))
# Create the a4aGr object as before but now we also include the vcov argument for the variance-covariance matrix
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=p["linf"], k=p["k"], t0=p["t0"], units=c("cm","ano-1","ano")), vcov=mm)

# transform idx and stk
aIdx <- l2a(rfTrawl.idx, vbObj)
aStk <- l2a(rfLen.stk, vbObj)

# check 
sum(catch.n(rfTrawl.idx), na.rm=T)==sum(catch.n(aIdx), na.rm=T)
sum(catch.n(rfTrawl.idx)*catch.wt(rfTrawl.idx), na.rm=T)==sum(catch.n(aIdx)*catch.wt(aIdx), na.rm=T)

sum(catch.n(rfLen.stk), na.rm=T)==sum(catch.n(aStk), na.rm=T)
sum(catch.n(rfLen.stk)*catch.wt(rfLen.stk), na.rm=T)==sum(catch.n(aStk)*catch.wt(aStk), na.rm=T)

