data(ple4)
data(ple4.indices)

# creating idx 1
bioidx <- FLIndex(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])))
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("startf","endf")] <- c(0,0)

# creating idx 2
stkn <- stock.n(ple4)*exp(-z(ple4))*0.5
bioidx2 <- FLIndex(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])))
index(bioidx2) <- quantSums(stkn*stock.wt(ple4))*0.001
index(bioidx2) <- index(bioidx2)*exp(rnorm(index(bioidx2), sd=0.1))
range(bioidx2)[c("startf","endf")] <- c(0,1)

# FLIndices
idxs <- ple4.indices
idxs$b1 <- bioidx
idxs$b2 <- bioidx2

# fitting the model
fit <- sca(ple4, idxs, qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3), ~1, ~1))
res <- residuals(fit, ple4, idxs)
plot(res)

sca(ple4, idxs[-1], qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3), ~1, ~1)[-1])
sca(ple4, idxs[-2], qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3), ~1, ~1)[-2])
sca(ple4, idxs[-3], qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3), ~1, ~1)[-3])
sca(ple4, idxs[-4], qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3), ~1, ~1)[-4])
sca(ple4, idxs[-5], qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3), ~1, ~1)[-5])


