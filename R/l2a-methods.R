#' @title Method to convert age to length 
#' @name l2a 
#' @rdname l2a 
#' @aliases l2a l2a-methods l2a,FLQuant,a4aGr-method
#' @param object a \code{FLQuant} object
#' @param model a \code{a4aGr} object
#' @param stat the aggregation statistic, must be \"mean\" or \"sum\"
#' @return a \code{FLQuant} object
#' @examples
#' # red fish
#' # M=0.05; Linf=58.5, k=0.086
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(50, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
#' vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
#' data(rfLen)
#' cth <- catch.n(rfLen.stk) 
#' # both with iter=1
#' cthA1 <- l2a(cth, vbObj)
#' # both with iter=n
#' cthA2 <- l2a(propagate(cth,10), mvrnorm(10, vbObj))
#' # mod: iter=n, data: iter=1
#' cthA3 <- l2a(cth, mvrnorm(10, vbObj))
#' # mod: iter=1, data: iter=n
#' cthA4 <- l2a(propagate(cth,10), vbObj)
#' # converting a stock object
#' rfAge.stk <- l2a(rfLen.stk, vbObj)
#' rfAge.stk <- l2a(rfLen.stk, mvrnorm(10, vbObj))
#' rfAge.stk <- l2a(propagate(rfLen.stk, 10), vbObj)
#' # converting a index object
#' rfAge.idx <- l2a(rfTrawl.idx, vbObj)
#' rfAge.idx <- l2a(rfTrawl.idx, mvrnorm(10, vbObj))
#' rfAge.idx <- l2a(propagate(rfTrawl.idx, 10), vbObj)

# l2a
setGeneric("l2a", function(object, model, ...) standardGeneric("l2a"))
setMethod("l2a", c("FLQuant", "a4aGr"), function(object, model, stat="sum", weights=FLQuant(1, dimnames=dimnames(object)), ...){
	# constants
	cat("Converting lengths to ages ...\n")
	dnms <- dimnames(object)
	if(!all.equal(dnms, dimnames(weights))) stop("Weights must have the same dimensions as the data.")
	len <- as.numeric(dnms[[1]])
	mit <- niters(model)
	qit <- length(dnms[[6]]) 
	if(mit>1 & qit>1 & mit!=qit) stop("Can not operate with diferent iterations in each object.")

	# model
	mod <- FLModelSim(model=grInvMod(model), params=params(model), vcov=vcov(model), distr=distr(model))
	age <- floor(predict(mod, len=len))

	# aggregates all lengths above linf
	age <- apply(age, 2, function(x){
		x[which.max(x[!is.infinite(x)]):length(x)] <- max(x[!is.infinite(x)], na.rm=T)
		x
	})

    # Blow up object to have same iters as age 
    if (qit < mit){
        object <- propagate(object,iter=mit)
    }

	# slicing and aggregating
    qname <- names(dimnames(object))[1]
    names(dimnames(age))[1] <- qname
    agem <- melt(age, value.name="age")
    odf <- as.data.frame(object)
    # Force to be the same - necessary for join
    odf$len <- as.numeric(odf$len)
    agem$len <- as.numeric(agem$len)

    # Join the length based data with age-length key
    if (mit==dim(object)[6]){
        odf <- left_join(odf,agem, by=c("iter",qname))
    }
    # If age matrix has only 1 iter
    if (mit==1 & (dim(object)[6] > 1)){
        odf <- left_join(odf,agem[,c(qname,"age")], by=c(qname))
    }

    
	if(stat=="sum"){
        out <- summarise(group_by(odf, age, year, unit, season, area, iter), data = sum(data, na.rm = TRUE))
    }

	if(stat=="mean"){
        # Blow up weights to match same iters as growth model
        if (qit < mit){
            weights <- propagate(weights,iter=mit)
        }
        wtdf <- as.data.frame(weights)
        names(wtdf)[names(wtdf)=="data"] <- "weight"
        odf <- left_join(odf,wtdf, by=c(qname,"year","unit","season","area","iter"))
        out <- summarise(group_by(odf, age, year, unit, season, area, iter), data = sum(data * (weight/sum(weight, na.rm = TRUE)), na.rm = TRUE))
    }

    # Add missing ages as otherwise not contiguous
    all_ages <- min(out$age):max(out$age)
    missing_ages <- all_ages[!(all_ages %in% unique(out$age))]
    extra_rows <- expand.grid(age=missing_ages, year = unique(out$year), unit=unique(out$unit), season=unique(out$season), area=unique(out$area), iter=unique(out$iter), data=0)
    out <- rbind(out,extra_rows)

    # Turn the out dataframe into an FLQuant
    # using as* takes ages
    # flq <- as.FLQuant(outdf)
    # Use reshape2 instead
    outa <- acast(out, age~year~unit~season~area~iter, value.var = "data")
    names(dimnames(outa)) <- c("age","year","unit","season","area","iter")
    flq <- as.FLQuant(outa)
    flq[is.na(flq)] <- 0.0 # hack, where are these NAs coming from?

	units(flq) <- units(object)
    return(flq)
})


#' @rdname l2a 
#' @aliases l2a,FLStockLen,a4aGr-method
setMethod("l2a", c("FLStockLen", "a4aGr"), function(object, model, plusgroup="missing", ...){
	warning("Individual weights, M and maturity will be averaged accross lengths, everything else will be summed. If this is not what you want, you'll have to deal with these slots by hand.")
	args <- list(...)
	args$model <- model
	args$FUN <- "l2a"	
	dnms <- dimnames(catch.n(object))
	
	# trick to get the relevant slots
	lst <- qapply(object, function(x) list(x))
	lst <- lapply(lst, "[[", 1)	
	nms <- qapply(object, quant)

	# first the sums	
	isLen <- lapply(nms, "==", "len")
	# aggregated are not lenghts
	isLen[c("catch", "stock", "landings", "discards")][] <- FALSE
	# whatelse
	isLen[c("catch.wt", "stock.wt", "landings.wt", "discards.wt", "mat", "m", "m.spwn", "harvest.spwn", "harvest")][] <- FALSE
	lst[unlist(isLen)] <- lapply(lst[unlist(isLen)], function(x){
		args$object <- x
		do.call("l2a", args)
	})
	
	# now the means	
	isLen <- lapply(nms, "==", "len")
	# aggregated are not lenghts
	isLen[c("catch", "stock", "landings", "discards")][] <- FALSE
	# whatelse
	isLen[c("catch.n", "stock.n", "landings.n", "discards.n", "catch.wt", "stock.wt", "landings.wt", "discards.wt")][] <- FALSE
	args$stat <- "mean"
	lst[unlist(isLen)] <- lapply(lst[unlist(isLen)], function(x){
		args$object <- x
		do.call("l2a", args)
	})
	# now the weighted means
	#==============================
	# WARNING: IT'S NOT WORKING YET	
	#==============================
	isLen <- lapply(nms, "==", "len")
	# aggregated are not lenghts
	isLen[c("catch", "stock", "landings", "discards")][] <- FALSE
	# whatelse
	isLen[c("catch.n", "stock.n", "landings.n", "discards.n", "mat", "m", "m.spwn", "harvest.spwn", "harvest")][] <- FALSE
	args$stat <- "mean"
	lst[unlist(isLen)] <- lapply(lst[unlist(isLen)], function(x){
		args$object <- x
		do.call("l2a", args)
	})

	# set the plus group on the first non continuous age
	isLen[] <- TRUE
	isLen[c("catch", "stock", "landings", "discards")][] <- FALSE
	lst00 <- lst[unlist(isLen)]

	if(missing(plusgroup)){
		lst0 <- lapply(lst00, function(x){
			v <- as.numeric(dimnames(x)[[1]])
			min(v[-length(v)][(v[-length(v)]-v[-1])< -1])
		}) 
		plusgroup <- min(unlist(lst0))
	}
	
	# spagethi ...	
	lst1 <- lapply(lst00, function(x){
		v <- as.numeric(dimnames(x)[[1]])
		min(v):max(v)	
	}) 
	dnms$age <- sort(unique(unlist(lst1)))
	dnms$len <- NULL
	flq <- FLQuant(dimnames=dnms)
	lst00$flq <- flq
	lst00 <- mcf(lst00)
	lst[unlist(isLen)] <- lst00[-match("flq", names(lst00))]

	# create the stock object
	lst <- lapply(lst, "quant<-", "age")
	lst$name <- object@name
	lst$desc <- object@desc
	stk <- do.call("FLStock", lst)
	# set the plus group on the first non continuous age
	stk <- setPlusGroup(stk, plusgroup, na.rm=T)
	# set harvest units
	units(harvest(stk)) <- units(object@harvest)
	})

#' @rdname l2a 
#' @aliases l2a,FLIndex,a4aGr-method
setMethod("l2a", c("FLIndex", "a4aGr"), function(object, model, ...){
	warning("Catch in numbers will be summed accross lenghths, everything else will be averaged. If this is not what you want, you'll have to deal with these slots by hand.")
	args <- list(...)
	args$model <- model
	args$FUN <- "l2a"	
	dnms <- dimnames(index(object))
	
	# trick to get the relevant slots
	lst <- qapply(object, function(x) list(x))
	lst <- lapply(lst, "[[", 1)	
	nms <- qapply(object, quant)

	# first the sums	
	isLen <- lapply(nms, "==", "len")
	# aggregated are not lenghts
	isLen[c("effort")][] <- FALSE
	# whatelse 
	isLen[c("index", "index.var", "catch.wt", "sel.pattern", "index.q")][] <- FALSE
	lst[unlist(isLen)] <- lapply(lst[unlist(isLen)], function(x){
		args$object <- x
		do.call("l2a", args)
	})

	# now the means	
	isLen <- lapply(nms, "==", "len")
	# aggregated are not lenghts
	isLen[c("effort")][] <- FALSE
	# whatelse
	isLen[c("catch.n", "catch.wt")][] <- FALSE
	args$stat <- "mean"
	lst[unlist(isLen)] <- lapply(lst[unlist(isLen)], function(x){
		args$object <- x
		do.call("l2a", args)
	})

	# now the weighted means	
	isLen <- lapply(nms, "==", "len")
	# aggregated are not lenghts
	isLen[c("effort")][] <- FALSE
	# whatelse
	isLen[c("catch.n", "index", "index.var", "sel.pattern", "index.q")][] <- FALSE
	args$stat <- "mean"
	lst[unlist(isLen)] <- lapply(lst[unlist(isLen)], function(x){
		args$object <- x
		args$weights <- catch.n(object)
		do.call("l2a", args)
	})

	# create the index object
	lst <- lapply(lst, "quant<-", "age")
	lst$name <- object@name
	lst$desc <- object@desc
	idx <- do.call("FLIndex", lst)
	idx@range[c("startf","endf")] <- object@range[c("startf","endf")]
	idx
})




