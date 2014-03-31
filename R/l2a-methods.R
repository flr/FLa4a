#' @title Method to convert age to length 
#' @name l2a 
#' @rdname l2a 
#' @aliases l2a l2a-methods l2a,FLQuant,a4aGr-method
#' @param object an \code{FLQuant}, or \code{FLStockLen} object. 
#' @param model a \code{a4aGr} object
#' @param stat the aggregation statistic, must be \"mean\" or \"sum\". Only used if object is an \code{FLQuant}.
#' @param plusgroup the plusgroup of the stock. Only used if the object is a \code{FLStockLen}.
#' @return an age based \class{FLQuant}, \class{FLStock}
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
	#mod <- FLModelSim(model=grInvMod(model), params=params(model), vcov=vcov(model), distr=distr(model))
	#age <- floor(predict(mod, len=len))
	age <- floor(predict(model, len=len))
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


    # Make the stock piece by piece to avoid memory problems
    catch.n <- l2a(catch.n(object), model, stat="sum",...)
    stk <- FLStock(catch.n=catch.n) 
    cat("Processing sum slots\n")
    sum_slots_names <- c("discards.n","landings.n","stock.n")
    for(slot_counter in sum_slots_names){
        slot(stk,slot_counter) <- l2a(slot(object,slot_counter), model, stat="sum",...)
        gc()
    }
    cat("Processing mean slots\n")
    mean_slots_names <- c("catch.wt","discards.wt","landings.wt","stock.wt","m","mat","harvest.spwn","m.spwn","harvest")
    for(slot_counter in mean_slots_names){
        slot(stk,slot_counter) <- l2a(slot(object,slot_counter), model, stat="mean",...)
        gc()
    }

    cat("Washing up\n")
	stk@name <- object@name
	stk@desc <- object@desc
	units(harvest(stk)) <- units(object@harvest)

    # Calc landings, discards, catch and stock
    # Correct weights as we are using average weights and totals in
    # age based object should match the length based object
    landings.wt(stk) <- sweep(landings.wt(stk),2:6, computeLandings(object) / computeLandings(stk), "*")
    landings(stk) <- computeLandings(stk)
    discards.wt(stk) <- sweep(discards.wt(stk),2:6, computeDiscards(object) / computeDiscards(stk), "*")
    discards(stk) <- computeDiscards(stk)
    catch.wt(stk) <- sweep(catch.wt(stk),2:6, computeCatch(object) / computeCatch(stk), "*")
    catch(stk) <- computeCatch(stk)
    stock.wt(stk) <- sweep(stock.wt(stk),2:6, computeStock(object) / computeStock(stk), "*")
    stock(stk) <- computeStock(stk)


	# set the plus group on the first non continuous age
    if(!missing(plusgroup)){
        stk <- setPlusGroup(stk, plusgroup, na.rm=T)
    }
    return(stk)
})

#' @rdname l2a 
#' @aliases l2a,FLIndex,a4aGr-method
setMethod("l2a", c("FLIndex", "a4aGr"), function(object, model, ...){
	warning("Catch in numbers will be summed accross lenghths, everything else will be averaged. If this is not what you want, you'll have to deal with these slots by hand.")
	args <- list(...)

    catch.n <- l2a(catch.n(object), model, stat="sum",...)
    idx <- FLIndex(catch.n=catch.n) 
    mean_slots_names <- c("index","catch.wt","index.var","sel.pattern","index.q")
    cat("Processing mean slots\n")
    for(slot_counter in mean_slots_names){
        slot(idx,slot_counter) <- l2a(slot(object,slot_counter), model, stat="mean",...)
        gc()
    }
    # Weighted means?

	idx@name <- object@name
	idx@desc <- object@desc
	idx@range[c("startf","endf")] <- object@range[c("startf","endf")]
	return(idx)
})




