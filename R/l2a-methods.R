#' @title Method to convert age to length 
#' @name l2a 
#' @rdname l2a 
#' @aliases l2a l2a-methods l2a,FLQuant,a4aGr-method
#' @param object an \code{FLQuant}, or \code{FLStockLen} object. 
#' @param model a \code{a4aGr} object
#' @param stat the aggregation statistic, must be \"mean\" or \"sum\". Only used if object is an \code{FLQuant}.
#' @param min_age minimum age of the FLQuant (all younger ages are set to this age). Only used if the object is a \code{FLQuant}.
#' @param max_age maximum age of the FLQuant (all older ages are set to this age). Only used if the object is a \code{FLQuant} or an \code{FLIndex}.
#' @param plusgroup the plusgroup of the stock. Only used if the object is a \code{FLStockLen}.
#' @return an age based \code{FLQuant}, \code{FLStock}
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
setMethod("l2a", c("FLQuant", "a4aGr"),
	function(object, model, stat="sum", max_age=NA, min_age=0,
		weights=FLQuant(1, dimnames=dimnames(object))) {
    
	# constants
	#cat("Converting lengths to ages ...\n")
	dnms <- dimnames(object)
	if(!all.equal(dnms, dimnames(weights))) stop("Weights must have the same dimensions as the data.")
	len <- as.numeric(dnms[[1]])
	mit <- niters(model)
	qit <- length(dnms[[6]]) 
	if(mit>1 & qit>1 & mit!=qit) stop("Can not operate with diferent iterations in each object.")
	age <- floor(predict(model, len=len))
	# aggregates all lengths above linf
	age <- apply(age, 2, function(x){
		x[which.max(x[!is.infinite(x)]):length(x)] <- max(x[!is.infinite(x)], na.rm=T)
		x
	})
    # set min and max ages
    age[age < min_age] <- min_age
    if (!is.na(max_age)){
        age[age > max_age] <- max_age
    }
	ages <- seq(min(age, na.rm=TRUE),max(age, na.rm=TRUE))
	# FUN to use
    if(stat=="sum"){
        FUN <- colSums
    }
    if(stat=="mean"){
        FUN <- colMeans
    }
    res <- array(NA, dim=c(length(ages), dim(object)[-c(1,6)], max(qit,mit)))
    if (mit >= qit){
        length_array <- array(rep(c(object), dim(age)[2]), dim=c(dim(object)[1:5], mit)) # propagate obj from 1 iter to n iters
        lak_array <- aperm(array(rep(c(age), times=prod(dim(object)[2:5])), dim=c(dim(object)[1], dim(age)[2],dim(object)[2:5])),c(1,3,4,5,6,2))
    }
    if (mit < qit){
        length_array <- array(c(object), dim=dim(object)) # make object with no dimnames
        lak_array <- aperm(array(rep(c(age), times=prod(dim(object)[2:6])), dim=c(dim(object)[1], dim(object)[6],dim(object)[2:5])),c(1,3,4,5,6,2))
    }
    s <- array(NA, dim=c(dim(object)[1:5], max(qit,mit)))
    unique_ages <- unique(c(age))
    #for(i in seq(length(unique_ages))) {
    for(i in unique_ages) {
        s[,,,,,] <- NA 
        idx <- lak_array == i
        s[idx] <- length_array[idx]
        res[i-min(ages)+1,,,,,] <- do.call(FUN,list(s,na.rm=TRUE))
    }
	
    # fill NAs with 0
    res[is.na(res)] <- 0
    out <-  FLQuant(res, dimnames=c(list(age=ages), dimnames(object)[-c(1, 6)]), units=units(object))
	return(out)
})


#' @rdname l2a 
#' @aliases l2a,FLStockLen,a4aGr-method
setMethod("l2a", c("FLStockLen", "a4aGr"), function(object, model, plusgroup=NA, ...){
	warning("Individual weights, M and maturity will be (weighted) averaged accross lengths, harvest is not computed and everything else will be summed.\n If this is not what you want, you'll have to deal with these slots by hand.")

    # Make the stock piece by piece to avoid memory problems
    #cat("Processing sum slots\n")
    catch.n <- l2a(catch.n(object), model, stat="sum", max_age=plusgroup,...)
    stk <- FLStock(catch.n=catch.n)
    qsize <- prod(dim(catch.n))
    # check which slots need slicing
    sum_slots_names <- c("discards.n","landings.n","stock.n")
	# if there's no discards landings=catches
	sum_slots_names <- sum_slots_names[c(rep(sum(is.na(discards.n(object)))!=qsize,2),sum(is.na(stock.n(object)))!=qsize)]
	if(!is.empty(sum_slots_names)){
		for(slot_counter in sum_slots_names){
		    slot(stk,slot_counter) <- l2a(slot(object,slot_counter), model, stat="sum", max_age=plusgroup,...)
		    gc()
		}
	}

    #cat("Processing mean slots\n")
    #mean_slots_names <- c("catch.wt","discards.wt","landings.wt","stock.wt","m","mat","harvest.spwn","m.spwn")
    mean_slots_names <- c("m","mat","harvest.spwn","m.spwn")
    for(slot_counter in mean_slots_names){
        slot(stk,slot_counter) <- l2a(slot(object,slot_counter), model, stat="mean", max_age=plusgroup,...)
        gc()
    }

    #cat("Processing weighted mean slots\n")
    weighted_means_slots_names <- c("catch","discards","landings","stock")
    for(slot_counter in weighted_means_slots_names){
        total_slice <- l2a(slot(object,paste(slot_counter,".wt",sep="")) * slot(object,paste(slot_counter,".n",sep="")), model, stat="mean", max_age=plusgroup,...)
        slot(stk,paste(slot_counter,".wt",sep="")) <- total_slice / slot(stk,paste(slot_counter,".n",sep=""))
        gc()
    }

    #cat("Washing up\n")
	stk@name <- object@name
	stk@desc <- object@desc
	units(harvest(stk)) <- units(object@harvest)

    # Calc landings, discards, catch and stock
    landings(stk) <- computeLandings(stk)
    discards(stk) <- computeDiscards(stk)
    catch(stk) <- computeCatch(stk)
    stock(stk) <- computeStock(stk)

	# set the plus group on the first non continuous age
    if(!is.na(plusgroup)){
        stk <- setPlusGroup(stk, plusgroup, na.rm=T)
    }
    return(stk)
})

#' @rdname l2a 
#' @aliases l2a,FLIndex,a4aGr-method
setMethod("l2a", c("FLIndex", "a4aGr"), function(object, model, ...){
	warning("Catch in numbers will be summed accross lengths, everything else will be averaged. If this is not what you want, you'll have to deal with these slots by hand.")
	args <- list(...)

    catch.n <- l2a(catch.n(object), model, stat="sum",...)
    idx <- FLIndex(catch.n=catch.n) 
    mean_slots_names <- c("index","index.var","sel.pattern","index.q")
    #cat("Processing mean slots\n")
    for(slot_counter in mean_slots_names){
        slot(idx,slot_counter) <- l2a(slot(object,slot_counter), model, stat="mean",...)
        gc()
    }

    #cat("Processing weighted mean slots\n")
    weighted_means_slots_names <- c("catch")
    for(slot_counter in weighted_means_slots_names){
        total_slice <- l2a(slot(object,paste(slot_counter,".wt",sep="")) * slot(object,paste(slot_counter,".n",sep="")), model, stat="mean",...)
        slot(idx,paste(slot_counter,".wt",sep="")) <- total_slice / slot(idx,paste(slot_counter,".n",sep=""))
        gc()
    }


	idx@name <- object@name
	idx@desc <- object@desc
	idx@range[c("startf","endf")] <- object@range[c("startf","endf")]
	return(idx)
})




