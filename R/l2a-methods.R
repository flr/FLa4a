
#' @title Method to convert length-based data to age-based
#' @description Method to convert length-based data to age-based
#' @details
#' A deterministic slicing method converts the length-based data to age-based data, using the supplied growth model (the \code{a4aGr} object).
#' Each length-based observation is allocated to a corresponding age, based on the growth model, and aggregated accordingly (either the sum or the mean).
#' There should be 1 or n iterations in both the object being sliced and the growth model.
#' This means that although the slicing method is deterministic, the length-based data is sliced by each iteration of the growth parameters, thereby propagating the uncertainty in the biological growth parameters (representing process uncertainty) through to the age-based data.
#' 
#' @name l2a 
#' @rdname l2a 
#' @aliases l2a l2a-methods
#' @param object an \code{FLQuant}, or \code{FLStockLen} object. 
#' @param model an \code{a4aGr} object
#' @param halfwidth the halfwidths of the length classes; a single numeric or vector the size of the number of length classes; not used if object is an \code{FLStockLen}
#' @param stat the aggregation statistic, which must be \code{mean} or \code{sum}; only used if object is an \code{FLQuant}.
#' @param max_age the maximum age in the returned \code{FLQuant}; all ages above this are set to \code{max_age}; only used if object is an \code{FLQuant}
#' @param plusgroup the plusgroup of the stock; only used if the object is an \code{FLStockLen}.
#' @template dots
#' @return an age based \code{FLQuant}, \code{FLStock}
#' @examples
#' # Southern hake
#' # Variance-covariance matrix for parameters
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(2310, 0.13, 0.84)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(-7.22,-6.28,0.08)
#' # Make the von Bertalanffy growth model
#' md <- ~linf*(1-exp(-k*(t-t0)))
#' imd <- ~t0-1/k*log(1-len/linf)
#' prs <- FLPar(linf=130, k=0.164, t0=-0.092, units=c("cm","yr-1","yr"))
#' vbObj <- a4aGr(grMod=md, grInvMod=imd, params=prs, vcov=mm, distr="norm")
#' # Make a triangle copula for simulating process error
#' linf <- list(a=104.5, b=155.5, c=130) 
#' k <- list(a=0.132, b=0.196, c=0.164)
#' t0 <- list(a=-0.184, b=0, c=-0.092)
#' tri_pars <- list(linf = linf, k = k, t0 = t0)
#' # Simulate 10 iterations from it
#' vbObj_tri <- mvrtriangle(10, vbObj, paramMargins=tri_pars)
#' data(southernHakeLen)
#' # Extract the catch numbers at length from stock object
#' cth <- catch.n(shake_len) 
#' # Slice the data using the unsimulated growth object 
#' # so the stock and the growth object have 1 iteration
#' cthA1 <- l2a(cth, vbObj)
#' # Slice with 1 iteration in stock and multiple in growth object
#' cthA2 <- l2a(cth, vbObj_tri)
#' # Result is age based catch with multiple iterations
#' # mod: iter=1, data: iter=n
#' cthA3 <- l2a(propagate(cth,10), vbObj)
#' # both with iter=n
#' cthA4 <- l2a(propagate(cth,10), vbObj_tri)
#' # converting a stock object
#' shake_age <- l2a(shake_len, vbObj)
#' shake_age <- l2a(shake_len, vbObj_tri)
#' shake_age <- l2a(propagate(shake_len, 10), vbObj)
#' shake_age <- l2a(propagate(shake_len, 10), vbObj_tri)
#' # converting a index object
#' index_pt_age <- l2a(index_pt_len, vbObj)
#' index_pt_age <- l2a(index_pt_len, mvrnorm(10, vbObj))
#' index_pt_age <- l2a(propagate(index_pt_len, 10), vbObj)
# l2a
setGeneric("l2a", function(object, model, ...) standardGeneric("l2a"))
#' @rdname l2a 
setMethod("l2a", c("FLQuant", "a4aGr"),
	function(object, model, halfwidth= c(diff(as.numeric(dimnames(object)[[1]])), tail(diff(as.numeric(dimnames(object)[[1]])),1))/2 , stat="sum", max_age=NA) {

# EJ: don't know what weights are doing ...
#	function(object, model, halfwidth= c(diff(as.numeric(dimnames(object)[[1]])), tail(diff(as.numeric(dimnames(object)[[1]])),1))/2 , stat="sum", weights=FLQuant(1, dimnames=dimnames(object)), max_age=NA) {
	# constants
	dnms <- dimnames(object)
#	if(!all.equal(dnms, dimnames(weights))) stop("Weights must have the same dimensions as the data.")
	len <- as.numeric(dnms[[1]]) + halfwidth 
	mit <- niters(model) # iters in the growth model
	qit <- length(dnms[[6]])  # iters in the FLQuant
	if(mit>1 & qit>1 & mit!=qit) stop("Iterations in model and object must be either 1 or X.")
	age <- floor(predict(model, len=len))
	# aggregates all lengths above linf
	age <- apply(age, 2, function(x){
		x[which.max(x[!is.infinite(x)]):length(x)] <- max(x[!is.infinite(x)], na.rm=T)
		x
	})
    # Truncate age range to save time and memory
    if (!(is.na(max_age))){
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
    if (mit > qit){
        # propagate obj from 1 iter to n iters (same as growth model)
        length_array <- array(rep(c(object), dim(age)[2]), dim=c(dim(object)[1:5], mit)) 
        lak_array <- aperm(array(rep(c(age), times=prod(dim(object)[2:5])), dim=c(dim(object)[1], dim(age)[2],dim(object)[2:5])),c(1,3,4,5,6,2))
    }
    if (mit < qit){
        length_array <- array(c(object), dim=dim(object)) # make object with no dimnames
        lak_array <- aperm(array(rep(c(age), times=prod(dim(object)[2:6])), dim=c(dim(object)[1], dim(object)[6],dim(object)[2:5])),c(1,3,4,5,6,2))
    }
    if (mit==qit){
        length_array <- array(c(object), dim=dim(object)) # make object with no dimnames
        lak_array <- aperm(array(rep(c(age), times=prod(dim(object)[2:5])), dim=c(dim(object)[1], dim(age)[2],dim(object)[2:5])),c(1,3,4,5,6,2))
    }
    s <- array(NA, dim=c(dim(object)[1:5], max(qit,mit)))
    unique_ages <- unique(c(age))
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
setMethod("l2a", c("FLStockLen", "a4aGr"), function(object, model, plusgroup=NA, ...){
	warning("Individual weights, M and maturity will be (weighted) averaged accross lengths,\n harvest is not computed and everything else will be summed.\n If this is not what you want, you'll have to deal with these slots by hand.")
    # Use the catch.n slot to build the resulting FLStock
    catch.n <- suppressWarnings(l2a(catch.n(object), model, halfwidth=halfwidth(object), stat="sum", max_age=plusgroup+1, ...))
    stk <- FLStock(catch.n=catch.n)
    # Abundance slots - sum these up
    sum_slots_names <- c("discards.n","landings.n","stock.n")
    for(slot_counter in sum_slots_names){
        # Only slice if there are some values in there
        if(!all(is.na(slot(object,slot_counter)))){
            slot(stk,slot_counter) <- suppressWarnings(l2a(slot(object,slot_counter), model, halfwidth=halfwidth(object), stat="sum", max_age=plusgroup+1, ...))
        }
    }
    # Weight slots - weighted means
    weighted_means_slots_names <- c("catch","discards","landings","stock")
    for(slot_counter in weighted_means_slots_names){
        total_quant <- slot(object,paste(slot_counter,".wt",sep="")) * slot(object,paste(slot_counter,".n",sep=""))
        # Only slice if not empty (either wt or n can be NA, e.g. stock.n before assessment)
        if(!all(is.na(total_quant))){
            total_slice <- suppressWarnings(l2a(total_quant, model, halfwidth=halfwidth(object), stat="sum", max_age=plusgroup+1, ...))
            slot(stk,paste(slot_counter,".wt",sep="")) <- total_slice / slot(stk,paste(slot_counter,".n",sep=""))
            # Replace any NAs with zeros
            slot(stk,paste(slot_counter,".wt",sep=""))[is.na(slot(stk,paste(slot_counter,".wt",sep="")))] <- 0
        }
    }

    # Other slots - mean
    mean_slots_names <- c("m","mat","harvest.spwn","m.spwn")
    for(slot_counter in mean_slots_names){
        if(!all(is.na(slot(object,slot_counter)))){
            slot(stk,slot_counter) <- suppressWarnings(l2a(slot(object,slot_counter), model, halfwidth=halfwidth(object), stat="mean", max_age=plusgroup+1, ...))
        }
    }

    # Check for ages < 0; report problem and trim
    if(range(stk)["min"] < 0){
        warning("Some ages are less than 0, indicating a mismatch between input data lengths\n and growth parameters (possibly t0)")
        warning("Trimming age range to a minimum of 0")
        stk <- trim(stk, age=0:range(stk)["max"])
    }
    # washing up
	stk@name <- object@name
	stk@desc <- object@desc
	units(harvest(stk)) <- units(object@harvest)
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
setMethod("l2a", c("FLIndex", "a4aGr"), function(object, model, ...){
    # Slots are treated differently
    # Sum: index, index.var, catch.n
    # Weighted sum: catch.wt, index.q - weighted by catch.n
    # Ignored: sel.pattern
	args <- list(...)

    # Start with index - most likely to not be empty
    index <- suppressWarnings(l2a(index(object), model, stat="sum", ...))
    idx <- FLIndex(index=index)


    sum_slots_names <- c("catch.n","index.var")
    for(slot_counter in sum_slots_names){
        # Only slice if there are some values in there
        if(all(!is.na(slot(object,slot_counter)))){
            slot(idx,slot_counter) <- suppressWarnings(l2a(slot(object,slot_counter), model, stat="sum", ...))
        }
    }

    weighted_means_slots_names <- c("catch.wt","sel.pattern","index.q")
    for(slot_counter in weighted_means_slots_names){
        total_quant <- slot(object,slot_counter) * slot(object,"catch.n")
        # Only slice if not empty (either wt or n can be NA, e.g. stock.n before assessment)
        if(all(!is.na(total_quant))){
            total_slice <- suppressWarnings(l2a(total_quant, model, stat="sum", ...))
            slot(idx, slot_counter) <- total_slice / slot(idx,"catch.n")
        }
    }

    # Check for ages < 0; report problem and trim
    if(range(idx)["min"] < 0){
        warning("Some ages are less than 0, indicating a mismatch between input data lengths\n and growth parameters (possibly t0)")
        warning("Trimming age range to a minimum of 0")
        idx <- trim(idx, age=0:range(idx)["max"])
    }

    # Washing up
	idx@name <- object@name
	idx@desc <- object@desc
	idx@distribution <- object@distribution
	idx@type <- object@type
	effort(idx)[] <- effort(object)
	idx@range[c("startf","endf")] <- object@range[c("startf","endf")]
	return(idx)
})




