
#' @title S4 class \code{a4aFitCatchDiagn}
#' @description The \code{a4aFitCatchDiagn} class extends \code{FLQuants} to store information to run diagnostics on aggregated catch estimated  by the a4a stock assessment fit.
#' @docType class
#' @name a4aFitCatchDiagn-class
#' @rdname a4aFitCatchDiagn-class
#' @aliases a4aFitCatchDiagn-class

setClass("a4aFitCatchDiagn", contain="FLQuants")
# VALIDITY CHECK NAMES


#' @rdname a4aFitCatchDiagn-class
#' @aliases a4aFitCatchDiagn a4aFitCatchDiagn-methods
#' @template bothargs
#' @param stock \code{FLStock} object used to fit the model
#' @param indices \code{FLIndices} object used to fit the model
#' @examples
#' data(ple4)
#' data(ple4.index)
#' fit <- sca(ple4, ple4.index)
#' flqs <- computeCatchDiagnostics(fit, ple4)
setGeneric("computeCatchDiagnostics", function(object, ...) standardGeneric("computeCatchDiagnostics"))

#' @rdname a4aFit-class
#' @aliases clock clock-methods
setMethod("computeCatchDiagnostics", signature(object="a4aFit"), function(object, stock, ...) {
	args <- list(...)
	it <- 500
	pred <- predict(object)
	if(dims(object)$iter==1) fits <- simulate(object, it) else it <- dims(object)$iter
	stk_new <- stock + object
	cth_est <- catch(stk_new)
	cth_obs <- catch(stock)
	# sim with observation error
	flqoe <- quantSums(rlnorm(it, log(catch.n(object)), pred$vmodel$catch)*catch.wt(stock))
	# sim with estimation error
	flqee <- catch(stock + fits)
	# both
	flqe <- quantSums(rlnorm(it, log(catch.n(object)), sqrt(pred$vmodel$catch^2 + iterVars(log(catch.n(fits)))))*catch.wt(stock))
	# standardized residuals
	resstd <- stdlogres(cth_obs, cth_est)
	# pearson residuals
	sdlog <- sqrt(iterVars(log(flqoe)))
	resprs <- stdlogres(cth_obs, cth_est, sdlog=sdlog)
	# deviances
	devs <- log(cth_obs/cth_est)
	# out
	flqs <- FLQuants(list(obs = cth_obs, est = cth_est, oe=flqoe, ee=flqee, oee=flqe, resstd=resstd, resprs=resprs, resraw=devs))
	new("a4aFitCatchDiagn", flqs)
  }
)

#' @title Plot of aggregated catch standardized log residuals
#' @name plot of catch residuals
#' @docType methods
#' @rdname plot-catch
#' @aliases plot,a4aFitCatchDiagn,missing-method
#' @description Method to produce scatterplots of aggregated catch residuals
#' @param x an \code{a4aFit} object with the model fit
#' @param y the \code{FLStock} object used to fit the model
#' @param ... additional argument list that might never be used
#' @return a \code{plot} with stardardized log residuals
#' @examples
#' data(ple4)
#' data(ple4.index)
#' fit <- sca(ple4, ple4.index)
#' flqs <- computeCatchDiagnostics(fit, ple4)
#' plot(flqs)

setMethod("plot", c("a4aFitCatchDiagn", "missing"), function(x, y=missing, ...){

	#----------------------------------------------------------------
	# absolute catches
	#----------------------------------------------------------------
	
	# build datasets catch
	probs <- c(0.10, 0.50, 0.90)
	d1 <- as.data.frame(quantile(x$oe, probs))%>%mutate(test="Observation error")%>%
	  pivot_wider(names_from = iter, values_from = data)%>%rename(median=`50%`)
	d2 <- as.data.frame(quantile(x$ee, probs))%>%mutate(test="Estimation error")%>%
	  pivot_wider(names_from = iter, values_from = data)%>%rename(median=`50%`)
	d3 <- as.data.frame(quantile(x$oee, probs))%>%mutate(test="Prediction error")%>%
	  pivot_wider(names_from = iter, values_from = data)%>%rename(median=`50%`)
	obs <- as.data.frame(x$obs)%>%select(year,data)%>%rename(observed=data)
	plot_data<-rbind(d1,d2,d3)%>%left_join(obs,by="year")%>%
	  select(year,test,`10%`,`90%`,median,observed)%>%
	  pivot_longer(!c("year","test","10%","90%"),names_to="CI 80%",values_to="Catch (t)")
	
	# plot catch
	plot_catch<-ggplot(plot_data,aes(x=year,y=`Catch (t)`,linetype=`CI 80%`))+
	  geom_line()+
	  theme_bw()+
	  geom_ribbon(aes(ymin=`10%`, ymax=`90%`), fill = "grey",alpha = 0.3)+
	  facet_wrap(~factor(test,levels=c("Observation error","Prediction error","Estimation error")),ncol=1)+
	  labs(x = "Year")+
	  theme(text=element_text(size=18))
	
	#----------------------------------------------------------------
	# residuals
	#----------------------------------------------------------------
	
	
	# build datasets residuals
	res1 <- as.data.frame(x$resstd)%>%mutate(test="Standardized residuals")%>%
	  pivot_wider(names_from = iter, values_from = data)
	res2 <- as.data.frame(x$resprs)%>%mutate(test="Pearson residuals")%>%
	  pivot_wider(names_from = iter, values_from = data)
	res3 <- as.data.frame(x$resraw, probs)%>%mutate(test="Raw residuals (deviances)")%>%
	  pivot_wider(names_from = iter, values_from = data)
	plot_data_res<-rbind(res1,res2,res3)%>%
	  select(year,test,`1`)
	
	# plot residuals
	plot_res<-ggplot(plot_data_res, aes(x = year, y = `1`)) +
	  geom_segment(aes(x = year, xend = year, y = 0, yend = `1`), color = "black") +
	  theme_bw()+geom_hline(yintercept=0, linetype='dashed')+
	  geom_point(size = 2) +facet_wrap(~test,ncol=1)+
	  labs(y = "Residuals",x = "Year")+
	  theme(text=element_text(size=18))
	
	#join plot
	final_plot<-ggarrange(plot_catch,plot_res,common.legend = T, ncol=2,legend="bottom")
	annotate_figure(final_plot,top = text_grob("Aggregated catch",size=18))

})

