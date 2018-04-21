#' This function computes the predicted trajectories based on etimated FECs and sparse observatios using the SOAP method 
#' @param betalist a list of all the FECs' coefficients
#' @param ylist a list of observed longtitudinal data; each element in the list is a numeric vector for one subject's observed data
#' @param tlist a list of observed time points
#' @param spline_basis the B-spline basis for the fitted FEC
#' @param nminus the total number of FECs in the betalist -  nminus is the number of FECs used for predicting the individual's trajectory
#' @export predict_SOAP
#' @import fda
#' @import dplyr
#' @import lsei
#' @examples
#' \dontrun{
#' 
#' library(fda)
#' data(simulate_data)
#' observed = simulate_data$observed
#' timepoints =simulate_data$timepoints
#' spline.basis=create.bspline.basis(rangeval=c(0,364),nbasis=15,norder=4)
#' beta0 = rep(observed%>%do.call(c,.)%>%mean,spline.basis$nbasis)
#' FEC1  = firstFEC(ylist=observed, tlist=timepoints,spline_basis=spline.basis, gamma=1e6,threshold=1e-5)
#' previous_beta = list()
#' previous_beta[[1]] = FEC1$beta
#' FEC2 = FECconditional(ylist=observed, tlist=timepoints,beta_previous=previous_beta,spline_basis=spline.basis, gamma=1e6,threshold=1e-5)
#' previous_beta[[2]] = FEC2$beta
#' FEC3 = FECconditional(ylist=observed, tlist=timepoints,beta_previous=previous_beta,spline_basis=spline.basis, gamma=1e6,threshold=1e-5)
#' previous_beta[[3]] = FEC3$beta
#' 
#' betas = do.call(cbind, previous_beta)
#' colnames(betas)  =c("FEC1","FEC2","FEC3")
#' fecs = fd(betas, spline.basis)
#' 
#' library(ggplot2)
#' fdagg(fecs)
#' 
#' predict_y = predict_SOAP(previous_beta,ylist=observed, tlist=timepoints, spline_basis=spline.basis,nminus=2)
#' 
#' i=6
#' plot(predict_y$predict[i],ylim=range(observed[[i]]))
#' lines(simulate_data$yfds[i])
#' points(x=timepoints[[i]],y=observed[[i]])
#' }

predict_SOAP = function(betalist,ylist, tlist,spline_basis,testt, nminus=2){
if(missing(testt)) testt = NULL
cfits = lapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	xmat = lapply(1:length(betalist), function(i){
		pc_fit  = fd(betalist[[i]], spline_basis)
		eval.fd(timei, pc_fit)%>%as.numeric
	})%>%do.call(cbind,.)
	index_pc = min(length(ylist[[x]])-nminus, length(betalist))
	if(index_pc<=0) index_pc=1
	cfit = lm(ylist[[x]]~0+xmat[,1:index_pc])%>%coef%>%as.numeric
	if(length(cfit)<length(betalist)) cfit = c(cfit, rep(0,length(betalist)-length(cfit)))
	cfit
})%>%do.call(rbind,.)

yfits = lapply(1:nrow(cfits), function(i){
	cfit = cfits[i,]
	rowSums(mapply("*",cfit, betalist,SIMPLIFY=TRUE))
}
)%>%do.call(cbind,.)

yfitsfd = fd(yfits, spline_basis)


residuals = sapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	resid = ylist[[x]] - as.numeric(eval.fd(timei,yfitsfd[x]))
	if (all(resid==0)) {
		return(NA)
	} else {
	return(mean(resid^2))
	}
	
})
 sigmahat = mean(residuals[!is.na(residuals)])
 testd = list();testd$time = ttest
if (!is.null(testt)) {
	ytestt = eval.fd(testt, yfitsfd) 
	testd$y = ytestt
}
list(predicted=  yfitsfd,testd=testd,sigmahat = sigmahat )

}