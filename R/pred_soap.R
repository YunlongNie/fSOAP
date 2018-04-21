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
#'library(fda)
#'library(cFuSIM)
#'data(bike_cFuSIM)
#'timepts = bike$timepts
#'norder=4 ## cubic B-spline
#'nbasis=norder+length(timepts)-2; 
#'spline_basis=create.bspline.basis(rangeval=c(1,24),nbasis#'norder,timepts)
#'wull = bike$temp
#'xfds=  Data2fd(y=wull%>%t, argvals=bike$timepts)
#'y = bike$y
#'train_sample = 1:length(y)
#'y = y[train_sample]
#'xfd = xfds[train_sample]
#'res_c = cFuSIM_index(y, xfd, spline_basis)
#'beta_fd = fd(res_c$coefBeta, res_c$basisBeta)
#'plot(beta_fd,ylab="index function", xlab='time')
#'fdagg(beta_fd)
#'score_fit = (res_c$score_fit)
#' pred_y = localpoly.reg(score_fit, y, degree.pol = 1, kernel.type = "gaussian",bandwidth = "CV",deriv=0,points#'score_fit)
#'plot(x=score_fit, y=y)
#'lines(pred_y$predicted[order(score_fit)],x=score_fit[order#'score_fit)],col=4)
#' }
predict_SOAP = function(betalist,ylist, tlist,spline_basis,nminus=2){

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
list(predicted=  yfitsfd,sigmahat = sigmahat)

}