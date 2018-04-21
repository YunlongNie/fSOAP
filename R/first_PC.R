#' This function computes the first FEC using the SOAP method
#' @param ylist a list of observed longtitudinal data; each element in the list is a numeric vector for one subject's observed data
#' @param tlist a list of observed time points
#' @param spline_basis the B-spline basis for the fitted FEC
#' @param gamma a positive number, the smoothing parameter 
#' @export firstFEC
#' @import fda
#' @import dplyr
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


firstFEC  = function(ylist,tlist,spline_basis, gamma=0,threshold=1e-5,maxit=100){

	initial_beta = rep(ylist%>%do.call(c,.)%>%mean,spline_basis$nbasis)

	thresh = 1
	it = 1
	initial_beta_before = rep(0, length(initial_beta))
	value = -1e14
	pc_fit = fd(initial_beta, spline_basis)
	pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
	initial_beta = coef(pc_fit)%>%as.numeric
	R = inprod(spline_basis, spline_basis,2,2)

while (thresh >  threshold|it<10){
	if (it >maxit) {cat('Reach maximum iteration\b');break }
	initial_beta_before = initial_beta
	value_before = value

	alpha_fit = function(subj_index){
		(timei =tlist[[subj_index]])
		xmati = eval.fd(timei, pc_fit)
		lm(ylist[[subj_index]]~0+xmati)%>%coef%>%as.numeric
	}

	sfit= sapply(1:length(ylist),alpha_fit)

	residual_fit = function(subj_index){
		(timei =tlist[[subj_index]])
		xmati = eval.fd(timei, pc_fit)
		(lm(ylist[[subj_index]]~xmati+0)%>%residuals)^2%>%mean
	}

	rfit= lapply(1:length(ylist),residual_fit)%>%do.call(c,.)
	value = mean(rfit^2)+gamma*inprod(pc_fit,pc_fit,2,2)
	if(abs(value_before - value/value_before) < threshold) break

	yem = do.call(c,ylist)

	xalpha = lapply(1:length(ylist),function(subj_index){

		(timei =tlist[[subj_index]])
		xmati = eval.basis(timei, spline_basis)
		xmati*sfit[subj_index]

	})%>%do.call(rbind,.)
	initial_beta = solve(t(xalpha%>%as.matrix)%*%(as.matrix(xalpha)) + gamma*R, t(xalpha%>%as.matrix)%*%as.matrix(yem))
	pc_fit = fd(initial_beta, spline_basis)
	pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
	initial_beta = coef(pc_fit)%>%as.numeric
	thresh =  max(abs(initial_beta_before-initial_beta))
	it = it+1
}

pc_fit = fd(initial_beta, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
initial_beta = coef(pc_fit)%>%as.numeric
coveraged = (thresh<threshold)
options = list()
options$coveraged=coveraged
options$gamma = gamma
options$thresh=thresh
options$iterations= it
options$basis = spline_basis

return(list(beta=initial_beta, fec = pc_fit,fec_score = sfit, options=options))
}


