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


