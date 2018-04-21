#' This function computes the next FEC conditional on previous fitted FEC using the SOAP method
#' @param ylist a list of observed longtitudinal data; each element in the list is a numeric vector for one subject's observed data
#' @param tlist a list of observed time points
#' @param beta_previous a list of previously estimated FECs' coefficients
#' @param spline_basis the B-spline basis for the fitted FEC
#' @param gamma a positive number, the smoothing parameter 
#' @export FECconditional 
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

FECconditional  = function(ylist, tlist, beta_previous, spline_basis, threshold=1e-4,gamma=0){
beta3 = runif(spline_basis$nbasis,-2,2)
# if(missing(pc_index)) 
pc_index= length(beta_previous)+2
thresh = 1
it = 1
R = inprod(spline_basis, spline_basis,2,2)
E = inprod(spline_basis, spline_basis,0,0)

pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric

pc_fits_previous = lapply(1:length(beta_previous), function(x){
pc_fit1 = fd(beta_previous[[x]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
pc_fit1
})
value = -1e14

ylist2 =ylist
ylist2[which(sapply(ylist,length)<=pc_index)]=NULL
tlist2 = tlist
tlist2[which(sapply(ylist,length)<=pc_index)]=NULL

beta3_before = rep(0, length(beta_previous[[1]]))
while (thresh >  threshold|it<10){
	
	beta3_before = beta3
	value_before = value

alpha_fit = function(subj_index){
	(timei =tlist2[[subj_index]])
	
	xmat_previous = lapply(pc_fits_previous, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit)
	lm(ylist2[[subj_index]]~xmat_previous+xmati+0)%>%coef%>%as.numeric
}

sfit= lapply(1:length(ylist2),alpha_fit)%>%do.call(rbind,.)

residual_fit = function(subj_index){
	(timei =tlist2[[subj_index]])
	xmat_previous = lapply(pc_fits_previous, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit)
	lm(ylist2[[subj_index]]~xmat_previous+xmati+0)%>%residuals
}

rfit= lapply(1:length(ylist2),residual_fit)%>%do.call(c,.)

value = as.numeric(mean(rfit^2)+gamma*inprod(pc_fit,pc_fit,2,2))
if(abs(value_before - value/value_before) < threshold) break
N = sapply(ylist2,length)%>%sum

yem = lapply(1:length(ylist2), function(i){
(timei =tlist2[[i]])
yfits_previous = lapply(1:length(pc_fits_previous), function(x){
	sfit[i,x]*eval.fd(timei,pc_fits_previous[[x]])%>%as.numeric
	})%>%do.call(rbind,.)%>%colSums

(ylist2[[i]] - yfits_previous)/sqrt(N)
})%>%do.call(c,.)


xalpha = lapply(1:length(ylist2),function(subj_index){

	(timei =tlist2[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit[subj_index,ncol(sfit)])/sqrt(N)

})%>%do.call(rbind,.)


A = xalpha
qmat = 2*(t(A)%*%A+gamma*R)
pmat = as.numeric(-2*t(yem)%*%A)
betamat=  do.call(rbind,beta_previous)
cmat = betamat%*%E
beta3  = qp(qmat, pmat, cmat, d=rep(0,length(beta_previous)))
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric

thresh =  max(abs(beta3_before-beta3))
loss_value= mean((yem - A%*%beta3)^2)+as.numeric(gamma*t(beta3)%*%R%*%beta3)
#print(value)
it = it+1
# if(it%%10==0) {print(it);print(thresh)}
}
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric
coveraged = (thresh<threshold)
options = list()
options$coveraged=coveraged
options$gamma = gamma
options$thresh=thresh
options$iterations= it
options$basis = spline_basis
options$beta_previous = beta_previous
return(list(beta=beta3,fec = pc_fit, fec_score = sfit, options=options))
}
