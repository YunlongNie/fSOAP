

seed_set= 10
library(fda);library(dplyr)
tempdat = daily$tempav
all_matrix=tempdat;all_matrix = sweep(all_matrix,1,rowMeans(all_matrix))
timepts=0:364;
norder=4 ## cubic B-spline
nbasis=norder+length(timepts)-2;
spline.basis=create.bspline.basis(rangeval=c(0,364),nbasis=10,norder=4)
D2Lfd <- int2Lfd(m=2)
D2fdPar <- fdPar(spline.basis, D2Lfd, 1e1)
all.fd =  Data2fd(y=all_matrix, argvals=timepts)
fpca = pca.fd(all.fd,2,D2fdPar)
pc1 = fpca$harmonics[1]
pc2 = fpca$harmonics[2]

set.seed(seed_set)
nsim = 300
nobs = sample(x=1:5,size=nsim, replace=TRUE)
s=(rnorm(n=nsim,sd=30))
# s=(rgamma(n=nsim,shape=1,rate=0.03))
# s = s-mean(s)
coefsy= matrix(rep(s, each =pc1$basis$nbasis), pc1$basis$nbasis)*(matrix(rep(coef(pc1), length(s)), length(s),byrow=TRUE)%>%t)
s2=(rnorm(n=nsim,sd=10))
# s2=(rgamma(n=nsim,shape=1,rate=0.1))
# s2 = s2-mean(s2)
coefsy2= matrix(rep(s2, each =pc2$basis$nbasis), pc2$basis$nbasis)*(matrix(rep(coef(pc2), length(s)), length(s),byrow=TRUE)%>%t)
yfds = fd(coefsy+coefsy2, pc2$basis)

timepoints = lapply(1:length(nobs), function(x) runif(nobs[x],0,364)%>%sort)

true_y = lapply(1:length(nobs), function(x) eval.fd(timepoints[[x]], yfds[x])%>%as.numeric)

sigma = 0.05
observed = lapply(1:length(nobs), function(x) true_y[[x]]+rnorm(length(true_y[[x]]), sd=sigma ))

simulate_data=list()
simulate_data$timepoints=timepoints
simulate_data$true_y = true_y
simulate_data$observed = observed
simulate_data$yfds = yfds
simulate_data$pc1 = pc1
simulate_data$pc2 = pc2


save(simulate_data, file="/Users/joha/Dropbox/Rpackages/SOAP_source/simulate_data.rda")

load("/Users/joha/Dropbox/Rpackages/SOAP_source/simulate_data.rda")


source('/project/6006512/ynie01/predict_last_obs_EM/funs_SOAP.R')


spline_basis = spline.basis
spline.basis=create.bspline.basis(rangeval=c(0,364),nbasis=15,norder=4)

observed%>%do.call(c,.)%>%mean
beta0 = rep(observed%>%do.call(c,.)%>%mean,spline.basis$nbasis)

FEC1  = firstFEC(ylist=observed, tlist=timepoints,spline_basis=spline.basis, gamma=1e6,threshold=1e-5)

plot(FEC1$fec,col=2)

previous_beta = list()
previous_beta[[1]] = FEC1$beta


pc2s = third_FPC_conditional(ylist=observed, tlist=timepoints,beta_previous=previous_beta,spline_basis=spline.basis, gamma=1e6,,threshold=1e-5)
 plot(pc2s$pc_fit)
 lines(pc2,col=2)

previous_beta[[2]] = pc2s$beta

pc3s = third_FPC_conditional( ylist=observed, tlist=timepoints, pc_index=4, gamma=1e6,betalist =previous_beta,threshold=1e-4)
previous_beta[[3]] = pc3s$beta

 # pc4s = third_FPC_conditional(runif(spline_basis$nbasis,-2,2), observed=observed, timepoints=timepoints, pc_index=5, gamma=1e4,betalist =previous_beta,threshold=1e-6)
 # previous_beta[[4]] = pc4s$beta

# compute AIC


betalist0=previous_beta
betalist=list()
AIC = c()
sig=c()
for (i in 1:length(previous_beta)){

betalist[1:i] = betalist0[1:i]

ylist = observed
tlist = timepoints
cfits = lapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	xmat = lapply(1:length(betalist), function(i){
		pc_fit  = fd(betalist[[i]], spline_basis)
		eval.fd(timei, pc_fit)%>%as.numeric
	})%>%do.call(cbind,.)
	index_pc = min(length(ylist[[x]])-1, length(betalist)) # to chose AIC we use -1 
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

residuals = lapply(1:length(observed), function(x){
	timei = timepoints[[x]]
	resid = observed[[x]] - as.numeric(eval.fd(timei,yfitsfd[x]))
	if (length(timei)<=-2) {
		return(NA)
	} else {
	return(resid^2)
	}
	
})

	
observed2 =observed
observed2[which(sapply(observed,length)<=1)]=NULL
timepoints2 = timepoints
timepoints2[which(sapply(observed,length)<=1)]=NULL

N = sapply(observed2,length)%>%sum
n = length(observed2)
# print(N)
# print(n)
sigmahat = mean(residuals%>%do.call(c,.),na.rm=TRUE)
sig = c(sig,sigmahat)
AIC  = c(AIC,N*log(sigmahat) + N  + 2*n*i)
}
print(AIC)
print(sig)
(k_select =which.min(AIC))


## predict 

# previous_beta = list()
# previous_beta[[1]] = pc1s$beta
# previous_beta[[2]] = pc2s$beta
# previous_beta[[3]] = pc3s$beta
# k_select=2
previous_beta0 = previous_beta[1:k_select]

res = pred_SOAP(previous_beta0,observed_test, timepoints, spline_basis,nminus=1)

error_soap1 = sapply(1:length(observed_test), function(x){
	inprod(yfds_test[x] - (res$yfd_fit)[x],yfds_test[x] - (res$yfd_fit)[x])
})
res = pred_SOAP(previous_beta0,observed_test, timepoints, spline_basis,nminus=2)

error_soap2 = sapply(1:length(observed_test), function(x){
	inprod(yfds_test[x] - (res$yfd_fit)[x],yfds_test[x] - (res$yfd_fit)[x])
})

print(summary(error_soap1[sapply(timepoints,length)>1]))
print(summary(error_soap2[sapply(timepoints,length)>1]))
print(summary(error_pace[sapply(timepoints,length)>1]))
dir.create('normal/')
save(error_soap1,error_soap2,error_pace,seed_set,res_pace, AIC,timepoints, pc1s,pc2s,pc3s, file=sprintf("normal/res%s_normal.Rdata", seed_set))
}

t2=Sys.time()
t2-t1
