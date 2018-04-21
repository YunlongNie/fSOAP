

library(fda)
data(simulate_data)
observed = simulate_data$observed
timepoints =simulate_data$timepoints
spline.basis=create.bspline.basis(rangeval=c(0,364),nbasis=15,norder=4)
beta0 = rep(observed%>%do.call(c,.)%>%mean,spline.basis$nbasis)
FEC1  = firstFEC(ylist=observed, tlist=timepoints,spline_basis=spline.basis, gamma=1e6,threshold=1e-5)
previous_beta = list()
previous_beta[[1]] = FEC1$beta
FEC2 = FECconditional(ylist=observed, tlist=timepoints,beta_previous=previous_beta,spline_basis=spline.basis, gamma=1e6,threshold=1e-5)
previous_beta[[2]] = FEC2$beta
FEC3 = FECconditional(ylist=observed, tlist=timepoints,beta_previous=previous_beta,spline_basis=spline.basis, gamma=1e6,threshold=1e-5)
previous_beta[[3]] = FEC3$beta

betas = do.call(cbind, previous_beta)
colnames(betas)  =c("FEC1","FEC2","FEC3")
fecs = fd(betas, spline.basis)

library(ggplot2)
fdagg(fecs)

predict_y = predict_SOAP(previous_beta,ylist=observed, tlist=timepoints, spline_basis=spline.basis,nminus=2)

i=6
plot(predict_y$predict[i],ylim=range(observed[[i]]))
lines(simulate_data$yfds[i])
points(x=timepoints[[i]],y=observed[[i]])
