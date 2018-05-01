# fSOAP
## install this packge by runing the following code: 

required_packages = c('dplyr','fda','ggplot2','lsei','devtools')
(missed_packages = setdiff(required_packages,rownames(installed.packages())))

if(length(missed_packages)){
	sapply(missed_packages, install.packages)
} 

devtools::install_github("YunlongNie/fSOAP")


# you can find a demo example in the demo folder on this github repository

