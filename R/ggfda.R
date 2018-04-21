#' This function plots the FECs 
#' @param fdobj FECs
#' @export fdagg
#' @import fda
#' @import dplyr
#' @examples
#' \dontrun{
#' }
fdagg = function(fdobj, ylab='Value', xlab="Time", addzero = FALSE,title="", gridout = TRUE)
{
range= fdobj$basis$rangeval

t = seq(range[1], range[2], len=1e3)

names_y = eval.fd(t, fdobj)%>%colnames

data=  data.frame(Time=t, eval.fd(t, fdobj))
names(data) = c('Time',names_y)

plotdata = reshape2::melt(data, id=1)

plot = ggplot2::ggplot(plotdata)+geom_line(aes(x=Time, y=value,group=factor(variable), color=factor(variable)))+theme_bw()+ylab(ylab)+xlab(xlab)+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
if (addzero) plot = plot+geom_hline(yintercept=0, linetype=2, col=2)
if (gridout) plot = plot+facet_wrap(~variable, scales="free")
return(plot)
}
