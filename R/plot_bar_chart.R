#' @name plot_bar_chart
#' @aliases plot_bar_chart
#' @title Plot the bar chart of residence times
#' @description This function plots the bar chart of residence times with the ones that correspond to the sites of interest highlighted by a change in colour.
#' @usage plot_bar_chart(Name, R, first='n', number_sites=-1, colours= c("red","darkgreen"))
#'
#' @details Plots the bar chart of ranked residency times of non-overlapping circles. The percent drops are the relative differences between these consecutive bars and where the maximum percent drop occurs is identified using a change in colour between the bars. This plot can be used to see the percent drops visually and if the maximum percent drop can easily be identified by sight.
#'
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param R radius value to use
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param number_sites number of sites to manually show the results for
#' @param colours list of the bars' colours
#' 
#' @return Bar chart plot of the the residence times of non-overlapping circles
#' @export plot_bar_chart
#' @import graphics
#'
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Plots
#'
#' @seealso See also \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the coordinates of the centres of the sites and non-overlapping circles from the csv files produced by \code{\link{Alt_Alg}}.
#' 
#' @examples \donttest{##Find the current working directory
#' wd = getwd()
#' ##Set the working directory as the temporary one
#' setwd(tempdir())
#' ##Load the data
#' data(OU_14)
#' t=unlist(OU_14["t"])
#' X=unlist(OU_14["X"])
#' Y=unlist(OU_14["Y"])
#'
#' ##Calculate the residence time with a radius of 0.3 and not including the first circle
#' Alt_Alg("OU_14",t,X,Y,0.3,first='y',save='y')
#'
#' ##Plot the bar chart of ranked non-overlapping residence times
#' plot_bar_chart("OU_14",0.3,first='y')
#'
#' ##It is possible to choose manually where the cut off between sites and none sites should be
#' plot_bar_chart("OU_14",0.3,first='y', number_sites=4)
#'
#' ##The colours can also be changed
#' plot_bar_chart("OU_14",0.3,first='y', colours = c("darkgreen","red"))
#' 
#' ##Reset the original working directory
#' setwd(wd)}
plot_bar_chart <- function(Name,R,first='n',number_sites=-1,colours=c("red","darkgreen")){
  N_no_overlap = Sites(Name,R, first = first, number_sites=number_sites)[[2]] #the number of non-overlapping circles
  Number_identified_sites = Sites(Name,R, first=first, number_sites=number_sites)[[8]] #the number of identified sites
  psi_sort_no_overlap2 = Sites(Name,R, first=first, number_sites=number_sites)[[9]] #the sorted (descending order) and non-overlapping residence times
  
  vals = 0:N_no_overlap #the ranks of the circles
  breaks = c(-Inf, Number_identified_sites, Inf) #where the colour change ocurs
  cols = colours[findInterval(vals, vec=breaks)] #the colours used for the bars
  #Plot the bar chart
  barplot(psi_sort_no_overlap2, main=paste(Name,"\n circle residence times with no overlap \n (R=",toString(R),'m)',sep=''), xlab="Rank", ylab=paste("Residence time ",expression(psi),sep=''),col=cols, names.arg=seq(1,N_no_overlap))
}



