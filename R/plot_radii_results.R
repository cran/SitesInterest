#' @name plot_radii_results
#' @aliases plot_radii_results
#' @title Plots the results from using different radii
#' @description Plots the results (maximum percent drop and number of identified sites) from already calculated residence times for different radii values.
#' @usage plot_radii_results(Name, Radii, first = 'n', legend_loc = "topright")
#'
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param Radii set of radius values
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param legend_loc location of the legend
#'
#' @details After finding the residence times for various radii values using \code{\link{Alt_Alg}}, this function plots the number of identified sites of interest and their corresponding maximum percent drops for each radii value. This plot could be used to identify which radius value to use.
#'
#' @return Plot of results from various radius values
#' @export plot_radii_results
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Plots
#'
#' @seealso See also \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the number of identified sites of interest and corresponding maximum percent drop from the csv files produced by \code{\link{Alt_Alg}}.
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
#' ##Run the algorithm for multiple radii values
#' Radii=seq(0.2,1.0,0.1)
#' for (R in Radii){
#'   Alt_Alg("OU_14",t,X,Y,R,first='y',save='y')}
#'
#' ##Plot the results (i.e. maximum percent drop and number of sites identified) from applying
#' ##the algorithm with various radius values
#' plot_radii_results("OU_14",Radii,first='y')
#'
#' ##The location of the legend can be changed
#' plot_radii_results("OU_14",Radii,first='y', legend_loc="bottomleft")
#' 
#' ##Reset the original working directory
#' setwd(wd)}
plot_radii_results <- function(Name,Radii,first='n',legend_loc="topright"){
  N_Radii = length(Radii) #the number of radii values
  Output_Matrix = matrix(0,N_Radii,3) #the matrix of results for each radii value
  for (i in seq(1,N_Radii)){
    Output_Matrix[i,1]=Radii[i]
    Output_Matrix[i,2]=unlist(Sites(Name,Radii[i], first = first)[7])
    Output_Matrix[i,3]=unlist(Sites(Name,Radii[i], first = first)[8])}

  Max_percent_drop = Output_Matrix[,2] #the maximum percent drops
  Number_of_sites = Output_Matrix[,3] #the number of identified sites
  par(mar = c(5, 4, 4, 4) + 0.3, xpd=TRUE)
  ##Plot the maximum percent drops
  plot(Radii,Max_percent_drop*100, type = "o", xlab = "Radii", ylab = "Maximum percent drop (%)", main=(paste(Name," maximum percent drop and number of sites against radii",sep='')), col="blue", pch=1)
  par(new = TRUE) #use the same x-axis
  ##Plot the number of sites
  plot(Radii,Number_of_sites, type = "o", axes = FALSE, xlab = "", bty = "n", ylab = "", col="red", pch=1)
  axis(side=4, at = pretty(range(Number_of_sites)))
  mtext("Number of sites", side=4, line=3)

  ##Add the legend
  legend(legend_loc, legend=c("Maximum percent drop","Number of sites"), pch=c(1,1), col=c("blue","red"), lty = c(1,1))
}
