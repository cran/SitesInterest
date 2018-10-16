#' @name plot_bars_and_hoops
#' @aliases plot_bars_and_hoops
#' @title Plots both the bar chart and hoops plot in one figure
#' @description The bar chart is the plot of the ranked non-overlapping residence times and the hoops plot shows where along the trajectory these different circles are located.
#' @usage plot_bars_and_hoops(x, y, Name, R, colours = c('orange','darkgreen','red'), lwds = 
#' c(2,2,2), number_sites=-1, first='n')
#'
#' @param x array of the x-coordinates describing the trajectory
#' @param y array of the y-coordinates describing the trajectory
#' @param Name the name of the data used to save the files and also in the title of the plot
#' @param R radius value to use
#' @param colours list of colours used for the bars and hoops
#' @param lwds list of the widths of the circles
#' @param number_sites number of sites to manually show the results for
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#'
#' @details The first subplot (\code{\link{plot_bar_chart}}) shows the residence times of circles that are left after overlapping ones are removed and identifies where the maximum percent drop occurs by a change in the colour of bars.
#' 
#' The second subplot (\code{\link{plot.hoops}}) includes the animal's trajectory as well as the three different types of circles overlaid (sites, non-overlapping circles and removed circles). The colour and thickness of these circles can be changed.
#'
#' @return Plot of the bar chart of residence times and site positions
#' @export plot_bars_and_hoops
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Plots
#'
#' @seealso See also \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the coordinates of the centres of the sites and non-overlapping circles from the csv files produced by \code{\link{Alt_Alg}}. See \code{\link{plot_bar_chart}} and \code{\link{plot.hoops}} for how to plot each separately and also the different presentation options.
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
#' class(X) <- "bars_and_hoops"
#' class(Y) <- "bars_and_hoops"
#'
#' ##Calculate the residence time with a radius of 0.3 and not including the first circle
#' Alt_Alg("OU_14",t,X,Y,0.3,first='y',save='y')
#'
#' ##Plot the bar chart of ranked non-overlapping residence times and the plot showing
#' ##the positions of these circles
#' plot_bars_and_hoops(X,Y,"OU_14",0.3,first='y')
#' 
#' ##Reset the original working directory
#' setwd(wd)}
plot_bars_and_hoops <- function(x,y,Name,R,colours = c('orange','darkgreen','red'),lwds = c(2,2,2),number_sites=-1,first='n'){
  par(mfrow=c(1,2)) #the subplots are positioned side by side
  plot_bar_chart(Name,R,number_sites=number_sites,first=first,colours=c(colours[3],colours[2])) #plot the bar chart
  class(x) = "hoops"
  class(y) = "hoops"
  plot(x,y,Name,R,colours = colours,lwds = lwds, number_sites=number_sites,first=first) #plot the site positions
  par(mfrow=c(1,1)) #any plots afterwards will not be subplots
}


