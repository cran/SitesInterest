#' @name plot.hoops
#' @aliases plot.hoops
#' @title Plots the trajectory with sites of interest and other hoops
#' @description Plots the trajectory of the animal with the different hoops (sites of interest, removed circles and other circles) overlaid.
#' @usage \method{plot}{hoops}(x, y, Name, R, first = 'n', colours = c('orange','darkgreen','red'), 
#' lwds = c(2,2,2), number_sites=-1, ...)
#' 
#' @method plot hoops
#'
#' @param x array of the x-coordinates describing the trajectory
#' @param y array of the y-coordinates describing the trajectory
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param R radius value to use
#' @param colours list of the hoops' colours
#' @param lwds list of hoop widths
#' @param number_sites number of sites to manually show the results for
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param ... additional arguments to \link[graphics]{plot}
#'
#' @details This function plots the trajectory of the animal with the identified sites of interest, removed circles and other circles plotted on top. The colours of the three types of the circles and their widths can also be defined. This can be used to see visually where the sites are along the trajectory.
#'
#' @return Plot of the identified sites' positions
#' @export 
#' @exportClass hoops
#' @importFrom plotrix draw.circle
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
#' class(X) = "hoops"
#' class(Y) = "hoops"
#'
#' ##Calculate the residence time with a radius of 0.3 and not including the first circle
#' Alt_Alg("OU_14",t,X,Y,0.3,first='y',save='y')
#'
#' ##Plot the positions of the identified sites as well as the non-overlapping circles
#' plot(X,Y,"OU_14",0.3)
#'
#' ##The colours for the hoops can be changed
#' plot(X,Y,"OU_14",0.3,first='y',colours=c('tan','chocolate','maroon'))
#'
#' ##The thickness of hoops can also be changed
#' plot(X,Y,"OU_14",0.3,first='y',lwds=c(0.5,2,3.5))
#'
#' ##It is also possible to manually choose the number of sites
#' plot(X,Y,"OU_14",0.3,first='y',number_sites=4)
#' 
#' ##Reset the original working directory
#' setwd(wd)}
plot.hoops <- function(x,y,Name,R,first='n',colours = c('orange','darkgreen','red'),lwds = c(2,2,2),number_sites=-1,...){
  x=unclass(x)
  y=unclass(y)
  
  ##Import the already calculated residence times
  df = read.csv(paste(Name,"_UD_alt_R",toString(R),".csv",sep=''))
  
  X_list = df[2]
  Y_list = df[3]
  
  X_centers = unlist(X_list) #the x coordinates of circle centers
  Y_centers = unlist(Y_list) #the y coordinates of circle centers
  N_centers = length(X_centers) #the total number of circles
  
  Number_identified_sites = Sites(Name,R, first = first,number_sites=number_sites)[[8]] #the number of sites
  sites_index = Sites(Name,R, first = first, number_sites=number_sites)[1] #the index among all circles of the sites
  N_no_overlap = Sites(Name,R, first = first, number_sites=number_sites)[[2]] #the number of non-overlapping circles
  X_no_overlap = Sites(Name,R, first = first, number_sites=number_sites)[[3]] #the x coordinates of non-overlapping circles
  Y_no_overlap = Sites(Name,R, first = first, number_sites=number_sites)[[4]] #the y coordinates of non-overlapping circles
  X_sites = Sites(Name,R, first = first, number_sites=number_sites)[[5]] #the x coordinates of sites
  Y_sites = Sites(Name,R, first = first, number_sites=number_sites)[[6]] #the y coordinates of sites
  
  ##Plot the path
  plot(x,y,type="l", main=paste(Name,'\n identified sites \n (R=',toString(R),'m)',sep=''), col="blue", asp = 1)
  ##Plot all circles
  for (i in seq(1,N_centers)){
    plotrix::draw.circle(X_centers[i],Y_centers[i],R,border=colours[1],col=NA, lwd = lwds[1])}
  ##Plot the non-overlapping circles
  for (i in seq(1,N_no_overlap)){
    plotrix::draw.circle(X_no_overlap[i],Y_no_overlap[i],R,border=colours[2],col=NA, lwd = lwds[2])}
  ##Plot the sites
  for (i in seq(1,Number_identified_sites)){
    plotrix::draw.circle(X_sites[i],Y_sites[i],R,border=colours[3],col=NA, lwd = lwds[3])}
}


