#' @name plot.displacement
#' @aliases plot.displacement
#' @title Plots the displacement from a particular point over a chosen window
#' @description The displacement from a particular point to every other point along the animal's trajectory is calculated and this is then plotted over a defined window.
#' @usage \method{plot}{displacement}(x, y, Name, t, R, t_a, t_1, t_2, ...)
#'
#' @method plot displacement
#'
#' @param x array of the x-coordinates describing the trajectory
#' @param y array of the y-coordinates describing the trajectory
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param t an array of the times that the positions are recorded at
#' @param R radius value to use
#' @param t_a starting time to calculate the displacement from
#' @param t_1 start of interval to view the displacement over
#' @param t_2 end of interval to view the displacement over
#' @param ... additional arguments to \link[graphics]{plot}
#' 
#' @details The displacement from a particular time point (\code{t_a}) to every other point along the trajectory is calculated and this is then plotted over a defined window [\code{t_1}, \code{t_2}]. A line representing the radius is also drawn to see when the trajectory enters or leaves the circle centred at the point (\code{t_a}). This plot can be used to see how far away the animal moves after leaving a particular circle.
#'
#' @return Plot of the displacement
#' @export 
#' @exportClass displacement
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Plots
#'
#' 
#' @examples ##Load the data
#' data(OU_14)
#' t=unlist(OU_14["t"])
#' X=unlist(OU_14["X"])
#' Y=unlist(OU_14["Y"])
#' 
#' class(X) = "displacement"
#' class(Y) = "displacement"
#'
#' ##Plot the displacement from the starting point (t=0) for t=0 to t=2.9999
#' plot(X, Y, "OU14", t, 0.3, 0, 0, 2.9999)
plot.displacement <- function(x, y, Name, t, R, t_a, t_1, t_2, ...)
{m1 = match(t_1,t) #find the index of the starting point of the domain
m2 = match(t_2,t) #find the index of the end point of the domain
t_a_index = match(t_a,t) #find the index of the point from which the displacement is measured
N=length(t) #the number of time points

D = integer(N) #array of the displacements
#the displacement is calculated using the Euclidean metric
for (i in seq(1,N)){
  D[i] = sqrt((x[i] - x[t_a_index])**2 + (y[i] - y[t_a_index])**2)}

#plot the displacement over the time window
plot(t[m1:m2], D[m1:m2], col='blue', main=paste(Name,' displacement plot (R=',R,")",sep=''), xlab='t', ylab='Displacement', type='l')
#plot a red line representing the circle's edge
lines(c(t[m1], t[m2]), c(R,R), col='red')}


