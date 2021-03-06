% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.schematic.R
\name{plot.schematic}
\alias{plot.schematic}
\title{Plots a schematic representation of the movement trajectory}
\usage{
\method{plot}{schematic}(x, y, Name, R, first = 'n', number_sites = -1, len_arrow = 0, lwd_arrow = 0.5, 
lwd_r = 1, text_size = 1, legend_loc = "topright", ...)
}
\arguments{
\item{x}{array of the x-coordinates describing the trajectory}

\item{y}{array of the y-coordinates describing the trajectory}

\item{Name}{name of the data, which is used for any saved files and plot titles}

\item{R}{radius value to use}

\item{first}{if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site}

\item{number_sites}{number of sites to manually show the results for}

\item{len_arrow}{length of the arrows}

\item{lwd_arrow}{thickness of the arrows}

\item{lwd_r}{width of the hoops}

\item{text_size}{size of the labels}

\item{legend_loc}{location of the legend}

\item{...}{additional arguments to \link[graphics]{plot}}
}
\value{
Plot of the schematic representation
}
\description{
Plots a schematic representation of the movement from site to site, but does not include returns to the same site.
}
\details{
This function plots a schematic representation of the animal's movements from site to site, with arrows indicating the direction of movement. This plot simplifies the movements of the animal so that movement is only direct from site to site. Several features of the plot can be defined including; the size of the arrow heads, thickness of the arrows and hoops, text size of the arrow labels and the location of the legend.
}
\examples{
\donttest{##Find the current working directory
wd = getwd()
##Set the working directory as the temporary one
setwd(tempdir())
##Load the data
data(OU_14)
t=unlist(OU_14["t"])
X=unlist(OU_14["X"])
Y=unlist(OU_14["Y"])

class(X) = "schematic"
class(Y) = "schematic"

##Calculate the residence time with a radius of 0.3 and not including the first circle
Alt_Alg("OU_14",t,X,Y,0.3,first='y',save='y')

##Plot the schematic representation of movements between sites
plot(X,Y,"OU_14",0.3,first='y')

##There is also the option to make changes to:
##the length of the arrow head
plot(X,Y,"OU_14",0.3,first='y',len_arrow=0.25)

##the thickness of the arrow
plot(X,Y,"OU_14",0.3,first='y',lwd_arrow=2)

##the thickness of the hoops
plot(X,Y,"OU_14",0.3,first='y',lwd_r=2)

##the size of the arrow labels
plot(X,Y,"OU_14",0.3,first='y',text_size=2)

##the location of the legend
plot(X,Y,"OU_14",0.3,first='y',legend_loc="bottomleft")

##Reset the original working directory
setwd(wd)}
}
\references{
Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
}
\seealso{
See also \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the coordinates of the centres of the sites from the csv files produced by \code{\link{Alt_Alg}}. \code{\link{print_site_visits}} prints the order of sites visited, the length of time for each visit and the amount of time spent between site visits.
}
\author{
Rhys Munden <rdmunden1@sheffield.ac.uk>
}
\keyword{Plots}
