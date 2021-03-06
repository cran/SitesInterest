% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_sites_pos.R
\name{print_sites_pos}
\alias{print_sites_pos}
\alias{print.sites_pos}
\title{Prints the positions of the identified sites}
\usage{
print_sites_pos(Name,R, first = 'n', number_sites=-1, save="n")
}
\arguments{
\item{Name}{name of the data, which is used for any saved files and plot titles}

\item{R}{radius value to use}

\item{first}{if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site}

\item{number_sites}{number of sites to manually show the results for}

\item{save}{if \code{'y'}, the results will be saved as a csv file}
}
\value{
Prints a summary of the site positions.
}
\description{
Prints  the positions of the identified sites of interest, which are defined by the radius and circle centre.
}
\details{
For a given radius value the already identified sites' positions are clearly printed as a table, where their positions are given by the x and y coordinates of the circle centres. There is also the option of saving the data in a csv file, with the title `\emph{Name}'_sites_R`\emph{R}'.csv.
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

##Calculate the residence time with a radius of 0.3 and not including the first circle
Alt_Alg("OU_14",t,X,Y,0.3,first='y',save='y')

##Print the coordinates of the centres of all identified sites
print_sites_pos("OU_14",0.3,first='y')

##There is also the option of saving the results as a csv file
print_sites_pos("OU_14",0.3,first='y')

##Reset the original working directory
setwd(wd)}
}
\references{
Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
}
\seealso{
See \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the coordinates of the site's centres and number of identified sites of interest from the csv files produced by \code{\link{Alt_Alg}}.
}
\author{
Rhys Munden <rdmunden1@sheffield.ac.uk>
}
\keyword{Printed_summary}
