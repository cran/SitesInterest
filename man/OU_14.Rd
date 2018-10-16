\name{OU_14}
\alias{OU_14}
\docType{data}
\title{Trajectory of an Ornstein-Uhlenbeck (OU) simulation.
}
\description{The trajectory of a switching OU simulation with 5 points of attraction.
}
\usage{data("OU_14")}
\format{
A data frame of 30000 rows and 3 columns.
}
\details{A switching OU simulation is one where the simulated object will move toward a point of attraction, where the strength of attraction is proportional to the distance from the point. The point of attraction changes at particular points in time. There is also an element of (Gaussian) randomness in the movement. In this particular simulation the points of attraction were (8,3), (5,9), (9,4), (8,2) and	(6,9). The long term standard deviation about these points is 15.689, however it should be noted that the simulation was stopped before reaching this.

\itemize{
\item{t.}{ simulated time recordings.}
\item{X.}{ simulated x-coordinates.}
\item{Y.}{ simulated y-coordinates.}}


}
\source{
Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
}

\examples{
##Load the data
data(OU_14)

}
\keyword{datasets}
