#' @name Alt_Alg_discont
#' @aliases Alt_Alg_discont
#' @title Application of Alt_Alg to discontinuous data
#' @description Calculates the residence times from discontinuous data, which can then be used to identify sites of interest. This could also be used to identify sites for a group of animals, by treating each animal's trajectory as one segment of a discontinuous set.
#' @usage Alt_Alg_discont(Overall_name,Names, t_all, X_all, Y_all, R, s = 10, m = 500, save = 'n')
#'
#' @param Overall_name name for the set of separate trajectories
#' @param Names list of names for each trajectory
#' @param t_all list of arrays, one for each trajectory of times when the positions were recorded
#' @param X_all list of arrays, one for each trajectory of x-coordinates
#' @param Y_all list of arrays, one for each trajectory of y-coordinates
#' @param R radius value to use
#' @param s number of time steps between checks for entrances and exits
#' @param m estimate of the maximum number of crossings across all circles
#' @param save if \code{'y'}, save the files
#'
#' @details This function is used specifically with discontinuous data to calculate the residence times for the trajectories. It works in the same way as for \code{\link{Alt_Alg}}, by linking together \code{\link{Alt_Alg_mini}} and \code{\link{combining}}.
#'
#' @export Alt_Alg_discont
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Discontinuous
#'
#' @seealso  See also \code{\link{Alt_Alg}} for the algorithim used on continuous data. \code{\link{Alt_Alg_mini}} is used to calculate the residence times for a particular set of circles and a particular trajectory, then using \code{\link{combining}} all the residence times for the same circles are summed.
#'
#' @examples ##Find the current working directory
#' wd = getwd()
#' ##Set the working directory as the temporary one
#' setwd(tempdir())
#' ##Load the data
#' data(OU_14)
#' t=unlist(OU_14["t"])
#' X=unlist(OU_14["X"])
#' Y=unlist(OU_14["Y"])
#'
#' ##Number of path sections
#' n=5
#' ##Number of recorded locations
#' N = length(t)
#'
#' ##A list of arrays of the time recoding for the 3 of the trajectory segments
#' t_all = list(t[seq(1,floor(N/n))], t[seq(floor(N/n)*2,floor(N/n)*3)], 
#' t[seq(floor(N/n)*4,floor(N/n)*5)])
#'
#' ##A list of arrays of the x-coordinates for the 3 of the trajectory segments
#' X_all = list(X[seq(1,floor(N/n))], X[seq(floor(N/n)*2,floor(N/n)*3)], 
#' X[seq(floor(N/n)*4,floor(N/n)*5)])
#'
#' ##A list of arrays of the y-coordinates for the 3 of the trajectory segments
#' Y_all = list(Y[seq(1,floor(N/n))], Y[seq(floor(N/n)*2,floor(N/n)*3)], 
#' Y[seq(floor(N/n)*4,floor(N/n)*5)])
#'
#' ##The calculation of the residence time for discontibuous data
#' Alt_Alg_discont("OU_14_discont",c("OU_14.1","OU_14.3","OU14.5"),t_all,X_all,Y_all,0.3,save='y')
#' 
#' ##Reset the original working directory
#' setwd(wd)
Alt_Alg_discont <- function(Overall_name,Names, t_all, X_all, Y_all, R, s = 10, m = 500, save = 'n'){
  N_names = length(Names)
  for (i in seq(1,N_names)){
    Alt_Alg(Names[i],unlist(t_all[i]),unlist(X_all[i]),unlist(Y_all[i]),R,s,m,save=save)
  }

  Circle_names = Names
  Path_names = Names

  for (Circles_name in Circle_names){
    df = read.csv(paste(Circles_name,"_UD_alt_R",R,".csv",sep=''))
    t_centers_list = df[1]
    X_centers_list = df[2]
    Y_centers_list = df[3]
    t_centers = unlist(t_centers_list)
    X_centers = unlist(X_centers_list)
    Y_centers = unlist(Y_centers_list)
    for (Path_name in Path_names){
      if (Circles_name != Path_name){
        index = match(Path_name,Names)
        Alt_Alg_mini(Circles_name, t_centers, X_centers, Y_centers, Path_name, unlist(t_all[index]), unlist(X_all[index]), unlist(Y_all[index]), R, s, m, save=save)}
    }
  }

  combining(Overall_name, Circle_names, Path_names, R)}
