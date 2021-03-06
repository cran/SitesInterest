#' @name combining
#' @aliases combining
#' @title Combines the residence times for the same set of circles and different trajectories
#' @description Combines the residence times for the same set of circles and different trajectories and also for the entire set.
#' @usage combining(Overall_name, Circle_names, Path_names, R, save = 'n')
#'
#' @param Overall_name name for the set of separate trajectories
#' @param Circle_names the list of names for each set of circles
#' @param Path_names the list of names for each trajectory
#' @param R radius value to use
#' @param save if \code{'y'}, save the files 
#'
#' @details This function combines the already calculated residence times for each of the path segments as well as combining the residence times of all path segments across the same set of circles. The results are stored in a csv file named `\emph{Overall_name}'_combined_UD_alt_R`\emph{R}'.csv.
#'
#' @export combining
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Discontinuous
#'
#' @seealso See also \code{\link{Alt_Alg_mini}} for how to calculate the residence times for a particular set of circles and a particular trajectory. \code{\link{Alt_Alg_discont}} is used to calculate the residence times across a whole set of discontinuous trajectories.
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
#' #Number of recorded locations
#' N = length(t)
#'
#' ##A list of arrays of the time recoding for the 3 of the path segments
#' t_all = list(t[seq(1,floor(N/n))], t[seq(floor(N/n)*2,floor(N/n)*3)], 
#' t[seq(floor(N/n)*4,floor(N/n)*5)])
#' ##A list of arrays of the x-coordinates for the 3 of the path segments
#' X_all = list(X[seq(1,floor(N/n))], X[seq(floor(N/n)*2,floor(N/n)*3)], 
#' X[seq(floor(N/n)*4,floor(N/n)*5)])
#' ##A list of arrays of the y-coordinates for the 3 of the path segments
#' Y_all = list(Y[seq(1,floor(N/n))], Y[seq(floor(N/n)*2,floor(N/n)*3)], 
#' Y[seq(floor(N/n)*4,floor(N/n)*5)])
#'
#' ##The names of each path segment
#' Names = c("OU_14.1","OU_14.3","OU_14.5")
#'
#' ##Calculates the residence time for each path segment individually
#' Alt_Alg("OU_14.1",unlist(t_all[1]),unlist(X_all[1]),unlist(Y_all[1]),0.3,first='y',save='y')
#' Alt_Alg("OU_14.3",unlist(t_all[2]),unlist(X_all[2]),unlist(Y_all[2]),0.3,first='y',save='y')
#' Alt_Alg("OU_14.5",unlist(t_all[3]),unlist(X_all[3]),unlist(Y_all[3]),0.3,first='y',save='y')
#'
#' Circle_names = Names
#' Path_names = Names
#'
#' ##Calculate the residence time for each set of circles and each path segment
#' for (Circles_name in Circle_names){
#'   df = read.csv(paste(Circles_name,"_UD_alt_R",0.3,".csv",sep=''))
#'   t_centers_list = df[1]
#'   X_centers_list = df[2]
#'   Y_centers_list = df[3]
#'   t_centers = unlist(t_centers_list)
#'   X_centers = unlist(X_centers_list)
#'   Y_centers = unlist(Y_centers_list)
#'   for (Path_name in Path_names){
#'   if (Circles_name != Path_name){
#'        index = match(Path_name,Names)
#'        Alt_Alg_mini(Circles_name, t_centers, X_centers, Y_centers, Path_name, 
#'        unlist(t_all[index]),unlist(X_all[index]),unlist(Y_all[index]),0.3,s=10,m=500,save='y')}}}
#'
#' ##Combine all the residence times for the same circles
#' combining("OU_14_discont", Circle_names, Path_names, 0.3,save='y')
#' 
#' ##Reset the original working directory
#' setwd(wd)
combining <-function(Overall_name, Circle_names, Path_names, R, save = 'n'){
  matrix = matrix(, nrow = 0, ncol = 5)
  for (circle_name in Circle_names){
    df = read.csv(paste(circle_name,"_UD_alt_R",R,".csv",sep=''))

    t_centers_list = df[1]
    X_centers_list = df[2]
    Y_centers_list = df[3]
    zeta_list = df[6]
    psi_list = df[7]

    t_centers = unlist(t_centers_list)
    X_centers = unlist(X_centers_list)
    Y_centers = unlist(Y_centers_list)
    zeta = unlist(zeta_list)
    psi = unlist(psi_list)

    N_centers = length(t_centers)
    Matrix = matrix(0,N_centers,5)
    Matrix[,1] = t_centers
    Matrix[,2] = X_centers
    Matrix[,3] = Y_centers
    Matrix[,4] = zeta
    Matrix[,5] = psi
    
    for (path_name in Path_names){
      if (circle_name != path_name){
        df = read.csv(paste(circle_name,"_multi_",path_name,"_UD_alt_R",R,".csv",sep=''))
        zeta_list = df[4]
        psi_list = df[5]
        zeta = unlist(zeta_list)
        psi = unlist(psi_list)
        Matrix[,4] = Matrix[,4] + zeta
        Matrix[,5] = Matrix[,5] + psi
        if (save == 'y'){
          write.table(Matrix, file = paste(circle_name,"_combined_UD_alt_R",R,".csv",sep=''), sep=",", row.names=FALSE, col.names=c("t","X","Y","zeta","psi"))}}}
    matrix = rbind(matrix,Matrix)
    }
  if (save == 'y'){
    write.table(matrix, file = paste(Overall_name,"_combined_UD_alt_R",R,".csv",sep=''), sep=",", row.names=FALSE, col.names=c("t","X","Y","zeta","psi"))
    }}
  