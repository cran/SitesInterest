#' @name Alt_Alg_mini
#' @aliases Alt_Alg_mini
#' @title Calculates the residence times for circles taken from one trajectory, but applied to another
#' @description The circles found from applying \code{\link{Alt_Alg}} to one trajectory are used to find the residence times of another trajectory passing through these circles.
#' @usage Alt_Alg_mini(Circles_name, t_centers, X_centers, Y_centers, Path_name, t, X, Y, R, 
#' s = 10, m = 500, save = 'n')
#'
#' @param Circles_name name of the trajectory used to find the circles
#' @param t_centers array of times when the positions were recorded
#' @param X_centers array of the x-coordinates of the circles' centres
#' @param Y_centers array of the y-coordinates of the circles' centres
#' @param Path_name name of the trajectory that the residence times are found from
#' @param t array of the times that the positions are recorded at
#' @param X array of the x-coordinates describing the trajectory
#' @param Y array of the y-coordinates describing the trajectory
#' @param R radius value to use
#' @param s number of time steps between checks for entrances and exits
#' @param m estimate of the maximum number of crossings across all circles
#' @param save if \code{'y'}, save the files
#'
#' @details This functions works in a similar way to \code{\link{Alt_Alg}}, but the circles are found from one trajectory and are applied to another. The results are stored in a csv file `\emph{Circles_name}'_multi_`\emph{Path_name}'_UD_alt_R`\emph{R}'.csv and the crossing times are stored in `\emph{Circles_name}'_multi_`\emph{Path_name}'_M_alt_R`\emph{R}'.csv 
#'
#' @export Alt_Alg_mini
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Discontinuous
#'
#' @seealso See also \code{\link{Alt_Alg}} for how to apply the algorithm to continuous data. \code{\link{Alt_Alg_discont}} perform \code{\link{Alt_Alg_mini}} on all trjectories, then \code{\link{combining}} combines the results from each application of \code{\link{Alt_Alg_mini}}. 
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
#' ##A list of arrays of the time recoding for the 3 of the path segments
#' t_all = list(t[seq(1,floor(N/n))], t[seq(floor(N/n)*2,floor(N/n)*3)], 
#' t[seq(floor(N/n)*4,floor(N/n)*5)])
#'
#' ##A list of arrays of the x-coordinates for the 3 of the path segments
#' X_all = list(X[seq(1,floor(N/n))], X[seq(floor(N/n)*2,floor(N/n)*3)], 
#' X[seq(floor(N/n)*4,floor(N/n)*5)])
#'
#' ##A list of arrays of the y-coordinates for the 3 of the path segments
#' Y_all = list(Y[seq(1,floor(N/n))], Y[seq(floor(N/n)*2,floor(N/n)*3)], 
#' Y[seq(floor(N/n)*4,floor(N/n)*5)])
#'
#' ##Calculates the residence time for one particular path segment
#' Alt_Alg("OU_14.1",unlist(t_all[1]),unlist(X_all[1]),unlist(Y_all[1]),0.3,first='y',save='y')
#'
#' ##Load the data of the circles found from Alt_Alg
#' df = read.csv(paste("OU_14.1","_UD_alt_R",0.3,".csv",sep=''))
#' t_centers = unlist(df[1])
#' X_centers = unlist(df[2])
#' Y_centers = unlist(df[3])
#'
#' ##Calculates the residence time from path segment 3, using circles from path segment 1
#' Alt_Alg_mini("OU14.1", t_centers, X_centers, Y_centers, "OU_14.3", unlist(t_all[2]), 
#' unlist(X_all[2]), unlist(Y_all[2]), 0.3,save='y')
#' 
#' ##Reset the original working directory
#' setwd(wd)
Alt_Alg_mini <-function(Circles_name, t_centers, X_centers, Y_centers, Path_name, t, X, Y, R, s=10, m=500, save = 'n')
{N=length(t)
T_c=0 #time threshold between visits
N_centers = length(X_centers)
M = matrix(0,N_centers,m)##mth passage across circle perimeter
##run through the rest of the algorithm in two parts
for (q in seq(1,N_centers)){
  l_f=1
  k_f=2
  if (sqrt((X_centers[q]-X[1])**2 + (Y_centers[q]-Y[1])**2) <= R){
    M[q,1] = t[1]
    l_f=l_f+2
  }
  #Calculate ranges
  if (N%%s == 0){
    end = N-s
  }else{
    end = N-N%%s
  }
  for (i in seq(1,end,s)){
    #exits
    if (sqrt((X[i]-X_centers[q])**2+(Y[i]-Y_centers[q])**2) <= R & sqrt((X[i+s]-X_centers[q])**2+(Y[i+s]-Y_centers[q])**2) >= R){
      for (j in seq(i,i+s+1)){
        if (sqrt((X[j]-X_centers[q])**2+(Y[j]-Y_centers[q])**2) <= R & sqrt((X[j+1]-X_centers[q])**2+(Y[j+1]-Y_centers[q])**2) >= R){
          if (R - sqrt((X[j]-X_centers[q])**2+(Y[j]-Y_centers[q])**2) <= sqrt((X[j+1]-X_centers[q])**2+(Y[j+1]-Y_centers[q])**2) - R){
            M[q,k_f]=t[j]
            k_f=k_f+2
            break
          }else{
            M[q,k_f]=t[j+1]
            k_f=k_f+2
            break
          }
        }
      }
    }
    #entries
    if (sqrt((X[i]-X_centers[q])**2+(Y[i]-Y_centers[q])**2) >= R & sqrt((X[i+s]-X_centers[q])**2+(Y[i+s]-Y_centers[q])**2) <= R){
      for (j in seq(i,i+s+1)){
        if (sqrt((X[j]-X_centers[q])**2+(Y[j]-Y_centers[q])**2) >= R & sqrt((X[j+1]-X_centers[q])**2+(Y[j+1]-Y_centers[q])**2) <= R){
          if (sqrt((X[j]-X_centers[q])**2+(Y[j]-Y_centers[q])**2) - R <= R - sqrt((X[j+1]-X_centers[q])**2+(Y[j+1]-Y_centers[q])**2)){
            M[q,l_f]=t[j]
            l_f=l_f+2
            break
          }else{
            M[q,l_f]=t[j+1]
            l_f=l_f+2
            break
          }
        }
      }
    }
  }
  #solve problem at the edges
  for (j in seq(end,N-1)){
    #exit
    if (sqrt((X[j]-X_centers[q])**2+(Y[j]-Y_centers[q])**2) <= R & sqrt((X[j+1]-X_centers[q])**2+(Y[j+1]-Y_centers[q])**2) >= R){
      if (R - sqrt((X[j]-X_centers[q])**2+(Y[j]-Y_centers[q])**2) <= sqrt((X[j+1]-X_centers[q])**2+(Y[j+1]-Y_centers[q])**2) - R){
        M[q,k_f]=t[j]
        k_f=k_f+2
      }else{
        M[q,k_f]=t[j+1]
        k_f=k_f+2
      }
    }
    #entry
    if (sqrt((X[j]-X_centers[q])**2+(Y[j]-Y_centers[q])**2) >= R & sqrt((X[j+1]-X_centers[q])**2+(Y[j+1]-Y_centers[q])**2) <= R){
      if (sqrt((X[j]-X_centers[q])**2+(Y[j]-Y_centers[q])**2) - R <= R - sqrt((X[j+1]-X_centers[q])**2+(Y[j+1]-Y_centers[q])**2)){
        M[q,l_f]=t[j]
        l_f=l_f+2
      }else{
        M[q,l_f]=t[j+1]
        l_f=l_f+2
      }
    }
  }
  if (sqrt((X_centers[q]-X[N])**2 + (Y_centers[q]-Y[N])**2) <= R){
    M[q,k_f] = t[N]
  }
}
#Here I calculate the number of centers
C=N_centers
##Number of visits
zeta = integer(C)##number of visits
for (i in seq(1,C)){
  if (M[i,1] != 0){
    zeta[i] = 1
    for (j in seq(1,m/2)){
      if (M[i,2*j] != 0){
        if (M[i,2*j]-M[i,2*j-1] >= T_c){
          zeta[i] = zeta[i] + 1
        }
      }
    }
  }
}

psi = integer(C)##residence time
##We will calculate the residence time
for (i in seq(1,C)){
  for (k in seq(1,m/2-1)){
    psi[i]=psi[i]+M[i,2*k]-M[i,2*k-1]
  }
}

##Number of visits
zeta = integer(C)##number of visits
for (i in seq(1,C)){
  if (M[i,1] != 0){
    zeta[i] = 1
    for (j in seq(1,m/2)){
      if (M[i,2*j] != 0){
        if (M[i,2*j+1]-M[i,2*j] >= T_c){zeta[i] = zeta[i] + 1}
      }
    }
  }
}



m2=2*max(zeta)
M2 = matrix(0,C,m2)
##Export the data##
a = matrix(1,C,5)
a[,1]=t_centers
a[,2]=X_centers
a[,3]=Y_centers
for (i in seq(1,C)){
  a[i,4]=zeta[i]
  a[i,5]=psi[i]
}

for (i in seq(1,C)){
  for (j_b in seq(0,m2)){
    M2[i,j_b]=M[i,j_b]
  }
}

if (save == 'y'){
  write.table(a, file = paste(Circles_name,"_multi_",Path_name,"_UD_alt_R",R,".csv",sep=''), sep=",", row.names=FALSE, col.names=c("t","X","Y","zeta","psi"))
  write.table(M2, file = paste(Circles_name,"_multi_",Path_name,"_M_alt_R",R,".csv",sep=''), sep=",", row.names=FALSE, col.names=seq(1,m2))}}
