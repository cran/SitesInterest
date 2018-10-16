#' @name Sites
#' @aliases Sites
#' @title Finds the number of sites of interest from already calculated residence times
#' @description Finds the number of sites of interest as well as other information from already calculated residence times.
#' @usage Sites(Name, R, first = 'n', number_sites = -1)
#'
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param R radius value to use
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param number_sites number of sites to manually show the results for
#'
#' @details This function finds all the necessary information from the results of the already applied algorithm \code{\link{Alt_Alg}}. It returns information relating to both the non-overlapping circles and also the identified sites, which is then used by other functions. The information is extracted from the already saved csv files.
#'
#' @return sites_index - array of indices of the sites of interest among all the circles
#' @return N_no_overlap - number of non-overlapping circles
#' @return X_no_overlap - x-coordinates of non-overlapping circles
#' @return Y_no_overlap - y-coordinates of non-overlapping circles
#' @return X_sites - x-coordinates of identified sites of interest
#' @return Y_sites - y-coordinates of identified sites of interest
#' @return max_percent_drop - maximum percent drop
#' @return number_identified_sites - number of identified sites
#' @return psi_sort_no_overlap2 - ordered list of non-overlapping residence times
#' @export Sites
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#'
#' @seealso See also \code{\link{Alt_Alg}} to find the residence times. 
#' 
#' @examples \donttest{##Find the current working directory
#' wd = getwd()
#' ##Set the working directory as the temporary one
#' setwd(tempdir())
#' ##Set the working directory as the temporary one
#' setwd(tempdir())
#' ##Load the data
#' data(OU_14)
#' t=unlist(OU_14["t"])
#' X=unlist(OU_14["X"])
#' Y=unlist(OU_14["Y"])
#'
#' ##Calculate the residence time with a radius of 0.3 and not including the first circle
#' Alt_Alg("OU_14",t,X,Y,0.3,first='y',save='y')
#'
#' ##Calculate all the necessary information to be used elsewhere
#' Sites("OU_14",0.3,first='y')
#' 
#' ##Reset the original working directory
#' setwd(wd)}
Sites <- function(Name,R, first = 'n',number_sites=-1)
{##Import the already calculated residence times
df = read.csv(paste(Name,"_UD_alt_R",toString(R),".csv",sep=''))

t_list = df[1]
X_list = df[2]
Y_list = df[3]
beta_list = df[4]
phi_list = df[5]
zeta_list = df[6]
psi_list = df[7]

t_centers = unlist(t_list) #the time points that circles are centered on
X_centers = unlist(X_list) #the x-coordinates of circle centers
Y_centers = unlist(Y_list) #the y-coordinates of circle centers
beta = unlist(beta_list) #the number of backward visits
phi = unlist(phi_list) #the number of forward visits
zeta = unlist(zeta_list) #the total number of visits
psi = unlist(psi_list) #the residence times

N_centers=length(t_centers) #the number of circles
if (N_centers!= 1){
  M=cbind(t_centers,X_centers,Y_centers,beta,phi,zeta,psi) #matrix of all the information
  M_sort = M[order(M[,7],M[,1]),] #sort by residence time and then time time stamp if they're the same
  #We now have the following sorted information which will be delted from to get the non-overlapping information
  t_no_overlap = M_sort[,1]
  X_no_overlap = M_sort[,2]
  Y_no_overlap = M_sort[,3]
  beta_no_overlap = M_sort[,4]
  phi_no_overlap = M_sort[,5]
  zeta_no_overlap = M_sort[,6]
  psi_sort_no_overlap = M_sort[,7]
  ##Here we then delete overlapping circles
  for (i in seq(N_centers-1,1,-1)){
    for (j in seq(length(X_no_overlap),i+1,-1)){
      if (sqrt((X_no_overlap[i]-X_no_overlap[j])**2+(Y_no_overlap[i]-Y_no_overlap[j])**2)<=2*R){
        psi_sort_no_overlap = psi_sort_no_overlap[-i]
        t_no_overlap = t_no_overlap[-i]
        X_no_overlap = X_no_overlap[-i]
        Y_no_overlap = Y_no_overlap[-i]
        beta_no_overlap = beta_no_overlap[-i]
        phi_no_overlap = phi_no_overlap[-i]
        zeta_no_overlap = zeta_no_overlap[-i]
        break}}}}else{
          t_no_overlap = t_centers
          X_no_overlap = X_centers
          Y_no_overlap = Y_centers
          beta_no_overlap = beta
          phi_no_overlap = phi
          zeta_no_overlap = zeta
          psi_sort_no_overlap = psi
        }

psi_sort_no_overlap2 = sort(psi_sort_no_overlap, decreasing = TRUE) #the sorted non-overlapping residence times in descending order

N_no_overlap = length(psi_sort_no_overlap2) #the number of non-overlapping circles

for (i in seq(N_no_overlap,1,-1)){
  if (psi_sort_no_overlap[i] == 0){
    psi_sort_no_overlap[i] = psi_sort_no_overlap[-i]
  }
}



if (number_sites==-1){ #check if the user has manually defined the number of sites
if (N_no_overlap == 1){ #if there is only 1 non-overlapping circle then the mpd is 1 and number of sites is 1
  max_percent_drop = 1
  number_identified_sites=1}else{
    Percent_drop = integer(N_no_overlap-1)
    for (i in seq(1,N_no_overlap-1)){
      Percent_drop[i] = 1 -(psi_sort_no_overlap2[i+1]/psi_sort_no_overlap2[i]) #calculate the percentage drop between consecutive non-overlapping circles
      max_percent_drop = max(Percent_drop)
      number_identified_sites=match(max_percent_drop,Percent_drop)}}
if (first == 'y'){ #if instead the user does not want to include the first circle
  if (number_identified_sites == N_no_overlap - 1 & t_no_overlap[1] == t_centers[1]){ #this happens if the only non-overlapping circle is the first one
    Percent_drop = integer(N_no_overlap-2)
    for (i in seq(1,N_no_overlap-2)){
      Percent_drop[i] = 1 -(psi_sort_no_overlap2[i+1]/psi_sort_no_overlap2[i])
      max_percent_drop = max(Percent_drop)
      number_identified_sites=match(max_percent_drop,Percent_drop)}
  }
}}else{
    number_identified_sites = number_sites #the user defined number of sites
    Percent_drop = integer(N_no_overlap-1)
    for (i in seq(1,N_no_overlap-1)){
      Percent_drop[i] = 1 -(psi_sort_no_overlap2[i+1]/psi_sort_no_overlap2[i])} #calculate the percentage drop between consecutive non-overlapping circles
    max_percent_drop = Percent_drop[number_identified_sites]
  }

t_sites = t_no_overlap[seq(N_no_overlap-number_identified_sites+1,N_no_overlap)] #the time points the sites are centered on
X_sites = X_no_overlap[seq(N_no_overlap-number_identified_sites+1,N_no_overlap)] #the x-coordinates of site centers
Y_sites = Y_no_overlap[seq(N_no_overlap-number_identified_sites+1,N_no_overlap)] #the y-coordinates of site centers
sites_index = match(t_sites[seq(1,number_identified_sites)],t_centers)

##All the results are contained in a list
results = list(sites_index,N_no_overlap,X_no_overlap,Y_no_overlap,X_sites,Y_sites,max_percent_drop,number_identified_sites,psi_sort_no_overlap2)
##The list entries are given names
names(results) = list('sites_index','N_no_overlap','X_no_overlap','Y_no_overlap','X_sites','Y_sites','percent_drop','Number_identified_sites','psi_sort_no_overlap')
##The list is returned
  return (results)
}

