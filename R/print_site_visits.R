#' @name print_site_visits
#' @aliases print_site_visits
#' @title Prints the site visitation results
#' @description Prints a summary of the site positions and the time spent at and in between each site.
#' @usage print_site_visits(Name, X, Y, R, first = 'n', number_sites = -1, save = 'n')
#'
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param X array of the x-coordinates describing the trajectory
#' @param Y array of the y-coordinates describing the trajectory
#' @param R radius value to use
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param save if \code{'y'}, the results will be saved as a csv file
#' @param number_sites number of sites to manually show the results for
#'
#' @details A summary is printed including the position of the site centres in the order that they are visited. The time spent at each site for each visit is also included as well as the time spent in between sites. There is also the option of saving the data in a csv file, with the title `\emph{Name}'_visits_R`\emph{R}'.csv.
#'
#' @return A summary table of the sites visited.
#' @export print_site_visits
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Printed_summary
#'
#' @seealso See also \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the coordinates of the centres of the sites and non-overlapping circles from the csv files produced by \code{\link{Alt_Alg}}. Then \code{\link{plot.schematic}} gives a visual representation of the order of sites visited.
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
#' ##Calculate the residence time with a radius of 0.3 and not including the first circle
#' Alt_Alg("OU_14",t,X,Y,0.3,first='y',save='y')
#'
#' ##Prints a summary of the site visitation results
#' print_site_visits("OU_14",X,Y,0.3,first='y')
#'
#' ##There is also the option of saving the results as a csv file
#' print_site_visits("OU_14",X,Y,0.3,first='y',save='y')
#' 
#' ##Reset the original working directory
#' setwd(wd)}
print_site_visits <- function(Name,X,Y,R, first = 'n',number_sites=-1,save="n"){
  N_t = length(X) #the number of positions

  #Import the already calculated data
  df = read.csv(paste(Name,"_UD_alt_R",toString(R),".csv",sep=''))

  X_list = df[2]
  Y_list = df[3]

  X_centers = unlist(X_list) #the x-coordinates of the circle centers
  Y_centers = unlist(Y_list) #the y-coordinates of the circle centers

  Number_identified_sites = Sites(Name,R, first = first,number_sites=number_sites)[[8]] #the number of identified sites
  sites_index = Sites(Name,R, first = first,number_sites=number_sites)[[1]] #the index of the sites among all circles
  N_no_overlap = Sites(Name,R, first = first,number_sites=number_sites)[[2]] #the number of non-overlapping circles
  X_no_overlap = Sites(Name,R, first = first,number_sites=number_sites)[[3]] #the x-coordinates of non-overlapping circles
  Y_no_overlap = Sites(Name,R, first = first,number_sites=number_sites)[[4]] #the y-coordinates of non-overlapping circles

  F = read.csv(paste(Name,"_F_alt_R",toString(R),".csv",sep=''), header = FALSE) #the matrix of forward passages entrannce and exit times
  F_col = ncol(F) #the number of columns
  ##This matrix is then reduced so that it only features the ones corresponding to sites
  F_reduced = matrix(0,Number_identified_sites,F_col)
  for (i in seq(1,Number_identified_sites)){
    for (j in seq(1,F_col)){
      F_reduced[i,j] = F[sites_index[i]+1,j]}}


  B = read.csv(paste(Name,"_B_alt_R",toString(R),".csv",sep=''), header = FALSE) #the matrix of backward passages entrannce and exit times
  B_col = ncol(B) #the number of columns
  ##This matrix is then reduced so that it only features the ones corresponding to sites
  B_reduced = matrix(0,Number_identified_sites,B_col)
  for (i in seq(1,Number_identified_sites)){
    for (j in seq(1,B_col)){
      B_reduced[i,j] = B[sites_index[i]+1,j]}}

  df = read.csv(paste(Name,"_UD_alt_R",toString(R),".csv",sep=''))
  zeta_list = df[6]
  zeta = unlist(zeta_list) #the number of visits

  Total_site_visits = 0 #the total number of visits to all sites
  for (i in seq(1,Number_identified_sites)){
    Total_site_visits = Total_site_visits + zeta[sites_index[i]]}

  Ordered_sites = integer(Total_site_visits*2) #array representing the order sites are visited in

  ##Find the backward crossing times for each site
  crossing_times1 = B_reduced[1,]
  if (Number_identified_sites >= 2){
    for (i in seq(2,Number_identified_sites)){
      crossing_times1 = c(crossing_times1,B_reduced[i,])}}

  ##Find the forward crossing times for each site
  crossing_times2 = list(F_reduced[1,])
  if (Number_identified_sites >= 2){
    for (i in seq(2,Number_identified_sites)){
      crossing_times2 = c(crossing_times2,F_reduced[i,])}}

  ##Combine the crossing times of forward and backward visits
  crossing_times = c(crossing_times1,crossing_times2)
  crossing_times = unlist(crossing_times)
  crossing_times = sort(crossing_times) #sort the crossing times
  crossing_times = crossing_times[-seq(1,length(crossing_times)-Total_site_visits*2)] #remove all the zeros

  visit_durations = integer(Total_site_visits)
  for (i in seq(1,Total_site_visits)){
    visit_durations[i] = crossing_times[2*i] - crossing_times[i]
  }
  time_between = integer(Total_site_visits)
  for (i in seq(1,Total_site_visits-1)){
    time_between[i] = crossing_times[2*i+1] - crossing_times[2*i]
  }

  Matrix = cbind(F_reduced,B_reduced) #combine the forward and backward passages matrices

  ##Find the sites corresponding to each passage in the ordered crossing times
  for (i in seq(1,Total_site_visits*2)){
    if (length(which(Matrix == crossing_times[i])[1]) == 1){
      Ordered_sites[i]=sites_index[which(Matrix == crossing_times[i], arr.ind = TRUE)[1]]}
    else{
      if (i%%2 == 1){
        Ordered_sites[i]=Ordered_sites[i-1]}}}

  for (i in seq((Total_site_visits*2-1),0,-1)){
    if (length(which(Matrix != crossing_times[i])[0]) == 1){
      if (Ordered_sites[i] == 0 & i%%2 == 0){
        j=i
        while (Ordered_sites[j] == 0){
          Ordered_sites[j] = Ordered_sites[j+1]
          j=j-1}}}}

  #Calculate the number of site visits not including imediate revists to the same site
  k=1
  for (i in seq(2,Total_site_visits*2)){
    if (Ordered_sites[i] != Ordered_sites[i-1]){
      k=k+1}}

  Ordered_sites3=integer(Total_site_visits)
  for (i in seq(1,Total_site_visits)){
    Ordered_sites3[i] = Ordered_sites[2*i]
  }

  Ordered_sites2=integer(k) #order of site visits not including imediate revisits to the same site

  k=1
  for (i in seq(2,Total_site_visits*2-2,2)){
    if (Ordered_sites[i] != Ordered_sites[i+2]){
      Ordered_sites2[k]=Ordered_sites[i]
      k=k+1}}

  Ordered_sites2[length(Ordered_sites2)] = Ordered_sites[Total_site_visits*2-1]

  Ordered_site_index = Ordered_sites2
  for (i in seq(length(Ordered_sites2),1)){
    if (match(Ordered_sites2[i],Ordered_sites2) != i)
      Ordered_site_index = Ordered_site_index[-i]
  }

  Ordered_sites4=integer(Total_site_visits)
  Ordered_sites4[1]=1
  l=1
  if (Number_identified_sites >= 2){
    for (i in seq(2,Total_site_visits)){
      if (!(match(Ordered_sites3[i],Ordered_site_index) %in%  Ordered_sites4)){
        l=l+1
        Ordered_sites4[i] = l
      }
      else{
        Ordered_sites4[i] = Ordered_sites4[match(Ordered_sites3[i],Ordered_sites3)]
      }
    }}

  X_sites_path = integer(length(Ordered_sites3)) #the x-coordinates of the ordered sites
  Y_sites_path = integer(length(Ordered_sites3)) #the y-coordinates of the ordered sites

  for (i in seq(1, length(Ordered_sites3))){
    X_sites_path[i] = X_centers[Ordered_sites3[i]]
    Y_sites_path[i] = Y_centers[Ordered_sites3[i]]}

  #The x and y coordinates together of the ordered sites
  visits_table = matrix(c(X_sites_path,Y_sites_path,visit_durations,time_between),nrow=length(Ordered_sites3),ncol=4)

  row_names = c("Site 1")
  if (Number_identified_sites >= 2){
    for (i in seq(2,length(Ordered_sites4))){
      row_names = c(row_names, c(paste("Site ",toString(Ordered_sites4[i]),sep='')))
    }}

  dimnames(visits_table) = list( row_names, c("X coordinate", "Y coordinate","Visit duration","Time to next visit"))

  cat("Call:\n")
  cat(paste("The order of visited sites of interest of ",Name," with a Radius of ",R,".\n",sep=''))
  cat("\nSite coordinates:\n")
  print(visits_table)

  if (save == "y"){
    table = matrix(0,Total_site_visits,4)
    for (i in seq(1,Total_site_visits)){
    table[i,1] = X_sites_path[i]
    table[i,2] = Y_sites_path[i]
    table[i,3] = visit_durations[i]
    table[i,4] = time_between[i]}
    write.table(table, file = paste(Name,"_visits_R",R,".csv"),sep=',', row.names=FALSE, col.names=c("X_sites", "Y_sites","Visit_duration","Time_to_next_visit"))
  }
}

