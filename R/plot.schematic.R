#' @name plot.schematic
#' @aliases plot.schematic
#' @title Plots a schematic representation of the movement trajectory
#' @description Plots a schematic representation of the movement from site to site, but does not include returns to the same site.
#' @usage \method{plot}{schematic}(x, y, Name, R, first = 'n', number_sites = -1, len_arrow = 0, lwd_arrow = 0.5, 
#' lwd_r = 1, text_size = 1, legend_loc = "topright", ...)
#' 
#' @method plot schematic
#'
#' @param x array of the x-coordinates describing the trajectory
#' @param y array of the y-coordinates describing the trajectory
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param R radius value to use
#' @param lwd_arrow thickness of the arrows
#' @param len_arrow length of the arrows
#' @param lwd_r width of the hoops
#' @param text_size size of the labels
#' @param legend_loc location of the legend
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param number_sites number of sites to manually show the results for
#' @param ... additional arguments to \link[graphics]{plot}
#'
#' @details This function plots a schematic representation of the animal's movements from site to site, with arrows indicating the direction of movement. This plot simplifies the movements of the animal so that movement is only direct from site to site. Several features of the plot can be defined including; the size of the arrow heads, thickness of the arrows and hoops, text size of the arrow labels and the location of the legend.
#'
#' @return Plot of the schematic representation
#' @export 
#' @exportClass schematic
#' @importFrom plotrix draw.circle
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Plots
#'
#' @seealso See also \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the coordinates of the centres of the sites from the csv files produced by \code{\link{Alt_Alg}}. \code{\link{print_site_visits}} prints the order of sites visited, the length of time for each visit and the amount of time spent between site visits. 
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
#' class(X) = "schematic"
#' class(Y) = "schematic"
#'
#' ##Calculate the residence time with a radius of 0.3 and not including the first circle
#' Alt_Alg("OU_14",t,X,Y,0.3,first='y',save='y')
#'
#' ##Plot the schematic representation of movements between sites
#' plot(X,Y,"OU_14",0.3,first='y')
#'
#' ##There is also the option to make changes to:
#' ##the length of the arrow head
#' plot(X,Y,"OU_14",0.3,first='y',len_arrow=0.25)
#'
#' ##the thickness of the arrow
#' plot(X,Y,"OU_14",0.3,first='y',lwd_arrow=2)
#'
#' ##the thickness of the hoops
#' plot(X,Y,"OU_14",0.3,first='y',lwd_r=2)
#'
#' ##the size of the arrow labels
#' plot(X,Y,"OU_14",0.3,first='y',text_size=2)
#'
#' ##the location of the legend
#' plot(X,Y,"OU_14",0.3,first='y',legend_loc="bottomleft")
#' 
#' ##Reset the original working directory
#' setwd(wd)}
plot.schematic <- function(x,y,Name,R, first = 'n',number_sites=-1,len_arrow=0,lwd_arrow=0.5,lwd_r=1,text_size=1,legend_loc="topright",...){
  N_t = length(x) #the number of positions

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

  Ordered_sites2=integer(k) #order of site visits not including imediate revisits to the same site

  k=1
  for (i in seq(2,Total_site_visits*2-2,2)){
    if (Ordered_sites[i] != Ordered_sites[i+2]){
      Ordered_sites2[k]=Ordered_sites[i]
      k=k+1}}

  Ordered_sites2[length(Ordered_sites2)] = Ordered_sites[Total_site_visits*2-1]

  X_sites_path = integer(length(Ordered_sites2)) #the x-coordinates of the ordered sites
  Y_sites_path = integer(length(Ordered_sites2)) #the y-coordinates of the ordered sites

  for (i in seq(1, length(Ordered_sites2))){
    X_sites_path[i] = X_centers[Ordered_sites2[i]]
    Y_sites_path[i] = Y_centers[Ordered_sites2[i]]}

  ##Calculate the distances between sites
  Distances = integer(length(Ordered_sites2)-1)
  for (i in seq(1, length(Ordered_sites2)-1)){
    Distances[i] = sqrt((X_sites_path[i]-X_sites_path[i+1])**2 + (Y_sites_path[i]-Y_sites_path[i+1])**2)
  }
  #If no arrow head length is defined we take it to be a 10th of the minimum distance between sites
  if (len_arrow == 0){
    len_arrow = 0.1*min(Distances)
  }

  label = seq(1,length(Ordered_sites2)) #the labels for site visits representing the order

  #A function for drawing the arrows pointing from A and B
  drawArrow <- function(A, B){
    arrows(A[1], A[2], x1=0.5*(B[1] + A[1]), y1=0.5*(B[2] + A[2]), col="blue",lwd=lwd_arrow,length=len_arrow)}

  #The x and y coordinates together of the ordered sites
  Sites_coordinates = list(c(X_sites_path[1],Y_sites_path[1]))
  for (i in seq(2,length(Ordered_sites2))){
    Sites_coordinates=c(Sites_coordinates,list(c(X_sites_path[i],Y_sites_path[i])))}

  ##Plot Markov Chain Schematic
  #plot a line between the starting point and first site
  plot(c(x[1],X_sites_path[1]),c(y[1],Y_sites_path[1]),type="l",col="blue", xlim = c(min(x)-R,max(x)+R), ylim = c(min(y)-R,max(y)+R), xlab="X", ylab="Y", main=paste(Name,' schematic (R=',R,')',sep=''), asp=1)
  #plot lines between sites representing transitions between the two sites
  lines(c(X_sites_path[length(Ordered_sites2)],x[N_t]),c(Y_sites_path[length(Ordered_sites2)],y[N_t]),type="l",col="blue")
  #plot a line between the final site and end point
  lines(X_sites_path,Y_sites_path,col="blue")
  #identify site centers with red dots
  points(X_sites_path,Y_sites_path, col="black", pch=21, bg="red")
  #identify the start point with a green dot
  points(x[1],y[1], col="black", pch=21, bg="darkgreen")
  #identify the end point with a yellow dot
  points(x[N_t],y[N_t], col="black", pch=21, bg="yellow")
  #draw circles representing the edges of sites in red
  for (i in seq(N_no_overlap-Number_identified_sites+1,N_no_overlap)){
    plotrix::draw.circle(X_no_overlap[i],Y_no_overlap[i],as.double(R),border='red',col=NA, lwd = lwd_r)}

  #draw arrow between start point and first site if they aren't the same
  if (x[1] != X_sites_path[1] & y[1] != Y_sites_path[1]){
    drawArrow(c(x[1],y[1]),c(X_sites_path[1],Y_sites_path[1]))}
  #draw arrow between the final site and end point  if the aren't the same
  if (x[N_t] != X_sites_path[length(Ordered_sites2)] & y[N_t] != Y_sites_path[length(Ordered_sites2)]){
    drawArrow(c(X_sites_path[length(Ordered_sites2)],Y_sites_path[length(Ordered_sites2)]),c(x[N_t],y[N_t]))}

  #draw arrows between sites to represent transitions
  if (length(Ordered_sites2)-1 >=1){
    for (i in seq(1,length(Ordered_sites2)-1)){
      drawArrow(unlist(Sites_coordinates[i]),unlist(Sites_coordinates[i+1]))}}

  #labels the arrows representing the order of transitions.
  for (i in seq(1,length(label))){
    text((X_sites_path[i+1] + X_sites_path[i])/2, (Y_sites_path[i+1] + Y_sites_path[i])/2, toString(label[i]),col="magenta",cex=text_size)}

  #add the legend
  legend(legend_loc, c("Sites of interest","Start point","End point"),pch=c(21,21,21),col=c("black","black","black"),pt.bg=c("red","darkgreen","yellow"))
}


