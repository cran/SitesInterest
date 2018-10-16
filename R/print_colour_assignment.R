#' @name print_colour_assignment
#' @aliases print_colour_assignment
#' @title Prints a summary table of results from the different criteria
#' @description Prints a summary of the results from using the stability criterion and threshold criterion as well as the colour assigned.
#' @usage print_colour_assignment(Name, Threshold, Radii, first = 'n', number_sites = -1)
#'
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param Threshold  threshold value
#' @param Radii set of radius values
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param number_sites number of sites to manually show the results for
#'
#' @details This function prints a clear summary of the results from using both the stability and threshold criteria, including the maximum percent drop and number of identified sites. The threshold criterion requires that the maximum percent drop must be greater than a given value (\code{Threshold}) and also be a local maximum. The stability criterion requires again that the maximum percent drop be a local maximum and also that the radius values either side, result in the same number of sites identified.
#' 
#' The colour assigned is also printed, where Green is assigned if the two criteria result in the same radius, Amber is assigned if the number of sites are the same for the two criteria, but not the same radius and Red is assigned if the number of sites are different. This gives a qualitative level of confidence in the results produced.
#'
#' @return Summary table of the colour assignment part of the algorithm.
#' 
#' @export print_colour_assignment
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Printed_summary
#'
#' @seealso See also \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the maximum percent drop and number of sites of interest from the csv files produced by \code{\link{Alt_Alg}}.
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
#' ##Run the algorithm for multiple radius values
#' Radii=seq(0.2,1.0,0.1)
#' for (R in Radii){
#'   Alt_Alg("OU_14",t,X,Y,R,first='y',save='y')}
#'
#' ##Print a summary table of the results from the two criteria and the colour assigned
#' print_colour_assignment("OU_14",65,Radii,first='y')
#' 
#' ##Reset the original working directory
#' setwd(wd)}
print_colour_assignment <- function(Name,Threshold,Radii,first='n',number_sites=-1)
{ Threshold=Threshold/100 #calculate the threshold value as a decimal
  N_Radii = length(Radii) #the number of radii values
  Output_Matrix = matrix(0,N_Radii,3) #matrix to fill with results from the diffent mpd and number of sites for each radius
  #the matrix is then filled out
  for (i in seq(1,N_Radii)){
    Output_Matrix[i,1]=Radii[i]
    Output_Matrix[i,2]=unlist(Sites(Name,Radii[i], first = first, number_sites=number_sites)[7])
    Output_Matrix[i,3]=unlist(Sites(Name,Radii[i], first = first, number_sites=number_sites)[8])}

  Max_percent_drop = Output_Matrix[,2]
  Number_of_sites = Output_Matrix[,3]

  N = length(Radii)

  ##Threshold criteria
  Thresh_max_percent_drop=0 #the maximum percent drop found using the threshold criteria
  for (k in seq(1,N-1)){
    if (Max_percent_drop[k]>=Max_percent_drop[k+1]){ #check that it is a local maximum
      if(Max_percent_drop[k] >= Threshold){ #check that it exceeds the threshold
        Thresh_max_percent_drop = Max_percent_drop[k]
        Thresh_radius = Radii[k] #the radius associated with the above mapd
        Thresh_number_of_sites = Number_of_sites[k] #the number of sites identified using the above radius
        break}}}

  ##First Max and stable
  Stable_max_percent_drop=0 #the maximum percent drop found using the stability criteria
  for (k in seq(2,N-1)){
    if (Max_percent_drop[k-1]<=Max_percent_drop[k] & Max_percent_drop[k]>=Max_percent_drop[k+1]){ #check that it is a local maximum
      if (Number_of_sites[k-1] == Number_of_sites[k] & Number_of_sites[k+1] == Number_of_sites[k]){ #check that it is stable
        Stable_max_percent_drop = Max_percent_drop[k]
        Stable_radius = Radii[k] #the radius associated with the above mpd
        Stable_number_of_sites = Number_of_sites[k] #the number of sites identified using the above radius
        break}}}

  ##Assign a colour
  if (Stable_max_percent_drop != 0){
  if (Stable_number_of_sites != Thresh_number_of_sites){
    colour_assignment = "Red"}else{ #red is assigned if the number of sites are different
      if (Stable_radius != Thresh_radius){
        colour_assignment = "Amber"}else{ #amber is assigned if the number of sites are the same, but the radii are different
          colour_assignment = "Green"}} #green is assigned if booth criteria agree

  #The results are then printed
  A = matrix(c(Thresh_radius,Thresh_max_percent_drop*100,as.integer(Thresh_number_of_sites),Stable_radius,Stable_max_percent_drop*100,as.integer(Stable_number_of_sites)),nrow=3,ncol=2)
  dimnames(A) = list(c("Radius","Maximum percent drop","Number of sites"), c("Threshold", "Stability"))
  cat("Call:\n")
  cat(paste("The results of ",Name," with a threshold value of ",Threshold*100,"%.\n",sep=''))
  cat("\nResults table:\n")
  print(A)
  cat("\nColour category assigned:\n")
  cat(colour_assignment)}else{
    print("No stable radius was found.\n")
  }
}


