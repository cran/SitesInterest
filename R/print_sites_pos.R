#' @name print_sites_pos
#' @aliases print.sites_pos
#' @title Prints the positions of the identified sites
#' @description Prints  the positions of the identified sites of interest, which are defined by the radius and circle centre.
#' @usage print_sites_pos(Name,R, first = 'n', number_sites=-1, save="n")
#'
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param R radius value to use
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param number_sites number of sites to manually show the results for
#' @param save if \code{'y'}, the results will be saved as a csv file
#'
#' @details For a given radius value the already identified sites' positions are clearly printed as a table, where their positions are given by the x and y coordinates of the circle centres. There is also the option of saving the data in a csv file, with the title `\emph{Name}'_sites_R`\emph{R}'.csv.
#'
#' @return Prints a summary of the site positions.
#' @export print_sites_pos
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#' @keywords Printed_summary
#'
#' @seealso See \code{\link{Alt_Alg}} to find the residence times. \code{\link{Sites}} can be used to find the coordinates of the site's centres and number of identified sites of interest from the csv files produced by \code{\link{Alt_Alg}}.
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
#' ##Print the coordinates of the centres of all identified sites
#' print_sites_pos("OU_14",0.3,first='y')
#'
#' ##There is also the option of saving the results as a csv file
#' print_sites_pos("OU_14",0.3,first='y')
#' 
#' ##Reset the original working directory
#' setwd(wd)}
print_sites_pos <- function(Name,R, first = 'n', number_sites=-1,save="n")
  {X_sites = unlist(Sites(Name,R, first = first,number_sites=number_sites)[5]) #the x-coordinates of site centers
   Y_sites = unlist(Sites(Name,R, first = first,number_sites=number_sites)[6]) #the y-coordinates of site centers
   Number_identified_sites = Sites(Name,R, first = first,number_sites=number_sites)[[8]] #the number of sites

  A = matrix(c(X_sites,Y_sites),nrow=Number_identified_sites,ncol=2) #matrix of results
  ##Here we name the rows
  row_names = c("Site 1")
  if (Number_identified_sites >= 2){
    for (i in seq(2,Number_identified_sites)){
      row_names = c(row_names, c(paste("Site ",toString(i),sep='')))
    }}
  ##Here we assign the names
  dimnames(A) = list( row_names, c("X coordinate", "Y coordinate"))
  ##The site center coordinates are then printed
  cat("Call:\n")
  cat(paste("The identified sites of interest of ",Name," with a Radius of ",R,".\n",sep=''))
  cat("\nResults table:\n")
  print(A)

  if (save == "y"){
    table = matrix(0,Number_identified_sites,2)
    for (i in seq(1,Number_identified_sites) ){
      table[i,1] = X_sites[i]
      table[i,2] = Y_sites[i]}
      write.table(table, file = paste(Name,"_sites_R",R,".csv"),sep=',', row.names=FALSE, col.names=c("X coordinate", "Y coordinate"))}
  }

