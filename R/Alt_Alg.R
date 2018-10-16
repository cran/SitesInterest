#' @name Alt_Alg
#' @aliases Alt_Alg
#' @title Calculates the residence times and identifies the sites of interest
#' @description From positional data, this function finds the residence times as well as the number of visits. It also identifies which areas are sites of interest.
#' @usage Alt_Alg(Name, t, X, Y, R, s = 10, m = 500, first = 'n', save = 'n')
#'
#' @param Name name of the data, which is used for any saved files and plot titles
#' @param t array of the times that the positions are recorded at
#' @param X array of the x-coordinates describing the trajectory
#' @param Y array of the y-coordinates describing the trajectory
#' @param R radius value to use
#' @param s number of time steps between checks for entrances and exits
#' @param m estimate of the maximum number of crossings across all circles
#' @param first if \code{'y'}, the algorithm will look for the second greatest maximum percent drop if the first results in the first circle being the only non-identified site
#' @param save if \code{'y'}, save the files
#'
#' @details This is the main part of the algorithm. Given the trajectory of an animal, the algorithm centres a circle at the start of the trajectory. The next circle is centred at the point when the animal first leaves the previous circle. This is then repeated all along the trajectory until the whole trajectory is covered. The algorithm then finds the times of when the animal has crossed the perimeter of these circles, the results of which are saved in two csv files (comma separated values), with titles `\emph{Name}'_F_alt_R`\emph{R}'.csv and `\emph{Name}'_B_alt_R`\emph{R}'.csv for the forward and backward crossings respectively. A forward crossing is one that occurs after the time point that the circle is centred on and backward crossings occur before this.
#' 
#' The residence times and number of visits for each circle are calculated, which is also stored in a csv file with the title `\emph{Name}'_UD_alt_R`\emph{R}'.csv. The function orders all the circles in descending order of residence time and then moving down the list, any that overlap with one of a higher residence time are removed. The relative difference between consecutive residence times (percent drops) are calculated and the maximum is identified (\code{maximum_percent_drop}). The number of identified sites (\code{Number_identified_sites}) is the number of circles before this maximum percent drop. This can be seen visually, using \code{\link{plot_bar_chart}}.
#'
#' @return max_percent_drop - maximum percent drop, which is the relative difference between the bars, indicating the difference between sites and non-site circles
#' @return number_identified_sites - number of circles identified as sites of interest
#' @export Alt_Alg
#'
#' @import utils
#'
#' @references Munden, R., Borger , L., Wilson, R.P., Redcliffe, J., Loison, A., Garel, M. and Potts, J.P. in review. Making sense of ultra-high-resolution movement data: an algorithm for inferring sites of interest.
#' @author Rhys Munden <rdmunden1@sheffield.ac.uk>
#'
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
#' ##Reset the original working directory
#' setwd(wd)}
Alt_Alg <- function(Name, t, X, Y, R, s = 10, m = 500, first = 'n', save = 'n')
{N = length(t)
T_c = 0 #time threshold between visits
F = matrix(0,N,m)##forward mth passage across circle perimeter
B = matrix(0,N,m)##backward mth passage across circle perimeter
P=integer(N) #array identifying the circle centers
P[1]=1#array identifying the circle centers

##Check the R value is large enough
d = integer(N)
for (i in seq(1,N-1)){
  d[i] = sqrt((X[i]-X[i+1])**2 + (Y[i]-Y[i+1])**2)}
R_min = max(d)
if (R < R_min){
  print(paste('Your value for R is too small. It should be greater than ',R_min))}

##Choose the first point
for( i in seq(1,N-s,s)){
  if (sqrt((X[i]-X[1])**2+(Y[i]-Y[1])**2) <= R & sqrt((X[i+s]-X[1])**2+(Y[i+s]-Y[1])**2) >= R){
    for (j in seq(i,i+s)){
      if (sqrt((X[j]-X[1])**2+(Y[j]-Y[1])**2) <= R & sqrt((X[j+1]-X[1])**2+(Y[j+1]-Y[1])**2) >= R){
        if (R - sqrt((X[j]-X[1])**2+(Y[j]-Y[1])**2) <= sqrt((X[j+1]-X[1])**2+(Y[j+1]-Y[1])**2) - R){
          F[1,1]=t[j]
          P[2]=j
          break}else{
            F[1,1]=t[j+1]
            P[2]=j+1
            break}}}
    break}}

##find all the forward passages from the first point
if (P[2] != 0){
k=3
l=2
if (P[2]%%s == 1){
  if (sqrt((X[P[2]+s]-X[1])**2+(Y[P[2]+s]-Y[1])**2) >= R){
  Range = seq(P[2]+s, N-2*s+1,s)}else{
    if (sqrt((X[P[2]]-X[1])**2+(Y[P[2]]-Y[1])**2) >= R & sqrt((X[P[2]+1]-X[1])**2+(Y[P[2]+1]-Y[1])**2) <= R){
      Range = seq(P[2], N-2*s+1,s)}else{
    Range = seq(P[2]+1, N-2*s+1,s)}
  }
}else{
  Range = seq(ceiling(P[2]/s)*s+1, N-2*s+1, s)}
for(i in Range){
  if (sqrt((X[i]-X[1])**2+(Y[i]-Y[1])**2) <= R & sqrt((X[i+s]-X[1])**2+(Y[i+s]-Y[1])**2) >= R){
    for (j in seq(i,i+s)){
      if (sqrt((X[j]-X[1])**2+(Y[j]-Y[1])**2) <= R & sqrt((X[j+1]-X[1])**2+(Y[j+1]-Y[1])**2) >= R){
        if (R - sqrt((X[j]-X[1])**2+(Y[j]-Y[1])**2) <= sqrt((X[j+1]-X[1])**2+(Y[j+1]-Y[1])**2) - R){
          F[1,k]=t[j]
          k=k+2
          break}else{
            F[1,k]=t[j+1]
            k=k+2
            break}}}}
  if (sqrt((X[i+s]-X[1])**2+(Y[i+s]-Y[1])**2) <= R & sqrt((X[i]-X[1])**2+(Y[i]-Y[1])**2) >= R){
    for( j in seq(i,i+s)){
      if (sqrt((X[j]-X[1])**2+(Y[j]-Y[1])**2) >= R & sqrt((X[j+1]-X[1])**2+(Y[j+1]-Y[1])**2) <= R){
        if (sqrt((X[j]-X[1])**2+(Y[j]-Y[1])**2) - R <= R - sqrt((X[j+1]-X[1])**2+(Y[j+1]-Y[1])**2)){
          F[1,l]=t[j]
          l=l+2
          break}else{
            F[1,l]=t[j+1]
            l=l+2
            break}}}}}

for (j in seq(N-s+1, N-1)){
    if (sqrt((X[j]-X[P[1]])**2+(Y[j]-Y[P[1]])**2) <= R & sqrt((X[j+1]-X[P[1]])**2+(Y[j+1]-Y[P[1]])**2) >= R){
      if (R - sqrt((X[j]-X[P[1]])**2+(Y[j]-Y[P[1]])**2) <= sqrt((X[j+1]-X[P[1]])**2+(Y[j+1]-Y[P[1]])**2) - R){
        F[1,k]=t[j]
        k=k+2
      }else{
        F[1,k]=t[j+1]
        k=k+2
      }}
  if (sqrt((X[j]-X[P[1]])**2+(Y[j]-Y[P[1]])**2) >= R & sqrt((X[j+1]-X[P[1]])**2+(Y[j+1]-Y[P[1]])**2) <= R){
    if (sqrt((X[j]-X[P[1]])**2+(Y[j]-Y[P[1]])**2) - R <= R - sqrt((X[j+1]-X[P[1]])**2+(Y[j+1]-Y[P[1]])**2)){
      F[1,l]=t[j]
      l=l+2
    }else{
      F[1,l]=t[j+1]
      l=l+2
    }}}
if (sqrt((X[P[1]]-X[N])**2 + (Y[P[1]]-Y[N])**2) <= R){
  F[1,k] = t[N]}}else{
    F[1,1] = t[N]
  }

B[1,1] = t[1]

##run through the rest of the algorithm in two parts
for( q in seq(2,N)){
  if (P[q] != 0){
    l_b=1
    k_b=2
    l_f=2
    k_f=1
    #Backward passages
    if (P[q]%%s == 0){
      max_lim = 2*s
    }else{
      max_lim = P[q]%%s+s
    }
    if (P[q] >= max_lim){
      for (i in seq(P[q],max_lim,-s)){
        #exits
        if (sqrt((X[i-s]-X[P[q]])**2+(Y[i-s]-Y[(P[q])])**2) <= R & sqrt((X[i]-X[P[q]])**2+(Y[i]-Y[P[q]])**2) >= R){
          for (j in seq(i,i-s-1,-1)){
            if (sqrt((X[j-1]-X[P[q]])**2+(Y[j-1]-Y[P[q]])**2) <= R & sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) >= R){
              if (R - sqrt((X[j-1]-X[P[q]])**2+(Y[j-1]-Y[P[q]])**2) <= sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) - R){
                B[q,k_b]=t[j-1]
                k_b=k_b+2
                break}else{
                  B[q,k_b]=t[j]
                  k_b=k_b+2
                  break}}}}
        #entries
        if (sqrt((X[i-s]-X[P[q]])**2+(Y[i-s]-Y[P[q]])**2) >= R & sqrt((X[i]-X[P[q]])**2+(Y[i]-Y[P[q]])**2) <= R){
          for (j in seq(i,i-s-1,-1)){
            if (sqrt((X[j-1]-X[P[q]])**2+(Y[j-1]-Y[P[q]])**2) >= R & sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) <= R){
              if (sqrt((X[j-1]-X[P[q]])**2+(Y[j-1]-Y[P[q]])**2) - R <= R - sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2)){
                B[q,l_b]=t[j-1]
                l_b=l_b+2
                break}else{
                  B[q,l_b]=t[j]
                  l_b=l_b+2
                  break}}}}}}
    #Solve problems at edges
    if (max_lim-s >= 2){
      for (j in seq(max_lim-s,2,-1)){
        if (sqrt((X[j-1]-X[P[q]])**2+(Y[j-1]-Y[P[q]])**2) <= R & sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) >= R){
          if (R - sqrt((X[j-1]-X[P[q]])**2+(Y[j-1]-Y[P[q]])**2) <= sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) - R){
            B[q,k_b]=t[j-1]
            k_b=k_b+2
          }else{
            B[q,k_b]=t[j]
            k_b=k_b+2
          }}
        if (sqrt((X[j-1]-X[P[q]])**2+(Y[j-1]-Y[P[q]])**2) >= R & sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) <= R){
          if (sqrt((X[j-1]-X[P[q]])**2+(Y[j-1]-Y[P[q]])**2) - R <= R - sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2)){
            B[q,l_b]=t[j-1]
            l_b=l_b+2
          }else{
            B[q,l_b]=t[j]
            l_b=l_b+2
          }}}}
    if (sqrt((X[P[q]]-X[1])**2 + (Y[P[q]]-Y[1])**2) <= R){
      B[q,l_b] = t[1]}

    #Forward passages
    if (P[q]+1 < N-s){
      for (i in seq(P[q],P[q]+floor((N-P[q]-s)/s)*s,s)){
        #exits
        if (sqrt((X[i]-X[P[q]])**2+(Y[i]-Y[P[q]])**2) <= R & sqrt((X[i+s]-X[P[q]])**2+(Y[i+s]-Y[P[q]])**2) >= R){
          for (j in seq(i,i+s+1)){
            if (sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) <= R & sqrt((X[j+1]-X[P[q]])**2+(Y[j+1]-Y[P[q]])**2) >= R){
              if (R - sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) <= sqrt((X[j+1]-X[P[q]])**2+(Y[j+1]-Y[P[q]])**2) - R){
                F[q,k_f]=t[j]
                if (k_f == 1){
                  P[q+1]=j}
                k_f=k_f+2
                break}else{
                  F[q,k_f]=t[j+1]
                  if (k_f == 1){
                    P[q+1]=j+1}
                  k_f=k_f+2
                  break}}}}
        #entries
        if (sqrt((X[i]-X[P[q]])**2+(Y[i]-Y[P[q]])**2) >= R & sqrt((X[i+s]-X[P[q]])**2+(Y[i+s]-Y[P[q]])**2) <= R){
          for (j in seq(i,i+s+1)){
            if (sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) >= R & sqrt((X[j+1]-X[P[q]])**2+(Y[j+1]-Y[P[q]])**2) <= R){
              if (sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) - R <= R - sqrt((X[j+1]-X[P[q]])**2+(Y[j+1]-Y[P[q]])**2)){
                F[q,l_f]=t[j]
                l_f=l_f+2
                break}else{
                  F[q,l_f]=t[j+1]
                  l_f=l_f+2
                  break}}}}}}
      #Solve problems at edges
      if (P[q]+floor((N-P[q]-s)/s)*s+s <= N-1){
        for (j in seq(P[q]+floor((N-P[q]-s)/s)*s+s, N-1)){
          if (sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) <= R & sqrt((X[j+1]-X[P[q]])**2+(Y[j+1]-Y[P[q]])**2) >= R){
            if (R - sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) <= sqrt((X[j+1]-X[P[q]])**2+(Y[j+1]-Y[P[q]])**2) - R){
              F[q,k_f]=t[j]
              k_f=k_f+2
            }else{
              F[q,k_f]=t[j+1]
              k_f=k_f+2
            }}
        if (sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) >= R & sqrt((X[j+1]-X[P[q]])**2+(Y[j+1]-Y[P[q]])**2) <= R){
          if (sqrt((X[j]-X[P[q]])**2+(Y[j]-Y[P[q]])**2) - R <= R - sqrt((X[j+1]-X[P[q]])**2+(Y[j+1]-Y[P[q]])**2)){
            F[q,l_f]=t[j]
            l_f=l_f+2
          }else{
            F[q,l_f]=t[j+1]
            l_f=l_f+2
          }}}}
    if (sqrt((X[P[q]]-X[N])**2 + (Y[P[q]]-Y[N])**2) <= R){
      F[q,k_f] = t[N]}}}

#Here I calculate the number of centers
C=0
for (i in seq(1,N)){
  if (P[i] != 0){
    C=C+1}}

#F[C,k_f]=t[N]

##Time spent inside the circle
psi = integer(C)##residence time
##We will calculate the residence time
for (i in seq(1,C)){
  f=0
  b=0
  for (j in seq(1,floor(m/2)-1)){
    f=f+F[i,2*j+1]-F[i,2*j]}
  for (k in seq(1,floor(m/2)-1)){
    b=b+B[i,2*k]-B[i,2*k+1]}
  psi[i]=F[i,1]-B[i,1]+f+b}

##Number of visits
zeta = integer(C)##number of visits
phi = integer(C)##number of forward visits
beta = integer(C)##number of backward visits

for (i in seq(1,C)){
  for (j in seq(1,floor(m/2))){
    if (F[i,2*j] != 0){
      if (F[i,2*j]-F[i,2*j-1] >= T_c){
        phi[i] = phi[i] + 1}}}}

for (i in seq(1,C)){
  for (j in seq(1,floor(m/2))){
    if (B[i,2*j] != 0){
      if (B[i,2*j-1]-B[i,2*j] >= T_c){
        beta[i] = beta[i] +1}}}}

zeta = beta + phi + 1

m1=2*max(phi)+1
m2=2*max(beta)+1

F2 = matrix(0,C,m1)
B2 = matrix(0,C,m2)

##Export the data##
a = matrix(0,C,7)
for (i in seq(1,C)){
  a[i,1]=t[P[i]]
  a[i,2]=X[P[i]]
  a[i,3]=Y[P[i]]
  a[i,4]=beta[i]
  a[i,5]=phi[i]
  a[i,6]=zeta[i]
  a[i,7]=psi[i]}
for (i in seq(1,C)){
  for (j in seq(1,m1)){
    F2[i,j]=F[i,j]}
  for (j in seq(1,m2)){
    B2[i,j]=B[i,j]}}

if (save == 'y'){
  write.table(a, file = paste(Name,"_UD_alt_R",R,".csv",sep=''), sep=",", row.names=FALSE, col.names=c("t","X","Y","beta","phi","zeta","psi"))
  write.table(F2, file = paste(Name,"_F_alt_R",R,".csv",sep=''), sep=",", row.names=FALSE, col.names=seq(1,m1))
  write.table(B2, file = paste(Name,"_B_alt_R",R,".csv",sep=''), sep=",", row.names=FALSE, col.names=seq(1,m2))
} 

t_centers=integer(C)
X_centers=integer(C)
Y_centers=integer(C)
for (j in seq(1,C)){
  t_centers[j]=t[P[j]]
  X_centers[j]=X[P[j]]
  Y_centers[j]=Y[P[j]]}

####Site Fidelity####
N_centers=length(t_centers)
if (N_centers!= 1){
  M=cbind(round(t_centers,4),X_centers,Y_centers,beta,phi,zeta,round(psi,4))
  M_sort = M[order(M[,7],M[,1]),]
  t_no_overlap = M_sort[,1]
  X_no_overlap = M_sort[,2]
  Y_no_overlap = M_sort[,3]
  beta_no_overlap = M_sort[,4]
  phi_no_overlap = M_sort[,5]
  zeta_no_overlap = M_sort[,6]
  psi_sort_no_overlap = M_sort[,7]
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

psi_sort_no_overlap2 = sort(psi_sort_no_overlap, decreasing = TRUE)
N_no_overlap = length(psi_sort_no_overlap2)
if (N_no_overlap == 1){
  max_percent_drop = 1
  number_identified_sites=1}else{
    Percent_drop = integer(N_no_overlap-1)
    for (i in seq(1,N_no_overlap-1)){
      Percent_drop[i] = 1 -(psi_sort_no_overlap2[i+1]/psi_sort_no_overlap2[i])
      max_percent_drop = max(Percent_drop)
      number_identified_sites=match(max_percent_drop,Percent_drop)}}
if (first == 'y'){
  if (number_identified_sites == N_no_overlap - 1 & t_no_overlap[1] == t[1]){
    Percent_drop = integer(N_no_overlap-2)
    for (i in seq(1,N_no_overlap-2)){
      Percent_drop[i] = 1 -(psi_sort_no_overlap2[i+1]/psi_sort_no_overlap2[i])
      max_percent_drop = max(Percent_drop)
      number_identified_sites=match(max_percent_drop,Percent_drop)}
  }
}
return(c(max_percent_drop,
         number_identified_sites))}


