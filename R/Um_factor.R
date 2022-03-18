## This file is part of RMaCzek

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


## Krzysztof Bartoszek documentation and made non-hidden
#'@title Calculate the Um factor
#'@description The function calculates the Um factor associated with an ordering of the rows and columns of a distance matrix. Lower values indicate a better grouping of similar objects. This was the original objective function proposed in the MaCzek program for producing Czekanowski's Diagram.
#'@param distMatrix a 'dist' object, matrix of distances between observations.
#'@param order a vector, if NULL, then the value of the factor is calculate for the distance matrix as is, otherwise the rows and columns are reordered according to the vector order.
#'@param matrix_conversion_coefficient numeric, value to be added to the distances, so that a division by 0 error is not thrown.
#'@param inverse_um logical, if TRUE, then the negative is returned. Default TRUE as the function is called in the genetic algorithm maximization procedures.
#'@export
#'@return The function returns a numeric value equalling the Um_factor.
#'@examples
#'# Set data ####
#'x<-mtcars
#'
#'mD<-stats::dist(scale(x))
#'mCz<-czek_matrix(x)
#'Um_factor(mD)
#'Um_factor(mD,order=attr(mCz,"order"))
## ========================================================================

Um_factor<-function(distMatrix,
                order=NULL,
                matrix_conversion_coefficient=1,
                inverse_um=TRUE){
##                ,scale_factor=FALSE){ ## Krzysztof Bartoszek uncommented

  # Change the class to matrix
##  if(class(distMatrix)!="matrix"){
  if(!is.matrix(distMatrix)){
    distMatrix<-as.matrix(distMatrix)
  }

  if(!is.null(order)){
    distMatrix<-distMatrix[order,order]
  }



  #Dimensions of the distance matrix
  d <- dim(distMatrix)
## Krzysztof Bartoszek ====================================================
  n<-d[1]

#  mDtmp<-matrix(1:n^2,n,n)
# mDtmp<-apply(mDtmp,c(1,2),function(ij){
#    i<-(ij-1)%/%n+1
#    j<-(ij-1)%%n+1
#    (i-j)^2
# })
# 
#  mDtmp<-mDtmp/(distMatrix+1)
  mDtmp<-(outer(1:n,1:n,function(x,y){(x-y)^2}))/(distMatrix+1)
  um_value<-2*sum(mDtmp[upper.tri(mDtmp,diag=FALSE)])/n^2
## ========================================================================  

  if(inverse_um){
    um_value<-(-um_value)
  }

## Krzysztof Bartosek made commented
##  if(scale_factor=="n^4"){
##    um_value<-um_value/nrow(distMatrix)^2
##  }
## ========================================================================

  return(um_value)

}
