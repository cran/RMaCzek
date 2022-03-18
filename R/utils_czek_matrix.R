## This file is part of RMaCzek

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


#'@title Prints information concerning Czekanowski's Diagram
#'@description This is a function that prints out information on a Czekanowski's Diagram.
#'@param x  a matrix with class czek_matrix.
#'@param print_raw logical, if TRUE print out raw, as if base::print() was called, in particular this prints out the matrix itself, if FALSE (default) print out a summary. Furthermore, with print_raw=TRUE the attributes "levels", "partition_boundaries" and "n_classes" defining the diagram will be printed out.
#'@param ... specifies further parameters that can be passed on to the print function.
#'@export
#'@return The function returns a Czekanowski's Diagram.
#'@examples
#'# Set data ####
#'x<-czek_matrix(mtcars)
#'
#'
#'# Standard print ############
#'print(x)
#'print.czek_matrix(x)
#'# Print out the raw object ############
#'print(x,print_raw=TRUE)
#'print.czek_matrix(x,print_raw=TRUE)

print.czek_matrix<-function(x,
			    print_raw=FALSE,
                           ...){
    if (print_raw){
    	print(.as_matrix_czek_matrix(x,b_remove_attributes=FALSE))
    }else{
	if (!is.null(rownames(x))){v_orderednames<-rownames(x)[attr(x,"order")]}
	else{
	    if (!is.null(colnames(x))){v_orderednames<-colnames(x)[attr(x,"order")]}
	    else{v_orderednames<-attr(x,"order")}
	}
	v_output<-paste(v_orderednames,collapse=", ")
	v_output<-paste0(v_output," ","\n")
	if (is.element("criterion_value",names(attributes(x)))){
	    crit_value<-attr(x,"criterion_value")
	    y<-crit_value[[1]]
	    z<-names(crit_value)
	    if ((!is.null(y))&&(!is.na(y))&&(!is.null(z))){
	    	v_output<-paste0(v_output,z," criterion value: ",y,"\n")
	    }
	}
	if (is.element("Um",names(attributes(x)))){
	    crit_value<-attr(x,"Um")
	    y<-crit_value[[1]]
	    if ((!is.null(y))&&(!is.na(y))){
	    	v_output<-paste0(v_output,"Um factor: ",y,"\n")
	    }
	}
	if (is.element("Path_length",names(attributes(x)))){
	    crit_value<-attr(x,"Path_length")
	    y<-crit_value[[1]]
	    if ((!is.null(y))&&(!is.na(y))){
	    	v_output<-paste0(v_output,"Hamiltonian path length factor: ",y,"\n")
	    }
	}	
	cat(v_output)
    }
}

#'@title Prints information concerning Czekanowski's Diagram
#'@description This is a function that prints out information on a Czekanowski's Diagram, but when the actual distances were saved.
#'@param x  a matrix with class czek_matrix_dist.
#'@param print_raw logical, if TRUE print out raw, as if base::print() was called, in particular this prints out the matrix itself, if FALSE (default) print out a summary. Furthermore, with print_raw=TRUE the attributes "levels", "partition_boundaries" and "n_classes" defining the diagram will be printed out.
#'@param ... specifies further parameters that can be passed on to the print function.
#'@export
#'@return The function returns a Czekanowski's Diagram.
#'@examples
#'# Set data ####
#'x<-czek_matrix(mtcars,as_dist=TRUE)
#'
#'
#'# Standard print ############
#'print(x)
#'print.czek_matrix(x)
#'# Print out the raw object ############
#'print(x,print_raw=TRUE)
#'print.czek_matrix(x,print_raw=TRUE)

print.czek_matrix_dist<-function(x,
			    print_raw=FALSE,
                           ...){
    print.czek_matrix(x,print_raw=print_raw)
}


#'@title Manually reorder  Czekanowski's Diagram
#'@description This is a function that allows the user to manully reorder Czekanowski's Diagram and recalculates all the factors.
#'@param x  a matrix with class czek_matrix, czek_matrix_dist or data matrix/data.frame or dist object.
#'@param v_neworder a numeric vector with the new ordering.
#'@param ... specifies further parameters that will be passed to the czek_matrix function.
#'@export
#'@return The function returns a Czekanowski's Diagram with the new order and recalculated factors.
#'@examples
#'# Set data ####
#' x<-mtcars
#'
#'# Calculate Czekanowski's diagram
#' czkm<-czek_matrix(x)
#' czkm_dist<-czek_matrix(x,as_dist=TRUE)
#'# new ordering
#' neworder<-attr(czkm,"order")
#' neworder[1:2]<-neworder[2:1]
#'# reorder the diagram
#' #if the output was Czekanowski's diagram without the distances
#' #remembered, then the original data has to be passed so that 
#' #factors can be recalculated.
#' new_czkm<-manual_reorder(czkm,v_neworder=neworder,orig_data=x)
#' new_czkm_dist<-manual_reorder(czkm_dist,v_neworder=neworder)
#' #we can also pass the original data directly
#' new_czkm<-manual_reorder(x,v_neworder=neworder)
#' #and this is equivalent to calling
#' czkm<-czek_matrix(x,order=neworder)
#' #up to the value of the "criterion_value" attribute
#' #which in the second case can be lost, as no information is passed
#' #on which one was originally used, while in the first case it might
#' #be impossible to recalculate-only criteria values from seriate are supported
#' #if a user has a custom seriation function, then they need to recalculate this
#' #value themselves

manual_reorder <- function(x, v_neworder, ...){
  UseMethod("manual_reorder")
}

manual_reorder.czek_matrix<-function(x,v_neworder,orig_data,...){
    x_dist<-czek_matrix(orig_data,order="Identity",as_dist=TRUE,...)
    x_dist<-manual_reorder(x_dist,v_neworder)
    attr(x,"Um")<-attr(x_dist,"Um")
    attr(x,"Path_length")<-attr(x_dist,"Path_length")
    attr(x,"criterion_value")<-attr(x_dist,"criterion_value")
    attr(x,"order")<-v_neworder
    x
}

manual_reorder.czek_matrix_dist<-function(x,v_neworder,...){
    m<-.as_matrix_czek_matrix(x)
    #attr(m,"Um")<-NULL;attr(m,"Path_length")<-NULL;attr(m,"order")<-NULL;attr(m,"criterion_value")<-NULL;
    attr(x,"Path_length")<-seriation::criterion(as.dist(m),order=seriation::ser_permutation(v_neworder),method="Path_length")
    attr(x,"Um")<-Um_factor(m,order=v_neworder,inverse_um=FALSE)
    if (is.element("criterion_value",names(attributes(x)))){
    	crit_value<-attr(x,"criterion_value")
	y<-crit_value[[1]]
	z<-names(crit_value)
	if ((!is.null(y))&&(!is.na(y))&&(!is.null(z))){
	    new_crit_value<-NA
	    if (as.character(z) %in% c("LS","Gradient_raw","Gradient_weighted","Path_length","Neumann_stress","2SUM","BAR","Inertia",orm="2SUM")){
		new_crit_value<-seriation::criterion(as.dist(m),order=seriation::ser_permutation(v_neworder),method=z)
	    }
	    names(new_crit_value)<-z
	    attr(x,"criterion_value")<-new_crit_value
	}else{if(is.null(z)){attr(x,"criterion_value")<-NA}}
    }else{attr(x,"criterion_value")<-NA}
    attr(x,"order")<-v_neworder
    x
}

manual_reorder.data.frame<-function(x,v_neworder,...){
    czkm<-czek_matrix(x,order="Identity",...)
    manual_reorder(czkm,v_neworder,orig_data=x,...)    
}

manual_reorder.matrix<-function(x,v_neworder,...){
    czkm<-czek_matrix(x,order="Identity",...)
    manual_reorder(czkm,v_neworder,orig_data=x,...)    
}

manual_reorder.dist<-function(x,v_neworder,...){
    czkm<-czek_matrix(x,order="Identity",...)
    czkm_res<-NA
    if (inherits(czkm,"czek_matrix")){
        czkm_res<-manual_reorder(czkm,v_neworder,orig_data=x,...)
    }else if(inherits(czkm,"czek_matrix_dist")){
    	czkm_res<-manual_reorder(czkm,v_neworder,...)
    }
    czkm_res
}

.as_matrix_czek_matrix<-function(czkm,b_remove_attributes=TRUE){
    if (b_remove_attributes){
	attr(czkm,"Um")<-NULL;attr(czkm,"Path_length")<-NULL;attr(czkm,"order")<-NULL;attr(czkm,"criterion_value")<-NULL;
	attr(czkm,"levels")<-NULL;attr(czkm,"partition_boundaries")<-NULL;attr(czkm,"n_classes")<-NULL
    }
    czkm<-as.matrix(czkm)
    class(czkm)<-"matrix"
    czkm
}