#'@title Preprocess data to produce Czekanowski's Diagram.
#'@description This is a function that divided the values inside a distance matrix into classes. The output can be used in the plot function to produce a Czekanowski's Diagram.
#'@param x  a numeric matrix, data frame or a 'dist' object.
#'@param order specifies which seriation method should be applied. The standard setting is the seriation method OLO. If NA or NULL, then no seriation is done and the original ordering is saved.  The user may provide their own ordering, through a number vector of indices. Also in this case no rearrangement will be done.
#'@param n_classes specifies how many classes the distances should be divided into. The standard setting is 5 classes.
#'@param interval_breaks specifies the partition boundaries for the distances. As a standard setting, each class represents an equal amount of distances. If the interval, breaks are positive and sum up to 1, then it is assumed that they specify percantages of the distances in each interval. Otherwise if provided as a numeric vector not summing up to 1, they specify the eact boundaries for the symbols representing distance groups.
#'@param monitor specifies if the distribution of the distances should be visualized. The standard setting is that the distribution will not be visualized. TRUE and "cumulativ_plot" is available.
#'@param distfun specifies which distance function should be used. Standard setting is the dist function which uses the Euclidean distance.
#'@param scale_data specifies if the data set should be scaled. The standard setting is that the data will be scaled.
#'@param focal_obj Numbers or names of objects (rows if x is a dataset and not 'dist' object) that are not to take part in the reordering procedure. These observations will be placed as last rows and columns of the output matrix. See Details.
#'@param as_dist If TRUE, then the distance matrix of x is returned, with object ordering, instead of the matrix with the levels assigned in place of the original distances. 
#'@param original_diagram If TRUE, then the returned matrix corresponds as close as possible to the original method proposed by Czekanowski (1909). The levels are column specific and not matrix specific. See Details
#'@param column_order_stat_grouping If original_diagram is TRUE, then here one can pass the partition boundaries for the ranking in each column.
#'@param ... specifies further parameters that can be passed on to the seriate function in the seriation package.
#'@export
#'@return The function returns a matrix with class czek_matrix. The return from the function is expected to be passed to the plot function. If as_dist is passed as TRUE, then a czek_matrix_dist object is returned and this is not suitable for the plotting. As an attribute of the output the optimized criterion value is returned. However, this is a guess based on seriation::seriate()'s and seriation::criterion()'s manuals. If something else was optimized, e.g. due to user's parameters, then this will be wrong. If unable to guess, then NA saved in the attribute.
#'@examples
#'# Set data ####
#'x<-mtcars
#'
#'
#'# Different type of input that give same result ############
#'czek_matrix(x)
#'czek_matrix(stats::dist(scale(x)))
#'
#'
#'# Change seriation method ############
#'#seriation::show_seriation_methods("dist")
#'czek_matrix(x,order = "GW")
#'czek_matrix(x,order = "ga")
#'czek_matrix(x,order = sample(1:nrow(x)))
#'
#'
#'# Change number of classes ############
#'czek_matrix(x,n_classes = 3)
#'
#'
#'# Change the partition boundaries ############
#'czek_matrix(x,interval_breaks = c(0.1,0.4,0.5)) #10%, 40% and 50%
#'czek_matrix(x,interval_breaks = c(0,1,4,6,8.48)) #[0,1] (1,4] (4,6] (6,8.48]
#'czek_matrix(x,interval_breaks = "equal_width_between_classes") #[0,1.7] (1.7,3.39]  (3.39,5.09] (5.09,6.78] (6.78,8.48]
#'
#'
#'# Change number of classes ############
#'czek_matrix(x,monitor = TRUE)
#'czek_matrix(x,monitor = "cumulativ_plot")
#'
#'
#'# Change distance function ############
#'czek_matrix(x,distfun = function(x) stats::dist(x,method = "manhattan"))
#'
#'
#'# Change dont scale the data ############
#'czek_matrix(x,scale_data = FALSE)
#'czek_matrix(stats::dist(x))
#'
#'
#'# Change additinal settings to the seriation method ############
#'czek_matrix(x,order="ga",control=list(popSize=200,
#'                                      suggestions=c("SPIN_STS","QAP_2SUM")))
#'
#'# Create matrix as originally described by Czekanowski (1909), with each column
#'# assigned levels according to how the order statistics of the  distances in it
#'# are grouped. The grouping below is the one used by Czekanowski (1909).
#'czek_matrix(x,original_diagram=TRUE,column_order_stat_grouping=c(3,4,5,6))
#'
#'# Create matrix with two focal object that will not influence seriation
#'czek_matrix(x,focal_obj=c("Merc 280","Merc 450SL"))
#'# Same results but with object indices
#'czek_res<-czek_matrix(x,focal_obj=c(10,13))
#'# we now place the two objects in a new place
#'attr(czek_res,"order")<-attr(czek_res,"order")[c(1:10,31,11:20,32,21:30)]
#'# and then correct the values of the different criteria so that they
#'# are consistent with the new ordering
#'attr(czek_res,"Path_length")<-seriation::criterion(stats::dist(scale(x)),order=seriation::ser_permutation(attr(czek_res, "order")),method="Path_length")
#'# Here we need to know what criterion was used for the seriation procedure
#'# If the seriation package was used, then see the manual for seriation::seriate()
#'# seriation::criterion().
#'# If the genetic algorithm shipped with RMaCzek was used, then it was the Um factor.
#'attr(czek_res,"criterion_value")<-seriation::criterion(stats::dist(scale(x)),order=seriation::ser_permutation(attr(czek_res, "order")),method="Path_length")
#'attr(czek_res,"Um")<-RMaCzek::Um_factor(stats::dist(scale(x)),order= attr(czek_res, "order"),inverse_um=FALSE)
#'@details
#'In his original paper Czekanowski (1909) did not have as the output a symmetric matrix
#'where each distance was assigned a level (symbol) depending in which numeric interval it was in.
#'Instead having the desired ordering, the following procedure was applied to each column.
#'The three smallest distances (in each column) obtain level (symbol) 1, the fourth smallest
#'level (symbol) 2, fifth smallest level (symbol) 3, sixth smallest level (symbol) 4 and
#'all the bigger distances the fifth symbol which was originally just a blank cell in the output
#'matrix. Here, we give the user more flexibility. In column_order_stat_grouping one may specify
#'how the order statistics should be grouped in each column. See Example.
#'The user may also choose some observations not to influence the ordering procedure.
#'This could be useful if e.g. a single observation is meant to be assigned to a cluster
#'and for some reason the clusters (that are to be read of from the ordering) should 
#'not be influenced by this observatin. One can pass such observations using the
#'focal_obs parameter.
#'A hopefully useful property is that the ordering inside the czek_matrix (and hence
#'of the diagram when one calls plot) can be manually changed. One merely manipulates
#'the order attribute as desired. However in such a case one should remember that
#'the Path_length, criterion_value and Um attributes will have incorrect values
#'and should be corrected (see Examples).

czek_matrix <- function(x,
                        order="OLO",
                        n_classes = 5,
                        interval_breaks=NULL,
                        monitor=FALSE,
                        distfun=dist,
                        scale_data=TRUE,
                        focal_obj=NULL,
                        as_dist=FALSE,
                        original_diagram=FALSE,
                        column_order_stat_grouping=NULL,
                        ...){

  # If not of class dist, make the data to class dist ####
##  if(class(x)!="dist"){
  if(all(class(x)!="dist")){

    # Scale data
    if(scale_data){
      x<-scale(x)
    }

    # Calculate a distance matrix
    x<-distfun(x)

  }
  
  ## Krzysztof Bartoszek
  if (!is.null(focal_obj)){
  ## by now x is a distance matrix
    x<-as.matrix(x)
    if (is.character(focal_obj)){
	focal_obj<-which(is.element(colnames(x),focal_obj))	
    }
    if ((length(focal_obj)>0)&&(is.numeric(focal_obj))&&(min(focal_obj)>0)&&(max(focal_obj)<(ncol(x)+1))){
	##mfocals<-x[,focal_obj,drop=FALSE]
	if (length(unique(focal_obj))==ncol(x)){
	    order<-NA
	    v_order_nofocal<-NA
	    x_org<-NA
	    focal_obj<-NULL
	    warning("Parameter focal_obj spans all the objects! Treating as if no ordering was requested!")
	}else{
	    if (length(focal_obj)!=length(unique(focal_obj))){
		focal_obj<-unique(focal_obj)
		warning("Repeated columns in focal_obj parameter, removing them!")
	    }
	    v_order_nofocal<-(1:ncol(x))[-focal_obj]
	    x_org<-x
	    x<-x[-focal_obj,-focal_obj,drop=FALSE]
	}
    }else{
	warning("There is some problem with the focal_obj parameter, ignoring it.")
	## mfocals<-NA
	v_order_nofocal<-NA
	x_org<-NA
	focal_obj<-NULL
    }    
    x<-as.dist(x)
  }
  ## =========================================
  res_seriate<-NA    
  # Seriation part ####
  # If the user have specified the order
  if((all(is.numeric(order)))&&(length(order)==attr(x,"Size"))){ ##Krzysztof Bartoszek added any() and also checking for length
    #Just for conviction
    new_order<-order
    ##order<-NA
  }
  # If standard settings is used
  else if (any(class(order[1])=="character")){
##    if (!.installed("seriation")) stop("Package 'seriation' needs to be installed!")
    # If standard settings is used
    if (order[1]=="ga"){
      .register_seriate_ga()
      order<-".seriate_ga"
    }
    new_order<-seriation::get_order(seriation::seriate(x,method=order[1],...))
  }
  ## Krzysztof Bartoszek
  # If the user dont want to change the order
  else if ((is.na(order[1]))||(is.null(order[1]))) {
    new_order<-1:attr(x,"Size")
    order<-NA
  } 
  ## By default do not change the order
  ##==========================================
  else {
    new_order<-1:attr(x,"Size")
    order<-NA
  }

  if (!is.null(focal_obj)){
    x<-x_org
    new_order<-c(v_order_nofocal[new_order],focal_obj)
  }
  # Change the class to matrix ####
  x<-as.matrix(x)

  if (!as_dist){ ## Krzysztof Bartoszek
    if (!original_diagram){ ## Krzysztof Bartoszek
	# Find the partition bounderies ####
	# If ther dont have specified the interval breaks
	if(is.null(interval_breaks)) {
	    # If NOT the user have specified the intervals
	    # Given 5 classes: 20% of the distance in class 1,..., 20% of the distance in class 5
	    interval_breaks<-stats::quantile(x[upper.tri(x)], probs=seq(0,1,length.out = n_classes+1), na.rm=TRUE)
	    interval_breaks[1]<-0
	 }
	 # If the user want equal width brettwen classes
	 else if("equal_width_between_classes"%in%interval_breaks){
	    interval_breaks<-max(x[upper.tri(x)])/n_classes*(0:n_classes)
	    probs<-stats::ecdf(x[upper.tri(x)])(interval_breaks)
	    names(interval_breaks)<-paste(round(probs,7)*100,"%",sep="")
	 }
	 # If interval_breaks is specified in procent
	 else if ((all(interval_breaks>=0))&&(sum(interval_breaks)==1)){ ##Krzysztof Bartoszek condition for all positive
	    probs<-c(0,cumsum(interval_breaks))
	    interval_breaks<-stats::quantile(x[upper.tri(x)], probs=probs, na.rm=TRUE)
	    interval_breaks[1]<-0
	}
	#If the user have specified the intervals
	else{
	    interval_breaks[1]<-0
	    interval_breaks[length(interval_breaks)]<-max(x)

	    probs<-ecdf(x[upper.tri(x)])(interval_breaks)
	    names(interval_breaks)<-paste(round(probs,7)*100,"%",sep="")
	}
	# Split the distances into classes ####

	# Make the partition boundaries
	cut_the_values <- cut(x, interval_breaks, include.lowest = TRUE)

	# Make the matrix that we want to plot later on
	czek_matrix<- matrix(as.numeric(cut_the_values),ncol = ncol(x))
	l_cut_the_values<-levels(cut_the_values) ## Krzysztof Bartoszek
    }else{
	### Krzysztof Bartoszek
	## Czekanowski's original diagram with grouping ranks in each column.
	b_default_ord_grouping<-FALSE
	if (is.null(column_order_stat_grouping)){
	    column_order_stat_grouping<-c(3,4,5,6) ## Czekanowski (1909)'s original order statistiscs grouping
	    b_default_ord_grouping<-TRUE
	}
	if (max(column_order_stat_grouping) > (ncol(x)-1)){
	    if (!b_default_ord_grouping){warning("Provided or default order statistics grouping is not compatible with matrix size. Correcting but please consider providing a compatible one.")}
	    v_OKgroups<-which(column_order_stat_grouping<ncol(x))
	    nOKgroups<-length(v_OKgroups)
	    if (nOKgroups==0){
		if (ncol(x)==0){stop("")}
		else if (ncol(x)==1){column_order_stat_grouping<-1}
		else if (ncol(x)==2){column_order_stat_grouping<-1}
		else if (ncol(x)==3){column_order_stat_grouping<-2}
		else if (ncol(x)==4){column_order_stat_grouping<-c(2,3)}
		else if (ncol(x)==5){column_order_stat_grouping<-c(2,3,4)}
		else if (ncol(x)==6){column_order_stat_grouping<-c(2,3,4,5)}
		else if (ncol(x)>7){column_order_stat_grouping<-c(3,4,5,6)}	
		else{stop("Something is wrong with provided order statistics groupings!")}
	    }
	    else if ((nOKgroups>0)&&(nOKgroups<5)){
		column_order_stat_grouping<-column_order_stat_grouping[v_OKgroups]
		while ((max(column_order_stat_grouping)<ncol(x)-1)&&(nOKgroups<6)){
		    column_order_stat_grouping<-c(column_order_stat_grouping,max(column_order_stat_grouping)+1)
		    nOKgroups<-nOKgroups+1
		}
	    }	    
	    else if (nOKgroups>4){
		column_order_stat_grouping<-column_order_stat_grouping[v_OKgroups]
	    }else{
		stop("Something is wrong with provided order statistics groupings!")
	    }
	}
	if ((length(column_order_stat_grouping)!=length(unique(column_order_stat_grouping)))||(!all(column_order_stat_grouping==cummax(column_order_stat_grouping)))){
	    column_order_stat_grouping<-unique(column_order_stat_grouping)
	    column_order_stat_grouping<-sort(column_order_stat_grouping)
	    warning("Parameter column_order_stat_grouping was not strictly increasing. Correcting but please consider supplying a strictly increasing one.")
	}
	czek_matrix<-apply(x,2,function(x_col,column_order_stat_grouping){
	    cut_the_values<-cut(base::rank(x_col),c(0,column_order_stat_grouping,length(x_col)),include_lowest=TRUE)
	    czek_x_col<-as.numeric(cut_the_values)
	    czek_x_col
	},column_order_stat_grouping=column_order_stat_grouping)
	cut_the_values<-cut(1:ncol(x),c(0,column_order_stat_grouping,ncol(x)),include_lowest=TRUE)
	tmp_cut_the_values<-cut_the_values
	l_cut_the_values<-levels(cut_the_values)
	l_cut_the_values[1]<-paste0("[1,",column_order_stat_grouping[1],"]")
	interval_breaks<-"Column specific, levels are the cut values of the order statistics for each column."
	### =================================================
    }
    # attr information to the matrix with classes ####
    # Add the partition boundaries to the matrix
    attr(czek_matrix, "levels") <- l_cut_the_values ## Krzysztof Bartoszek
    attr(czek_matrix, "partition_boundaries")<-interval_breaks
    attr(czek_matrix, "n_classes")<-length(levels(czek_matrix))
  }else{
    czek_matrix<-x
  }
  attr(czek_matrix, "order")<-new_order
  ## Krzysztof Bartoszek 
  ## seriation::criterion takes order directly for order parameter NOT rank
  attr(czek_matrix, "criterion_value")<-NA
  if (!is.na(order[1])){
    criterion_method<-NA
    ## Attempt as best as good based on ?seriation::seriate and ?seriation::criterion to guess the optimized criterion for every method
    criterion_method<-switch(as.character(order[1]),
	    ARSA="LS",BBURCG="Gradient_raw", BBWRCG="Gradient_weighted",
	    GW="Path_length",GW_average="Path_length",GW_complete="Path_length",GW_single="Path_length", GW_ward="Path_length",
	    HC=NA,HC_average=NA,HC_complete=NA,HC_single=NA,HC_ward=NA,Identity=NA,
	    MDS="Neumann_stress",MDS_angle="Neumann_stress",MDS_metric="Neumann_stress",MDS_nonmetric="Neumann_stress",
	    OLO="Path_length",OLO_average="Path_length",OLO_complete="Path_length",OLO_single="Path_length",OLO_ward="Path_length",
	    QAP_2SUM="2SUM",QAP_BAR="BAR",QAP_Inertia="Inertia",QAP_LS="LS",
	    R2E=NA,Random=NA,SA="Gradient_raw",Spectral="2SUM",Spectral_norm="2SUM",SPIN_NH=NA,SPIN_STS=NA,TSP="Path_length",VAT=NA,NA)
    if (!is.na(criterion_method)){
	attr(czek_matrix, "criterion_value")<-seriation::criterion(as.dist(x),order=seriation::ser_permutation(new_order),method=criterion_method)
    }
  }
  ## =================================================
  attr(czek_matrix, "Path_length")<-seriation::criterion(as.dist(x),order=seriation::ser_permutation(new_order),method="Path_length")
  attr(czek_matrix, "Um")<-Um_factor(x,order= attr(czek_matrix, "order"),inverse_um=FALSE) ## Krzysztof Bartoszek

  # Add row/col names to the color matrix
  rownames(czek_matrix)<-rownames(x)
  colnames(czek_matrix)<-colnames(x)

  if ((!as_dist)&&(!original_diagram)){ ## Krzysztof Bartoszek
	# Monitor ####
	if(monitor%in%c(TRUE,"cumulativ_plot")){

	    if(monitor==TRUE)
    		monitor<-"plot"

	    cum_probs<-as.numeric(gsub(pattern = "%",replacement = "",x = names(interval_breaks)))
	    plot_values<-cum_probs[-1]
	    my_title<-"Cumulative distribution of distances in each class"

	    if(monitor=="plot"){
    	    probs<-c()
    	    for(i in 2:(length(cum_probs))){
    		probs[i-1]<-cum_probs[i]-cum_probs[i-1]
    	    }
    	    plot_values<-probs
    	    my_title<-"Distribution of distances in classes"
	    }

	    names(plot_values)<-levels(cut_the_values)
	    graphics::barplot(plot_values,
        	main=my_title,
        	col=c("grey30"),
        	xlab = "Classes of distances",
        	ylim = c(0,100),yaxt="n")

	    graphics::axis(2, at = seq(0, 100, by = 10),
        	labels=paste(seq(0, 100, by = 10),"%",sep=""),
        	las=2)

	    graphics::box(col = "black")
	}
     # Change class and return ####
    # Change class
	class(czek_matrix)<-"czek_matrix"
    }else{
  ## Krzysztof Bartoszek
  ## the user wants a distance matrix to be returned
	if (as_dist){class(czek_matrix)<-"czek_matrix_dist"}
	else{class(czek_matrix)<-"czek_matrix"}
  ## ========================================================
  }
  # Return
  return(czek_matrix)
}

