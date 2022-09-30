## This file is part of RMaCzek

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.register_seriate_ga<-function(){
## KB edited for changes in seriation  
  if(!.my_check_is_ga_registered()){
 ## if(is.null(try(seriation::show_seriation_methods("dist")$dist_seriate_ga))){
##============================================    
    seriation::set_seriation_method(kind="dist",
                                    name=".seriate_ga",
                                    definition=.seriate_ga,
                                    description="Genetic Algorithm for permutation")

  }
}

## KB added for changes in seriation
.my_check_is_ga_registered<-function(){
    try(is.element(".seriate_ga",seriation::list_seriation_methods()$dist))
}
## ==============================================
