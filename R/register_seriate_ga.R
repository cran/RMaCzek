.register_seriate_ga<-function(){
  if(is.null(try(seriation::show_seriation_methods("dist")$dist_seriate_ga))){
    seriation::set_seriation_method(kind="dist",
                                    name=".seriate_ga",
                                    definition=.seriate_ga,
                                    description="Genetic Algorithm for permutation")

  }
}
