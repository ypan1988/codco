# DESIGN SPACE CLASS

DesignSpace <- R6::R6Class("DesignSpace",
                 public = list(
                   weights = NULL,
                   initialize = function(...,
                                         weights=NULL) {
                     samp.size <- c()
                     for (item in list(...)) {
                       if(!is(item,"Design"))stop(paste0(item," is not a Design"))
                     }
                     if(!is.null(weights)){
                       if(length(weights)!=length(list(...)))stop("weights not same length as designs")
                       self$weights <- weights
                     } else {
                       self$weights <- rep(1/length(list(...)),length(list(...)))
                     }
                     for (item in list(...)) {
                       self$add(item)
                     }

                   },
                   add = function(x) {
                     private$designs <- append(private$designs, list(x))
                     self$weights <- rep(1/length(private$designs),length(private$designs))
                     invisible(self)
                   },
                   remove = function(index) {
                     if (private$length() == 0) return(NULL)
                     private$designs <- private$designs[-index]
                   },
                   print = function(){
                     cat(paste0("Design space with ",self$n()," design(s): \n"))
                     for(i in 1:length(private$designs)){
                       cat(paste0("=========================================================\nDESIGN ",i,"(weight ",self$weights[i],"):\n"))
                       print(private$designs[[i]])
                     }
                   },
                   power = function(par,
                                    value,
                                    alpha=0.05){
                     pwr <- c()
                     if(self$n()>1 && length(par)==1)par <- rep(par,self$n())
                     if(self$n()>1 && length(value)==1)value <- rep(value,self$n())
                     for(i in 1:self$n()){
                       pwr[i] <- private$designs[[i]]$power(par[i],value[i],alpha)
                     }
                     return(data.frame(Design = 1:self$n(),power = pwr))
                   },
                   n = function(){
                     length(private$designs)
                   },
                   optimal = function(m,
                                      C,
                                      rm_cols=NULL,
                                      keep=TRUE,
                                      verbose=TRUE){
                     if(keep&verbose)message("linked design objects will be overwritten with the new design")
                     #initialise from random starting index
                     n <- private$designs[[1]]$mean_function$n()
                     idx_in <- sort(sample(1:n,m,replace=FALSE))
                     idx_out <- grad_robust_step(
                       idx_in,
                       C_list = list(matrix(C,ncol=1)),
                       X_list = private$genXlist(),
                       sig_list = private$genSlist(),
                       rm_cols = rm_cols,
                       trace=verbose
                     )
                     idx_out <- sort(idx_out)
                     if(keep){
                       for(i in 1:self$n()){
                         private$designs[[i]]$subset_rows(idx_out)
                         ncol <- 1:ncol(private$designs[[i]]$mean_function$X)
                         private$designs[[i]]$subset_cols(ncol[-rm_cols[[i]]])
                       }
                     }
                   },
                   show = function(i){
                     return(private$designs[[i]])
                   }
                 ),
                 private = list(
                   designs = list(),
                   genXlist = function(){
                     X_list <- list()
                     for(i in 1:self$n()){
                       X_list[[i]] <- as.matrix(private$designs[[i]]$mean_function$X)
                     }
                     return(X_list)
                   },
                   genSlist = function(){
                     S_list <- list()
                     for(i in 1:self$n()){
                       S_list[[i]] <- as.matrix(private$designs[[i]]$Sigma)
                     }
                     return(S_list)
                   }
                 )
)








