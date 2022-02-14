# MEAN FUNCTION CLASS
MeanFunction <- R6::R6Class("MeanFunction",
                        public = list(
                          formula = NULL,
                          data = NULL,
                          family = NULL,
                          parameters = NULL,
                          randomiser = NULL,
                          treat_var = NULL,
                          X = NULL,
                          n=function(){
                            nrow(self$data)
                          },
                          check = function(verbose=TRUE){
                            if(private$hash != private$hash_do()){
                              if(verbose)message("changes found, updating")
                              self$generate()
                            }},
                          initialize = function(formula,
                                                data,
                                                family,
                                                parameters ,
                                                verbose = FALSE,
                                                random_function=NULL,
                                                treat_par = NULL
                          ){
                            if(any(missing(formula),missing(family),missing(parameters))){
                              cat("not all inputs set. call generate() when set")
                            } else {
                              self$formula <- as.formula(formula, env=.GlobalEnv)
                              self$family <- family
                              self$parameters <- parameters

                              if(!is(data,"data.frame"))stop("data must be data frame")
                              # self$n <- nrow(data)
                              if(!is.null(random_function)){
                                if(is.null(treat_var)){
                                  stop("provide name of treatment variable treat_var")
                                }
                                #test random function
                                # check it produces a varia
                                self$randomise <- random_function
                              }
                              self$data <- data
                              self$generate(verbose=verbose)
                            }},
                          generate = function(verbose = FALSE){

                            if(length(self$formula)==3)stop("formula should not have dependent variable.")
                            #check if all parameters in data
                            if(any(!all.vars(self$formula)%in%colnames(self$data)))stop("variables not in data frame")

                            private$genTerms()
                            private$genX()
                            private$hash <- private$hash_do()
                          },
                          print = function(){
                            cat("Mean Function")
                            print(self$family)
                            cat("Formula:")
                            print(self$formula)
                            # cat("Data:\n")
                            # print(head(self$data))
                          },
                          colnames = function(names = NULL){
                            if(is.null(names)){
                              print(colnames(self$data))
                            } else {
                              colnames(self$data) <- names
                            }
                          },
                          subset_rows = function(index){
                            self$X <- self$X[index,]
                            self$data <- self$data[index,]
                          },
                          subset_cols = function(index){
                            self$X <- self$X[,index]
                          }
                        ),
                        private = list(
                          mod_string = NULL,
                          form = NULL,
                          Xb = NULL,
                          funs = NULL,
                          vars = NULL,
                          hash = NULL,
                          hash_do = function(){
                            digest::digest(c(self$formula,self$data,self$family,
                                             digest::digest(as.character(self$parameters),serialize = FALSE),
                                             self$randomiser))
                          },
                          genTerms = function(){
                            mf1 <- self$formula[[2]]
                            checkTerm <- TRUE
                            iter <- 0
                            funs <- list()
                            vars <- list()
                            while(checkTerm){
                              iter <- iter + 1
                              checkTerm <- I(length(mf1)>1 && (mf1[[1]]=="+"|mf1[[1]]=="-"))
                              if(checkTerm){
                                vars[[iter]] <- all.vars(mf1[[3]])
                                if(length(mf1[[3]])==1){
                                  funs[[iter]] <- "identity"
                                } else {
                                  funs[[iter]] <- as.character(mf1[[3]][[1]])
                                }
                                if(length(vars[[iter]])==0){
                                  vars[[iter]] <- funs[[iter]] <- "RMINT"
                                }
                                mf1 <- mf1[[2]]
                              } else {
                                vars[[iter]] <- all.vars(mf1)
                                if(length(mf1)==1){
                                  funs[[iter]] <- "identity"
                                } else {
                                  funs[[iter]] <- as.character(mf1[[1]])
                                }
                              }
                            }

                            private$funs <- rev(funs)
                            private$vars <- rev(vars)

                          },
                          genX = function(){
                            # generate model matrix X, including linearisation of non-linear terms,
                            X <- matrix(1,nrow=self$n(),ncol=1)
                            for(i in 1:length(private$funs)){
                              if(private$funs[[i]]=="RMINT")next
                              X <- cbind(X,
                                         do.call(paste0("d",private$funs[[i]]),list(list(
                                           data = as.matrix(self$data[,private$vars[[i]]]),
                                           pars = self$parameters[[i]]
                                         ))))
                            }
                            if(any(private$funs=="RMINT"))X <- X[,-1]
                            if(ncol(X)!=length(unlist(self$parameters)))stop("variables != parameters")
                            private$Xb <- X %*% matrix(unlist(self$parameters),ncol=1)
                            self$X <- Matrix::Matrix(X)

                          }
                        ))


###
# first order derivative functions

didentity <- function(x){
  return(x$data)
}

dfexp <- function(x){
  m <- as.matrix(x$data) %*% matrix(x$pars[2:length(x$pars)],ncol=1)
  X <- matrix(exp(m),ncol=1)
  for(i in 1:ncol(x$data)){
    X <- cbind(X,x$data[,i]*x$pars[i+1]*exp(m))
  }
  return(X)
}

dfactor <- function(x){
  matrix(model.matrix(~factor(a)-1,data.frame(a=x$data)),nrow=length(x$data))
}

dpexp <- function(x){

}

dlog <- function(x){

}
