
Covariance <- R6::R6Class("Covariance",
                      public = list(
                        data=NULL,
                        formula = NULL,
                        parameters = NULL,
                        Z = NULL,
                        D = NULL,
                        n= function(){
                          nrow(self$Z)
                        },
                        initialize = function(formula=NULL,
                                              data = NULL,
                                              parameters= NULL,
                                              verbose=TRUE){
                          if(any(is.null(data),is.null(formula),is.null(parameters))){
                            cat("not all attributes set. call check() when all attributes set.")
                          } else {
                            self$data = data
                            self$formula = as.formula(formula, env=.GlobalEnv)
                            self$parameters = parameters

                            private$cov_form()
                          }
                        },
                        check = function(verbose=TRUE){
                          new_hash <- private$hash_do()
                          if(private$hash[1] != new_hash[1]){
                            if(verbose)message("changes found, updating Z")
                            private$cov_form()
                          } else if(private$hash[2] != new_hash[2]){
                            if(verbose)message("changes found, updating D")
                            private$genD()
                          }

                          invisible(self)
                        },
                        print = function(){
                          # MAKE CLEARER ABOUT FUNCTIONS AND PARAMETERS

                          cat("Covariance\n")
                          print(self$formula)
                          cat("Parameters:")
                          print(unlist(self$parameters))
                          #print(head(self$data))
                        },
                        subset = function(index){
                          self$data <- self$data[index,]
                          self$check()
                        },
                        sampleD = function(parameters){
                          return(private$genD(update = FALSE,
                                              new_pars = parameters))
                        }
                      ),
                      private = list(
                        details = list(), # update to contain details of covariance for printing later
                        Distlist = NULL,
                        flist = NULL,
                        flistvars = NULL,
                        Zlist = NULL,
                        cov_functions = function(arg,x,pars){
                          if(arg=="exponential"){
                            return(pars[1]*exp(-x/pars[2]))
                          }
                          if(arg == "indicator"){
                            return(pars[1]*I(x==0))
                          }
                          if(arg == "exp_power"){
                            return(pars[1]^(x))
                          }
                        },
                        hash = NULL,
                        hash_do = function(){
                          c(digest::digest(c(self$data)),digest::digest(c(self$formula,self$parameters)))
                        },
                        cov_form = function(){
                          #1. split into independent components that can be combined in block diagonal form
                          flist <- list()
                          flistvars <- list()
                          formExt <- TRUE
                          count <- 0
                          form0 <- self$formula[[2]]
                          while(formExt){
                            count <- count + 1
                            formExt <- I(form0[[1]]=="+")
                            if(formExt){
                              flist[[count]] <- form0[[3]][[2]]
                              form0 <- form0[[2]]
                            } else {
                              flist[[count]] <- form0[[2]]
                            }
                            if("+"%in%all.names(flist[[count]][[3]]))stop("covariance only products")
                            rhsvars <- all.vars(flist[[count]][[3]])
                            funlist <- list()
                            rhsvargroup <- rep(NA,length(rhsvars))
                            formMult <- TRUE
                            countMult <- 0
                            form3 <- flist[[count]][[3]]
                            while(formMult){
                              countMult <- countMult + 1
                              formMult <- I(form3[[1]] == "*")
                              if(formMult){
                                funlist[[countMult]] <- form3[[3]][[1]]
                                rhsvargroup[match(all.vars(form3[[3]]),rhsvars)] <- countMult
                                form3 <- form3[[2]]
                              } else {
                                funlist[[countMult]] <- form3[[1]]
                                rhsvargroup[match(all.vars(form3),rhsvars)] <- countMult
                              }

                            }
                            flistvars[[count]] <- list(lhs=all.vars(flist[[count]][[2]]),
                                                       rhs = rhsvars,
                                                       funs = funlist,
                                                       groups = rhsvargroup)

                          }

                          # build each Z matrix and cbind
                          Zlist <- list()
                          Distlist <- list()
                          for(i in 1:length(flist)){
                            data_nodup <- self$data[!duplicated(self$data[,flistvars[[i]]$rhs]),flistvars[[i]]$rhs]
                            if(!is(data_nodup,"data.frame")){
                              data_nodup <- data.frame(data_nodup)
                              colnames(data_nodup) <- flistvars[[i]]$rhs
                            }
                            zdim2 <- nrow(data_nodup)
                            Zlist[[i]] <- match_rows(self$data,data_nodup,by=flistvars[[i]]$rhs)
                            if(length(flistvars[[i]]$lhs)>0){
                              ZlistNew <- list()
                              for(j in 1:length(flistvars[[i]]$lhs)){
                                ZlistNew[[j]] <- Zlist[[i]]*df[,flistvars[[i]]$lhs[j]]
                              }
                              Zlist[[i]] <- Reduce(cbind,ZlistNew)
                            }
                            Dist1 <- list()
                            for(j in 1:length(unique(flistvars[[i]]$groups))){
                              Dist1[[j]] <- as.matrix(dist(data_nodup[,flistvars[[i]]$rhs[flistvars[[i]]$groups==j]],upper = TRUE, diag=TRUE))
                            }
                            Distlist[[i]] <- Dist1
                          }
                          Z <- Reduce(cbind,Zlist)
                          Z <- Matrix::Matrix(Z)
                          if(ncol(Z)>nrow(Z))warning("Model underidentified")
                          self$Z <- Z
                          private$Distlist <- Distlist
                          private$flistvars <- flistvars
                          private$flist <- flist
                          for(i in 1:length(Zlist))Zlist[[i]] <- Matrix::Matrix(Zlist[[i]])
                          private$Zlist <- Zlist
                          private$genD()

                        },
                        genD = function(update=TRUE,
                                        new_pars=NULL){
                          D.sublist <- list()
                          for(d in 1:length(private$flist)){
                            ngroup <- length(unique(private$flistvars[[d]]$groups))
                            D.sublist[[d]] <- matrix(1, nrow=ncol(private$Zlist[[d]]),
                                                     ncol=ncol(private$Zlist[[d]]))
                            if(update){
                               pars <-rev(self$parameters)[[d]][1:ngroup]
                            } else {
                              pars <- rev(new_pars)[[d]]
                            }
                            for(j in 1:ngroup){
                              D.sublist[[d]] <- D.sublist[[d]]*
                                do.call(as.character(private$flistvars[[d]]$funs[[j]]),
                                        list(list(data=private$Distlist[[d]][[j]],
                                                  pars=pars[[j]])))
                            }
                            #remove AsIs class
                            class(D.sublist[[d]]) <- "matrix"

                            if(length(private$flistvars[[d]]$lhs)>0){
                              #fix cov parameter here
                              Dmatlist <- list(D.sublist[[d]])
                              for(j in 2:length(private$flistvars[[d]]$lhs)){
                                Dmatlist[[j]] <- diag(zdim2)*self$parameters[[d]][[ngroup-1+j]]
                              }
                              D.sublist[[d]] <- do.call("blockmat",Dmatlist)
                            }
                          }

                          #finally bdiag to combine them all
                          D <- do.call(Matrix::bdiag,D.sublist)

                          #add warning if number of re > n


                          if(update){
                            self$D <- D
                            private$hash <- private$hash_do()
                          } else {
                            return(D)
                          }

                        }
                      ))

