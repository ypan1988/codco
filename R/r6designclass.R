# DESIGN CLASS
Design <- R6::R6Class("Design",
                  public = list(
                    covariance = NULL,
                    mean_function = NULL,
                    exp_cond = NULL,
                    Sigma = NULL,
                    var_par = NULL,
                    initialize = function(covariance,
                                          mean.function,
                                          var_par = NULL,
                                          verbose=TRUE){
                      if(is(covariance,"R6")){
                        if(is(covariance,"Covariance")){
                          self$covariance <- covariance
                        } else {
                          stop("covariance should be Covariance class or list of appropriate arguments")
                        }
                      } else if(is(covariance,"list")){
                        self$covariance <- Covariance$new(
                          formula= covariance$formula,
                          data = covariance$data,
                          parameters = covariance$parameters,
                          verbose = verbose
                        )
                      }

                      if(is(mean.function,"R6")){
                        if(is(mean.function,"MeanFunction")){
                          self$mean_function <- mean.function
                        } else {
                          stop("mean.function should be MeanFunction class or list of appropriate arguments")
                        }
                      } else if(is(mean.function,"list")){
                        self$mean_function <- MeanFunction$new(
                          formula = mean.function$formula,
                          data = mean.function$data,
                          family = mean.function$family,
                          parameters = mean.function$parameters,
                          random_function = mean.function$random_function,
                          treat_par = mean.function$treat_par,
                          verbose = verbose
                        )
                      }



                      self$var_par <- var_par

                      self$generate()
                      private$hash <- private$hash_do()
                    },
                    print = function(){
                      cat("\n----------------------------------------\n")
                      print(self$mean_function)
                      cat("\n----------------------------------------\n")
                      print(self$covariance)
                      cat("\n----------------------------------------\n")
                    },
                    n = function(){
                      self$mean_function$n()
                    },
                    generate = function(){
                      # add check for var par with gaussian family

                      private$genW(family = self$mean_function$family,
                                   Xb = self$mean_function$.__enclos_env__$private$Xb,
                                   var_par = self$var_par)
                      private$genS(D = self$covariance$D,
                                   Z = self$covariance$Z,
                                   W = private$W)
                    },
                    power = function(par,
                                     value,
                                     alpha=0.05,
                                     type="linear",
                                     method=NULL,
                                     iter=10,
                                     skip.check = FALSE,
                                     parallel=TRUE){
                      if(!skip.check)self$check(verbose=FALSE)
                      if(missing(par)|missing(value))stop("parameter missing")
                      if(type=="linear"){
                        old_par <- self$mean_function$parameters[[par]]
                        self$mean_function$parameters[[par]] <- value
                        self$check(verbose=FALSE)
                        M <- private$information_matrix()
                        v0 <- solve(M)[par,par]
                        pwr <- pnorm(value/(sqrt(v0)) - qnorm(1-alpha/2))
                        self$mean_function$parameters[[par]] <- old_par
                      } else if(type=="sim"){
                        if(parallel){
                          cl <- parallel::makeCluster(parallel::detectCores()-2)
                          parallel::clusterEvalQ(cl,require(Matrix))
                          parallel::clusterEvalQ(cl,require(lme4))
                          parallel::clusterEvalQ(cl,require(minqa))
                         ests <- pbapply::pbreplicate(iter, private$lme_est(par=par,
                                                                            value=value),
                                                      cl=cl)
                         parallel::stopCluster(cl)
                        } else {
                          ests <- pbapply::pbreplicate(iter, private$lme_est(par=par,
                                                                             value=value))
                        }
                         ests <- apply(ests,2,function(x)x$b/x$se)
                         ests <- ests[par,]
                         ests <- (1-pnorm(abs(ests)))*2
                         pwr <- mean(I(ests < alpha),na.rm=TRUE)
                         if(any(is.na(ests)))message(paste0("failure to converge in ",length(which(is.na(ests)))," models"))
                      }
                      return(pwr)
                    },
                    subset_rows = function(index){
                      self$mean_function$subset_rows(index)
                      self$covariance$subset(index)
                    },
                    subset_cols = function(index){
                      self$mean_function$subset_cols(index)
                    },
                    plot = function(x,
                                    y,
                                    z=NULL,
                                    treat){
                      if(is.null(z)){
                        ggplot(data=self$covariance$data,aes(x=.data[[x]],y=.data[[y]]))+
                          geom_count(aes(color=self$mean_function$data[,treat]))+
                          theme_bw()+
                          theme(panel.grid=element_blank())+
                          scale_color_viridis_c(name=treat)+
                          scale_size_area()
                      } else {
                        ggplot(data=self$covariance$data,aes(x=.data[[x]],y=.data[[y]]))+
                          geom_count(aes(color=self$mean_function$data[,treat]))+
                          facet_wrap(~.data[[z]])+
                          theme_bw()+
                          theme(panel.grid=element_blank())+
                          scale_color_viridis_c(name=treat)+
                          scale_size_area()
                      }},
                    sim_data = function(type = "y",
                                        par=NULL,
                                        value=NULL){
                      if(!is.null(par)){
                        if(is.null(value))stop("set parameter value")
                        orig_par <- self$mean_function$parameters[par]
                        self$mean_function$parameters[par] <- value
                      }
                      re <- MASS::mvrnorm(n=1,mu=rep(0,nrow(self$covariance$D)),Sigma = self$covariance$D)
                      mu <- c(drop(self$mean_function$.__enclos_env__$private$Xb)) + as.matrix(self$covariance$Z%*%re)
                   
                      f <- self$mean_function$family
                      if(f[1]=="poisson"){
                        if(f[2]=="log"){
                          y <- rpois(self$n(),exp(mu))
                        }
                        if(f[2]=="identity"){
                          y <- rpois(self$n(),mu)
                        }
                      }

                      if(f[1]=="binomial"){
                        if(f[2]=="logit"){
                          y <- rbinom(self$n(),1,exp(mu)/(1+exp(mu)))
                        }
                        if(f[2]=="log"){
                          y <- rbinom(self$n(),1,exp(mu))
                        }
                        if(f[2]=="identity"){
                          y <- rbinom(self$n(),1,mu)
                        }
                        if(f[2]=="probit"){
                          y <- rbinom(self$n(),1,pnorm(mu))
                        }
                      }

                      if(f[1]=="gaussian"){
                        if(f[2]=="identity"){
                          if(is.null(self$var_par))stop("For gaussian(link='identity') provide var_par")
                          y <- rnorm(self$n(),mu,self$var_par)
                        }
                        if(f[2]=="log"){
                          if(is.null(self$var_par))stop("For gaussian(link='log') provide var_par")
                          #CHECK THIS IS RIGHT
                          y <- rnorm(self$n(),exp(mu),self$var_par)
                        }
                      }

                      if(f[1]=="gamma"){
                        if(f[2]=="inverse"){
                          if(is.null(self$var_par))stop("For gamma(link='inverse') provide var_par")
                          #CHECK THIS IS RIGHT
                          y <- rgamma(self$n(),shape = 1/(mu*self$var_par),rate = 1/self$var_par)
                        }
                      }
                      if(!missing(par))self$mean_function$parameters[par]<-orig_par
                      if(type=="y")return(y)
                      if(type=="data.frame")return(cbind(y,self$data$data,self$covariance$location))
                    },
                    check = function(verbose=TRUE){
                      self$covariance$check(verbose=verbose)
                      self$mean_function$check(verbose = verbose)
                      if(private$hash != private$hash_do()){
                        self$generate()
                      }
                    },
                    apv = function(prior,
                                   var,
                                   prior.fun,
                                   iter,
                                   verbose=TRUE){
                      if(verbose)message("Monte Carlo integration")
                      samps <- pbapply::pbreplicate(iter,self$posterior(prior,var,do.call(prior.fun,list())))
                      summary(samps)
                    },
                    posterior = function(prior,
                                         var,
                                         parameters){
                      #move to private and set this as Monte Carlo integration
                      #can just request a function that outputs a new set of covariance parameters
                      R <- solve(Matrix::Matrix(diag(prior)))
                      S <- private$genS(self$covariance$sampleD(parameters),self$covariance$Z,private$W,update=FALSE)
                      M <- R + Matrix::crossprod(self$mean_function$X,solve(S))%*%self$mean_function$X
                      M <- solve(M)
                      M[var,var]
                    }
                  ),
                  private = list(
                    W = NULL,
                    Xb = NULL,
                    logit = function(x){
                      exp(x)/(1+exp(x))
                    },
                    genW = function(family,
                                    Xb,
                                    var_par=NULL){
                      # assume random effects value is at zero
                      f <- family
                      Xb <- c(Xb)
                      if(!f[1]%in%c("poisson","binomial","gaussian","gamma"))stop("family must be one of Poisson, Binomial, Gaussian, Gamma")

                      if(f[1]=="poisson"){
                        if(f[2]=="log"){
                          W <- diag(1/(exp(Xb)))
                        }
                        if(f[2]=="identity"){
                          W <- diag(exp(Xb))
                        }
                      }

                      if(f[1]=="binomial"){
                        if(f[2]=="logit"){
                          W <- diag(1/(private$logit(Xb)*(1-private$logit(Xb))))
                        }
                        if(f[2]=="log"){
                          W <- diag((1-private$logit(Xb))/(private$logit(Xb)))
                        }
                        if(f[2]=="identity"){
                          W <- diag((private$logit(Xb)*(1-private$logit(Xb))))
                        }
                        if(f[2]=="probit"){
                          W <- diag((pnorm(Xb)*(1-pnorm(Xb)))/(dnorm(Xb)))
                        }
                      }

                      if(f[1]=="gaussian"){
                        if(f[2]=="identity"){
                          if(is.null(var_par))stop("For gaussian(link='identity') provide var_par")
                          W <- var_par*diag(length(Xb))
                        }
                        if(f[2]=="log"){
                          if(is.null(var_par))stop("For gaussian(link='log') provide var_par")
                          W <- diag(var_par/exp(Xb))
                        }
                      }

                      if(f[1]=="gamma"){
                        if(f[2]=="inverse"){
                          if(is.null(var_par))stop("For gamma(link='inverse') provide var_par")
                          W <- var_par*diag(length(Xb))
                        }
                      }
                      private$W <- Matrix::Matrix(W)
                    },
                    genS = function(D,Z,W,update=TRUE){
                      if(is(D,"numeric")){
                        S <- W + D * Matrix::tcrossprod(Z)
                      } else {
                        S <- W + Z %*% Matrix::tcrossprod(D,Z)
                      }
                      if(update){
                        self$Sigma <- Matrix::Matrix(S)
                        private$hash <- private$hash_do()
                      } else {
                        return(S)
                      }

                    },
                    hash = NULL,
                    hash_do = function(){
                      digest::digest(c(self$covariance$.__enclos_env__$private$hash,
                                       self$mean_function$.__enclos_env__$private$hash))
                    },
                    lme_est = function(par=NULL,
                                       value=NULL){
                      lambda <- Matrix::t(Matrix::chol(self$covariance$D))#,nrow=nrow(self$covariance$D)))
                      #lambda <- Matrix::Matrix(lambda)
                      
                      parInds <- list(covar = c(1,2),
                                      fixef = c(3:(2+ncol(self$mean_function$X))))
                      devfunList <- list(Lind = seq_along(lambda@x),
                                         pp = lme4::merPredD$new(
                                           X = self$mean_function$X,
                                           Zt = Matrix::t(self$covariance$Z),
                                           Lambdat = lambda,
                                           Lind = seq_along(lambda@x),
                                           theta = as.double(lambda@x),
                                           n = self$n()),
                                         resp = lme4::glmResp$new(
                                           y = self$sim_data(par=par,
                                                             value=value),
                                           family = self$mean_function$family,
                                           weights = rep(1, self$n())),
                                         lp0 = NULL,
                                         baseOffset = rep(0, self$n()),
                                         tolPwrss = 1e-6,
                                         maxit = 30,
                                         GQmat = lme4::GHrule(1),
                                         compDev = TRUE,
                                         fac = NULL,
                                         verbose = TRUE,
                                         parInds = parInds)
                      #NEED TO GENERALISE TO ALL COVARIANCES
                      cov2 <- self$covariance$clone(deep=TRUE)
                      

                      updateTheta <- function(pars){
                        cov2$parameters <- relist(cov2$parameters,
                                                  value = pars)[[1]]
                        #cov2$parameters <- list(list(pars[1],pars[2]))
                        cov2$check(verbose=FALSE)
                        newD <- cov2$D
                        cholD <- tryCatch(Matrix::chol(newD),error=function(e)NULL)
                        if(is.null(cholD)){
                          newD <- Matrix::nearPD(matrix(newD,nrow=nrow(cov2$D)))
                          cholD <- tryCatch(Matrix::chol(newD),error=function(e)print("help!"))
                        }
                        Matrix::t(cholD)@x
                      }
                      
                      # CHANGE FUNCTION OPTS BELOW
                      devfun <- function(pars) {
                        pp$setTheta(as.double(updateTheta(pars[parInds$covar])))
                        spars <- as.numeric(pars[parInds$fixef])
                        offset <- if (length(spars)==0) baseOffset else baseOffset + pp$X %*% spars
                        resp$setOffset(offset)
                        p <- lme4::glmerLaplaceHandle(pp$ptr(), resp$ptr(), 0, 1e-6, 30, TRUE)
                        resp$updateWts()
                        p
                      }

                      devfunEnv <- new.env()
                      environment(devfun) <- list2env(devfunList, envir = devfunEnv)
                      environment(devfun)$lp0 <- environment(devfun)$pp$linPred(1)
                      
                      n.cov.par <- length(unlist(cov2$parameters))

                      opt <- tryCatch(minqa::bobyqa(par = c(rep(0.1,n.cov.par),rep(0,ncol(self$mean_function$X))),
                                           fn = devfun,
                                           lower = c(rep(1e-5,n.cov.par),rep(-Inf,ncol(self$mean_function$X)))),error=function(e)NULL)
                      if(is.null(opt)){
                        return(data.frame(b=rep(NA,ncol(self$mean_function$X)),se=rep(NA,ncol(self$mean_function$X))))
                      } else {
                        cov2$parameters <- relist(cov2$parameters,
                                                  value = opt$par[seq_len(n.cov.par)])[[1]]#list(list(opt$par[1],opt$par[2]))
                        cov2$check(verbose=FALSE)
                        se <- sqrt(diag(solve(Matrix::crossprod(self$mean_function$X,cov2$Z)%*%cov2$D%*%Matrix::crossprod(cov2$Z,self$mean_function$X))))
                        b <- opt$par[parInds$fixef]

                        return(data.frame(b,se))
                      }

                    },
                    information_matrix = function(){
                      Matrix::crossprod(self$mean_function$X,solve(self$Sigma))%*%self$mean_function$X
                    }
                  ))


relist <- function(lst,value,p=0){
  if(is(lst,"list")){
    for(i in 1:length(lst)){
      out <- Recall(lst[[i]],value,p=p)
      lst[[i]] <- out[[1]]
      p <- out[[2]]
    }
  } else {
    for(i in 1:length(lst)){
      lst[i] <- value[p+1]
      p <- p + 1
    }
  }
  return(list(lst,p))
}



