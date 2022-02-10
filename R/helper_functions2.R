#FUNCTIONS TO GENERATE SPECIFIC DESIGN TYPES
# PLUS OTHER AUXILLIARY FUNCTIONS

# stepped wedge x
# cluster trial
# cluster cross over
# cohort study
# spatio-temporal sampling
# nelder split plot x

stepped_wedge <- function(J,
                          M,
                          nper=1,
                          beta=c(rep(0,J+1),0),
                          icc,
                          cac = NULL,
                          iac = NULL,
                          var = 1,
                          family = gaussian()){
  if(missing(icc))stop("icc must be set as a minimum")

  ndesigns <- length(icc) * ifelse(!is.null(cac[1]),length(cac),1) *
    ifelse(!is.null(iac[1]),length(iac),1)

  if(!is.null(cac[1]) && !is.na(cac[1])){
    wp_var <- icc[1]*var*(1-cac[1])
    bp_var <- icc[1]*var[1]*cac[1]
  } else {
    bp_var <- icc[1]*var
  }
  if(!is.null(iac[1]) && !is.na(iac[1])){
    ind_var <- var*(1-icc[1])*iac[1]
    sigma <- var*(1-icc[1])*(1-iac[1])
  } else {
    sigma <- var*(1-icc[1])
  }

  t <- J+1

  df <- nelder(formula(paste0("~ (J(",J,") > ind(",M,")) * t(",t,")")))
  df <- df[order(df$J,df$t),]
  int <- c()
  for(i in 1:J){
    int <- c(int,rep(c(rep(0,t-(i)),rep(1,i)),1))
  }
  df$int <- rep(int,each=M)

  if(is.null(cac[1]) || is.na(cac[1])){
   if(is.null(iac[1]) || is.na(iac[1])){
     f1 <- "~(1|gr(J))"
     pars <- list(list(bp_var))
   } else {
     f1 <- "~(1|gr(J)) + (1|gr(ind))"
     pars <- list(list(bp_var),list(ind_var))
   }
  } else {
    if(is.null(iac[1]) || is.na(iac[1])){
      f1 <- "~ (1|gr(J)) + (1|gr(J)*gr(t))"
      pars <- list(list(bp_var),list(1,wp_var))
    } else {
      f1 <- "~ (1|gr(J)) + (1|gr(J)*gr(t)) + (1|gr(ind))"
      pars <- list(list(bp_var),list(1,wp_var),list(ind_var))
    }
  }

  d1 <- Design$new(
    covariance = list(
      data=df,
      formula = f1,
      parameters = pars
    ),
    mean.function = list(
      formula = "~ factor(t) + int - 1",
      data = df,
      family = family,
      parameters = as.list(c(rep(0,t+1)))

    ),
    var_par = sigma
  )

  if(ndesigns>1){
    ds1 <- DesignSpace$new(d1)
    if(is.null(cac))cac <- NA
    if(is.null(iac))iac <- NA
    dsvalues <- expand.grid(icc=icc,cac=cac,iac=iac)

    for(i in 1:(ndesigns-1)){
      assign(paste0("d",i+1),
             Design$new(
               covariance = d1$covariance$clone(deep=TRUE),
               mean.function = d1$mean_function$clone(),
               var_par = 1
             ))

      if(!is.null(dsvalues$cac[i+1]) && !is.na(dsvalues$cac[i+1])){
        wp_var <- dsvalues$icc[i+1]*var*(1-dsvalues$cac[i+1])
        bp_var <- dsvalues$icc[i+1]*var*dsvalues$cac[i+1]
      } else {
        bp_var <- dsvalues$icc[i+1]*var
      }
      if(!is.null(dsvalues$iac[i+1]) && !is.na(dsvalues$iac[i+1])){
        ind_var <- var*(1-dsvalues$icc[i+1])*dsvalues$iac[i+1]
        sigma <- var*(1-dsvalues$icc[i+1])*(1-dsvalues$iac[i+1])
      } else {
        sigma <- var*(1-dsvalues$icc[i+1])
      }

      if(is.null(dsvalues$cac[i+1]) || is.na(dsvalues$cac[i+1])){
        if(is.null(dsvalues$iac[i+1]) || is.na(dsvalues$iac[i+1])){
          f1 <- "~(1|gr(J))"
          pars <- list(list(bp_var))
        } else {
          f1 <- "~(1|gr(J)) + (1|gr(ind))"
          pars <- list(list(bp_var),list(ind_var))
        }
      } else {
        if(is.null(dsvalues$iac[i+1]) || is.na(dsvalues$iac[i+1])){
          f1 <- "~ (1|gr(J)) + (1|gr(J)*gr(t))"
          pars <- list(list(bp_var),list(1,wp_var))
        } else {
          f1 <- "~ (1|gr(J)) + (1|gr(J)*gr(t)) + (1|gr(ind))"
          pars <- list(list(bp_var),list(1,wp_var),list(ind_var))
        }
      }


      ds1$add(get(paste0("d",i+1)))
      ds1$.__enclos_env__$private$designs[[i+1]]$var_par <- sigma
      ds1$.__enclos_env__$private$designs[[i+1]]$covariance <- Covariance$new(
        data=df,
        formula = f1,
        parameters = pars
      )

      # $parameters <- pars
      # ds1$.__enclos_env__$private$designs[[i+1]]$covariance$formula <- f1

    }
    return(ds1)
  } else {
    return(d1)
  }
}

# COHORT STUDY

# NELDER FUNCTION

#GENERATE CRISS CROSS
cross_df <- function(df1,df2){
  cnames <- c(colnames(df1),colnames(df2))
  df1 <- as.data.frame(df1[rep(1:nrow(df1),each=nrow(df2)),])
  df3 <- cbind(df1,df2[1:nrow(df2),])
  colnames(df3) <- cnames
  return(df3)
}

#GENERATE NESTING
nest_df <- function(df1,df2){
  df3 <- cbind(df1[rep(1:nrow(df1),each=nrow(df2)),],df2)
  colnames(df3)[1:ncol(df1)] <- colnames(df1)
  if(ncol(df1)>1)ids <- Reduce(paste0,df3[,1:ncol(df1)]) else ids <- df3[,1]
  df3[,(ncol(df1)+1):(ncol(df1)+ncol(df2))] <- apply(as.data.frame(df3[,(ncol(df1)+1):(ncol(df1)+ncol(df2))]),2,
                                                     function(i)as.numeric(as.factor(paste0(ids,i))))
  colnames(df3[,(ncol(df1)+1):(ncol(df1)+ncol(df2))]) <- colnames(df2)
  return(df3)
}

nelder <- function(formula){
  if(formula[[1]]=="~")formula <- formula[[2]]
  f1l <- formula[[2]]
  f1r <- formula[[3]]

  if(as.character(f1l[[1]])%in%c("*",">")){
    df1 <- Recall(f1l)
  } else if(as.character(f1l[[1]])%in%c("(")){
    df1 <- Recall(f1l[[2]])
  } else {
    df1 <- data.frame(a = seq(1,f1l[[2]]))
    colnames(df1) <- as.character(f1l[[1]])
  }

  if(as.character(f1r[[1]])%in%c("(")){
    df2 <- Recall(f1r[[2]])
  } else {
    df2 <- data.frame(a = seq(1,f1r[[2]]))
    colnames(df2) <- as.character(f1r[[1]])
  }

  if(formula[[1]] == "*"){
    df <- cross_df(df1,df2)
  } else if(formula[[1]] == ">"){
    df <- nest_df(df1,df2)
  }
  rownames(df) <- NULL
  return(df)
}

cycles <- function(a){
  fa <- a
  for(i in 1:(length(a)-1)){
    a <- c(a[2:length(a)],a[1])
    fa <- c(fa,a)
  }
  fa
}


# function to identify group membership
match_rows <- function(x,target,by){
  if(ncol(target)==1){
    tstr <- target[,by]
    xstr <- x[,by]
  } else {
    xstr <- Reduce(paste0,as.data.frame(apply(x[,by],2,function(i)paste0(i,".0000."))))
    tstr <- Reduce(paste0,as.data.frame(apply(target[,by],2,function(i)paste0(i,".0000."))))
  }
  Z <- matrix(0,nrow=length(xstr),ncol=length(tstr))
  mat <- lapply(tstr,function(i)which(xstr==i))
  for(i in 1:length(mat))Z[mat[[i]],i] <- 1
  return(Z)
}

## model non-linear functons below

fexp <- function(x){
  if(length(x$pars)!=2)stop("two parameters required for fexp")
  x$pars[1]*exp(-x$pars[2]*x$data)
}

pexp <- function(x){
  x$pars[1]^x$data
}

gr <- function(x){
  I(x$data==0)*x$pars[1]
}

## create block matrix

blockmat <- function(...){
  matlist <- list(...)
  n <- length(matlist)
  N <- 0:(n-1)
  rlist <- list()
  for(i in 1:n){
    N <- (N+1)%%n
    N[N==0] <- n
    rlist[[i]] <- Reduce(cbind,matlist[N])
  }
  Reduce(rbind,rlist)
}


# SPATIO-TEMPORAL SAMPLING
# build design for boundary and grid creation using rts2 code

# PARALLEL
