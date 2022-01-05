# require(Rcpp)
# require(RcppArmadillo)
#
# sourceCpp("rem_add_funcs.cpp")

# R versions of functions

choose_swap <- function(idx_in,A,sig,u){
  idx_out <- 1:nrow(sig)
  idx_out <- idx_out[!idx_out%in%idx_in]
  val_out <- sapply(1:length(idx_in),function(i)remove_one(A,i-1,u[idx_in]))
  rm1A <- remove_one_mat(A,which.max(val_out)-1)
  idx_in <- idx_in[-which.max(val_out)]
  val_in <- sapply(idx_out,function(i)add_one(rm1A,sig[i,i],sig[idx_in,i],u[c(idx_in,i)]))
  swap_idx <- idx_out[which.max(val_in)]
  newA <- add_one_mat(rm1A,sig[swap_idx,swap_idx],sig[idx_in,swap_idx])
  idx_in <- c(idx_in,swap_idx)
  val <- obj_fun(newA,u[idx_in])
  return(list(val,idx_in,newA))
}

choose_swap_robust <- function(idx_in,A_list,sig_list,u_list,weights){
  idx_out <- 1:nrow(sig_list[[1]])
  idx_out <- idx_out[!idx_out%in%idx_in]

  val_out_mat <- matrix(NA,nrow=length(idx_in),ncol=length(A_list))
  val_in_mat <- matrix(NA,nrow=length(idx_out),ncol=length(A_list))

  for(idx in 1:length(A_list)){
    val_out_mat[,idx] <- sapply(1:length(idx_in),function(i)
      remove_one(A_list[[idx]],i-1,u_list[[idx]][idx_in]))
  }
  val_out <- as.numeric(val_out_mat %*% weights)

  rm1A <- list()
  for(idx in 1:length(A_list)){
    rm1A[[idx]] <- remove_one_mat(A_list[[idx]],which.max(val_out)-1)
  }

  idx_in <- idx_in[-which.max(val_out)]

  for(idx in 1:length(A_list)){
    val_in_mat[,idx] <- sapply(idx_out,function(i)add_one(rm1A[[idx]],
                                                          sig_list[[idx]][i,i],
                                                          sig_list[[idx]][idx_in,i],
                                                          u_list[[idx]][c(idx_in,i)]))
  }
  val_in <- as.numeric(val_in_mat %*% weights)
  swap_idx <- idx_out[which.max(val_in)]

  newA <- list()
  for(idx in 1:length(A_list)){
    newA[[idx]] <- add_one_mat(rm1A[[idx]],sig_list[[idx]][swap_idx,swap_idx],
                               sig_list[[idx]][idx_in,swap_idx])
  }
  idx_in <- c(idx_in,swap_idx)
  val <- val_in[which.max(val_in)]
  return(list(val,idx_in,newA))
}

grad <- function(idx_in,A,sig,u,tol=1e-9, trace = TRUE){
  new_val <- obj_fun(A,u[idx_in])
  diff <- 1
  i <- 0
  while(diff > tol){
    val <- new_val
    i <- i + 1
    out <- choose_swap(idx_in,A,sig,u)
    new_val <- out[[1]]
    diff <- new_val - val
    if(diff>0){
      A <- out[[3]]
      idx_in <- out[[2]]
    }
    if (trace) {
      cat("\nIter: ",i)
      cat(" ",diff)
    }
  }
  return(idx_in)
}

#i've only changed this function, and added a couple of functions in gd_search.cpp
grad_robust <- function(idx_in,
                        C_list,
                        X_list,
                        sig_list,
                        w=NULL,
                        trace = TRUE,
                        rm_cols = NULL){

  if(is.null(w))w <- rep(1/length(sig_list),length(sig_list))
  if(sum(w)!=1)w <- w/sum(w)
  if(!is(w,"matrix"))w <- matrix(w,ncol=1)
  if(!all(unlist(lapply(C_list,function(x)is(x,"matrix")))))stop("All C_list must be matrices")
  if(!all(unlist(lapply(sig_list,function(x)is(x,"matrix")))))stop("All sig_list must be matrices")
  if(!all(unlist(lapply(X_list,function(x)is(x,"matrix")))))stop("All X_list must be matrices")
  if((length(C_list)!=length(X_list))|length(X_list)!=length(sig_list))stop("Lists must be same length")
  
  # added this function to give the user the option to remove columns from particular
  # designs quickly if the algorithm previously stopped and said to remove
  if(!is.null(rm_cols))
    {
    if(!is(rm_cols,"list"))stop("rm_cols should be a list")
    idx_original <- list()
    zero_idx <- c()
    idx_original <- 1:nrow(X_list[[1]])
    
    # find all the entries with non-zero values of the given columns in each design
    for(i in 1:length(rm_cols))
    {
      if(!is.null(rm_cols[[i]])){
        for(j in 1:length(rm_cols[[i]]))
        {
          zero_idx <- c(zero_idx,which(X_list[[i]][,rm_cols[[i]][j]]!=0))
        }
      }
    }
    zero_idx <- sort(unique(zero_idx))
    idx_original <- idx_original[-zero_idx]
    idx_in <- match(idx_in,idx_original)
    
    if(trace)message(paste0("removing ",length(zero_idx)," observations"))
    
    #update the matrices
    for(i in 1:length(rm_cols))
      {
      X_list[[i]] <- X_list[[i]][-zero_idx,-rm_cols[[i]]]
      C_list[[i]] <- matrix(C_list[[i]][-rm_cols[[i]]],ncol=1)
      sig_list[[i]] <- sig_list[[i]][-zero_idx,-zero_idx]
      
    }
    
    if(any(is.na(idx_in)))
    {
      if(trace)message("generating new random starting point")
      idx_in <- sample(1:nrow(X_list[[1]]),length(idx_in),replace=FALSE)
    }
  }
  
  #MAIN BODY OF THE FUNCTION
  
  # we need to calculate the M matrices for all the designs and store them
  # M is calculated for idx_in design rather than full design
  A_list <- list()
  u_list <- list()
  M_list <- list()
  for(i in 1:length(sig_list))
  {
    A_list[[i]] <- solve(sig_list[[i]][idx_in,idx_in])
    M_list[[i]] <- gen_m(X_list[[i]][idx_in,],A_list[[i]])
    cM <- t(C_list[[i]]) %*% solve(M_list[[i]])
    # print(cM%*%C_list[[i]])
    u_list[[i]] <- cM %*% t(X_list[[i]])
  }
  
  # the objective function here is now c^T M^-1 c - i've implemented c_obj_func in gd_search.cpp
  new_val_vec <- matrix(sapply(1:length(A_list),function(i)c_obj_fun(M_list[[i]], C_list[[i]])),nrow=1)
  new_val <- as.numeric(new_val_vec %*% w)
  
  diff <- -1
  i <- 0
  # we now need diff to be negative
  while(diff < 0)
    {
    #this code prints the M matrix if needed for debugging
    # for(j in 1:length(X_list)){
    #   M <- M_fun(C_list[[j]],X_list[[j]][sort(idx_in),],sig_list[[j]][sort(idx_in),sort(idx_in)])
    #   print(solve(M))
    #   print(cM%*%M%*%t(cM))
    # }
    
    val <- new_val
    i <- i + 1
    out <- choose_swap_robust(idx_in,A_list,sig_list,u_list, w)
    # new_val <- out[[1]]
    
    #we have to now recalculate all the lists of matrices for the new design proposed by the swap
    A_list <- out[[3]]
    for(j in 1:length(sig_list))
    {
      
      #checking matrix rank is relatively expensive. our problem is mostly due to removing all observations
      # from categorical variables, so within the loop we will just check if any of the columns have zero sum
      # and then stop the function. I've also added an option to the function to remove columns, see above
      if(any(colSums(X_list[[j]][out[[2]],])==0))
        {
        zero_col <- which(colSums(X_list[[j]][out[[2]],])==0)
        if(all(X_list[[j]][out[[2]],zero_col]==0))
          {
          stop(paste0("non-full rank information matrix. column ",zero_col," is not part of design ",j))
        }
      }
      
      
      M_list[[j]] <- gen_m(X_list[[j]][out[[2]],],A_list[[j]])
      cM <- t(C_list[[j]]) %*% solve(M_list[[j]])
      u_list[[j]] <- cM %*% t(X_list[[j]])
    }
    
   
    
    #calculate values for the new design - this is changed to new objective function
    new_val_vec <- matrix(sapply(1:length(A_list),function(j)c_obj_fun(M_list[[j]], C_list[[j]])),nrow=1)
    new_val <- as.numeric(new_val_vec %*% w)
    
    if (trace) {
      cat("\nIter: ",i)
      cat(" ",diff)
    }
    diff <- new_val - val
    # we are now looking for the smallest value rather than largest so diff<0
    if(diff<0){
      # A_list <- out[[3]]
      idx_in <- out[[2]]
    }

  }

  #check if matrix is full rank
  # algorithm should still work if not full rank due to factor covariates, 
  # but obviously conflicts with the original model so should warn user
  # and in some cases it won't
  # for(j in 1:length(X_list))
  # {
  #   r1 <- Matrix::rankMatrix(X_list[[j]][idx_in,])
  #   r2 <- ncol(X_list[[j]])
  #   if(r1[1]!=r2)message("solution does not have full rank, check model before using design.")
  # }
  
  ## if columns were removed then return the index to the indexing of the original X matrix
  if(!is.null(rm_cols))
  {
    idx_in <- idx_original[idx_in]
  }
  
  #return variance
  if(trace){
    var_vals <- c()
    for(i in 1:length(X_list))
    {
      var_vals[i] <- new_val_vec[i] #tryCatch(c_obj_fun(C_list[[i]],X_list[[i]][sort(idx_in),],sig_list[[i]][sort(idx_in),sort(idx_in)]),error=function(i){NA})
    }
    cat("\nVariance for individual model(s):\n")
    print(var_vals)
    if(length(A_list)>1){
      cat("\n Weighted average variance:\n")
      print(sum(var_vals*c(w)))
    }
  }
  
  
  
  return(idx_in)
}

optim_fun <- function(C,X,S){
  if(!any(is(C,"matrix"),is(X,"matrix"),is(S,"matrix")))stop("C, X, S must be matrices")
  M <- t(X) %*% solve(S) %*% X
  val <- t(C) %*% solve(M) %*% C
  return(val)
}

optim_fun2 <- function(C,X,S){
  if(!any(is(C,"matrix"),is(X,"matrix"),is(S,"matrix")))stop("C, X, S must be matrices")
  M <- t(X) %*% solve(S) %*% X
  val <- diag(solve(M))[c(C) != 0]
  return(val)
}

M_fun <- function(C,X,S){
  M <- t(X) %*% solve(S) %*% X
  return(M)
}


# sourceCpp("gd_search.cpp")
grad_robust2 <- function(idx_in, C_list, X_list, sig_list, w=NULL, tol=1e-9,
                        trace = TRUE){
  if(is.null(w))w <- rep(1/length(sig_list),length(sig_list))
  if(sum(w)!=1)w <- w/sum(w)
  if(!is(w,"matrix"))w <- matrix(w,ncol=1)
  if(!all(unlist(lapply(C_list,function(x)is(x,"matrix")))))stop("All C_list must be matrices")
  if(!all(unlist(lapply(sig_list,function(x)is(x,"matrix")))))stop("All sig_list must be matrices")
  if(!all(unlist(lapply(X_list,function(x)is(x,"matrix")))))stop("All X_list must be matrices")
  if((length(C_list)!=length(X_list))|length(X_list)!=length(sig_list))stop("Lists must be same length")

  A_list <- list()
  u_list <- list()
  for(i in 1:length(sig_list)){
    A_list[[i]] <- solve(sig_list[[i]][idx_in,idx_in])
    M <- t(X_list[[i]]) %*% solve(sig_list[[i]]) %*% X_list[[i]]
    cM <- t(C_list[[i]]) %*% solve(M)
    u_list[[i]] <- cM %*% t(X_list[[i]])
  }

  idx_in <- GradRobust(length(sig_list), idx_in-1, do.call(rbind, A_list),
            do.call(rbind, sig_list), do.call(cbind, u_list), w, tol, trace)
  idx_in <- idx_in + 1

  #return variance
  if(trace){
    var_vals <- c()
    for(i in 1:length(X_list)){
      var_vals[i] <- optim_fun(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in])
    }
    cat("\nVariance for individual model(s):\n")
    print(var_vals)
    if(length(A_list)>1){
      cat("\n Weighted average variance:\n")
      print(sum(var_vals*c(w)))
    }
  }
  return(idx_in)
}

# For a given m find the optimal power vector
max_var <- function(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE){

  # randomly generate starting position
  d <- sample(c(rep(1,m),rep(0,nrow(X_list[[1]])-m)),nrow(X_list[[1]]))
  idx_in <- which(d==1)

  if (length(idx_in) != nrow(X_list[[1]]))
  idx_in <- grad_robust2(idx_in, C_list, X_list, sig_list, w, 1e-9, trace)

  v0 <- c()
  for(i in 1:length(X_list)){
    v0 <- c(v0,optim_fun2(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in]))
  }

  v0
}


# For a given m find the optimal power vector
max_power <- function(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE){

  # randomly generate starting position
  d <- sample(c(rep(1,m),rep(0,nrow(X_list[[1]])-m)),nrow(X_list[[1]]))
  idx_in <- which(d==1)
  #idx_in <- (1:m)*round(nrow(X_list[[1]])/m,0)

  if (length(idx_in) != nrow(X_list[[1]]))
  idx_in <- grad_robust2(idx_in, C_list, X_list, sig_list, w, 1e-9, trace)

  v0 <- c()
  for(i in 1:length(X_list)){
    v0 <- c(v0,optim_fun2(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in]))
  }

  pow <- pnorm(sqrt(theta[unlist(C_list)!=0]/sqrt(v0)) - qnorm(1-alpha/2))

  pow
}

sample_size <- function(theta, alpha, pwr_target, m, C_list, X_list, sig_list, w) {
  iter <- 0
  pwr_new <- max_power(theta, alpha, m, C_list, X_list, sig_list, w)
  while (!all(pwr_new - pwr_target > 0) & m < nrow(X_list[[1]])) {
    iter <- iter + 1
    m <- m + 1
    pwr_new <- max_power(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE)
    cat("\nm = ", m)
    cat("\ntarget: ", pwr_target)
    cat("  minpwr: ", min(pwr_new))
  }
  return(m)
}

sample_size2 <- function(theta, alpha, pwr_target, C_list, X_list, sig_list, w) {

  cat("\nTarget power = ", pwr_target)

  lo <- max(unlist(lapply(C_list,function(i)length(unlist(i)))))*3
  hi <- nrow(X_list[[1]])
  pwr_new_lo <-NULL
  while(is.null(pwr_new_lo)){
    cat("\nlo = ", lo)
    pwr_new_lo <- tryCatch(
      max_power(theta, alpha, lo, C_list, X_list, sig_list, w, trace = FALSE),
      error=function(i)NULL)
    lo <- lo+10
  }
  pwr_new_hi <- max_power(theta, alpha, hi, C_list, X_list, sig_list, w, trace = FALSE)

  cat("\nmin power = ", min(pwr_new_lo))
  cat("\nmax power = ", min(pwr_new_hi))

  if (min(pwr_new_hi) < pwr_target | min(pwr_new_lo) > pwr_target)
    stop("\ntarget power is not in range of ", min(pwr_new_lo) , " and ", min(pwr_new_hi))

  v_hi    <- max_var(theta, alpha, hi, C_list, X_list, sig_list, w, trace = FALSE)
  v_target<- (theta[unlist(C_list)!=0]/(qnorm(pwr_target) + qnorm(1-alpha/2))^2)^2
  guess   <- round(max(v_hi / v_target * hi))
  pwr_new_guess <- max_power(theta, alpha, guess, C_list, X_list, sig_list, w, trace = FALSE)

  cat("\ninitial guess = ", guess, " with power = ", min(pwr_new_guess))

  if (min(pwr_new_guess) < pwr_target) lo <- guess
  if (min(pwr_new_guess) > pwr_target) hi <- guess

  while (lo <= hi) {
    mid <- lo + round((hi - lo) / 2)
    cat("\nlo = ", lo)
    cat("  hi = ", hi)
    cat(" mid = ", mid)
    pwr_new <- max_power(theta, alpha, mid, C_list, X_list, sig_list, w, trace = FALSE)
    if (pwr_target < min(pwr_new)) hi = mid - 1
    if (pwr_target > min(pwr_new)) lo = mid + 1
    cat("\ntarget: ", pwr_target)
    cat("  minpwr: ", min(pwr_new))
  }

  return(mid)
}
