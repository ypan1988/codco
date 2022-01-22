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
  #print(varval - val_in)
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

# check if psd using chol
check_psd <- function(M){
  cholM <- tryCatch(chol(M),error=function(e)NA)
  if(!is(cholM,"matrix"))
  {
    return(FALSE)
  } else {
    return(TRUE)
  }
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
      M_list[[j]] <- gen_m(X_list[[j]][out[[2]],],A_list[[j]])

      # check positive semi-definite before proceeding
      if(!check_psd(M_list[[j]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[j]])==0),
                                             " of design ",j," is not part of an optimal design."))

      cM <- t(C_list[[j]]) %*% solve(M_list[[j]])
      u_list[[j]] <- cM %*% t(X_list[[j]])
    }

    #calculate values for the new design - this is changed to new objective function
    new_val_vec <- matrix(sapply(1:length(A_list),function(j)c_obj_fun(M_list[[j]], C_list[[j]])),nrow=1)
    new_val <- as.numeric(new_val_vec %*% w)

    diff <- new_val - val
    if (trace) {
      cat("\nIter: ",i)
      cat(" ",diff)
    }

    # we are now looking for the smallest value rather than largest so diff<0
    if(diff<0){
      # A_list <- out[[3]]
      idx_in <- out[[2]]
    }

  }

  # as a check, see if the next swap would mean that M is not positive semi-definite, in which
  # case the algorithm has terminated early and the solution is not in this design
  # this code is copied from the choose_swap_robust function, so could be streamlined into its own function perhaps

  idx_test <- out[[2]]
  val_out_mat <- matrix(NA,nrow=length(idx_test),ncol=length(A_list))
  for(idx in 1:length(A_list)){
    val_out_mat[,idx] <- sapply(1:length(idx_test),function(i)
      remove_one(A_list[[idx]],i-1,u_list[[idx]][idx_test]))
  }
  val_out <- as.numeric(val_out_mat %*% w)
  idx_test <- idx_test[-which.max(val_out)]
  rm1A <- list()
  for(idx in 1:length(A_list)){
    rm1A[[idx]] <- remove_one_mat(A_list[[idx]],which.max(val_out)-1)
  }
  for(j in 1:length(sig_list))
  {
    M_list[[j]] <- gen_m(X_list[[j]][idx_test,],rm1A[[j]])

    if(!check_psd(M_list[[j]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[j]])==0),
                                           " of design ",j," is not part of an optimal design."))
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

#Updated version of algorithm below, plus reorder_obs() function below it
grad_robust_step <- function(idx_in,
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

  #define sample size
  n <- length(idx_in)
  N <- nrow(X_list[[1]])

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
      set.seed(888)
      idx_in <- sample(1:nrow(X_list[[1]]),n,replace=FALSE)
    }
  }

  #MAIN BODY OF THE FUNCTION

  # we need to calculate the M matrices for all the designs and store them
  # M is calculated for idx_in design rather than full design
  A_list <- list()
  A_list_sub <- list()
  u_list <- list()
  M_list <- list()
  rm1A <- list()
  for(i in 1:length(sig_list))
  {
    A_list[[i]] <- solve(sig_list[[i]][idx_in,idx_in])
    M_list[[i]] <- gen_m(X_list[[i]][idx_in,],A_list[[i]])
    u_list[[i]] <- (t(C_list[[i]]) %*% solve(M_list[[i]])) %*% t(X_list[[i]])
  }

  # the objective function here is now c^T M^-1 c - i've implemented c_obj_func in gd_search.cpp
  new_val_vec <- matrix(sapply(1:length(A_list),function(i)c_obj_fun(M_list[[i]], C_list[[i]])),nrow=1)
  new_val <- as.numeric(new_val_vec %*% w)

  diff <- -1
  i <- 0

  #1. order the list of observations to include and exclude
  # idx_in are the observations in the design and idx_out are the observations not in the design
  idx_list <- reorder_obs(idx_in,
                          A_list,
                          sig_list,
                          u_list,
                          w)
  idx_in <- idx_list[[1]]
  idx_out <- idx_list[[2]]
  A_list <- idx_list[[3]]

  while(diff < 0)
  {
    val <- new_val
    i <- i + 1

    best_val_vec <- new_val_vec
    #2. make swap at top of list

    #grad <- c()
    for(idx in 1:length(A_list)){
      rm1A[[idx]] <- remove_one_mat(A_list[[idx]],0)
      A_list_sub[[idx]] <- add_one_mat(rm1A[[idx]],
                                       sig_list[[idx]][idx_out[1],idx_out[1]],
                                       sig_list[[idx]][idx_in[2:n],idx_out[1]])
      M_list[[idx]] <- gen_m(X_list[[idx]][c(idx_in[2:n],idx_out[1]),],A_list_sub[[idx]])

      # check positive semi-definite before proceeding
      if(!check_psd(M_list[[idx]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[idx]])==0),
                                               " of design ",idx," is not part of an optimal design."))

    }

    #calculate values for the new design - this is changed to new objective function
    new_val_vec <- matrix(sapply(1:length(A_list),function(j)c_obj_fun(M_list[[j]], C_list[[j]])),nrow=1)
    new_val <- as.numeric(new_val_vec %*% w)

    diff <- new_val - val

    if (trace) {
      cat("\nIter: ",i)
      cat(" ",val)
    }

    if(diff<0){
      swap_idx <- idx_in[1]
      idx_in <- c(idx_in[2:n],idx_out[1])
      idx_out <- c(idx_out[2:n],swap_idx)
      A_list <- A_list_sub
      for(idx in 1:length(A_list))
      {
        u_list[[idx]] <- (t(C_list[[idx]]) %*% solve(M_list[[idx]])) %*% t(X_list[[idx]])
      }
    } else {
      #if no easy swaps can be made, reorder the list and try the top one again
      i <- i + 1
      idx_list <- reorder_obs(idx_in,
                              A_list,
                              sig_list,
                              u_list,
                              w)
      idx_in <- idx_list[[1]]
      idx_out <- idx_list[[2]]
      A_list <- idx_list[[3]]

      for(idx in 1:length(A_list))
      {
        rm1A[[idx]] <- remove_one_mat(A_list[[idx]],0)
        A_list_sub[[idx]] <- add_one_mat(rm1A[[idx]],
                                         sig_list[[idx]][idx_out[1],idx_out[1]],
                                         sig_list[[idx]][idx_in[2:n],idx_out[1]])

        M_list[[idx]] <- gen_m(X_list[[idx]][c(idx_in[2:n],idx_out[1]),],A_list_sub[[idx]])


        # check positive semi-definite before proceeding
        if(!check_psd(M_list[[idx]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[idx]])==0),
                                                 " of design ",idx," is not part of an optimal design."))
      }

      #calculate values for the new design - this is changed to new objective function
      new_val_vec <- matrix(sapply(1:length(A_list),function(j)c_obj_fun(M_list[[j]], C_list[[j]])),nrow=1)
      new_val <- as.numeric(new_val_vec %*% w)

      diff <- new_val - val

      # we are now looking for the smallest value rather than largest so diff<0
      if (trace) {
        cat("\nIter (reorder): ",i)
        cat(" ",val)
      }

      if(diff<0){
        swap_idx <- idx_in[1]
        idx_in <- c(idx_in[2:n],idx_out[1])
        idx_out <- c(idx_out[2:n],swap_idx)
        A_list <- A_list_sub
        for(idx in 1:length(A_list))
        {
          u_list[[idx]] <- (t(C_list[[idx]]) %*% solve(M_list[[idx]])) %*% t(X_list[[idx]])
        }
      } else {
        # if reordering doesn't find a better solution then check all the neighbours
        #check optimality and see if there are any neighbours that improve the solution
        # I have got it to go through in order of the ordered observations, and then break the loop if
        # it reaches one with a positive gradient because we know the ones after that will be worse

        if(trace)cat("\nChecking optimality...\n")
        psd_check <- c()
        val_in_mat <- c()
        for(obs in 1:n)
        {
          if(trace)cat(paste0("\rChecking neighbour block: ",obs," of ",n))
          rm1A <- list()
          for(idx in 1:length(A_list))
          {
            rm1A[[idx]] <- remove_one_mat(A_list[[idx]],obs-1)
          }
          for(obs_j in 1:(N-n))
          {
            for(idx in 1:length(A_list))
            {
              val_in_mat[idx] <- add_one(rm1A[[idx]],
                                     sig_list[[idx]][idx_out[obs_j],idx_out[obs_j]],
                                     sig_list[[idx]][idx_in[-obs],idx_out[obs_j]],
                                     u_list[[idx]][c(idx_in[-obs],idx_out[obs_j])])
            }
            val_in <- sum(val_in_mat*w)
            if(val - val_in < 0){
              for(idx in 1:length(A_list))
              {
                A_list_sub[[idx]] <- add_one_mat(rm1A[[idx]],
                                                 sig_list[[idx]][idx_out[obs_j],idx_out[obs_j]],
                                                 sig_list[[idx]][idx_in[-obs],idx_out[obs_j]])

                M_list[[idx]] <- gen_m(X_list[[idx]][c(idx_in[-obs],idx_out[obs_j]),],A_list_sub[[idx]])
                psd_check[idx] <- check_psd(M_list[[idx]])
              }

              if(any(!psd_check))next
              new_val_vec <- matrix(sapply(1:length(A_list),function(k)c_obj_fun(M_list[[k]], C_list[[k]])),nrow=1)
              new_val <- as.numeric(new_val_vec %*% w)
              if(new_val - val < 0 )
              {
                diff <- new_val - val
                if(trace)cat("\nImprovement found: ",new_val)
                swap_idx <- idx_in[obs]
                idx_in <- c(idx_in[-obs],idx_out[obs_j])
                idx_out <- c(idx_out[-obs_j],swap_idx)
                A_list <- A_list_sub
                for(idx in 1:length(A_list))
                {
                  u_list[[idx]] <- (t(C_list[[idx]]) %*% solve(M_list[[idx]])) %*% t(X_list[[idx]])
                }
                break
              }
            } else {
              break
            }
          }
          if(diff < 0)break
        }




      }

    }

  }

  # as a check, see if the next swap would mean that M is not positive semi-definite, in which
  # case the algorithm has terminated early and the solution is not in this design
  # this code is copied from the choose_swap_robust function, so could be streamlined into its own function perhaps

  idx_test <- c(idx_in[2:n],idx_out[1])
  val_out_mat <- matrix(NA,nrow=length(idx_test),ncol=length(A_list))
  for(idx in 1:length(A_list)){
    val_out_mat[,idx] <- sapply(1:length(idx_test),function(i)
      remove_one(A_list[[idx]],i-1,u_list[[idx]][idx_test]))
  }
  val_out <- as.numeric(val_out_mat %*% w)
  idx_test <- idx_test[-which.max(val_out)]
  rm1A <- list()
  for(idx in 1:length(A_list)){
    rm1A[[idx]] <- remove_one_mat(A_list[[idx]],which.max(val_out)-1)
  }
  for(j in 1:length(sig_list))
  {
    M_list[[j]] <- gen_m(X_list[[j]][idx_test,],rm1A[[j]])

    if(!check_psd(M_list[[j]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[j]])==0),
                                           " of design ",j," is not part of an optimal design."))
  }


  ## if columns were removed then return the index to the indexing of the original X matrix
  if(!is.null(rm_cols))
  {
    idx_in <- idx_original[idx_in]
  }

  #return variance
  if(trace){
    cat("\nVariance for individual model(s):\n")
    print(c(best_val_vec))
    if(length(A_list)>1){
      cat("\n Weighted average variance:\n")
      print(sum(best_val_vec*c(w)))
    }
  }

  return(idx_in)
}

#function to reorder the observations
reorder_obs <- function(idx_in,
                        A_list,
                        sig_list,
                        u_list,
                        weights){

  idx_out <- 1:nrow(sig_list[[1]])
  idx_out <- idx_out[!idx_out%in%idx_in]
  n <- length(idx_in)

  val_out_mat <- matrix(NA,nrow=length(idx_in),ncol=length(A_list))
  val_in_mat <- matrix(NA,nrow=length(idx_out),ncol=length(A_list))

  for(idx in 1:length(A_list)){
    val_out_mat[,idx] <- sapply(1:length(idx_in),function(i)
      remove_one(A_list[[idx]],i-1,u_list[[idx]][idx_in]))
  }
  val_out <- as.numeric(val_out_mat %*% weights)
  #print(head(data.frame(val=sort(val_out,decreasing = TRUE),idx=idx_in[order(val_out,decreasing = TRUE)])))
  idx_in <- idx_in[order(val_out, decreasing = TRUE)]

  #rm1A <- list()

  for(idx in 1:length(A_list)){
    A_list[[idx]] <- A_list[[idx]][order(val_out, decreasing = TRUE),order(val_out, decreasing = TRUE)]
    #rm1A[[idx]] <- remove_one_mat(A_list[[idx]],0)
    val_in_mat[,idx] <- sapply(idx_out,function(i)add_one(A_list[[idx]],
                                                          sig_list[[idx]][i,i],
                                                          sig_list[[idx]][idx_in,i],
                                                          u_list[[idx]][c(idx_in,i)]))
  }

  val_in <- as.numeric(val_in_mat %*% weights)
  #print(head(data.frame(val=sort(val_in,decreasing = TRUE),idx=idx_out[order(val_in,decreasing = TRUE)])))
  idx_out <- idx_out[order(val_in,decreasing = TRUE)]

  return(list(idx_in,idx_out,A_list))
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
grad_robust2 <-function(idx_in, C_list, X_list, sig_list, w=NULL, tol=1e-9,
                        trace = TRUE, rm_cols = NULL){
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
  for(i in 1:length(sig_list)){
    A_list[[i]] <- solve(sig_list[[i]][idx_in,idx_in])
    #M <- t(X_list[[i]]) %*% solve(sig_list[[i]]) %*% X_list[[i]]
    #cM <- t(C_list[[i]]) %*% solve(M)
    M_list[[i]] <- gen_m(X_list[[i]][idx_in,],A_list[[i]])
    cM <- t(C_list[[i]]) %*% solve(M_list[[i]])
    # print(cM%*%C_list[[i]])
    u_list[[i]] <- cM %*% t(X_list[[i]])
  }

  out_list <- GradRobustAlg1(idx_in -1, do.call(rbind, C_list), do.call(rbind, X_list), do.call(rbind, sig_list), weights = w)
  idx_in <- out_list[["idx_in"]] + 1
  idx_out <- out_list[["idx_out"]] + 1
  best_val_vec <- out_list[["best_val_vec"]]

  # idx_in <- GradRobust(length(sig_list), idx_in-1, do.call(rbind, A_list),
  #                      do.call(rbind, M_list),do.call(rbind, C_list),do.call(rbind, X_list),
  #                      do.call(rbind, sig_list), do.call(cbind, u_list), w, tol, trace)
  # idx_in <- idx_in + 1

  # as a check, see if the next swap would mean that M is not positive semi-definite, in which
  # case the algorithm has terminated early and the solution is not in this design
  # this code is copied from the choose_swap_robust function, so could be streamlined into its own function perhaps

  # idx_test <- out[[2]]
  # val_out_mat <- matrix(NA,nrow=length(idx_test),ncol=length(A_list))
  # for(idx in 1:length(A_list)){
  #   val_out_mat[,idx] <- sapply(1:length(idx_test),function(i)
  #     remove_one(A_list[[idx]],i-1,u_list[[idx]][idx_test]))
  # }
  # val_out <- as.numeric(val_out_mat %*% w)
  # idx_test <- idx_test[-which.max(val_out)]
  # rm1A <- list()
  # for(idx in 1:length(A_list)){
  #   rm1A[[idx]] <- remove_one_mat(A_list[[idx]],which.max(val_out)-1)
  # }
  # for(j in 1:length(sig_list))
  # {
  #   M_list[[j]] <- gen_m(X_list[[j]][idx_test,],rm1A[[j]])
  #
  #   if(!check_psd(M_list[[j]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[j]])==0),
  #                                          " of design ",j," is not part of an optimal design."))
  # }

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

grad_robust2_step <-function(idx_in, C_list, X_list, sig_list, w=NULL, tol=1e-9,
                             trace = TRUE, rm_cols = NULL){
  if(is.null(w))w <- rep(1/length(sig_list),length(sig_list))
  if(sum(w)!=1)w <- w/sum(w)
  if(!is(w,"matrix"))w <- matrix(w,ncol=1)
  if(!all(unlist(lapply(C_list,function(x)is(x,"matrix")))))stop("All C_list must be matrices")
  if(!all(unlist(lapply(sig_list,function(x)is(x,"matrix")))))stop("All sig_list must be matrices")
  if(!all(unlist(lapply(X_list,function(x)is(x,"matrix")))))stop("All X_list must be matrices")
  if((length(C_list)!=length(X_list))|length(X_list)!=length(sig_list))stop("Lists must be same length")

  #define sample size
  n <- length(idx_in)
  N <- nrow(X_list[[1]])

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
      set.seed(888)
      idx_in <- sample(1:nrow(X_list[[1]]),n,replace=FALSE)
    }
  }

  #MAIN BODY OF THE FUNCTION

  # we need to calculate the M matrices for all the designs and store them
  # M is calculated for idx_in design rather than full design
  A_list <- list()
  A_list_sub <- list()
  u_list <- list()
  M_list <- list()
  rm1A <- list()
  for(i in 1:length(sig_list)){
    A_list[[i]] <- solve(sig_list[[i]][idx_in,idx_in])
    M_list[[i]] <- gen_m(X_list[[i]][idx_in,],A_list[[i]])
    u_list[[i]] <- t(C_list[[i]]) %*% solve(M_list[[i]]) %*% t(X_list[[i]])
  }

  out_list <- GradRobustStep(idx_in -1, do.call(rbind, C_list), do.call(rbind, X_list), do.call(rbind, sig_list), weights = w)
  idx_in <- out_list[["idx_in"]] + 1
  idx_out <- out_list[["idx_out"]] + 1
  best_val_vec <- out_list[["best_val_vec"]]

  # as a check, see if the next swap would mean that M is not positive semi-definite, in which
  # case the algorithm has terminated early and the solution is not in this design
  # this code is copied from the choose_swap_robust function, so could be streamlined into its own function perhaps

  idx_test <- c(idx_in[2:n],idx_out[1])
  val_out_mat <- matrix(NA,nrow=length(idx_test),ncol=length(A_list))
  for(idx in 1:length(A_list)){
    val_out_mat[,idx] <- sapply(1:length(idx_test),function(i)
      remove_one(A_list[[idx]],i-1,u_list[[idx]][idx_test]))
  }
  val_out <- as.numeric(val_out_mat %*% w)
  idx_test <- idx_test[-which.max(val_out)]
  rm1A <- list()
  for(idx in 1:length(A_list)){
    rm1A[[idx]] <- remove_one_mat(A_list[[idx]],which.max(val_out)-1)
  }
  for(j in 1:length(sig_list))
  {
    M_list[[j]] <- gen_m(X_list[[j]][idx_test,],rm1A[[j]])

    if(!check_psd(M_list[[j]]))stop(paste0("M not positive semi-definite. Column ",which(colSums(M_list[[j]])==0),
                                           " of design ",j," is not part of an optimal design."))
  }

  ## if columns were removed then return the index to the indexing of the original X matrix
  if(!is.null(rm_cols))
  {
    idx_in <- idx_original[idx_in]
  }

  #return variance
  if(trace){
    cat("\nVariance for individual model(s):\n")
    print(c(best_val_vec))
    if(length(A_list)>1){
      cat("\n Weighted average variance:\n")
      print(sum(best_val_vec*c(w)))
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
