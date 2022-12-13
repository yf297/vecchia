#' Decompresses a sparse representation of a lower triangular matrix
#'
#' @param Linv The Vecchia approximation in compressed format
#' @param NNarray A matrix of indices of nonzero entries
#' @return Li The Vecchia approximation in non-compressed format
#' @export
expandm <- function(Linv, NNarray){
  n <- nrow(NNarray)
  m <- ncol(NNarray)-1
  Li <- matrix(0,n,n)

  for( i in 1:n){
    b  = min(i,m+1)
    ind <-
      Li[i,NNarray[i, !is.na(NNarray[i,]) ]] <- Linv[i, !is.na(NNarray[i,]) ]
  }
  return(Li)
}

new_ch <- function(S11, L11, k){
  new_ch <- rbind(cbind(L11, rep(0,k-1)), rep(0,k))
  r <- solve(new_ch[1:(k-1),1:(k-1)], S11[1:(k-1), (k)])
  new_ch[k,] <- c(r, (S11[k, k] - sum(r*r))^(1/2) )
  new_ch
}


remove_v <- function(win_j, array){
  array <-  array[-which(array==win_j)]
}



#' Uses a dynamic algorithm to construct the Vecchia approximation
#'
#' @param covfun The covariance function used to create covmat
#' @param covparms The covariance parameters used to create covmat
#' @param covmat The covariance matrix
#' @param NNarray2 A matrix of indices to search when doing dynamic search
#'  @return A list containing
#' \itemize{
#'   \item NNarray A matrix of indices of nonzero entries
#'   \item Linv - The Vecchia Approximation
#' }
#' @export
FSPAI <- function(covfun, covparms,covmat, NNarray2, locs, m){

  n <- nrow(NNarray2)
  m2 <-  ncol(NNarray2)-1

  NNarray <- matrix(NA,n,m+1)
  NNarray[,1] <- 1:n

  tau <- rep(0,n)
  Linv <- matrix(NA,n,m+1)
  Linv[,1] <- (covmat[,1])^(-1/2)

  for(i in 2:n){

    b = min(i,m+1)
    b2 = min(i,m2+1)
    array <- NNarray2[i,2:b2]

    k = 1
    ind <- 1:k

    for(j in array){
      tau[j] <- ( sum( covmat[j,NNarray[i,ind]] * Linv[i,ind] )^2 )/covmat[j,j]
    }

    if(sum(tau)==0)break

    win_j <- which.max(tau)

    NNarray[i,k+1] <- win_j
    array <- remove_v(win_j, array)

    loc_sub <- locs[NNarray[i,c(2:(k+1),1)],]
    sub_cov_all <- covfun(covparms, loc_sub)

    S11 <- sub_cov_all[1:k,1:k]
    S12 <- sub_cov_all[(k+1),1:k]
    L11 = matrix((covmat[i,i])^(1/2),1,1)
    ls <- solve(t(L11),solve(L11, S12))
    lii <- as.numeric(1/sqrt(sub_cov_all[(k+1),(k+1)] - S12%*%ls))
    Linv[i,1:(k+1)] =  c(lii,-lii*ls)

    tau <- rep(0,n)

    for(k in 2:(b-1)){
      ind <- 1:k

      for(j in array){
        tau[j] <- ( sum( covmat[j,NNarray[i,ind]] * Linv[i,ind] )^2 )/covmat[j,j]
      }

      if(sum(tau)==0)break

      win_j <- which.max(tau)
      NNarray[i,k+1] <- win_j
      array <- remove_v(win_j, array)

      loc_sub <- locs[NNarray[i,c(2:(k+1),1)],]
      sub_cov_all <- covfun(covparms, loc_sub)

      S11 <- sub_cov_all[1:k,1:k]
      S12 <- sub_cov_all[(k+1),1:k]
      L11 <- new_ch(S11,L11,k)
      ls <- backsolve(t(L11),forwardsolve(L11,S12))
      lii <- as.numeric(1/sqrt(sub_cov_all[(k+1),(k+1)] - S12%*%ls))
      Linv[i,1:(k+1)] =  c(lii,-lii*ls)

      tau <- rep(0,n)
    }
  }

  return(list(NNarray = NNarray, Linv = Linv))

}
