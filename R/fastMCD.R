#' fastMCD
#' Estimate location and scatter using FAST-MCD algorithm
#' @param X A 2D matrix to estimate location and scatter from
#' @param h An integer optionally specifying number of observations to use
#' @return A list of estimated location (center) and scatter (cov)
#' @export
#' @importFrom stats mahalanobis median qchisq cov
#' @import utils
#' @examples
#' set.seed(51234)
#' S <- matrix(runif(5^2), 5)
#' S <- t(S) %*% S
#' X <- MASS::mvrnorm(700, mu = rep(0, 5), Sigma = S) # generate random matrix
#' outliers <- MASS::mvrnorm(140, mu = rep(5, 5), Sigma = S)
#' X[seq(1, 700, 5), ] <- outliers # set 20% of observations to be outliers 
#' res <- fastMCD(X) # estimate location and scatter
fastMCD <- function(X, h = 0){
    n <- nrow(X); p <- ncol(X)
                                        # set h
    if (!h) {
        h <- (n + p + 1) / 2
    }
                                        # call helper for big or small data
    if (h == n) {
        mu <- colMeans(X)
        sigma2 <- cov(X)
    } else {
        if (n <= 600) {
            res <- smallMCD(X, h)
        } else {
            res <- bigMCD(X, h, p, n)
        }
                                        # reweight estimates
        d2 <- mahalanobis(X, center = res$T, cov = res$S)
        S_mcd <- median(d2) * res$S / qchisq(p = 0.5, df = p)
        w <- mahalanobis(X, center = res$T, cov = S_mcd)
        w <- w <= qchisq(p = 0.975, df = p)
        T <- colSums(X * w) / sum(w)
        S <- (t(X) %*% (X * w)) / (sum(w) - 1)
    }
    return(list(center = T, cov = S))
}

#' Obtain unweighted estimates for data with <= 600 observations
#' @param X A 2D matrix to estimate location and scatter from
#' @param h An integer specifying number of observations to use
#' @return A list of estimated location (center) and scatter (cov)
smallMCD <- function(X, h){
                                        # Sample intial subsets
    H_all <- sapply(1:500, function(i) {
        draw_h(X, h)
    })
                                        # Iterate for 10 best subsets
    H10 <- pick10(X, H_all, h)$H
    res10 <- apply(H10, 2, function(H) {
        S <- cov(X[H, ])
        T <- colMeans(X[H, ])
        step_it(X, T, S, h)
    })
    S_det <- sapply(res10, function(x) {
        x$det_S
        })
    res <- res10[[which.min(S_det)]]
    return(list(T = colMeans(X[res$H,]), S = cov(X[res$H,])))
}

#' Obtain unweighted estimates for data with > 600 observations
#' @param X A 2D matrix to estimate location and scatter from
#' @param h An integer specifying number of observations to use
#' @param p An integer specifying the number of columns in X
#' @param n An integer specifying the number of total observations
#' @return A list of estimated location (center) and scatter (cov)
bigMCD <- function(X, h, p, n){
    k <- min(5, ceiling(n / 300))
    n_merge  <- min(1500, n)
    i_start <- rep(n_merge %% k, k)
    n_sub <- floor(n_merge / k)
    h_sub <- floor(n_sub * h / n)
    h_merge <- floor(n_merge * h / n)
                                        # Set indices to partition data
    if(!(n_merge %% k)){
        i_start <- n_sub * (0:(k - 1)) + 1
        i_end <- i_start - 1 + n_merge / k
    } else{
        i_start[seq(n_merge %% k)] <- seq(n_merge %% k) - 1
        i_start <- n_sub * (0:(k - 1)) + i_start + 1
        i_end <- c(i_start[2:k] - 1, n_merge)            
    }
                                        # Permute data and obtain sub-partitions
    samp <- sample(n)
    ind_H <- lapply(1:k, function(i) {
        samp[i_start[i]:i_end[i]]
        })
    X_sub <- lapply(1:k, function(i) {
        X[ind_H[[i]], ]
    })
                                        # Pick 10 best in each sub-partition
    X_merge <- do.call('rbind', X_sub)
    T_sub <- S_sub <- NULL
    for (i in 1:k){
        H_all <- sapply(1:(500 / k), function(ind) {
            draw_h(X_sub[[i]], h_sub)
        })
        res10 <- pick10(X_sub[[1]], H_all, h_sub)
        T_sub <- cbind(T_sub, res10$T)
        S_sub <- append(S_sub, res10$S)
    }
                                        # Perform 2 C-steps and select 10 best
    T_merge <- S_merge <- det_merge <- NULL
    for (i in 1:ncol(T_sub)){
        res <- step_it(X_merge, T_sub[,i], S_sub[[i]], h_merge, 2)
        T_merge <- cbind(T_merge, res$T)
        S_merge <- append(S_merge, list(res$S))
        det_merge <- c(det_merge, res$det_S)
    }
    T_merge <- T_merge[,head(order(det_merge), 10)]
    S_merge <- S_merge[head(order(det_merge), 10)]
                                        # C-step until convergence for best subsets
    res <- NULL
    for (i in 1:10){
        res <- append(res, list(step_it(X, T_merge[, i], S_merge[[i]], h)))
    }
    imin <- which.min(sapply(1:10, function(i) res[[i]]$det_S))
    return(list(S = res[[imin]]$S, T = res[[imin]]$T))
}

#' Choose the 10 best estimates after iterating twice through initial sets
#' @param X A 2D matrix
#' @param H_all A 2D matrix where each row specifies a subset of observations
#' @param h An integer specifying number of observations to use
#' @return A list of best sets (H), scatter (S) and location (T)
pick10 <- function(X, H_all, h){
    res <- apply(H_all, 2, function(H) {
        S <- cov(X[H, ])
        T <- colMeans(X[H, ]) 
        step_it(X, T, S, h, 2)
    })
    H <- sapply(res, function(x) {
        x$H
        })
    det_S <- sapply(res, function(x) {
        x$det_S
        })
    S <- lapply(res, function(x) {
        x$S
        })
    T <- sapply(res, function(x) {
        x$T
        })
    return(list(H = H[, head(order(det_S), 10)], S = S[head(order(det_S), 10)], T = T[, head(order(det_S), 10)]))
}

#' Iterate through C-step
#' @param X A 2D matrix
#' @param T A vector of the initial location estimate
#' @param S A vector of the initial scatter estimate
#' @param h An integer specifying the number of observations to use
#' @param it An optional integer specifying the number of C-steps to perform.
#' With it = 0, C-step will be performed until convergence
#' @return A list of set (H), scatter (S) and location (T)
step_it <- function(X, T, S, h, it = 0){
   if(!it){
        det_old <- Inf
        det_new <- det(S)
        while(det_old > det_new & det_new > 0){
            res <- cstep(X, T, S, h)
            S <- res$S
            T <- res$T
            det_old <- det_new
            det_new <- res$det_S
        }
    } else{
        for (i in 1:it){
            res <- cstep(X, T, S, h)
            S <- res$S
            T <- res$T            
        }        
    }
    return(res)
}

#' Perform single iteration of C-step
#' @param X A 2D matrix
#' @param T A vector of the initial location estimate
#' @param S A vector of the initial scatter estimate
#' @param h An integer specifying the number of observations to use
#' @return A list of set (H), scatter (S) and location (T)
cstep <- function(X, T, S, h){
    d2 <- mahalanobis(X, center = T, cov = S)
    H2 <- head(order(d2), h)
    S2 <- cov(X[H2, ])
    T2 <- colMeans(X[H2, ])
    return(list(H = H2, S = S2, det_S = det(S2), T = T2))
}

#' Randomly draw a subset of observations
#' @param X A 2D matrix
#' @param h An integer specifying the number of observations to use
#' @return A vector representing an h-length subset of X
draw_h <- function(X, h){
    p <- ncol(X); n <- nrow(X)
    j <- sample(1:n, p + 1)
    T <- colMeans(X[j, ])
    S <- cov(X[j,])
    while(det(S) == 0 & length(j) < n - 1){
        j <- c(j, sample((1:n)[-j], 1))
        T <- colMeans(X[j, ])
        S <- cov(X[j, ])      
    }
    d2 <- mahalanobis(X, center = T, cov = S)
    H1 <- head(order(d2), h)
    return(H1)
}
