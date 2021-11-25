fastMCD <- function(X, h = 0){
    n <- nrow(X); p <- ncol(X)
    if(!h){
        h <- (n + p + 1) / 2
    }
    if(h == n){
        mu <- colMeans(X)
        sigma2 <- cov(X)
    } else{
        res <- if (n <= 600) {smallMCD(X, h)} else {bigMCD(X,h)}
        d <- d2(X, res$T, res$S)
        S_mcd <- median(d) * res$S / qchisq(0.5, p)
        w <- d2(X, res$T, S_mcd) <= qchisq(0.975, p)
        mu <- colSums(X * w) / sum(w)
        sigma2 <- (t(X) %*% (X * w)) / (sum(w) - 1)
    }
    return(list(mu = mu, sigma2 = sigma2))
}

smallMCD <- function(X, h){
    H_all <- sapply(1:500, function(i) {draw_h(X, h)})
    H10 <- pick10(X, H_all, h)$H
    res10 <- apply(H10, 2, function(H) {
        S <- cov(X[H,])
        T <- colMeans(X[H,])
        step_it(X, T, S, h)
    })
    S_det <- sapply(res10, function(x) x$det_S)
    res <- res10[[which.min(S_det)]]
    return(list(T = colMeans(X[res$H,]), S = cov(X[res$H,])))
}

bigMCD <- function(X, h, p){
    k <- min(5, ceiling(n / 300))
    n_merge  <- min(1500, n)
    i_start <- rep(n_merge %% k, k)
    n_sub <- floor(n_merge / k)
    h_sub <- floor(n_sub * h / n)
    h_merge <- floor(n_merge * h / n)
    if(!(n_merge %% k)){
        i_start <- n_sub * (0:(k - 1)) + 1
        i_end <- i_start - 1 + n_merge / k
    } else{
        i_start[seq(n_merge %% k)] <- seq(n_merge %% k) - 1
        i_start <- n_sub * (0:(k - 1)) + i_start + 1
        i_end <- c(i_start[2:k] - 1, n_merge)            
    }
    samp <- sample(n)
    ind_H <- lapply(1:k, function(i) samp[i_start[i]:i_end[i]])
    X_sub <- lapply(1:k, function(i) {X[ind_H[[i]],]})
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
    T_merge <- S_merge <- det_merge <- NULL
    for (i in 1:ncol(T_sub)){
        res <- step_it(X_merge, T_sub[,i], S_sub[[i]], h_merge, 2)
        T_merge <- cbind(T_merge, res$T)
        S_merge <- append(S_merge, list(res$S))
        det_merge <- c(det_merge, res$det_S)
    }
    T_merge <- T_merge[,head(order(det_merge), 10)]
    S_merge <- S_merge[head(order(det_merge), 10)]
    res <- NULL
    for (i in 1:10){
        res <- append(res, list(step_it(X, T_merge[,i], S_merge[[i]], h)))
    }
    imin <- which.min(sapply(1:10, function(i) res[[i]]$det_S))
    return(list(S = res[[imin]]$S, T = res[[imin]]$T))
}

pick10 <- function(X, H_all, h){
    res <- apply(H_all, 2, function(H) {
        S <- cov(X[H,])
        T <- colMeans(X[H,]) 
        step_it(X, T, S, h, 2)
    })
    H <- sapply(res, function(x) x$H)
    det_S <- sapply(res, function(x) x$det_S)
    S <- lapply(res, function(x) x$S)
    T <- sapply(res, function(x) x$T)
    return(list(H = H[, head(order(det_S), 10)], S = S[head(order(det_S), 10)], T = T[, head(order(det_S), 10)]))
}

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

cstep <- function(X, T, S, h){
    d <- d2(X, T, S)
    H2 <- head(order(d), h)
    S2 <- cov(X[H2,])
    T2 <- colMeans(X[H2,])
    return(list(H = H2, S = S2, det_S = det(S2), T = T2))
}

draw_h <- function(X, h){
    p <- ncol(X); n <- nrow(X)
    j <- sample(1:n, p + 1)
    mu <- colMeans(X[j,])
    s <- cov(X[j,])
    while(det(s) == 0 & length(j) < n - 1){
        j <- c(j, sample((1:n)[-j], 1))
        mu <- colMeans(X[j,])
        s <- cov(X[j,])      
    }
    d <- d1(X, mu, s)
    H1 <- head(order(d), h)
    return(H1)
}


d2 <- function(X, T, S){
    S_inv <- solve(S)
    d <- apply(X, 1, function(x) {
        (x - T) %*% S_inv %*% (x - T)
    })    
}
