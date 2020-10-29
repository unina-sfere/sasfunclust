spfclust <-
  function(x=NULL,curve=NULL,timeindex=NULL,data=NULL, q = 5, 
           der = 2, lambda_s = 1e1, lambda_l = 1e1, gamma_ada=1,
           K = 2, tol = 10^-7, maxit = 15, pert =  
             NA, grid = seq(0, 1, length = 100),par_LQA=list(eps_diff = 1e-06,MAX_iter_LQA=200,eps_LQA = 1e-05),
           CK=log(20),hard = F, plot= F,trace=F,perc_rankpapp=NA,mu_fd_vero=NA,lambda_s_ini=NA,lambda_s_pen=NA,init="kmeans",breaks=NA,varcon="non-diagonal")
  {
    # This is the main function to implement the FClust procedure.
    
    if (is.null(data))
      data <- list(x=x,curve=curve,timeindex=timeindex)
    
    if(K==1)lambda_l=0
    # Initialize the parameters
    initfit <- fclustinit(data = data, pert = pert, grid = grid,  q = q, K = K,der=der,gamma_ada=gamma_ada,mu_fd_vero=mu_fd_vero,lambda_s_ini=lambda_s_ini,lambda_s_pen=lambda_s_pen,init=init,breaks=breaks,varcon=varcon)
    
    parameters <- initfit$parameters
    vars <- initfit$vars
    S <- initfit$S
    FullS <- initfit$FullS
    W <- initfit$W
    AW_vec<-initfit$AW_vec
    P_tot<- initfit$P_tot
    P <-initfit$P
    basis<- initfit$basis
    sigma.old <- 0
    sigma.new <- parameters$sigma
    ind <- 1
    if(plot){
      cluster.mean <- matrix(0,K,dim(FullS)[1])
      for(k in 1:K)
        cluster.mean[k,] <- FullS %*% (t(parameters$mu)[,k])
      plot(grid,grid,ylim=range(cluster.mean),type='n',ylab="Cluster Mean")
      for(k in 1:K)
        lines(grid,cluster.mean[k,], col = 4, lwd = 2)
    }
    # Main loop. Iterates between M and E steps and stops when
    # sigma  has converged.
    lk_old=0
    lk_new=loglik (parameters, data, vars, S,W,AW_vec,P_tot,lambda_s,lambda_l,CK)
    if (trace)    print(paste("Iteration", ind,": Sigma = ",sigma.new," loglk = ",lk_new[1]," ploglk = ",lk_new[2]))
    while(abs(lk_old - lk_new) > tol & (ind <= maxit)) {
      parameters <- fclustMstep(parameters, data, vars, S, tol, hard,lambda_s,lambda_l,W,AW_vec,P_tot,par_LQA,CK,perc_rankpapp,varcon=varcon)
      vars <- fclustEstep(parameters, data, vars, S, hard)
      lk_old<-lk_new
      lk_i<-loglik (parameters, data, vars, S,W,AW_vec,P_tot,lambda_s,lambda_l,CK)
      lk_new<--lk_i[2]
      sigma.old <- sigma.new
      sigma.new <- parameters$sigma[1]
      
      if (trace)
        # print(paste("Iteration", ind,": Sigma = ",sigma.new," loglk = ",lk_i[1]," ploglk = ",lk_i[2]))
        print(paste("Iteration", ind,": Sigma = ",sigma.new," loglk = ",lk_i[1]," ploglk = ",lk_i[2]))
      #Plot cluster mean curves.
      if(plot){
        cluster.mean <- matrix(0,K,dim(FullS)[1])
        for(k in 1:K)
          cluster.mean[k,] <- FullS %*% (t(parameters$mu)[,k])
        plot(grid,grid,ylim=range(cluster.mean),type='n',ylab="Cluster Mean")
        for(k in 1:K)
          lines(grid,cluster.mean[k,], col = 4, lwd = 2)
      }
      ind <- ind + 1
    }
    list(data=data,parameters = parameters, vars = vars, FullS = FullS,grid=grid,S=S,
         W=W,AW_vec=AW_vec,P_tot=P_tot,P=P,lambda_s=lambda_s,lambda_l=lambda_l,CK=CK,basis=basis)
  }

fclustinit <-
  function(data, pert = 0, grid = seq(0.01, 1, length = 100),  q = 5,K = K,der=der,gamma_ada=gamma_ada,mu_fd_vero=mu_fd_vero,lambda_s_ini=lambda_s_ini,lambda_s_pen=lambda_s_pen,init=init,breaks=breaks,varcon=varcon){
    S <- FullS <- NULL
    # This function initializes all the parameters.
    # Produce spline basis matrix
    basis<-if(is.na(breaks))create.bspline.basis(c(grid[1],grid[length(grid)]),nbasis = q) else create.bspline.basis(c(grid[1],grid[length(grid)]),nbasis = q,breaks = breaks )
    FullS<-eval.basis(grid,basis)
    W<-eval.penalty(basis,der)
    S <- FullS[data$timeindex,  ]
    order<-q-length(basis$params)
    breaks<-basis$params
    
    ext_break<-c(rep(grid[1],order),breaks,rep(grid[length(grid)],order))
    weights_vec<-rep(diff(ext_break,lag=order)/order,each=K)
    
    # Get pairwise matrix
    if(K!=1)
    {
      P<-matrix(0,((K-1)^2+(K-1))/2,K)
      ind<-c(1,K-1)
      for (ii in 1:(K-1)) {
        
        P[ind[1]:ind[2],ii]<-rep(1,length(ind[1]:ind[2]))
        if(length(ind[1]:ind[2])==1)
          aa<--1
        else
          aa<-diag(rep(-1,length(ind[1]:ind[2])))
        
        P[ind[1]:ind[2],(ii+1):K]<-aa
        ind<-ind+c((K-1)-(ii-1),(K-1)-(ii))
        ind<-c(min(((K-1)^2+(K-1))/2,ind[1]),min(((K-1)^2+(K-1))/2,ind[2]))
      }
      
      P_tot<-matrix(0,q*((K-1)^2+(K-1))/2,K*q)
      for (ii in 1:K) {
        P_tot[, ((ii - 1) * q + 1):((ii) * q)]<-kronecker(diag(q), P[, ii])
      }
      P_tot<-Matrix(P_tot,sparse = TRUE)
    }
    else{
      P<-P_tot<-diag(q)
    }
    
    # Weight approximation L1 penalty
    order<-q-length(basis$params)
    breaks<-basis$params
    
    ext_break<-c(rep(grid[1],order),breaks,rep(grid[length(grid)],order))
    weights_vec<-rep(diff(ext_break,lag=order)/order,each=((K-1)^2+(K-1))/2)
    N <- length(table(data$curve))
    
    # Compute initial estimate of basis coefficients.
    if(!is.na(pert)){
      # print(paste("pert=",pert))
      points <- matrix(0,N,sum(q))
      for (i in 1:N){
        Si <- S[data$curve==i,]
        xi <- data$x[data$curve==i]
        points[i,] <- solve(t(Si) %*% Si + pert * diag(q)) %*% t(Si) %*%xi
      }
    }
    else{
      # print(paste("pert=NA"))
      d<-sapply(1:N,function(i)length(which(data$curve==i)))
      e<-lapply(1:N, function(i)data$timeindex[data$curve==i])
      
      if(length(unique(d))==1&length(unique(e))==1){
        basis_start<-basis
        # plot(basis_start)
        # points(grid,rep(0,length(grid)))
        grid_i <- grid[data$timeindex[data$curve==1]]
        X<-sapply(1:N,function(i) data$x[data$curve==i])
        loglam         = seq(-10, 10, 0.25)
        Gcvsave        = numeric()
        for(i in 1:length(loglam)){
          fdPari  = fdPar(basis_start, Lfdobj=2, 10^loglam[i])
          Sm.i    = smooth.basis(grid_i, X, fdPari)
          Gcvsave[i] = sum(Sm.i$gcv)
          
        }
        lambda_s=if(is.na(lambda_s_ini))10^loglam[which.min(Gcvsave)]else lambda_s_ini
        # print(lambda_s)
        fdPari  = fdPar(basis_start, Lfdobj=2,lambda_s)
        points<-t(smooth.basis(grid_i, X, fdPari)$fd$coefs)
        # print(lambda_s)
        # plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
        #      ylab=expression(GCV(lambda)), lwd=2 )
        # plot(smooth.basis(grid_i, X, fdPari)$fd)
      }
      else{
        basis_start<-basis
        loglam         = seq(-10, 10, 0.25)
        Gcvsave        = matrix(0,N,length(loglam))
        points <- matrix(0,N,sum(q))
        for (i in 1:N){
          print(i)
          xi <- data$x[data$curve==i]
          grid_i <- grid[data$timeindex[data$curve==i]]
          
          for(l in 1:length(loglam)){
            fdPari  = fdPar(basis_start, Lfdobj=2, 10^loglam[l])
            Sm.i    = smooth.basis(grid_i, xi, fdPari)
            Gcvsave[i,l] = sum(Sm.i$gcv)
          }
          lambda_s[i]=10^loglam[which.min(Gcvsave[i,])]
          fdPari  = fdPar(basis_start, Lfdobj=2,lambda_s[i])
          points[i,]<-t(smooth.basis(grid_i, xi, fdPari)$fd$coefs)
          
          # plot(loglam,Gcvsave[i,], 'o', las=1, xlab=expression(log[10](lambda)),
          #      ylab=expression(GCV(lambda)), lwd=2 )
          # plot(smooth.basis(grid_i, xi, fdPari)$fd)
        }
      }
    }
    
    # Use k-means to get initial cluster memberships from points.
    if(K > 1){
      if(init=="kmeans")class <-  kmeans(points, K, 1000,nstart = 10000)$cluster 
      else if(init=="model-based") class <- Mclust(points,K,verbose=FALSE,modelNames = "EII")$classification
      else if(init=="hierarchical")class <-  hcut(points,k =K )$cluster
    }
    else {
      class <- rep(1, N)
    }
    # matplot(t(points),type="l",col=class)
    # Initial estimates for the posterior probs of cluster membership.
    piigivej <- matrix(0, N, K)
    piigivej[col(piigivej) == class] <- 1
    pi=apply(piigivej,2,mean)
    # 
    # # Calculate coefficeints for cluster means.
    classmean <- matrix(0,K,q)
    for (k in 1:K)
      classmean[k,] <- apply(as.matrix(points[class==k,]),2,mean)
    
    
    #  X<-t(sapply(1:N,function(i) data$x[data$curve==i]))
    #  mu_eval<-matrix(0,dim(X)[2],K)
    # new_X<-sapply(1:K,function(i)if(is.null(dim(X[class==i,]))) X[class==i,] else apply(X[class==i,],2,mean))
    #  if(!is.na(mu_fd_vero))new_X<-eval.fd(grid,mu_fd_vero)
    #  basis<- create.bspline.basis(c(grid[1],grid[length(grid)]), nbasis=q)
    #  loglam         = seq(-10, 10, 0.25)
    #  Gcvsave        = numeric()
    #  for(i in 1:length(loglam)){
    #    fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
    #    Sm.i    = smooth.basis(grid, new_X, fdPari)
    #    Gcvsave[i] = sum(Sm.i$gcv)
    # 
    #  }
    #  lambda_s=if(is.na(lambda_s_pen))10^loglam[which.min(Gcvsave)]else lambda_s_pen
    #  # print(lambda_s)
    #  fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
    #  mu_fd<-smooth.basis(grid, new_X, fdPari)$fd
    #  mu_start=t(mu_fd$coefs)
    
    
    
    # Initialize mu.
    mu<-mu_start<-classmean
    gamma<-array(0, c(N, K, sum(q)))
    gprod <- NULL
    if(K==1){
      for(i in 1:N){
        gamma[i,,]<-t(points[i,] - t(mu))
        gprod <- cbind(gprod, (gamma[i,  ,]) %*% t(gamma[i, , ]))
      }
    }
    else{
      for(i in 1:N){
        gamma[i,,]<-t(points[i,] - t(mu))
        gprod <- cbind(gprod, t(gamma[i,  ,]) %*% (gamma[i, , ]))
      }
    }
    N <- dim(gamma)[1]
    ind <- matrix(rep(c(rep(c(1, rep(0, q - 1)), N), 0), q)[1:(N*q^2)], N * q, q)
    Gamma <-gprod %*% ind/N
    if(varcon=="diagonal")Gamma <-diag(diag( Gamma))
    if(varcon=="equal")Gamma=diag(q)*sum(diag( Gamma))/(q)
    #   
    #   X_fd<-fd(points,basis)
    # X_eval<-eval.fd(grid,X_fd)
    
    #  print(points[10,])
    # print( gamma[10,,])
    # matplot(t(gamma[10,,]),type="l")
    
    #  X<-t(sapply(1:N,function(i) data$x[data$curve==i]))
    #  mu_eval<-matrix(0,dim(X)[2],K)
    # new_X<-sapply(1:K,function(i)if(is.null(dim(X[class==i,]))) X[class==i,] else apply(X[class==i,],2,mean))
    #  if(!is.na(mu_fd_vero))new_X<-eval.fd(grid,mu_fd_vero)
    #  basis<- create.bspline.basis(c(grid[1],grid[length(grid)]), nbasis=q)
    #  loglam         = seq(-10, 10, 0.25)
    #  Gcvsave        = numeric()
    #  for(i in 1:length(loglam)){
    #    fdPari  = fdPar(basis, Lfdobj=2, 10^loglam[i])
    #    Sm.i    = smooth.basis(grid, new_X, fdPari)
    #    Gcvsave[i] = sum(Sm.i$gcv)
    # 
    #  }
    #  lambda_s=if(is.na(lambda_s_pen))10^loglam[which.min(Gcvsave)]else lambda_s_pen
    #  # print(lambda_s)
    #  fdPari  = fdPar(basis, Lfdobj=2,lambda_s)
    #  mu_fd<-smooth.basis(grid, new_X, fdPari)$fd
    #  mu_start=t(mu_fd$coefs)
    # 
    # plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
    #      ylab=expression(GCV(lambda)), lwd=2 )
    
    gcov <- matrix(0, sum(q), N * sum(q))
    
    # Get weight matrix
    if(K!=1){
      AW<-matrix(0,((K-1)^2+(K-1))/2,q)
      for (ii in 1:q) {
        
        AW[,ii]<-1/(abs(P%*%mu_start[,ii])^gamma_ada)
      }
      AW_vec<-vec(AW)*(weights_vec)
    }
    else{
      AW_vec<-rep(1,q)
    }
    # matplot(t(classmean),type="l")
    n <- length(data$curve)
    n_i<-sapply(1:N,function(ii)length(which(data$curve==ii)))
    sigma=as.numeric(get_sigma( data$x, data$curve, data$time,  as.matrix(S),  piigivej,  gcov,n_i,gamma,mu))/n
    
    list(S = S, W = W, AW_vec=AW_vec,P_tot=P_tot,P=P, FullS = FullS, parameters = list(mu=mu_start,sigma=sigma,pi=pi,Gamma=Gamma), vars = list(gamma = gamma,piigivej = piigivej,
                                                                                                                                               gprod = gprod, gcov = gcov ),basis=basis)
  }

fclustMstep <-
  function(parameters, data, vars, S, tol,  hard,lambda_s,lambda_l,W,AW_vec,P_tot,par_LQA,CK,perc_rankpapp,varcon=varcon)
  {
    # This function implements the M step of the EM algorithm.
    K <- dim(parameters$mu)[1]
    mu<-parameters$mu
    gamma <- vars$gamma
    gcov <- vars$gcov
    curve <- data$curve
    piigivej <- vars$piigivej
    N <- dim(gamma)[1]
    K <- dim(mu)[1]
    n <- length(curve)
    q <- dim(S)[2]
    # Compute pi.
    if(hard)
      parameters$pi <- rep(1/K, K)
    else parameters$pi <- (apply(vars$piigivej, 2, mean)*N+CK)/(N+K*CK)
    # Compute rank p estimate of Gamma
    ind <- matrix(rep(c(rep(c(1, rep(0, q - 1)), N), 0), q)[1:(N*q^2)], N * q, q)
    if(!is.na(perc_rankpapp)){
      gsvd <- svd(vars$gprod %*% ind/N)
      p<-which(cumsum( gsvd$d)/sum( gsvd$d)>=perc_rankpapp)[1]
      gsvd$d[ - (1:p)] <- 0
      parameters$Gamma <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u)
    }
    else{
      parameters$Gamma <-vars$gprod %*% ind/N
    }
    if(varcon=="diagonal")parameters$Gamma <-diag(diag( parameters$Gamma))
    if(varcon=="equal")parameters$Gamma=diag(q)*sum(diag(parameters$Gamma))/(q)
    # if(varcon=="diagonal-diff"){
    #   clus<-apply(vars$piigivej,1,which.max)
    #   diag_big_mat<-do.call(bdiag ,lapply(1:N,function(ii)diag(q)*clus[ii]))
    #   mat_list_w<-list()
    #  for (ii in 1:K) {
    #    diag_big_mat_pr<-diag_big_mat
    #    diag_big_mat_pr[which(diag_big_mat_pr!=ii)]=0
    #    diag_big_mat_pr[which(diag_big_mat_pr==ii)]=1
    #    ni<-length(which(clus==ii))
    #    mat_list_w[[ii]]<-vars$gprod%*%diag_big_mat_pr %*% ind/(ni*q)
    #    
    #  } 
    # }
    # parameters$Gamma<-tr(  parameters$Gamma)/q*diag(q)
    # Gamma[row(Gamma) - col(Gamma) == 1] <- Gamma[row(Gamma) - col(Gamma) == -1] <- parameters$Gamma[row(Gamma) - col(Gamma) == 1]
    # parameters$Gamma=Gamma
    # This loop iteratively LDA to get mu
    # when they have converged.
    W_star<-matrix(0,q*K,q*K)
    for (i in 1:K) {
      W_star[((i-1)*q+1):((i)*q),((i-1)*q+1):((i)*q)]<-W
    }
    W_star<-Matrix(W_star, sparse = TRUE)
    x <- data$x
    n_i<-sapply(1:N,function(ii)length(which(curve==ii)))
    numden<-get_numden( data$x, data$curve, data$time,  as.matrix(S),  piigivej,  gcov,n_i,gamma)
    VY<-Matrix(numden[[1]],sparse = TRUE)
    S.den<-Matrix(numden[[2]],sparse=TRUE)
    z_int=1
    diff_inter <- 100
    mu_old=mu
    
    while(diff_inter>par_LQA$eps_LQA && z_int<=par_LQA$MAX_iter_LQA) {
      
      mu_vec<-vec(t(mu))
      diff_mat<-abs(P_tot%*%mu_vec)
      diff_mat[diff_mat < par_LQA$eps_diff]<-par_LQA$eps_diff
      V_l<-Matrix(diag(as.numeric(AW_vec/(2*diff_mat))),sparse = TRUE)
      mu_vec<-solve(S.den*(1/parameters$sigma)+lambda_s*2*W_star+2*lambda_l*t(P_tot)%*%V_l%*%P_tot)%*%VY*(1/parameters$sigma)
      mu<-matrix(mu_vec,K,q,byrow = T)
      diff_inter<-sum(abs(mu-mu_old))/(sum(abs(mu_old)))
      mu_old<-mu
      z_int=z_int+1
    }
    
    # Calculate sigma 
    sigma<-get_sigma( data$x, data$curve, data$time,  as.matrix(S),  piigivej,  gcov,n_i,gamma,mu)
    sigma=as.numeric(sigma)/n
    parameters$mu <- mu
    parameters$sigma <- sigma
    parameters
  }

fclustEstep <-
  function(parameters, data, vars, S, hard)
  {
    # This function performs the E step of the EM algorithm by
    # calculating the expected values of gamma and gamma %*% t(gamma)
    # given the current parameter estimates.
    
    N <- dim(vars$gamma)[1]
    n_i<-sapply(1:N,function(ii)length(which(data$curve==ii)))
    parameters$sigma=as.matrix(parameters$sigma)
    parameters$pi=as.matrix(parameters$pi)
    aa<-get_Estep(parameters, data, vars, S, hard,n_i)
    vars$gamma=aa[[1]]
    vars$piigivej=aa[[2]]
    vars$gprod=aa[[3]]
    vars$gcov=aa[[4]]
    vars
  }


fclustconst <-
  function(data, parameters, vars, S)
  {
    # This function enforces the constraint (7) from the paper on the
    # parameters. This means that the alphas can be interpreted as the
    # number of standard deviations that the groups are apart etc.
    par <- parameters
    A <- t(S) %*% solve(par$sigma * diag(dim(S)[1]) + S %*% par$Gamma %*%
                          t(S)) %*% S
    svdA <- svd(A)
    sqrtA <- diag(sqrt(svdA$d)) %*% t(svdA$u)
    negsqrtA <- svdA$u %*% diag(1/sqrt(svdA$d))
    finalsvd <- svd(sqrtA %*% par$Lambda)
    par$Lambda <- negsqrtA %*% finalsvd$u
    if(dim(par$Lambda)[2] > 1)
      par$alpha <- t(diag(finalsvd$d) %*% t(finalsvd$v) %*% t(par$alpha))
    else par$alpha <- t(finalsvd$d * t(finalsvd$v) %*% t(par$alpha))
    meanalpha <- apply(par$alpha, 2, mean)
    par$alpha <- t(t(par$alpha) - meanalpha)
    par$lambda.zero <- par$lambda.zero + par$Lambda %*% meanalpha
    list(parameters = par, vars = vars)
  }

nummax <-
  function(X)
  {
    ind <- rep(1, dim(X)[1])
    m <- X[, 1]
    if(dim(X)[2] > 1)
      for(i in 2:dim(X)[2]) {
        test <- X[, i] > m
        ind[test] <- i
        m[test] <- X[test, i]
      }
    list(ind = ind, max = m)
  }

fclust.pred <-
  function(fit,data=NULL,reweight=F)
  {
    # This function produces the alpha hats used to provide a low
    # dimensional pictorial respresentation of each curve. It also
    # produces a class prediction for each curve. It takes as
    # input the fit from fldafit (for predictions on the original data)
    # or the fit and a new data set (for predictions on new data).
    if (is.null(data))
      data <- fit$data
    FullS <- fit$FullS
    par <- fit$parameters
    curve <- data$curve
    time <- data$time
    N <- length(table(curve))
    h <- dim(par$alpha)[2]
    alpha.hat <- matrix(0, N, h)
    K <- dim(fit$par$alpha)[1]
    distance <- matrix(0, N, K)
    Calpha <- array(0, c(N, h, h))
    for(i in 1:N) {
      Sij <- FullS[time[curve == i],  ]
      xij <- data$x[curve == i]
      n <- length(xij)
      Sigma <- par$sigma * diag(n) + Sij %*% par$Gamma %*% t(Sij)
      # Calculate covariance for each alpha hat.
      InvCalpha <- t(par$Lambda) %*% t(Sij) %*% solve(Sigma) %*% Sij %*%
        par$Lambda
      Calpha[i,  ,  ] <- solve(InvCalpha)
      # Calculate each alpha hat.
      alpha.hat[i,  ] <- Calpha[i,  ,  ] %*% t(par$Lambda) %*% t(
        Sij) %*% solve(Sigma) %*% (xij - Sij %*% par$lambda.zero)
      # Calculate the matrix of distances, relative to the
      # appropriate metric of each curve from each class centroid. 
      for (k in 1:K){
        y <- as.vector(alpha.hat[i,])-fit$par$alpha[k,]
        distance[i,k] <- t(y)%*%InvCalpha %*%y}}
    # Calculate final class predictions for each curve.
    class.pred <- rep(1, N)
    log.pi <- log(fit$par$pi)
    if (!reweight)
      log.pi <- rep(0,K)
    probs <- t(exp(log.pi-t(distance)/2))
    probs <- probs/apply(probs,1,sum)
    m <- probs[,1]
    if(K != 1)
      for(k in 2:K) {
        test <- (probs[, k] > m)
        class.pred[test] <- k
        m[test] <- probs[test, k]
      }
    list(Calpha = Calpha, alpha.hat = alpha.hat, class.pred = class.pred,
         distance = distance, m = m,probs=probs)
  }

fclust.curvepred <-
  function(fit, data=NULL, index=NULL, tau = 0.95, tau1 = 0.975)
  {
    if (is.null(data))
      data <- fit$data
    if (is.null(index))
      index <- 1:length(table(data$curve))
    tau2 <- tau/tau1
    sigma <- fit$par$sigma
    Gamma <- fit$par$Gamma
    mu<-fit$par$mu
    # Lambda <- fit$par$Lambda
    # alpha <- fit$par$alpha
    # lambda.zero <- as.vector(fit$par$lambda.zero)
    S <- fit$FullS
    N <- length(index)
    upci <-lowci <- uppi <- lowpi <- gpred <- matrix(0,N,nrow(S))
    etapred <- matrix(0,N,ncol(S))
    ind <- 1
    
    for (i in index){
      y <- data$x[data$curve == i]
      Si <- S[data$time[data$curve == i],  ]
      ni <- dim(Si)[1]
      invvar <- diag(1/rep(sigma, ni))
      covx <- Si %*% Gamma %*% t(Si) + solve(invvar)
      centx <- data$x[data$curve == i] - Si %*% t(mu)
      d <- exp( - diag(t(centx) %*% solve(covx) %*% centx)/2) * fit$par$pi
      pi <- d/sum(d)
      K <- length(pi)
      mu.pi <- t(mu* pi) %*% rep(1, K) 
      cov <- (Gamma - Gamma %*% t(Si) %*% solve(diag(ni) + Si %*% Gamma %*%
                                                  t(Si)/sigma) %*% Si %*% Gamma/sigma)/sigma
      etapred[ind,] <- mu.pi + cov %*% t(Si) %*% (y - Si %*% mu.pi)
      ord <- order( - pi)
      numb <- sum(cumsum(pi[ord]) <= tau1) + 1
      v <- diag(S %*% (cov * sigma) %*% t(S))
      pse <- sqrt(v + sigma)
      se <- sqrt(v)
      lower.p <- upper.p <- lower.c <- upper.c <- matrix(0, nrow(S), numb)
      for(j in 1:numb) {
        mean <- S %*% (t(mu)[,ord[j]] + cov %*% t(Si) %*% (y - Si %*% t(mu)[,ord[j]]))
        upper.p[, j] <- mean + qnorm(tau2) * pse
        lower.p[, j] <- mean - qnorm(tau2) * pse
        upper.c[, j] <- mean + qnorm(tau2) * se
        lower.c[, j] <- mean - qnorm(tau2) * se
      }
      upci[ind,] <- nummax(upper.c)$max
      lowci[ind,] <-  - nummax( - lower.c)$max
      uppi[ind,] <- nummax(upper.p)$max
      lowpi[ind,] <-  - nummax( - lower.p)$max
      gpred[ind,] <- as.vector(S %*%etapred[ind,])
      ind <- ind+1
    }
    meancurves <- S%*%t(mu)
    list(etapred = etapred, gpred = gpred,  upci = upci,lowci = lowci,  uppi = uppi, lowpi = lowpi,index=index,grid=fit$grid,data=data,meancurves=meancurves)
  }

fclust.discrim <-
  function(fit,absvalue=F){
    S <- fit$FullS
    sigma <- fit$par$sigma
    n <- nrow(S)
    Gamma <- fit$par$Gamma
    Sigma <- S%*%Gamma%*%t(S)+sigma*diag(n)
    Lambda <- fit$par$Lambda
    discrim <- solve(Sigma)%*%S%*%Lambda
    if (absvalue)
      discrim <- abs(discrim)
    n <- ncol(discrim)
    nrows <- ceiling(sqrt(n))
    par(mfrow=c(nrows,nrows))
    for (i in 1:n){
      plot(fit$grid,discrim[,i],ylab=paste("Discriminant Function ",i),xlab="Time",type='n')
      lines(fit$grid,discrim[,i],lwd=3)
      abline(0,0)}}

fclust.plotcurves <-
  function(object=NULL,fit=NULL,index=NULL,ci=T,pi=T,clustermean=F){
    if (is.null(object))
      object <- fclust.curvepred(fit)
    if (is.null(index))
      index <- 1:length(table(object$data$curve))
    r <- ceiling(sqrt(length(index)))
    # par(mfrow=c(r,r))
    for (i in index){
      grid <- object$grid
      upci <- object$upci[i,]
      uppi <- object$uppi[i,]
      lowci <- object$lowci[i,]
      lowpi <- object$lowpi[i,]
      gpred <- object$gpred[i,]
      meancurves <- (object$mean)
      if (clustermean)
        yrange <- c(min(c(lowpi,meancurves)),max(c(uppi,meancurves)))
      else
        yrange <- c(min(lowpi),max(uppi))
      plot(grid,grid,ylim=yrange,ylab="Predictions",xlab="Time",type='n',
           main=paste("Curve ",i))
      if (clustermean)
        for (k  in 1:ncol(meancurves))
          lines(grid,meancurves[,k],col=6,lty=2,lwd=2)
      if (ci){
        lines(grid,upci,col=3)
        lines(grid,lowci,col=3)}
      if (pi){
        lines(grid,uppi,col=4)
        lines(grid,lowpi,col=4)}
      lines(grid,gpred,col=2,lwd=2)
      lines(grid[object$data$time[object$data$curve==i]],object$data$x[object$data$curve==i],lwd=2)
      points(grid[object$data$time[object$data$curve==i]],object$data$x[object$data$curve==i],pch=19,cex=1.5)
    }
  }
loglik <- function(parameters, data, vars, FullS,W=NA,AW_vec=NA,P_tot=NA,lambda_s=NA,lambda_l=NA,CK=NA,CLC=FALSE, perc_rankpapp2=NA){
  
  gamma <- vars$gamma
  gcov <- vars$gcov
  curve <- data$curve
  pi <- parameters$pi
  
  S <- FullS[data$timeindex,  ]
  
  N<-length(unique(data$curve))
  K <- dim(vars$gamma)[2]
  q <- dim(vars$gamma)[3]
  
  Gamma <- parameters$Gamma
  if(!is.na(perc_rankpapp2)){
    print(det(Gamma))
    print(2)
    gsvd <- svd(Gamma)
    p<-which(cumsum( gsvd$d)/sum( gsvd$d)>=perc_rankpapp2)[1]
    gsvd$d[ - (1:p)] <- 0
    Gamma <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u)
    # Gamma=diag(diag(Gamma))
  }
  if(is.null(parameters$mu))mu<-t(matrix(parameters$lambda.zero,q,K)+parameters$Lambda%*%t(parameters$alpha))
  else mu<-parameters$mu
  
  loglk <- 0
  # print(det(Gamma))
  for(i in 1:N){
    y <- data$x[data$curve == i]
    Si <- S[data$curve == i,  ]
    ni <- dim(Si)[1]
    invvar <- diag(1/rep(parameters$sigma, ni))
    covx <- Si %*% Gamma %*% t(Si) +solve(invvar)
    covx_inv<-chol2inv(chol(covx))
    covx_det<-det(covx)
    # covx<-diag(diag(covx))
     temp<-sum(sapply(1:K,function(ll)pi[ll]*(2*base::pi)^(-ni/2)*covx_det^(-1/2)*exp(-(1/2)*t(y-Si%*%t(parameters$mu)[,ll])%*% covx_inv%*%(y-Si%*%t(parameters$mu)[,ll]))))
    #  print(temp)
    # temp<-sapply(1:K,function(ll)pi[ll]*mvtnorm::dmvnorm(t(y), mean=Si%*%t(mu)[,ll], covx))
    # print(sum(temp))
    loglk=loglk+log(sum(temp))
    # print(log(sum(temp)))
    
  }
  
  if(!is.na(lambda_l)){
    p_l=lambda_l*t(AW_vec)%*%abs(P_tot%*%vec(t(parameters$mu)))
    p_s=lambda_s*sum(sapply(1:K,function(ll)t(parameters$mu)[,ll]%*%W%*%t(parameters$mu)[,ll]))
    p_pi=CK*if(is.na(sum(sapply(1:K,function(ll)log(pi[ll]))))) is.na(sum(sapply(1:K,function(ll)log(pi[ll])))) else 0
    # print(paste(p_s," ",p_l, " ",p_pi))
    
    ploglk<-loglk-p_l-p_s+p_pi
    if(CLC==TRUE){
      EN<--sum(sapply(1:K,function(ii)sum(sapply(1:N,function(ll)vars$piigivej[ll,ii]*log(vars$piigivej[ll,ii]+10^-200)))))
      loglk=-2*loglk+2*EN
    }
    out<-round(c(loglk,ploglk[1,1]),digits = 2)
    
    
  }
  else{
    out<-loglk
  }
  return(out)
}
classify <- function(mod, data_new=NA){
  
  
  parameters<-mod$parameters
  vars<-mod$vars
  if(is.na(data_new))
    data=mod$data
  
  gamma <- vars$gamma
  gcov <- vars$gcov
  curve <- data$curve
  pi <- parameters$pi
  S <- mod$S
  N<-length(unique(data$curve))
  K <- dim(vars$gamma)[2]
  q <- dim(vars$gamma)[3]
  Gamma <- parameters$Gamma
  
  
  po_pr<-matrix(0,N,K)
  for(i in 1:N){
    y <- data$x[data$curve == i]
    Si <- S[data$time[data$curve == i],  ]
    ni <- dim(Si)[1]
    invvar <- diag(1/rep(parameters$sigma, ni))
    covx<- Si %*% Gamma %*% t(Si) + solve(invvar)
    covx_inv<-chol2inv(chol(covx))
    covx_det<-det(covx)
    # covx<-diag(diag(covx))
    temp<-sum(sapply(1:K,function(ll)pi[ll]*(2*base::pi)^(-ni/2)*covx_det^(-1/2)*exp(-(1/2)*t(y-Si%*%t(parameters$mu)[,ll])%*% covx_inv%*%(y-Si%*%t(parameters$mu)[,ll]))))
    # temp<-sapply(1:K,function(ll)pi[ll]*mvtnorm::dmvnorm(t(y), mean=Si%*%t(parameters$mu)[,ll], covx))
    po_pr[i,]=temp/sum(temp)
    
  }
  if(is.na(data_new))po_pr<-vars$piigivej
  classes<-apply(po_pr,1,which.max)
  out<-list(classes=classes,po_pr=po_pr)
  return(out)
}



spfclust_cv<-function(data,num_cluster=6,lambda_l_g=10^seq(-1,2),lambda_s_g=10^seq(-5,-3),gamma_ada=c(1),K_fold=5,test=NA,grid=NA,grid_2=NA,k1=1,k2=0,k3=1,maxit=15,q=30,
                      CK=log(20),init="kmeans",lambda_s_ini=NA,perc_rankpapp=NA,breaks=NA,varcon="full",...){
  
  
  N<-length(unique(data$curve))
  comb_list<-expand.grid(num_cluster,lambda_s_g,lambda_l_g,gamma_ada)
  # if(num_cluster[1]==1){
  # comb_list_new<-comb_list[-which(comb_list[,1]==1),]
  # comb_list_par<-expand.grid(1,lambda_s_g,0)
  # comb_list<-rbind(comb_list_par,comb_list_new)
  # }
  if(is.na(test)){
    parr_fun<-function(ii){
      
      parameters<-as.numeric(comb_list[ii,])
      
      num_clusters_i<-parameters[1]
      lambda_s_i<-parameters[2]
      lambda_l_i<-parameters[3]
      gamma_ada_i<-parameters[4]
      ran_seq<-sample(seq(1, N), N, replace=FALSE)
      split_vec<-split(ran_seq,cut(seq(1,N),breaks=K_fold,labels=FALSE))
      l_i<-zeros_vec<-numeric()
      data_fold<-data_i<-list()
      for(lll in 1:K_fold){
        
        ind_fold<-as.numeric(unlist(split_vec[-lll]))
        ind_i<-split_vec[[lll]]
        
        data_fold$x<-vec(sapply(1:length(ind_fold),function(ll)data$x[data$curve==ind_fold[ll]]))
        data_fold$timeindex<-unlist(lapply(1:length(ind_fold),function(ll)data$timeindex[data$curve==ind_fold[ll]]))
        data_fold$curve<-unlist(lapply(1:length(ind_fold),function(ll)rep(ll,length(data$timeindex[data$curve==ind_fold[ll]]))))
        data_i$x<-vec(sapply(1:length(ind_i),function(ll)data$x[data$curve==ind_i[ll]]))
        data_i$timeindex<-unlist(lapply(1:length(ind_i),function(ll)data$timeindex[data$curve==ind_i[ll]]))
        data_i$curve<-unlist(lapply(1:length(ind_i),function(ll)rep(ll,length(data$timeindex[data$curve==ind_i[ll]]))))
        
        mod<-spfclust_ss(data=data_fold,lambda_l = lambda_l_i,lambda_s =lambda_s_i,gamma_ada=gamma_ada_i,K=num_clusters_i,grid=grid,maxit=maxit,q=q,CK=CK,perc_rankpapp=perc_rankpapp,
                         init=init,lambda_s_ini=lambda_s_ini,breaks=breaks,varcon=varcon,...)
        l_i[lll]<-loglik(parameters = mod[[1]]$parameters,data = data_i,vars = mod[[1]]$vars, FullS = mod[[1]]$FullS,W = mod[[1]]$W,AW_vec = mod[[1]]$AW_vec,P_tot = mod[[1]]$P_tot,
                         lambda_s = mod[[1]]$lambda_s,CK=mod[[1]]$CK,perc_rankpapp2 = NA)[1]#/length(ind_i)
        # print(l_i[lll])
        zeros_vec[lll]<-get_zero(mod[[1]])
        # print(get_measure( clus_true[ind_i],classify(mod$mod,data =data_i )$classes))
        # print(Misclass(clus_true[ind_i],classify(mod$mod,data =data_i )$classes,best = TRUE))
      }
      mean<-mean(l_i)
      sd<-sd(l_i)/sqrt(K_fold)#
      zeros<-mean(zeros_vec)
      
      out<-list(mean=mean,
                sd=sd,
                zeros=zeros)
      return(out)
      
      
    }
  }
  else{
    
    if(!all(grid==grid_2))
      stop("Not equal grids between training and test set \n")
    
    parr_fun<-function(ii){
      
      parameters<-as.numeric(comb_list[ii,])
      num_clusters_i<-parameters[1]
      lambda_s_i<-parameters[2]
      lambda_l_i<-parameters[3]
      gamma_ada_i<-parameters[4]
      mod<-spfclust_ss(data=data,lambda_l = lambda_l_i,lambda_s =lambda_s_i,gamma_ada=gamma_ada_i,maxit=maxit,q=q,CK=CK,K=num_clusters_i,grid=grid,...)
      l_i<-loglik(mod[[1]]$parameters,test,mod[[1]]$vars, mod[[1]]$S,mod[[1]]$W,mod[[1]]$AW_vec,mod[[1]]$P_tot,mod[[1]]$lambda_s,mod[[1]]$lambda_l,CK=mod[[1]]$CK)[1]
      print(l_i)
      zeros<-get_zero(mod[[1]])
      mean<-l_i
      sd<-0
      # if(ii==1){
      #   plot(y=log10(lambda_l_g),x=log10(lambda_l_g),ylim=c(l_i-20000,l_i+20000),type='n',ylab="CV")
      #   }
      # points(x=log10(lambda_l_i),y=mean)
      out<-list(mean=mean,
                sd=sd,
                zeros=zeros)
      return(out)
      
      
    }
  }
  
  
  
  # if(is.na(test)){
  #   closeAllConnections()
  # cl<-parallel::makeCluster(12, setup_strategy = "sequential")
  # 
  # try(clusterEvalQ(cl, {
  #   
  #   library(MASS)
  #   library(fda)
  #   library(matrixcalc)
  #   library(Rcpp)
  #   Rcpp::sourceCpp('../fun_sc_cpp.cpp')
  #   source("../spfclust_sparse_par.R")
  #   source("../functions_competitors.R")}))
  # cat(22)
  # try(clusterExport(cl, c("comb_list","N","data","q","grid","maxit","K_fold","CK"),envir = environment()))
  # 
  # vec_par<-try(parLapply(cl, seq(1,length(comb_list[,1])),parr_fun))
  # print("closing clusters")
  # stopCluster(cl)
  # closeAllConnections()
  # }
  # else{
  if(is.na(test)){
    closeAllConnections()
    cl<-parallel::makeCluster(12, setup_strategy = "sequential")
    
    try(clusterEvalQ(cl, {
      
      library(MASS)
      library(fda)
      library(matrixcalc)
      library(Rcpp)
      library(mclust)
      source("../functions_competitors.R")
     Rcpp::sourceCpp('../fun_sc_cpp.cpp')
     source("../spfclust_sparse_def.R")
    }))
    cat(22)
    try(clusterExport(cl, c("comb_list","N","data","q","grid","maxit","K_fold","CK","init","lambda_s_ini","breaks","perc_rankpapp","varcon"),envir = environment()))
    
    vec_par<-try(parLapply(cl, seq(1,length(comb_list[,1])),parr_fun))
    print("closing clusters")
    stopCluster(cl)
    closeAllConnections()
  }
  else{
    cores <- detectCores()
    cat(cores)
    vec_par<-mclapply(seq(1,length(comb_list[,1])),parr_fun,mc.cores = 10)
  }
  # }
  par<-sapply(vec_par,"[[",1)
  sds<-sapply(vec_par,"[[",2)
  zeros<-sapply(vec_par,"[[",3)
  l_opt<-as.numeric(comb_list[max(which(par==max(par))),])
  
  
  # if(num_cluster[1]==1){
  #   comb_list<-expand.grid(num_cluster,lambda_s_g,lambda_l_g)
  #   par<-numeric()
  #   comb_list<-expand.grid(num_cluster,lambda_s_g,lambda_l_g)
  #   par<-
  #   comb_list[which(comb_list[,1]==1),]<-expand.grid(1,lambda_s_g,0)
  #   comb_list_par<-expand.grid(1,lambda_s_g,0)
  #   comb_list<-rbind(comb_list_par,comb_list_new)
  # }
  
  ksdrule<-get_ksdrule(par,sds,comb_list,k1,k2,k3)
  
  num_clusters_opt<-ksdrule[1]
  lambda_s_opt<-ksdrule[2]
  lambda_l_opt<-ksdrule[3]
  gamma_opt<-ksdrule[4]
  # 
  # mod_opt<-spfclust(data=data,trace=T,lambda_l = lambda_l_opt,lambda_s =lambda_s_opt,q=40,plot = TRUE,K=num_clusters_opt,...)
  # 
  out<-list(mod_opt=NA,
            num_clusters_opt=num_clusters_opt,
            lambda_l_opt=lambda_l_opt,
            lambda_s_opt=lambda_s_opt,
            gamma_opt=gamma_opt,
            comb_list=comb_list,
            CV=par,
            CV_sd=sds,
            zeros=zeros,ks=c(k1,k2,k3))
  
  return(out)
}



spfclust_ss<-function(data,...){
  
  
  
  
  mod<-spfclust(data=data,...)
  nbasis= dim(mod$S)[2]
  grid=mod$grid
  basis<- mod$basis
  mean_fd<-fd(t(mod$parameters$mu),basis)
  clus<-classify(mod)
  # 
  # mod_opt<-spfclust(data=data,trace=T,lambda_l = lambda_l_opt,lambda_s =lambda_s_opt,q=40,plot = TRUE,K=num_clusters_opt,...)
  # 
  out<-list(mod=mod,
            mean_fd=mean_fd,
            clus=clus)
  
  return(out)
}

get_ksdrule<-function(par,sds,comb_list,k1,k2,k3){
  
  lambda_s_g=unique(comb_list[,2])
  lambda_l_g=unique(comb_list[,3])
  gamma_ada=unique(comb_list[,4])
  
  
  new_comb_list3<-list()
  max_vec_l<-numeric()
  for (kkk in 1:length(gamma_ada)) {
    
    kk=1
    max_vec_nc<-sd_vec_nc<-numeric()
    new_comb_list<-matrix(0,length(lambda_s_g)*length(lambda_l_g),4)
    
    for (jj in 1:length(lambda_l_g)) {
      for (ii in 1:length(lambda_s_g)) {
        indexes<-which(comb_list[,2]==lambda_s_g[ii]&comb_list[,3]==lambda_l_g[jj]&comb_list[,4]==gamma_ada[kkk])
        par_index<-par[indexes]
        sd_index<-sds[indexes]
        max<-which.max(par_index)
        onese<-which(par_index[1:(max)]>=par_index[max]-k1*sd_index[max])[1]
        
        max_vec_nc[kk]<-par_index[onese]
        sd_vec_nc[kk]<-sd_index[onese]
        
        new_comb_list[kk,]<-as.numeric(comb_list[indexes[onese],])
        kk=kk+1
        
      }
    }
    kk=1
    max_vec_s<-sd_vec_s<-numeric()
    new_comb_list2<-matrix(0,length(lambda_l_g),4)
    for (ii in 1:length(lambda_l_g)) {
      indexes<-which(new_comb_list[,3]==lambda_l_g[ii]&new_comb_list[,4]==gamma_ada[kkk])
      par_index<-max_vec_nc[indexes]
      sd_index<-sd_vec_nc[indexes]
      max<-which.max(par_index)
      
      onese<-max(which(par_index>=par_index[max]-k2*sd_index[max]))
      max_vec_s[kk]<-par_index[onese]
      sd_vec_s[kk]<-sd_index[onese]
      new_comb_list2[kk,]<-as.numeric(new_comb_list[indexes[onese],])
      kk=kk+1
    }
    
    
    
    par_index<-max_vec_s
    sd_index<-sd_vec_s
    max<-which.max(par_index)
    if(k3*sd_index[max]>0.5*abs(max(par_index)-min(par_index)))lim=0.5*abs(max(par_index)-min(par_index)) else lim=k3*sd_index[max]
    onese<-max(which(par_index>=par_index[max]-lim))
    max_vec_l[kkk]<-par_index[onese]
    sd_vec_l<-sd_index[onese]
    new_comb_list3[[kkk]]<-as.numeric(new_comb_list2[onese,])
    
  }
  ind_max<-which.max(max_vec_l)
  
  num_clusters_opt<-new_comb_list3[[ind_max]][1]
  lambda_s_opt<-new_comb_list3[[ind_max]][2]
  lambda_l_opt<-new_comb_list3[[ind_max]][3]
  gamma_opt<-new_comb_list3[[ind_max]][4]
  
  return(c(num_clusters_opt=num_clusters_opt,
           lambda_s_opt=lambda_s_opt,
           lambda_l_opt=lambda_l_opt,
           gamma_opt=gamma_opt))
}

spfclust_fk<-function(mod_cv_grid,K_fixed,...){
  
  if(!K_fixed%in%unique(mod_cv_grid$comb_list[,1])) stop("K_fixed not evaluated",call. = FALSE)
  ind<-which(mod_cv_grid$comb_list[,1]==K_fixed)
  comb_list=mod_cv_grid$comb_list[ind,]
  par=mod_cv_grid$CV[ind]
  sds=mod_cv_grid$CV_sd[ind]
  zeros=mod_cv_grid$zeros[ind]
  k1=mod_cv_grid$ks[1]
  k2=mod_cv_grid$ks[2]
  k3=mod_cv_grid$ks[3]
  ksdrule<-get_ksdrule(par,sds,comb_list,k1,k2,k3)
  num_clusters_opt<-ksdrule[1]
  lambda_s_opt<-ksdrule[2]
  lambda_l_opt<-ksdrule[3]
  gamma_opt<-ksdrule[4]
  # 
  # mod_opt<-spfclust(data=data,trace=T,lambda_l = lambda_l_opt,lambda_s =lambda_s_opt,q=40,plot = TRUE,K=num_clusters_opt,...)
  # 
  out<-list(mod_opt=NA,
            num_clusters_opt=num_clusters_opt,
            lambda_l_opt=lambda_l_opt,
            lambda_s_opt=lambda_s_opt,
            gamma_opt=gamma_opt,
            comb_list=comb_list,
            CV=par,
            CV_sd=sds,
            zeros=zeros,ks=c(k1,k2,k3))
  
  return(out)
}
