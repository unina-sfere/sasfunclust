
#' @export
sasfclust <-
  function(x=NULL,curve=NULL,timeindex=NULL,data=NULL, q = 30,
            lambda_s = 1e1, lambda_l = 1e1, G = 2, tol = 10^-7, maxit = 50,
           grid = seq(0, 1, length = 100),par_LQA=list(eps_diff = 1e-06,MAX_iter_LQA=200,eps_LQA = 1e-05),
            plot= F,trace=F,init="kmeans",varcon="diagonal",lambda_s_ini=NA)
  {
    der=2
    gamma_ada=1
    CK=0
    hard = F
    perc_rankpapp=pert=NA
    if(G==1)lambda_l=0
    # Initialize the parameters
    initfit <- sasfclustinit(data = data, pert = pert, grid = grid,  q = q, G = G,der=der,gamma_ada=gamma_ada,lambda_s_ini=lambda_s_ini,init=init,varcon=varcon)

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
      cluster.mean <- matrix(0,G,dim(FullS)[1])
      for(k in 1:G)
        cluster.mean[k,] <- FullS %*% (t(parameters$mu)[,k])
      plot(grid,grid,ylim=range(cluster.mean),type='n',ylab="Cluster Mean")
      for(k in 1:G)
        graphics::lines(grid,cluster.mean[k,], col = 4, lwd = 2)
    }
    # Main loop. Iterates between M and E steps and stops when the stopping condition is met.
    lk_old=0
    lk_new=loglik (parameters, data, vars, S,W,AW_vec,P_tot,lambda_s,lambda_l)
    if (trace)    print(paste("Iteration", 0,": Sigma = ",sigma.new," loglk = ",lk_new[1]," ploglk = ",lk_new[2]))
    lk_new=lk_new[2]
    while(abs(lk_old - lk_new) > tol & (ind <= maxit)) {
      parameters <- sasfclustMstep(parameters, data, vars, S, tol, hard,lambda_s,lambda_l,W,AW_vec,P_tot,par_LQA,CK,perc_rankpapp,varcon=varcon)
      vars <- sasfclustEstep(parameters, data, vars, S, hard)
      lk_old<-lk_new
      lk_i<-loglik (parameters, data, vars, S,W,AW_vec,P_tot,lambda_s,lambda_l)
      lk_new<--lk_i[2]
      sigma.old <- sigma.new
      sigma.new <- parameters$sigma[1]
      # print(ind)
      if (trace)
        print(paste("Iteration", ind,": Sigma = ",sigma.new," loglk = ",lk_i[1]," ploglk = ",lk_i[2]))
      #Plot cluster mean curves.
      if(plot){
        cluster.mean <- matrix(0,G,dim(FullS)[1])
        for(k in 1:G)
          cluster.mean[k,] <- FullS %*% (t(parameters$mu)[,k])
        plot(grid,grid,ylim=range(cluster.mean),type='n',ylab="Cluster Mean")
        for(k in 1:G)
          graphics::lines(grid,cluster.mean[k,], col = 4, lwd = 2)
      }
      ind <- ind + 1

    }
    mod=list(data=data,parameters = parameters, vars = vars, FullS = FullS,grid=grid,S=S,
         W=W,AW_vec=AW_vec,P_tot=P_tot,P=P,lambda_s=lambda_s,lambda_l=lambda_l,CK=CK,basis=basis)
    mean_fd<-fda::fd(t(parameters$mu),basis)
    clus<-classify(mod)

    out<-list(mod=mod,
              mean_fd=mean_fd,
              clus=clus)

    return(out)
    }

sasfclustinit <-
  function(data, pert = 0, grid = seq(0.01, 1, length = 100),  q = 5,G = G,der=der,gamma_ada=gamma_ada,lambda_s_ini=lambda_s_ini,init=init,varcon=varcon){
    S <- FullS <- NULL
    # This function initializes all the parameters.
    # Produce spline basis matrix
    basis<-fda::create.bspline.basis(c(grid[1],grid[length(grid)]),nbasis = q)
    FullS<-fda::eval.basis(grid,basis)
    W<-fda::eval.penalty(basis,2)
    S <- FullS[data$timeindex,  ]
    order<-q-length(basis$params)
    breaks<-basis$params

    ext_break<-c(rep(grid[1],order),breaks,rep(grid[length(grid)],order))
    weights_vec<-rep(diff(ext_break,lag=order)/order,each=G)

    # Get pairwise matrix
    if(G!=1)
    {
      P<-matrix(0,((G-1)^2+(G-1))/2,G)
      ind<-c(1,G-1)
      for (ii in 1:(G-1)) {

        P[ind[1]:ind[2],ii]<-rep(1,length(ind[1]:ind[2]))
        if(length(ind[1]:ind[2])==1)
          aa<--1
        else
          aa<-diag(rep(-1,length(ind[1]:ind[2])))

        P[ind[1]:ind[2],(ii+1):G]<-aa
        ind<-ind+c((G-1)-(ii-1),(G-1)-(ii))
        ind<-c(min(((G-1)^2+(G-1))/2,ind[1]),min(((G-1)^2+(G-1))/2,ind[2]))
      }

      P_tot<-matrix(0,q*((G-1)^2+(G-1))/2,G*q)
      for (ii in 1:G) {
        P_tot[, ((ii - 1) * q + 1):((ii) * q)]<-kronecker(diag(q), P[, ii])
      }
      P_tot<-Matrix::Matrix(P_tot,sparse = TRUE)
    }
    else{
      P<-P_tot<-diag(q)
    }

    # Weight approximation L1 penalty
    order<-q-length(basis$params)
    breaks<-basis$params

    ext_break<-c(rep(grid[1],order),breaks,rep(grid[length(grid)],order))
    weights_vec<-rep(diff(ext_break,lag=order)/order,each=((G-1)^2+(G-1))/2)
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
        loglam         = seq(-10, 6, 0.25)
        Gcvsave        = numeric()
        for(i in 1:length(loglam)){
          fdPari  = fda::fdPar(basis_start, Lfdobj=2, 10^loglam[i])
          Sm.i    = fda::smooth.basis(grid_i, X, fdPari)
          Gcvsave[i] = sum(Sm.i$gcv)

        }
        lambda_s=if(is.na(lambda_s_ini))10^loglam[which.min(Gcvsave)]else lambda_s_ini
        # print(lambda_s)
        fdPari  = fda::fdPar(basis_start, Lfdobj=2,lambda_s)
        points<-t(fda::smooth.basis(grid_i, X, fdPari)$fd$coefs)
        # print(lambda_s)
        # plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
        #      ylab=expression(GCV(lambda)), lwd=2 )
        # plot(smooth.basis(grid_i, X, fdPari)$fd)
      }
      else{
        basis_start<-basis
        loglam         = seq(-10, 6, 0.25)
        Gcvsave        = matrix(0,N,length(loglam))
        points <- matrix(0,N,sum(q))
        for (i in 1:N){
          print(i)
          xi <- data$x[data$curve==i]
          grid_i <- grid[data$timeindex[data$curve==i]]

          for(l in 1:length(loglam)){
            fdPari  = fda::fdPar(basis_start, Lfdobj=2, 10^loglam[l])
            Sm.i    = fda::smooth.basis(grid_i, xi, fdPari)
            Gcvsave[i,l] = sum(Sm.i$gcv)
          }
          lambda_s[i]=10^loglam[which.min(Gcvsave[i,])]
          fdPari  = fda::fdPar(basis_start, Lfdobj=2,lambda_s[i])
          points[i,]<-t(fda::smooth.basis(grid_i, xi, fdPari)$fd$coefs)

          # plot(loglam,Gcvsave[i,], 'o', las=1, xlab=expression(log[10](lambda)),
          #      ylab=expression(GCV(lambda)), lwd=2 )
          # plot(smooth.basis(grid_i, xi, fdPari)$fd)
        }
      }
    }

    # Use k-means to get initial cluster memberships from points.
    if(G > 1){
      if(init=="kmeans")class <-  stats::kmeans(points, G, 1000,nstart = 10000)$cluster
      else if(init=="model-based") class <- mclust::Mclust(points,G,verbose=FALSE,modelNames = "EII")$classification
      else if(init=="hierarchical")class <-  factoextra::hcut(points,k =G )$cluster
    }
    else {
      class <- rep(1, N)
    }
    # matplot(t(points),type="l",col=class)
    # Initial estimates for the posterior probs of cluster membership.
    piigivej <- matrix(0, N, G)
    piigivej[col(piigivej) == class] <- 1
    pi=apply(piigivej,2,mean)
    #
    # # Calculate coefficeints for cluster means.
    classmean <- matrix(0,G,q)
    for (k in 1:G)
      classmean[k,] <- apply(as.matrix(points[class==k,]),2,mean)


    #  X<-t(sapply(1:N,function(i) data$x[data$curve==i]))
    #  mu_eval<-matrix(0,dim(X)[2],G)
    # new_X<-sapply(1:G,function(i)if(is.null(dim(X[class==i,]))) X[class==i,] else apply(X[class==i,],2,mean))
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
    gamma<-array(0, c(N, G, sum(q)))
    gprod <- NULL
    if(G==1){
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
    #  mu_eval<-matrix(0,dim(X)[2],G)
    # new_X<-sapply(1:G,function(i)if(is.null(dim(X[class==i,]))) X[class==i,] else apply(X[class==i,],2,mean))
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
    if(G!=1){
      AW<-matrix(0,((G-1)^2+(G-1))/2,q)
      for (ii in 1:q) {

        AW[,ii]<-1/(abs(P%*%mu_start[,ii])^gamma_ada)
      }
      AW_vec<-matrixcalc::vec(AW)*(weights_vec)
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

sasfclustMstep <-
  function(parameters, data, vars, S, tol,  hard,lambda_s,lambda_l,W,AW_vec,P_tot,par_LQA,CK,perc_rankpapp,varcon=varcon)
  {
    # This function implements the M step of the EM algorithm.
    G <- dim(parameters$mu)[1]
    mu<-parameters$mu
    gamma <- vars$gamma
    gcov <- vars$gcov
    curve <- data$curve
    piigivej <- vars$piigivej
    N <- dim(gamma)[1]
    G <- dim(mu)[1]
    n <- length(curve)
    q <- dim(S)[2]
    # Compute pi.
    if(hard)
      parameters$pi <- rep(1/G, G)
    else parameters$pi <- (apply(vars$piigivej, 2, mean)*N+CK)/(N+G*CK)
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

    # This loop iteratively LDA to get mu
    # when they have converged.
    W_star<-matrix(0,q*G,q*G)
    for (i in 1:G) {
      W_star[((i-1)*q+1):((i)*q),((i-1)*q+1):((i)*q)]<-W
    }
    W_star<-Matrix::Matrix(W_star, sparse = TRUE)
    x <- data$x
    n_i<-sapply(1:N,function(ii)length(which(curve==ii)))
    numden<-get_numden( data$x, data$curve, data$time,  as.matrix(S),  piigivej,  gcov,n_i,gamma)
    VY<-Matrix::Matrix(numden[[1]],sparse = TRUE)
    S.den<-Matrix::Matrix(numden[[2]],sparse=TRUE)

    z_int=1
    diff_inter <- 100
    mu_old=mu

    while(diff_inter>par_LQA$eps_LQA && z_int<=par_LQA$MAX_iter_LQA) {

      mu_vec<-matrixcalc::vec(t(mu))
      diff_mat<-abs(P_tot%*%mu_vec)
      diff_mat[diff_mat < par_LQA$eps_diff]<-par_LQA$eps_diff
      V_l<-Matrix::Matrix(diag(as.numeric(AW_vec/(2*diff_mat))),sparse = TRUE)
      mu_vec<-solve(S.den*(1/parameters$sigma)+lambda_s*2*W_star+2*lambda_l*Matrix::t(P_tot)%*%V_l%*%P_tot)%*%VY*(1/parameters$sigma)
      mu<-matrix(mu_vec,G,q,byrow = T)
      diff_inter<-sum(abs(mu-mu_old))/(sum(abs(mu_old)))
      mu_old<-mu
      z_int=z_int+1
    }


    sigma<-get_sigma( data$x, data$curve, data$time,  as.matrix(S),  piigivej,  gcov,n_i,gamma,mu)
    sigma=as.numeric(sigma)/n
    parameters$mu <- mu
    parameters$sigma <- sigma

    parameters
  }

sasfclustEstep <-
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



#' @export
sasfclust_cv<-function(data,num_cluster=6,lambda_l_g=10^seq(-1,2),lambda_s_g=10^seq(-5,-3),K_fold=5,test=NA,grid=NA,grid_2=NA,k1=1,k2=0,k3=1,maxit=50,q=30,
                     init="kmeans",lambda_s_ini=NA,varcon="diagonal",...){

  N<-length(unique(data$curve))
  comb_list<-expand.grid(num_cluster,lambda_s_g,lambda_l_g)
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

     ran_seq<-sample(seq(1, N), N, replace=FALSE)
     split_vec<-split(ran_seq,cut(seq(1,N),breaks=K_fold,labels=FALSE))
     l_i<-zeros_vec<-numeric()
     data_fold<-data_i<-list()
     for(lll in 1:K_fold){

       ind_fold<-as.numeric(unlist(split_vec[-lll]))
       ind_i<-split_vec[[lll]]

       data_fold$x<-matrixcalc::vec(sapply(1:length(ind_fold),function(ll)data$x[data$curve==ind_fold[ll]]))
       data_fold$timeindex<-unlist(lapply(1:length(ind_fold),function(ll)data$timeindex[data$curve==ind_fold[ll]]))
       data_fold$curve<-unlist(lapply(1:length(ind_fold),function(ll)rep(ll,length(data$timeindex[data$curve==ind_fold[ll]]))))
       data_i$x<-matrixcalc::vec(sapply(1:length(ind_i),function(ll)data$x[data$curve==ind_i[ll]]))
       data_i$timeindex<-unlist(lapply(1:length(ind_i),function(ll)data$timeindex[data$curve==ind_i[ll]]))
       data_i$curve<-unlist(lapply(1:length(ind_i),function(ll)rep(ll,length(data$timeindex[data$curve==ind_i[ll]]))))

       mod<-sasfclust(data=data_fold,lambda_l = lambda_l_i,lambda_s =lambda_s_i,G=num_clusters_i,grid=grid,maxit=maxit,q=q,
                      init=init,lambda_s_ini=lambda_s_ini,varcon=varcon,...)
       l_i[lll]<-loglik(parameters = mod[[1]]$parameters,data = data_i,vars = mod[[1]]$vars, FullS = mod[[1]]$FullS,W = mod[[1]]$W,AW_vec = mod[[1]]$AW_vec,P_tot = mod[[1]]$P_tot,
                        lambda_s = mod[[1]]$lambda_s)[1]#/length(ind_i)
       # print(l_i[lll])
       zeros_vec[lll]<-get_zero(mod[[1]])
       rm(mod)
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

      mod<-sasfclust(data=data,lambda_l = lambda_l_i,lambda_s =lambda_s_i,maxit=maxit,q=q,G=num_clusters_i,grid=grid,...)
      l_i<-loglik(mod[[1]]$parameters,test,mod[[1]]$vars, mod[[1]]$S,mod[[1]]$W,mod[[1]]$AW_vec,mod[[1]]$P_tot,mod[[1]]$lambda_s,mod[[1]]$lambda_l)[1]
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


  cores <- parallel::detectCores()
  #cat(cores)
  vec_par<-parallel::mclapply(seq(1,length(comb_list[,1])),parr_fun,mc.cores = cores)

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
  # mod_opt<-spfclust(data=data,trace=T,lambda_l = lambda_l_opt,lambda_s =lambda_s_opt,q=40,plot = TRUE,G=num_clusters_opt,...)
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



#' @export
simulate_data<-function(scenario,n_i=50,nbasis=30,length_tot=50,sd=1,sd2_basis=1) {


  grid<-seq(0,1,length.out = length_tot)
  domain<-c(0,1)


  X_basis<-fda::create.bspline.basis(domain,norder = 4,nbasis = nbasis)
  mean_list<-list()
  if(scenario=="Scenario I"){
    mean_list[[1]]<-c(rep(1.5,nbasis/6),rep(0,nbasis*5/6))
    mean_list[[2]]<-c(rep(-1.5,nbasis/6),rep(0,nbasis*5/6))
    # mean_list[[4]]<-c(rep(-2.5,nbasis/3),rep(-1.5,nbasis/3),rep(0,nbasis/3))

  }
  if(scenario=="Scenario II"){
    mean_list[[1]]<-c(rep(3,nbasis/6),rep(1.5,nbasis/6),rep(0,nbasis/6),rep(0,nbasis/2))
    mean_list[[2]]<-c(rep(0,nbasis/6),rep(1.5,nbasis/6),rep(0,nbasis/6),rep(0,nbasis/2))
    mean_list[[3]]<-c(rep(0,nbasis/6),rep(-1.5,nbasis/6),rep(0,nbasis/6),rep(0,nbasis/2))
    # mean_list[[4]]<-c(rep(-2.5,nbasis/3),rep(-1.5,nbasis/3),rep(0,nbasis/3))

  }
  if(scenario=="Scenario III"){
    mean_list[[1]]<-c(rep(1.5,nbasis/6),rep(3,nbasis/6),rep(1.5,nbasis/6),rep(0,nbasis/2))
    mean_list[[2]]<-c(rep(1.5,nbasis/6),rep(0,nbasis/6),rep(1.5,nbasis/6),rep(0,nbasis/2))
    mean_list[[3]]<-c(rep(-1.5,nbasis/6),rep(0,nbasis/6),rep(-1.5,nbasis/6),rep(0,nbasis/2))
    mean_list[[4]]<-c(rep(-1.5,nbasis/6),rep(-3,nbasis/6),rep(-1.5,nbasis/6),rep(0,nbasis/2))
    # mean_list[[4]]<-c(rep(-2.5,nbasis/3),rep(-1.5,nbasis/3),rep(0,nbasis/3))

  }

  if(scenario=="Scenario IV"){
    mean_list[[1]]<-c(rep(1.5,nbasis/3),rep(3,nbasis/3),rep(1.5,nbasis/3))
    mean_list[[2]]<-c(rep(1.5,nbasis/3),rep(0,nbasis/3),rep(1.5,nbasis/3))
    mean_list[[3]]<-c(rep(-1.5,nbasis/3),rep(0,nbasis/3),rep(-1.5,nbasis/3))
    mean_list[[4]]<-c(rep(-1.5,nbasis/3),rep(-3,nbasis/3),rep(-1.5,nbasis/3))
    # mean_list[[4]]<-c(rep(-2.5,nbasis/3),rep(-1.5,nbasis/3),rep(0,nbasis/3))

  }
  mu_fd<-fda::fd(t(do.call("rbind",mean_list)),X_basis)
  if(length(mean_list)==1)X_coef<-t(MASS::mvrnorm(n_i, mean_list[[1]],diag(nbasis)*sd2_basis))
  else{
    X_coef<-t(MASS::mvrnorm(n_i, mean_list[[1]],diag(nbasis)*sd2_basis))
    for (ii in 2:length(mean_list)) {
      X_coef<-cbind(X_coef,t(MASS::mvrnorm(n_i, mean_list[[ii]],diag(nbasis)*sd2_basis)))

    }
  }

  X_fd<-fda::fd(X_coef,X_basis)
  X<-fda::eval.fd(grid,X_fd)
  X<-X+matrix(stats::rnorm(dim(X)[1]*dim(X)[2],0,sd),dim(X)[1],dim(X)[2])
  n_obs<-dim(X)[2]

  vec<-list(x=matrixcalc::vec(X),timeindex=rep(1:length(grid),n_obs),curve=rep(1:n_obs,each=length(grid)))

  out<-list(X=X,
            X_fd=X_fd,
            mu_fd=mu_fd,
            grid=grid,
            vec=vec)

  return(out)
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
  # mod_opt<-spfclust(data=data,trace=T,lambda_l = lambda_l_opt,lambda_s =lambda_s_opt,q=40,plot = TRUE,G=num_clusters_opt,...)
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
