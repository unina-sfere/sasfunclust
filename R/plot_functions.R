#' @export
plot<-function(mod,m1=NULL,m2=NULL,m3=NULL,lambda_s_min=NULL,...){

  if(mod$class=="sasfclust_cv"){
    comb_list_i<-mod$comb_list
    CV_i<-mod$CV
    sd_i<-mod$CV_sd
    zeros_i<-mod$zeros

    if(!is.null(lambda_s_min)){
      comb_list_i<-mod$comb_list[which(mod$comb_list[,2]>=lambda_s_min),1:3]
      CV_i<-mod$CV[which(mod$comb_list[,2]>=lambda_s_min)]
      sd_i<-mod$CV_sd[which(mod$comb_list[,2]>=lambda_s_min)]
      zeros_i<-mod$zeros[which(mod$comb_list[,2]>=lambda_s_min)]
    }

    x<-seq(1,length(comb_list_i[,1]))
    labels<-lapply(1:length(comb_list_i[,1]),function(ii){a<-as.character(signif(comb_list_i[ii,],digits = 1));    paste(a[1],a[2],a[3])})
    layout(matrix(rbind(c(1,1),c(2,3)),2,2))
    base::plot(CV_i,pch=16,cex=0.5,col=2,type="l",xaxt="n",xlab="",ylab="CV",ylim=c(min(CV_i)-max(sd_i),max(CV_i)+max(sd_i)),...)
    points(CV_i,pch=16,cex=0.5,col=2)
    segments(x-0.1,CV_i+sd_i,x+0.1)
    segments(x-0.1,CV_i-sd_i,x+0.1)
    text(x=x, y=par()$usr[3]-0.00001*(par()$usr[4]-par()$usr[3]),
         labels=labels, srt=90, adj=1, xpd=TRUE)
    text(x=x, y=par()$usr[4]+0.04*(par()$usr[4]-par()$usr[3]),
         labels=as.character(round(zeros_i*100)), srt=90, adj=1, xpd=TRUE)
    abline(v=which(CV_i==max(CV_i)))
    abline(h=max(CV_i))
    lamb_s<-unique(comb_list_i[,2])
    lamb_L<-unique(comb_list_i[,3])
    num_cluster<-unique(comb_list_i[,1])
    par<-CV_i
    sds<-sd_i
    zeros<-zeros_i
    comb_list<-comb_list_i
    if(is.null(m1)) m1<-mod$ms[1]
    if(is.null(m2)) m2<-mod$ms[2]
    if(is.null(m3)) m3<-mod$ms[3]

    kk=1
    max_vec_nc<-sd_vec_nc<-zero_vec<-numeric()
    new_comb_list<-matrix(0,length(lamb_s)*length(lamb_L),3)

    for (jj in 1:length(lamb_L)) {
      for (ii in 1:length(lamb_s)) {
        indexes<-which(comb_list[,2]==lamb_s[ii]&comb_list[,3]==lamb_L[jj])
        par_index<-par[indexes]
        sd_index<-sds[indexes]
        zero_index<-zeros[indexes]
        max<-which.max(par_index)
        if(m1*sd_index[max]>0.5*abs(max(par_index)-min(par_index)))lim=0.5*abs(max(par_index)-min(par_index)) else lim=m1*sd_index[max]
        onese<-which(par_index[1:(max)]>=par_index[max]-lim)[1]
        max_vec_nc[kk]<-par_index[onese]
        sd_vec_nc[kk]<-sd_index[onese]
        zero_vec[kk]<-zero_index[onese]
        new_comb_list[kk,]<-as.numeric(comb_list[indexes[onese],])
        kk=kk+1

      }
    }

    labels<-lapply(1:length(new_comb_list[,1]),function(ii){a<-as.character(signif(new_comb_list[ii,],digits = 1));    paste(a[1],a[2],a[3])})
    base::plot(max_vec_nc,pch=16,cex=0.5,col=2,type="l",xaxt="n",xlab="",ylab="CV fixed G",ylim=c(min(max_vec_nc)-max(sd_vec_nc),max(max_vec_nc)+max(sd_vec_nc)))
    points(max_vec_nc,pch=16,cex=0.5,col=2)
    segments(x-0.1,max_vec_nc+sd_vec_nc,x+0.1)
    segments(x-0.1,max_vec_nc-sd_vec_nc,x+0.1)
    text(x=x, y=par()$usr[3]-0.00001*(par()$usr[4]-par()$usr[3]),
         labels=labels, srt=90, adj=1, xpd=TRUE)
    text(x=x, y=par()$usr[4]+0.1*(par()$usr[4]-par()$usr[3]),
         labels=as.character(round(zero_vec*100)), srt=90, adj=1, xpd=TRUE)

    kk=1
    max_vec_s<-sd_vec_s<-zero_vec2<-numeric()
    new_comb_list2<-matrix(0,length(lamb_L),3)
    for (ii in 1:length(lamb_L)) {
      indexes<-which(new_comb_list[,3]==lamb_L[ii])
      par_index<-max_vec_nc[indexes]
      sd_index<-sd_vec_nc[indexes]
      zero_index<-zero_vec[indexes]
      max<-which.max(par_index)
      onese<-max(which(par_index>=par_index[max]-m2*sd_index[max]))
      max_vec_s[kk]<-par_index[onese]
      sd_vec_s[kk]<-sd_index[onese]
      zero_vec2[kk]<-zero_index[onese]
      new_comb_list2[kk,]<-as.numeric(new_comb_list[indexes[onese],])
      kk=kk+1
    }

    labels_L<-lapply(1:length(new_comb_list2[,1]),function(ii){a<-as.character(signif(new_comb_list2[ii,],digits = 1));    paste(a[1],a[2],a[3])})
    base::plot(max_vec_s,pch=16,cex=0.5,col=2,type="l",xaxt="n",xlab="",ylab="CV fixed G and lambda_s",ylim=c(min(max_vec_s)-max(sd_vec_s),max(max_vec_s)+max(sd_vec_s)))
    points(max_vec_s,pch=16,cex=0.5,col=2)
    segments(x-0.1,max_vec_s+sd_vec_s,x+0.1)
    segments(x-0.1,max_vec_s-sd_vec_s,x+0.1)
    text(x=x, y=par()$usr[3]-0.00001*(par()$usr[4]-par()$usr[3]),
         labels=labels_L, srt=90, adj=1, xpd=TRUE)
    text(x=x, y=par()$usr[4]+0.1*(par()$usr[4]-par()$usr[3]),
         labels=as.character(round(zero_vec2*100)), srt=90, adj=1, xpd=TRUE)

  }
  else if(mod$class=="sasfclust"){
    G<-dim(mod$mean_fd$coefs)[2]
    range<-mod$mean_fd$basis$rangeval
    grid_eval<-seq(range[1],range[2],length.out = 500)
    eval_mu<- fda::eval.fd(grid_eval,mod$mean_fd)
    par(mfrow=c(1,2))

    graphics::matplot(grid_eval,eval_mu,ylab = "",xlab="",lty=1:G,type="l",xlim=range,ylim=c(min(mod$mod$data$x),max(mod$mod$data$x)))
    abline(h=0,col=adjustcolor("grey", alpha = 0.6))
    graphics::title("Cluster means")
    legend("topright",legend = paste0("Cluster ",1:G),lty=1:G,col=1:G)

    base::plot(0,type="n",xlim=range,ylim=c(min(mod$mod$data$x),max(mod$mod$data$x)),xlab="",ylab="")
    graphics::title("Classified observations")
    for(ii in 1:length(unique(mod$mod$data$curve))){
      lines(mod$mod$grid[mod$mod$data$timeindex[which(mod$mod$data$curve==ii)]],mod$mod$data$x[which(mod$mod$data$curve==ii)],col=mod$clus[[1]][ii],lty=mod$clus[[1]][ii])
    }
    legend("topright",legend = paste0("Cluster ",1:G),lty=1:G,col=1:G)


  }





}

