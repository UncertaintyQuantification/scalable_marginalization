library(deSolve) 
library(spatPomp)
library(FastGaSP)
library(RobustGaSP)
library(nloptr)
library(Matrix)
library(Rcpp)
library(RcppEigen)

sourceCpp(file='src/functions_Feb26_2022.cpp')

negative_LJ_truncated_kernel <- function(r){
  # r is the vector of Euclidean distances
  sigma = 1
  p =8
  q =2
  epsilon = 1
  r_trunc = 0.95*sigma
  f = rep(0, length(r))
  
  index_trunc = which(r<r_trunc)
  
  if(length(index_trunc)>0){
    rinv = sigma/r[-index_trunc]
    rinv_trunc=sigma/r_trunc
    f[-index_trunc] = p*q*epsilon*(rinv^(q+2) -rinv^(p+2))/((p-q)*sigma^2)
    
    f1=p*q*epsilon*(rinv_trunc^(q+2) -rinv_trunc^(p+2))/((p-q)*sigma^2)
    f2=p*q*epsilon/((p-q)*sigma^3)*((rinv_trunc)^(p+3)*(p+2)-(rinv_trunc)^(q+3)*(q+2))
    a=f2/(-12*f1*r_trunc^11)
    b=f1*exp(a*r_trunc^12)
    f[index_trunc]=b*exp(-a*r[index_trunc]^12)
  }
  else{
    rinv = sigma/r
    f = p*q*epsilon*(rinv^(q+2) -rinv^(p+2))/((p-q)*sigma^2)
  }
  
  return(f)
}



negative_LJ_kernel<- function(r){
  # r is the vector of Euclidean distances
  sigma = 1
  p =8
  q =2
  epsilon = 1
  f = rep(0, length(r))
  
  rinv = sigma/r
  f = p*q*epsilon*(rinv^(q+2) -rinv^(p+2))/((p-q)*sigma^2)
  
  return(f)
}

norm_vec <- function(v){
  sqrt(sum(v^2))
}


###modified version to another order
RHSfn <- function(t,x,N){
  # The RHS function of the ODE
  # \dot X (t) = RHSfn(t,X)
  # Input
  #   x   - the state vector DN x 1;
  #        [x_1(1), ..., x_N(1), ... , x_1(D),..., x_N(D)]^T
  #   t   - time
  #
  # Output
  #   dx  -  right hand side of the ode  DN x 1
  #          [x_1'(1),..., x_N'(1), ... , x_1'(D),..., x_N'(D)]^T
  #
  #x <- as.vector(x)
  D = length(x)/N
  y = matrix(x,N,D)  #[x_1...x_N]  d x N
  aa= matrix(0,N,D)
  ##change this to current order?
  for(i in 1:N){
    temp  = (t(y)- matrix(y[i,],D,N))    # d x N [x_1-x_i,..., x_N-x_i]
    
    DD  = apply(temp,2,norm_vec)      # 1 x N
    DD[DD==0] = 1                     #this step is to avoid 0 x inf= NAN
    aa[i,] = (temp)%*%negative_LJ_truncated_kernel(DD)    # d x 1 not normalized
    
  }
  dx = as.vector(aa)  # dN x 1
  return( list(dx) )
  
}



Generate_training_data <- function(sysInfo,obsInfo,solverInfo){
  # sysInfo = c(N, d)
  # obsInfo <- list(time_vec, M)
  # 
  
  ## basic setting of the system
  N = sysInfo[1]         # number of agents
  D = sysInfo[2]        # dim of state vectors
  # myODE = RHSfn(t,x,N,sysInfo.phi{1});
  time_vec = obsInfo[[1]]
  M = obsInfo[[2]]
  #use_derivative = obsInfo[[3]]
  
  xpath_train = array(0, c(D*N, length(time_vec)+1, M))
  dxpath_train = array(0, c(D*N, length(time_vec)+1, M))
  
  # test_time_vec= [time_vec(1:end-1) linspace(obsInfo.time_vec(end),sysInfo. T_f,200)];                                                                     % final time the system will reach steady state
  # xpath_test = zeros(d*N,length(test_time_vec),obsInfo.M);
  # dxpath_test = zeros(d*N,length(test_time_vec),obsInfo.M);
  
  for(i in 1:M){    # trajectories with random initial conditions for learning interaction kernel
    if(obsInfo[[3]]=='normal'){
      x0 = rnorm(D*N)     # random initial condition
    }else if(obsInfo[[3]]=='uniform'){
      x0 = 5*runif(D*N)     # random initial condition
      
    }else if(obsInfo[[3]]=='log_uniform'){
      x0 = log(10^{-3})+(log(5)-log(10^{-3}) )*runif(D*N)     # random initial condition
      x0=exp(x0)
    }else if(obsInfo[[3]]=='small_distance'){
      x0 = log(10^{-3})+(log(10^{-1})-log(10^{-3}) )*runif(D*N)     # random initial condition
      x0=exp(x0)
    }
    output = ode(x0, times=c(0, time_vec), func=RHSfn, parms=N, method =solverInfo) ##let me just try euler first
    
    output = t(as.matrix(output[,-1])) ##keep initial value
    
    xpath_train[,,i] = output
    
    for(j in 1: (length(time_vec)+1)){
      dxpath_train[,j,i] = (RHSfn(t=time_vec[j], x=xpath_train[,j,i], N=N))[[1]]
    }
  }
  res = as.list(1:2)
  res[[1]] = xpath_train
  res[[2]] = dxpath_train
  
  return(res)
}




S=30
n_repeat=5
record_time_CG=matrix(NA,S,n_repeat)
record_time_full=matrix(NA,S,n_repeat)

record_rmse_CG_full=matrix(NA,S,n_repeat)
record_rmse_full=matrix(NA,S,n_repeat)
record_rmse_CG=matrix(NA,S,n_repeat)


for(i_S in 1:S){
  for(j_repeat in 1:n_repeat){
    set.seed(i_S+S*j_repeat)
    
    N=i_S*20
    D=2
    L_sim=1      # time steps
    L=0      ##how many data to be use, initial as 0
    #h = 0.05
    h = 5*10^{-4}
    
    M = 1    # the number of trajectories
    
    ###test euler first, work on other approach later
    sysInfo = c(N, D) 
    solverInfo='euler' ##just for demonstration purposes
    obsInfo = as.list(1:3)
    obsInfo[[1]] = seq(from=h, by=h, length.out=L_sim)
    obsInfo[[2]] = M
    obsInfo[[3]] = 'log_uniform'
    
    res = Generate_training_data(sysInfo,obsInfo,solverInfo)
    input_pos = res[[1]][,,1] ##contains zero
    
    v_all = (res[[2]][,,1])[,1:(L+1)]
    
    d_all = matrix(NA,N*(L+1),N)
    A_all = matrix(NA,D*N*(L+1),N)
    for(t in 1:(L+1)){
      A=matrix(NA,D*N,N)
      for(i in 1:N){
        A[1:N,i]=-(input_pos[1:N,t]-input_pos[i,t])#input_pos[i,t]-input_pos[1:N,t]
        A[(N+1):(2*N),i]=-(input_pos[N+(1:N),t]-input_pos[N+i,t])#input_pos[N+i,t]-input_pos[N+(1:N),t]
      }
      A_all[((t-1)*N*D+1):(t*N*D),]=A
      
      
      d=matrix(NA,N,N) ##ith column is the distance of particle i
      for(i in 1:N){
        d[,i]=((input_pos[i,t]-input_pos[1:N,t])^2+(input_pos[i+N,t]-input_pos[N+(1:N),t])^2)^(1/2)
      }
      
      d_all[((t-1)*N+1):(t*N),] = d
    }
    
    

    
    ###done with simulation
    beta=.2
    nu  = 10^{-5}
    param = log(c(beta, nu))
    
    
    record_time_CG[i_S,j_repeat]=system.time(
      for(i_here in 1:1){
    
        sort_d_all=sort(d_all,index.return=T)
        sort_d_all$x=sort_d_all$x[-(1:(N*(L+1)))]
        sort_d_all$ix=sort_d_all$ix[-(1:(N*(L+1)))] ##delete zero
        unique_d_all_sorted=unique(sort_d_all$x)
        n_unique_d_all=length(unique_d_all_sorted)
        #n_unique_d_all-N*(N-1)/2*(L+1)
        
        Pr_Pc = Get_Pr_Pc(sort_d_all$ix,unique_d_all_sorted, L=L,n=N,d=d_all)
        permutation_A=Pr_Pc[[1]]
        p_r = Pr_Pc[[2]]
        ##
        delta_x_all = unique_d_all_sorted[-1] - unique_d_all_sorted[-n_unique_d_all]
        Q_K_2 = matrix(1.0, N*(N-1)*(L+1)/2, 1)
        
        
        rho = exp(-beta*delta_x_all)
        phi_L_T_inv = Get_phi_L_T_inv(rho, n_unique_d_all)
        D_input = phi_L_T_inv[[1]]
        E = phi_L_T_inv[[2]]
        
        gg <- c(0, rho)
        sqrt_Q = c(1,sqrt(1-rho^2))
        
        tol=10^{-6}*N^2
        

        ##not normalized
        res_LJ_multi_step = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=(v_all), gg,
                                                      Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
                                                      param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp", 
                                                      L,n=N,D ,tol = tol, maxIte = 1000)
        
        #res_LJ_multi_step[[2]]
        
        testing_n =100
        testing_d=as.numeric(seq(5/testing_n,5,5/(testing_n)))

        r0=abs(outer(testing_d,unique_d_all_sorted,'-'))
        r = exp(-beta*r0)
        pred_mean_CG_LJ_multi_step = r%*%t_sparse_A_times_x_multi(A_all, p_r, res_LJ_multi_step[[1]], L, n=N, D)
       

    })[3]
    
    testing_output=negative_LJ_truncated_kernel(testing_d)
    
    record_rmse_full[i_S,j_repeat]=sqrt(mean( (pred_mean_CG_LJ_multi_step-testing_output)^2))
    
    #plot(testing_d,testing_output,type='p',col='black',pch=19,main="CG_LJ_multi_steps")
    #lines(testing_d,pred_mean_CG_LJ_multi_step,type='p',col='blue',pch=20)
    
    
    if(i_S<=8){
      record_time_full[i_S,j_repeat]=system.time(
        for(i in 1:1){
          A_full=matrix(0,D*N,n_unique_d_all) ##put zero for now, we only need to record nonzero
          
          for(i in 1:N){
            #A_full[c(i,i+N),permutation_A[i,-i] ]=-t(cbind( (A[i,((1:N)[-i])]), (A[i+N,((1:N)[-i])])))
            A_full[c(i,i+N),permutation_A[i,-i] ]=t(cbind( (A[i,((1:N)[-i])]), (A[i+N,((1:N)[-i])])))
            
          }
          # 
          d0=abs(outer(unique_d_all_sorted,unique_d_all_sorted,'-'))
          
          phi_full=exp(-beta*d0)
          
          R=A_full%*%phi_full%*%t(A_full)
          
          r=matrix(NA,testing_n,N*D)
          
          for(i_testing in 1:testing_n){
            d_star=abs(outer(testing_d[i_testing],unique_d_all_sorted,'-'))
            

            r[i_testing,]=exp(-beta*d_star)%*%t(A_full) ##need sparse coding
          }
          
          v_all=as.matrix(v_all)          

          pred_mean_full=(r)%*%(solve(R+nu*diag(D*N))%*%v_all[,1]) ##times N to get back
          
        }
      )[3]
      
      record_rmse_full[i_S,j_repeat]=sqrt(mean( (pred_mean_full-testing_output)^2))    
      record_rmse_CG_full[i_S,j_repeat]=sqrt(mean( (pred_mean_CG_LJ_multi_step-pred_mean_full)^2))
      
    }
  }

}



time_full= rowMeans(record_time_full)

time_CG= rowMeans(record_time_CG)

pdf('time_full_GP_CG_GP.pdf',height=4,width=4)
plot(D*20*(1:S),time_full,type='p',col='red',pch=17,xlab='N',ylab='Time(s)',mgp=c(2.5,1,0))
lines(D*20*(1:S),time_CG,type='p',col='blue',pch=19,xlab='N',ylab='Time(s)',mgp=c(2.5,1,0))
legend("topright",
       pch=c(17,19),
       col=c('red','blue'),
       legend=c('Full GP', 'CG GP'))
dev.off()
save(D,S,time_full,time_CG,record_rmse_CG_full,record_rmse_full,record_rmse_CG,file='time_comparison_particle_interaction.RData')

  

        