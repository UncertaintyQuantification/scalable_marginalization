library(deSolve) 
library(spatPomp)
library(FastGaSP)
library(RobustGaSP)
library(nloptr)
library(Matrix)
library(Rcpp)
library(RcppEigen)

sourceCpp(file='src/functions_Feb26_2022.cpp')



OD_influence_kernel <- function(r){
  n = length(r)
  influence = rep(0, n)
  cutoff = 1/(sqrt(2))
  support = 1
  delta = 0.05
  for(i in 1:n){
    if(r[i] < (cutoff - delta)){
      influence[i] = 0.4
    } 
    if( (cutoff - delta) <= r[i] & r[i] < (cutoff + delta)){
      y_1 = 0.4
      y_2 = 1
      influence[i] = (y_2 - y_1)/(-2)*(cos(pi/(2 * delta)*(r[i] - (cutoff - delta))) - 1) + y_1
    } 
    if( (cutoff + delta) <= r[i] & r[i] < (support - delta)){
      influence[i] = 1
    } 
    if( (support - delta) <= r[i] & r[i] < (support + delta)){
      y_1 = 1
      y_2 = 0
      influence[i] = (y_2 - y_1)/(-2)*(cos(pi/(2 * delta)*(r[i] - (support - delta))) - 1) + y_1
    }
  }
  return(influence)
}



norm_vec <- function(v){
  sqrt(sum(v^2))
}



RHSfn_OD <- function(t,x,N){
  D = length(x)/N
  y = matrix(x,N,D)  #[x_1...x_N]  d x N
  aa= matrix(0,N,D)
  for(i in 1:N){
    #temp  = -(t(y)- matrix(y[i,],D,N))
    temp  = (t(y)- matrix(y[i,],D,N))
    DD  = apply(temp,2,norm_vec)      # 1 x N
    DD[DD==0] = 1                     #this step is to avoid 0 x inf= NAN
    aa[i,] = temp%*%OD_influence_kernel(DD)/N   # d x 1
  }
  dx = as.vector(aa)  # dN x 1 
  return( list(dx) )
  
}



RHSfn_pred<-function(t,x,params){
  N=params[1]
  weights=params[2:length(params)]
  
  D = length(x)/N
  y = matrix(x,N,D)  #[x_1...x_N]  d x N
  aa= matrix(0,N,D)
  
  
  
  for(i in 1:N){
    temp  = (t(y)- matrix(y[i,],D,N))    # d x N [x_1-x_i,..., x_N-x_i]
    
    DD  = apply(temp,2,norm_vec)      # 1 x N
    DD[DD==0] = 1                     #this step is to avoid 0 x inf= NAN
    
    r0=abs(outer(DD,unique_d_all_sorted,'-'))
    r_here = exp(-beta*r0)
    pred_mean_CG_OD_here = r_here%*%weights
    aa[i,] = (temp)%*%pred_mean_CG_OD_here/N
  }
  
  dx = as.vector(aa)  # dN x 1
  return( list(dx) )
  
  
}



Generate_OD_data <- function(sysInfo,obsInfo,solverInfo){
  # sysInfo = c(N, d)
  # obsInfo <- list(time_vec, M)
  N = sysInfo[1]         # number of agents
  D = sysInfo[2]        # dim of state vectors
  time_vec = obsInfo[[1]]
  M = obsInfo[[2]]

  xpath_train = array(0, c(D*N, length(time_vec)+1, M))
  dxpath_train = array(0, c(D*N, length(time_vec)+1, M))
  
  for(i in 1:M){    # trajectories with random initial conditions for learning interaction kernel
    if(obsInfo[[3]]=='normal'){
      x0 = rnorm(D*N,mean=0,sd=0.25)     # random initial condition
    }else if(obsInfo[[3]]=='uniform'){
      x0 = 1.5*runif(D*N)     # random initial condition
      
    }else if(obsInfo[[3]]=='log_uniform'){
      x0 = log(10^{-3})+(log(1.5)-log(10^{-3}) )*runif(D*N)     # random initial condition
      x0=exp(x0)
    }else if(obsInfo[[3]]=='test_design'){
      x0 = runif(D*N)     # random initial condition
    }
    
    output = ode(x0, times=c(0, time_vec), func=RHSfn_OD, parms=N, method =solverInfo)
    output = t(as.matrix(output[,-1])) ##keep initial value
    
    xpath_train[,,i] = output
    
    for(j in 1: (length(time_vec)+1)){
      dxpath_train[,j,i] = (RHSfn_OD(t=time_vec[j], x=xpath_train[,j,i], N=N))[[1]]
    }
  } 
  res = as.list(1:2)
  res[[1]] = xpath_train
  res[[2]] = dxpath_train
  
  return(res)
  
}


######################################################################################################
num_time_type=2    ##L=1 or L=10
num_particle_type=2 ##N=50 or N=200
num_design_type=3  

num_repetition=10  ##repetition for each simulation   

record_rmse_phi_CG=array(NA,c(num_time_type,num_particle_type,num_design_type))
record_rmse_phi_CG_full=array(NA,c(num_time_type,num_particle_type,num_design_type))


testing_n =1000

#record_phi_pred=matrix(NA,testing_n,S)

#record_phi_truth=matrix(NA,N,S)

obsInfo_all=c('uniform','normal','log_uniform')
N_Info_all=c(50,200)
L_Info_all=c(1,10)


##held-out testing 
testing_d=as.numeric(seq(1.5/testing_n,1.5,1.5/(testing_n)))
testing_output=OD_influence_kernel(testing_d)

##record the estimation of phi 
record_est_phi_CG=matrix(NA,num_time_type*num_particle_type*num_design_type*num_repetition,testing_n)

##let's focus on RMSE for now
count_record_est_phi=0;

for(i_L in 1:num_time_type){
  print(i_L)
  for(j_N in 1:num_particle_type){
    print(j_N)
    
    for(k_design in 1:num_design_type){
      
      sum_squared_error_phi_CG=0 
      sum_squared_error_phi_CG_full=0
      
      for(l_repetition in 1:num_repetition){
        count_record_est_phi=count_record_est_phi+1
        
        seed_num=(j_N-1)*num_particle_type*num_design_type*num_repetition+ (j_N-1)*num_design_type*num_repetition+(k_design-1)*num_repetition+l_repetition
        
        set.seed(seed_num)
        
        ###1.  generating data 
        L_sim=L_Info_all[i_L]      # time steps>1
        L=L_Info_all[i_L]-1        ##how many data to be use, only initial as 0
        
        N=N_Info_all[j_N]
        
        
        D=2
        #h = 0.05
        h = 0.1
        
        M = 1    # the number of trajectories
        
        ###test euler first, work on other approach later
        sysInfo = c(N, D) 
        solverInfo='euler'  ##just for demonstration purposes
        obsInfo = as.list(1:3)
        obsInfo[[1]] = seq(from=h, by=h, length.out=L_sim)
        obsInfo[[2]] = M
        obsInfo[[3]] = obsInfo_all[k_design]
        
        res = Generate_OD_data(sysInfo,obsInfo,solverInfo)
        input_pos = res[[1]][,,1] ##contains zero
        
        v_all = (res[[2]][,,1])[,1:L_sim]
        
        
        ###2.start prediction
        
        d_all = matrix(NA,N*(L+1),N)
        A_all = matrix(NA,D*N*(L+1),N)
        #A_r_all = matrix(0,D*N,(L+1)*N*(N-1)/2)
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
            #v3[c(i,i+N),1]= c(A[1:N,i]%*%OD_truncated_kernel(d[,i]),A[N+(1:N),i]%*%OD_truncated_kernel(d[,i]))
          }
          
          d_all[((t-1)*N+1):(t*N),] = d
        }
        
        
        #plot(input_pos[1:N,1],input_pos[N+1:N,1],pch=19,xlim=c(-2,7),ylim=c(-2,7))
        #lines(input_pos[1:N,2],input_pos[N+1:N,2],pch=18,col='lightblue',type='p')
        
        beta=.2
        nu  = 10^{-5}
        param = log(c(beta, nu))
        
        

        sort_d_all=sort(d_all,index.return=T)
        sort_d_all$x=sort_d_all$x[-(1:(N*(L+1)))]
        sort_d_all$ix=sort_d_all$ix[-(1:(N*(L+1)))] ##delete zero
        unique_d_all_sorted=unique(sort_d_all$x)
        n_unique_d_all=length(unique_d_all_sorted)
        n_unique_d_all-N*(N-1)/2*(L+1)
        
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
        res_OD_multi_step = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=(v_all*N), gg,
                                                      Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
                                                      param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp", 
                                                      L,n=N,D ,tol = tol, maxIte = 1000)
        

        
        r0=abs(outer(testing_d,unique_d_all_sorted,'-'))
        r = exp(-beta*r0)
        A_t_tilde_R_inv_y=t_sparse_A_times_x_multi(A_all, p_r, res_OD_multi_step[[1]], L, n=N, D)
        
        pred_mean_CG_OD_multi_step = r%*%A_t_tilde_R_inv_y
        sum_squared_error_phi_CG=sum_squared_error_phi_CG+sum( (pred_mean_CG_OD_multi_step-testing_output)^2 )
        
        
        ##record the full CG
        record_est_phi_CG[count_record_est_phi,]=pred_mean_CG_OD_multi_step
        
        if((i_L==1)&&(j_N==1)){
          A_full=matrix(0,D*N,n_unique_d_all) ##put zero for now, we only need to record nonzero
          
          for(i in 1:N){
            A_full[c(i,i+N),permutation_A[i,-i] ]=t(cbind( (A[i,((1:N)[-i])]), (A[i+N,((1:N)[-i])])))
          }
          
          d0=abs(outer(unique_d_all_sorted,unique_d_all_sorted,'-'))
          
          phi_full=exp(-beta*d0)
          
          R=A_full%*%phi_full%*%t(A_full)
          
          r_full=matrix(NA,testing_n,N*D)
          
          for(i_testing in 1:testing_n){
            d_star=abs(outer(testing_d[i_testing],unique_d_all_sorted,'-'))
            
            #A_i=cbind((A[1:N,i]),(A[N+(1:N),i]))
            
            r_full[i_testing,]=exp(-beta*d_star)%*%t(A_full) ##need sparse coding
          }
          
          v_all=as.matrix(v_all)          
          
          pred_mean_full=(r_full)%*%(solve(R+nu*diag(D*N))%*%v_all[,1])*N ##times N to get back
          
          sum_squared_error_phi_CG_full=sum_squared_error_phi_CG_full+sum( (pred_mean_CG_OD_multi_step-pred_mean_full)^2 )
        }
        
        
        
      }
      
      record_rmse_phi_CG[i_L,j_N,k_design]=sqrt(sum_squared_error_phi_CG/(length(testing_output)*num_repetition))
      if((i_L==1)&&(j_N==1)){
        record_rmse_phi_CG_full[i_L,j_N,k_design]=sqrt(sum_squared_error_phi_CG_full/(length(testing_output)*num_repetition))
        print(record_rmse_phi_CG_full[i_L,j_N,k_design])
        
      }
      print(record_rmse_phi_CG[i_L,j_N,k_design])
    }
    
    
  }
}





record_rmse_phi_CG/sd(testing_output)



####interval 

LB95=as.list(1:2)
UB95=as.list(1:2)

for(i_L in 1:2 ){
  j_N=1
  k_design=3
  l_repetition=1
  
  
  seed_num=(j_N-1)*num_particle_type*num_design_type*num_repetition+ (j_N-1)*num_design_type*num_repetition+(k_design-1)*num_repetition+l_repetition
  
  
  set.seed(seed_num)
  
  ###1.  generating data 
  L_sim=L_Info_all[i_L]      # time steps>1
  L=L_Info_all[i_L]-1        ##how many data to be use, only initial as 0
  
  N=N_Info_all[j_N]
  
  
  D=2
  h = 0.1
  
  M = 1    # the number of trajectories
  
  ###test euler first, work on other approach later
  sysInfo = c(N, D) 
  solverInfo='euler'  ##just for demonstration purposes
  obsInfo = as.list(1:3)
  obsInfo[[1]] = seq(from=h, by=h, length.out=L_sim)
  obsInfo[[2]] = M
  obsInfo[[3]] = obsInfo_all[k_design]
  
  res = Generate_OD_data(sysInfo,obsInfo,solverInfo)
  input_pos = res[[1]][,,1] ##contains zero
  
  v_all = (res[[2]][,,1])[,1:L_sim]
  
  
  ###2.start prediction
  
  d_all = matrix(NA,N*(L+1),N)
  A_all = matrix(NA,D*N*(L+1),N)
  #A_r_all = matrix(0,D*N,(L+1)*N*(N-1)/2)
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
      #v3[c(i,i+N),1]= c(A[1:N,i]%*%OD_truncated_kernel(d[,i]),A[N+(1:N),i]%*%OD_truncated_kernel(d[,i]))
    }
    
    d_all[((t-1)*N+1):(t*N),] = d
  }
  
  
  #plot(input_pos[1:N,1],input_pos[N+1:N,1],pch=19,xlim=c(-2,7),ylim=c(-2,7))
  #lines(input_pos[1:N,2],input_pos[N+1:N,2],pch=18,col='lightblue',type='p')
  
  beta=.2
  nu  = 10^{-5}
  param = log(c(beta, nu))
  
  
  sort_d_all=sort(d_all,index.return=T)
  sort_d_all$x=sort_d_all$x[-(1:(N*(L+1)))]
  sort_d_all$ix=sort_d_all$ix[-(1:(N*(L+1)))] ##delete zero
  unique_d_all_sorted=unique(sort_d_all$x)
  n_unique_d_all=length(unique_d_all_sorted)
  n_unique_d_all-N*(N-1)/2*(L+1)
  
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
  res_OD_multi_step = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=(v_all*N), gg,
                                                Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
                                                param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp", 
                                                L,n=N,D ,tol = tol, maxIte = 1000)
  
  #res_OD_multi_step[[2]]
  
  
  r0=abs(outer(testing_d,unique_d_all_sorted,'-'))
  r = exp(-beta*r0)
  A_t_tilde_R_inv_y=t_sparse_A_times_x_multi(A_all, p_r, res_OD_multi_step[[1]], L, n=N, D)
  
  pred_mean_CG_OD_multi_step = r%*%A_t_tilde_R_inv_y

  
  c_star=rep(NA,testing_n)
  system.time(
    for(i in 1:testing_n ){
      print(i)

      Ur=sparse_A_times_x_multi( A_all,  P=permutation_A, r[i,], L, N, D)

      tol=10^{-6}*(N*L)^2
      res_r = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=Ur, gg,
                                       Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
                                       param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp",
                                       L,n=N,D ,tol = tol, maxIte = 1000)
  
  
  
      
      r_R_inv_r = r[i,]%*%t_sparse_A_times_x_multi(A_all,p_r,res_r[[1]],L, n=N, D)
      
      c_star[i]=1-r_R_inv_r
    }
  )
  v_all=as.matrix(v_all)
  sigma_2_est = as.vector(v_all*N)%*%res_OD_multi_step[[1]]/length(v_all) ##why times N?
  

  LB95[[i_L]]=    pred_mean_CG_OD_multi_step+sqrt(as.numeric(sigma_2_est)*c_star)*qnorm(0.025)
  UB95[[i_L]]=    pred_mean_CG_OD_multi_step+sqrt(as.numeric(sigma_2_est)*c_star)*qnorm(0.975)
}

pdf('truncated_OD_n_200_L_1.pdf',height=4.3,width=4.3)
plot(testing_d,testing_output,type='l',col='black',
     lty=1,xlab='d',ylab=expression(phi),mgp=c(2.5,1,0),main='OD, L=1',ylim=c(min(LB95[[1]]),max(UB95[[1]])))
polygon( c(testing_d,rev(testing_d)), c(LB95[[1]],rev(UB95[[1]])), col = "grey80", border = F)
lines(testing_d,testing_output,type='l',col='black',lty=1)


lines(testing_d,record_est_phi_CG[1,],type='l',col='green',lty=2)
lines(testing_d,record_est_phi_CG[num_repetition+1,],type='l',col='orange',lty=2)
lines(testing_d,record_est_phi_CG[2*num_repetition+1,],type='l',col='blue',lty=2)
lines(testing_d,record_est_phi_CG[3*num_repetition+1,],type='l',col='green',lty=3)
lines(testing_d,record_est_phi_CG[4*num_repetition+1,],type='l',col='orange',lty=3)
lines(testing_d,record_est_phi_CG[5*num_repetition+1,],type='l',col='blue',lty=3)

dev.off()


pdf('truncated_OD_n_200_L_10.pdf',height=4.3,width=4.3)
plot(testing_d,testing_output,type='l',col='black',
     lty=1,xlab='d',ylab=expression(phi),mgp=c(2.5,1,0),main='OD, L=10',ylim=c(min(LB95[[2]]),max(UB95[[2]])))
polygon( c(testing_d,rev(testing_d)), c(LB95[[2]],rev(UB95[[2]])), col = "grey80", border = F)
lines(testing_d,testing_output,type='l',col='black',lty=1)


lines(testing_d,record_est_phi_CG[6*num_repetition+1,],type='l',col='green',lty=2)
lines(testing_d,record_est_phi_CG[7*num_repetition+1,],type='l',col='orange',lty=2)
lines(testing_d,record_est_phi_CG[8*num_repetition+1,],type='l',col='blue',lty=2)
lines(testing_d,record_est_phi_CG[9*num_repetition+1,],type='l',col='green',lty=3)
lines(testing_d,record_est_phi_CG[10*num_repetition+1,],type='l',col='orange',lty=3)
lines(testing_d,record_est_phi_CG[11*num_repetition+1,],type='l',col='blue',lty=3)

dev.off()





##3. predicting velocity, this only for plotting one figure, we don't need to show comparison
##            
set.seed(1)

L_sim=L_Info_all[2]      # time steps>1
L=L_Info_all[2]-1        ##how many data to be use, only initial as 0

N=N_Info_all[1]
k_design=3

D=1
#h = 0.05
h = 0.1

M = 1    # the number of trajectories

###test euler first, work on other approach later
sysInfo = c(N, D) 
solverInfo='euler' ##just for demonstration purposes
obsInfo = as.list(1:3)
obsInfo[[1]] = seq(from=h, by=h, length.out=L_sim)
obsInfo[[2]] = M
obsInfo[[3]] = obsInfo_all[k_design]

res = Generate_OD_data(sysInfo,obsInfo,solverInfo)
input_pos = res[[1]][,,1] ##contains zero

v_all = (res[[2]][,,1])[,1:L_sim]


d_all = matrix(NA,N*(L+1),N)
A_all = matrix(NA,D*N*(L+1),N)
#A_r_all = matrix(0,D*N,(L+1)*N*(N-1)/2)
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
    #v3[c(i,i+N),1]= c(A[1:N,i]%*%OD_truncated_kernel(d[,i]),A[N+(1:N),i]%*%OD_truncated_kernel(d[,i]))
  }
  
  d_all[((t-1)*N+1):(t*N),] = d
}


beta=0.2
nu = 10^{-5}
param = log(c(beta, nu^2))

sort_d_all=sort(d_all,index.return=T)
sort_d_all$x=sort_d_all$x[-(1:(N*(L+1)))]
sort_d_all$ix=sort_d_all$ix[-(1:(N*(L+1)))] ##delete zero
unique_d_all_sorted=unique(sort_d_all$x)
n_unique_d_all=length(unique_d_all_sorted)
n_unique_d_all-N*(N-1)/2*(L+1)

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
res_OD_multi_step = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=(v_all*N), gg,
                                              Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
                                              param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp", 
                                              L,n=N,D ,tol = tol, maxIte = 1000)

#res_OD_multi_step[[2]]


r0=abs(outer(testing_d,unique_d_all_sorted,'-'))
r = exp(-beta*r0)

###this is the only thing we need
A_t_tilde_R_inv_y=t_sparse_A_times_x_multi(A_all, p_r, res_OD_multi_step[[1]], L, n=N, D)
weights=A_t_tilde_R_inv_y

###generating testing and prediction 
test_obsInfo = as.list(1:3)

L_testing=100
test_obsInfo[[1]] = seq(from=h, by=h, length.out=L_testing)
test_obsInfo[[2]] = M
test_obsInfo[[3]] = 'test_design'

test_res = Generate_OD_data(sysInfo=sysInfo,obsInfo=test_obsInfo,solverInfo=solverInfo)
test_input_pos = test_res[[1]][,,1] ##contains zero




x0=test_input_pos[,1] ### the initial velocity
params=c(N,weights)
time_vec = test_obsInfo[[1]]

####this could be very slow, we use RHSfn_pred based on our method 
predict_input_pos = ode(x0, times=c(0, time_vec), func=RHSfn_pred, parms=params, method =solverInfo) ##let me just try euler first

predict_input_pos = t(as.matrix(predict_input_pos[,-1])) ##keep initial value


##second step and last step
predict_input_pos[,2]-test_input_pos[,2]

predict_input_pos[,L_testing+1]-test_input_pos[,L_testing+1]



pdf('prediction_trajectory_truncated_OD_n_200.pdf',height=4.3,width=4.3)

particle_num=5
plot(t(test_input_pos[c(particle_num,particle_num+N),]),type='l',
     xlim=c(0,1),ylim=c(0.2,1),lty=2,col='black',
     xlab=expression(x[1]),ylab=expression(x[2]),main='OD')
lines(t(predict_input_pos[c(particle_num,particle_num+N),]),type='l',lty=4,col='blue')
lines(t(test_input_pos[c(particle_num,particle_num+N),1]),pch=20,type='p',cex=.5)

abs_v=sqrt((test_input_pos[c(particle_num),L_testing+1]-test_input_pos[c(particle_num),L_testing])^2+ (test_input_pos[c(particle_num+N),L_testing+1]-test_input_pos[c(particle_num+N),L_testing])^2)

for(i_num in 2:10){
  particle_num=(i_num-1)*5+5
  
  lines(t(test_input_pos[c(particle_num,particle_num+N),]),type='l',lty=2,col='black')
  lines(t(predict_input_pos[c(particle_num,particle_num+N),]),type='l',lty=4,col='blue')
  lines(t(test_input_pos[c(particle_num,particle_num+N),1]),pch=20,type='p',cex=.5)
  abs_v=sqrt((test_input_pos[c(particle_num),L_testing+1]-test_input_pos[c(particle_num),L_testing])^2+ (test_input_pos[c(particle_num+N),L_testing+1]-test_input_pos[c(particle_num+N),L_testing])^2)

}
legend("bottomleft",lty=c(2,4),
       col=c('black','blue'),
       legend=c( 'Truth','Prediction')
)
dev.off()


