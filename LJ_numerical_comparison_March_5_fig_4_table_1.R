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



negLJ_kernel<- function(r){
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

RHSfn_pred<-function(t,x,params){
  N=params[1]
  weights=params[2:length(params)]
  
  D = length(x)/N
  y = matrix(x,N,D)  #[x_1...x_N]  d x N
  #bb= matrix(0,N,D)
  aa= matrix(0,N,D)
  
  
  
  ##change this to current order
  for(i in 1:N){
    temp  = (t(y)- matrix(y[i,],D,N))    # d x N [x_1-x_i,..., x_N-x_i]
    
    DD  = apply(temp,2,norm_vec)      # 1 x N
    DD[DD==0] = 1                     #this step is to avoid 0 x inf= NAN
    
    r0=abs(outer(DD,unique_d_all_sorted,'-'))
    r_here = exp(-beta*r0)
    pred_mean_CG_LJ_here = r_here%*%weights
    aa[i,] = (temp)%*%pred_mean_CG_LJ_here
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
  time_vec = obsInfo[[1]]
  M = obsInfo[[2]]
  xpath_train = array(0, c(D*N, length(time_vec)+1, M))
  dxpath_train = array(0, c(D*N, length(time_vec)+1, M))
  

  for(i in 1:M){    # trajectories with random initial conditions for learning interaction kernel
    if(obsInfo[[3]]=='normal'){
      x0 = rnorm(D*N)     # random initial condition
    }else if(obsInfo[[3]]=='uniform'){
      x0 = 5*runif(D*N)     # random initial condition
      
    }else if(obsInfo[[3]]=='log_uniform'){
      x0 = log(10^{-3})+(log(5)-log(10^{-3}) )*runif(D*N)     # random initial condition
      x0=exp(x0)
    }else if(obsInfo[[3]]=='test_design'){
      x0 = 1*runif(D*N)     # random initial condition
    }
    # could still be large
    #x0 =5*runif(d*N)     # random initial condition for truncated LJ
    
    # x0 = as.numeric(seq(5/(d*N),5,5/(d*N)))    # designed initial condition
    #output = ode(x0, times=c(0, time_vec), func=RHSfn, parms=N, method = "rk4")
    output = ode(x0, times=c(0, time_vec), func=RHSfn, parms=N, method =solverInfo) ##let me just try euler first
    
    #output = ( t(as.matrix(output[,-1])) )[,-1]
    output = t(as.matrix(output[,-1])) ##keep initial value
    
    xpath_train[,,i] = output
    
    for(j in 1: (length(time_vec)+1)){
      #dxpath_train[,j,i] = (RHSfn(t=time_vec[j], x=xpath_train[,j,i], N=N))[[1]]
      dxpath_train[,j,i] = (RHSfn(t=time_vec[j], x=xpath_train[,j,i], N=N))[[1]]
    }
  }
  res = as.list(1:2)
  res[[1]] = xpath_train
  res[[2]] = dxpath_train
  
  return(res)
}





##############three design 
###let's do e.g. N=50 and N=200 particles, then 1 step and 10 steps,
#based on these 3 designs
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
testing_d=as.numeric(seq(5/testing_n,5,5/(testing_n)))
testing_output=negative_LJ_truncated_kernel(testing_d)

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
            h = 5*10^{-4}
            
            M = 1    # the number of trajectories
            
            ###test euler first, work on other approach later
            sysInfo = c(N, D) 
            solverInfo='euler' ##just for demonstration purposes
            obsInfo = as.list(1:3)
            obsInfo[[1]] = seq(from=h, by=h, length.out=L_sim)
            obsInfo[[2]] = M
            obsInfo[[3]] = obsInfo_all[k_design]
            
            res = Generate_training_data(sysInfo,obsInfo,solverInfo)
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
                #v3[c(i,i+N),1]= c(A[1:N,i]%*%LJ_truncated_kernel(d[,i]),A[N+(1:N),i]%*%LJ_truncated_kernel(d[,i]))
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
            
      
            r0=abs(outer(testing_d,unique_d_all_sorted,'-'))
            r = exp(-beta*r0)
            A_t_tilde_R_inv_y=t_sparse_A_times_x_multi(A_all, p_r, res_LJ_multi_step[[1]], L, n=N, D)
            
            pred_mean_CG_LJ_multi_step = r%*%A_t_tilde_R_inv_y
            sum_squared_error_phi_CG=sum_squared_error_phi_CG+sum( (pred_mean_CG_LJ_multi_step-testing_output)^2 )
            
            ##record the full CG
            record_est_phi_CG[count_record_est_phi,]=pred_mean_CG_LJ_multi_step
          
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
                  
                  #pred_mean_full=(r)%*%(solve(R+nu*diag(D*N))%*%v_all[,1])*N ##times N to get back
                  #not normalized
                  pred_mean_full=(r_full)%*%(solve(R+nu*diag(D*N))%*%v_all[,1]) ##times N to get back
                  
                  sum_squared_error_phi_CG_full=sum_squared_error_phi_CG_full+sum( (pred_mean_CG_LJ_multi_step-pred_mean_full)^2 )
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

save(record_rmse_phi_CG,record_rmse_phi_CG_full,testing_d,record_est_phi_CG,testing_output,file='numerical_comparison_particle_interaction_LJ.RData')


#plot(pred_mean_CG_LJ_multi_step,ylim=c(-12,1))
#lines(testing_output)



### to compute the uncertainty for L>1
# c_star=rep(NA,testing_n)
# system.time(
#   for(i in 1:testing_n ){
#     print(i)
#     Ur=sparse_A_times_x( A_all,   P=permutation_A, r[i,])
#     A_t_tilde_R_inv_y=t_sparse_A_times_x_multi(A_all, p_r, res_LJ_multi_step[[1]], L, n=N, D)

#     ##Ur=t_sparse_A_times_x(A,r_index,r_star[i,],N)
#     tol=10^{-6}*(N*L)^2
#     res_r = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=Ur, gg,
#                                      Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
#                                      param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp",
#                                      L,n=N,D ,tol = tol, maxIte = 1000)
# 
# 
#     r_R_inv_r = r[i,]%*%t_sparse_A_times_x(A,p_r,res_r[[1]],N)
#     c_star[i]=1-r_R_inv_r
#   }
# )
# v_all=as.matrix(v_all)
# sigma_2_est = t(v_all[,1])%*%res_LJ_multi_step[[1]]/length(v_all[,1])
# 
# LB95=rep(NA,testing_n)
# UB95=rep(NA,testing_n)
# 
# LB95=    record_phi_pred[,3]+sqrt(as.numeric(sigma_2_est)*c_star)*qnorm(0.025)
# UB95=    record_phi_pred[,3]+sqrt(as.numeric(sigma_2_est)*c_star)*qnorm(0.975)





record_rmse_phi_CG/sd(testing_output)



####interval 
#count_record_est_phi=count_record_est_phi+1
##record 2
LB95=as.list(1:2)
UB95=as.list(1:2)

for(i_L in 1:2 ){
  #i_L=1  ####
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
  
  res = Generate_training_data(sysInfo,obsInfo,solverInfo)
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
      #v3[c(i,i+N),1]= c(A[1:N,i]%*%LJ_truncated_kernel(d[,i]),A[N+(1:N),i]%*%LJ_truncated_kernel(d[,i]))
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
  res_LJ_multi_step = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=(v_all), gg,
                                                Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
                                                param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp", 
                                                L,n=N,D ,tol = tol, maxIte = 1000)
  
  #res_LJ_multi_step[[2]]
  
  
  r0=abs(outer(testing_d,unique_d_all_sorted,'-'))
  r = exp(-beta*r0)
  A_t_tilde_R_inv_y=t_sparse_A_times_x_multi(A_all, p_r, res_LJ_multi_step[[1]], L, n=N, D)
  
  pred_mean_CG_LJ_multi_step = r%*%A_t_tilde_R_inv_y
  #sum_squared_error_phi_CG=sum_squared_error_phi_CG+sum( (pred_mean_CG_LJ_multi_step-testing_output)^2 )
  
  

  c_star=rep(NA,testing_n)
  system.time(
    for(i in 1:testing_n ){
      z#print(i)
      #Ur=sparse_A_times_x( A_all,   P=permutation_A, r[i,])
      
      Ur=sparse_A_times_x_multi( A_all,  P=permutation_A, r[i,], L, N, D)
      #Ur=as.matrix(Ur)
      
      ##Ur=t_sparse_A_times_x(A,r_index,r_star[i,],N)
      tol=10^{-6}*(N*L)^2
      res_r = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=Ur, gg,
                                        Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
                                        param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp",
                                        L,n=N,D ,tol = tol, maxIte = 1000)
      
      
      
      
      #r_R_inv_r = r[i,]%*%t_sparse_A_times_x(A_all,p_r,res_r[[1]],N)
      r_R_inv_r = r[i,]%*%t_sparse_A_times_x_multi(A_all,p_r,res_r[[1]],L, n=N, D)
      
      c_star[i]=1-r_R_inv_r
    }
  )
  v_all=as.matrix(v_all)
  #sigma_2_est = t(v_all[,1])%*%res_LJ_multi_step[[1]]/length(v_all[,1])
  #sigma_2_est = as.vector(v_all*N)%*%res_LJ_multi_step[[1]]/length(v_all) ##OD
  sigma_2_est = as.vector(v_all)%*%res_LJ_multi_step[[1]]/length(v_all) ##LJ
  
  #LB95=rep(NA,testing_n)
  #UB95=rep(NA,testing_n)
  
  LB95[[i_L]]=    pred_mean_CG_LJ_multi_step+sqrt(as.numeric(sigma_2_est)*c_star)*qnorm(0.025)
  UB95[[i_L]]=    pred_mean_CG_LJ_multi_step+sqrt(as.numeric(sigma_2_est)*c_star)*qnorm(0.975)
}
###

pdf('truncated_LJ_n_200_L_1.pdf',height=4.3,width=4.3)
plot(testing_d,testing_output,type='l',col='black',
     lty=1,xlab='d',ylab=expression(phi),mgp=c(2.5,1,0),main='Truncated LJ, L=1',ylim=c(min(LB95[[1]]),max(UB95[[1]])))
polygon( c(testing_d,rev(testing_d)), c(LB95[[1]],rev(UB95[[1]])), col = "grey80", border = F)
lines(testing_d,testing_output,type='l',col='black',lty=1)


lines(testing_d,record_est_phi_CG[1,],type='l',col='green',lty=2)
lines(testing_d,record_est_phi_CG[num_repetition+1,],type='l',col='orange',lty=2)
lines(testing_d,record_est_phi_CG[2*num_repetition+1,],type='l',col='blue',lty=2)
lines(testing_d,record_est_phi_CG[3*num_repetition+1,],type='l',col='green',lty=3)
lines(testing_d,record_est_phi_CG[4*num_repetition+1,],type='l',col='orange',lty=3)
lines(testing_d,record_est_phi_CG[5*num_repetition+1,],type='l',col='blue',lty=3)


legend("bottomright",lty=c(2,2,2,3,3,3,1),
       col=c('green','orange','blue','green','orange','blue','black'),
       legend=c( 'uniform, n=50','normal, n=50','log-uniform, n=50',
                 'uniform, n=200','normal, n=200','log-uniform, n=200',
                 'Truth')
)
dev.off()


pdf('truncated_LJ_n_200_L_10.pdf',height=4.3,width=4.3)
plot(testing_d,testing_output,type='l',col='black',
     lty=1,xlab='d',ylab=expression(phi),mgp=c(2.5,1,0),main='Truncated LJ, L=10',ylim=c(min(LB95[[2]]),max(UB95[[2]])))
polygon( c(testing_d,rev(testing_d)), c(LB95[[2]],rev(UB95[[2]])), col = "grey80", border = F)
lines(testing_d,testing_output,type='l',col='black',lty=1)

lines(testing_d,record_est_phi_CG[6*num_repetition+1,],type='l',col='green',lty=2)
lines(testing_d,record_est_phi_CG[7*num_repetition+1,],type='l',col='orange',lty=2)
lines(testing_d,record_est_phi_CG[8*num_repetition+1,],type='l',col='blue',lty=2)
lines(testing_d,record_est_phi_CG[9*num_repetition+1,],type='l',col='green',lty=3)
lines(testing_d,record_est_phi_CG[10*num_repetition+1,],type='l',col='orange',lty=3)
lines(testing_d,record_est_phi_CG[11*num_repetition+1,],type='l',col='blue',lty=3)

#legend("bottomright",lty=c(2,2,2,3,3,3,1),
#       col=c('green','orange','blue','green','orange','blue','black'),
#       legend=c( 'uniform, n=50','normal, n=50','log-uniform, n=50',
#                 'uniform, n=200','normal, n=200','log-uniform, n=200',
#                 'Truth'))
dev.off()




##3. predicting velocity, this only for plotting one figure, we don't need to show comparison
##            
set.seed(1)

L_sim=L_Info_all[2]      # time steps>1
L=L_Info_all[2]-1        ##how many data to be use, only initial as 0

N=N_Info_all[1]
k_design=3

D=2
#h = 0.05
h = 5*10^{-4}

M = 1    # the number of trajectories

###test euler first, work on other approach later
sysInfo = c(N, D) 
solverInfo='euler' ##just for demonstration purposes
obsInfo = as.list(1:3)
obsInfo[[1]] = seq(from=h, by=h, length.out=L_sim)
obsInfo[[2]] = M
obsInfo[[3]] = obsInfo_all[k_design]

res = Generate_training_data(sysInfo,obsInfo,solverInfo)
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
res_LJ_multi_step = CG_sparse_KF_Thomas_multi(A=A_all,diag_phi_L_T_inv=D_input,off_phi_L_T_inv=E,b=(v_all), gg,
                                              Q_K_2, sqrt_Q,  P = permutation_A, P_r = p_r,
                                              param=param, delta_x=delta_x_all,have_noise=F, kernel_type="exp", 
                                              L,n=N,D ,tol = tol, maxIte = 1000)



r0=abs(outer(testing_d,unique_d_all_sorted,'-'))
r = exp(-beta*r0)

###this is the only thing we need
A_t_tilde_R_inv_y=t_sparse_A_times_x_multi(A_all, p_r, res_LJ_multi_step[[1]], L, n=N, D)
weights=A_t_tilde_R_inv_y

###generating testing and prediction 
test_obsInfo = as.list(1:3)

L_testing=200
test_obsInfo[[1]] = seq(from=h, by=h, length.out=L_testing)
test_obsInfo[[2]] = M
test_obsInfo[[3]] = 'test_design'

test_res = Generate_training_data(sysInfo=sysInfo,obsInfo=test_obsInfo,solverInfo=solverInfo)
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


pdf('prediction_trajectory_truncated_LJ_n_200.pdf',height=4.3,width=4.3)

particle_num=5
plot(t(test_input_pos[c(particle_num,particle_num+N),]),type='l',
     xlim=c(-0.55,1.55),ylim=c(-1,2),lty=2,col='black',
     xlab=expression(x[1]),ylab=expression(x[2]),main='Truncated LJ')
lines(t(predict_input_pos[c(particle_num,particle_num+N),]),type='l',
     xlim=c(-0.55,1.55),ylim=c(-1,2),lty=4,col='blue')
lines(t(test_input_pos[c(particle_num,particle_num+N),1]),pch=20,type='p',cex=.5)

abs_v=sqrt((test_input_pos[c(particle_num),L_testing+1]-test_input_pos[c(particle_num),L_testing])^2+ (test_input_pos[c(particle_num+N),L_testing+1]-test_input_pos[c(particle_num+N),L_testing])^2)
arrows(x0=test_input_pos[c(particle_num),L_testing],y0=test_input_pos[c(particle_num+N),L_testing],
       x1=test_input_pos[c(particle_num),L_testing]+(test_input_pos[c(particle_num),L_testing+1]-test_input_pos[c(particle_num),L_testing])/abs_v*0.02,
       y1=test_input_pos[c(particle_num+N),L_testing]+(test_input_pos[c(particle_num+N),L_testing+1]-test_input_pos[c(particle_num+N),L_testing])/abs_v*0.02,
       length=0.03
)

for(i_num in 2:10){
  particle_num=(i_num-1)*5+5

  lines(t(test_input_pos[c(particle_num,particle_num+N),]),type='l',lty=2,col='black')
  lines(t(predict_input_pos[c(particle_num,particle_num+N),]),type='l',lty=4,col='blue')
  lines(t(test_input_pos[c(particle_num,particle_num+N),1]),pch=20,type='p',cex=.5)
  abs_v=sqrt((test_input_pos[c(particle_num),L_testing+1]-test_input_pos[c(particle_num),L_testing])^2+ (test_input_pos[c(particle_num+N),L_testing+1]-test_input_pos[c(particle_num+N),L_testing])^2)
  arrows(x0=test_input_pos[c(particle_num),L_testing],y0=test_input_pos[c(particle_num+N),L_testing],
         x1=test_input_pos[c(particle_num),L_testing]+(test_input_pos[c(particle_num),L_testing+1]-test_input_pos[c(particle_num),L_testing])/abs_v*0.02,
         y1=test_input_pos[c(particle_num+N),L_testing]+(test_input_pos[c(particle_num+N),L_testing+1]-test_input_pos[c(particle_num+N),L_testing])/abs_v*0.02,
         length=0.03
  )
  
}
legend("bottomleft",lty=c(2,4),
       col=c('black','blue'),
       legend=c( 'Truth','Prediction')
)
dev.off()
