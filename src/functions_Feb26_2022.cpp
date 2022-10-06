
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 

// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath> 
#include "ctools.h"
// [[Rcpp::depends(RcppEigen)]] 

using namespace Rcpp;
using namespace std;
using namespace Eigen; 



////L_times_z for exponential kernel
// [[Rcpp::export]] 
VectorXd L_times_z_CG(const Eigen::VectorXd gg,  const Eigen::VectorXd Q_K_2,  const Eigen::VectorXd sqrt_Q, 
                      const Eigen::VectorXd z){
  
  int n=z.rows();
  
  Eigen::MatrixXd K=Q_K_2;
  
  Eigen::VectorXd m=Eigen::VectorXd::Zero(n);
  double a;
  

  Eigen::VectorXd output=Eigen::VectorXd::Zero(n);
  for(int t=0;t<n;t++){
    a=gg[t]*m[0];
    output[t]=a+z[t]*sqrt_Q[t];
    m[0]=output[t];
  }
  
  return output;
}


// [[Rcpp::export]] 
List Get_phi_L_T_inv(const Eigen::VectorXd rho, const int n_unique){
  
  double temp = 1.0;
  List res(2);
  VectorXd diag=Eigen::VectorXd::Zero(n_unique);
  VectorXd off_diag=Eigen::VectorXd::Zero(n_unique-1);
  diag[0] = 1.0;
  for(int i=0;i<n_unique-1;i++){
    temp = 1.0/(1-pow(rho[i],2.0));
    diag[i+1] = sqrt(temp);
    temp = 1.0-pow(rho[i],2.0);
    off_diag[i] = -rho[i]/sqrt(temp);
  }
  res[1] = off_diag;
  res[0] = diag;
  return res; 
}


//Feb 16, 2022
// [[Rcpp::export]] 
MatrixXd Get_permutation_A(const VectorXi sort_d_all_ix, 
                           const VectorXd unique_d_sorted, 
                           int n, const VectorXd d){ 
  //n is the number of interacting particles here
  int index_count=0;  
  int index_here;
  MatrixXd permutation_A=MatrixXd::Zero(n,n); 
  int j_here;
  int i_here;
  for(int i=0; i< (pow(n,2)-n);i++){
    
    index_here=sort_d_all_ix[i]-1;
    
    if(d[index_here]>unique_d_sorted[index_count]){
      index_count=index_count+1;
    }
    
    
    if(d[index_here]==unique_d_sorted[index_count]){
      
      j_here=floor((sort_d_all_ix[i]-1)/n)+1;
      i_here=sort_d_all_ix[i]-(j_here-1)*n;
      
      permutation_A(i_here-1,j_here-1)=index_count+1; //this is the index to R
    }
    
    
  }
  
  return permutation_A; 
}

// [[Rcpp::export]] 
List Get_Pr_Pc(const VectorXi sort_d_all_ix, const VectorXd unique_d_sorted, const int L,
               int n, const VectorXd d){ 
  //n is the number of interacting particles here
  List res(2);
  int index_count=0;  
  int index_here;
  MatrixXd P_c=MatrixXd::Zero(n*(L+1),n); 
  int j_here;
  int i_here;
  MatrixXd P_r=MatrixXd::Zero(n*(n-1)*(L+1)/2,2);
  for(int i=0; i< (n*(n-1)*(L+1));i++){
    index_here=sort_d_all_ix[i]-1;
    
    if(d[index_here]>unique_d_sorted[index_count]){
      index_count=index_count+1;
    }
    
    
    if(d[index_here]==unique_d_sorted[index_count]){
      
      j_here=floor((sort_d_all_ix[i]-1)/(n*(L+1)))+1;
      i_here=sort_d_all_ix[i]-(j_here-1)*n*(L+1);
      
      P_c(i_here-1,j_here-1)=index_count+1; //this is the index to R
      if((i_here-floor((i_here-1)/n)*n)>j_here){
        P_r(index_count,0)=i_here;
      }else{
        P_r(index_count,1)=i_here;
      }
    }
  }
  res[0]=P_c;
  res[1]=P_r;//P_r.rightCols(2);
  return res; 
}

// [[Rcpp::export]] 
List Get_permutation_A_r_index(const VectorXi sort_d_all_ix, 
                               const VectorXd unique_d_sorted, 
                               int n, const VectorXd d){ 
  //n is the number of interacting particles here
  List res(2);
  int index_count=0;  
  int index_here;
  MatrixXd permutation_A=MatrixXd::Zero(n,n); 
  int j_here;
  int i_here;
  MatrixXd r_index=MatrixXd::Zero((pow(n,2)-n)/2,2+1);
  for(int i=0; i< (pow(n,2)-n);i++){
    
    index_here=sort_d_all_ix[i]-1;
    
    if(d[index_here]>unique_d_sorted[index_count]){
      index_count=index_count+1;
    }
    
    
    if(d[index_here]==unique_d_sorted[index_count]){
      
      j_here=floor((sort_d_all_ix[i]-1)/n)+1;
      i_here=sort_d_all_ix[i]-(j_here-1)*n;
      
      permutation_A(i_here-1,j_here-1)=index_count+1; //this is the index to R
      r_index(index_count,0)+=1;
      r_index(index_count,r_index(index_count,0))=i_here;
      ////A_r_index[[index_count]] = c(A_r_index[[index_count]],i_here)
      
    }
    
    
  }
  res[0]=permutation_A;
  res[1]=r_index.rightCols(2);
  return res; 
}

////Implement Thomas algorithm 
// [[Rcpp::export]] 
VectorXd Thomas_for_CG(const Eigen::MatrixXd A, const Eigen::VectorXd Y){
  int n=A.rows();
  //double w=0;
  Eigen::VectorXd X=Eigen::VectorXd::Zero(n); 
  Eigen::VectorXd y=Y; 
  Eigen::MatrixXd a=A;
  //for(int i=1; i<n;i++){
    //w = a(i,i-1)/a(i-1,i-1);
    //a(i,i) = a(i,i) - w*a(i-1,i);
    //y[i] = y[i] - w*y[i-1];
 // }
  X[n-1] = y[n-1]/a(n-1,n-1);
  for(int i=n-2; i>=0;i--){
    X[i] = (y[i] - a(i,i+1)*X[i+1])/a(i,i);
  }
  
  return X;
}


////Implement Thomas algorithm 
// [[Rcpp::export]] 
VectorXd Fast_Thomas_for_CG(const Eigen::VectorXd D, const Eigen::VectorXd E, const Eigen::VectorXd Y){
  int n=Y.size();
  //double w=0;
  Eigen::VectorXd X=Eigen::VectorXd::Zero(n); 
  Eigen::VectorXd d=D; 
  Eigen::VectorXd e=E;
  Eigen::VectorXd y=Y;

  X[n-1] = y[n-1]/d[n-1];
  for(int i=n-2; i>=0;i--){
    X[i] = (y[i] - e[i]*X[i+1])/d[i];
  }
  
  return X;
}


// [[Rcpp::export]] 
VectorXd sparse_A_times_x(const Eigen::MatrixXd A, const Eigen::MatrixXd P, const Eigen::VectorXd Y){
//// now this only works for 2D particles
  int n=A.rows();
  int N=P.rows();
  int index=0;
  Eigen::VectorXd res=Eigen::VectorXd::Zero(n); 
 //Eigen::VectorXd X=Eigen::VectorXd::Zero(n);
 //  Eigen::VectorXd r=Eigen::VectorXd::Zero(n-2); 
 // Eigen::VectorXd y=Y; 
 // Eigen::MatrixXd a=A;
  double res_1=0.0;
  double res_2=0.0;
  
  for(int i=0; i<N;i++){
    res_1=0.0;
    res_2=0.0;
    for(int j=0; j<N;j++){
      if(i!=j){
        index = P(i,j)-1;
        res_1 = res_1 + A(j,i)*Y[index];
        res_2 = res_2 + A(j+N,i)*Y[index];
      }
    }
    res[i] = res_1;
    res[i+N] = res_2;
  }
  return res;
}


// [[Rcpp::export]] 
VectorXd sparse_A_times_x_multi(const Eigen::MatrixXd A, const Eigen::MatrixXd P, const Eigen::VectorXd Y, 
                                const int L, int n, int D){
  //// now this only works for 2D particles
  int n_tilde=A.rows();
  int N=P.rows();
  int index=0;
  int t = 0;
  Eigen::VectorXd res=Eigen::VectorXd::Zero(n_tilde); 
  
  double res_1=0.0;
  double res_2=0.0;
  
  for(int i=0; i<N;i++){
    res_1=0.0;
    res_2=0.0;
    t = floor(i/n);
    for(int j=0; j<n;j++){
      if((i-t*n)!=j){
        index = P(i,j)-1;
        res_1 = res_1 + A((i-t*n+D*n*t),j)*Y[index];
        res_2 = res_2 + A((i-t*n+D*n*t+n),j)*Y[index];
      }
    }
    res[i+t*n] = res_1;
    res[i+t*n+n] = res_2;
  }
  return res;
}


// [[Rcpp::export]] 
VectorXd t_sparse_A_times_x(const Eigen::MatrixXd A, const Eigen::MatrixXd r_index, 
                            const Eigen::VectorXd Y, int n){
  int R = r_index.rows(); //R = N(N-1)/2
  Eigen::VectorXd res=Eigen::VectorXd::Zero(R); 
  for(int i=0; i<R;i++){
    res[i]=A((r_index(i,0)-1),(r_index(i,1)-1))*(Y[(r_index(i,1)-1)]-Y[(r_index(i,0)-1)]);
    res[i]+=A((r_index(i,0)-1+n),(r_index(i,1)-1))*(Y[(r_index(i,1)-1+n)]-Y[(r_index(i,0)-1+n)]);
  }
  return res;
}


// [[Rcpp::export]] 
VectorXd t_sparse_A_times_x_multi(const Eigen::MatrixXd A, const Eigen::MatrixXd P_r, 
                                  const Eigen::VectorXd Y, const int L, int n, int D){
  
  int R = P_r.rows(); //R = n(n-1)(L+1)/2
  int t = 0;
  Eigen::VectorXd res=Eigen::VectorXd::Zero(R); 
  for(int i=0; i<R;i++){
    t = floor((P_r(i,0)-1)/n);
    res[i]=A((P_r(i,0)-t*n+D*n*t-1),(P_r(i,1)-t*n-1))*(-Y[(P_r(i,1)+t*n-1)]+Y[(P_r(i,0)+t*n-1)]);
    res[i]+=A((P_r(i,0)-t*n+D*n*t-1+n),(P_r(i,1)-t*n-1))*(-Y[(P_r(i,1)+t*n-1+n)]+Y[(P_r(i,0)+t*n-1+n)]);
  }
  return res;
}








////CG with Kalman Filter and Thomas for sparse A
// [[Rcpp::export]]
List CG_sparse_KF_Thomas_multi(Eigen::MatrixXd A, Eigen::VectorXd diag_phi_L_T_inv, 
                               Eigen::VectorXd off_phi_L_T_inv, VectorXd b, const Eigen::VectorXd gg,  
                               const Eigen::VectorXd Q_K_2,  const Eigen::VectorXd sqrt_Q,
                               Eigen::MatrixXd P, Eigen::MatrixXd P_r,
                               VectorXd param, VectorXd delta_x, bool have_noise, String kernel_type,
                               const int L, int n, int D, float tol=1e-6, int maxIte = 1000){
  List return_list(3);
  int R = A.rows();
  int C = R;
  VectorXd x=VectorXd::Zero(C);
  
  VectorXd r = b;//-A*x;
  VectorXd p = r;
  double rs_old = (r.transpose() * r).value();
  double rs_new=1.0;
  double rs_ratio;
  VectorXd Bp = VectorXd::Zero(R);
  double alpha;
  int ite = 0;
  VectorXd resid = VectorXd::Zero(maxIte);
  
  while((ite < maxIte) && (rs_new > tol)){
    //Ap = A*p;
    //Bp = A.transpose()*p;
    Bp = t_sparse_A_times_x_multi(A, P_r, p, L, n, D);
    //Bp = L.transpose()*Bp;
    Bp = Fast_Thomas_for_CG(diag_phi_L_T_inv,off_phi_L_T_inv,Bp);
    //Bp = L_times_z(param, have_noise = false, delta_x, Bp, kernel_type = "exp");
    Bp = L_times_z_CG(gg,Q_K_2, sqrt_Q, Bp);
    //Bp = A*Bp;
    Bp = sparse_A_times_x_multi(A, P, Bp,L, n,D);
    Bp += exp(param[1])*p;
    
    alpha = rs_old / (p.transpose() * Bp).value();
    x += alpha*p;
    r -= alpha*Bp;
    rs_new = (r.transpose() * r).value();
    rs_ratio = rs_new / rs_old;
    p = r + rs_ratio * p;
    rs_old = rs_new;
    resid[ite] = rs_new; //mean((res2[[1]]-z)^2)
    ite++;
  }
  if (ite >= maxIte){
    Rcout << "cg did not converge." << endl;
  }
  
  return_list[0]=x;
  return_list[1]=ite;
  return_list[2]=resid.head(ite);
  return return_list;
}






////General CG
// [[Rcpp::export]]
List ConjGrad(Eigen::MatrixXd A, VectorXd b, float tol=1e-6, int maxIte = 1000){
  List return_list(3);
  int C = A.cols();
  int R = A.rows();
  VectorXd x=VectorXd::Zero(C);
  
  VectorXd r = b-A*x;
  VectorXd p = r;
  double rs_old = (r.transpose() * r).value();
  double rs_new=1.0;
  double rs_ratio;
  VectorXd Ap = VectorXd::Zero(R);
  double alpha;
  int ite = 0;
  VectorXd resid = VectorXd::Zero(maxIte);
  
  while((ite < maxIte) && (rs_new > tol)){
    Ap = A*p;
    alpha = rs_old / (p.transpose() * Ap).value();
    x += alpha*p;
    // if (ite%50 == 0){
    //   r = b-A*x;
    // }else{
    //   r -= alpha*Ap;
    // }
    r -= alpha*Ap;
    rs_new = (r.transpose() * r).value();
    rs_ratio = rs_new / rs_old;
    p = r + rs_ratio * p;
    rs_old = rs_new;
    resid[ite] = rs_new; //mean((res2[[1]]-z)^2)
    ite++;
  }
  if (ite >= maxIte){
    Rcout << "cg did not converge." << endl;
  }
  
  return_list[0]=x;
  return_list[1]=ite;
  return_list[2]=resid.head(ite);
  return return_list;
}

