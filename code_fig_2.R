grlee12 <- function(x)
{
  ##########################################################################
  #
  # GRAMACY & LEE (2012) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  
  term1 <- sin(10*pi*x) / (2*x)
  term2 <- (x-1)^4
  
  y <- term1 + term2
  return(y)
}


library(lhs)
library(RobustGaSP)
library(FastGaSP)

library(plot3D)

M=25
time_record=matrix(0,M,2)
colnames(time_record)=c('fastgasp','fullgasp')
  

p=1
record_pred_mean=as.list(1:M)
record_input_output=as.list(1:M)

for(i_M in 1:M){
  print(i_M)
  set.seed(i_M)
  n=200*i_M
  #input <- 0.5+2*seq(0,1,1/(n-1)) ##equally spaced
  #input=sort(0.5+2*maximinLHS(n=n, k=p))  # maximin lhd sample, sort it for now
  input=sort(0.5+2*runif(n))  # maximin lhd sample, sort it for now
  
  output=grlee12(input)+0.2*rnorm(n)
  n_testing=1000
  testing_input=0.5+2*(seq(0,1,1/(n_testing-1) )) ##need to sort, need to remove this 
  testing_output_noise_free=grlee12(testing_input)
  testing_output=testing_output_noise_free+0.2*rnorm(n_testing)
  #
  gamma=0.5
  eta=10^{-4} ## sigma_2_0/sigma_2

  time_record[i_M,1]=system.time(
    for(i in 1:5){
      fgasp.model=fgasp(input, output,kernel_type='matern_5_2') ##only implement expo or matern_5_2 for this version
      pred_mean_fastgasp=predict(param=c(log(1/gamma),log(eta)),object=fgasp.model, testing_input=testing_input, var_data=F) 
    }
  )[3]/5  ##repeat 5 time and compute average

  time_record[i_M,2]=system.time(
    for(i in 1:5){
      R0=abs(outer(input,input,'-'))
      R=matern_5_2_funct(R0,beta_i=1/gamma)
      R_tilde=R+eta*diag(n)
      r0=abs(outer(input,(testing_input),'-'))
      r=matern_5_2_funct(r0,beta_i=1/gamma)
      chol_R_tilde=t(chol(R_tilde))
      pred_mean=t(r)%*%(backsolve(t(chol_R_tilde), forwardsolve(chol_R_tilde,output)))
    }
  )[3]/5   ##repeat 5 time and compute average
  
  record_pred_mean[[i_M]]=matrix(NA,n_testing,4)
  record_pred_mean[[i_M]][,1]=pred_mean_fastgasp@mean
  record_pred_mean[[i_M]][,2]=pred_mean
  record_pred_mean[[i_M]][,3]=testing_output_noise_free
  record_pred_mean[[i_M]][,4]=testing_output
  
  record_input_output[[i_M]]=matrix(NA,n,2)
  record_input_output[[i_M]][,1]=input
  record_input_output[[i_M]][,2]=output
  
}
# image2D(matrix(pred_model$mean,100,100),testing_input[1:100,1],testing_input[seq(1:100)*100,2],
#         zlim=zlim,main='Prediction, n=24',xlab=expression(x[1]),ylab=expression(x[2]),
#         cex.axis=1.45, cex.lab=1.45,cex.main=1.45,cex.clab=1.45)
# points(input[,1],input[,2],cex=1,lwd=2)

#pred_mean2=t(r)%*%(solve(R_tilde)%*%output) ##may code it in C++

pdf(file='time_comparison_KF.pdf',height=4.2,width=4,2)
plot(200*(1:M),time_record[,2],col='red',type='p',pch=17,xlab='N',ylab='time (s)',
     main='Computational cost',mgp=c(2.5,1,0))
lines(200*(1:M),time_record[,1],col='blue',type='p',pch=19)
legend('topleft', legend=c('GP by direct computation','GP by KF'), pch=c(17,19),col=c('red','blue'))
dev.off()

pdf(file='pred_comparison_KF.pdf',height=4.2,width=4.2)

i_M=5
plot(testing_input,record_pred_mean[[i_M]][,3],col='black',
     type='l',lty=1,xlab='x',ylab='y',ylim=c(-1.1,5.3),
     main='Prediction, N=1000',mgp=c(2.5,1,0))
lines(record_input_output[[i_M]][,1], record_input_output[[i_M]][,2],col='black',type='p',pch=20,cex=.15)

lines(testing_input,record_pred_mean[[i_M]][,1],col='red',type='l',lty=2,lwd=2)
lines(testing_input,record_pred_mean[[i_M]][,2],col='blue',type='l',lty=3,lwd=2)
legend('topleft', legend=c('Observations','Truth', 'GP by direct computation','GP by KF'),
       pch=c(20,NA,NA,NA),lty=c(NA,1,2,3),col=c('black','black','red','blue'),lwd=c(NA,1,2,2))
dev.off()

save(testing_input,record_pred_mean,time_record,M,record_input_output,file='fig_2_comp.RData')


#pred_mean_fastgasp@mean-(testing_output_noise_free)
#sqrt(mean( (pred_mean-pred_mean_fastgasp@mean)^2))
sqrt(mean( (record_pred_mean[[i_M]][,1]-record_pred_mean[[i_M]][,2])^2))


# 
# 
# 
# time_record_fastgasp=system.time(
#   for(ii in 1:1){
#     
#     fgasp.model=fgasp(x, y,kernel_type='exp') ##only implement expo or matern_5_2 for this version
#     
#     ###first term
#     log_det_S2=Get_log_det_S2( c(log(1/gamma),log(sigma_2_0/sigma_2)),fgasp.model@have_noise,fgasp.model@delta_x,
#                                fgasp.model@output,fgasp.model@kernel_type)
#     
#     log_lik_fastgasp=-n/2*log(2*pi*sigma_2)-log_det_S2[[1]]/2-log_det_S2[[2]]/(2*sigma_2)
#   }
# )
# 
# 
# 
# ##robust gasp
# model=rgasp(design=input,response=output)
# pred_model=predict(model,testing_input)
# pred_model$mean -testing_output
# 
# 
# 
# library(RobustGaSP)
# library(plot3D)
# set.seed(1)
# n=12
# p=2
# LB=c(-5,0)
# UB=c(10,15)
# range=UB-LB
# 
# input <- maximinLHS(n=n, k=p)  # maximin lhd sample
# 
# 
#   for(i in 1:p){
#     input[,i]=LB[i]+range[i]*input[,i]
#   }
# 
# output=matrix(0,n,1)
# for(i in 1:n){
#   output[i]=branin(input[i,])
# }
# 
# model=rgasp(design=input,response=output)
# 
# testing_input=matrix(NA,100^2,2)
# testing_input[,1]=as.vector(matrix(seq(0,1,1/99),100,100))
# testing_input[,2]=as.vector(t(matrix(seq(0,1,1/99),100,100)))
# 
# for(i in 1:p){
#   testing_input[,i]=LB[i]+range[i]*testing_input[,i]
# }
# 
# testing_output=matrix(0,100^2,1)
# for(i in 1:(100^2) ){
#   testing_output[i]=branin(testing_input[i,])
# }
# 
# 
# pred_model=predict(model,testing_input)
# 
# sqrt(mean((pred_model$mean-testing_output)^2))
# sd(testing_output)
# 
# 
# testing_output_mat=matrix(testing_output,100,100)
# 
# input_12=input
# pred_model_mean_12=pred_model$mean
# 
# ###24
# set.seed(1)
# n=24
# p=2
# LB=c(-5,0)
# UB=c(10,15)
# range=UB-LB
# 
# input <- maximinLHS(n=n, k=p)  # maximin lhd sample
# 
# 
# for(i in 1:p){
#   input[,i]=LB[i]+range[i]*input[,i]
# }
# 
# output=matrix(0,n,1)
# for(i in 1:n){
#   output[i]=branin(input[i,])
# }
# 
# model=rgasp(design=input,response=output)
# 
# testing_input=matrix(NA,100^2,2)
# testing_input[,1]=as.vector(matrix(seq(0,1,1/99),100,100))
# testing_input[,2]=as.vector(t(matrix(seq(0,1,1/99),100,100)))
# 
# for(i in 1:p){
#   testing_input[,i]=LB[i]+range[i]*testing_input[,i]
# }
# 
# testing_output=matrix(0,100^2,1)
# for(i in 1:(100^2) ){
#   testing_output[i]=branin(testing_input[i,])
# }
# 
# 
# pred_model=predict(model,testing_input)
# 
# sqrt(mean((pred_model$mean-testing_output)^2))
# sd(testing_output)
# 
# 
# testing_output_mat=matrix(testing_output,100,100)
# 
# 
# zlim=c(-20,315)
# 
# 
# pdf('Truth_Branin.pdf',height=5,width=5)
# image2D(matrix(testing_output,100,100),testing_input[1:100,1],
#         testing_input[seq(1:100)*100,2],zlim=zlim,main='Truth',
#         xlab=expression(x[1]),ylab=expression(x[2]),colkey = list(plot = F, side = 4),
#         cex.axis=1.45, cex.lab=1.45,cex.main=1.45)
# dev.off()
# pdf('Prediction_Branin_N_12.pdf',height=5,width=5)
# 
# image2D(matrix(pred_model_mean_12,100,100),testing_input[1:100,1],testing_input[seq(1:100)*100,2],
#         zlim=zlim,main='Prediction, N=12',xlab=expression(x[1]),ylab=expression(x[2]),colkey = list(plot = FALSE, side = 4),
#         cex.axis=1.45, cex.lab=1.45,cex.main=1.45,cex.clab=1.45)
# points(input_12[,1],input_12[,2],cex=1,lwd=2)
# dev.off()
# pdf('Prediction_Branin_N_24.pdf',height=5,width=5)
# 
# image2D(matrix(pred_model$mean,100,100),testing_input[1:100,1],testing_input[seq(1:100)*100,2],
#         zlim=zlim,main='Prediction, N=24',xlab=expression(x[1]),ylab=expression(x[2]),colkey = list(plot = T, side = 4),
#         cex.axis=1.45, cex.lab=1.45,cex.main=1.45,cex.clab=1.45)
# points(input[,1],input[,2],cex=1,lwd=2)
# dev.off()

# pdf('Truth_Branin.pdf',height=5,width=5)
# image2D(matrix(testing_output,100,100),testing_input[1:100,1],
#         testing_input[seq(1:100)*100,2],zlim=zlim,main='Truth',
#         xlab=expression(x[1]),ylab=expression(x[2]),colkey = list(plot = FALSE, side = 3),
#         cex.axis=1.45, cex.lab=1.45,cex.main=1.45)
# dev.off()
# pdf('Prediction_Branin_n_24.pdf',height=5,width=5)
# image2D(matrix(pred_model$mean,100,100),testing_input[1:100,1],testing_input[seq(1:100)*100,2],
#         zlim=zlim,main='Prediction, n=24',xlab=expression(x[1]),ylab=expression(x[2]),
#         cex.axis=1.45, cex.lab=1.45,cex.main=1.45,cex.clab=1.45)
# points(input[,1],input[,2],cex=1,lwd=2)
# dev.off()
# 
