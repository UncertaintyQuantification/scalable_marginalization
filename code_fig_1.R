
library(lhs)
library(RobustGaSP)
library(plot3D)

branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi))
{
  ##########################################################################
  #
  # BRANIN FUNCTION
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
  #
  # INPUTS:
  #
  # xx = c(x1, x2)
  # a = constant (optional), with default value 1
  # b = constant (optional), with default value 5.1/(4*pi^2)
  # c = constant (optional), with default value 5/pi
  # r = constant (optional), with default value 6
  # s = constant (optional), with default value 10
  # t = constant (optional), with default value 1/(8*pi)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  return(y)
}


set.seed(1)
n=12
p=2
LB=c(-5,0)
UB=c(10,15)
range=UB-LB

input <- maximinLHS(n=n, k=p)  # maximin lhd sample


  for(i in 1:p){
    input[,i]=LB[i]+range[i]*input[,i]
  }

output=matrix(0,n,1)
for(i in 1:n){
  output[i]=branin(input[i,])
}

model=rgasp(design=input,response=output)

testing_input=matrix(NA,100^2,2)
testing_input[,1]=as.vector(matrix(seq(0,1,1/99),100,100))
testing_input[,2]=as.vector(t(matrix(seq(0,1,1/99),100,100)))

for(i in 1:p){
  testing_input[,i]=LB[i]+range[i]*testing_input[,i]
}

testing_output=matrix(0,100^2,1)
for(i in 1:(100^2) ){
  testing_output[i]=branin(testing_input[i,])
}


pred_model=predict(model,testing_input)

sqrt(mean((pred_model$mean-testing_output)^2))
sd(testing_output)


testing_output_mat=matrix(testing_output,100,100)

input_12=input
pred_model_mean_12=pred_model$mean

###24
set.seed(1)
n=24
p=2
LB=c(-5,0)
UB=c(10,15)
range=UB-LB

input <- maximinLHS(n=n, k=p)  # maximin lhd sample


for(i in 1:p){
  input[,i]=LB[i]+range[i]*input[,i]
}

output=matrix(0,n,1)
for(i in 1:n){
  output[i]=branin(input[i,])
}

model=rgasp(design=input,response=output)

testing_input=matrix(NA,100^2,2)
testing_input[,1]=as.vector(matrix(seq(0,1,1/99),100,100))
testing_input[,2]=as.vector(t(matrix(seq(0,1,1/99),100,100)))

for(i in 1:p){
  testing_input[,i]=LB[i]+range[i]*testing_input[,i]
}

testing_output=matrix(0,100^2,1)
for(i in 1:(100^2) ){
  testing_output[i]=branin(testing_input[i,])
}


pred_model=predict(model,testing_input)

sqrt(mean((pred_model$mean-testing_output)^2))
sd(testing_output)


testing_output_mat=matrix(testing_output,100,100)


zlim=c(-20,315)


pdf('Truth_Branin.pdf',height=5,width=5)
image2D(matrix(testing_output,100,100),testing_input[1:100,1],
        testing_input[seq(1:100)*100,2],zlim=zlim,main='Truth',
        xlab=expression(x[1]),ylab=expression(x[2]),colkey = list(plot = F, side = 4),
        cex.axis=1.45, cex.lab=1.45,cex.main=1.45)
dev.off()
# pdf('Prediction_Branin_N_12.pdf',height=5,width=5)
# 
# image2D(matrix(pred_model_mean_12,100,100),testing_input[1:100,1],testing_input[seq(1:100)*100,2],
#         zlim=zlim,main='Prediction, N=12',xlab=expression(x[1]),ylab=expression(x[2]),colkey = list(plot = FALSE, side = 4),
#         cex.axis=1.45, cex.lab=1.45,cex.main=1.45,cex.clab=1.45)
# points(input_12[,1],input_12[,2],cex=1,lwd=2)
# dev.off()
pdf('Prediction_Branin_N_24.pdf',height=5,width=5)

image2D(matrix(pred_model$mean,100,100),testing_input[1:100,1],testing_input[seq(1:100)*100,2],
        zlim=zlim,main='Prediction, N=24',xlab=expression(x[1]),ylab=expression(x[2]),colkey = list(plot = T, side = 4),
        cex.axis=1.45, cex.lab=1.45,cex.main=1.45,cex.clab=1.45)
points(input[,1],input[,2],cex=1,lwd=2)
dev.off()

