rm(list = ls())

setwd("/Users/oscar/Desktop/community detection")
output_directory = "/Users/oscar/Desktop/community detection/simulations/"

library(glmgen)
#library(genlasso)
library(blockmodels)
source("utils.R")
source("utilsV2.R")
source("gfl2d_admm.R")
source("funs.R")
library(NMI)
library(cvxclustr)



n_grid = c(500,750,1000,1250,1500)
NMC =  100
K_grid =   c(12)


ind_K = 1

a  =  list()
a[[1]] = sort(runif(K_grid[1]))
a[[2]] =  sort(runif(K_grid[2]))
a[[3]] =  sort(runif(K_grid[3]))
lambda_grid = 10^seq(-1.2,.2,length= 20)#10^seq(-2,log(50)/log(10),length = 30)

MSE_1 =   array(0,c(length(K_grid),length(n_grid),NMC  ))
MSE_2 =   array(0,c(length(K_grid),length(n_grid),NMC  ))
MSE_3 =   array(0,c(length(K_grid),length(n_grid),NMC  ))


average_MSE_1 =   array(0,c(length(K_grid),length(n_grid)))
average_MSE_2 =   array(0,c(length(K_grid),length(n_grid) ))
average_MSE_3 =   array(0,c(length(K_grid),length(n_grid) ))


for(ind_K in 1:length(K_grid) )
{ 
  K =  K_grid[ind_K]
  b  =  matrix(0,K,K)
  diag(b) =  runif(K,0.2,1)
  aux = diag(b)
  
  b = matrix( min(aux)*runif(K*K) ,K,K)
  diag(b)  =  aux
  
  #matrix( runif(K*K) ,K,K)
  print("K")
  print(K)
  
  for(ind_n in 1:length(n_grid))
  {
    n =  n_grid[ind_n]
    print("n")
    print(n)
    
    for(iter  in 1:NMC)
    {
      print("iter")
      
      xi  =  sort(runif(n))
      
      P =  f0_1(xi,a[[ind_K]],b,K)
      P =  as.matrix(forceSymmetric(P))
      A =  1*(matrix(runif(n*n),n,n)<P )
      A =  as.matrix(forceSymmetric(A))
      diag(A) = 0
      diag(P) = 0
      
      ######################
      A_name  =  paste("A_ind1",'n',toString(n),'K',toString(K),"iter",iter,sep='')
      A_name = paste(output_directory,A_name,".txt",sep="")
      #write.table(A,A_name,col.names = FALSE,row.names = FALSE)
      
      P_name  =  paste("P_ind1",'n',toString(n),'K',toString(K),"iter",iter,sep='')
      P_name = paste(output_directory,P_name,".txt",sep="")
      #write.table(P,P_name,col.names = FALSE,row.names = FALSE)
      
      ###   S-FL
      degree  = colSums(A)  
      tour = order(degree)
      tour_inv = 1:n
      tour_inv[tour] = 1:n
      
      temp = TwoD_BernFL_pred(A[tour,tour],lambda_grid,n)
      theta_hat_tour = temp$best_sol
      #image(theta_hat_tour) 
      theta_hat = theta_hat_tour[tour_inv,tour_inv]
      diag(theta_hat) = 0
      MSE_1[ind_K,ind_n,iter] =  mean(((theta_hat + t(theta_hat))/2    -P )^2)
      
      ######################## L1 distance
      system.time({tour2  = nn_TSP(A, n, n,sample(1:n,1))})
      tour_inv2 = 1:n
      tour_inv2[tour2] = 1:n
      temp = TwoD_BernFL_pred(A[tour2,tour2],lambda_grid,n)
      theta_hat_tour2 = temp$best_sol
      #image(theta_hat_tour2) 
      theta_hat2 = theta_hat_tour2[tour_inv2,tour_inv2]
      diag(theta_hat2) = 0
      MSE_2[ind_K,ind_n,iter] =   mean(((theta_hat2 + t(theta_hat2))/2    -P )^2)
      
      ######################## inner product based
      inner =  matrix(0,n,n)
      indices =  sample(1:n,min(n,800))
      
      for(i in 1:n)
      {
        inner[i,] = apply(A[,indices],1,function(t){sum(t*A[i,indices])  })
      }
      
      system.time({tour3  = new_nn_TSP(A, inner, n, n,sample(1:n,1))})
      tour_inv3 = 1:n
      tour_inv3[tour3] = 1:n
      temp = TwoD_BernFL_pred(A[tour3,tour3],lambda_grid,n)
      theta_hat_tour3 = temp$best_sol
      #image(theta_hat_tour3) 
      theta_hat3 = theta_hat_tour3[tour_inv3,tour_inv3]
      diag(theta_hat3) = 0
      MSE_3[ind_K,ind_n,iter] =   mean(((theta_hat3 + t(theta_hat3))/2    -P )^2)
      
      
      
      
      
    }
    average_MSE_1[ind_K,ind_n]  = mean(MSE_1[ind_K,ind_n,]) 
    average_MSE_2[ind_K,ind_n]  = mean(MSE_2[ind_K,ind_n,]) 
    average_MSE_3[ind_K,ind_n]  = mean(MSE_3[ind_K,ind_n,]) 
    
    print("S-FL")
    print(average_MSE_1[ind_K,ind_n])
    print("L1-FL")
    print(average_MSE_2[ind_K,ind_n])
    print("NN-FL")
    print(average_MSE_3[ind_K,ind_n])
  }###############3
  
}##  close function




name1  =  paste("average_MSE_1_ind1",'n',toString(n),'K',toString(K),sep='')
name1 = paste(output_directory,name1,".txt",sep="")
write.table(average_MSE_1,name1,col.names = FALSE,row.names = FALSE)

name1  =  paste("MSE_1_ind1",'n',toString(n),'K',toString(K),sep='')
name1 = paste(output_directory,name1,".txt",sep="")
write.table(drop(MSE_1),name1,col.names = FALSE,row.names = FALSE)


name2  =  paste("average_MSE_2_ind1",'n',toString(n),'K',toString(K),sep='')
name2 = paste(output_directory,name2,".txt",sep="")
write.table(average_MSE_2,name2,col.names = FALSE,row.names = FALSE)

name2  =  paste("MSE_2_ind1",'n',toString(n),'K',toString(K),sep='')
name2 = paste(output_directory,name2,".txt",sep="")
write.table(drop(MSE_2),name2,col.names = FALSE,row.names = FALSE)

name3  =  paste("MSE_3_ind1",'n',toString(n),'K',toString(K),sep='')
name3 = paste(output_directory,name3,".txt",sep="")
write.table(drop(MSE_3),name3,col.names = FALSE,row.names = FALSE)

