
####################################################################################
sotft_thresholding = function(y,lambda)
{
  return( sign(y)*pmax(abs(y)-lambda,0))
}
####################################################################################

################################################
library(Matrix)
g = function(x)
{
  x1 = x[,1]
  x2 = x[,2]
  #   
  #   a=  x2*x1 + 3*x2^2 
  #   # a = abs(3*x1 - 2*x2) ## good
  # #   
  #   for(i in 1:length(a))
  #   { 
  #     if(x1[i]<1/2 && x2[i]<1/2) 
  #     #{a[i] = -1.2}
  #     {a[i] =-2 - x1[i] - 2*sin(2*pi*(x2[i]))} 
  #   }
  ################################################
  #    a =  abs(2*x1-x2)^.2
  # 
  #     for(i in 1:length(a))
  #     { 
  #       if( sqrt((x1[i] -.5)^2 + (x2[i]-.5)^2)< .2) 
  #       #{a[i] = -1.2}
  #       {a[i] = abs(2*x1[i]^2+x1[i]*x2[i] - x2[i])^.4    + 1} 
  #     }
  #    )
  
  ################################################
  a = rep(0,length(x1))
  
  for(i in 1:length(a))
  { 
    if(  1/2*(x1[i] + x2[i]-1/2)^2 + 8/10*(x2[i]-x1[i])^2 < 1/32     ) 
    {a[i] = -1}
    if(  1/2*(x1[i] -x2[i])^2 + 10/8*(x2[i]+ x1[i] - 14/8)^2 < 1/32     ) 
    {a[i] = 1}
    #if( .2*(x1[i]-1/4)^2 + .8*(x2[i]-1/4)^2 < .3 ) 
    #{a[i] =  1}
    #if(x1[i] > 1/2 && x2[i]<1/2) 
    #{a[i] =  1}
    #if(x1[i] > 1/2 && x2[i]>1/2) 
    #{a[i] = -1}
  }
  
  return( a )
  
  a = rep(0,dim(x)[1])
  #   for(i in 1:dim(x)[1])
  #   {
  #     x1 = x[i,1]
  #     x2 = x[i,2]
  #     if(sin(2*pi*x1 )/10 + 1/2  < x2)
  #     {
  #       a[i] = 1.5 
  #     }
  #   }
  #   return(a)
  #   x1 = x[,1]
  #   x2 = x[,2]
  #   
  #   a = rep(0,dim(x)[1])
  #   
  #   for(i in 1:dim(x)[1])
  #   {
  #     if ((x1[i]-.5)^2 +  (x2[i] - 0.5)^2  < .2^2)
  #     {
  #       a[i] = exp(x1[i]) + exp(x2[i])  
  #     }
  #     
  #   }
  #   #3 b*x2 + c - A*sin(w*x1+1) + 2*A*sin(w*x2+1) + (x1+2)^2
  #   #return(2*sin(w*x1+1)  +   2*cos(w*x2+1)+ x2*x2 - x1^3)
  #   # return( sin(w*x1+1)  + abs(x2)+  2*cos(w*x2+1)+ x2*x1 - x1^3 + 3.4*x2^2+ sqrt(x1^2+abs(x2)))
  #   return(  x2*x1 - x1^3 + 3.4*x2^2)
}

f =  function(x)
{
  
  if(length(x)==2)
  {
    out = 0
    temp = (x[1]-.5)^2 +  (x[2]-.5)^2
    if(temp < .15^2)
    {
      out =  1
    }
    return(out)
  }
  
  out = rep(0,dim(x)[1])
  for(i in 1:dim(x)[1])
  {
    temp = (x[i,1]-.5)^2 +  (x[i,2]-.5)^2
    if(temp < .25)
    {
      out[i] =  1
    }
  }
  return(out)
}
################################################
################################################
################################################

tsp_npr = function(y_tour,w,prop= .25,lambda_grid = 2^seq(-1,10,length = 200))
{
  n = length(y_tour)
  #prop = .25
  indices = sample(2:(n-1),floor(n*prop))
  y_training = y_tour[-indices]
  y_test  =  y_tour[indices]
  
  loc = seq(0,1,length = n)
  MSE = rep(0,length(lambda_grid))
  for(i in 1:length(lambda_grid))
  {
    temp = FlsaDp(y_training, lambda2  = lambda_grid[i], obsWt = w)
    #temp = FlsaDp(y_tour, lambda2  = lambda_grid[i])#, obsWt = w)
    prediction =  approx(x= loc[-indices] ,y= temp, xout = loc[indices]) 
    MSE[i] = mean((prediction$y  - y_test)^2)
    #MSE[i] = mean((f0[tour] - temp)^2)
  }
  min(MSE)
  i =which.min(MSE)
  f_hat = FlsaDp(y_tour, lambda2 = lambda_grid[i], obsWt = w)
  # plot(f_hat)
  return(f_hat) 
}


################################################
################################################
################################################

tsp_npr_degrees_of_freedom = function(y_tour,w,prop= .25,lambda_grid = 2^seq(-1,10,length = 200))
{
  n = length(y_tour)
  
  sigma = sqrt(.5/(n-1)*sum((diff(y_tour))^2))
  #prop = .25
  #indices = sample(2:(n-1),floor(n*prop))
  #y_training = y_tour[-indices]
  #y_test  =  y_tour[indices]
  
  loc = seq(0,1,length = n)
  aic = rep(0,length(lambda_grid))
  temp = matrix(0,length(lambda_grid),n)
  for(i in 1:length(lambda_grid))
  {
    temp[i,] = FlsaDp(y_tour, lambda2  = lambda_grid[i], obsWt = w)
    aic[i] = sum((temp[i,] - y_tour)^2) + 2*sigma*sigma* length(which(abs(diff(temp[i,]))>10^-4)) 
    #temp = FlsaDp(y_tour, lambda2  = lambda_grid[i])#, obsWt = w)
    #prediction =  approx(x= loc[-indices] ,y= temp, xout = loc[indices]) 
    #MSE[i] = mean((prediction$y  - y_test)^2)
    #MSE[i] = mean((f0[tour] - temp)^2)
  }
  #min(MSE)
  i =which.min(aic)
  f_hat = temp[i,]
  #FlsaDp(y_tour, lambda2 = lambda_grid[i], obsWt = w)
  # plot(f_hat)
  return(f_hat) 
}

#############################################################################
#############################################################################
########################################################

fused_lasso_admm = function(y,Delta,lambda)
{
  n = length(y)
  rho = 1
  z = rep(0,dim(Delta)[1])
  u = rep(0,dim(Delta)[1])
  u_prev =  u
  z_prev = z
  beta =  rep(0,n)
  beta_prev = beta
  r =  rep(0,length(z))
  s =  rep(0, n)
  
  max_iter = 400
  
  I = sparseMatrix(i = 1:n,j = 1:n,x = rep(1,n),dims = c(n,n))
  
  tDelta = t(Delta)
  tDeltaDelta = t(Delta)%*% Delta
  
  for(iter in 1:max_iter)
  {
    # print(iter)
    ### update beta
    beta =  solve( I + rho*tDeltaDelta , y + drop(tDelta %*% (z-u)),sparse = TRUE ) 
    ### update z
    z  = sotft_thresholding(  drop(Delta%*%beta) + u,lambda/rho)
    ### update u
    r = drop(Delta %*% beta) -  z
    u =  u + r
    s = rho*drop(tDelta %*%(z-z_prev) )
    
    ### update rho
    norm_s =   sum(s*s)
    norm_r =   sum(r*r) 
    
    if(norm_r > 10*norm_s)
    {
      rho =  2*rho
      u = u/2
    }
    if(norm_s > 10*norm_r)
    {
      rho =  rho/2
      u =  2*u
    }
    u_prev  =u
    s_prev = s
    z_prev = z
    beta_prev = beta
    
    ################################################
    ###  stopping criteria
    
    if(mean(r^2)< 10^-6  && mean(s^2)< 10^-6)
    {
      break;
    }
  }
  
  return(beta)
}

################################################

library(Rcpp)
cppFunction('
            
            NumericVector nn_tour(NumericVector f, int n,int init)
            {
            NumericVector state(n,0.0);
            state[init-1] = 1;
            NumericVector order(n,0.0); 
            order[0] = init;
            
            int next;
            float dist;
            float max_dist = pow(10.0,15);
            float temp1;   
            
            for(int i = 1; i < n; i++)
            {
            next = 1;
            dist = max_dist;
            int k = order[i-1];
            
            for(int j = 1; j<(n+1);j++)
            {
            if( state[j-1] < .5)
            {
            
            temp1 = f[k-1 ]  - f[j-1];
            
            if(temp1 < 0)
            {
            temp1 =  - temp1;            
            }
            
            if( temp1 < dist)
            {
            dist = temp1;
            next = j;
            }
            }
            }
            order[i] = next;  
            state[next-1] = 1; 
            }
            return order;  
            }')



# aa = nn_tour(svm.pred,n,tour_temp[1])
# sum(abs(diff(svm.pred[aa])))
# sum(abs(diff(svm.pred[tour_temp])))
# aa[1:10]
# tour_temp[1:10]
# 
# abs(svm.pred[aa[1]] -svm.pred[aa[2]] )
# abs(svm.pred[tour_temp[1]] -svm.pred[tour_temp[2]] )

########################################################################
######################################################################
####################################################################



g_new = function(x,indice)
{
  x1 = x[,1]
  x2 = x[,2]
  
  if(indice==1)
  {
    a = rep(0,length(x1))
    
    for(i in 1:length(a))
    {
      if( (x1[i]-6/8)^2 + (x2[i]-6/8)^2 < 1/32 )
      {
        a[i] =1
      }
      if( (x1[i]-2/8)^2 + (x2[i]-3/8)^2 < 1/16 )
      {
        a[i] =-1
      }
      #  if(  1/2*(x1[i] + x2[i]-1/2)^2 + 8/10*(x2[i]-x1[i])^2 < 1/32     ) 
      #  {a[i] = -1}
      #  if(  1/2*(x1[i] -x2[i])^2 + 10/8*(x2[i]+ x1[i] - 14/8)^2 < 1/32     ) 
      #  {a[i] = 1}
    }
    return( a )
  }
  
  a = rep(0,length(x1))
  if(indice == 2)
  {
    #       a = abs(3*x1 - 2*x2^2) ## good
    #       for(i in 1:length(a))
    #       { 
    #         if(x1[i]<1/2 && x2[i]<1/2) 
    #         #{a[i] = -1.2}
    #         {a[i] =-1}
    #       }
    ind = which(x1 < cos(x2*4*pi) ) 
    a[ind ] = 1
    
    #     ind =   which(x2 < cos(x1*4*pi)) 
    #     a[ind ] = -1    
    
    return(a)
  }
  if(indice ==3)
  {
    
    return( x1^2 + x2^2 - 3*x1*x2 +  .5*x1  )
  }
  
  
  a =  rep(0,length(x1))
  if(indice==4)
  {
    
    ind = which(x1 < sin(x2*2*pi)) 
    #      ind1 =  which( (x1-.5)^2 +  (x2-.5)^2 > 1/70 )
    #      ind2 =  which( (x1-.5)^2 +  (x2-.5)^2 < 1/6 )
    #      ind = intersect(ind1,ind2)
    #      a[ind ] = -1
    #      
    #      ind = which((x1-.5)^2 +  (x2-.5)^2 < 1/64)
    #      
    a[ind] = 1
    
    return(a)
  }
  
}

#####################################

astp_func =  function(x,n)
{
  D =  matrix(0,n,n)
  for(i in 1:n)
  {
    temp = sqrt((x[,1] - x[i,1])^2 + (x[,2] - x[i,2])^2)
    D[i,] =   temp
    D[,i] =   temp  
  }
  atsp <- ATSP(D)
  return(atsp)
}

#############################################

TSP1_est = function(x,y,N,n,astp,df = 1)
{
  
  f_hat2  = rep(0,n)
  f_hat22  = rep(0,n)
  
  tour_temp_inv = rep(0,n)
  tour_temp_inv2 = rep(0,n)
  
  lambda_grid = 2^seq(-3,9,length = 250)
  
  save_f_hat = matrix(0,N,n)
  
  for(i in 1:N)
  {
    #as.numeric(solve_TSP(atsp) )#
    tour_temp = as.numeric(solve_TSP(atsp,method ="cheapest_insertion"))
    tour_temp_inv[tour_temp] = 1:n 
    
    w = rep(1,n)
    if(df==1)
    {
      aux = tsp_npr_degrees_of_freedom(y[tour_temp],w,prop = .25,lambda_grid )
    }
    else
    {
      aux = tsp_npr(y[tour_temp],w,prop = .25,lambda_grid ) 
    }
    aux = aux[tour_temp_inv]
    
    save_f_hat[i,] = aux
    
  }###
  
  if(N == 1)
  {
    return(as.vector(save_f_hat))
  }
  
  return(list(est = apply(save_f_hat[1:N,], 2,mean),raw = save_f_hat) )
}
###########################################################
####################################################

find_closest_point =  function(x,x_test)
{
  n = dim(x)[1]
  n_test = dim(x_test)[1]
  closet_ind = matrix(0,n_test,2)
  
  for( i in 1:n_test)
  {
    d = apply(x,1,function(t){ sum((t-x_test[i,])^2) })
    closet_ind[i,1] = which.min(d)
    d[closet_ind[i,1]] = 10^5
    closet_ind[i,2] = which.min(d)
  }
  return(closet_ind)
}
######################################################
######################################################
######################################################
######################################################

L2_distance = function(x1,x2)
{
  if(is.matrix(x2))
  {
    n2 =  dim(x2)[1]
    dist = rep(0,n2)
    
    for(i in 1:n2)
    {
      dist[i] =  sqrt(sum( (x1 -x2[i,] )^2 ))
    }
  }
  if(is.matrix(x2)==0)
  {
    dist =  sqrt(sum( (x1 -x2)^2 ))
  }
  
  return(dist)
}

#.5*sqrt(100)


######################################################
######################################################
######################################################


f0_func = function(x,ind)
{
  p = dim(x)[2]
  n = dim(x)[1]
  
  eval = rep(0,n)
  
  if(ind ==0)
  {
    mu1 =  rep(1/2,p)
    for( i in 1:n)
    {
      eval[i] =  L2_distance(x[i,],mu1)^2
    }### loop for i
    
  } ###  case 1  
  
  
  if(ind ==1)
  {
    mu1 =  rep(1/2,p)
    for( i in 1:n)
    {
      eval[i] =  L2_distance(x[i,],mu1)^.75
    }### loop for i
    
  } ###  case 1  
  
  
  if(ind ==2)
  {
    mu1 =  rep(1/4,p)
    mu2 =  rep(3/4,p)
    for( i in 1:n)
    {
      temp1 = L2_distance(x[i,],mu1)^.75
      temp2 = L2_distance(x[i,],mu2)^.75
      
      if(temp1 < temp2)
      {
        eval[i] =  1 
      }
      else{
        eval[i] =  -1
      }
    }### loop for i
  } ###  case 1
  
  
  if(ind ==3)
  {
    mu1 =  rep(1/4,p)
    mu2 =  rep(3/4+.1,p)
    for( i in 1:n)
    {
      temp1 = L2_distance(x[i,],mu1)^.75
      temp2 = L2_distance(x[i,],mu2)^.75
      
      if(temp1 < temp2)
      {
        eval[i] =  1 
      }
      else{
        eval[i] =  -1
      }
    }### loop for i
  } ###  case 1
  
  if(ind ==4)
  {
    mu1 =  rep(1/4+.1,p)
    mu2 =  rep(3/4-.05,p)
    mu3 =  rep(1/2,p)
    for( i in 1:n)
    {
      temp1 = L2_distance(x[i,],mu1)
      temp2 = L2_distance(x[i,],mu2)
      temp3 = L2_distance(x[i,],mu3)
      
      ind   = which.min(c(temp1,temp2,temp3))
      eval[i] =  ind - 2 
      
    }### loop for i
  } ###  case 1
  
  if(ind ==5)
  {
    mu1 =  rep(1/4-.1,p)
    mu2 =  rep(3/4,p)
    for( i in 1:n)
    {
      temp1 = L2_distance(x[i,],mu1)
      temp2 = L2_distance(x[i,],mu2)
      
      if(temp1 < temp2)
      {
        eval[i] = -1
      }
      else{
        eval[i]  = temp2^2  
      }
      
    }### loop for i
  } ###  case 1
  
  
  
  if(ind ==6)
  {
    mu2 = rep(0.5,p)
    mu3 = rep(0.5,p)
    mu1 = rep(0.5,p)
    
    # mu1[1:floor(p/3)] = .5 
    mu2[ 1:floor(p/2)]  = .6
    mu3[(1+floor(p/2)):p ] =  .6
    
    for( i in 1:n)
    {
      temp1 = L2_distance(x[i,],mu1)
      temp2 = L2_distance(x[i,],mu2)
      temp3 = L2_distance(x[i,],mu3)
      
      if(temp1 < min(temp2,temp3))
      {
        eval[i] = 0
      }
      if(temp2 < min(temp1,temp3))
      {
        eval[i]  = -2 
      }
      if(temp3 < min(temp1,temp2))
      {
        eval[i]  = 2 
      }
      
    }### loop for i
  } ###  case 1
  
  
  
  if(ind ==7)
  {
    mu1 = rep(0.5,p)
    mu2 = rep(0.5,p)
    mu3 = rep(0.5,p)
    mu4 = rep(0.5,p)
    
    mu1[1:floor(p/2)] = .4 
    mu2[ 1:floor(p/2)]  = .6
    mu3[(1+floor(p/2)):p ] =  .4
    mu4[(1+floor(p/2)):p ] =  .6
    
    for( i in 1:n)
    {
      temp1 = L2_distance(x[i,],mu1)
      temp2 = L2_distance(x[i,],mu2)
      temp3 = L2_distance(x[i,],mu3)
      temp4 = L2_distance(x[i,],mu4)
      
      if(temp1 < min(c(temp2,temp3,temp4)))
      {
        eval[i] = -2
      }
      if(temp2 < min(c(temp1,temp3,temp4)))
      {
        eval[i]  = 2 
      }
      if(temp3 < min(c(temp2,temp1,temp4)))
      {
        eval[i]  = 1 
      }
      if(temp4 < min(c(temp2,temp1,temp3)))
      {
        eval[i]  = -1 
      }
      
    }### loop for i
  } ###  case 1
  
  
  return(eval)
}

# #####
# L2_distance(mu1,x)
# 
#  R = 1/2
#  exp(log(pi)*p/2 + p*log(R) - lgamma(p/2 +1))
# # 
#  R = 1/2
#  exp(log(pi)*2/2 + 2*log(R) - lgamma(2/2 +1))
# pi * R^2 
find_ordering_NN = function(X1)
{
  
  tour_length = 0  
  
  n1 =  dim(X1)[1];
  n_max = 2000;
  i1 =  min(n1, 1+ floor(runif(1)*n1));
  
  factor =  rep(1,n1);
  factor[i1] = 10^10;
  order =  rep(0,n1);
  order[1] = i1;
  
  for(j in 2:n1)
  {
    indices = which(factor <10^10);
    n_new = min(n_max,length(indices));
    D = rep(0,n_new);
    subset_indices = sample(1:length(indices),n_new);
    
    for(k in 1:n_new)
    {
      D[k] =  sqrt(sum( (X1[i1,]-X1[indices[subset_indices[k]],])^2 ))
    }#######
    ind = which( (D+1)*factor[indices[subset_indices]] == min((D+1)*factor[indices[subset_indices]])  )
    tour_length =   tour_length + sqrt(D[ind])
    ind = indices[subset_indices[ind[1]]];
    order[j] = ind;
    i1 = ind;
    factor[ind] = 10^10;
    
  }#####
  
  return(list(order = order, cost =  tour_length))
}##########



#########################################################
#####################################################
#####################################################3
##############################################
##############################################

TSP_FL = function(x,y,N,n,df = 1, Num_lambda = 370)
{
  p = dim(x)[2]
  f_hat2  = rep(0,n)
  f_hat22  = rep(0,n)
  
  tour_temp_inv = rep(0,n)
  tour_temp_inv2 = rep(0,n)
  
  lambda_grid =  2^seq(-3,9,length =  Num_lambda)
  
  save_f_hat = matrix(0,N,n)
  
  cost = rep(0,N)
  
  for(i in 1:N)
  {
    #as.numeric(solve_TSP(atsp) )#
    temp = nn_TSP(x, p, n,sample(1:n,1))
    #find_ordering_NN(x)
    tour_temp = temp#$order
    #cost[i] = temp$cost
    #tour_temp = as.numeric(solve_TSP(atsp,method ="nn"))
    tour_temp_inv[tour_temp] = 1:n 
    
    w = rep(0,n)
    # w = exp(-.1*sum(diff(x[tour_temp,])^2))
    
    if(df==1)
    {
      aux = tsp_npr_degrees_of_freedom(y[tour_temp],w = NULL,prop = .25,lambda_grid )
    }
    else
    {
      aux = tsp_npr(y[tour_temp],w,prop = .25,lambda_grid ) 
    }
    aux = aux[tour_temp_inv]
    
    save_f_hat[i,] = aux
    
  }###
  
  if(N == 1)
  {
    return(list(est = as.vector(save_f_hat)))
  }
  
  # apply(save_f_hat[1:N,], 2,mean)
  return(list(est  = apply(save_f_hat[1:N,], 2,mean),raw = save_f_hat) )
}


#################################################################
#################################################################
#

Knn_TSP_FL  = function(x,y,k)
{
  n =  length(y)
  
  y_hat =  rep(0,n)
  
  for(i in 1:n)
  {
    temp = L2_distance(x[i,],x)
    temp  = order(temp,decreasing = FALSE)
    indices = temp[2:k]
    y_hat[i] = mean(y[indices])
  }
  # mean((y_hat - f0)^2)
  tour =  order(y_hat)
  
  lambda_grid = 2^seq(-2,10,length = 30)
  
  w = rep(0,n)
  aux = tsp_npr_degrees_of_freedom(y[tour],w,prop = .25,lambda_grid )
  tour_inv =  1:n
  tour_inv[tour] = 1:n
  theta_hat =  aux[tour_inv] 
  
  #    mean((theta_hat - f0)^2)
  
  return(theta_hat)   
}
# 
#  i
# temp[1:k]
#  L2_distance(x[i,],mu1)< log(p)/2.2
# f0[temp[1:k]]
#  L2_distance(x[temp[5],],mu1)< log(p)/2.2 
#  L2_distance(x[i,],x[temp[1:k],])
# # exp(log(.30)/p)/2
# max(abs((x[i,]-mu1)))
# max(abs((x[temp[2],]-mu1)))
# f0[temp[1]]
# f0[temp[2]]


#################################################################
#################################################################
##################################################################
#################################################################
##################################################################
#################################################################
#
cppFunction('
            NumericVector nn_TSP(NumericMatrix f,int p, int n,int init)
            {
            NumericVector state(n,0.0);
            state[init-1] = 1;
            NumericVector order(n,0.0); 
            order[0] = init;
            
            int next;
            float dist;
            float max_dist = pow(10.0,15);
            float temp1;   
            
            for(int i = 1; i < n; i++)
            {
            next = 1;
            dist = max_dist;
            int k = order[i-1];
            
            for(int j = 1; j<(n+1);j++)
            {
            if( state[j-1] < .5)
            {
            temp1 = 0.0;
            for(int s = 0; s<p;s++)
            {
            temp1 = temp1 +  (f(k-1,s)  - f(j-1,s))*(f(k-1,s)  - f(j-1,s));
            }
            if( temp1 < dist)
            {
            dist = temp1;
            next = j;
            }
            }
            }
            order[i] = next;  
            state[next-1] = 1; 
            }
            return order;  
            }')


#tour1 =  NumericVector nn_TSP(x, p, n,sample(1:n,1))
#(f[k-1,s]  - f[j-1,s])*(f[k-1,s]  - f[j-1,s]);

########################################################
########################################################
########################################################
########################################################

cppFunction('
            NumericVector insertion_TSP(NumericMatrix f,int p, int n,int init)
            {
            NumericVector state(n,0.0);
            state[init-1] = 1;
            NumericVector order(n,0.0); 
            order[0] = init;
            
            int next;
            int u;
            int optimal_j;
            int optimal_ind;
            float dist;
            float max_dist = pow(10.0,15);
            float temp1;   
            float temp2; 
            float temp3; 
            float total_cost; 
            total_cost  = 0.0;
            double optimal_cost;
            
            
            for(int i = 1; i < n; i++)
            {
            optimal_cost = pow(10.0,10);  
            optimal_j = 1;
            optimal_ind = 1;
            
            for(int j = 1; j<(n+1);j++)
            {
            if( state[j-1] < .5)
            {
            for(int r = 0; r<i; r++)
            {
            temp1 = total_cost; 
            if( r+1 < i)
            {
            temp2 = 0.0;
            for(int s = 0; s<p;s++)
            {
            temp2 = temp2 +  (f(order[r]-1,s)  - f(order[r+1],s))*(f(order[r]-1,s)  - f(order[r+1],s));
            }   
            temp1 =  temp1 - pow(temp2,0.5);
            
            temp2 = 0.0;
            
            for(int s = 0; s<p;s++)
            {
            temp2 = temp2 +  (f(order[r]-1,s)  - f(j-1,s))*(f(order[r]-1,s)  - f(j-1,s));
            }  
            temp1 =  temp1 + pow(temp2,0.5);
            
            temp2 = 0.0;
            
            for(int s = 0; s<p;s++)
            {
            temp2 = temp2 +  (f(order[r+1]-1,s)  - f(j-1,s))*(f(order[r+1]-1,s)  - f(j-1,s));
            } 
            temp1 =  temp1 + pow(temp2,0.5);
            }
            if(r+1 == i)
            {
            temp2 = 0.0;
            for(int s = 0; s<p;s++)
            {
            temp2 = temp2 +  (f(order[r]-1,s)  - f(j-1,s))*(f(order[r]-1,s)  - f(j-1,s));
            } 
            temp1 =  temp1 + pow(temp2,0.5);                           
            }
            
            if( temp1 < optimal_cost)
            {
            optimal_cost = temp1;
            optimal_ind = r;
            optimal_j = j;  
            }
            }
            }
            }
            
            
            total_cost =  optimal_cost;
            
            for(int r = (optimal_ind +1) ; r<(i+1); r++)
            {  
            u  =    i  + optimal_ind + 1 - r;
            order[u] =  order[u-1];     
            }
            order[optimal_ind +1] =  optimal_j;
            state[optimal_j-1] = 1;
            }
            return order;  
            }')


########################################################
########################################################
########################################################
########################################################



insertion_TSP_FL = function(x,y,N,n,df = 1)
{
  p = dim(x)[2]
  f_hat2  = rep(0,n)
  f_hat22  = rep(0,n)
  
  tour_temp_inv = rep(0,n)
  tour_temp_inv2 = rep(0,n)
  
  lambda_grid =  2^seq(-3,9,length = 250)
  
  save_f_hat = matrix(0,N,n)
  
  cost = rep(0,N)
  
  for(i in 1:N)
  {
    #as.numeric(solve_TSP(atsp) )#
    temp = insertion_TSP(x, p, n,sample(1:n,1))
    #find_ordering_NN(x)
    tour_temp = temp#$order
    #cost[i] = temp$cost
    #tour_temp = as.numeric(solve_TSP(atsp,method ="nn"))
    tour_temp_inv[tour_temp] = 1:n 
    
    w = rep(0,n)
    # w = exp(-.1*sum(diff(x[tour_temp,])^2))
    
    if(df==1)
    {
      aux = tsp_npr_degrees_of_freedom(y[tour_temp],w = NULL,prop = .25,lambda_grid )
    }
    else
    {
      aux = tsp_npr(y[tour_temp],w,prop = .25,lambda_grid ) 
    }
    aux = aux[tour_temp_inv]
    
    save_f_hat[i,] = aux
    
  }###
  
  if(N == 1)
  {
    return(as.vector(save_f_hat))
  }
  
  # apply(save_f_hat[1:N,], 2,mean)
  return(list(est  = apply(save_f_hat[1:N,], 2,mean),raw = save_f_hat) )
}


###############################################################
########################################################



Linfinity_distance = function(x1,x2)
{
  if(is.matrix(x2))
  {
    n2 =  dim(x2)[1]
    dist = rep(0,n2)
    
    for(i in 1:n2)
    {
      dist[i] =  max(abs( x1 -x2[i,] ))
    }
  }
  if(is.matrix(x2)==0)
  {
    dist =  max(abs( x1 -x2))
  }
  
  return(dist)
}

########################################################################
########################################################################
######################################################
#############################################

graph_consturction = function(x,n,p)
{
  delta =  (n)^{-1/p}
  delta = 1/ceiling(1/delta)
  loc = seq(0,1,by = delta)  
  loc = loc[which(loc!=0)]
  loc = loc[which(loc!=1)]
  loc = c(0,loc,1)
  mids =  loc[1:(length(loc)-1)] + diff(loc)/2
  len_mids = length(mids)
  
  p_prime =  ceiling(log(n)/log(1/delta))
  
  indices = rep(1,p_prime)
  count  = 1
  x_grid = matrix(0,length(mids)^p_prime,p_prime)
  x_grid[1,] = rep(mids[1],p_prime)
  
  while(min(indices)<len_mids)
  {
    temp =  which(indices < len_mids)
    
    if(length(temp) > 0 )
    {
      aux = which(indices == len_mids)
      if(length(aux)==0)
      {
        indices[temp[length(temp )]] = indices[temp[length(temp )]]  +1
        count = count + 1
        x_grid[count,] =  mids[indices]
        #rbind(x_grid,mids[indices])
      }
      if(length(aux)>0)
      {
        aux2  = which(aux > temp[length(temp )])
        if( length(aux2)>0 )
        {
          indices[temp[length(temp )]] = indices[temp[length(temp )]]  +1
          indices[aux[aux2]] = 1
          count = count + 1
          x_grid[count,] =  mids[indices]
          #x_grid = rbind(x_grid,mids[indices])            
        }
        else
        {
          indices[temp[length(temp )]] = indices[temp[length(temp )]]  +1
          # x_grid = rbind(x_grid,mids[indices])
          count = count + 1
          x_grid[count,] =  mids[indices]
        }
      }
    }
  }
  colnames(x_grid) = NULL
  rownames(x_grid) = NULL
  #mids
  #x_grid
  # system.time({ ind_closest = find_closest_point_V2(x_grid,x) })
  x_grid_new = cbind(x_grid,matrix(0.5,dim(x_grid)[1],p-p_prime) )
  ind_closest = find_closest_point_V2(x_grid_new,x)
  # ind_closest = ind_closest[,1]
  #x[2,]
  #x_grid[ind_closest[2],]
  
  #   n_grid =  dim(x_grid)[1]
  #   #adj =  matrix(0,n_grid,n_grid)
  #   i_ind =  c()
  #   j_ind = c()
  #   for(i in 1:n_grid)
  #   {
  #     temp   =  Linfinity_distance( x_grid[i,],x_grid)
  #     #temp2 =  order(temp)
  #     ind = which(temp == delta) 
  #     i_ind =  c(i_ind,rep(i,length(ind)))
  #     j_ind =  c(j_ind,ind)
  #     #adj[i,ind] = 1
  #   }  
  
  n_grid =  dim(x_grid)[1]
  num_edges = 0
  for(i in 1:n_grid)
  {
    temp   =  L2_distance( x_grid[i,],x_grid)
    temp[i] =  10
    #temp2 =  order(temp)
    ind = which(temp < delta + 10^-6) 
    num_edges  = num_edges  +   length(ind)
    #i_ind =  c(i_ind,rep(i,length(ind)))
    #j_ind =  c(j_ind,ind)
    #adj[i,ind] = 1
  } 
  i_ind =  rep(0,num_edges)
  j_ind =   rep(0,num_edges)
  #  weight =  rep(0,nmum_edges)
  count = 0
  for(i in 1:n_grid)
  {
    temp   =  L2_distance( x_grid[i,],x_grid)
    temp[i] = 10
    #temp2 =  order(temp)
    ind = which(temp < delta + 10^-6)  
    i_ind[(count+1):(count + length(ind))] =  rep(i,length(ind))
    j_ind[(count+1):(count + length(ind))] =  ind
    count =  count + length(ind)
    #j_ind =  c(j_ind,ind)
    #adj[i,ind] = 1
  } 
  adj =  sparseMatrix(i = i_ind, j = j_ind,  x = rep(1,length(i_ind)), dims  = c(n_grid,n_grid))
  
  #   for(i in 1:n_grid)
  #   {
  #     temp   =  Linfinity_distance( x_grid[i,],x_grid)
  #     #temp2 =  order(temp)
  #     ind = which(temp == delta)  
  #     adj[i,ind] = 1
  #   }
  G = graph_from_adjacency_matrix(adj,"undirected")
  
  return(list(G = G,n_grid = n_grid,x_grid = x_grid_new,adj = adj))
}


#####################################################3
##############################################
##############################################

GDFS_FL = function(x,y,N,n,df = 1)
{
  p = dim(x)[2]
  
  #system.time({temp = graph_consturction(x,n,p)})
  temp = graph_consturction(x,4*n,p)
  G = temp$G
  n_grid =  temp$n_grid
  x_grid =  temp$x_grid
  adj  =   temp$adj
  
  
  ind_closest  = find_closest_point(x_grid,x)
  ind_closest = ind_closest[,1]
  
  plot(f0_func(x_grid[ind_closest,],2) - f0_func(x,2) )
  
  ptm <- proc.time()
  lambda_grid =  2^seq(-3,9,length = 370)#350
  
  save_f_hat = matrix(0,N,n)
  
  cost = rep(0,N)
  
  for(i in 1:N)
  {
    perm  = sample(1:n_grid,n_grid) 
    perm_inv = 1:n_grid
    perm_inv[perm] = 1:n_grid
    
    adj_perm =  adj[perm,perm]
    G_perm = graph_from_adjacency_matrix(adj_perm,"undirected")
    temp=  dfs(G_perm,sample(1:n_grid,1) )
    dfs_tour  = as.numeric(temp$order)
    dfs_tour =  perm[dfs_tour]
    sum(abs(diff(f0_func(x_grid[dfs_tour,],2))))
    dfs_tour[1:2]
    #temp=  dfs(G,sample(1:n_grid,1) )
    #dfs_tour  = as.numeric(temp$order)
    dfs_tour_inv  =  1:n_grid
    dfs_tour_inv[dfs_tour] = 1:n_grid 
    
    tour_temp =   order(dfs_tour_inv[ind_closest])
    tour_temp_inv =  1:n
    tour_temp_inv[tour_temp] = 1:n
    f0_tilde = f0_func(x_grid,ind)
    sum(abs(diff(f0[tour_temp])))
    sum(abs(diff(f0_tilde)))
    sum(abs(diff(f0_tilde[dfs_tour])))
    sum(abs(diff(f0[sample(1:n,n)])))
    sum(abs(diff(f0[nn_TSP(x, p, n,sample(1:n,1))]))) 
    w = rep(0,n)
    # w = exp(-.1*sum(diff(x[tour_temp,])^2))
    
    
    cost[i] =   sum(sqrt(rowSums((diff(x[tour_temp,])^2))))
    
    if(df==1)
    {
      aux = tsp_npr_degrees_of_freedom(y[tour_temp],w = NULL,prop = .25,lambda_grid )
    }
    else
    {
      aux = tsp_npr(y[tour_temp],w,prop = .25,lambda_grid ) 
    }
    aux = aux[tour_temp_inv]
    
    save_f_hat[i,] = aux
    
  }###
  proc.time() - ptm
  
  if(N == 1)
  {
    return(list(est = as.vector(save_f_hat)))
  }
  
  # apply(save_f_hat[1:N,], 2,mean)
  #apply(save_f_hat[1:N,], 2,mean)
  #save_f_hat[which.min(cost),]
  return(list(est  = apply(save_f_hat[1:N,], 2,mean),raw = save_f_hat) )
}



# adj = matrix(0,10,10)
# adj[1,c(2,3)] =  1
# adj[2,c(1,4,5)] =1
# adj[3,c(1,10)] =1
# adj[4,c(2,6,7)] =1
# adj[5,c(2,8,9)] =1
# adj[6,c(4)] =1
# adj[7,c(4)] =1
# adj[8,c(5)] =1
# adj[9,c(5)] =1
# adj[10,c(3)] =1
# 
# adj
# perm = sample(1:10,10)#c(4,1,3,2)
# perm_inv =  1:10
# perm_inv[perm] = 1:10
# perm
# perm_inv
# adj[perm,perm]
# G = graph_from_adjacency_matrix(adj,"undirected")
# plot(G)
# 
# adj_perm =  adj[perm,perm]
# G_perm = graph_from_adjacency_matrix(adj_perm,"undirected")
# temp=  dfs(G_perm,sample(1:10,1) )
# dfs_tour  = as.numeric(temp$order)
# dfs_tour =  perm[dfs_tour]
# dfs_tour 
#init =  sample(1:10,1)
#init
#dfs(G,root = init)$order

###########################################################
####################################################

find_closest_point_V2 =  function(x,x_test)
{
  n = dim(x)[1]
  n_test = dim(x_test)[1]
  closet_ind = rep(0,n_test)
  
  for( i in 1:n_test)
  {
    d = apply(x,1,function(t){ sum(abs(t-x_test[i,])) })
    closet_ind[i] = which.min(d)
  }
  return(closet_ind)
}
######################################################
######################################################
######################################################
######################################################
#   while(min(indices)<len_mids)
#   {
#     temp =  which(indices < len_mids)
#     
#     if(length(temp) > 0 )
#     {
#       aux = which(indices == len_mids)
#       if(length(aux)==0)
#       {
#         indices[temp[length(temp )]] = indices[temp[length(temp )]]  +1
#         x_grid = rbind(x_grid,mids[indices])
#       }
#       if(length(aux)>0)
#       {
#         aux2  = which(aux > temp[length(temp )])
#         if( length(aux2)>0 )
#         {
#           indices[temp[length(temp )]] = indices[temp[length(temp )]]  +1
#           indices[aux[aux2]] = 1
#           x_grid = rbind(x_grid,mids[indices])            
#         }
#         else
#         {
#           indices[temp[length(temp )]] = indices[temp[length(temp )]]  +1
#           x_grid = rbind(x_grid,mids[indices])
#         }
#       }
#     }
#   }
dfs_order_2d =  function(beta0)
{
  n =  dim(beta0)[1]*dim(beta0)[2]
  
  i = c()
  j = c()
  for(ind in 1:dim(beta0)[2])
  {
    i  = c(i, (ind-1)*dim(beta0)[1] + 1:(dim(beta0)[1] - 1) )  
    j  = c(j, (ind-1)*dim(beta0)[1] + 2:dim(beta0)[1] )  
  }
  for(ind in 1:dim(beta0)[1])
  {
    i  = c(i, dim(beta0)[1]*(0:(dim(beta0)[2]-2))  + ind )  
    j  = c(j, dim(beta0)[1]*(1:(dim(beta0)[2]-1))  + ind  )  
  }
  # x = rep(1,length(i))
  #A = sparseMatrix(i= i,j=j,x=x,dims=c(n,n))
  
  x = rep(1,length(i))
  A = sparseMatrix(i= i,j=j,x=x,dims=c(n,n))
  indices = sample(1:n,n, replace = FALSE, prob = NULL)
  gr = graph_from_adjacency_matrix(A[indices,indices],mode="undirected")
  DFS_gr  = dfs(gr,sample(1:n, 1, replace = FALSE, prob = NULL))
  dfs_order = indices[DFS_gr$order]
  
  return(dfs_order)
}



new_create_adjacency =  function(V, w, n, method = "ama") 
{
  if (!is.null(method) && !(method %in% c("ama", "admm"))) 
    stop("method must be 'ama', 'admm', or NULL.")
  differences <- apply(V, 2, FUN = function(x) {
    # norm(as.matrix(x), "f")
    mean((x)^2)
  })
  connected_ix <- which(differences <2^-4 )
  if (method == "ama") {
    ix <- vec2tri(which(w > 0), n)
  }
  else {
    ix <- vec2tri(1:(n * (n - 1)/2), n)
  }
  i <- ix[connected_ix, 1]
  j <- ix[connected_ix, 2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j - 1) * n + i] <- 1
  return(A)
}


########################################################################
########################################################################
########################################################################
########################################################################

TwoD_BernFL_bic = function(M_tour,lambda_grid,n)
{
  ## M_tour =  M[tour,tour]
  df =rep(0,length(lambda_grid))
  BIC = rep(0,length(lambda_grid))
  minus_lik = rep(0,length(lambda_grid))
  
  D = getD2dSparse(n,n)
  #sum(abs(  D%*% as.vector(theta0[tour,tour])   ))/n^2
  D = summary(D)
  ind_i = D$i
  ind_j = D$j
  D= 0 
  
  aux = ind_j[order(ind_i)]
  ind_i = aux[2*(1:(length(aux)/2))]
  ind_j = aux[2*(1:(length(aux)/2)) - 1]
  
  best_sol = matrix(0,n,n)
  best_BIC = Inf
  
  for(j in 1:length(lambda_grid))
  {
    # temp = fusedlasso_logit_admm2d(M[tour,tour], N = matrix(1.1,n,n), lambda = lambda_grid[j], iter_max = 500, rel_tol = 1e-3, inflate=1.5) 
    if(j==1)
    {
      temp = fusedlasso_weightedl2_admm2d(M_tour,lambda = lambda_grid[j])
    }
    
    if(j >1)
    {
      temp = fusedlasso_weightedl2_admm2d(M_tour,lambda = lambda_grid[j],x_init= temp$x)
    }
    #   temp = temp$x
    #   temp = as.vector(temp)
    #   temp[which(temp>1)] = 1
    #   temp[which(temp<0)] = 0  
    #   temp = matrix(temp,n(,n) 
    theta_hat =  as.vector(temp$x)
    
    nz_ind = which(abs(theta_hat[ind_i] - theta_hat[ind_j] )<10^-2)
    new_ind_i =  ind_i[nz_ind]
    new_ind_j =  ind_j[nz_ind]
    A  = sparseMatrix(i = new_ind_i, j = new_ind_j,x = rep(1,length(new_ind_i)), dims = c(n^2,n^2),symmetric = TRUE)
    G = graph_from_adjacency_matrix(A,mode = "undirected")  
    comp = components(G)
    comp = comp$csize
    df[j] = length(comp)
    theta_hat_mat = matrix(theta_hat,n,n) 
    # sigma_hat = mean((M[tour,tour] - theta_hat_mat)^2)
    # path[j,,] = theta_hat
    #    MSE[j] = mean((theta_hat_mat - theta0[tour,tour])^2)
    BIC[j] =  -2*sum(log(dbinom(M_tour,1,theta_hat_mat ))) + log(n^2)*df[j] 
    minus_lik[j]  = -sum(log(dbinom(M_tour,1,theta_hat_mat )))   
    
    if(BIC[j]< best_BIC)
    {
      best_BIC =  BIC[j]
      best_sol = theta_hat_mat
    }
    print("j")
    print(j)
  }### close for lambda
  
  return(list(best_sol = best_sol,BIC= BIC) )
  
  
}### close function


########################################################################
########################################################################
########################################################################
########################################################################

TwoD_BernFL_pred = function(M_tour,lambda_grid,n)
{
  ## M_tour =  M[tour,tour]
  df =rep(0,length(lambda_grid))
  pred = rep(0,length(lambda_grid))
  minus_lik = rep(0,length(lambda_grid))
  
  non_Missing = 1*(matrix(runif(n*n),n,n)<.8)
  M_obs = non_Missing*M_tour 
  Missing  = 1-non_Missing 
  
  
  #D = getD2dSparse(n,n)
  #sum(abs(  D%*% as.vector(theta0[tour,tour])   ))/n^2
  #D = summary(D)
  #ind_i = D$i
  #ind_j = D$j
  #D= 0 
  
  #aux = ind_j[order(ind_i)]
  #ind_i = aux[2*(1:(length(aux)/2))]
  #ind_j = aux[2*(1:(length(aux)/2)) - 1]
  
  best_sol = matrix(0,n,n)
  best_pred = Inf
  indices = which(Missing>0)
  
  for(j in 1:length(lambda_grid))
  {
    # temp = fusedlasso_logit_admm2d(M[tour,tour], N = matrix(1.1,n,n), lambda = lambda_grid[j], iter_max = 500, rel_tol = 1e-3, inflate=1.5) 
    if(j==1)
    {
      temp = fusedlasso_weightedl2_admm2d(M_obs,lambda = lambda_grid[j])
    }
    
    if(j >1)
    {
      temp = fusedlasso_weightedl2_admm2d(M_obs,lambda = lambda_grid[j],x_init= temp$x)
    }
    #   temp = temp$x
    #   temp = as.vector(temp)
    #   temp[which(temp>1)] = 1
    #   temp[which(temp<0)] = 0  
    #   temp = matrix(temp,n(,n) 
    theta_hat =  as.vector(temp$x)
    
    theta_hat_mat = matrix(theta_hat,n,n) 
    # sigma_hat = mean((M[tour,tour] - theta_hat_mat)^2)
    # path[j,,] = theta_hat
    #    MSE[j] = mean((theta_hat_mat - theta0[tour,tour])^2)
    pred[j] = -sum(log(dbinom(M_tour[indices],1,theta_hat_mat[indices] ))) 
    #mean((Missing*theta_hat_mat - Missing*M_tour)^2)
    #  -2*sum(log(dbinom(M_tour,1,theta_hat_mat ))) + log(n^2)*df[j] 
    #minus_lik[j]  = -sum(log(dbinom(M_tour,1,theta_hat_mat )))   
    
    if(pred[j]< best_pred)
    {
      best_pred =  pred[j]
      best_sol = theta_hat_mat
    }
    # print("j")
    #print(j)
  }### close for lambda
  
  ind = which.min(pred)
  temp = fusedlasso_weightedl2_admm2d(M_tour,lambda = lambda_grid[ind],x_init= temp$x)
  best_sol = temp$x
  
  return(list(best_sol = best_sol,pred=pred) )
  
  
}### close function
########################################################################





TwoD_BernFL_pred_new = function(M_tour,lambda_grid)
{
  n = dim(M_tour)[1]
  p = dim(M_tour)[2]  
  ## M_tour =  M[tour,tour]
  df =rep(0,length(lambda_grid))
  pred = rep(0,length(lambda_grid))
  minus_lik = rep(0,length(lambda_grid))
  
  non_Missing = 1*(matrix(runif(n*p),n,p)<.8)
  M_obs = non_Missing*M_tour 
  Missing  = 1-non_Missing 
  
  
  #D = getD2dSparse(n,n)
  #sum(abs(  D%*% as.vector(theta0[tour,tour])   ))/n^2
  #D = summary(D)
  #ind_i = D$i
  #ind_j = D$j
  #D= 0 
  
  #aux = ind_j[order(ind_i)]
  #ind_i = aux[2*(1:(length(aux)/2))]
  #ind_j = aux[2*(1:(length(aux)/2)) - 1]
  
  best_sol = matrix(0,n,n)
  best_pred = Inf
  indices = which(Missing>0)
  
  for(j in 1:length(lambda_grid))
  {
    # temp = fusedlasso_logit_admm2d(M[tour,tour], N = matrix(1.1,n,n), lambda = lambda_grid[j], iter_max = 500, rel_tol = 1e-3, inflate=1.5) 
    if(j==1)
    {
      temp = fusedlasso_weightedl2_admm2d(M_obs,lambda = lambda_grid[j])
    }
    
    if(j >1)
    {
      temp = fusedlasso_weightedl2_admm2d(M_obs,lambda = lambda_grid[j],x_init= temp$x)
    }
    #   temp = temp$x
    #   temp = as.vector(temp)
    #   temp[which(temp>1)] = 1
    #   temp[which(temp<0)] = 0  
    #   temp = matrix(temp,n(,n) 
    theta_hat =  as.vector(temp$x)
    
    theta_hat_mat = matrix(theta_hat,n,p) 
    # sigma_hat = mean((M[tour,tour] - theta_hat_mat)^2)
    # path[j,,] = theta_hat
    #    MSE[j] = mean((theta_hat_mat - theta0[tour,tour])^2)
    pred[j] = -sum(log(dbinom(M_tour[indices],1,theta_hat_mat[indices] ))) 
    #mean((Missing*theta_hat_mat - Missing*M_tour)^2)
    #  -2*sum(log(dbinom(M_tour,1,theta_hat_mat ))) + log(n^2)*df[j] 
    #minus_lik[j]  = -sum(log(dbinom(M_tour,1,theta_hat_mat )))   
    
    if(pred[j]< best_pred)
    {
      best_pred =  pred[j]
      best_sol = theta_hat_mat
    }
    print("j")
    print(j)
  }### close for lambda
  
  ind = which.min(pred)
  temp = fusedlasso_weightedl2_admm2d(M_tour,lambda = lambda_grid[ind],x_init= temp$x)
  best_sol = temp$x
  
  return(list(best_sol = best_sol,pred=pred) )
  
  
}### close function


cppFunction('
            NumericVector new_nn_TSP(NumericMatrix f,NumericMatrix inner,int p, int n,int init)
            {
            NumericVector state(n,0.0);
            state[init-1] = 1;
            NumericVector order(n,0.0); 
            order[0] = init;
            
            int next;
            float dist;
            float dist2;
            float max_dist = pow(10.0,15);
            float temp1;   
            float temp2;
            
            for(int i = 1; i < n; i++)
            {
            next = 1;
            dist = max_dist;
            int k = order[i-1];
            
            for(int j = 1; j<(n+1);j++)
            {
            if( state[j-1] < .5)
            {
            temp1 = 0.0;
            for(int s = 0; s<p;s++)
            {
            temp2 = 0.0; 
            temp2 = inner(k-1,s) - inner(j-1,s);
            if(temp2 <0)
            {
            temp2 =  -temp2;
            }
            
            if(j-1 == s)
            {
            temp2 =  -max_dist;   
            }
            if(k-1 == s)
            {
            temp2 =  -max_dist;   
            }
            if(temp2 >  temp1)
            {
            temp1 =  temp2;  
            }
            }
            if( temp1 < dist)
            {
            dist = temp1;
            next = j;
            }
            }
            }
            order[i] = next;  
            state[next-1] = 1; 
            }
            return order;  
            }')
