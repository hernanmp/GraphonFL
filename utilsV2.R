

f0_1  =  function(xi,a,b,K)
{
  n =  length(xi)
  Z =  rep(0,length(n))
  
  #   for(i in 1:n)
  #   {
  #     cum_sum = 0
  #     
  #     for(k in 1:K)
  #     {
  #       cum_sum = cum_sum + a[k]
  #       if(cum_sum > xi[i])
  #       {
  #         j = k
  #         break;
  #       }
  #     }
  #     Z[i] = j
  #   }
  breakpoints =  rep(0,K)
  for(k in 1:K)
  {
    ind =  which(xi<a[k])
    if(k==1)
    {
      breakpoints[k] = ind[length(ind)] 
    }
    #####
    else
    {
      breakpoints[k] = breakpoints[k] +  ind[length(ind)]
    }
    #####################
  }
  P = matrix(0,n,n)
  for( k1 in 1:K)
  {
    for( k2 in 1:K)
    {
      ############
      if(k1 ==1)
      {
        I1 =  1:breakpoints[1]
      }
      if(k1 > 1  &&  k1 < K )
      {
        I1 =  (breakpoints[k1-1]+1):breakpoints[k1]
      }
      if(k1 == K)
      {
        I1 =  (breakpoints[k1-1]+1):n
      }
      if(k2 ==1)
      {
        I2  =  1:breakpoints[1]
      }
      if(k2 > 1  &&  k2 < K )
      {
        I2 =  (breakpoints[k2-1]+1):breakpoints[k2]
      }
      if(k2 == K)
      {
        I2 =  (breakpoints[k2-1]+1):n
      }
      ##################3
      P[I1,I2] =  b[k1,k2]
    }######
  } 
  ##########################################
  
  
  return(P)
}

###########################################################




f0_2  =  function(xi,a,b,K)
{
  n =  length(xi)
  Z =  rep(0,length(n))
  
  #   for(i in 1:n)
  #   {
  #     cum_sum = 0
  #     
  #     for(k in 1:K)
  #     {
  #       cum_sum = cum_sum + a[k]
  #       if(cum_sum > xi[i])
  #       {
  #         j = k
  #         break;
  #       }
  #     }
  #     Z[i] = j
  #   }
  breakpoints =  rep(0,K)
  for(k in 1:K)
  {
    ind =  which(xi<a[k])
    if(k==1)
    {
      breakpoints[k] = ind[length(ind)] 
    }
    #####
    else
    {
      breakpoints[k] = breakpoints[k] +  ind[length(ind)]
    }
    #####################
  }
  P = matrix(0,n,n)
  for( k1 in 1:K)
  {
    for( k2 in 1:K)
    {
      ############
      if(k1 ==1)
      {
        I1 =  1:breakpoints[1]
      }
      if(k1 > 1  &&  k1 < K )
      {
        I1 =  (breakpoints[k1-1]+1):breakpoints[k1]
      }
      if(k1 == K)
      {
        I1 =  (breakpoints[k1-1]+1):n
      }
      if(k2 ==1)
      {
        I2  =  1:breakpoints[1]
      }
      if(k2 > 1  &&  k2 < K )
      {
        I2 =  (breakpoints[k2-1]+1):breakpoints[k2]
      }
      if(k2 == K)
      {
        I2 =  (breakpoints[k2-1]+1):n
      }
      ##################3
      #P[I1,I2] =  b[k1,k2]
      
      for(ii in 1:length(I1))
      {
        for(jj in 1:length(I2))
        {
          i = I1[ii]
          j = I2[jj]
          P[i,j] = (b[k1]*sqrt(xi[i])  + (1-b[k1])*sqrt(xi[j]) )/2 
          P[i,j] =  P[i,j] + (b[k2]*sqrt(xi[j])  + (1-b[k2])*sqrt(xi[i]) )/2 
        }
      }
      
      
      
      
    }######
  } 
  ##########################################
  
  
  return(P)
}
###########################################################




f0_3  =  function(xi,a,b,K)
{
  n =  length(xi)
  Z =  rep(0,length(n))
  

  ##########################################
  P =  matrix(0.05,n,n)
  for(i in 1:n)
  {
     for(j in 1:n)
     {
        for(k in 1:length(a))
        {
          temp  = -  xi[j] + a[k]
          if(xi[i]< temp)
          {
            P[i,j] =  (1 - k%%2 )/2 + 0.05
            break;
          }
        }
     }
  }
  
  return(P)
}

###########################################################




f0_4  =  function(xi,a,b,K)
{
  n =  length(xi)
  Z =  rep(0,length(n))
  
  
  ##########################################
  P =  matrix(0,n,n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      for(k in 1:length(a))
      {
      #  temp  = (cos(xi[i]*pi*a[k]) )/2
        if(xi[j]< a[k])
        {
          P[i,j]   =  xi[i]^2  +  xi[j]^2  -  2*(xi[i]*xi[j])# min(0.5,(cos(xi[i]*6*pi*a[k]))/2)
          P[i,j] =  (xi[i] + xi[j]+P[i,j] )/3 #max(P[i,j] ,0.05)
            #min(.5*(a[k] +  xi[i]^2  +  xi[j]^2  -  2*(xi[i]*xi[j])),.6) #+ cos(xi[j]*xi[i]*2*pi)
            #cos((xi[i]+xi[j])*pi*2) + sin((xi[i]+xi[j])*pi*2)
            #min(0.5,(cos(xi[i]*pi*a[k]))/2)
          break;
        }
      }
    }
  }
  image(P)
  
  return(P)
}


