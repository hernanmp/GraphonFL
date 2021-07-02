#source("~/Rpackages/convextools/genlasso_ADMM_cholcache.R")
sourceCpp('fl_dp.cpp')

flogit = function(x) log(x/(1-x))
ilogit = function(x) 1/{1+exp(-x)}

softthresh = function(x, lambda) {
	return(sign(x)*pmax(0, abs(x) - lambda))
}


fusedlasso_weightedl2_admm2d = function(y, w=NULL, lambda = 1, iter_max = 500, rel_tol = 1e-3, inflate=1.5, x_init = NULL) {
	# solves the fused lasso on a 2D grid graph under weighted l2 loss
	# loss = sum_i w_i (y_i - x_i)^2 where x is the variable of optimization
	# y: matrix of observations
	# N: matrix of counts
	# lambda: penalty
	# inflate: governs how rapidly we rescale the dual variables to ensure primal/dual residuals of same magnitude
	n = nrow(y)
	m = ncol(y)
	
	# Initialize
	x = x_init 
	  #matrix(mean(y), nrow=n, ncol=m)
	
	if(is.null(x_init) == 1)
	{
	  x = matrix(mean(y), nrow=n, ncol=m)
	}

	z = matrix(0, nrow=n, ncol=m) # slack variable for likelihood
	u = matrix(0, nrow=n, ncol=m) # scaled dual variable for constraint
	a = 2*lambda
	if(missing(w)) w = matrix(1, nrow=n, ncol=m)

	primal_trace = NULL
	dual_trace = NULL
	converged = FALSE
	counter = 0
	while(!converged & counter < iter_max) {

		# Update rows
		this_y = (w*y + a*(z - u)) / (w+a)
		for(i in 1:n) {
			x[i,] = fl_dp_weight(this_y[i,], a + w[i,], lambda)
		}
		
		# Update columns
		z_new = z
		for(i in 1:m) {
			z_new[,i] = fl_dp(x[,i] + u[,i], lambda)
		}
		dual_residual = a*(z_new - z)
		z = z_new
		
		# Update scaled dual variable
		primal_residual = x - z
		u = u + primal_residual
		
		# check convergence
		primal_resnorm = sqrt(mean(primal_residual^2))
		dual_resnorm = sqrt(mean(dual_residual^2))
		primal_check = primal_resnorm / (rel_tol + max(sqrt(mean(x^2)), sqrt(mean(z^2))))
		dual_check = dual_resnorm / (rel_tol + sqrt(mean(u^2)))
		
		if(dual_resnorm < rel_tol && primal_resnorm < rel_tol) {
			converged=TRUE
		}
		
 		# Update step-size parameter based on norm of primal and dual residuals
		if(primal_resnorm > 5*dual_resnorm) {
			a = inflate*a
			u = u/inflate
		} else if(dual_resnorm > 5*primal_resnorm) {
			a = a/inflate
			u = inflate*u
		}
		primal_trace = c(primal_trace, primal_resnorm)
		dual_trace = c(dual_trace, dual_resnorm)
		counter = counter+1
	}
	list(x=x, z=z, u=u, primal_trace = primal_trace, dual_trace= dual_trace, counter=counter);
}





fusedlasso_logit_admm2d = function(y, N, lambda = 1, iter_max = 500, rel_tol = 1e-3, inflate=1.5) {
	# y: matrix of successes
	# N: matrix of counts
	# lambda: penalty
	# inflate: governs how rapidly we rescale the dual variables to ensure primal/dual residuals of same magnitude
	n = nrow(y)
	m = ncol(y)
	
	# Do some pre-smoothing to get the initial guess
	prob = (y + 1)/(N + 2)
	x = flogit(prob)
	z = matrix(0, nrow=n, ncol=m) # slack variable for likelihood
	u = matrix(0, nrow=n, ncol=m) # scaled dual variable for constraint
	w = matrix(1, nrow=n, ncol=m) # weights
	a = 2*lambda

	primal_trace = NULL
	dual_trace = NULL
	converged = FALSE
	counter = 0
	while(!converged & counter < iter_max) {
		# Update weights and pseudo data
		w = N*prob*(1-prob) + 1e-12  # tiny regularization term so we don't get zero weights
		r = x - (N*prob - y)/w
		
		# Update rows
		this_y = (w*r + a*(z - u)) / (w+a)
		for(i in 1:n) {
			x[i,] = fl_dp_weight(this_y[i,], a + w[i,], lambda)
		}
		
		# Update columns
		z_new = z
		for(i in 1:m) {
			z_new[,i] = fl_dp(x[,i] + u[,i], lambda)
		}
		dual_residual = a*(z_new - z)
		z = z_new
		
		# Update scaled dual variable
		primal_residual = x - z
		u = u + primal_residual
		prob = ilogit(x)
		
		# check convergence
		primal_resnorm = sqrt(mean(primal_residual^2))
		dual_resnorm = sqrt(mean(dual_residual^2))
		primal_check = primal_resnorm / (rel_tol + max(sqrt(mean(x^2)), sqrt(mean(z^2))))
		dual_check = dual_resnorm / (rel_tol + sqrt(mean(u^2)))
		
		if(dual_resnorm < rel_tol && primal_resnorm < rel_tol) {
			converged=TRUE
		}
		
 		# Update step-size parameter based on norm of primal and dual residuals
		if(primal_resnorm > 5*dual_resnorm) {
			a = inflate*a
			u = u/inflate
		} else if(dual_resnorm > 5*primal_resnorm) {
			a = a/inflate
			u = inflate*u
		}
		primal_trace = c(primal_trace, primal_resnorm)
		dual_trace = c(dual_trace, dual_resnorm)
		counter = counter+1
	}
	list(x=x, z=z, u=u, primal_trace = primal_trace, dual_trace= dual_trace, counter=counter);
}






