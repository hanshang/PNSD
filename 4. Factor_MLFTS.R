##################################
# select the number of components
##################################

# eigenvalue: estimated eigenvalues
# index: 1 to T
# sample_size: T
# no_pop: N

select_K_new <- function(eigenvalue, index, sample_size, no_pop)
{
    return((eigenvalue[index])/sample_size + index * (max(sample_size, no_pop)^(-0.5)))
}

#################################
# Factor decomposition for HDFTS
#################################

# data: an array of dimension (T by N by p)
 
factor_model_fun <- function(data)
{
    sample_size = dim(data)[1]
    no_pop = dim(data)[2]
    D_val = dim(data)[3]

    Delta = matrix(0, sample_size, sample_size)
    for(ij in 1:no_pop)
    {
        temp = matrix(NA, sample_size, sample_size)
        for(t in 1:sample_size)
        {
            for(s in 1:sample_size)
            {
                temp[t,s] = matrix(data[t,ij,], nrow = 1) %*% matrix(data[s,ij,], ncol = 1)
            }
        }
        Delta = Delta + temp/D_val
        rm(temp)
    }
    rm(t); rm(s); rm(ij)

    Delta_mat = Delta/no_pop
    Delta_mat_eigen = eigen(Delta_mat)

    # new way of selecting the number of components

    K_val = vector("numeric", sample_size)
    for(ik in 1:sample_size)
    {
        K_val[ik] = select_K_new(eigenvalue = Delta_mat_eigen$values, index = ik, sample_size = sample_size,
                                 no_pop = no_pop)
        rm(ik)
    }
    q_val_est = which.min(K_val) - 1
    if(q_val_est == 0)
    {
        warning("The number of components is zero.")
    }
    
    ###############################
    # estimate the eigen-dimension
    ###############################

    # factors
    
    Delta_eigen_vector = as.matrix(Delta_mat_eigen$vectors[,1:q_val_est]) * sqrt(sample_size)

  	# factor loadings
	
    factor_loading = array(NA, dim = c(q_val_est, D_val, no_pop))
    for(ij in 1:no_pop)
    {
        factor_loading[,,ij] = (crossprod(Delta_eigen_vector, data[,ij,]))/sample_size
        rm(ij)
    }

    factor_decomp = error = array(NA, dim = c(sample_size, no_pop, D_val), dimnames = list(1:sample_size, 1:no_pop, 1:D_val))
    for(ij in 1:no_pop)
    {
        factor_decomp[,ij,] = (Delta_eigen_vector %*% factor_loading[,,ij])
        error[,ij,] = data[,ij,] - factor_decomp[,ij,]
        rm(ij)
    }
    
    # forecast scores
    
    if(q_val_est == 1)
    {
        # univariate time-series forecasting of the extracted factor
      
        factor_fore = forecast(ets(ts(Delta_eigen_vector, start = 1, end = length(Delta_eigen_vector))), h = 1)$mean
    }
    else
    {
        warning("There is more than one factor!")
    }
    
    # multiply forecasted factors by factor loadings for each PM size
    
    factor_decomp_fore = matrix(NA, D_val, no_pop)
    for(ij in 1:D_val)
    {
        factor_decomp_fore[ij,] = factor_fore %*% factor_loading[,ij,]
        rm(ij)
    }
    
    # model residuals (HDFTS) by MLFTS
    
    MLFTS_residual <- matrix(unlist(MLFTS_model(data_input = aperm(error, c(3, 1, 2)),
                          aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,no_pop)
    
    fore_val = factor_decomp_fore + MLFTS_residual
    return(list(factors = Delta_eigen_vector, factor_loading = factor_loading, 
                factor_decomp = factor_decomp, error = error, q_val_est = q_val_est,
                fore_val = fore_val))
}    

#######################
# factor model + MLFTS
#######################

# data_set: data array

factor_MLFTS_fun <- function(data_set)
{
    result = array(NA, dim = c(dim(data_set)[3], dim(data_set)[1], n_test))
    factor_ncomp = vector("numeric", n_test)
    for(iw in 1:n_test)
    {
        n_train = (408 * 3/4 + iw - 1)
    
        # factor model decomposition
        
        factor_decomp = factor_model_fun(data = aperm(data_set[,1:n_train,], c(2, 1, 3)))
        
        # number of factors is often 1
        
        factor_ncomp[iw] = factor_decomp$q_val_est
        result[,,iw] = factor_decomp$fore_val
        
        print(iw); rm(iw); rm(n_train); rm(factor_decomp); rm(factor_fore); rm(factor_decomp_fore);
        rm(MLFTS_residual)
    }
    return(list(result = result, ncomp = factor_ncomp))
}

factor_MLFTS_Mon_fore = factor_MLFTS_fun(data_set = PNSD_log_Mon_408)
factor_MLFTS_Tue_fore = factor_MLFTS_fun(data_set = PNSD_log_Tue_408)
factor_MLFTS_Wed_fore = factor_MLFTS_fun(data_set = PNSD_log_Wed_408)
factor_MLFTS_Thu_fore = factor_MLFTS_fun(data_set = PNSD_log_Thu_408)
factor_MLFTS_Fri_fore = factor_MLFTS_fun(data_set = PNSD_log_Fri_408)
factor_MLFTS_Sat_fore = factor_MLFTS_fun(data_set = PNSD_log_Sat_408)
factor_MLFTS_Sun_fore = factor_MLFTS_fun(data_set = PNSD_log_Sun_408)

# retained number of factors

all(factor_MLFTS_Mon_fore$ncomp == 1) # TRUE
all(factor_MLFTS_Tue_fore$ncomp == 1) # TRUE
all(factor_MLFTS_Wed_fore$ncomp == 1) # TRUE
all(factor_MLFTS_Thu_fore$ncomp == 1) # TRUE
all(factor_MLFTS_Fri_fore$ncomp == 1) # TRUE
all(factor_MLFTS_Sat_fore$ncomp == 1) # TRUE
all(factor_MLFTS_Sun_fore$ncomp == 1) # TRUE 

#################
# compute errors
#################

# MAPE

mape <- ftsa:::mape

# MAPE aggregated across 51 sizes

mape_h1_Mon_factor = mape_h1_Tue_factor = mape_h1_Wed_factor = mape_h1_Thu_factor = 
mape_h1_Fri_factor = mape_h1_Sat_factor = mape_h1_Sun_factor = matrix(NA, 24, n_test)
for(ik in 1:n_test)
{
    for(ij in 1:24)
    {
        mape_h1_Mon_factor[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Mon_fore$result) - 1)[,ij,ik], true = PNSD_Mon_test[ij,ik,])
        mape_h1_Tue_factor[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Tue_fore$result) - 1)[,ij,ik], true = PNSD_Tue_test[ij,ik,])
        mape_h1_Wed_factor[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Wed_fore$result) - 1)[,ij,ik], true = PNSD_Wed_test[ij,ik,])
        mape_h1_Thu_factor[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Thu_fore$result) - 1)[,ij,ik], true = PNSD_Thu_test[ij,ik,])
        mape_h1_Fri_factor[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Fri_fore$result) - 1)[,ij,ik], true = PNSD_Fri_test[ij,ik,])
        mape_h1_Sat_factor[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Sat_fore$result) - 1)[,ij,ik], true = PNSD_Sat_test[ij,ik,])
        mape_h1_Sun_factor[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Sun_fore$result) - 1)[,ij,ik], true = PNSD_Sun_test[ij,ik,])
        rm(ij)
    }
    rm(ik)
}

mape_h1_Mon_factor_mean = apply(mape_h1_Mon_factor, 1, mean)
mape_h1_Tue_factor_mean = apply(mape_h1_Tue_factor, 1, mean)
mape_h1_Wed_factor_mean = apply(mape_h1_Wed_factor, 1, mean)
mape_h1_Thu_factor_mean = apply(mape_h1_Thu_factor, 1, mean)
mape_h1_Fri_factor_mean = apply(mape_h1_Fri_factor, 1, mean)
mape_h1_Sat_factor_mean = apply(mape_h1_Sat_factor, 1, mean)
mape_h1_Sun_factor_mean = apply(mape_h1_Sun_factor, 1, mean)

length(which(mape_h1_Mon_factor_mean < mape_h1_Mon_TS_mean)) # 22
length(which(mape_h1_Tue_factor_mean < mape_h1_Tue_TS_mean)) # 24  
length(which(mape_h1_Wed_factor_mean < mape_h1_Wed_TS_mean)) # 24
length(which(mape_h1_Thu_factor_mean < mape_h1_Thu_TS_mean)) # 23
length(which(mape_h1_Fri_factor_mean < mape_h1_Fri_TS_mean)) # 15
length(which(mape_h1_Sat_factor_mean < mape_h1_Sat_TS_mean)) # 7
length(which(mape_h1_Sun_factor_mean < mape_h1_Sun_TS_mean)) # 12


mape_h1_days_factor_mean = rbind(mape_h1_Mon_factor_mean, mape_h1_Tue_factor_mean, mape_h1_Wed_factor_mean,
                                 mape_h1_Thu_factor_mean, mape_h1_Fri_factor_mean, mape_h1_Sat_factor_mean,
                                 mape_h1_Sun_factor_mean)
rownames(mape_h1_days_factor_mean) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(mape_h1_days_factor_mean) = 1:24

# MAPE aggregated across 24 hours

mape_h1_Mon_factor_size = mape_h1_Tue_factor_size = mape_h1_Wed_factor_size = mape_h1_Thu_factor_size = 
mape_h1_Fri_factor_size = mape_h1_Sat_factor_size = mape_h1_Sun_factor_size = matrix(NA, 51, n_test)
for(ik in 1:n_test)
{
    for(ij in 1:51)
    {
        mape_h1_Mon_factor_size[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Mon_fore$result) - 1)[ij,,ik], true = PNSD_Mon_test[,ik,ij])
        mape_h1_Tue_factor_size[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Tue_fore$result) - 1)[ij,,ik], true = PNSD_Tue_test[,ik,ij])
        mape_h1_Wed_factor_size[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Wed_fore$result) - 1)[ij,,ik], true = PNSD_Wed_test[,ik,ij])
        mape_h1_Thu_factor_size[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Thu_fore$result) - 1)[ij,,ik], true = PNSD_Thu_test[,ik,ij])
        mape_h1_Fri_factor_size[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Fri_fore$result) - 1)[ij,,ik], true = PNSD_Fri_test[,ik,ij])
        mape_h1_Sat_factor_size[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Sat_fore$result) - 1)[ij,,ik], true = PNSD_Sat_test[,ik,ij])
        mape_h1_Sun_factor_size[ij,ik] = ftsa:::mape(forecast = (10^(factor_MLFTS_Sun_fore$result) - 1)[ij,,ik], true = PNSD_Sun_test[,ik,ij])
        rm(ij)
    }
    rm(ik)
}

mape_h1_Mon_factor_size_mean = apply(mape_h1_Mon_factor_size, 1, mean)
mape_h1_Tue_factor_size_mean = apply(mape_h1_Tue_factor_size, 1, mean)
mape_h1_Wed_factor_size_mean = apply(mape_h1_Wed_factor_size, 1, mean)
mape_h1_Thu_factor_size_mean = apply(mape_h1_Thu_factor_size, 1, mean)
mape_h1_Fri_factor_size_mean = apply(mape_h1_Fri_factor_size, 1, mean)
mape_h1_Sat_factor_size_mean = apply(mape_h1_Sat_factor_size, 1, mean)
mape_h1_Sun_factor_size_mean = apply(mape_h1_Sun_factor_size, 1, mean)

length(which(mape_h1_Mon_factor_size_mean < mape_h1_Mon_TS_size_mean)) # 36
length(which(mape_h1_Tue_factor_size_mean < mape_h1_Tue_TS_size_mean)) # 29
length(which(mape_h1_Wed_factor_size_mean < mape_h1_Wed_TS_size_mean)) # 30
length(which(mape_h1_Thu_factor_size_mean < mape_h1_Thu_TS_size_mean)) # 28
length(which(mape_h1_Fri_factor_size_mean < mape_h1_Fri_TS_size_mean)) # 31
length(which(mape_h1_Sat_factor_size_mean < mape_h1_Sat_TS_size_mean)) # 26
length(which(mape_h1_Sun_factor_size_mean < mape_h1_Sun_TS_size_mean)) # 29

mape_h1_days_factor_size_mean = rbind(mape_h1_Mon_factor_size_mean,
                                      mape_h1_Tue_factor_size_mean,
                                      mape_h1_Wed_factor_size_mean,
                                      mape_h1_Thu_factor_size_mean,
                                      mape_h1_Fri_factor_size_mean,
                                      mape_h1_Sat_factor_size_mean,
                                      mape_h1_Sun_factor_size_mean)
colnames(mape_h1_days_factor_size_mean) = x_grid
rownames(mape_h1_days_factor_size_mean) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

