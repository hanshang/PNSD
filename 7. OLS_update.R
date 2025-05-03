# data: data_array 24 x T x 51
# newdata: newly arrived data 
# holdoutdata: holdout data
# order: number of retained functional principal components
# fmethod: ARIMA or ETS

ols_update <- function(data, newdata, holdoutdata, order, fmethod)
{
    n_hour = dim(data)[1]
    n_day = dim(data)[2]
    n_PM = dim(data)[3]
    
    # convert data array to data matrix via stacking
    
    data_set_mat = matrix(NA, (n_hour * n_PM), n_day)
    for(iw in 1:n_PM)
    {
        data_set_mat[((iw - 1) * n_hour + 1):(iw * n_hour), ] = data[,1:n_day,iw]
        rm(iw)
    }
    
    ftsmobject = ftsm(fts(1:(n_hour * n_PM), data_set_mat), order = order)
    coef = ftsmobject$coeff
    base = ftsmobject$basis
    
    base_array = array(NA, dim = c(n_hour, n_PM, (order + 1)))
    for(ik in 1:n_PM)
    {
        base_array[,ik,] = base[((ik - 1) * 24 + 1):(ik * 24),]
        rm(ik)
    }
    
    n2 = dim(newdata)[1]
    base1 = base_array[1:n2,,]
    base2 = array(base_array[(n2 + 1):n_hour,,], dim = c(length((n2 + 1):n_hour), n_PM, (order + 1)))
    
    ols_mat = matrix(NA, (order + 1), n_PM)
    for(ik in 1:n_PM)
    {
        ols_mat[,ik] = ginv(t(base1[,ik,]) %*% base1[,ik,]) %*% (t(base1[,ik,]) %*% newdata[,ik])
        rm(ik)
    }   
    
    update_period = (nrow(newdata) + 1):n_hour
    ols_forecasts = matrix(NA, length(update_period), n_PM)
    for(ik in 1:n_PM)
    {
        ols_forecasts[,ik] = base2[,ik,] %*% matrix(ols_mat[,ik], (order + 1), 1)
        rm(ik)
    }
    
    err = vector("numeric", n_PM)
    for(ik in 1:n_PM)
    {
        err[ik] = ftsa:::mape(forecast = (10^(ols_forecasts[,ik]) - 1), true = (10^(matrix(holdoutdata, length(update_period), n_PM)[,ik]) - 1))
        rm(ik)      
    }
    rm(ftsmobject); rm(coef); rm(base); rm(base1); rm(base2); rm(n2); rm(update_period)
    return(list(fore = t(10^(ols_forecasts) - 1), err = err))
}

##############
# evaluation 
##############

ols_update_test <- function(data_set, observed_period, fore_method)
{
    n_obs = dim(data_set)[2]
    ols_h1_err = matrix(NA, dim(data_set)[3], n_test)
    for(iw in 1:n_test)
    {
        n_train = (n_obs * 3/4 + iw - 1)
        ols_h1_err[,iw] = ols_update(data = data_set[,1:n_train,], 
                                     newdata =  data_set[1:observed_period,(n_train + 1),], 
                                     holdoutdata = data_set[(observed_period + 1):24,(n_train + 1),], 
                                     order = 6, fmethod = fore_method)$err
        rm(iw); rm(n_train)
    }
    return(ols_h1_err)
}

# forecasts for Mon-Sun

ols_log_Mon_list <- foreach(iwk = 2:23) %dopar% ols_update_test(data_set = PNSD_log_Mon_408, 
                                                                observed_period = iwk,
                                                                fore_method = "ets")

ols_log_Tue_list <- foreach(iwk = 2:23) %dopar% ols_update_test(data_set = PNSD_log_Tue_408, 
                                                                observed_period = iwk,
                                                                fore_method = "ets")

ols_log_Wed_list <- foreach(iwk = 2:23) %dopar% ols_update_test(data_set = PNSD_log_Wed_408, 
                                                                observed_period = iwk,
                                                                fore_method = "ets")

ols_log_Thu_list <- foreach(iwk = 2:23) %dopar% ols_update_test(data_set = PNSD_log_Thu_408, 
                                                                observed_period = iwk,
                                                                fore_method = "ets")

ols_log_Fri_list <- foreach(iwk = 2:23) %dopar% ols_update_test(data_set = PNSD_log_Fri_408, 
                                                                observed_period = iwk,
                                                                fore_method = "ets")

ols_log_Sat_list <- foreach(iwk = 2:23) %dopar% ols_update_test(data_set = PNSD_log_Sat_408, 
                                                                observed_period = iwk,
                                                                fore_method = "ets")

ols_log_Sun_list <- foreach(iwk = 2:23) %dopar% ols_update_test(data_set = PNSD_log_Sun_408, 
                                                                observed_period = iwk,
                                                                fore_method = "ets")

ols_log_Mon_array = ols_log_Tue_array = 
ols_log_Wed_array = ols_log_Thu_array = 
ols_log_Fri_array = ols_log_Sat_array = 
ols_log_Sun_array = array(NA, dim = c(51, 102, 22), dimnames = list(x_grid, 1:102, 3:24))
for(iwk in 1:22)
{
    ols_log_Mon_array[,,iwk] = ols_log_Mon_list[[iwk]]
    ols_log_Tue_array[,,iwk] = ols_log_Tue_list[[iwk]]
    ols_log_Wed_array[,,iwk] = ols_log_Wed_list[[iwk]]
    ols_log_Thu_array[,,iwk] = ols_log_Thu_list[[iwk]]
    ols_log_Fri_array[,,iwk] = ols_log_Fri_list[[iwk]]
    ols_log_Sat_array[,,iwk] = ols_log_Sat_list[[iwk]]
    ols_log_Sun_array[,,iwk] = ols_log_Sun_list[[iwk]]
    rm(iwk)
}

# 22 x 51

ols_log_Mon_array_mean = t(apply(ols_log_Mon_array, c(1, 3), mean))
ols_log_Tue_array_mean = t(apply(ols_log_Tue_array, c(1, 3), mean))
ols_log_Wed_array_mean = t(apply(ols_log_Wed_array, c(1, 3), mean))
ols_log_Thu_array_mean = t(apply(ols_log_Thu_array, c(1, 3), mean))
ols_log_Fri_array_mean = t(apply(ols_log_Fri_array, c(1, 3), mean))
ols_log_Sat_array_mean = t(apply(ols_log_Sat_array, c(1, 3), mean))
ols_log_Sun_array_mean = t(apply(ols_log_Sun_array, c(1, 3), mean))

# by hour

ols_log_Mon_array_mean_hour = apply(ols_log_Mon_array_mean, 1, mean)
ols_log_Tue_array_mean_hour = apply(ols_log_Tue_array_mean, 1, mean)
ols_log_Wed_array_mean_hour = apply(ols_log_Wed_array_mean, 1, mean)
ols_log_Thu_array_mean_hour = apply(ols_log_Thu_array_mean, 1, mean)
ols_log_Fri_array_mean_hour = apply(ols_log_Fri_array_mean, 1, mean)
ols_log_Sat_array_mean_hour = apply(ols_log_Sat_array_mean, 1, mean)
ols_log_Sun_array_mean_hour = apply(ols_log_Sun_array_mean, 1, mean)

ols_days_array_mean_hour = rbind(ols_log_Mon_array_mean_hour,
                                 ols_log_Mon_array_mean_hour,
                                 ols_log_Mon_array_mean_hour,
                                 ols_log_Mon_array_mean_hour,
                                 ols_log_Mon_array_mean_hour,
                                 ols_log_Mon_array_mean_hour,
                                 ols_log_Mon_array_mean_hour)
colnames(ols_days_array_mean_hour) = 3:24
rownames(ols_days_array_mean_hour) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  
# by size

ols_log_Mon_array_mean_size = apply(ols_log_Mon_array_mean, 2, mean)
ols_log_Tue_array_mean_size = apply(ols_log_Tue_array_mean, 2, mean)
ols_log_Wed_array_mean_size = apply(ols_log_Wed_array_mean, 2, mean)
ols_log_Thu_array_mean_size = apply(ols_log_Thu_array_mean, 2, mean)
ols_log_Fri_array_mean_size = apply(ols_log_Fri_array_mean, 2, mean)
ols_log_Sat_array_mean_size = apply(ols_log_Sat_array_mean, 2, mean)
ols_log_Sun_array_mean_size = apply(ols_log_Sun_array_mean, 2, mean)

ols_days_array_mean_size = rbind(ols_log_Mon_array_mean_size, 
                                 ols_log_Tue_array_mean_size, 
                                 ols_log_Wed_array_mean_size, 
                                 ols_log_Thu_array_mean_size, 
                                 ols_log_Fri_array_mean_size, 
                                 ols_log_Sat_array_mean_size, 
                                 ols_log_Sun_array_mean_size)
colnames(ols_days_array_mean_size) = x_grid
rownames(ols_days_array_mean_size) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

