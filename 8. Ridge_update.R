# lambda: tuning parameter
# data: data_array 24 x T x 51
# newdata: newly arrived data 
# holdoutdata: holdout data
# order: number of retained functional principal components
# fmethod: ARIMA or ETS

ridge_update <- function(lambda, data, newdata, holdoutdata, order, fmethod)
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
    
    ridge_mat = matrix(NA, (order + 1), n_PM)
    for(ik in 1:n_PM)
    {
        I = diag(dim(base)[2])
        ridge_mat[,ik] = ginv(t(base1[,ik,]) %*% base1[,ik,] + lambda * I) %*% (t(base1[,ik,]) %*% newdata[,ik])
        rm(ik); rm(I)
    }   
    
    update_period = (nrow(newdata) + 1):n_hour
    ridge_forecasts = matrix(NA, length(update_period), n_PM)
    for(ik in 1:n_PM)
    {
        ridge_forecasts[,ik] = base2[,ik,] %*% matrix(ridge_mat[,ik], (order + 1), 1)
        rm(ik)
    }
    
    err = vector("numeric", n_PM)
    for(ik in 1:n_PM)
    {
        err[ik] = ftsa:::mape(forecast = (10^(ridge_forecasts[,ik]) - 1), 
                              true = (10^(matrix(holdoutdata, length(update_period), n_PM)[,ik]) - 1))
        rm(ik)      
    }
    rm(data)
    return(list(fore = t(10^(ridge_forecasts) - 1), err = err))
}

# validation set: 205:306

ridge_update_validation <- function(lambda_val, data_set, observed_period, fore_method)
{
    ridge_h1_err = matrix(NA, dim(data_set)[3], n_test)
    for(iw in 1:n_test)
    {
        n_train = (408 * 1/2 + iw - 1)
        ridge_h1_err[,iw] = ridge_update(lambda = lambda_val, data = data_set[,1:n_train,], 
                                     newdata =  data_set[1:observed_period,(n_train + 1),], 
                                     holdoutdata = data_set[(observed_period + 1):24,(n_train + 1),], 
                                     order = 6, fmethod = fore_method)$err
        rm(iw)
    }
    return(round(mean(ridge_h1_err), 4))
}

# ridge forecasts (Mon - Sun)

registerDoMC(22)
ridge_lambda_log_Mon_list = foreach(iwk = 2:23) %dopar% optimise(f = ridge_update_validation, 
                                    interval = c(0, 10^3), data_set = PNSD_log_Mon_408, observed_period = iwk)

ridge_lambda_log_Tue_list = foreach(iwk = 2:23) %dopar% optimise(f = ridge_update_validation, 
                                    interval = c(0, 10^3), data_set = PNSD_log_Tue_408, observed_period = iwk)

ridge_lambda_log_Wed_list = foreach(iwk = 2:23) %dopar% optimise(f = ridge_update_validation, 
                                    interval = c(0, 10^3), data_set = PNSD_log_Wed_408, observed_period = iwk)

ridge_lambda_log_Thu_list = foreach(iwk = 2:23) %dopar% optimise(f = ridge_update_validation, 
                                    interval = c(0, 10^3), data_set = PNSD_log_Thu_408, observed_period = iwk)

ridge_lambda_log_Fri_list = foreach(iwk = 2:23) %dopar% optimise(f = ridge_update_validation, 
                                    interval = c(0, 10^3), data_set = PNSD_log_Fri_408, observed_period = iwk)

ridge_lambda_log_Sat_list = foreach(iwk = 2:23) %dopar% optimise(f = ridge_update_validation, 
                                    interval = c(0, 10^3), data_set = PNSD_log_Sat_408, observed_period = iwk)

ridge_lambda_log_Sun_list = foreach(iwk = 2:23) %dopar% optimise(f = ridge_update_validation, 
                                    interval = c(0, 10^3), data_set = PNSD_log_Sun_408, observed_period = iwk)


ridge_lambda_log_Mon_para = ridge_lambda_log_Tue_para = ridge_lambda_log_Wed_para = 
ridge_lambda_log_Thu_para = ridge_lambda_log_Fri_para = ridge_lambda_log_Sat_para = ridge_lambda_log_Sun_para = 
ridge_lambda_log_Mon_obj = ridge_lambda_log_Tue_obj = ridge_lambda_log_Wed_obj = 
ridge_lambda_log_Thu_obj = ridge_lambda_log_Fri_obj = ridge_lambda_log_Sat_obj = ridge_lambda_log_Sun_obj = vector("numeric", 22)
for(iwk in 1:22)
{
    ridge_lambda_log_Mon_para[iwk] = ridge_lambda_log_Mon_list[[iwk]]$minimum
    ridge_lambda_log_Tue_para[iwk] = ridge_lambda_log_Tue_list[[iwk]]$minimum
    ridge_lambda_log_Wed_para[iwk] = ridge_lambda_log_Wed_list[[iwk]]$minimum
    ridge_lambda_log_Thu_para[iwk] = ridge_lambda_log_Thu_list[[iwk]]$minimum
    ridge_lambda_log_Fri_para[iwk] = ridge_lambda_log_Fri_list[[iwk]]$minimum
    ridge_lambda_log_Sat_para[iwk] = ridge_lambda_log_Sat_list[[iwk]]$minimum
    ridge_lambda_log_Sun_para[iwk] = ridge_lambda_log_Sun_list[[iwk]]$minimum
    
    ridge_lambda_log_Mon_obj[iwk] = ridge_lambda_log_Mon_list[[iwk]]$objective
    ridge_lambda_log_Tue_obj[iwk] = ridge_lambda_log_Tue_list[[iwk]]$objective
    ridge_lambda_log_Wed_obj[iwk] = ridge_lambda_log_Wed_list[[iwk]]$objective
    ridge_lambda_log_Thu_obj[iwk] = ridge_lambda_log_Thu_list[[iwk]]$objective
    ridge_lambda_log_Fri_obj[iwk] = ridge_lambda_log_Fri_list[[iwk]]$objective
    ridge_lambda_log_Sat_obj[iwk] = ridge_lambda_log_Sat_list[[iwk]]$objective
    ridge_lambda_log_Sun_obj[iwk] = ridge_lambda_log_Sun_list[[iwk]]$objective
    rm(iwk)
}

ridge_lambda_days_para = cbind(ridge_lambda_log_Mon_para, ridge_lambda_log_Tue_para, ridge_lambda_log_Wed_para,
                               ridge_lambda_log_Thu_para, ridge_lambda_log_Fri_para, ridge_lambda_log_Sat_para,
                               ridge_lambda_log_Sun_para)
rownames(ridge_lambda_days_para) = 3:24
colnames(ridge_lambda_days_para) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

savefig("Fig_3a", width = 12, height = 10, toplines = 0.8, type = "png")
plot(3:24, ridge_lambda_days_para[,1], col = 1, lty = 1, pch = 1, ylim = c(0, 11), xlab = "Intraday period", ylab = "Tuning parameter", main = "Ridge estimator")
points(3:24, ridge_lambda_days_para[,2], col = 2, lty = 2, pch = 2)
points(3:24, ridge_lambda_days_para[,3], col = 3, lty = 3, pch = 3)
points(3:24, ridge_lambda_days_para[,4], col = 4, lty = 4, pch = 4)
points(3:24, ridge_lambda_days_para[,5], col = 5, lty = 5, pch = 5)
points(3:24, ridge_lambda_days_para[,6], col = 6, lty = 6, pch = 6)
points(3:24, ridge_lambda_days_para[,7], col = 7, lty = 7, pch = 7)
legend("topleft", c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"), col = 1:7, pch = 1:7, ncol = 2, cex = 0.8)
dev.off()

##############
# evaluation 
##############

ridge_update_test <- function(lambda_val, data_set, observed_period, fore_method)
{
    ridge_h1_err = matrix(NA, dim(data_set)[3], n_test)
    for(iw in 1:n_test)
    {
        n_train = (408 * 3/4 + iw - 1)
        ridge_h1_err[,iw] = ridge_update(lambda = lambda_val, data = data_set[,1:n_train,], 
                                     newdata =  data_set[1:observed_period,(n_train + 1),], 
                                     holdoutdata = data_set[(observed_period + 1):24,(n_train + 1),], 
                                     order = 6, fmethod = fore_method)
        rm(iw)
    }
    return(ridge_h1_err)
}

# ridge estimator (Mon - Sun)

ridge_log_Mon_list <- foreach(iwk = 2:23) %dopar% ridge_update_test(lambda_val = ridge_lambda_log_Mon_para[iwk-1], 
                                                                data_set = PNSD_log_Mon_408, observed_period = iwk)

ridge_log_Tue_list <- foreach(iwk = 2:23) %dopar% ridge_update_test(lambda_val = ridge_lambda_log_Tue_para[iwk-1], 
                                                                data_set = PNSD_log_Tue_408, observed_period = iwk)

ridge_log_Wed_list <- foreach(iwk = 2:23) %dopar% ridge_update_test(lambda_val = ridge_lambda_log_Wed_para[iwk-1], 
                                                                data_set = PNSD_log_Wed_408, observed_period = iwk)

ridge_log_Thu_list <- foreach(iwk = 2:23) %dopar% ridge_update_test(lambda_val = ridge_lambda_log_Thu_para[iwk-1], 
                                                                data_set = PNSD_log_Thu_408, observed_period = iwk)

ridge_log_Fri_list <- foreach(iwk = 2:23) %dopar% ridge_update_test(lambda_val = ridge_lambda_log_Fri_para[iwk-1], 
                                                                data_set = PNSD_log_Fri_408, observed_period = iwk)

ridge_log_Sat_list <- foreach(iwk = 2:23) %dopar% ridge_update_test(lambda_val = ridge_lambda_log_Sat_para[iwk-1], 
                                                                data_set = PNSD_log_Sat_408, observed_period = iwk)

ridge_log_Sun_list <- foreach(iwk = 2:23) %dopar% ridge_update_test(lambda_val = ridge_lambda_log_Sun_para[iwk-1], 
                                                                data_set = PNSD_log_Sun_408, observed_period = iwk)

ridge_log_Mon_array = ridge_log_Tue_array = ridge_log_Wed_array = ridge_log_Thu_array = 
ridge_log_Fri_array = ridge_log_Sat_array = ridge_log_Sun_array = array(NA, dim = c(22, 51, 102), dimnames = list(3:24, x_grid, 1:102))
for(iwk in 1:22)
{
    ridge_log_Mon_array[iwk,,] = ridge_log_Mon_list[[iwk]]   
    ridge_log_Tue_array[iwk,,] = ridge_log_Tue_list[[iwk]]   
    ridge_log_Wed_array[iwk,,] = ridge_log_Wed_list[[iwk]]   
    ridge_log_Thu_array[iwk,,] = ridge_log_Thu_list[[iwk]]   
    ridge_log_Fri_array[iwk,,] = ridge_log_Fri_list[[iwk]]   
    ridge_log_Sat_array[iwk,,] = ridge_log_Sat_list[[iwk]]   
    ridge_log_Sun_array[iwk,,] = ridge_log_Sun_list[[iwk]]   
    rm(iwk)
}

ridge_log_Mon_array_mean = apply(ridge_log_Mon_array, c(1, 2), mean)
ridge_log_Tue_array_mean = apply(ridge_log_Tue_array, c(1, 2), mean)
ridge_log_Wed_array_mean = apply(ridge_log_Wed_array, c(1, 2), mean)
ridge_log_Thu_array_mean = apply(ridge_log_Thu_array, c(1, 2), mean)
ridge_log_Fri_array_mean = apply(ridge_log_Fri_array, c(1, 2), mean)
ridge_log_Sat_array_mean = apply(ridge_log_Sat_array, c(1, 2), mean)
ridge_log_Sun_array_mean = apply(ridge_log_Sun_array, c(1, 2), mean)

# by hour

ridge_log_Mon_array_mean_hour = apply(ridge_log_Mon_array_mean, 1, mean)
ridge_log_Tue_array_mean_hour = apply(ridge_log_Tue_array_mean, 1, mean)
ridge_log_Wed_array_mean_hour = apply(ridge_log_Wed_array_mean, 1, mean)
ridge_log_Thu_array_mean_hour = apply(ridge_log_Thu_array_mean, 1, mean)
ridge_log_Fri_array_mean_hour = apply(ridge_log_Fri_array_mean, 1, mean)
ridge_log_Sat_array_mean_hour = apply(ridge_log_Sat_array_mean, 1, mean)
ridge_log_Sun_array_mean_hour = apply(ridge_log_Sun_array_mean, 1, mean)

ridge_days_array_mean_hour = rbind(ridge_log_Mon_array_mean_hour,
                                   ridge_log_Tue_array_mean_hour,
                                   ridge_log_Wed_array_mean_hour,
                                   ridge_log_Thu_array_mean_hour,
                                   ridge_log_Fri_array_mean_hour,
                                   ridge_log_Sat_array_mean_hour,
                                   ridge_log_Sun_array_mean_hour)
colnames(ridge_days_array_mean_hour) = 3:24
rownames(ridge_days_array_mean_hour) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  
# by size

ridge_log_Mon_array_mean_size = apply(ridge_log_Mon_array_mean, 2, mean)
ridge_log_Tue_array_mean_size = apply(ridge_log_Tue_array_mean, 2, mean)
ridge_log_Wed_array_mean_size = apply(ridge_log_Wed_array_mean, 2, mean)
ridge_log_Thu_array_mean_size = apply(ridge_log_Thu_array_mean, 2, mean)
ridge_log_Fri_array_mean_size = apply(ridge_log_Fri_array_mean, 2, mean)
ridge_log_Sat_array_mean_size = apply(ridge_log_Sat_array_mean, 2, mean)
ridge_log_Sun_array_mean_size = apply(ridge_log_Sun_array_mean, 2, mean)

ridge_days_array_mean_size = rbind(ridge_log_Mon_array_mean_size,
                                   ridge_log_Tue_array_mean_size,
                                   ridge_log_Wed_array_mean_size,
                                   ridge_log_Thu_array_mean_size,
                                   ridge_log_Fri_array_mean_size,
                                   ridge_log_Sat_array_mean_size,
                                   ridge_log_Sun_array_mean_size)
colnames(ridge_days_array_mean_size) = x_grid
rownames(ridge_days_array_mean_size) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

