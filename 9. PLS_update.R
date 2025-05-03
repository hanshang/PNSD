# lambda: tuning parameter
# data: data_array 24 x T x 51
# newdata: newly arrived data 
# holdoutdata: holdout data
# order: number of retained functional principal components
# fmethod: ARIMA or ETS

PLS_update <- function(lambda, data, newdata, holdoutdata, order, fmethod)
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
        base_array[,ik,] = base[((ik - 1) * n_hour + 1):(ik * n_hour),]
        rm(ik)
    }
    
    n2 = dim(newdata)[1]
    base1 = base_array[1:n2,,]
    base2 = array(base_array[(n2 + 1):n_hour,,], dim = c(length((n2 + 1):n_hour), n_PM, (order + 1)))
    
    fore = matrix(NA, (order + 1), 1)
    if(fmethod == "arima") 
    {
        for(ik in 1:(order + 1)) 
        {
            fore[ik, ] = forecast(auto.arima(coef[, ik]), h = 1)$mean
            rm(ik)
        }
    }
    else if(fmethod == "ets") 
    {
        for(ik in 1:(order + 1)) 
        {
            fore[ik, ] = forecast(ets(coef[, ik]), h = 1)$mean
            rm(ik)
        }
    }
    
    pls_mat = matrix(NA, (order + 1), n_PM)
    for(ik in 1:n_PM)
    {
        I = diag(dim(base)[2])
        pls_mat[,ik] = ginv(t(base1[,ik,]) %*% base1[,ik,] + lambda * I) %*% (t(base1[,ik,]) %*% newdata[,ik] + lambda * fore)
        rm(ik); rm(I)
    }   
    
    update_period = (nrow(newdata) + 1):n_hour
    pls_forecasts = matrix(NA, length(update_period), n_PM)
    for(ik in 1:n_PM)
    {
        pls_forecasts[,ik] = base2[,ik,] %*% matrix(pls_mat[,ik], (order + 1), 1)
        rm(ik)
    }
    
    err = vector("numeric", n_PM)
    for(ik in 1:n_PM)
    {
        err[ik] = ftsa:::mape(forecast = (10^(pls_forecasts[,ik]) - 1), true = (10^(matrix(holdoutdata, length(update_period), n_PM)[,ik]) - 1))
        rm(ik)
    }
    rm(data)
    return(err)
}

# validation set: 205:306

PLS_update_validation <- function(lambda_val, data_set, observed_period, fore_method)
{
    PLS_h1_err = matrix(NA, dim(data_set)[3], n_test)
    for(iw in 1:n_test)
    {
        n_train = (408 * 1/2 + iw - 1)
        PLS_h1_err[,iw] = PLS_update(lambda = lambda_val, data = data_set[,1:n_train,], 
                         newdata =  data_set[1:observed_period,(n_train + 1),], 
                         holdoutdata = data_set[(observed_period + 1):24,(n_train + 1),], 
                         order = 6, fmethod = fore_method)
        rm(iw); rm(n_train)
    }
    return(round(mean(PLS_h1_err), 4))
}

# PLS forecast (Mon - Sun)

registerDoMC(22)
PLS_lambda_log_Mon_list = foreach(iwk = 2:23) %dopar% optimise(f = PLS_update_validation, 
          interval = c(0, 10^3), data_set = PNSD_log_Mon_408, observed_period = iwk, fore_method = "arima")

PLS_lambda_log_Tue_list = foreach(iwk = 2:23) %dopar% optimise(f = PLS_update_validation, 
          interval = c(0, 10^3), data_set = PNSD_log_Tue_408, observed_period = iwk, fore_method = "arima")

PLS_lambda_log_Wed_list = foreach(iwk = 2:23) %dopar% optimise(f = PLS_update_validation, 
          interval = c(0, 10^3), data_set = PNSD_log_Wed_408, observed_period = iwk, fore_method = "arima")

PLS_lambda_log_Thu_list = foreach(iwk = 2:23) %dopar% optimise(f = PLS_update_validation, 
          interval = c(0, 10^3), data_set = PNSD_log_Thu_408, observed_period = iwk, fore_method = "arima")

PLS_lambda_log_Fri_list = foreach(iwk = 2:23) %dopar% optimise(f = PLS_update_validation, 
          interval = c(0, 10^3), data_set = PNSD_log_Fri_408, observed_period = iwk, fore_method = "arima")

PLS_lambda_log_Sat_list = foreach(iwk = 2:23) %dopar% optimise(f = PLS_update_validation, 
          interval = c(0, 10^3), data_set = PNSD_log_Sat_408, observed_period = iwk, fore_method = "arima")

PLS_lambda_log_Sun_list = foreach(iwk = 2:23) %dopar% optimise(f = PLS_update_validation, 
          interval = c(0, 10^3), data_set = PNSD_log_Sun_408, observed_period = iwk, fore_method = "arima")

PLS_lambda_log_Mon_para = PLS_lambda_log_Tue_para = PLS_lambda_log_Wed_para = 
PLS_lambda_log_Thu_para = PLS_lambda_log_Fri_para = PLS_lambda_log_Sat_para = PLS_lambda_log_Sun_para = 
PLS_lambda_log_Mon_obj = PLS_lambda_log_Tue_obj = PLS_lambda_log_Wed_obj = 
PLS_lambda_log_Thu_obj = PLS_lambda_log_Fri_obj = PLS_lambda_log_Sat_obj = PLS_lambda_log_Sun_obj = vector("numeric", 22)
for(iwk in 1:22)
{
    PLS_lambda_log_Mon_para[iwk] = PLS_lambda_log_Mon_list[[iwk]]$minimum
    PLS_lambda_log_Tue_para[iwk] = PLS_lambda_log_Tue_list[[iwk]]$minimum
    PLS_lambda_log_Wed_para[iwk] = PLS_lambda_log_Wed_list[[iwk]]$minimum
    PLS_lambda_log_Thu_para[iwk] = PLS_lambda_log_Thu_list[[iwk]]$minimum
    PLS_lambda_log_Fri_para[iwk] = PLS_lambda_log_Fri_list[[iwk]]$minimum
    PLS_lambda_log_Sat_para[iwk] = PLS_lambda_log_Sat_list[[iwk]]$minimum
    PLS_lambda_log_Sun_para[iwk] = PLS_lambda_log_Sun_list[[iwk]]$minimum
    
    PLS_lambda_log_Mon_obj[iwk] = PLS_lambda_log_Mon_list[[iwk]]$objective
    PLS_lambda_log_Tue_obj[iwk] = PLS_lambda_log_Tue_list[[iwk]]$objective
    PLS_lambda_log_Wed_obj[iwk] = PLS_lambda_log_Wed_list[[iwk]]$objective
    PLS_lambda_log_Thu_obj[iwk] = PLS_lambda_log_Thu_list[[iwk]]$objective
    PLS_lambda_log_Fri_obj[iwk] = PLS_lambda_log_Fri_list[[iwk]]$objective
    PLS_lambda_log_Sat_obj[iwk] = PLS_lambda_log_Sat_list[[iwk]]$objective
    PLS_lambda_log_Sun_obj[iwk] = PLS_lambda_log_Sun_list[[iwk]]$objective
    rm(iwk)
}

PLS_lambda_days_para = cbind(PLS_lambda_log_Mon_para, PLS_lambda_log_Tue_para, PLS_lambda_log_Wed_para, 
                             PLS_lambda_log_Thu_para, PLS_lambda_log_Fri_para, PLS_lambda_log_Sat_para,
                             PLS_lambda_log_Sun_para)
rownames(PLS_lambda_days_para) = 3:24
colnames(PLS_lambda_days_para) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  

savefig("Fig_3b", width = 12, height = 10, toplines = 0.8, type = "png")
plot(3:24, PLS_lambda_days_para[,1], col = 1, lty = 1, pch = 1, ylim = c(0, 162), xlab = "Intraday period", 
     ylab = "", main = "PLS estimator")
points(3:24, PLS_lambda_days_para[,2], col = 2, lty = 2, pch = 2)
points(3:24, PLS_lambda_days_para[,3], col = 3, lty = 3, pch = 3)
points(3:24, PLS_lambda_days_para[,4], col = 4, lty = 4, pch = 4)
points(3:24, PLS_lambda_days_para[,5], col = 5, lty = 5, pch = 5)
points(3:24, PLS_lambda_days_para[,6], col = 6, lty = 6, pch = 6)
points(3:24, PLS_lambda_days_para[,7], col = 7, lty = 7, pch = 7)
dev.off()

##############
# evaluation 
##############

PLS_update_test <- function(lambda_val, data_set, observed_period, fore_method)
{
    PLS_h1_err = matrix(NA, dim(data_set)[3], n_test)
    for(iw in 1:n_test)
    {
        n_train = (408 * 3/4 + iw - 1)
        PLS_h1_err[,iw] = PLS_update(lambda = lambda_val, data = data_set[,1:n_train,], 
                                     newdata =  data_set[1:observed_period,(n_train + 1),], 
                                     holdoutdata = data_set[(observed_period + 1):24,(n_train + 1),], 
                                     order = 6, fmethod = fore_method)
        rm(iw)
    }
    return(PLS_h1_err)
}

# PLS forecasts (Mon - Sun)

PLS_log_Mon_list <- foreach(iwk = 2:23) %dopar% PLS_update_test(lambda_val = PLS_lambda_log_Mon_para[iwk-1], 
                                                                data_set = PNSD_log_Mon_408, observed_period = iwk, fore_method = "arima")

PLS_log_Tue_list <- foreach(iwk = 2:23) %dopar% PLS_update_test(lambda_val = PLS_lambda_log_Tue_para[iwk-1], 
                                                                data_set = PNSD_log_Tue_408, observed_period = iwk, fore_method = "arima")

PLS_log_Wed_list <- foreach(iwk = 2:23) %dopar% PLS_update_test(lambda_val = PLS_lambda_log_Wed_para[iwk-1], 
                                                                data_set = PNSD_log_Wed_408, observed_period = iwk, fore_method = "arima")

PLS_log_Thu_list <- foreach(iwk = 2:23) %dopar% PLS_update_test(lambda_val = PLS_lambda_log_Thu_para[iwk-1], 
                                                                data_set = PNSD_log_Thu_408, observed_period = iwk, fore_method = "arima")

PLS_log_Fri_list <- foreach(iwk = 2:23) %dopar% PLS_update_test(lambda_val = PLS_lambda_log_Fri_para[iwk-1], 
                                                                data_set = PNSD_log_Fri_408, observed_period = iwk, fore_method = "arima")

PLS_log_Sat_list <- foreach(iwk = 2:23) %dopar% PLS_update_test(lambda_val = PLS_lambda_log_Sat_para[iwk-1], 
                                                                data_set = PNSD_log_Sat_408, observed_period = iwk, fore_method = "arima")

PLS_log_Sun_list <- foreach(iwk = 2:23) %dopar% PLS_update_test(lambda_val = PLS_lambda_log_Sun_para[iwk-1], 
                                                                data_set = PNSD_log_Sun_408, observed_period = iwk, fore_method = "arima")


PLS_log_Mon_array = PLS_log_Tue_array = PLS_log_Wed_array = PLS_log_Thu_array = 
PLS_log_Fri_array = PLS_log_Sat_array = PLS_log_Sun_array = array(NA, dim = c(22, 51, 102), dimnames = list(3:24, x_grid, 1:102))
for(iwk in 1:22)
{
    PLS_log_Mon_array[iwk,,] = PLS_log_Mon_list[[iwk]]   
    PLS_log_Tue_array[iwk,,] = PLS_log_Tue_list[[iwk]]   
    PLS_log_Wed_array[iwk,,] = PLS_log_Wed_list[[iwk]]   
    PLS_log_Thu_array[iwk,,] = PLS_log_Thu_list[[iwk]]   
    PLS_log_Fri_array[iwk,,] = PLS_log_Fri_list[[iwk]]   
    PLS_log_Sat_array[iwk,,] = PLS_log_Sat_list[[iwk]]   
    PLS_log_Sun_array[iwk,,] = PLS_log_Sun_list[[iwk]]   
    rm(iwk)
}

PLS_log_Mon_array_mean = apply(PLS_log_Mon_array, c(1, 2), mean)
PLS_log_Tue_array_mean = apply(PLS_log_Tue_array, c(1, 2), mean)
PLS_log_Wed_array_mean = apply(PLS_log_Wed_array, c(1, 2), mean)
PLS_log_Thu_array_mean = apply(PLS_log_Thu_array, c(1, 2), mean)
PLS_log_Fri_array_mean = apply(PLS_log_Fri_array, c(1, 2), mean)
PLS_log_Sat_array_mean = apply(PLS_log_Sat_array, c(1, 2), mean)
PLS_log_Sun_array_mean = apply(PLS_log_Sun_array, c(1, 2), mean)

# by hour

PLS_log_Mon_array_mean_hour = apply(PLS_log_Mon_array_mean, 1, mean)
PLS_log_Tue_array_mean_hour = apply(PLS_log_Tue_array_mean, 1, mean)
PLS_log_Wed_array_mean_hour = apply(PLS_log_Wed_array_mean, 1, mean)
PLS_log_Thu_array_mean_hour = apply(PLS_log_Thu_array_mean, 1, mean)
PLS_log_Fri_array_mean_hour = apply(PLS_log_Fri_array_mean, 1, mean)
PLS_log_Sat_array_mean_hour = apply(PLS_log_Sat_array_mean, 1, mean)
PLS_log_Sun_array_mean_hour = apply(PLS_log_Sun_array_mean, 1, mean)

PLS_days_array_mean_hour = rbind(PLS_log_Mon_array_mean_hour,
                                 PLS_log_Tue_array_mean_hour,
                                 PLS_log_Wed_array_mean_hour,
                                 PLS_log_Thu_array_mean_hour,
                                 PLS_log_Fri_array_mean_hour,
                                 PLS_log_Sat_array_mean_hour,
                                 PLS_log_Sun_array_mean_hour)
colnames(PLS_days_array_mean_hour) = 3:24
rownames(PLS_days_array_mean_hour) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  
# by size

PLS_log_Mon_array_mean_size = apply(PLS_log_Mon_array_mean, 2, mean)
PLS_log_Tue_array_mean_size = apply(PLS_log_Tue_array_mean, 2, mean)
PLS_log_Wed_array_mean_size = apply(PLS_log_Wed_array_mean, 2, mean)
PLS_log_Thu_array_mean_size = apply(PLS_log_Thu_array_mean, 2, mean)
PLS_log_Fri_array_mean_size = apply(PLS_log_Fri_array_mean, 2, mean)
PLS_log_Sat_array_mean_size = apply(PLS_log_Sat_array_mean, 2, mean)
PLS_log_Sun_array_mean_size = apply(PLS_log_Sun_array_mean, 2, mean)

PLS_days_array_mean_size = rbind(PLS_log_Mon_array_mean_size,
                                 PLS_log_Tue_array_mean_size,
                                 PLS_log_Wed_array_mean_size,
                                 PLS_log_Thu_array_mean_size,
                                 PLS_log_Fri_array_mean_size,
                                 PLS_log_Sat_array_mean_size,
                                 PLS_log_Sun_array_mean_size)
colnames(PLS_days_array_mean_size) = x_grid
rownames(PLS_days_array_mean_size) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

