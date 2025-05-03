##########################################
# Multilevel functional time series model
##########################################

# data_set: a list of p by n data matrix
# aux_var: an aggregated p by n data matrix
# ncomp_method: method for selecting the number of components
# fh: forecast horizon
# fore_method: univariate time-series forecasting method, such as exponential smoothing "ETS"

MLFTS_model <- function(data_input, aux_var, ncomp_method, fh, fore_method)
{
    n_age  = dim(data_input)[1]
    n_year = dim(data_input)[2]
    n_pop  = dim(data_input)[3]
    
    # clean the data
    
    data_set_array = array(NA, dim = c(n_age, n_year, n_pop))
    if(any(!is.finite(data_input)))
    {
        for(iw in 1:n_pop)
        {
            for(ij in 1:n_age)
            {
                data_set_array[ij,,iw] = na.interp(data_input[ij,,iw])
            }
        }
    }
    else
    {
        data_set_array = data_input
    }
    
    # compute the mean function
    
    mean_function_list = list()
    for(ik in 1:n_pop)
    {
        mean_function_list[[ik]] = rowMeans(data_set_array[,,ik], na.rm = TRUE)
        rm(ik)
    }
    
    data_set = array(NA, dim = c(n_age, n_year, n_pop))
    for(ik in 1:n_pop)
    {
        data_set[,,ik] = t(scale(t(data_set_array[,,ik]), center = TRUE, scale = FALSE))
        rm(ik)
    }
    
    if(missing(aux_var)|is.null(aux_var))
    {
        aggregate_data = apply(data_set, c(1, 2), mean)
    }
    else
    {
        aggregate_data = t(aux_var)
    }
    colnames(aggregate_data) = 1:n_year
    rownames(aggregate_data) = 1:n_age
    
    # 1st FPCA
    
    eigen_value_aggregate = eigen(cov(t(aggregate_data)))$values
    if(ncomp_method == "EVR")
    {
        ncomp_aggregate = select_K(tau = 10^-3, eigenvalue = eigen_value_aggregate)
    }
    else if(ncomp_method == "provide")
    {
        ncomp_aggregate = 6
    }
    ftsm_aggregate = ftsm(y = fts(1:n_age, aggregate_data), order = ncomp_aggregate)
    
    # calculate sum of lambda_k
    sum_lambda_k = sum(eigen_value_aggregate[1:ncomp_aggregate])
    
    # compute the residual trend
    data_residual = array(NA, dim = c(n_age, n_year, n_pop))
    for(iw in 1:n_pop)
    {
        data_residual[,,iw] = data_set[,,iw] - ftsm_aggregate$fitted$y
        colnames(data_residual[,,iw]) = 1:n_year
        rownames(data_residual[,,iw]) = 1:n_age
        rm(iw)
    }
    
    # 2nd FPCA
    
    if(ncomp_method == "EVR")
    {
        ncomp_resi = vector("numeric", n_pop)
        for(iw in 1:n_pop)
        {
            eigen_value_resi = eigen(cov(t(data_residual[,,iw])))$values
            ncomp_resi[iw] = select_K(tau = 10^-3, eigenvalue = eigen_value_resi)
        }
    }
    else if(ncomp_method == "provide")
    {
        ncomp_resi = rep(6, n_pop)
    }
    
    sum_lambda_l = vector("numeric", n_pop)
    for(iw in 1:n_pop)
    {
        eigen_value_resi = eigen(cov(t(data_residual[,,iw])))$values
        
        # calculate sum of lambda_l
        sum_lambda_l[iw] = sum(eigen_value_resi[1:(ncomp_resi[iw])])
    }
    
    ftsm_resi = list()
    for(iw in 1:n_pop)
    {
        y = fts(1:n_age, data_residual[,,iw])
        benchmark_ngrid = 500
        dum = try(ftsm(y, order = ncomp_resi[iw], ngrid = max(benchmark_ngrid, ncol(y$y))), silent = TRUE)
        while(any(class(dum) == "try-error"))
        {
            dum = try(ftsm(y, order = ncomp_resi[iw], ngrid = max(benchmark_ngrid + floor(runif(1, -150, 150)), ncol(y$y))), 
                      silent = TRUE)
        }
        ftsm_resi[[iw]] = dum
        rm(iw); rm(dum)
    }
    
    # within-cluster variability
    
    within_cluster_variability = vector("numeric", n_pop)
    for(iw in 1:n_pop)
    {
        within_cluster_variability[iw] = sum_lambda_k/(sum_lambda_k + sum_lambda_l[iw])
    }
    
    # reconstruction
    
    coef_fore = matrix(NA, ncomp_aggregate, fh)
    if(fore_method == "arima")
    {
        for(ik in 1:ncomp_aggregate)
        {
            coef_fore[ik,] = forecast(auto.arima(ftsm_aggregate$coeff[,ik+1]), h = fh)$mean
        }
    }
    else if(fore_method == "ets")
    {
        for(ik in 1:ncomp_aggregate)
        {
            coef_fore[ik,] = forecast(ets(ftsm_aggregate$coeff[,ik+1]), h = fh)$mean
        }
    }
    else
    {
        warning("Forecasting method can either be ARIMA or ETS.")
    }
    rownames(coef_fore) = 1:ncomp_aggregate
    colnames(coef_fore) = 1:fh
    
    if(ncomp_aggregate == 1)
    {
        aggregate_fore = as.matrix(ftsm_aggregate$basis[,2]) %*% matrix(coef_fore, nrow = 1)
    }
    else
    {
        aggregate_fore = ftsm_aggregate$basis[,2:(ncomp_aggregate+1)] %*% coef_fore
    }
    
    # residual forecasts
    
    coef_fore_resi_list = list()
    for(iw in 1:n_pop)
    {
        coef_fore_resi = matrix(NA, ncomp_resi[iw], fh)
        if(fore_method == "arima")
        {
            for(ik in 1:ncomp_resi[iw])
            {
                coef_fore_resi[ik,] = forecast(auto.arima(ftsm_resi[[iw]]$coeff[,ik+1]), h = fh)$mean
            }
        }
        else if(fore_method == "ets")
        {
            for(ik in 1:ncomp_resi[iw])
            {
                coef_fore_resi[ik,] = forecast(ets(ftsm_resi[[iw]]$coeff[,ik+1]), h = fh)$mean
            }
        }
        else
        {
            warning("Forecasting method can either be ARIMA or ETS.")
        }
        coef_fore_resi_list[[iw]] = coef_fore_resi
        rm(iw)
    }
    
    resi_fore = list()
    for(iw in 1:n_pop)
    {
        resi_fore[[iw]] = ftsm_resi[[iw]]$basis[,2:(ncomp_resi[iw] + 1)] %*% coef_fore_resi_list[[iw]]
        rm(iw)
    }
    
    final_fore = list()
    for(iw in 1:n_pop)
    {
        final_fore[[iw]] = mean_function_list[[iw]] + (aggregate_fore + resi_fore[[iw]])[,fh]
        rm(iw)
    }
    return(final_fore)
}

########################
# expanding window
# h = 1 (one-day-ahead)
########################

n_test = dim(PNSD_Mon_test)[2]

MLFTS_fit_h1_Mon = MLFTS_fit_h1_Tue = MLFTS_fit_h1_Wed = MLFTS_fit_h1_Thu = MLFTS_fit_h1_Fri = 
MLFTS_fit_h1_Sat = MLFTS_fit_h1_Sun = array(NA, dim = c(51, 24, n_test), 
                                            dimnames = list(x_grid, 1:24, 1:n_test))
for(iw in 1:n_test)
{
    n_train = (408 * 3/4 + iw - 1)
    MLFTS_fit_h1_Mon[,,iw] <- matrix(unlist(MLFTS_model(data_input = aperm(PNSD_log_Mon_408[,1:n_train,], c(3, 2, 1)),
                                                        aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
    
    MLFTS_fit_h1_Tue[,,iw] <- matrix(unlist(MLFTS_model(data_input = aperm(PNSD_log_Tue_408[,1:n_train,], c(3, 2, 1)),
                                                        aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
    
    MLFTS_fit_h1_Wed[,,iw] <- matrix(unlist(MLFTS_model(data_input = aperm(PNSD_log_Wed_408[,1:n_train,], c(3, 2, 1)),
                                                        aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
    
    MLFTS_fit_h1_Thu[,,iw] <- matrix(unlist(MLFTS_model(data_input = aperm(PNSD_log_Thu_408[,1:n_train,], c(3, 2, 1)),
                                                        aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
    
    MLFTS_fit_h1_Fri[,,iw] <- matrix(unlist(MLFTS_model(data_input = aperm(PNSD_log_Fri_408[,1:n_train,], c(3, 2, 1)),
                                                        aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
    
    MLFTS_fit_h1_Sat[,,iw] <- matrix(unlist(MLFTS_model(data_input = aperm(PNSD_log_Sat_408[,1:n_train,], c(3, 2, 1)),
                                                        aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
    
    MLFTS_fit_h1_Sun[,,iw] <- matrix(unlist(MLFTS_model(data_input = aperm(PNSD_log_Sun_408[,1:n_train,], c(3, 2, 1)),
                                                        aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
    print(n_train); rm(iw); rm(n_train)
}

# take the 10 power transformation back to the original scale

MLFTS_fore_h1_Mon = 10^(MLFTS_fit_h1_Mon) - 1
MLFTS_fore_h1_Tue = 10^(MLFTS_fit_h1_Tue) - 1
MLFTS_fore_h1_Wed = 10^(MLFTS_fit_h1_Wed) - 1
MLFTS_fore_h1_Thu = 10^(MLFTS_fit_h1_Thu) - 1
MLFTS_fore_h1_Fri = 10^(MLFTS_fit_h1_Fri) - 1
MLFTS_fore_h1_Sat = 10^(MLFTS_fit_h1_Sat) - 1
MLFTS_fore_h1_Sun = 10^(MLFTS_fit_h1_Sun) - 1

#################
# compute errors
#################

# MAPE

mape <- ftsa:::mape

# MAPE aggregated across 51 sizes

mape_h1_Mon_TS = mape_h1_Tue_TS = mape_h1_Wed_TS = mape_h1_Thu_TS = 
mape_h1_Fri_TS = mape_h1_Sat_TS = mape_h1_Sun_TS = matrix(NA, 24, n_test)
for(ik in 1:n_test)
{
    for(ij in 1:24)
    {
        mape_h1_Mon_TS[ij,ik] = mape(forecast = MLFTS_fore_h1_Mon[,ij,ik], true = PNSD_Mon_test[ij,ik,])
        mape_h1_Tue_TS[ij,ik] = mape(forecast = MLFTS_fore_h1_Tue[,ij,ik], true = PNSD_Tue_test[ij,ik,])
        mape_h1_Wed_TS[ij,ik] = mape(forecast = MLFTS_fore_h1_Wed[,ij,ik], true = PNSD_Wed_test[ij,ik,])
        mape_h1_Thu_TS[ij,ik] = mape(forecast = MLFTS_fore_h1_Thu[,ij,ik], true = PNSD_Thu_test[ij,ik,])
        mape_h1_Fri_TS[ij,ik] = mape(forecast = MLFTS_fore_h1_Fri[,ij,ik], true = PNSD_Fri_test[ij,ik,])
        mape_h1_Sat_TS[ij,ik] = mape(forecast = MLFTS_fore_h1_Sat[,ij,ik], true = PNSD_Sat_test[ij,ik,])
        mape_h1_Sun_TS[ij,ik] = mape(forecast = MLFTS_fore_h1_Sun[,ij,ik], true = PNSD_Sun_test[ij,ik,])
    }
}

mape_h1_Mon_TS_mean = apply(mape_h1_Mon_TS, 1, mean)
mape_h1_Tue_TS_mean = apply(mape_h1_Tue_TS, 1, mean)
mape_h1_Wed_TS_mean = apply(mape_h1_Wed_TS, 1, mean)
mape_h1_Thu_TS_mean = apply(mape_h1_Thu_TS, 1, mean)
mape_h1_Fri_TS_mean = apply(mape_h1_Fri_TS, 1, mean)
mape_h1_Sat_TS_mean = apply(mape_h1_Sat_TS, 1, mean)
mape_h1_Sun_TS_mean = apply(mape_h1_Sun_TS, 1, mean)

mape_h1_days_TS_mean = rbind(mape_h1_Mon_TS_mean, mape_h1_Tue_TS_mean, mape_h1_Wed_TS_mean,
                             mape_h1_Thu_TS_mean, mape_h1_Fri_TS_mean, mape_h1_Sat_TS_mean,
                             mape_h1_Sun_TS_mean)
rownames(mape_h1_days_TS_mean) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(mape_h1_days_TS_mean) = 1:24

# boxplot

boxplot(t(mape_h1_Mon_TS), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Mon_factor), col = rgb(1, 0, 0, 0.4), outline = FALSE, add = TRUE, boxwex = 0.5)
legend("topright", legend = c("MLFTS", "Factor model + MLFTS"),
       fill = c("gray", rgb(1, 0, 0, 0.4)), cex = 0.8)

boxplot(t(mape_h1_Tue_TS), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Tue_factor), col = rgb(1, 0, 0, 0.4), outline = FALSE, add = TRUE, boxwex = 0.5)

boxplot(t(mape_h1_Wed_TS), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Wed_factor), col = rgb(1, 0, 0, 0.4), outline = FALSE, add = TRUE, boxwex = 0.5)

boxplot(t(mape_h1_Thu_TS), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Fri_TS), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Sat_TS), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Sun_TS), notch = TRUE, outline = FALSE, ylab = "MAPE")

plot(1:24, mape_h1_Mon_TS_mean, ylab = "MAPE", xlab = "Hour", type = "l", ylim = c(50, 350))
lines(1:24, mape_h1_Tue_TS_mean, col = 2, lty = 2)
lines(1:24, mape_h1_Wed_TS_mean, col = 3, lty = 3)
lines(1:24, mape_h1_Thu_TS_mean, col = 4, lty = 4)
lines(1:24, mape_h1_Fri_TS_mean, col = 5, lty = 5)
lines(1:24, mape_h1_Sat_TS_mean, col = 6, lty = 6)
lines(1:24, mape_h1_Sun_TS_mean, col = 7, lty = 7)
legend("topright", c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"), col = 1:7, lty = 1:7, cex = 0.8, ncol = 2)


##################################
# MAPE aggregated across 24 hours
##################################

mape_h1_Mon_TS_size = mape_h1_Tue_TS_size = mape_h1_Wed_TS_size = mape_h1_Thu_TS_size = 
mape_h1_Fri_TS_size = mape_h1_Sat_TS_size = mape_h1_Sun_TS_size = matrix(NA, 51, n_test)
for(ik in 1:n_test)
{
    for(ij in 1:51)
    {
        mape_h1_Mon_TS_size[ij,ik] = mape(forecast = MLFTS_fore_h1_Mon[ij,,ik], true = PNSD_Mon_test[,ik,ij])
        mape_h1_Tue_TS_size[ij,ik] = mape(forecast = MLFTS_fore_h1_Tue[ij,,ik], true = PNSD_Tue_test[,ik,ij])
        mape_h1_Wed_TS_size[ij,ik] = mape(forecast = MLFTS_fore_h1_Wed[ij,,ik], true = PNSD_Wed_test[,ik,ij])
        mape_h1_Thu_TS_size[ij,ik] = mape(forecast = MLFTS_fore_h1_Thu[ij,,ik], true = PNSD_Thu_test[,ik,ij])
        mape_h1_Fri_TS_size[ij,ik] = mape(forecast = MLFTS_fore_h1_Fri[ij,,ik], true = PNSD_Fri_test[,ik,ij])
        mape_h1_Sat_TS_size[ij,ik] = mape(forecast = MLFTS_fore_h1_Sat[ij,,ik], true = PNSD_Sat_test[,ik,ij])
        mape_h1_Sun_TS_size[ij,ik] = mape(forecast = MLFTS_fore_h1_Sun[ij,,ik], true = PNSD_Sun_test[,ik,ij])
    }
}

mape_h1_Mon_TS_size_mean = apply(mape_h1_Mon_TS_size, 1, mean)
mape_h1_Tue_TS_size_mean = apply(mape_h1_Tue_TS_size, 1, mean)
mape_h1_Wed_TS_size_mean = apply(mape_h1_Wed_TS_size, 1, mean)
mape_h1_Thu_TS_size_mean = apply(mape_h1_Thu_TS_size, 1, mean)
mape_h1_Fri_TS_size_mean = apply(mape_h1_Fri_TS_size, 1, mean)
mape_h1_Sat_TS_size_mean = apply(mape_h1_Sat_TS_size, 1, mean)
mape_h1_Sun_TS_size_mean = apply(mape_h1_Sun_TS_size, 1, mean)


mape_h1_days_TS_size_mean = rbind(mape_h1_Mon_TS_size_mean,
                                  mape_h1_Tue_TS_size_mean,
                                  mape_h1_Wed_TS_size_mean,
                                  mape_h1_Thu_TS_size_mean,
                                  mape_h1_Fri_TS_size_mean,
                                  mape_h1_Sat_TS_size_mean,
                                  mape_h1_Sun_TS_size_mean)
colnames(mape_h1_days_TS_size_mean) = x_grid
rownames(mape_h1_days_TS_size_mean) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

boxplot(t(mape_h1_Mon_TS_size), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Tue_TS_size), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Wed_TS_size), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Thu_TS_size), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Fri_TS_size), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Sat_TS_size), notch = TRUE, outline = FALSE, ylab = "MAPE")
boxplot(t(mape_h1_Sun_TS_size), notch = TRUE, outline = FALSE, ylab = "MAPE")

plot(x_grid, mape_h1_Mon_TS_size_mean, ylab = "MAPE", xlab = "Hour", type = "l", ylim = c(50, 450))
lines(x_grid, mape_h1_Tue_TS_size_mean, col = 2, lty = 2)
lines(x_grid, mape_h1_Wed_TS_size_mean, col = 3, lty = 3)
lines(x_grid, mape_h1_Thu_TS_size_mean, col = 4, lty = 4)
lines(x_grid, mape_h1_Fri_TS_size_mean, col = 5, lty = 5)
lines(x_grid, mape_h1_Sat_TS_size_mean, col = 6, lty = 6)
lines(x_grid, mape_h1_Sun_TS_size_mean, col = 7, lty = 7)

