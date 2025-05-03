#####################
# conformal approach
# ridge
#####################

# data_set: 24 x T x 51
# observed_period: observation period
# fore_method: forecasting method
# level_sig: 1 - nominal coverage probability

conformal_PI_ridge <- function(data_set, observed_period, fore_method, level_sig)
{
    n_hour = dim(data_set)[1]
    n_PM = dim(data_set)[3]
    
    ridge_h1_fore = array(NA, dim = c(n_PM, (n_hour - observed_period), n_test), 
                          dimnames = list(x_grid, (1:24)[-(1:observed_period)], 1:102))
    for(iw in 1:n_test)
    {
        n_train = (408 * 1/2 + iw - 1)
        ridge_h1_fore[,,iw] = ridge_update(lambda = lambda_val, data = data_set[,1:n_train,], 
                                           newdata =  data_set[1:observed_period,(n_train + 1),], 
                                           holdoutdata = data_set[(observed_period + 1):24,(n_train + 1),], 
                                           order = 6, fmethod = fore_method)$fore
        rm(iw); rm(n_train)
    }
    
    # holdout data in the original scale
    
    holdout_set = array((10^(aperm(data_set[,((408*1/2+1):(408*3/4)),], c(3, 1, 2))) - 1)[,(1:24)[-(1:observed_period)],],
                        dim = c(n_PM, length((1:24)[-(1:observed_period)]), n_test))
    
    # difference between holdout data and forecasts (51 x 24 x n_test)
    
    resi_mat = holdout_set - ridge_h1_fore
    
    # compute quantiles
    
    quantile_resid = matrix(NA, n_PM, (n_hour - observed_period))
    for(iw in 1:n_PM)
    {
        quantile_resid[iw,] = apply(matrix(resi_mat[iw,,], length((1:24)[-(1:observed_period)]), n_test), 1, function(x) quantile(abs(x), probs = (1 - level_sig)))
        rm(iw)
    }
    rm(holdout_set); rm(ridge_h1_fore)
    
    # holdout data in the original scale for the testing period
    
    holdout_set = array((10^(aperm(data_set[,((408*3/4+1):408),], c(3, 1, 2))) - 1)[,(1:24)[-(1:observed_period)],],
                        dim = c(n_PM, length((1:24)[-(1:observed_period)]), n_test))
    
    # forecast via ridge regression 
    
    ridge_h1_fore = array(NA, dim = c(n_PM, (n_hour - observed_period), n_test), 
                          dimnames = list(x_grid, (1:24)[-(1:observed_period)], 1:102))
    for(iw in 1:n_test)
    {
        n_train = (408 * 3/4 + iw - 1)
        ridge_h1_fore[,,iw] = ridge_update(lambda = lambda_val, data = data_set[,1:n_train,], 
                                           newdata =  data_set[1:observed_period,(n_train + 1),], 
                                           holdoutdata = data_set[(observed_period + 1):24,(n_train + 1),], 
                                           order = 6, fmethod = fore_method)$fore
        rm(iw); rm(n_train)
    }
    
    # compute the lower and upper bounds
    
    test_ridge_h1_lb = test_ridge_h1_ub = array(NA, dim = dim(ridge_h1_fore), dimnames = dimnames(ridge_h1_fore))
    for(iw in 1:n_PM)
    {
        test_ridge_h1_lb[iw,,] = ridge_h1_fore[iw,,] - quantile_resid[iw,]
        test_ridge_h1_ub[iw,,] = ridge_h1_fore[iw,,] + quantile_resid[iw,]
        rm(iw)
    }
    
    # compute empirical coverage probability (ECP), coverage probability difference (CPD), interval score
    
    int_val = matrix(NA, n_PM, 3) 
    for(iw in 1:n_PM)
    {
        int_val[iw,] = interval_score(holdout = holdout_set[iw,,], 
                                      lb = test_ridge_h1_lb[iw,,],
                                      ub = test_ridge_h1_ub[iw,,], 
                                      alpha = level_sig)
        rm(iw)
    }
    colnames(int_val) = c("ECP", "CPD", "score")
    rm(ridge_h1_fore); rm(test_ridge_h1_lb); rm(test_ridge_h1_ub); rm(holdout_set); rm(quantile_resid)
    return(int_val)   
}

##############################
# level of significance = 0.2
##############################

registerDoMC(detectCores())
PNSD_Mon_ridge_conformal = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Mon_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.2)

PNSD_Tue_ridge_conformal = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Tue_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.2)

PNSD_Wed_ridge_conformal = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Wed_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.2)

PNSD_Thu_ridge_conformal = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Thu_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.2)

PNSD_Fri_ridge_conformal = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Fri_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.2)

PNSD_Sat_ridge_conformal = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Sat_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.2)

PNSD_Sun_ridge_conformal = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Sun_408,
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.2)


PNSD_Mon_ridge_conformal_array = PNSD_Tue_ridge_conformal_array = PNSD_Wed_ridge_conformal_array = 
PNSD_Thu_ridge_conformal_array = PNSD_Fri_ridge_conformal_array = PNSD_Sat_ridge_conformal_array = 
PNSD_Sun_ridge_conformal_array = array(NA, dim = c(51, 3, 22), dimnames = list(x_grid, c("ECP", "CPD", "score"), 3:24))
for(iwk in 1:22)
{
    PNSD_Mon_ridge_conformal_array[,,iwk] = PNSD_Mon_ridge_conformal[[iwk]]
    PNSD_Tue_ridge_conformal_array[,,iwk] = PNSD_Tue_ridge_conformal[[iwk]]
    PNSD_Wed_ridge_conformal_array[,,iwk] = PNSD_Wed_ridge_conformal[[iwk]]
    PNSD_Thu_ridge_conformal_array[,,iwk] = PNSD_Thu_ridge_conformal[[iwk]]
    PNSD_Fri_ridge_conformal_array[,,iwk] = PNSD_Fri_ridge_conformal[[iwk]]
    PNSD_Sat_ridge_conformal_array[,,iwk] = PNSD_Sat_ridge_conformal[[iwk]]
    PNSD_Sun_ridge_conformal_array[,,iwk] = PNSD_Sun_ridge_conformal[[iwk]]
    print(iwk); rm(iwk)
}

PNSD_days_ridge_conformal_array_CPD = rbind(apply(PNSD_Mon_ridge_conformal_array[,2,], 2, mean),
                                            apply(PNSD_Tue_ridge_conformal_array[,2,], 2, mean),
                                            apply(PNSD_Wed_ridge_conformal_array[,2,], 2, mean),
                                            apply(PNSD_Thu_ridge_conformal_array[,2,], 2, mean),
                                            apply(PNSD_Fri_ridge_conformal_array[,2,], 2, mean),
                                            apply(PNSD_Sat_ridge_conformal_array[,2,], 2, mean),
                                            apply(PNSD_Sun_ridge_conformal_array[,2,], 2, mean))

PNSD_days_ridge_conformal_array_score = rbind(apply(PNSD_Mon_ridge_conformal_array[,3,], 2, mean),
                                              apply(PNSD_Tue_ridge_conformal_array[,3,], 2, mean),
                                              apply(PNSD_Wed_ridge_conformal_array[,3,], 2, mean),
                                              apply(PNSD_Thu_ridge_conformal_array[,3,], 2, mean),
                                              apply(PNSD_Fri_ridge_conformal_array[,3,], 2, mean),
                                              apply(PNSD_Sat_ridge_conformal_array[,3,], 2, mean),
                                              apply(PNSD_Sun_ridge_conformal_array[,3,], 2, mean))
rownames(PNSD_days_ridge_conformal_array_CPD) = rownames(PNSD_days_ridge_conformal_array_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(PNSD_days_ridge_conformal_array_CPD) = colnames(PNSD_days_ridge_conformal_array_score) = 3:24

###############################
# level of significance = 0.05
###############################

registerDoMC(detectCores())
PNSD_Mon_ridge_conformal_0.05 = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Mon_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.05)

PNSD_Tue_ridge_conformal_0.05 = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Tue_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.05)

PNSD_Wed_ridge_conformal_0.05 = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Wed_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.05)

PNSD_Thu_ridge_conformal_0.05 = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Thu_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.05)

PNSD_Fri_ridge_conformal_0.05 = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Fri_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.05)

PNSD_Sat_ridge_conformal_0.05 = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Sat_408, 
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.05)

PNSD_Sun_ridge_conformal_0.05 = foreach(iwk = 1:22) %dopar% conformal_PI_ridge(data_set = PNSD_log_Sun_408,
                                                                          observed_period = (iwk + 1), 
                                                                          fore_method = "ets", level_sig = 0.05)


PNSD_Mon_ridge_conformal_array_0.05 = PNSD_Tue_ridge_conformal_array_0.05 = PNSD_Wed_ridge_conformal_array_0.05 = 
PNSD_Thu_ridge_conformal_array_0.05 = PNSD_Fri_ridge_conformal_array_0.05 = PNSD_Sat_ridge_conformal_array_0.05 = 
PNSD_Sun_ridge_conformal_array_0.05 = array(NA, dim = c(51, 3, 22), dimnames = list(x_grid, c("ECP", "CPD", "score"), 3:24))
for(iwk in 1:22)
{
    PNSD_Mon_ridge_conformal_array_0.05[,,iwk] = PNSD_Mon_ridge_conformal_0.05[[iwk]]
    PNSD_Tue_ridge_conformal_array_0.05[,,iwk] = PNSD_Tue_ridge_conformal_0.05[[iwk]]
    PNSD_Wed_ridge_conformal_array_0.05[,,iwk] = PNSD_Wed_ridge_conformal_0.05[[iwk]]
    PNSD_Thu_ridge_conformal_array_0.05[,,iwk] = PNSD_Thu_ridge_conformal_0.05[[iwk]]
    PNSD_Fri_ridge_conformal_array_0.05[,,iwk] = PNSD_Fri_ridge_conformal_0.05[[iwk]]
    PNSD_Sat_ridge_conformal_array_0.05[,,iwk] = PNSD_Sat_ridge_conformal_0.05[[iwk]]
    PNSD_Sun_ridge_conformal_array_0.05[,,iwk] = PNSD_Sun_ridge_conformal_0.05[[iwk]]
    print(iwk); rm(iwk)
}

PNSD_days_ridge_conformal_array_CPD_0.05 = rbind(apply(PNSD_Mon_ridge_conformal_array_0.05[,2,], 2, mean),
                                                 apply(PNSD_Tue_ridge_conformal_array_0.05[,2,], 2, mean),
                                                 apply(PNSD_Wed_ridge_conformal_array_0.05[,2,], 2, mean),
                                                 apply(PNSD_Thu_ridge_conformal_array_0.05[,2,], 2, mean),
                                                 apply(PNSD_Fri_ridge_conformal_array_0.05[,2,], 2, mean),
                                                 apply(PNSD_Sat_ridge_conformal_array_0.05[,2,], 2, mean),
                                                 apply(PNSD_Sun_ridge_conformal_array_0.05[,2,], 2, mean))

PNSD_days_ridge_conformal_array_score_0.05 = rbind(apply(PNSD_Mon_ridge_conformal_array_0.05[,3,], 2, mean),
                                                   apply(PNSD_Tue_ridge_conformal_array_0.05[,3,], 2, mean),
                                                   apply(PNSD_Wed_ridge_conformal_array_0.05[,3,], 2, mean),
                                                   apply(PNSD_Thu_ridge_conformal_array_0.05[,3,], 2, mean),
                                                   apply(PNSD_Fri_ridge_conformal_array_0.05[,3,], 2, mean),
                                                   apply(PNSD_Sat_ridge_conformal_array_0.05[,3,], 2, mean),
                                                   apply(PNSD_Sun_ridge_conformal_array_0.05[,3,], 2, mean))

rownames(PNSD_days_ridge_conformal_array_CPD_0.05) = rownames(PNSD_days_ridge_conformal_array_score_0.05) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(PNSD_days_ridge_conformal_array_CPD_0.05) = colnames(PNSD_days_ridge_conformal_array_score_0.05) = 3:24

