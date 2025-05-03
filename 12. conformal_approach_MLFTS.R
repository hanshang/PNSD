source("load_packages.R")

# data_set: 51 x 24 x n_test
# fore_model: forecasting model
# level_sig: level of significance == (1 - nominal coverage probability)

conformal_PI <- function(data_set, fore_model, level_sig)
{
    n_hour = dim(data_set)[1]
    n_PM = dim(data_set)[3]
    
    # holdout data in the original scale for the validation period (205:306)
    
    holdout_set = 10^(aperm(data_set[,((408*1/2+1):(408*3/4)),], c(3, 1, 2))) - 1
    
    validation_MLFTS_fit_h1 = array(NA, dim = c(n_PM, n_hour, n_test), dimnames = list(x_grid, 1:n_hour, 1:n_test))
    for(iw in 1:n_test)
    {
        n_train = (408 * 1/2 + iw - 1)
      
        if(fore_model == "MLFTS")
        {
            validation_MLFTS_fit_h1[,,iw] = matrix(unlist(MLFTS_model(data_input = aperm(data_set[,1:n_train,], c(3, 2, 1)),
                                                                      aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
        }
        else if(fore_model == "factor_MLFTS")
        {
            validation_MLFTS_fit_h1[,,iw] = factor_model_fun(data = aperm(data_set[,1:n_train,], c(2, 1, 3)))$fore_val
        }
        print(n_train); rm(n_train); rm(iw)
    }
    
    # difference between holdout data and forecasts (51 x 24 x n_test)
    
    resi_mat = holdout_set - (10^validation_MLFTS_fit_h1 - 1)
    
    # compute quantiles
    
    quantile_resid = matrix(NA, n_PM, n_hour)
    for(iw in 1:n_PM)
    {
        quantile_resid[iw,] = apply(resi_mat[iw,,], 1, function(x) quantile(abs(x), probs = (1 - level_sig)))
        rm(iw)
    }
    rm(holdout_set)
    
    # holdout data in the original scale for the testing period
    
    holdout_set = 10^(aperm(data_set[,((408*3/4+1):408),], c(3, 1, 2))) - 1
    
    # forecasts in the testing period
    
    test_MLFTS_fit_h1 = array(NA, dim = c(n_PM, n_hour, n_test), dimnames = list(x_grid, 1:n_hour, 1:n_test))
    for(iw in 1:n_test)
    {
        n_train = (408 * 3/4 + iw - 1)
      
        if(fore_model == "MLFTS")
        {
            test_MLFTS_fit_h1[,,iw] = matrix(unlist(MLFTS_model(data_input = aperm(data_set[,1:n_train,], c(3, 2, 1)),
                                                            aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)
        }
        else if(fore_model == "factor_MLFTS")
        {
            test_MLFTS_fit_h1[,,iw] = factor_model_fun(data = aperm(data_set[,1:n_train,], c(2, 1, 3)))$fore_val
        }
        print(n_train); rm(n_train); rm(iw)
    }
    
    # compute the lower and upper bounds
    
    test_MLFTS_fit_h1_lb = test_MLFTS_fit_h1_ub = array(NA, dim = dim(test_MLFTS_fit_h1), dimnames = dimnames(test_MLFTS_fit_h1))
    for(iw in 1:n_PM)
    {
        test_MLFTS_fit_h1_lb[iw,,] = (10^(test_MLFTS_fit_h1[iw,,]) - 1) - quantile_resid[iw,]
        test_MLFTS_fit_h1_ub[iw,,] = (10^(test_MLFTS_fit_h1[iw,,]) - 1) + quantile_resid[iw,]
        rm(iw)
    }
    
    # compute empirical coverage probability (ECP), coverage probability difference (CPD), interval score
    
    int_val = matrix(NA, n_PM, 3) 
    for(iw in 1:n_PM)
    {
        int_val[iw,] = interval_score(holdout = holdout_set[iw,,], 
                                      lb = test_MLFTS_fit_h1_lb[iw,,],
                                      ub = test_MLFTS_fit_h1_ub[iw,,], 
                                      alpha = level_sig)
        rm(iw)
    }
    colnames(int_val) = c("ECP", "CPD", "score")
    rm(test_MLFTS_fit_h1); rm(test_MLFTS_fit_h1_lb); rm(test_MLFTS_fit_h1_ub); rm(holdout_set)
    rm(quantile_resid); rm(resi_mat)
    return(int_val)   
}

##############################
# level of significance = 0.2
##############################

PNSD_Mon_conformal = conformal_PI(data_set = PNSD_log_Mon_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Tue_conformal = conformal_PI(data_set = PNSD_log_Tue_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Wed_conformal = conformal_PI(data_set = PNSD_log_Wed_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Thu_conformal = conformal_PI(data_set = PNSD_log_Thu_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Fri_conformal = conformal_PI(data_set = PNSD_log_Fri_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Sat_conformal = conformal_PI(data_set = PNSD_log_Sat_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Sun_conformal = conformal_PI(data_set = PNSD_log_Sun_408, fore_model = "MLFTS", level_sig = 0.2)

# ECP, CPD, score

PNSD_ECP_conformal_test = cbind(PNSD_Mon_conformal[,1], PNSD_Tue_conformal[,1],
                                PNSD_Wed_conformal[,1], PNSD_Thu_conformal[,1], 
                                PNSD_Fri_conformal[,1], PNSD_Sat_conformal[,1], 
                                PNSD_Sun_conformal[,1])

PNSD_CPD_conformal_test = cbind(PNSD_Mon_conformal[,2], PNSD_Tue_conformal[,2],
                                PNSD_Wed_conformal[,2], PNSD_Thu_conformal[,2], 
                                PNSD_Fri_conformal[,2], PNSD_Sat_conformal[,2], 
                                PNSD_Sun_conformal[,2])

PNSD_score_conformal_test = cbind(PNSD_Mon_conformal[,3], PNSD_Tue_conformal[,3],
                                  PNSD_Wed_conformal[,3], PNSD_Thu_conformal[,3], 
                                  PNSD_Fri_conformal[,3], PNSD_Sat_conformal[,3], 
                                  PNSD_Sun_conformal[,3])
rownames(PNSD_CPD_conformal_test) = rownames(PNSD_score_conformal_test) = x_grid
colnames(PNSD_CPD_conformal_test) = colnames(PNSD_score_conformal_test) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

###############################
# level of significance = 0.05
###############################

PNSD_Mon_conformal_0.95 = conformal_PI(data_set = PNSD_log_Mon_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Tue_conformal_0.95 = conformal_PI(data_set = PNSD_log_Tue_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Wed_conformal_0.95 = conformal_PI(data_set = PNSD_log_Wed_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Thu_conformal_0.95 = conformal_PI(data_set = PNSD_log_Thu_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Fri_conformal_0.95 = conformal_PI(data_set = PNSD_log_Fri_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Sat_conformal_0.95 = conformal_PI(data_set = PNSD_log_Sat_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Sun_conformal_0.95 = conformal_PI(data_set = PNSD_log_Sun_408, fore_model = "MLFTS", level_sig = 0.05)

# ECP, CPD, score

PNSD_ECP_conformal_test_0.95 = cbind(PNSD_Mon_conformal_0.95[,1], PNSD_Tue_conformal_0.95[,1],
                                     PNSD_Wed_conformal_0.95[,1], PNSD_Thu_conformal_0.95[,1], 
                                     PNSD_Fri_conformal_0.95[,1], PNSD_Sat_conformal_0.95[,1], 
                                     PNSD_Sun_conformal_0.95[,1])

PNSD_CPD_conformal_test_0.95 = cbind(PNSD_Mon_conformal_0.95[,2], PNSD_Tue_conformal_0.95[,2],
                                     PNSD_Wed_conformal_0.95[,2], PNSD_Thu_conformal_0.95[,2], 
                                     PNSD_Fri_conformal_0.95[,2], PNSD_Sat_conformal_0.95[,2], 
                                     PNSD_Sun_conformal_0.95[,2])

PNSD_score_conformal_test_0.95 = cbind(PNSD_Mon_conformal_0.95[,3], PNSD_Tue_conformal_0.95[,3],
                                       PNSD_Wed_conformal_0.95[,3], PNSD_Thu_conformal_0.95[,3], 
                                       PNSD_Fri_conformal_0.95[,3], PNSD_Sat_conformal_0.95[,3], 
                                       PNSD_Sun_conformal_0.95[,3])
colnames(PNSD_ECP_conformal_test_0.95) = colnames(PNSD_CPD_conformal_test_0.95) = 
colnames(PNSD_score_conformal_test_0.95) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
rownames(PNSD_ECP_conformal_test_0.95) = rownames(PNSD_CPD_conformal_test_0.95) = 
rownames(PNSD_score_conformal_test_0.95) = x_grid

