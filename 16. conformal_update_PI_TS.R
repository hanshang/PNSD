# data_set: an array
# fore_model: MLFTS or factor_MLFTS model
# level_sig: level of significance

sd_PI_test_TS_conformal <- function(data_set, fore_model, level_sig)
{
    n_hour = dim(data_set)[1]
    n_PM = dim(data_set)[3]

    # holdout data in the validation set (original scale)

    holdout_set = array((10^(aperm(data_set[,((408*1/2+1):(408*3/4)),], c(3, 1, 2))) - 1), dim = c(n_PM, n_hour, n_test))
    test_MLFTS_fit_h1 = array(NA, dim = c(n_PM, n_hour, n_test), dimnames = list(x_grid, 1:n_hour, 1:n_test))
    for(iw in 1:n_test)
    {
        n_train = (408 * 1/2 + iw - 1)
    
        if(fore_model == "MLFTS")
        {
            test_MLFTS_fit_h1[,,iw] = matrix(unlist(MLFTS_model(data_input = aperm(data_set[,1:n_train,], c(3, 2, 1)),
                                                                aux_var = NULL, ncomp_method = "provide", fh = 1, 
                                                                fore_method = "ets")),,24)
        }
        else if(fore_model == "factor_MLFTS")
        {
            test_MLFTS_fit_h1[,,iw] = factor_model_fun(data = aperm(data_set[,1:n_train,], c(2, 1, 3)))$fore_val
        } 
        rm(n_train); rm(iw)
    }

    # difference between holdout data and forecasts (51 x 24 x n_test)

    resi_mat = holdout_set - (10^test_MLFTS_fit_h1 - 1)
    rm(test_MLFTS_fit_h1); rm(holdout_set)

    test_MLFTS_fit_h1 = array(NA, dim = c(n_PM, n_hour, n_test), dimnames = list(x_grid, 1:n_hour, 1:n_test))
    for(iw in 1:n_test)
    {
        n_train = (408 * 3/4 + iw - 1)
    
        if(fore_model == "MLFTS")
        {
            test_MLFTS_fit_h1[,,iw] = matrix(unlist(MLFTS_model(data_input = aperm(data_set[,1:n_train,], c(3, 2, 1)),
                                                                aux_var = NULL, ncomp_method = "provide", fh = 1, 
                                                                fore_method = "ets")),,n_hour)
        }
        else if(fore_model == "factor_MLFTS")
        {
            test_MLFTS_fit_h1[,,iw] = factor_model_fun(data = aperm(data_set[,1:n_train,], c(2, 1, 3)))$fore_val
        } 
        rm(n_train); rm(iw)
    }

    # holdout data in the testing period (original scale)
    
    holdout_set = array((10^(aperm(data_set[,((408*3/4+1):408),], c(3, 1, 2))) - 1), dim = c(n_PM, n_hour, n_test))
    
    int_val = array(NA, dim = c(n_PM, 3, (n_hour - 2)), dimnames = list(x_grid, c("ECP", "CPD", "score"), 3:n_hour)) 
    for(iwk in 2:23)
    {
        observed_period = 1:iwk
  
        # compute quantiles

        quantile_resid = matrix(NA, n_PM, (n_hour - length(observed_period)))
        for(iw in 1:n_PM)
        {
            quantile_resid[iw,] = apply(matrix(resi_mat[iw,(1:24)[-observed_period],], length((1:24)[-observed_period]), n_test), 1, 
                                        function(x) quantile(abs(x), probs = (1 - level_sig)))
            rm(iw)
        }
    
        # compute the lower and upper bounds
    
        test_MLFTS_fit_h1_lb = test_MLFTS_fit_h1_ub = array(NA, dim = c(n_PM, (n_hour - length(observed_period)), n_test))
        for(iw in 1:n_PM)
        {
            test_MLFTS_fit_h1_lb[iw,,] = (10^(test_MLFTS_fit_h1[iw,(1:24)[-observed_period],]) - 1) - quantile_resid[iw,]
            test_MLFTS_fit_h1_ub[iw,,] = (10^(test_MLFTS_fit_h1[iw,(1:24)[-observed_period],]) - 1) + quantile_resid[iw,]
            rm(iw)
        }
    
        # compute CPD, interval score
    
        for(iw in 1:n_PM)
        {
            int_val[iw,,(iwk-1)] = interval_score(holdout = matrix(holdout_set[iw,(1:n_hour)[-observed_period],], 
                                                  nrow = length((1:n_hour)[-observed_period]), n_test), 
                                                  lb = matrix(test_MLFTS_fit_h1_lb[iw,,],,n_test),
                                                  ub = matrix(test_MLFTS_fit_h1_ub[iw,,],,n_test), alpha = level_sig)
            rm(iw)
        }
        print(iwk); rm(iwk); rm(observed_period); rm(test_MLFTS_fit_h1_lb); rm(test_MLFTS_fit_h1_ub); rm(quantile_resid)
    }
    return(int_val)
}

########################
## fore_model = "MLFTS"
########################

# level of significance = 0.2

PNSD_Mon_TS_conformal_dynamic = sd_PI_test_TS_conformal(data_set = PNSD_log_Mon_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Tue_TS_conformal_dynamic = sd_PI_test_TS_conformal(data_set = PNSD_log_Tue_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Wed_TS_conformal_dynamic = sd_PI_test_TS_conformal(data_set = PNSD_log_Wed_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Thu_TS_conformal_dynamic = sd_PI_test_TS_conformal(data_set = PNSD_log_Thu_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Fri_TS_conformal_dynamic = sd_PI_test_TS_conformal(data_set = PNSD_log_Fri_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Sat_TS_conformal_dynamic = sd_PI_test_TS_conformal(data_set = PNSD_log_Sat_408, fore_model = "MLFTS", level_sig = 0.2)
PNSD_Sun_TS_conformal_dynamic = sd_PI_test_TS_conformal(data_set = PNSD_log_Sun_408, fore_model = "MLFTS", level_sig = 0.2)

PNSD_days_TS_conformal_dynamic_CPD = rbind(apply(PNSD_Mon_TS_conformal_dynamic[,2,], 2, mean),
                                           apply(PNSD_Tue_TS_conformal_dynamic[,2,], 2, mean),
                                           apply(PNSD_Wed_TS_conformal_dynamic[,2,], 2, mean),
                                           apply(PNSD_Thu_TS_conformal_dynamic[,2,], 2, mean),
                                           apply(PNSD_Fri_TS_conformal_dynamic[,2,], 2, mean),
                                           apply(PNSD_Sat_TS_conformal_dynamic[,2,], 2, mean),
                                           apply(PNSD_Sun_TS_conformal_dynamic[,2,], 2, mean))

PNSD_days_TS_conformal_dynamic_score = rbind(apply(PNSD_Mon_TS_conformal_dynamic[,3,], 2, mean),
                                             apply(PNSD_Tue_TS_conformal_dynamic[,3,], 2, mean),
                                             apply(PNSD_Wed_TS_conformal_dynamic[,3,], 2, mean),
                                             apply(PNSD_Thu_TS_conformal_dynamic[,3,], 2, mean),
                                             apply(PNSD_Fri_TS_conformal_dynamic[,3,], 2, mean),
                                             apply(PNSD_Sat_TS_conformal_dynamic[,3,], 2, mean),
                                             apply(PNSD_Sun_TS_conformal_dynamic[,3,], 2, mean))

rownames(PNSD_days_TS_conformal_dynamic_CPD) = rownames(PNSD_days_TS_conformal_dynamic_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(PNSD_days_TS_conformal_dynamic_CPD) = colnames(PNSD_days_TS_conformal_dynamic_score) = 3:24

# level of significance = 0.05

PNSD_Mon_TS_conformal_dynamic_0.05 = sd_PI_test_TS_conformal(data_set = PNSD_log_Mon_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Tue_TS_conformal_dynamic_0.05 = sd_PI_test_TS_conformal(data_set = PNSD_log_Tue_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Wed_TS_conformal_dynamic_0.05 = sd_PI_test_TS_conformal(data_set = PNSD_log_Wed_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Thu_TS_conformal_dynamic_0.05 = sd_PI_test_TS_conformal(data_set = PNSD_log_Thu_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Fri_TS_conformal_dynamic_0.05 = sd_PI_test_TS_conformal(data_set = PNSD_log_Fri_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Sat_TS_conformal_dynamic_0.05 = sd_PI_test_TS_conformal(data_set = PNSD_log_Sat_408, fore_model = "MLFTS", level_sig = 0.05)
PNSD_Sun_TS_conformal_dynamic_0.05 = sd_PI_test_TS_conformal(data_set = PNSD_log_Sun_408, fore_model = "MLFTS", level_sig = 0.05)

PNSD_days_TS_conformal_dynamic_0.05_CPD = rbind(apply(PNSD_Mon_TS_conformal_dynamic_0.05[,2,], 2, mean),
                                                apply(PNSD_Tue_TS_conformal_dynamic_0.05[,2,], 2, mean),
                                                apply(PNSD_Wed_TS_conformal_dynamic_0.05[,2,], 2, mean),
                                                apply(PNSD_Thu_TS_conformal_dynamic_0.05[,2,], 2, mean),
                                                apply(PNSD_Fri_TS_conformal_dynamic_0.05[,2,], 2, mean),
                                                apply(PNSD_Sat_TS_conformal_dynamic_0.05[,2,], 2, mean),
                                                apply(PNSD_Sun_TS_conformal_dynamic_0.05[,2,], 2, mean))

PNSD_days_TS_conformal_dynamic_0.05_score = rbind(apply(PNSD_Mon_TS_conformal_dynamic_0.05[,3,], 2, mean),
                                                  apply(PNSD_Tue_TS_conformal_dynamic_0.05[,3,], 2, mean),
                                                  apply(PNSD_Wed_TS_conformal_dynamic_0.05[,3,], 2, mean),
                                                  apply(PNSD_Thu_TS_conformal_dynamic_0.05[,3,], 2, mean),
                                                  apply(PNSD_Fri_TS_conformal_dynamic_0.05[,3,], 2, mean),
                                                  apply(PNSD_Sat_TS_conformal_dynamic_0.05[,3,], 2, mean),
                                                  apply(PNSD_Sun_TS_conformal_dynamic_0.05[,3,], 2, mean))

rownames(PNSD_days_TS_conformal_dynamic_0.05_CPD) = rownames(PNSD_days_TS_conformal_dynamic_0.05_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(PNSD_days_TS_conformal_dynamic_0.05_CPD) = colnames(PNSD_days_TS_conformal_dynamic_0.05_score) = 3:24
  
###############################
## fore_model = "factor_MLFTS"
###############################

## level of significance = 0.2

PNSD_Mon_TS_conformal_dynamic_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Mon_408, 
                                                      fore_model = "factor_MLFTS", level_sig = 0.2)
PNSD_Tue_TS_conformal_dynamic_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Tue_408, 
                                                      fore_model = "factor_MLFTS", level_sig = 0.2)
PNSD_Wed_TS_conformal_dynamic_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Wed_408, 
                                                      fore_model = "factor_MLFTS", level_sig = 0.2)
PNSD_Thu_TS_conformal_dynamic_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Thu_408, 
                                                      fore_model = "factor_MLFTS", level_sig = 0.2)
PNSD_Fri_TS_conformal_dynamic_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Fri_408, 
                                                      fore_model = "factor_MLFTS", level_sig = 0.2)
PNSD_Sat_TS_conformal_dynamic_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Sat_408, 
                                                      fore_model = "factor_MLFTS", level_sig = 0.2)
PNSD_Sun_TS_conformal_dynamic_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Sun_408, 
                                                      fore_model = "factor_MLFTS", level_sig = 0.2)

# CPD

factor_PNSD_days_TS_conformal_dynamic_CPD = rbind(apply(PNSD_Mon_TS_conformal_dynamic_factor_MLFTS[,2,], 2, mean),
                                                  apply(PNSD_Tue_TS_conformal_dynamic_factor_MLFTS[,2,], 2, mean),
                                                  apply(PNSD_Wed_TS_conformal_dynamic_factor_MLFTS[,2,], 2, mean),
                                                  apply(PNSD_Thu_TS_conformal_dynamic_factor_MLFTS[,2,], 2, mean),
                                                  apply(PNSD_Fri_TS_conformal_dynamic_factor_MLFTS[,2,], 2, mean),
                                                  apply(PNSD_Sat_TS_conformal_dynamic_factor_MLFTS[,2,], 2, mean),
                                                  apply(PNSD_Sun_TS_conformal_dynamic_factor_MLFTS[,2,], 2, mean))

# interval score

factor_PNSD_days_TS_conformal_dynamic_score = rbind(apply(PNSD_Mon_TS_conformal_dynamic_factor_MLFTS[,3,], 2, mean),
                                                    apply(PNSD_Tue_TS_conformal_dynamic_factor_MLFTS[,3,], 2, mean),
                                                    apply(PNSD_Wed_TS_conformal_dynamic_factor_MLFTS[,3,], 2, mean),
                                                    apply(PNSD_Thu_TS_conformal_dynamic_factor_MLFTS[,3,], 2, mean),
                                                    apply(PNSD_Fri_TS_conformal_dynamic_factor_MLFTS[,3,], 2, mean),
                                                    apply(PNSD_Sat_TS_conformal_dynamic_factor_MLFTS[,3,], 2, mean),
                                                    apply(PNSD_Sun_TS_conformal_dynamic_factor_MLFTS[,3,], 2, mean))

rownames(factor_PNSD_days_TS_conformal_dynamic_CPD) = rownames(factor_PNSD_days_TS_conformal_dynamic_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(factor_PNSD_days_TS_conformal_dynamic_CPD) = colnames(factor_PNSD_days_TS_conformal_dynamic_score) = 3:24  

## level of significance = 0.05

PNSD_Mon_TS_conformal_dynamic_0.05_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Mon_408, 
                                                                     fore_model = "factor_MLFTS", level_sig = 0.05)
PNSD_Tue_TS_conformal_dynamic_0.05_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Tue_408, 
                                                                     fore_model = "factor_MLFTS", level_sig = 0.05)
PNSD_Wed_TS_conformal_dynamic_0.05_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Wed_408, 
                                                                     fore_model = "factor_MLFTS", level_sig = 0.05)
PNSD_Thu_TS_conformal_dynamic_0.05_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Thu_408, 
                                                                     fore_model = "factor_MLFTS", level_sig = 0.05)
PNSD_Fri_TS_conformal_dynamic_0.05_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Fri_408, 
                                                                     fore_model = "factor_MLFTS", level_sig = 0.05)
PNSD_Sat_TS_conformal_dynamic_0.05_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Sat_408, 
                                                                     fore_model = "factor_MLFTS", level_sig = 0.05)
PNSD_Sun_TS_conformal_dynamic_0.05_factor_MLFTS = sd_PI_test_TS_conformal(data_set = PNSD_log_Sun_408, 
                                                                     fore_model = "factor_MLFTS", level_sig = 0.05)

# CPD

factor_PNSD_days_TS_conformal_dynamic_0.05_CPD = rbind(apply(PNSD_Mon_TS_conformal_dynamic_0.05_factor_MLFTS[,2,], 2, mean),
                                                             apply(PNSD_Tue_TS_conformal_dynamic_0.05_factor_MLFTS[,2,], 2, mean),
                                                             apply(PNSD_Wed_TS_conformal_dynamic_0.05_factor_MLFTS[,2,], 2, mean),
                                                             apply(PNSD_Thu_TS_conformal_dynamic_0.05_factor_MLFTS[,2,], 2, mean),
                                                             apply(PNSD_Fri_TS_conformal_dynamic_0.05_factor_MLFTS[,2,], 2, mean),
                                                             apply(PNSD_Sat_TS_conformal_dynamic_0.05_factor_MLFTS[,2,], 2, mean),
                                                             apply(PNSD_Sun_TS_conformal_dynamic_0.05_factor_MLFTS[,2,], 2, mean))

# interval score

factor_PNSD_days_TS_conformal_dynamic_0.05_score = rbind(apply(PNSD_Mon_TS_conformal_dynamic_0.05_factor_MLFTS[,3,], 2, mean),
                                                               apply(PNSD_Tue_TS_conformal_dynamic_0.05_factor_MLFTS[,3,], 2, mean),
                                                               apply(PNSD_Wed_TS_conformal_dynamic_0.05_factor_MLFTS[,3,], 2, mean),
                                                               apply(PNSD_Thu_TS_conformal_dynamic_0.05_factor_MLFTS[,3,], 2, mean),
                                                               apply(PNSD_Fri_TS_conformal_dynamic_0.05_factor_MLFTS[,3,], 2, mean),
                                                               apply(PNSD_Sat_TS_conformal_dynamic_0.05_factor_MLFTS[,3,], 2, mean),
                                                               apply(PNSD_Sun_TS_conformal_dynamic_0.05_factor_MLFTS[,3,], 2, mean))

rownames(factor_PNSD_days_TS_conformal_dynamic_0.05_CPD) = rownames(factor_PNSD_days_TS_conformal_dynamic_0.05_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(factor_PNSD_days_TS_conformal_dynamic_0.05_CPD) = colnames(factor_PNSD_days_TS_conformal_dynamic_0.05_score) = 3:24

