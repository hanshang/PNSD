##################################################
# TS method for constructing prediction intervals
##################################################

# data_set: an array
# fore_model: MLFTS or factor_MLFTS
# tune_para_select: selected tuning parameter
# sd_val_mat: standard deviation
# level_sig: level of significance

sd_PI_test_TS <- function(data_set, fore_model, tune_para_select, sd_val_mat, level_sig)
{
    n_hour = dim(data_set)[1]
    n_PM = dim(data_set)[3]
    
    holdout_set = 10^(aperm(data_set[,((408*3/4+1):408),], c(3, 1, 2))) - 1
    
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
    
    int_val = array(NA, dim = c(n_PM, 3, 22), dimnames = list(x_grid, c("ECP", "CPD", "score"), 3:24)) 
    for(iwk in 2:23)
    {
        observed_period = 1:iwk
    
        test_MLFTS_fit_h1_lb = test_MLFTS_fit_h1_ub = array(NA, dim = c(n_PM, n_hour - length(observed_period), n_test), 
                                                            dimnames = list(x_grid, (1:n_hour)[-observed_period], 1:n_test))
        for(iw in 1:n_PM)
        {
            test_MLFTS_fit_h1_lb[iw,,] = (10^(test_MLFTS_fit_h1[iw,(1:24)[-observed_period],]) - 1) - tune_para_select[iw] * sd_val_mat[iw,(1:24)[-observed_period]]
            test_MLFTS_fit_h1_ub[iw,,] = (10^(test_MLFTS_fit_h1[iw,(1:24)[-observed_period],]) - 1) + tune_para_select[iw] * sd_val_mat[iw,(1:24)[-observed_period]]
            rm(iw)
        }
    
        # compute CPD, interval score
    
        for(iw in 1:n_PM)
        {
            int_val[iw,,(iwk-1)] = interval_score(holdout = matrix(holdout_set[iw,(1:n_hour)[-observed_period],], nrow = length((1:n_hour)[-observed_period]), n_test), 
                                          lb = matrix(test_MLFTS_fit_h1_lb[iw,,],,n_test),
                                          ub = matrix(test_MLFTS_fit_h1_ub[iw,,],,n_test), alpha = level_sig)
            rm(iw)
        }
        colnames(int_val) = c("ECP", "CPD", "score")
        print(iwk); rm(iwk); rm(observed_period)
    }
    return(int_val)   
}

##############################
# level of significance = 0.2
##############################

# MLFTS

PNSD_Mon_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Mon_408, fore_model = "MLFTS",
                                    tune_para_select = PNSD_Mon_sd_tune$tune_para_find,
                                    sd_val_mat = PNSD_Mon_sd_tune$sd_val, level_sig = 0.2)

PNSD_Tue_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Tue_408, fore_model = "MLFTS",
                                    tune_para_select = PNSD_Tue_sd_tune$tune_para_find,
                                    sd_val_mat = PNSD_Tue_sd_tune$sd_val, level_sig = 0.2)

PNSD_Wed_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Wed_408, fore_model = "MLFTS",
                                    tune_para_select = PNSD_Wed_sd_tune$tune_para_find,
                                    sd_val_mat = PNSD_Wed_sd_tune$sd_val, level_sig = 0.2)

PNSD_Thu_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Thu_408, fore_model = "MLFTS",
                                    tune_para_select = PNSD_Thu_sd_tune$tune_para_find,
                                    sd_val_mat = PNSD_Thu_sd_tune$sd_val, level_sig = 0.2)

PNSD_Fri_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Fri_408, fore_model = "MLFTS",
                                    tune_para_select = PNSD_Fri_sd_tune$tune_para_find,
                                    sd_val_mat = PNSD_Fri_sd_tune$sd_val, level_sig = 0.2)

PNSD_Sat_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Sat_408, fore_model = "MLFTS",
                                    tune_para_select = PNSD_Sat_sd_tune$tune_para_find,
                                    sd_val_mat = PNSD_Sat_sd_tune$sd_val, level_sig = 0.2)

PNSD_Sun_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Sun_408, fore_model = "MLFTS",
                                    tune_para_select = PNSD_Sun_sd_tune$tune_para_find,
                                    sd_val_mat = PNSD_Sun_sd_tune$sd_val, level_sig = 0.2)

# CPD

PNSD_days_sd_test_TS_CPD = rbind(apply(PNSD_Mon_sd_test_TS[,2,], 2, mean),
                                 apply(PNSD_Tue_sd_test_TS[,2,], 2, mean),
                                 apply(PNSD_Wed_sd_test_TS[,2,], 2, mean),
                                 apply(PNSD_Thu_sd_test_TS[,2,], 2, mean),
                                 apply(PNSD_Fri_sd_test_TS[,2,], 2, mean),
                                 apply(PNSD_Sat_sd_test_TS[,2,], 2, mean),
                                 apply(PNSD_Sun_sd_test_TS[,2,], 2, mean))

# interval score

PNSD_days_sd_test_TS_score = rbind(apply(PNSD_Mon_sd_test_TS[,3,], 2, mean),
                                   apply(PNSD_Tue_sd_test_TS[,3,], 2, mean),
                                   apply(PNSD_Wed_sd_test_TS[,3,], 2, mean),
                                   apply(PNSD_Thu_sd_test_TS[,3,], 2, mean),
                                   apply(PNSD_Fri_sd_test_TS[,3,], 2, mean),
                                   apply(PNSD_Sat_sd_test_TS[,3,], 2, mean),
                                   apply(PNSD_Sun_sd_test_TS[,3,], 2, mean))

rownames(PNSD_days_sd_test_TS_CPD) = rownames(PNSD_days_sd_test_TS_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(PNSD_days_sd_test_TS_CPD) = colnames(PNSD_days_sd_test_TS_score) = 3:24
  
# factor_MLFTS

factor_PNSD_Mon_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Mon_408, fore_model = "factor_MLFTS",
                                           tune_para_select = factor_PNSD_Mon_sd_tune$tune_para_find,
                                           sd_val_mat = factor_PNSD_Mon_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Tue_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Tue_408, fore_model = "factor_MLFTS",
                                           tune_para_select = factor_PNSD_Tue_sd_tune$tune_para_find,
                                           sd_val_mat = factor_PNSD_Tue_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Wed_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Wed_408, fore_model = "factor_MLFTS",
                                           tune_para_select = factor_PNSD_Wed_sd_tune$tune_para_find,
                                           sd_val_mat = factor_PNSD_Wed_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Thu_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Thu_408, fore_model = "factor_MLFTS",
                                           tune_para_select = factor_PNSD_Thu_sd_tune$tune_para_find,
                                           sd_val_mat = factor_PNSD_Thu_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Fri_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Fri_408, fore_model = "factor_MLFTS",
                                           tune_para_select = factor_PNSD_Fri_sd_tune$tune_para_find,
                                           sd_val_mat = factor_PNSD_Fri_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Sat_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Sat_408, fore_model = "factor_MLFTS",
                                           tune_para_select = factor_PNSD_Sat_sd_tune$tune_para_find,
                                           sd_val_mat = factor_PNSD_Sat_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Sun_sd_test_TS = sd_PI_test_TS(data_set = PNSD_log_Sun_408, fore_model = "factor_MLFTS",
                                           tune_para_select = factor_PNSD_Sun_sd_tune$tune_para_find,
                                           sd_val_mat = factor_PNSD_Sun_sd_tune$sd_val, level_sig = 0.2)

# CPD

factor_PNSD_days_sd_test_TS_CPD = rbind(apply(factor_PNSD_Mon_sd_test_TS[,2,], 2, mean),
                                        apply(factor_PNSD_Tue_sd_test_TS[,2,], 2, mean),
                                        apply(factor_PNSD_Wed_sd_test_TS[,2,], 2, mean),
                                        apply(factor_PNSD_Thu_sd_test_TS[,2,], 2, mean),
                                        apply(factor_PNSD_Fri_sd_test_TS[,2,], 2, mean),
                                        apply(factor_PNSD_Sat_sd_test_TS[,2,], 2, mean),
                                        apply(factor_PNSD_Sun_sd_test_TS[,2,], 2, mean))

# interval score

factor_PNSD_days_sd_test_TS_score = rbind(apply(factor_PNSD_Mon_sd_test_TS[,3,], 2, mean),
                                          apply(factor_PNSD_Tue_sd_test_TS[,3,], 2, mean),
                                          apply(factor_PNSD_Wed_sd_test_TS[,3,], 2, mean),
                                          apply(factor_PNSD_Thu_sd_test_TS[,3,], 2, mean),
                                          apply(factor_PNSD_Fri_sd_test_TS[,3,], 2, mean),
                                          apply(factor_PNSD_Sat_sd_test_TS[,3,], 2, mean),
                                          apply(factor_PNSD_Sun_sd_test_TS[,3,], 2, mean))

rownames(factor_PNSD_days_sd_test_TS_CPD) = rownames(factor_PNSD_days_sd_test_TS_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(factor_PNSD_days_sd_test_TS_CPD) = colnames(factor_PNSD_days_sd_test_TS_score) = 3:24

###############################
# level of significance = 0.05
###############################

# MLFTS

PNSD_Mon_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Mon_408, fore_model = "MLFTS",
                                         tune_para_select = PNSD_Mon_sd_tune_0.95$tune_para_find,
                                         sd_val_mat = PNSD_Mon_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Tue_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Tue_408, fore_model = "MLFTS",
                                         tune_para_select = PNSD_Tue_sd_tune_0.95$tune_para_find,
                                         sd_val_mat = PNSD_Tue_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Wed_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Wed_408, fore_model = "MLFTS",
                                         tune_para_select = PNSD_Wed_sd_tune_0.95$tune_para_find,
                                         sd_val_mat = PNSD_Wed_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Thu_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Thu_408, fore_model = "MLFTS",
                                         tune_para_select = PNSD_Thu_sd_tune_0.95$tune_para_find,
                                         sd_val_mat = PNSD_Thu_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Fri_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Fri_408, fore_model = "MLFTS",
                                         tune_para_select = PNSD_Fri_sd_tune_0.95$tune_para_find,
                                         sd_val_mat = PNSD_Fri_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Sat_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Sat_408, fore_model = "MLFTS",
                                         tune_para_select = PNSD_Sat_sd_tune_0.95$tune_para_find,
                                         sd_val_mat = PNSD_Sat_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Sun_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Sun_408, fore_model = "MLFTS",
                                         tune_para_select = PNSD_Sun_sd_tune$tune_para_find,
                                         sd_val_mat = PNSD_Sun_sd_tune$sd_val, level_sig = 0.05)

# CPD

PNSD_days_sd_test_TS_CPD_0.05 = rbind(apply(PNSD_Mon_sd_test_TS_0.05[,2,], 2, mean),
                                      apply(PNSD_Tue_sd_test_TS_0.05[,2,], 2, mean),
                                      apply(PNSD_Wed_sd_test_TS_0.05[,2,], 2, mean),
                                      apply(PNSD_Thu_sd_test_TS_0.05[,2,], 2, mean),
                                      apply(PNSD_Fri_sd_test_TS_0.05[,2,], 2, mean),
                                      apply(PNSD_Sat_sd_test_TS_0.05[,2,], 2, mean),
                                      apply(PNSD_Sun_sd_test_TS_0.05[,2,], 2, mean))

# interval score

PNSD_days_sd_test_TS_score_0.05 = rbind(apply(PNSD_Mon_sd_test_TS_0.05[,3,], 2, mean),
                                        apply(PNSD_Tue_sd_test_TS_0.05[,3,], 2, mean),
                                        apply(PNSD_Wed_sd_test_TS_0.05[,3,], 2, mean),
                                        apply(PNSD_Thu_sd_test_TS_0.05[,3,], 2, mean),
                                        apply(PNSD_Fri_sd_test_TS_0.05[,3,], 2, mean),
                                        apply(PNSD_Sat_sd_test_TS_0.05[,3,], 2, mean),
                                        apply(PNSD_Sun_sd_test_TS_0.05[,3,], 2, mean))

rownames(PNSD_days_sd_test_TS_CPD_0.05) = rownames(PNSD_days_sd_test_TS_score_0.05) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(PNSD_days_sd_test_TS_CPD_0.05) = colnames(PNSD_days_sd_test_TS_score_0.05) = 3:24

# factor_MLFTS

factor_PNSD_Mon_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Mon_408, fore_model = "factor_MLFTS",
                                                tune_para_select = factor_PNSD_Mon_sd_tune_0.95$tune_para_find,
                                                sd_val_mat = factor_PNSD_Mon_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Tue_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Tue_408, fore_model = "factor_MLFTS",
                                                tune_para_select = factor_PNSD_Tue_sd_tune_0.95$tune_para_find,
                                                sd_val_mat = factor_PNSD_Tue_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Wed_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Wed_408, fore_model = "factor_MLFTS",
                                                tune_para_select = factor_PNSD_Wed_sd_tune_0.95$tune_para_find,
                                                sd_val_mat = factor_PNSD_Wed_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Thu_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Thu_408, fore_model = "factor_MLFTS",
                                                tune_para_select = factor_PNSD_Thu_sd_tune_0.95$tune_para_find,
                                                sd_val_mat = factor_PNSD_Thu_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Fri_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Fri_408, fore_model = "factor_MLFTS",
                                                tune_para_select = factor_PNSD_Fri_sd_tune_0.95$tune_para_find,
                                                sd_val_mat = factor_PNSD_Fri_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Sat_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Sat_408, fore_model = "factor_MLFTS",
                                                tune_para_select = factor_PNSD_Sat_sd_tune_0.95$tune_para_find,
                                                sd_val_mat = factor_PNSD_Sat_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Sun_sd_test_TS_0.05 = sd_PI_test_TS(data_set = PNSD_log_Sun_408, fore_model = "factor_MLFTS",
                                                tune_para_select = factor_PNSD_Sun_sd_tune_0.95$tune_para_find,
                                                sd_val_mat = factor_PNSD_Sun_sd_tune_0.95$sd_val, level_sig = 0.05)

# CPD

factor_PNSD_days_sd_test_TS_0.05_CPD = rbind(apply(factor_PNSD_Mon_sd_test_TS_0.05[,2,], 2, mean),
                                             apply(factor_PNSD_Tue_sd_test_TS_0.05[,2,], 2, mean),
                                             apply(factor_PNSD_Wed_sd_test_TS_0.05[,2,], 2, mean),
                                             apply(factor_PNSD_Thu_sd_test_TS_0.05[,2,], 2, mean),
                                             apply(factor_PNSD_Fri_sd_test_TS_0.05[,2,], 2, mean),
                                             apply(factor_PNSD_Sat_sd_test_TS_0.05[,2,], 2, mean),
                                             apply(factor_PNSD_Sun_sd_test_TS_0.05[,2,], 2, mean))

# interval score

factor_PNSD_days_sd_test_TS_0.05_score = rbind(apply(factor_PNSD_Mon_sd_test_TS_0.05[,3,], 2, mean),
                                               apply(factor_PNSD_Tue_sd_test_TS_0.05[,3,], 2, mean),
                                               apply(factor_PNSD_Wed_sd_test_TS_0.05[,3,], 2, mean),
                                               apply(factor_PNSD_Thu_sd_test_TS_0.05[,3,], 2, mean),
                                               apply(factor_PNSD_Fri_sd_test_TS_0.05[,3,], 2, mean),
                                               apply(factor_PNSD_Sat_sd_test_TS_0.05[,3,], 2, mean),
                                               apply(factor_PNSD_Sun_sd_test_TS_0.05[,3,], 2, mean))

rownames(factor_PNSD_days_sd_test_TS_0.05_CPD) = rownames(factor_PNSD_days_sd_test_TS_0.05_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(factor_PNSD_days_sd_test_TS_0.05_CPD) = colnames(factor_PNSD_days_sd_test_TS_0.05_score) = 3:24

