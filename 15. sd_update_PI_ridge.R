##############
# sd approach
# ridge
##############

# lambda_val: lambda value
# data_set: 24 x T x 51
# observed_period: observation period
# fore_method: forecasting method
# ncp: nominal coverage probability

sd_PI_ridge <- function(lambda_val, data_set, observed_period, fore_method, ncp)
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
        rm(iw)
    }
    
    # holdout data in the original scale
    
    holdout_set = array((10^(aperm(data_set[,((408*1/2+1):(408*3/4)),], c(3, 1, 2))) - 1)[,(1:24)[-(1:observed_period)],],
                        dim = c(n_PM, length((1:24)[-(1:observed_period)]), n_test))
    
    # difference between holdout data and forecasts (51 x 24 x n_test)
    
    resi_mat = holdout_set - ridge_h1_fore
    
    # compute pointwise standard deviation
    
    sd_val = matrix(NA, n_PM, (n_hour - observed_period))
    for(iw in 1:n_PM)
    {
        sd_val[iw,] = apply(matrix(resi_mat[iw,,], length((1:24)[-(1:observed_period)]), n_test), 1, sd)
        rm(iw)
    }
    
    tune_para_find = obj_val_min = vector("numeric", n_PM)
    for(iw in 1:n_PM)
    {
        tune_para_find_val_1 = optimise(f = tune_para_find_function, interval = c(0, 1),
                                        resi_mat = matrix(resi_mat[iw,,], length((1:24)[-(1:observed_period)]), n_test), 
                                        sd_val_input = sd_val[iw,],
                                        alpha_level = ncp)
        
        tune_para_find_val_2 = optimise(f = tune_para_find_function, interval = c(0, 5),
                                        resi_mat = matrix(resi_mat[iw,,], length((1:24)[-(1:observed_period)]), n_test), 
                                        sd_val_input = sd_val[iw,],
                                        alpha_level = ncp)
        
        tune_para_find_val_3 = optimise(f = tune_para_find_function, interval = c(0, 10),
                                        resi_mat = matrix(resi_mat[iw,,], length((1:24)[-(1:observed_period)]), n_test), 
                                        sd_val_input = sd_val[iw,],
                                        alpha_level = ncp)
        
        tune_para_find_val_4 = optimise(f = tune_para_find_function, interval = c(0, 20),
                                        resi_mat = matrix(resi_mat[iw,,], length((1:24)[-(1:observed_period)]), n_test), 
                                        sd_val_input = sd_val[iw,],
                                        alpha_level = ncp)
        
        tune_para_find_val_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B",
                                     resi_mat = matrix(resi_mat[iw,,], length((1:24)[-(1:observed_period)]), n_test), 
                                     sd_val_input = sd_val[iw,],
                                     alpha_level = ncp)
        
        tune_para_find_val_6 = optim(par = 1, fn = tune_para_find_function,
                                     resi_mat = matrix(resi_mat[iw,,], length((1:24)[-(1:observed_period)]), n_test), 
                                     sd_val_input = sd_val[iw,],
                                     alpha_level = ncp)
        
        obj_val = c(tune_para_find_val_1$objective, 
                    tune_para_find_val_2$objective, 
                    tune_para_find_val_3$objective, 
                    tune_para_find_val_4$objective, 
                    tune_para_find_val_5$value, 
                    tune_para_find_val_6$value)
        obj_val_min[iw] = min(obj_val)
        
        tune_para_find[iw] = c(tune_para_find_val_1$minimum, 
                               tune_para_find_val_2$minimum,
                               tune_para_find_val_3$minimum,
                               tune_para_find_val_4$minimum,
                               tune_para_find_val_5$par,
                               tune_para_find_val_6$par)[which.min(obj_val)]
        rm(obj_val); rm(iw); rm(tune_para_find_val_1); rm(tune_para_find_val_2)
        rm(tune_para_find_val_3); rm(tune_para_find_val_4); rm(tune_para_find_val_5)
        rm(tune_para_find_val_6)
    }
    return(list(tune_para_find = tune_para_find, obj_val_min = obj_val_min, sd_val = sd_val))
}

##############################
# level of significance = 0.2
##############################

registerDoMC(detectCores())
PNSD_Mon_sd_ridge_tune = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Mon_para[iwk], 
                                                                 data_set = PNSD_log_Mon_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.8)

PNSD_Tue_sd_ridge_tune = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Tue_para[iwk], 
                                                                 data_set = PNSD_log_Tue_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.8)

PNSD_Wed_sd_ridge_tune = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Wed_para[iwk], 
                                                                 data_set = PNSD_log_Wed_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.8)

PNSD_Thu_sd_ridge_tune = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Thu_para[iwk], 
                                                                 data_set = PNSD_log_Thu_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.8)

PNSD_Fri_sd_ridge_tune = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Fri_para[iwk], 
                                                                 data_set = PNSD_log_Fri_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.8)

PNSD_Sat_sd_ridge_tune = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Sat_para[iwk], 
                                                                 data_set = PNSD_log_Sat_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.8)

PNSD_Sun_sd_ridge_tune = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Sun_para[iwk], 
                                                                 data_set = PNSD_log_Sun_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.8)

# tuning parameter

PNSD_Mon_sd_ridge_tune_para = PNSD_Tue_sd_ridge_tune_para = PNSD_Wed_sd_ridge_tune_para = 
PNSD_Thu_sd_ridge_tune_para = PNSD_Fri_sd_ridge_tune_para = PNSD_Sat_sd_ridge_tune_para = 
PNSD_Sun_sd_ridge_tune_para = matrix(NA, 51, 22)
for(iw in 1:22)
{
    PNSD_Mon_sd_ridge_tune_para[,iw] = PNSD_Mon_sd_ridge_tune[[iw]]$tune_para_find
    PNSD_Tue_sd_ridge_tune_para[,iw] = PNSD_Tue_sd_ridge_tune[[iw]]$tune_para_find
    PNSD_Wed_sd_ridge_tune_para[,iw] = PNSD_Wed_sd_ridge_tune[[iw]]$tune_para_find
    PNSD_Thu_sd_ridge_tune_para[,iw] = PNSD_Thu_sd_ridge_tune[[iw]]$tune_para_find
    PNSD_Fri_sd_ridge_tune_para[,iw] = PNSD_Fri_sd_ridge_tune[[iw]]$tune_para_find
    PNSD_Sat_sd_ridge_tune_para[,iw] = PNSD_Sat_sd_ridge_tune[[iw]]$tune_para_find
    PNSD_Sun_sd_ridge_tune_para[,iw] = PNSD_Sun_sd_ridge_tune[[iw]]$tune_para_find
    rm(iw)
}

colnames(PNSD_Mon_sd_ridge_tune_para) = colnames(PNSD_Tue_sd_ridge_tune_para) =
colnames(PNSD_Wed_sd_ridge_tune_para) = colnames(PNSD_Thu_sd_ridge_tune_para) =
colnames(PNSD_Fri_sd_ridge_tune_para) = colnames(PNSD_Sat_sd_ridge_tune_para) =
colnames(PNSD_Sun_sd_ridge_tune_para) = 3:24

rownames(PNSD_Mon_sd_ridge_tune_para) = rownames(PNSD_Tue_sd_ridge_tune_para) = 
rownames(PNSD_Wed_sd_ridge_tune_para) = rownames(PNSD_Thu_sd_ridge_tune_para) = 
rownames(PNSD_Fri_sd_ridge_tune_para) = rownames(PNSD_Sat_sd_ridge_tune_para) = 
rownames(PNSD_Sun_sd_ridge_tune_para) = x_grid

# objective

PNSD_Mon_sd_ridge_obj = PNSD_Tue_sd_ridge_obj = PNSD_Wed_sd_ridge_obj = 
PNSD_Thu_sd_ridge_obj = PNSD_Fri_sd_ridge_obj = PNSD_Sat_sd_ridge_obj = PNSD_Sun_sd_ridge_obj = matrix(NA, 51, 22)
for(iw in 1:22)
{
    PNSD_Mon_sd_ridge_obj[,iw] = PNSD_Mon_sd_ridge_tune[[1]]$obj_val_min
    PNSD_Tue_sd_ridge_obj[,iw] = PNSD_Tue_sd_ridge_tune[[1]]$obj_val_min
    PNSD_Wed_sd_ridge_obj[,iw] = PNSD_Wed_sd_ridge_tune[[1]]$obj_val_min
    PNSD_Thu_sd_ridge_obj[,iw] = PNSD_Thu_sd_ridge_tune[[1]]$obj_val_min
    PNSD_Fri_sd_ridge_obj[,iw] = PNSD_Fri_sd_ridge_tune[[1]]$obj_val_min
    PNSD_Sat_sd_ridge_obj[,iw] = PNSD_Sat_sd_ridge_tune[[1]]$obj_val_min
    PNSD_Sun_sd_ridge_obj[,iw] = PNSD_Sun_sd_ridge_tune[[1]]$obj_val_min
    print(iw); rm(iw)
}

colnames(PNSD_Mon_sd_ridge_obj) = colnames(PNSD_Tue_sd_ridge_obj) = colnames(PNSD_Wed_sd_ridge_obj) = 
colnames(PNSD_Thu_sd_ridge_obj) = colnames(PNSD_Fri_sd_ridge_obj) = colnames(PNSD_Sat_sd_ridge_obj) = 
colnames(PNSD_Sun_sd_ridge_obj) = 3:24

rownames(PNSD_Mon_sd_ridge_obj) = rownames(PNSD_Tue_sd_ridge_obj) = rownames(PNSD_Wed_sd_ridge_obj) = 
rownames(PNSD_Thu_sd_ridge_obj) = rownames(PNSD_Fri_sd_ridge_obj) = rownames(PNSD_Sat_sd_ridge_obj) = 
rownames(PNSD_Sun_sd_ridge_obj) = x_grid
  

###############################
# level of significance = 0.05
###############################

registerDoMC(detectCores())
PNSD_Mon_sd_ridge_tune_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Mon_para[iwk], 
                                                                 data_set = PNSD_log_Mon_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.95)

PNSD_Tue_sd_ridge_tune_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Tue_para[iwk], 
                                                                 data_set = PNSD_log_Tue_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.95)

PNSD_Wed_sd_ridge_tune_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Wed_para[iwk], 
                                                                 data_set = PNSD_log_Wed_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.95)

PNSD_Thu_sd_ridge_tune_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Thu_para[iwk], 
                                                                 data_set = PNSD_log_Thu_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.95)

PNSD_Fri_sd_ridge_tune_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Fri_para[iwk], 
                                                                 data_set = PNSD_log_Fri_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.95)

PNSD_Sat_sd_ridge_tune_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Sat_para[iwk], 
                                                                 data_set = PNSD_log_Sat_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.95)

PNSD_Sun_sd_ridge_tune_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge(lambda_val = ridge_lambda_log_Sun_para[iwk], 
                                                                 data_set = PNSD_log_Sun_408, 
                                                                 observed_period = (iwk + 1), 
                                                                 fore_method = "ets", ncp = 0.95)
    
# tuning parameter

PNSD_Mon_sd_ridge_tune_para_0.05 = PNSD_Tue_sd_ridge_tune_para_0.05 = PNSD_Wed_sd_ridge_tune_para_0.05 = 
PNSD_Thu_sd_ridge_tune_para_0.05 = PNSD_Fri_sd_ridge_tune_para_0.05 = PNSD_Sat_sd_ridge_tune_para_0.05 = 
PNSD_Sun_sd_ridge_tune_para_0.05 = matrix(NA, 51, 22)
for(iw in 1:22)
{
    PNSD_Mon_sd_ridge_tune_para_0.05[,iw] = PNSD_Mon_sd_ridge_tune_0.05[[iw]]$tune_para_find
    PNSD_Tue_sd_ridge_tune_para_0.05[,iw] = PNSD_Tue_sd_ridge_tune_0.05[[iw]]$tune_para_find
    PNSD_Wed_sd_ridge_tune_para_0.05[,iw] = PNSD_Wed_sd_ridge_tune_0.05[[iw]]$tune_para_find
    PNSD_Thu_sd_ridge_tune_para_0.05[,iw] = PNSD_Thu_sd_ridge_tune_0.05[[iw]]$tune_para_find
    PNSD_Fri_sd_ridge_tune_para_0.05[,iw] = PNSD_Fri_sd_ridge_tune_0.05[[iw]]$tune_para_find
    PNSD_Sat_sd_ridge_tune_para_0.05[,iw] = PNSD_Sat_sd_ridge_tune_0.05[[iw]]$tune_para_find
    PNSD_Sun_sd_ridge_tune_para_0.05[,iw] = PNSD_Sun_sd_ridge_tune_0.05[[iw]]$tune_para_find
    print(iw); rm(iw)
}

colnames(PNSD_Mon_sd_ridge_tune_para_0.05) = colnames(PNSD_Tue_sd_ridge_tune_para_0.05) =
colnames(PNSD_Wed_sd_ridge_tune_para_0.05) = colnames(PNSD_Thu_sd_ridge_tune_para_0.05) =
colnames(PNSD_Fri_sd_ridge_tune_para_0.05) = colnames(PNSD_Sat_sd_ridge_tune_para_0.05) =
colnames(PNSD_Sun_sd_ridge_tune_para_0.05) = 3:24

rownames(PNSD_Mon_sd_ridge_tune_para_0.05) = rownames(PNSD_Tue_sd_ridge_tune_para_0.05) = 
rownames(PNSD_Wed_sd_ridge_tune_para_0.05) = rownames(PNSD_Thu_sd_ridge_tune_para_0.05) = 
rownames(PNSD_Fri_sd_ridge_tune_para_0.05) = rownames(PNSD_Sat_sd_ridge_tune_para_0.05) = 
rownames(PNSD_Sun_sd_ridge_tune_para_0.05) = x_grid

# objective

PNSD_Mon_sd_ridge_obj_0.05 = PNSD_Tue_sd_ridge_obj_0.05 = PNSD_Wed_sd_ridge_obj_0.05 = 
PNSD_Thu_sd_ridge_obj_0.05 = PNSD_Fri_sd_ridge_obj_0.05 = PNSD_Sat_sd_ridge_obj_0.05 = PNSD_Sun_sd_ridge_obj_0.05 = matrix(NA, 51, 22)
for(iw in 1:22)
{
    PNSD_Mon_sd_ridge_obj_0.05[,iw] = PNSD_Mon_sd_ridge_tune_0.05[[1]]$obj_val_min
    PNSD_Tue_sd_ridge_obj_0.05[,iw] = PNSD_Tue_sd_ridge_tune_0.05[[1]]$obj_val_min
    PNSD_Wed_sd_ridge_obj_0.05[,iw] = PNSD_Wed_sd_ridge_tune_0.05[[1]]$obj_val_min
    PNSD_Thu_sd_ridge_obj_0.05[,iw] = PNSD_Thu_sd_ridge_tune_0.05[[1]]$obj_val_min
    PNSD_Fri_sd_ridge_obj_0.05[,iw] = PNSD_Fri_sd_ridge_tune_0.05[[1]]$obj_val_min
    PNSD_Sat_sd_ridge_obj_0.05[,iw] = PNSD_Sat_sd_ridge_tune_0.05[[1]]$obj_val_min
    PNSD_Sun_sd_ridge_obj_0.05[,iw] = PNSD_Sun_sd_ridge_tune_0.05[[1]]$obj_val_min
    print(iw); rm(iw)
}

colnames(PNSD_Mon_sd_ridge_obj_0.05) = colnames(PNSD_Tue_sd_ridge_obj_0.05) = colnames(PNSD_Wed_sd_ridge_obj_0.05) = 
colnames(PNSD_Thu_sd_ridge_obj_0.05) = colnames(PNSD_Fri_sd_ridge_obj_0.05) = colnames(PNSD_Sat_sd_ridge_obj_0.05) = 
colnames(PNSD_Sun_sd_ridge_obj_0.05) = 3:24

rownames(PNSD_Mon_sd_ridge_obj_0.05) = rownames(PNSD_Tue_sd_ridge_obj_0.05) = rownames(PNSD_Wed_sd_ridge_obj_0.05) = 
rownames(PNSD_Thu_sd_ridge_obj_0.05) = rownames(PNSD_Fri_sd_ridge_obj_0.05) = rownames(PNSD_Sat_sd_ridge_obj_0.05) = 
rownames(PNSD_Sun_sd_ridge_obj_0.05) = x_grid


##################################
# sd approach for constructing PI
##################################
    
# data_set: 24 by T by 51
# fore_model: forecasting model
# tune_para_select: selected tuning parameters
# sd_val_mat: computed pointwise standard deviation matrix 51 by 24
# level_sig: level of significance 
    
sd_PI_ridge_test <- function(data_set, lambda_val, tune_para_select, sd_val_mat, observed_period, level_sig)
{
    n_hour = dim(data_set)[1]
    n_PM = dim(data_set)[3]
    holdout_set = array((10^(aperm(data_set[,((408*3/4+1):408),], c(3, 1, 2))) - 1)[,(1:24)[-(1:observed_period)],],
                        dim = c(n_PM, length((1:24)[-(1:observed_period)]), n_test))
    
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
    
    test_ridge_h1_lb = test_ridge_h1_ub = array(NA, dim = dim(ridge_h1_fore), dimnames = dimnames(ridge_h1_fore))
    for(iw in 1:n_PM)
    {
        test_ridge_h1_lb[iw,,] = ridge_h1_fore[iw,,] - tune_para_select[iw] * sd_val_mat[iw,]
        test_ridge_h1_ub[iw,,] = ridge_h1_fore[iw,,] + tune_para_select[iw] * sd_val_mat[iw,]
        rm(iw)
    }
      
    # compute CPD, interval score
      
    int_val = matrix(NA, n_PM, 3) 
    for(iw in 1:n_PM)
    {
        int_val[iw,] = interval_score(holdout = holdout_set[iw,,], lb = test_ridge_h1_lb[iw,,],
                                      ub = test_ridge_h1_ub[iw,,], alpha = level_sig)
        rm(iw)
    }
    colnames(int_val) = c("ECP", "CPD", "score")
    return(int_val)   
}

##################
# level_sig = 0.2
##################

registerDoMC(detectCores())
PNSD_Mon_sd_ridge_PI_test = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Mon_408, 
                                                                      lambda_val = ridge_lambda_log_Mon_para[iwk], 
                                                                      tune_para_select = PNSD_Mon_sd_ridge_tune_para[,iwk],
                                                                      sd_val_mat = PNSD_Mon_sd_ridge_tune[[iwk]]$sd_val, 
                                                                      observed_period = (iwk + 1), level_sig = 0.2) 

PNSD_Tue_sd_ridge_PI_test = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Tue_408, 
                                                                         lambda_val = ridge_lambda_log_Tue_para[iwk], 
                                                                         tune_para_select = PNSD_Tue_sd_ridge_tune_para[,iwk],
                                                                         sd_val_mat = PNSD_Tue_sd_ridge_tune[[iwk]]$sd_val, 
                                                                         observed_period = (iwk + 1), level_sig = 0.2) 

PNSD_Wed_sd_ridge_PI_test = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Wed_408, 
                                                                         lambda_val = ridge_lambda_log_Wed_para[iwk], 
                                                                         tune_para_select = PNSD_Wed_sd_ridge_tune_para[,iwk],
                                                                         sd_val_mat = PNSD_Wed_sd_ridge_tune[[iwk]]$sd_val, 
                                                                         observed_period = (iwk + 1), level_sig = 0.2) 

PNSD_Thu_sd_ridge_PI_test = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Thu_408, 
                                                                         lambda_val = ridge_lambda_log_Thu_para[iwk], 
                                                                         tune_para_select = PNSD_Thu_sd_ridge_tune_para[,iwk],
                                                                         sd_val_mat = PNSD_Thu_sd_ridge_tune[[iwk]]$sd_val, 
                                                                         observed_period = (iwk + 1), level_sig = 0.2) 

PNSD_Fri_sd_ridge_PI_test = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Fri_408, 
                                                                         lambda_val = ridge_lambda_log_Fri_para[iwk], 
                                                                         tune_para_select = PNSD_Fri_sd_ridge_tune_para[,iwk],
                                                                         sd_val_mat = PNSD_Fri_sd_ridge_tune[[iwk]]$sd_val, 
                                                                         observed_period = (iwk + 1), level_sig = 0.2) 

PNSD_Sat_sd_ridge_PI_test = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Sat_408, 
                                                                         lambda_val = ridge_lambda_log_Sat_para[iwk], 
                                                                         tune_para_select = PNSD_Sat_sd_ridge_tune_para[,iwk],
                                                                         sd_val_mat = PNSD_Sat_sd_ridge_tune[[iwk]]$sd_val, 
                                                                         observed_period = (iwk + 1), level_sig = 0.2) 

PNSD_Sun_sd_ridge_PI_test = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Sun_408, 
                                                                         lambda_val = ridge_lambda_log_Sun_para[iwk], 
                                                                         tune_para_select = PNSD_Sun_sd_ridge_tune_para[,iwk],
                                                                         sd_val_mat = PNSD_Sun_sd_ridge_tune[[iwk]]$sd_val, 
                                                                         observed_period = (iwk + 1), level_sig = 0.2) 

PNSD_Mon_sd_ridge_PI_test_array = PNSD_Tue_sd_ridge_PI_test_array = 
PNSD_Wed_sd_ridge_PI_test_array = PNSD_Thu_sd_ridge_PI_test_array = 
PNSD_Fri_sd_ridge_PI_test_array = PNSD_Sat_sd_ridge_PI_test_array = 
PNSD_Sun_sd_ridge_PI_test_array = array(NA, dim = c(51, 3, 22))
for(iwk in 1:22)
{
    PNSD_Mon_sd_ridge_PI_test_array[,,iwk] = PNSD_Mon_sd_ridge_PI_test[[iwk]]
    PNSD_Tue_sd_ridge_PI_test_array[,,iwk] = PNSD_Tue_sd_ridge_PI_test[[iwk]]
    PNSD_Wed_sd_ridge_PI_test_array[,,iwk] = PNSD_Wed_sd_ridge_PI_test[[iwk]]
    PNSD_Thu_sd_ridge_PI_test_array[,,iwk] = PNSD_Thu_sd_ridge_PI_test[[iwk]]
    PNSD_Fri_sd_ridge_PI_test_array[,,iwk] = PNSD_Fri_sd_ridge_PI_test[[iwk]]
    PNSD_Sat_sd_ridge_PI_test_array[,,iwk] = PNSD_Sat_sd_ridge_PI_test[[iwk]]
    PNSD_Sun_sd_ridge_PI_test_array[,,iwk] = PNSD_Sun_sd_ridge_PI_test[[iwk]]
    rm(iwk)
}

PNSD_days_sd_ridge_PI_test_array_CPD = rbind(apply(PNSD_Mon_sd_ridge_PI_test_array[,2,], 2, mean),
                                             apply(PNSD_Tue_sd_ridge_PI_test_array[,2,], 2, mean),
                                             apply(PNSD_Wed_sd_ridge_PI_test_array[,2,], 2, mean),
                                             apply(PNSD_Thu_sd_ridge_PI_test_array[,2,], 2, mean),
                                             apply(PNSD_Fri_sd_ridge_PI_test_array[,2,], 2, mean),
                                             apply(PNSD_Sat_sd_ridge_PI_test_array[,2,], 2, mean),
                                             apply(PNSD_Sun_sd_ridge_PI_test_array[,2,], 2, mean))

PNSD_days_sd_ridge_PI_test_array_score = rbind(apply(PNSD_Mon_sd_ridge_PI_test_array[,3,], 2, mean),
                                               apply(PNSD_Tue_sd_ridge_PI_test_array[,3,], 2, mean),
                                               apply(PNSD_Wed_sd_ridge_PI_test_array[,3,], 2, mean),
                                               apply(PNSD_Thu_sd_ridge_PI_test_array[,3,], 2, mean),
                                               apply(PNSD_Fri_sd_ridge_PI_test_array[,3,], 2, mean),
                                               apply(PNSD_Sat_sd_ridge_PI_test_array[,3,], 2, mean),
                                               apply(PNSD_Sun_sd_ridge_PI_test_array[,3,], 2, mean))
rownames(PNSD_days_sd_ridge_PI_test_array_CPD) = rownames(PNSD_days_sd_ridge_PI_test_array_score) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(PNSD_days_sd_ridge_PI_test_array_CPD) = colnames(PNSD_days_sd_ridge_PI_test_array_score) = 3:24

###################
# level_sig = 0.05
###################

registerDoMC(detectCores())
PNSD_Mon_sd_ridge_PI_test_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Mon_408, 
                                                                         lambda_val = ridge_lambda_log_Mon_para[iwk], 
                                                                         tune_para_select = PNSD_Mon_sd_ridge_tune_para_0.05[,iwk],
                                                                         sd_val_mat = PNSD_Mon_sd_ridge_tune[[iwk]]$sd_val, 
                                                                         observed_period = (iwk + 1), level_sig = 0.05) 

PNSD_Tue_sd_ridge_PI_test_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Tue_408, 
                                                                              lambda_val = ridge_lambda_log_Tue_para[iwk], 
                                                                              tune_para_select = PNSD_Tue_sd_ridge_tune_para_0.05[,iwk],
                                                                              sd_val_mat = PNSD_Tue_sd_ridge_tune[[iwk]]$sd_val, 
                                                                              observed_period = (iwk + 1), level_sig = 0.05) 

PNSD_Wed_sd_ridge_PI_test_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Wed_408, 
                                                                              lambda_val = ridge_lambda_log_Wed_para[iwk], 
                                                                              tune_para_select = PNSD_Wed_sd_ridge_tune_para_0.05[,iwk],
                                                                              sd_val_mat = PNSD_Wed_sd_ridge_tune[[iwk]]$sd_val, 
                                                                              observed_period = (iwk + 1), level_sig = 0.05) 

PNSD_Thu_sd_ridge_PI_test_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Thu_408, 
                                                                              lambda_val = ridge_lambda_log_Thu_para[iwk], 
                                                                              tune_para_select = PNSD_Thu_sd_ridge_tune_para_0.05[,iwk],
                                                                              sd_val_mat = PNSD_Thu_sd_ridge_tune[[iwk]]$sd_val, 
                                                                              observed_period = (iwk + 1), level_sig = 0.05) 

PNSD_Fri_sd_ridge_PI_test_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Fri_408, 
                                                                              lambda_val = ridge_lambda_log_Fri_para[iwk], 
                                                                              tune_para_select = PNSD_Fri_sd_ridge_tune_para_0.05[,iwk],
                                                                              sd_val_mat = PNSD_Fri_sd_ridge_tune[[iwk]]$sd_val, 
                                                                              observed_period = (iwk + 1), level_sig = 0.05) 

PNSD_Sat_sd_ridge_PI_test_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Sat_408, 
                                                                              lambda_val = ridge_lambda_log_Sat_para[iwk], 
                                                                              tune_para_select = PNSD_Sat_sd_ridge_tune_para_0.05[,iwk],
                                                                              sd_val_mat = PNSD_Sat_sd_ridge_tune[[iwk]]$sd_val, 
                                                                              observed_period = (iwk + 1), level_sig = 0.05) 

PNSD_Sun_sd_ridge_PI_test_0.05 = foreach(iwk = 1:22) %dopar% sd_PI_ridge_test(data_set = PNSD_log_Sun_408, 
                                                                              lambda_val = ridge_lambda_log_Sun_para[iwk], 
                                                                              tune_para_select = PNSD_Sun_sd_ridge_tune_para_0.05[,iwk],
                                                                              sd_val_mat = PNSD_Sun_sd_ridge_tune[[iwk]]$sd_val, 
                                                                              observed_period = (iwk + 1), level_sig = 0.05) 

PNSD_Mon_sd_ridge_PI_test_array_0.05 = PNSD_Tue_sd_ridge_PI_test_array_0.05 = 
PNSD_Wed_sd_ridge_PI_test_array_0.05 = PNSD_Thu_sd_ridge_PI_test_array_0.05 = 
PNSD_Fri_sd_ridge_PI_test_array_0.05 = PNSD_Sat_sd_ridge_PI_test_array_0.05 = 
PNSD_Sun_sd_ridge_PI_test_array_0.05 = array(NA, dim = c(51, 3, 22))
for(iwk in 1:22)
{
    PNSD_Mon_sd_ridge_PI_test_array_0.05[,,iwk] = PNSD_Mon_sd_ridge_PI_test_0.05[[iwk]]
    PNSD_Tue_sd_ridge_PI_test_array_0.05[,,iwk] = PNSD_Tue_sd_ridge_PI_test_0.05[[iwk]]
    PNSD_Wed_sd_ridge_PI_test_array_0.05[,,iwk] = PNSD_Wed_sd_ridge_PI_test_0.05[[iwk]]
    PNSD_Thu_sd_ridge_PI_test_array_0.05[,,iwk] = PNSD_Thu_sd_ridge_PI_test_0.05[[iwk]]
    PNSD_Fri_sd_ridge_PI_test_array_0.05[,,iwk] = PNSD_Fri_sd_ridge_PI_test_0.05[[iwk]]
    PNSD_Sat_sd_ridge_PI_test_array_0.05[,,iwk] = PNSD_Sat_sd_ridge_PI_test_0.05[[iwk]]
    PNSD_Sun_sd_ridge_PI_test_array_0.05[,,iwk] = PNSD_Sun_sd_ridge_PI_test_0.05[[iwk]]
    print(iwk);rm(iwk)
}

PNSD_days_sd_ridge_PI_test_array_CPD_0.05 = rbind(apply(PNSD_Mon_sd_ridge_PI_test_array_0.05[,2,], 2, mean),
                                                  apply(PNSD_Tue_sd_ridge_PI_test_array_0.05[,2,], 2, mean),
                                                  apply(PNSD_Wed_sd_ridge_PI_test_array_0.05[,2,], 2, mean),
                                                  apply(PNSD_Thu_sd_ridge_PI_test_array_0.05[,2,], 2, mean),
                                                  apply(PNSD_Fri_sd_ridge_PI_test_array_0.05[,2,], 2, mean),
                                                  apply(PNSD_Sat_sd_ridge_PI_test_array_0.05[,2,], 2, mean),
                                                  apply(PNSD_Sun_sd_ridge_PI_test_array_0.05[,2,], 2, mean))

PNSD_days_sd_ridge_PI_test_array_score_0.05 = rbind(apply(PNSD_Mon_sd_ridge_PI_test_array_0.05[,3,], 2, mean),
                                                    apply(PNSD_Tue_sd_ridge_PI_test_array_0.05[,3,], 2, mean),
                                                    apply(PNSD_Wed_sd_ridge_PI_test_array_0.05[,3,], 2, mean),
                                                    apply(PNSD_Thu_sd_ridge_PI_test_array_0.05[,3,], 2, mean),
                                                    apply(PNSD_Fri_sd_ridge_PI_test_array_0.05[,3,], 2, mean),
                                                    apply(PNSD_Sat_sd_ridge_PI_test_array_0.05[,3,], 2, mean),
                                                    apply(PNSD_Sun_sd_ridge_PI_test_array_0.05[,3,], 2, mean))

rownames(PNSD_days_sd_ridge_PI_test_array_CPD_0.05) = rownames(PNSD_days_sd_ridge_PI_test_array_score_0.05) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
colnames(PNSD_days_sd_ridge_PI_test_array_CPD_0.05) = colnames(PNSD_days_sd_ridge_PI_test_array_score_0.05) = 3:24

