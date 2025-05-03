# tune_para: tuning parameter
# resi_mat: residual matrix
# sd_val_input: standard deviation
# alpha_level: nominal coverage probability

source("load_packages.R")

# objective function to be optimalized

tune_para_find_function <- function(tune_para, resi_mat, sd_val_input, alpha_level)
{
    n_age = nrow(resi_mat)
    ind = matrix(NA, n_age, ncol(resi_mat))
    for(iw in 1:ncol(resi_mat))
    {
        ind[,iw] = ifelse(between(resi_mat[,iw], -tune_para * sd_val_input, tune_para * sd_val_input), 1, 0)
        rm(iw)
    }
    ecp = sum(ind)/(n_age * ncol(resi_mat))
    rm(ind)
    return(abs(ecp - alpha_level))
}

# data_set: 51 x 24 x n_test
# fore_model: MLFTS alone or factor + MLFTS model
# ncp: nominal coverage probability

sd_PI <- function(data_set, fore_model, ncp)
{
    n_hour = dim(data_set)[1]
    n_days = dim(data_set)[2]
    n_PM = dim(data_set)[3]
    
    # holdout data in the original scale
    
    holdout_set = 10^(aperm(data_set[,((n_days*1/2+1):(n_days*3/4)),], c(3, 1, 2))) - 1
    
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
        rm(iw); rm(n_train)
    }
    
    # difference between holdout data and forecasts (51 x 24 x n_test)
    
    resi_mat = holdout_set - (10^validation_MLFTS_fit_h1 - 1)
    
    # compute pointwise standard deviation
    
    sd_val = matrix(NA, n_PM, n_hour)
    for(iw in 1:n_PM)
    {
        sd_val[iw,] = apply(resi_mat[iw,,], 1, sd)
        rm(iw)
    }

    # find the optimal tuning parameter
    
    tune_para_find = obj_val_min = vector("numeric", n_PM)
    for(iw in 1:n_PM)
    {
        tune_para_find_val_1 = optimise(f = tune_para_find_function, interval = c(0, 1),
                                        resi_mat = resi_mat[iw,,], sd_val_input = sd_val[iw,],
                                        alpha_level = ncp)
    
        tune_para_find_val_2 = optimise(f = tune_para_find_function, interval = c(0, 5),
                                        resi_mat = resi_mat[iw,,], sd_val_input = sd_val[iw,],
                                        alpha_level = ncp)
        
        tune_para_find_val_3 = optimise(f = tune_para_find_function, interval = c(0, 10),
                                        resi_mat = resi_mat[iw,,], sd_val_input = sd_val[iw,],
                                        alpha_level = ncp)
        
        tune_para_find_val_4 = optimise(f = tune_para_find_function, interval = c(0, 20),
                                        resi_mat = resi_mat[iw,,], sd_val_input = sd_val[iw,],
                                        alpha_level = ncp)
        
        tune_para_find_val_5 = optim(par = 1, fn = tune_para_find_function, lower = 0, method = "L-BFGS-B",
                                       resi_mat = resi_mat[iw,,], sd_val_input = sd_val[iw,],
                                       alpha_level = ncp)
        
        tune_para_find_val_6 = optim(par = 1, fn = tune_para_find_function,
                                       resi_mat = resi_mat[iw,,], sd_val_input = sd_val[iw,],
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

############
# ncp = 0.8
############

PNSD_Mon_sd_tune = sd_PI(data_set = PNSD_log_Mon_408, fore_model = "MLFTS", ncp = 0.8)
PNSD_Tue_sd_tune = sd_PI(data_set = PNSD_log_Tue_408, fore_model = "MLFTS", ncp = 0.8)
PNSD_Wed_sd_tune = sd_PI(data_set = PNSD_log_Wed_408, fore_model = "MLFTS", ncp = 0.8)
PNSD_Thu_sd_tune = sd_PI(data_set = PNSD_log_Thu_408, fore_model = "MLFTS", ncp = 0.8)
PNSD_Fri_sd_tune = sd_PI(data_set = PNSD_log_Fri_408, fore_model = "MLFTS", ncp = 0.8)
PNSD_Sat_sd_tune = sd_PI(data_set = PNSD_log_Sat_408, fore_model = "MLFTS", ncp = 0.8)
PNSD_Sun_sd_tune = sd_PI(data_set = PNSD_log_Sun_408, fore_model = "MLFTS", ncp = 0.8)

round(PNSD_Mon_sd_tune$obj_val_min, 4) 
round(PNSD_Tue_sd_tune$obj_val_min, 4)
round(PNSD_Wed_sd_tune$obj_val_min, 4)
round(PNSD_Thu_sd_tune$obj_val_min, 4)
round(PNSD_Fri_sd_tune$obj_val_min, 4)
round(PNSD_Sat_sd_tune$obj_val_min, 4)
round(PNSD_Sun_sd_tune$obj_val_min, 4)

round(PNSD_Mon_sd_tune$tune_para_find, 4)
round(PNSD_Tue_sd_tune$tune_para_find, 4)
round(PNSD_Wed_sd_tune$tune_para_find, 4)
round(PNSD_Thu_sd_tune$tune_para_find, 4)
round(PNSD_Fri_sd_tune$tune_para_find, 4)
round(PNSD_Sat_sd_tune$tune_para_find, 4)
round(PNSD_Sun_sd_tune$tune_para_find, 4)

PNSD_days_sd_tune = cbind(PNSD_Mon_sd_tune$tune_para_find,
                          PNSD_Tue_sd_tune$tune_para_find,
                          PNSD_Wed_sd_tune$tune_para_find,
                          PNSD_Thu_sd_tune$tune_para_find,
                          PNSD_Fri_sd_tune$tune_para_find,
                          PNSD_Sat_sd_tune$tune_para_find,
                          PNSD_Sun_sd_tune$tune_para_find)
rownames(PNSD_days_sd_tune) = x_grid
colnames(PNSD_days_sd_tune) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

# save figure

savefig("Fig_4a", width = 12, height = 10, toplines = 0.8)
plot(x_grid, PNSD_Mon_sd_tune$tune_para_find, col = 1, pch = 1, ylim = range(PNSD_days_sd_tune), 
     ylab = "Tuning parameter", xlab = "Particle size")
points(x_grid, PNSD_Tue_sd_tune$tune_para_find, col = 2, pch = 2)
points(x_grid, PNSD_Wed_sd_tune$tune_para_find, col = 3, pch = 3)
points(x_grid, PNSD_Thu_sd_tune$tune_para_find, col = 4, pch = 4)
points(x_grid, PNSD_Fri_sd_tune$tune_para_find, col = 5, pch = 5)
points(x_grid, PNSD_Sat_sd_tune$tune_para_find, col = 6, pch = 6)
points(x_grid, PNSD_Sun_sd_tune$tune_para_find, col = 7, pch = 7)
legend("topright", c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"), col = 1:7, pch = 1:7, cex = 0.8, ncol = 2)
dev.off()

#############
# ncp = 0.95
#############

PNSD_Mon_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Mon_408, fore_model = "MLFTS", ncp = 0.95)
PNSD_Tue_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Tue_408, fore_model = "MLFTS", ncp = 0.95)
PNSD_Wed_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Wed_408, fore_model = "MLFTS", ncp = 0.95)
PNSD_Thu_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Thu_408, fore_model = "MLFTS", ncp = 0.95)
PNSD_Fri_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Fri_408, fore_model = "MLFTS", ncp = 0.95)
PNSD_Sat_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Sat_408, fore_model = "MLFTS", ncp = 0.95)
PNSD_Sun_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Sun_408, fore_model = "MLFTS", ncp = 0.95)

round(PNSD_Mon_sd_tune_0.95$obj_val_min, 4) 
round(PNSD_Tue_sd_tune_0.95$obj_val_min, 4)
round(PNSD_Wed_sd_tune_0.95$obj_val_min, 4)
round(PNSD_Thu_sd_tune_0.95$obj_val_min, 4)
round(PNSD_Fri_sd_tune_0.95$obj_val_min, 4)
round(PNSD_Sat_sd_tune_0.95$obj_val_min, 4)
round(PNSD_Sun_sd_tune_0.95$obj_val_min, 4)

round(PNSD_Mon_sd_tune_0.95$tune_para_find, 4)
round(PNSD_Tue_sd_tune_0.95$tune_para_find, 4)
round(PNSD_Wed_sd_tune_0.95$tune_para_find, 4)
round(PNSD_Thu_sd_tune_0.95$tune_para_find, 4)
round(PNSD_Fri_sd_tune_0.95$tune_para_find, 4)
round(PNSD_Sat_sd_tune_0.95$tune_para_find, 4)
round(PNSD_Sun_sd_tune_0.95$tune_para_find, 4)

PNSD_days_sd_tune_0.95 = cbind(PNSD_Mon_sd_tune_0.95$tune_para_find,
                               PNSD_Tue_sd_tune_0.95$tune_para_find,
                               PNSD_Wed_sd_tune_0.95$tune_para_find,
                               PNSD_Thu_sd_tune_0.95$tune_para_find,
                               PNSD_Fri_sd_tune_0.95$tune_para_find,
                               PNSD_Sat_sd_tune_0.95$tune_para_find,
                               PNSD_Sun_sd_tune_0.95$tune_para_find)
rownames(PNSD_days_sd_tune_0.95) = x_grid
colnames(PNSD_days_sd_tune_0.95) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

# save figure

savefig("Fig_4b", width = 12, height = 10, toplines = 0.8)
plot(x_grid, PNSD_Mon_sd_tune_0.95$tune_para_find, xlab = "", ylab = "", type = "p", ylim = range(PNSD_days_sd_tune_0.95))
points(x_grid, PNSD_Tue_sd_tune_0.95$tune_para_find, col = 2, pch = 2)
points(x_grid, PNSD_Wed_sd_tune_0.95$tune_para_find, col = 3, pch = 3)
points(x_grid, PNSD_Thu_sd_tune_0.95$tune_para_find, col = 4, pch = 4)
points(x_grid, PNSD_Fri_sd_tune_0.95$tune_para_find, col = 5, pch = 5)
points(x_grid, PNSD_Sat_sd_tune_0.95$tune_para_find, col = 6, pch = 6)
points(x_grid, PNSD_Sun_sd_tune_0.95$tune_para_find, col = 7, pch = 7)
dev.off()

##################################
# sd approach for constructing PI
##################################

# data_set: 24 by T by 51
# fore_model: forecasting model
# tune_para_select: selected tuning parameters
# sd_val_mat: computed pointwise standard deviation matrix 51 by 24
# level_sig: level of significance 

sd_PI_test <- function(data_set, fore_model, tune_para_select, sd_val_mat, level_sig)
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

    test_MLFTS_fit_h1_lb = test_MLFTS_fit_h1_ub = array(NA, dim = dim(test_MLFTS_fit_h1), dimnames = dimnames(test_MLFTS_fit_h1))
    for(iw in 1:n_PM)
    {
        test_MLFTS_fit_h1_lb[iw,,] = (10^(test_MLFTS_fit_h1[iw,,]) - 1) - tune_para_select[iw] * sd_val_mat[iw,]
        test_MLFTS_fit_h1_ub[iw,,] = (10^(test_MLFTS_fit_h1[iw,,]) - 1) + tune_para_select[iw] * sd_val_mat[iw,]
        rm(iw)
    }

    # compute empirical coverage probability (ECP), coverage probability difference (CPD), interval score
    
    int_val = matrix(NA, n_PM, 3) 
    for(iw in 1:n_PM)
    {
        int_val[iw,] = interval_score(holdout = holdout_set[iw,,], lb = test_MLFTS_fit_h1_lb[iw,,],
                                      ub = test_MLFTS_fit_h1_ub[iw,,], alpha = level_sig)
        rm(iw)
    }
    colnames(int_val) = c("ECP", "CPD", "score")
    return(int_val)   
}

###################
## level_sig = 0.2
###################

PNSD_Mon_sd_test = sd_PI_test(data_set = PNSD_log_Mon_408, fore_model = "MLFTS",
                              tune_para_select = PNSD_Mon_sd_tune$tune_para_find,
                              sd_val_mat = PNSD_Mon_sd_tune$sd_val, level_sig = 0.2)

PNSD_Tue_sd_test = sd_PI_test(data_set = PNSD_log_Tue_408, fore_model = "MLFTS",
                              tune_para_select = PNSD_Tue_sd_tune$tune_para_find,
                              sd_val_mat = PNSD_Tue_sd_tune$sd_val, level_sig = 0.2)

PNSD_Wed_sd_test = sd_PI_test(data_set = PNSD_log_Wed_408, fore_model = "MLFTS",
                              tune_para_select = PNSD_Wed_sd_tune$tune_para_find,
                              sd_val_mat = PNSD_Wed_sd_tune$sd_val, level_sig = 0.2)

PNSD_Thu_sd_test = sd_PI_test(data_set = PNSD_log_Thu_408, fore_model = "MLFTS",
                              tune_para_select = PNSD_Thu_sd_tune$tune_para_find,
                              sd_val_mat = PNSD_Thu_sd_tune$sd_val, level_sig = 0.2)

PNSD_Fri_sd_test = sd_PI_test(data_set = PNSD_log_Fri_408, fore_model = "MLFTS",
                              tune_para_select = PNSD_Fri_sd_tune$tune_para_find,
                              sd_val_mat = PNSD_Fri_sd_tune$sd_val, level_sig = 0.2)

PNSD_Sat_sd_test = sd_PI_test(data_set = PNSD_log_Sat_408, fore_model = "MLFTS",
                              tune_para_select = PNSD_Sat_sd_tune$tune_para_find,
                              sd_val_mat = PNSD_Sat_sd_tune$sd_val, level_sig = 0.2)

PNSD_Sun_sd_test = sd_PI_test(data_set = PNSD_log_Sun_408, fore_model = "MLFTS",
                              tune_para_select = PNSD_Sun_sd_tune$tune_para_find,
                              sd_val_mat = PNSD_Sun_sd_tune$sd_val, level_sig = 0.2)

# collect ECP, CPD and score

PNSD_ECP_sd_test = cbind(PNSD_Mon_sd_test[,1], PNSD_Tue_sd_test[,1], PNSD_Wed_sd_test[,1],
                         PNSD_Thu_sd_test[,1], PNSD_Fri_sd_test[,1], PNSD_Sat_sd_test[,1], 
                         PNSD_Sun_sd_test[,1])

PNSD_CPD_sd_test = cbind(PNSD_Mon_sd_test[,2], PNSD_Tue_sd_test[,2], PNSD_Wed_sd_test[,2],
                         PNSD_Thu_sd_test[,2], PNSD_Fri_sd_test[,2], PNSD_Sat_sd_test[,2], 
                         PNSD_Sun_sd_test[,2])

PNSD_score_sd_test = cbind(PNSD_Mon_sd_test[,3], PNSD_Tue_sd_test[,3], PNSD_Wed_sd_test[,3],
                           PNSD_Thu_sd_test[,3], PNSD_Fri_sd_test[,3], PNSD_Sat_sd_test[,3], 
                           PNSD_Sun_sd_test[,3])
colnames(PNSD_ECP_sd_test) = colnames(PNSD_CPD_sd_test) = colnames(PNSD_score_sd_test) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
rownames(PNSD_ECP_sd_test) = rownames(PNSD_CPD_sd_test) = rownames(PNSD_score_sd_test) = x_grid

# graphical displays

boxplot(PNSD_ECP_sd_test,   notch = TRUE)
boxplot(PNSD_CPD_sd_test,   notch = TRUE)
boxplot(PNSD_score_sd_test, notch = TRUE)

####################
## level_sig = 0.05
####################

PNSD_Mon_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Mon_408, fore_model = "MLFTS",
                              tune_para_select = PNSD_Mon_sd_tune_0.95$tune_para_find,
                              sd_val_mat = PNSD_Mon_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Tue_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Tue_408, fore_model = "MLFTS",
                                   tune_para_select = PNSD_Tue_sd_tune_0.95$tune_para_find,
                                   sd_val_mat = PNSD_Tue_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Wed_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Wed_408, fore_model = "MLFTS",
                                   tune_para_select = PNSD_Wed_sd_tune_0.95$tune_para_find,
                                   sd_val_mat = PNSD_Wed_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Thu_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Thu_408, fore_model = "MLFTS",
                                   tune_para_select = PNSD_Thu_sd_tune_0.95$tune_para_find,
                                   sd_val_mat = PNSD_Thu_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Fri_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Fri_408, fore_model = "MLFTS",
                                   tune_para_select = PNSD_Fri_sd_tune_0.95$tune_para_find,
                                   sd_val_mat = PNSD_Fri_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Sat_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Sat_408, fore_model = "MLFTS",
                                   tune_para_select = PNSD_Sat_sd_tune_0.95$tune_para_find,
                                   sd_val_mat = PNSD_Sat_sd_tune_0.95$sd_val, level_sig = 0.05)

PNSD_Sun_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Sun_408, fore_model = "MLFTS",
                                   tune_para_select = PNSD_Sun_sd_tune_0.95$tune_para_find,
                                   sd_val_mat = PNSD_Sun_sd_tune_0.95$sd_val, level_sig = 0.05)

# collect ECP, CPD and score

PNSD_ECP_sd_test_0.95 = cbind(PNSD_Mon_sd_test_0.95[,1], PNSD_Tue_sd_test_0.95[,1], PNSD_Wed_sd_test_0.95[,1],
                              PNSD_Thu_sd_test_0.95[,1], PNSD_Fri_sd_test_0.95[,1], PNSD_Sat_sd_test_0.95[,1], 
                              PNSD_Sun_sd_test_0.95[,1])

PNSD_CPD_sd_test_0.95 = cbind(PNSD_Mon_sd_test_0.95[,2], PNSD_Tue_sd_test_0.95[,2], PNSD_Wed_sd_test_0.95[,2],
                              PNSD_Thu_sd_test_0.95[,2], PNSD_Fri_sd_test_0.95[,2], PNSD_Sat_sd_test_0.95[,2], 
                              PNSD_Sun_sd_test_0.95[,2])

PNSD_score_sd_test_0.95 = cbind(PNSD_Mon_sd_test_0.95[,3], PNSD_Tue_sd_test_0.95[,3], PNSD_Wed_sd_test_0.95[,3],
                                PNSD_Thu_sd_test_0.95[,3], PNSD_Fri_sd_test_0.95[,3], PNSD_Sat_sd_test_0.95[,3], 
                                PNSD_Sun_sd_test_0.95[,3])
colnames(PNSD_ECP_sd_test_0.95) = colnames(PNSD_CPD_sd_test_0.95) = colnames(PNSD_score_sd_test_0.95) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
rownames(PNSD_ECP_sd_test_0.95) = rownames(PNSD_CPD_sd_test_0.95) = rownames(PNSD_score_sd_test_0.95) = x_grid

