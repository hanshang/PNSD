############
# ncp = 0.8
############

factor_PNSD_Mon_sd_tune = sd_PI(data_set = PNSD_log_Mon_408, fore_model = "factor_MLFTS", ncp = 0.8)
factor_PNSD_Tue_sd_tune = sd_PI(data_set = PNSD_log_Tue_408, fore_model = "factor_MLFTS", ncp = 0.8)
factor_PNSD_Wed_sd_tune = sd_PI(data_set = PNSD_log_Wed_408, fore_model = "factor_MLFTS", ncp = 0.8)
factor_PNSD_Thu_sd_tune = sd_PI(data_set = PNSD_log_Thu_408, fore_model = "factor_MLFTS", ncp = 0.8)
factor_PNSD_Fri_sd_tune = sd_PI(data_set = PNSD_log_Fri_408, fore_model = "factor_MLFTS", ncp = 0.8)
factor_PNSD_Sat_sd_tune = sd_PI(data_set = PNSD_log_Sat_408, fore_model = "factor_MLFTS", ncp = 0.8)
factor_PNSD_Sun_sd_tune = sd_PI(data_set = PNSD_log_Sun_408, fore_model = "factor_MLFTS", ncp = 0.8)

# objective value

round(factor_PNSD_Mon_sd_tune$obj_val_min, 4) 
round(factor_PNSD_Tue_sd_tune$obj_val_min, 4)
round(factor_PNSD_Wed_sd_tune$obj_val_min, 4)
round(factor_PNSD_Thu_sd_tune$obj_val_min, 4)
round(factor_PNSD_Fri_sd_tune$obj_val_min, 4)
round(factor_PNSD_Sat_sd_tune$obj_val_min, 4)
round(factor_PNSD_Sun_sd_tune$obj_val_min, 4)

# selected tuning parameter

round(factor_PNSD_Mon_sd_tune$tune_para_find, 4)
round(factor_PNSD_Tue_sd_tune$tune_para_find, 4)
round(factor_PNSD_Wed_sd_tune$tune_para_find, 4)
round(factor_PNSD_Thu_sd_tune$tune_para_find, 4)
round(factor_PNSD_Fri_sd_tune$tune_para_find, 4)
round(factor_PNSD_Sat_sd_tune$tune_para_find, 4)
round(factor_PNSD_Sun_sd_tune$tune_para_find, 4)

#############
# ncp = 0.95
#############

factor_PNSD_Mon_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Mon_408, fore_model = "factor_MLFTS", ncp = 0.95)
factor_PNSD_Tue_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Tue_408, fore_model = "factor_MLFTS", ncp = 0.95)
factor_PNSD_Wed_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Wed_408, fore_model = "factor_MLFTS", ncp = 0.95)
factor_PNSD_Thu_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Thu_408, fore_model = "factor_MLFTS", ncp = 0.95)
factor_PNSD_Fri_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Fri_408, fore_model = "factor_MLFTS", ncp = 0.95)
factor_PNSD_Sat_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Sat_408, fore_model = "factor_MLFTS", ncp = 0.95)
factor_PNSD_Sun_sd_tune_0.95 = sd_PI(data_set = PNSD_log_Sun_408, fore_model = "factor_MLFTS", ncp = 0.95)

# objective value

round(factor_PNSD_Mon_sd_tune_0.95$obj_val_min, 4) 
round(factor_PNSD_Tue_sd_tune_0.95$obj_val_min, 4)
round(factor_PNSD_Wed_sd_tune_0.95$obj_val_min, 4)
round(factor_PNSD_Thu_sd_tune_0.95$obj_val_min, 4)
round(factor_PNSD_Fri_sd_tune_0.95$obj_val_min, 4)
round(factor_PNSD_Sat_sd_tune_0.95$obj_val_min, 4)
round(factor_PNSD_Sun_sd_tune_0.95$obj_val_min, 4)

# selected tuning parameter

round(factor_PNSD_Mon_sd_tune_0.95$tune_para_find, 4)
round(factor_PNSD_Tue_sd_tune_0.95$tune_para_find, 4)
round(factor_PNSD_Wed_sd_tune_0.95$tune_para_find, 4)
round(factor_PNSD_Thu_sd_tune_0.95$tune_para_find, 4)
round(factor_PNSD_Fri_sd_tune_0.95$tune_para_find, 4)
round(factor_PNSD_Sat_sd_tune_0.95$tune_para_find, 4)
round(factor_PNSD_Sun_sd_tune_0.95$tune_para_find, 4)

###################
## level_sig = 0.2
###################

factor_PNSD_Mon_sd_test = sd_PI_test(data_set = PNSD_log_Mon_408, fore_model = "factor_MLFTS",
                                     tune_para_select = factor_PNSD_Mon_sd_tune$tune_para_find,
                                     sd_val_mat = factor_PNSD_Mon_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Tue_sd_test = sd_PI_test(data_set = PNSD_log_Tue_408, fore_model = "factor_MLFTS",
                                     tune_para_select = factor_PNSD_Tue_sd_tune$tune_para_find,
                                     sd_val_mat = factor_PNSD_Tue_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Wed_sd_test = sd_PI_test(data_set = PNSD_log_Wed_408, fore_model = "factor_MLFTS",
                                     tune_para_select = factor_PNSD_Wed_sd_tune$tune_para_find,
                                     sd_val_mat = factor_PNSD_Wed_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Thu_sd_test = sd_PI_test(data_set = PNSD_log_Thu_408, fore_model = "factor_MLFTS",
                                     tune_para_select = factor_PNSD_Thu_sd_tune$tune_para_find,
                                     sd_val_mat = factor_PNSD_Thu_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Fri_sd_test = sd_PI_test(data_set = PNSD_log_Fri_408, fore_model = "factor_MLFTS",
                                     tune_para_select = factor_PNSD_Fri_sd_tune$tune_para_find,
                                     sd_val_mat = factor_PNSD_Fri_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Sat_sd_test = sd_PI_test(data_set = PNSD_log_Sat_408, fore_model = "factor_MLFTS",
                                     tune_para_select = factor_PNSD_Sat_sd_tune$tune_para_find,
                                     sd_val_mat = factor_PNSD_Sat_sd_tune$sd_val, level_sig = 0.2)

factor_PNSD_Sun_sd_test = sd_PI_test(data_set = PNSD_log_Sun_408, fore_model = "factor_MLFTS",
                                     tune_para_select = factor_PNSD_Sun_sd_tune$tune_para_find,
                                     sd_val_mat = factor_PNSD_Sun_sd_tune$sd_val, level_sig = 0.2)

# empirical coverage probability (ECP), coverage probability difference (CPD), interval score

factor_PNSD_ECP_sd_test = cbind(factor_PNSD_Mon_sd_test[,1], factor_PNSD_Tue_sd_test[,1],
                                factor_PNSD_Wed_sd_test[,1], factor_PNSD_Thu_sd_test[,1],
                                factor_PNSD_Fri_sd_test[,1], factor_PNSD_Sat_sd_test[,1],
                                factor_PNSD_Sun_sd_test[,1])

factor_PNSD_CPD_sd_test = cbind(factor_PNSD_Mon_sd_test[,2], factor_PNSD_Tue_sd_test[,2],
                                factor_PNSD_Wed_sd_test[,2], factor_PNSD_Thu_sd_test[,2],
                                factor_PNSD_Fri_sd_test[,2], factor_PNSD_Sat_sd_test[,2],
                                factor_PNSD_Sun_sd_test[,2])

factor_PNSD_score_sd_test = cbind(factor_PNSD_Mon_sd_test[,3], factor_PNSD_Tue_sd_test[,3],
                                  factor_PNSD_Wed_sd_test[,3], factor_PNSD_Thu_sd_test[,3],
                                  factor_PNSD_Fri_sd_test[,3], factor_PNSD_Sat_sd_test[,3],
                                  factor_PNSD_Sun_sd_test[,3])
rownames(factor_PNSD_CPD_sd_test) = rownames(factor_PNSD_score_sd_test) = x_grid
colnames(factor_PNSD_CPD_sd_test) = colnames(factor_PNSD_score_sd_test) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

####################
## level_sig = 0.05
####################

factor_PNSD_Mon_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Mon_408, fore_model = "factor_MLFTS",
                                          tune_para_select = factor_PNSD_Mon_sd_tune_0.95$tune_para_find,
                                          sd_val_mat = factor_PNSD_Mon_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Tue_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Tue_408, fore_model = "factor_MLFTS",
                                          tune_para_select = factor_PNSD_Tue_sd_tune_0.95$tune_para_find,
                                          sd_val_mat = factor_PNSD_Tue_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Wed_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Wed_408, fore_model = "factor_MLFTS",
                                          tune_para_select = factor_PNSD_Wed_sd_tune_0.95$tune_para_find,
                                          sd_val_mat = factor_PNSD_Wed_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Thu_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Thu_408, fore_model = "factor_MLFTS",
                                          tune_para_select = factor_PNSD_Thu_sd_tune_0.95$tune_para_find,
                                          sd_val_mat = factor_PNSD_Thu_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Fri_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Fri_408, fore_model = "factor_MLFTS",
                                          tune_para_select = factor_PNSD_Fri_sd_tune_0.95$tune_para_find,
                                          sd_val_mat = factor_PNSD_Fri_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Sat_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Sat_408, fore_model = "factor_MLFTS",
                                          tune_para_select = factor_PNSD_Sat_sd_tune_0.95$tune_para_find,
                                          sd_val_mat = factor_PNSD_Sat_sd_tune_0.95$sd_val, level_sig = 0.05)

factor_PNSD_Sun_sd_test_0.95 = sd_PI_test(data_set = PNSD_log_Sun_408, fore_model = "factor_MLFTS",
                                          tune_para_select = factor_PNSD_Sun_sd_tune_0.95$tune_para_find,
                                          sd_val_mat = factor_PNSD_Sun_sd_tune_0.95$sd_val, level_sig = 0.05)

# empirical coverage probability (ECP), coverage probability difference (CPD), interval score

factor_PNSD_ECP_sd_test_0.95 = cbind(factor_PNSD_Mon_sd_test_0.95[,1], factor_PNSD_Tue_sd_test_0.95[,1],
                                     factor_PNSD_Wed_sd_test_0.95[,1], factor_PNSD_Thu_sd_test_0.95[,1],
                                     factor_PNSD_Fri_sd_test_0.95[,1], factor_PNSD_Sat_sd_test_0.95[,1],
                                     factor_PNSD_Sun_sd_test_0.95[,1])

factor_PNSD_CPD_sd_test_0.95 = cbind(factor_PNSD_Mon_sd_test_0.95[,2], factor_PNSD_Tue_sd_test_0.95[,2],
                                     factor_PNSD_Wed_sd_test_0.95[,2], factor_PNSD_Thu_sd_test_0.95[,2],
                                     factor_PNSD_Fri_sd_test_0.95[,2], factor_PNSD_Sat_sd_test_0.95[,2],
                                     factor_PNSD_Sun_sd_test_0.95[,2])

factor_PNSD_score_sd_test_0.95 = cbind(factor_PNSD_Mon_sd_test_0.95[,3], factor_PNSD_Tue_sd_test_0.95[,3],
                                       factor_PNSD_Wed_sd_test_0.95[,3], factor_PNSD_Thu_sd_test_0.95[,3],
                                       factor_PNSD_Fri_sd_test_0.95[,3], factor_PNSD_Sat_sd_test_0.95[,3],
                                       factor_PNSD_Sun_sd_test_0.95[,3])
rownames(factor_PNSD_CPD_sd_test_0.95) = rownames(factor_PNSD_score_sd_test_0.95) = x_grid
colnames(factor_PNSD_CPD_sd_test_0.95) = colnames(factor_PNSD_score_sd_test_0.95) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

