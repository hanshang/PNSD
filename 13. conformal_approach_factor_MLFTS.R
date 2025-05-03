##############################
# level of significance = 0.2
##############################

factor_PNSD_Mon_conformal = conformal_PI(data_set = PNSD_log_Mon_408, fore_model = "factor_MLFTS", level_sig = 0.2)
factor_PNSD_Tue_conformal = conformal_PI(data_set = PNSD_log_Tue_408, fore_model = "factor_MLFTS", level_sig = 0.2)
factor_PNSD_Wed_conformal = conformal_PI(data_set = PNSD_log_Wed_408, fore_model = "factor_MLFTS", level_sig = 0.2)
factor_PNSD_Thu_conformal = conformal_PI(data_set = PNSD_log_Thu_408, fore_model = "factor_MLFTS", level_sig = 0.2)
factor_PNSD_Fri_conformal = conformal_PI(data_set = PNSD_log_Fri_408, fore_model = "factor_MLFTS", level_sig = 0.2)
factor_PNSD_Sat_conformal = conformal_PI(data_set = PNSD_log_Sat_408, fore_model = "factor_MLFTS", level_sig = 0.2)
factor_PNSD_Sun_conformal = conformal_PI(data_set = PNSD_log_Sun_408, fore_model = "factor_MLFTS", level_sig = 0.2)

# ECP, CPD, score

factor_PNSD_ECP_conformal_test = cbind(factor_PNSD_Mon_conformal[,1], factor_PNSD_Tue_conformal[,1],
                                       factor_PNSD_Wed_conformal[,1], factor_PNSD_Thu_conformal[,1], 
                                       factor_PNSD_Fri_conformal[,1], factor_PNSD_Sat_conformal[,1], 
                                       factor_PNSD_Sun_conformal[,1])

factor_PNSD_CPD_conformal_test = cbind(factor_PNSD_Mon_conformal[,2], factor_PNSD_Tue_conformal[,2],
                                       factor_PNSD_Wed_conformal[,2], factor_PNSD_Thu_conformal[,2], 
                                       factor_PNSD_Fri_conformal[,2], factor_PNSD_Sat_conformal[,2], 
                                       factor_PNSD_Sun_conformal[,2])

factor_PNSD_score_conformal_test = cbind(factor_PNSD_Mon_conformal[,3], factor_PNSD_Tue_conformal[,3],
                                         factor_PNSD_Wed_conformal[,3], factor_PNSD_Thu_conformal[,3], 
                                         factor_PNSD_Fri_conformal[,3], factor_PNSD_Sat_conformal[,3], 
                                         factor_PNSD_Sun_conformal[,3])
rownames(factor_PNSD_CPD_conformal_test) = rownames(factor_PNSD_score_conformal_test) = x_grid
colnames(factor_PNSD_CPD_conformal_test) = colnames(factor_PNSD_score_conformal_test) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

###############################
# level of significance = 0.05
###############################

factor_PNSD_Mon_conformal_0.95 = conformal_PI(data_set = PNSD_log_Mon_408, fore_model = "factor_MLFTS", level_sig = 0.05)
factor_PNSD_Tue_conformal_0.95 = conformal_PI(data_set = PNSD_log_Tue_408, fore_model = "factor_MLFTS", level_sig = 0.05)
factor_PNSD_Wed_conformal_0.95 = conformal_PI(data_set = PNSD_log_Wed_408, fore_model = "factor_MLFTS", level_sig = 0.05)
factor_PNSD_Thu_conformal_0.95 = conformal_PI(data_set = PNSD_log_Thu_408, fore_model = "factor_MLFTS", level_sig = 0.05)
factor_PNSD_Fri_conformal_0.95 = conformal_PI(data_set = PNSD_log_Fri_408, fore_model = "factor_MLFTS", level_sig = 0.05)
factor_PNSD_Sat_conformal_0.95 = conformal_PI(data_set = PNSD_log_Sat_408, fore_model = "factor_MLFTS", level_sig = 0.05)
factor_PNSD_Sun_conformal_0.95 = conformal_PI(data_set = PNSD_log_Sun_408, fore_model = "factor_MLFTS", level_sig = 0.05)

# ECP, CPD, score

factor_PNSD_ECP_conformal_test_0.95 = cbind(factor_PNSD_Mon_conformal_0.95[,1], factor_PNSD_Tue_conformal_0.95[,1],
                                            factor_PNSD_Wed_conformal_0.95[,1], factor_PNSD_Thu_conformal_0.95[,1], 
                                            factor_PNSD_Fri_conformal_0.95[,1], factor_PNSD_Sat_conformal_0.95[,1], 
                                            factor_PNSD_Sun_conformal_0.95[,1])

factor_PNSD_CPD_conformal_test_0.95 = cbind(factor_PNSD_Mon_conformal_0.95[,2], factor_PNSD_Tue_conformal_0.95[,2],
                                            factor_PNSD_Wed_conformal_0.95[,2], factor_PNSD_Thu_conformal_0.95[,2], 
                                            factor_PNSD_Fri_conformal_0.95[,2], factor_PNSD_Sat_conformal_0.95[,2], 
                                            factor_PNSD_Sun_conformal_0.95[,2])

factor_PNSD_score_conformal_test_0.95 = cbind(factor_PNSD_Mon_conformal_0.95[,3], factor_PNSD_Tue_conformal_0.95[,3],
                                              factor_PNSD_Wed_conformal_0.95[,3], factor_PNSD_Thu_conformal_0.95[,3], 
                                              factor_PNSD_Fri_conformal_0.95[,3], factor_PNSD_Sat_conformal_0.95[,3], 
                                              factor_PNSD_Sun_conformal_0.95[,3])
colnames(factor_PNSD_ECP_conformal_test_0.95) = colnames(factor_PNSD_CPD_conformal_test_0.95) = 
colnames(factor_PNSD_score_conformal_test_0.95) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
rownames(factor_PNSD_ECP_conformal_test_0.95) = rownames(factor_PNSD_CPD_conformal_test_0.95) = 
rownames(factor_PNSD_score_conformal_test_0.95) = x_grid

