#############
# evaluation 
#############

registerDoMC(10)
factor_BM_log_Mon_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Mon_408, newdata_period = 1:iwk, fore_method = "factor_MLFTS")
factor_BM_log_Tue_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Tue_408, newdata_period = 1:iwk, fore_method = "factor_MLFTS")
factor_BM_log_Wed_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Wed_408, newdata_period = 1:iwk, fore_method = "factor_MLFTS")
factor_BM_log_Thu_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Thu_408, newdata_period = 1:iwk, fore_method = "factor_MLFTS")
factor_BM_log_Fri_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Fri_408, newdata_period = 1:iwk, fore_method = "factor_MLFTS")
factor_BM_log_Sat_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Sat_408, newdata_period = 1:iwk, fore_method = "factor_MLFTS")
factor_BM_log_Sun_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Sun_408, newdata_period = 1:iwk, fore_method = "factor_MLFTS")

factor_BM_log_Mon = factor_BM_log_Tue = factor_BM_log_Wed = factor_BM_log_Thu = factor_BM_log_Fri = 
factor_BM_log_Sat = factor_BM_log_Sun = array(NA, dim = c(22, 51, n_test), dimnames = list(1:22, x_grid, 1:n_test))
for(iwk in 1:22)
{
    factor_BM_log_Mon[iwk,,] = factor_BM_log_Mon_list[[iwk]]
    factor_BM_log_Tue[iwk,,] = factor_BM_log_Tue_list[[iwk]]
    factor_BM_log_Wed[iwk,,] = factor_BM_log_Wed_list[[iwk]]
    factor_BM_log_Thu[iwk,,] = factor_BM_log_Thu_list[[iwk]]
    factor_BM_log_Fri[iwk,,] = factor_BM_log_Fri_list[[iwk]]
    factor_BM_log_Sat[iwk,,] = factor_BM_log_Sat_list[[iwk]]
    factor_BM_log_Sun[iwk,,] = factor_BM_log_Sun_list[[iwk]]
    rm(iwk)
}

# take averages

factor_BM_log_Mon_mean = apply(factor_BM_log_Mon, c(1, 2), mean)
factor_BM_log_Tue_mean = apply(factor_BM_log_Tue, c(1, 2), mean)
factor_BM_log_Wed_mean = apply(factor_BM_log_Wed, c(1, 2), mean)
factor_BM_log_Thu_mean = apply(factor_BM_log_Thu, c(1, 2), mean)
factor_BM_log_Fri_mean = apply(factor_BM_log_Fri, c(1, 2), mean)
factor_BM_log_Sat_mean = apply(factor_BM_log_Sat, c(1, 2), mean)
factor_BM_log_Sun_mean = apply(factor_BM_log_Sun, c(1, 2), mean)

# by hour

factor_BM_log_Mon_mean_hour = apply(factor_BM_log_Mon_mean, 1, mean)
factor_BM_log_Tue_mean_hour = apply(factor_BM_log_Tue_mean, 1, mean)
factor_BM_log_Wed_mean_hour = apply(factor_BM_log_Wed_mean, 1, mean)
factor_BM_log_Thu_mean_hour = apply(factor_BM_log_Thu_mean, 1, mean)
factor_BM_log_Fri_mean_hour = apply(factor_BM_log_Fri_mean, 1, mean)
factor_BM_log_Sat_mean_hour = apply(factor_BM_log_Sat_mean, 1, mean)
factor_BM_log_Sun_mean_hour = apply(factor_BM_log_Sun_mean, 1, mean)

factor_BM_days_mean_hour = rbind(factor_BM_log_Mon_mean_hour,
                                 factor_BM_log_Tue_mean_hour,
                                 factor_BM_log_Wed_mean_hour,
                                 factor_BM_log_Thu_mean_hour,
                                 factor_BM_log_Fri_mean_hour,
                                 factor_BM_log_Sat_mean_hour,
                                 factor_BM_log_Sun_mean_hour)
colnames(factor_BM_days_mean_hour) = 3:24
rownames(factor_BM_days_mean_hour) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun") 

# by size

factor_BM_log_Mon_mean_size = apply(factor_BM_log_Mon_mean, 2, mean)
factor_BM_log_Tue_mean_size = apply(factor_BM_log_Tue_mean, 2, mean)
factor_BM_log_Wed_mean_size = apply(factor_BM_log_Wed_mean, 2, mean)
factor_BM_log_Thu_mean_size = apply(factor_BM_log_Thu_mean, 2, mean)
factor_BM_log_Fri_mean_size = apply(factor_BM_log_Fri_mean, 2, mean)
factor_BM_log_Sat_mean_size = apply(factor_BM_log_Sat_mean, 2, mean)
factor_BM_log_Sun_mean_size = apply(factor_BM_log_Sun_mean, 2, mean)

factor_BM_days_mean_size = rbind(factor_BM_log_Mon_mean_size,
                                 factor_BM_log_Tue_mean_size,
                                 factor_BM_log_Wed_mean_size,
                                 factor_BM_log_Thu_mean_size,
                                 factor_BM_log_Fri_mean_size,
                                 factor_BM_log_Sat_mean_size,
                                 factor_BM_log_Sun_mean_size)
colnames(factor_BM_days_mean_size) = x_grid
rownames(factor_BM_days_mean_size) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun") 

########################  
# comparison with MLFTS
########################

factor_MLFTS_fore_h1_Mon = (10^(factor_MLFTS_Mon_fore$result) - 1)
factor_MLFTS_fore_h1_Tue = (10^(factor_MLFTS_Tue_fore$result) - 1)
factor_MLFTS_fore_h1_Wed = (10^(factor_MLFTS_Wed_fore$result) - 1)
factor_MLFTS_fore_h1_Thu = (10^(factor_MLFTS_Thu_fore$result) - 1)
factor_MLFTS_fore_h1_Fri = (10^(factor_MLFTS_Fri_fore$result) - 1)
factor_MLFTS_fore_h1_Sat = (10^(factor_MLFTS_Sat_fore$result) - 1)
factor_MLFTS_fore_h1_Sun = (10^(factor_MLFTS_Sun_fore$result) - 1)

# MAPE aggregated across 51 sizes

mape_h1_Mon_factor_TS_BM_period = mape_h1_Tue_factor_TS_BM_period = mape_h1_Wed_factor_TS_BM_period = mape_h1_Thu_factor_TS_BM_period = 
mape_h1_Fri_factor_TS_BM_period = mape_h1_Sat_factor_TS_BM_period = mape_h1_Sun_factor_TS_BM_period = array(NA, dim = c(22, n_test, 51),
                                                                                         dimnames = list(3:24, 1:n_test, x_grid))
for(ik in 1:n_test)
{
    for(ij in 3:24)
    {
        for(iw in 1:51)
        {
            mape_h1_Mon_factor_TS_BM_period[(ij-2),ik,iw] = mape(forecast = factor_MLFTS_fore_h1_Mon[iw,ij:24,ik], true = PNSD_Mon_test[ij:24,ik,iw])
            mape_h1_Tue_factor_TS_BM_period[(ij-2),ik,iw] = mape(forecast = factor_MLFTS_fore_h1_Tue[iw,ij:24,ik], true = PNSD_Tue_test[ij:24,ik,iw])
            mape_h1_Wed_factor_TS_BM_period[(ij-2),ik,iw] = mape(forecast = factor_MLFTS_fore_h1_Wed[iw,ij:24,ik], true = PNSD_Wed_test[ij:24,ik,iw])
            mape_h1_Thu_factor_TS_BM_period[(ij-2),ik,iw] = mape(forecast = factor_MLFTS_fore_h1_Thu[iw,ij:24,ik], true = PNSD_Thu_test[ij:24,ik,iw])
            mape_h1_Fri_factor_TS_BM_period[(ij-2),ik,iw] = mape(forecast = factor_MLFTS_fore_h1_Fri[iw,ij:24,ik], true = PNSD_Fri_test[ij:24,ik,iw])
            mape_h1_Sat_factor_TS_BM_period[(ij-2),ik,iw] = mape(forecast = factor_MLFTS_fore_h1_Sat[iw,ij:24,ik], true = PNSD_Sat_test[ij:24,ik,iw])
            mape_h1_Sun_factor_TS_BM_period[(ij-2),ik,iw] = mape(forecast = factor_MLFTS_fore_h1_Sun[iw,ij:24,ik], true = PNSD_Sun_test[ij:24,ik,iw])
            rm(iw)
        }
        rm(ij)
    }
    rm(ik)
}

# take average

factor_TS_BM_period_Mon_mean = apply(mape_h1_Mon_factor_TS_BM_period , c(1, 3), mean)
factor_TS_BM_period_Tue_mean = apply(mape_h1_Tue_factor_TS_BM_period , c(1, 3), mean)
factor_TS_BM_period_Wed_mean = apply(mape_h1_Wed_factor_TS_BM_period , c(1, 3), mean)
factor_TS_BM_period_Thu_mean = apply(mape_h1_Thu_factor_TS_BM_period , c(1, 3), mean)
factor_TS_BM_period_Fri_mean = apply(mape_h1_Fri_factor_TS_BM_period , c(1, 3), mean)
factor_TS_BM_period_Sat_mean = apply(mape_h1_Sat_factor_TS_BM_period , c(1, 3), mean)
factor_TS_BM_period_Sun_mean = apply(mape_h1_Sun_factor_TS_BM_period , c(1, 3), mean)

# by hour

factor_TS_BM_period_Mon_mean_hour = apply(factor_TS_BM_period_Mon_mean, 1, mean)
factor_TS_BM_period_Tue_mean_hour = apply(factor_TS_BM_period_Tue_mean, 1, mean)
factor_TS_BM_period_Wed_mean_hour = apply(factor_TS_BM_period_Wed_mean, 1, mean)
factor_TS_BM_period_Thu_mean_hour = apply(factor_TS_BM_period_Thu_mean, 1, mean)
factor_TS_BM_period_Fri_mean_hour = apply(factor_TS_BM_period_Fri_mean, 1, mean)
factor_TS_BM_period_Sat_mean_hour = apply(factor_TS_BM_period_Sat_mean, 1, mean)
factor_TS_BM_period_Sun_mean_hour = apply(factor_TS_BM_period_Sun_mean, 1, mean)

factor_TS_BM_period_days_mean_hour = rbind(factor_TS_BM_period_Mon_mean_hour,
                                           factor_TS_BM_period_Tue_mean_hour,
                                           factor_TS_BM_period_Wed_mean_hour,
                                           factor_TS_BM_period_Thu_mean_hour,
                                           factor_TS_BM_period_Fri_mean_hour,
                                           factor_TS_BM_period_Sat_mean_hour,
                                           factor_TS_BM_period_Sun_mean_hour)
colnames(factor_TS_BM_period_days_mean_hour) = 3:24
rownames(factor_TS_BM_period_days_mean_hour) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

# by size

factor_TS_BM_period_Mon_mean_size = apply(factor_TS_BM_period_Mon_mean, 2, mean)
factor_TS_BM_period_Tue_mean_size = apply(factor_TS_BM_period_Tue_mean, 2, mean)
factor_TS_BM_period_Wed_mean_size = apply(factor_TS_BM_period_Wed_mean, 2, mean)
factor_TS_BM_period_Thu_mean_size = apply(factor_TS_BM_period_Thu_mean, 2, mean)
factor_TS_BM_period_Fri_mean_size = apply(factor_TS_BM_period_Fri_mean, 2, mean)
factor_TS_BM_period_Sat_mean_size = apply(factor_TS_BM_period_Sat_mean, 2, mean)
factor_TS_BM_period_Sun_mean_size = apply(factor_TS_BM_period_Sun_mean, 2, mean)

factor_TS_BM_period_days_mean_size = rbind(factor_TS_BM_period_Mon_mean_size,
                                           factor_TS_BM_period_Tue_mean_size,
                                           factor_TS_BM_period_Wed_mean_size,
                                           factor_TS_BM_period_Thu_mean_size,
                                           factor_TS_BM_period_Fri_mean_size,
                                           factor_TS_BM_period_Sat_mean_size,
                                           factor_TS_BM_period_Sun_mean_size)
colnames(factor_TS_BM_period_days_mean_size) = x_grid
rownames(factor_TS_BM_period_days_mean_size) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")


length(which(factor_BM_log_Mon_mean < factor_TS_BM_period_Mon_mean)) # 1112
length(which(factor_BM_log_Tue_mean < factor_TS_BM_period_Tue_mean)) # 1116
length(which(factor_BM_log_Wed_mean < factor_TS_BM_period_Wed_mean)) # 1080
length(which(factor_BM_log_Thu_mean < factor_TS_BM_period_Thu_mean)) # 1119
length(which(factor_BM_log_Fri_mean < factor_TS_BM_period_Fri_mean)) # 1088
length(which(factor_BM_log_Sat_mean < factor_TS_BM_period_Sat_mean)) # 1040
length(which(factor_BM_log_Sun_mean < factor_TS_BM_period_Sun_mean)) # 1084

