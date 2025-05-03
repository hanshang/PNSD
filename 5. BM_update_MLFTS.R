####################################
# Dynamic updating via block moving
####################################

# data: data array of dimension 24 x T x 51
# newdata: data matrix new_period x 51
# holdoutdata: data matrix updating_period x 51
# fmethod: forecasting method

BM_update <- function(data, newdata, holdoutdata, fmethod)
{
    n_size = dim(data)[3]
    update_period = (nrow(newdata) + 1):(dim(data)[1])
    new_data_array = array(NA, dim = dim(data), dimnames = dimnames(data))
    for(ik in 1:n_size)
    {
        new_data_array[,,ik] = matrix(c(as.numeric(data[,,ik])[-c(1:nrow(newdata))], newdata[,ik]), 24, dim(data)[2])
        rm(ik)
    }
    
    if(fmethod == "MLFTS")
    {
        BM_fore = matrix(matrix(unlist(MLFTS_model(data_input = aperm(new_data_array, c(3, 2, 1)),
            aux_var = NULL, ncomp_method = "provide", fh = 1, fore_method = "ets")),,24)[,1:length(update_period)], n_size, length(update_period))
    }
    else if(fmethod == "factor_MLFTS")
    {
        BM_fore = matrix(factor_model_fun(data = aperm(new_data_array, c(2, 1, 3)))$fore_val[,1:length(update_period)], n_size, length(update_period))
    }
    else
    {
        warning("fmethod should be chosen from the list.")
    }
    
    # remember to transform back to the original scale
     
    BM_err = vector("numeric", n_size)
    for(ik in 1:n_size)
    {
        BM_err[ik] = ftsa:::mape(forecast = (10^(BM_fore[ik,]) - 1), true = (10^(holdoutdata[,ik]) - 1))
        rm(ik)
    }
    return(BM_err)
}

######################
# evaluation function
######################

# data_set: data_array
# newdata_period: period where newly arrived data become available
# fore_method: forecasting method

BM_update_eval <- function(data_set, newdata_period, fore_method)
{
    MLFTS_BM_h1_err = matrix(NA, dim(data_set)[3], n_test)
    for(iw in 1:n_test)
    {
        n_train = (408 * 3/4 + iw - 1)
  
        MLFTS_BM_h1_err[,iw] = BM_update(data = data_set[,1:n_train,], 
                        newdata = data_set[newdata_period,(n_train + 1),],
                        holdoutdata = matrix(data_set[(1:24)[-newdata_period],(n_train + 1),], length((1:24)[-newdata_period]), dim(data_set)[3]),
                        fmethod = fore_method)  
        print(iw); rm(iw); rm(n_train)
    }
    return(MLFTS_BM_h1_err)
}
  
# evaluation 

BM_log_Mon_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Mon_408, newdata_period = 1:iwk, fore_method = "MLFTS")
BM_log_Tue_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Tue_408, newdata_period = 1:iwk, fore_method = "MLFTS")
BM_log_Wed_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Wed_408, newdata_period = 1:iwk, fore_method = "MLFTS")
BM_log_Thu_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Thu_408, newdata_period = 1:iwk, fore_method = "MLFTS")
BM_log_Fri_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Fri_408, newdata_period = 1:iwk, fore_method = "MLFTS")
BM_log_Sat_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Sat_408, newdata_period = 1:iwk, fore_method = "MLFTS")
BM_log_Sun_list <- foreach(iwk = 2:23) %dopar% BM_update_eval(data_set = PNSD_log_Sun_408, newdata_period = 1:iwk, fore_method = "MLFTS")

BM_log_Mon = BM_log_Tue = BM_log_Wed = BM_log_Thu = BM_log_Fri = 
BM_log_Sat = BM_log_Sun = array(NA, dim = c(22, 51, n_test), dimnames = list(1:22, x_grid, 1:n_test))
for(iwk in 1:22)
{
    BM_log_Mon[iwk,,] = BM_log_Mon_list[[iwk]]
    BM_log_Tue[iwk,,] = BM_log_Tue_list[[iwk]]
    BM_log_Wed[iwk,,] = BM_log_Wed_list[[iwk]]
    BM_log_Thu[iwk,,] = BM_log_Thu_list[[iwk]]
    BM_log_Fri[iwk,,] = BM_log_Fri_list[[iwk]]
    BM_log_Sat[iwk,,] = BM_log_Sat_list[[iwk]]
    BM_log_Sun[iwk,,] = BM_log_Sun_list[[iwk]]
    rm(iwk)
}

# take averages

BM_log_Mon_mean = apply(BM_log_Mon, c(1, 2), mean)
BM_log_Tue_mean = apply(BM_log_Tue, c(1, 2), mean)
BM_log_Wed_mean = apply(BM_log_Wed, c(1, 2), mean)
BM_log_Thu_mean = apply(BM_log_Thu, c(1, 2), mean)
BM_log_Fri_mean = apply(BM_log_Fri, c(1, 2), mean)
BM_log_Sat_mean = apply(BM_log_Sat, c(1, 2), mean)
BM_log_Sun_mean = apply(BM_log_Sun, c(1, 2), mean)

# by hour

BM_log_Mon_mean_hour = apply(BM_log_Mon_mean, 1, mean)
BM_log_Tue_mean_hour = apply(BM_log_Tue_mean, 1, mean)
BM_log_Wed_mean_hour = apply(BM_log_Wed_mean, 1, mean)
BM_log_Thu_mean_hour = apply(BM_log_Thu_mean, 1, mean)
BM_log_Fri_mean_hour = apply(BM_log_Fri_mean, 1, mean)
BM_log_Sat_mean_hour = apply(BM_log_Sat_mean, 1, mean)
BM_log_Sun_mean_hour = apply(BM_log_Sun_mean, 1, mean)

BM_days_mean_hour = rbind(BM_log_Mon_mean_hour,
                          BM_log_Tue_mean_hour,
                          BM_log_Wed_mean_hour,
                          BM_log_Thu_mean_hour,
                          BM_log_Fri_mean_hour,
                          BM_log_Sat_mean_hour,
                          BM_log_Sun_mean_hour)
colnames(BM_days_mean_hour) = 3:24
rownames(BM_days_mean_hour) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun") 
  
# by size

BM_log_Mon_mean_size = apply(BM_log_Mon_mean, 2, mean)
BM_log_Tue_mean_size = apply(BM_log_Tue_mean, 2, mean)
BM_log_Wed_mean_size = apply(BM_log_Wed_mean, 2, mean)
BM_log_Thu_mean_size = apply(BM_log_Thu_mean, 2, mean)
BM_log_Fri_mean_size = apply(BM_log_Fri_mean, 2, mean)
BM_log_Sat_mean_size = apply(BM_log_Sat_mean, 2, mean)
BM_log_Sun_mean_size = apply(BM_log_Sun_mean, 2, mean)

BM_days_mean_size = rbind(BM_log_Mon_mean_size,
                          BM_log_Tue_mean_size,
                          BM_log_Wed_mean_size,
                          BM_log_Thu_mean_size,
                          BM_log_Fri_mean_size,
                          BM_log_Sat_mean_size,
                          BM_log_Sun_mean_size)
colnames(BM_days_mean_size) = x_grid
rownames(BM_days_mean_size) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun") 

# ridge vs BM

length(which(ridge_log_Mon_array_mean < BM_log_Mon_mean)) # 1114
length(which(ridge_log_Tue_array_mean < BM_log_Tue_mean)) # 1118
length(which(ridge_log_Wed_array_mean < BM_log_Wed_mean)) # 1110
length(which(ridge_log_Thu_array_mean < BM_log_Thu_mean)) # 1122
length(which(ridge_log_Fri_array_mean < BM_log_Fri_mean)) # 1062
length(which(ridge_log_Sat_array_mean < BM_log_Sat_mean)) # 928
length(which(ridge_log_Sun_array_mean < BM_log_Sun_mean)) # 1106

# PLS vs BM

length(which(PLS_log_Mon_array_mean < BM_log_Mon_mean)) # 971
length(which(PLS_log_Tue_array_mean < BM_log_Tue_mean)) # 996
length(which(PLS_log_Wed_array_mean < BM_log_Wed_mean)) # 1028
length(which(PLS_log_Thu_array_mean < BM_log_Thu_mean)) # 1042
length(which(PLS_log_Fri_array_mean < BM_log_Fri_mean)) # 1034
length(which(PLS_log_Sat_array_mean < BM_log_Sat_mean)) # 777
length(which(PLS_log_Sun_array_mean < BM_log_Sun_mean)) # 972

# mean ridge vs mean BM

all(rowMeans(ridge_log_Mon_array_mean) < rowMeans(BM_log_Mon_mean)) # TRUE
all(rowMeans(ridge_log_Tue_array_mean) < rowMeans(BM_log_Tue_mean)) # TRUE
all(rowMeans(ridge_log_Wed_array_mean) < rowMeans(BM_log_Wed_mean)) # TRUE
all(rowMeans(ridge_log_Thu_array_mean) < rowMeans(BM_log_Thu_mean)) # TRUE
all(rowMeans(ridge_log_Fri_array_mean) < rowMeans(BM_log_Fri_mean)) # TRUE
all(rowMeans(ridge_log_Sat_array_mean) < rowMeans(BM_log_Sat_mean)) # TRUE
all(rowMeans(ridge_log_Sun_array_mean) < rowMeans(BM_log_Sun_mean)) # TRUE

# mean PLS vs mean BM

all(rowMeans(PLS_log_Mon_array_mean) < rowMeans(BM_log_Mon_mean)) # TRUE
all(rowMeans(PLS_log_Tue_array_mean) < rowMeans(BM_log_Tue_mean)) # TRUE
all(rowMeans(PLS_log_Wed_array_mean) < rowMeans(BM_log_Wed_mean)) # TRUE
all(rowMeans(PLS_log_Thu_array_mean) < rowMeans(BM_log_Thu_mean)) # TRUE
all(rowMeans(PLS_log_Fri_array_mean) < rowMeans(BM_log_Fri_mean)) # TRUE
all(rowMeans(PLS_log_Sat_array_mean) < rowMeans(BM_log_Sat_mean)) # FALSE
all(rowMeans(PLS_log_Sun_array_mean) < rowMeans(BM_log_Sun_mean)) # TRUE

# mean ridge vs mean PLS

all(rowMeans(ridge_log_Mon_array_mean) < rowMeans(PLS_log_Mon_array_mean)) # FALSE
all(rowMeans(ridge_log_Tue_array_mean) < rowMeans(PLS_log_Tue_array_mean)) # FALSE
all(rowMeans(ridge_log_Wed_array_mean) < rowMeans(PLS_log_Wed_array_mean)) # FALSE
all(rowMeans(ridge_log_Thu_array_mean) < rowMeans(PLS_log_Thu_array_mean)) # FALSE
all(rowMeans(ridge_log_Fri_array_mean) < rowMeans(PLS_log_Fri_array_mean)) # FALSE
all(rowMeans(ridge_log_Sat_array_mean) < rowMeans(PLS_log_Sat_array_mean)) # FALSE
all(rowMeans(ridge_log_Sun_array_mean) < rowMeans(PLS_log_Sun_array_mean)) # TRUE

# mean values of ridge and PLS

mean(ridge_log_Mon_array_mean); mean(PLS_log_Mon_array_mean)
mean(ridge_log_Tue_array_mean); mean(PLS_log_Tue_array_mean)
mean(ridge_log_Wed_array_mean); mean(PLS_log_Wed_array_mean)
mean(ridge_log_Thu_array_mean); mean(PLS_log_Thu_array_mean)
mean(ridge_log_Fri_array_mean); mean(PLS_log_Fri_array_mean)
mean(ridge_log_Sat_array_mean); mean(PLS_log_Sat_array_mean)
mean(ridge_log_Sun_array_mean); mean(PLS_log_Sun_array_mean)


########################  
# comparison with MLFTS
########################

# MAPE aggregated across 51 sizes

mape_h1_Mon_TS_BM_period = mape_h1_Tue_TS_BM_period = mape_h1_Wed_TS_BM_period = mape_h1_Thu_TS_BM_period = 
mape_h1_Fri_TS_BM_period = mape_h1_Sat_TS_BM_period = mape_h1_Sun_TS_BM_period = array(NA, dim = c(22, n_test, 51),
                                                        dimnames = list(3:24, 1:n_test, x_grid))
for(ik in 1:n_test)
{
    for(ij in 3:24)
    {
        for(iw in 1:51)
        {
            mape_h1_Mon_TS_BM_period[(ij-2),ik,iw] = mape(forecast = MLFTS_fore_h1_Mon[iw,ij:24,ik], true = PNSD_Mon_test[ij:24,ik,iw])
            mape_h1_Tue_TS_BM_period[(ij-2),ik,iw] = mape(forecast = MLFTS_fore_h1_Tue[iw,ij:24,ik], true = PNSD_Tue_test[ij:24,ik,iw])
            mape_h1_Wed_TS_BM_period[(ij-2),ik,iw] = mape(forecast = MLFTS_fore_h1_Wed[iw,ij:24,ik], true = PNSD_Wed_test[ij:24,ik,iw])
            mape_h1_Thu_TS_BM_period[(ij-2),ik,iw] = mape(forecast = MLFTS_fore_h1_Thu[iw,ij:24,ik], true = PNSD_Thu_test[ij:24,ik,iw])
            mape_h1_Fri_TS_BM_period[(ij-2),ik,iw] = mape(forecast = MLFTS_fore_h1_Fri[iw,ij:24,ik], true = PNSD_Fri_test[ij:24,ik,iw])
            mape_h1_Sat_TS_BM_period[(ij-2),ik,iw] = mape(forecast = MLFTS_fore_h1_Sat[iw,ij:24,ik], true = PNSD_Sat_test[ij:24,ik,iw])
            mape_h1_Sun_TS_BM_period[(ij-2),ik,iw] = mape(forecast = MLFTS_fore_h1_Sun[iw,ij:24,ik], true = PNSD_Sun_test[ij:24,ik,iw])
            rm(iw)
        }
        rm(ij)
    }
    rm(ik)
}

# take average

TS_BM_period_Mon_mean = apply(mape_h1_Mon_TS_BM_period , c(1, 3), mean)
TS_BM_period_Tue_mean = apply(mape_h1_Tue_TS_BM_period , c(1, 3), mean)
TS_BM_period_Wed_mean = apply(mape_h1_Wed_TS_BM_period , c(1, 3), mean)
TS_BM_period_Thu_mean = apply(mape_h1_Thu_TS_BM_period , c(1, 3), mean)
TS_BM_period_Fri_mean = apply(mape_h1_Fri_TS_BM_period , c(1, 3), mean)
TS_BM_period_Sat_mean = apply(mape_h1_Sat_TS_BM_period , c(1, 3), mean)
TS_BM_period_Sun_mean = apply(mape_h1_Sun_TS_BM_period , c(1, 3), mean)

# by hour

TS_BM_period_Mon_mean_hour = apply(TS_BM_period_Mon_mean, 1, mean)
TS_BM_period_Tue_mean_hour = apply(TS_BM_period_Tue_mean, 1, mean)
TS_BM_period_Wed_mean_hour = apply(TS_BM_period_Wed_mean, 1, mean)
TS_BM_period_Thu_mean_hour = apply(TS_BM_period_Thu_mean, 1, mean)
TS_BM_period_Fri_mean_hour = apply(TS_BM_period_Fri_mean, 1, mean)
TS_BM_period_Sat_mean_hour = apply(TS_BM_period_Sat_mean, 1, mean)
TS_BM_period_Sun_mean_hour = apply(TS_BM_period_Sun_mean, 1, mean)

TS_BM_period_days_mean_hour = rbind(TS_BM_period_Mon_mean_hour,
                                    TS_BM_period_Tue_mean_hour,
                                    TS_BM_period_Wed_mean_hour,
                                    TS_BM_period_Thu_mean_hour,
                                    TS_BM_period_Fri_mean_hour,
                                    TS_BM_period_Sat_mean_hour,
                                    TS_BM_period_Sun_mean_hour)
colnames(TS_BM_period_days_mean_hour) = 3:24
rownames(TS_BM_period_days_mean_hour) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

# by size

TS_BM_period_Mon_mean_size = apply(TS_BM_period_Mon_mean, 2, mean)
TS_BM_period_Tue_mean_size = apply(TS_BM_period_Tue_mean, 2, mean)
TS_BM_period_Wed_mean_size = apply(TS_BM_period_Wed_mean, 2, mean)
TS_BM_period_Thu_mean_size = apply(TS_BM_period_Thu_mean, 2, mean)
TS_BM_period_Fri_mean_size = apply(TS_BM_period_Fri_mean, 2, mean)
TS_BM_period_Sat_mean_size = apply(TS_BM_period_Sat_mean, 2, mean)
TS_BM_period_Sun_mean_size = apply(TS_BM_period_Sun_mean, 2, mean)

TS_BM_period_days_mean_size = rbind(TS_BM_period_Mon_mean_size,
                                    TS_BM_period_Tue_mean_size,
                                    TS_BM_period_Wed_mean_size,
                                    TS_BM_period_Thu_mean_size,
                                    TS_BM_period_Fri_mean_size,
                                    TS_BM_period_Sat_mean_size,
                                    TS_BM_period_Sun_mean_size)
colnames(TS_BM_period_days_mean_size) = x_grid
rownames(TS_BM_period_days_mean_size) = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")


length(which(BM_log_Mon_mean < TS_BM_period_Mon_mean)) # 1121
length(which(BM_log_Tue_mean < TS_BM_period_Tue_mean)) # 1122
length(which(BM_log_Wed_mean < TS_BM_period_Wed_mean)) # 1092
length(which(BM_log_Thu_mean < TS_BM_period_Thu_mean)) # 1122
length(which(BM_log_Fri_mean < TS_BM_period_Fri_mean)) # 1122
length(which(BM_log_Sat_mean < TS_BM_period_Sat_mean)) # 1041
length(which(BM_log_Sun_mean < TS_BM_period_Sun_mean)) # 1092

# filled contour

filled.contour(3:24, x_grid, BM_log_Mon_mean)
filled.contour(3:24, x_grid, TS_BM_period_Mon_mean)

