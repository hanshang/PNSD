# load R packages

source("load_packages.R")

# load data set

load("~/Dropbox/Todos/Hanlin-Isra/Data/PNSD.RData")

# data and time

date_time = (PNSD$Data0)[,1]
date = unique(substr(date_time, 1, 10))
day_of_week <- weekdays(as.Date(date))
n_date = length(date)

indices_by_day <- split(seq_along(date), day_of_week)
class(indices_by_day) # list

# there are 2861 dates from Jan-1-2011 to 31-Oct-2018

x_grid = PNSD$x
PNSD_array = array(NA, dim = c(24, n_date, 51), dimnames = list(1:24, 1:n_date, x_grid))
for(ik in 1:51)
{
    PNSD_array[,,ik] = matrix(PNSD$Data[,ik], 24, n_date)
    rm(ik)
}

dim(PNSD_array)   # 24 2861 51
range(PNSD_array) # 0  232629.4

# replace 0 by NA

PNSD_array_zero_NA = array(NA, dim(PNSD_array), dimnames = dimnames(PNSD_array))
for(ik in 1:n_date)
{
    for(ij in 1:51)
    {
        data_dum = PNSD_array[,ik,ij]
        index = which(data_dum == 0)
        if(length(index) > 0)
        {
             PNSD_array_zero_NA[,ik,ij] = replace(data_dum, index, NA)
        }
        else
        {
            PNSD_array_zero_NA[,ik,ij] = data_dum
        }
        rm(data_dum); rm(ij)
    }
    rm(ik)
}

dates_NA = which(apply(PNSD_array_zero_NA, 2, function(x) any(is.na(x))))

# na.interp interpolation in the forecast package

PNSD_array_na_interp = array(NA, dim = dim(PNSD_array), dimnames = dimnames(PNSD_array))
for(ik in 1:n_date)
{
    for(ij in 1:51)
    {
        PNSD_array_na_interp[,ik,ij] = na.interp(PNSD_array_zero_NA[,ik,ij])
        rm(ij)
    }
    rm(ik)
}

any(is.na(PNSD_array_na_interp)) # FALSE
round(range(PNSD_array_na_interp), 4) # 0.0127 232629.4400

# rainbow plot

par(mfrow = c(1, 2))
plot(fts(x_grid, t(PNSD_array_na_interp[,1,])), xlab = "Size", ylab = "Count")
plot(fts(x_grid, t(PNSD_array_na_interp[,2,])), xlab = "Size", ylab = "Count")

# take transformation, square-root, cubic-root or log_10(count + 1)

PNSD_square = sqrt(PNSD_array_na_interp)
PNSD_cubic  = (PNSD_array_na_interp)^(1/3)
PNSD_log    = log(PNSD_array_na_interp + 1, base = 10)

# separate days

PNSD_log_Monday    = PNSD_log[,indices_by_day$Monday,]
PNSD_log_Tuesday   = PNSD_log[,indices_by_day$Tuesday,]
PNSD_log_Wednesday = PNSD_log[,indices_by_day$Wednesday,]
PNSD_log_Thursday  = PNSD_log[,indices_by_day$Thursday,]
PNSD_log_Friday    = PNSD_log[,indices_by_day$Friday,]
PNSD_log_Saturday  = PNSD_log[,indices_by_day$Saturday,]
PNSD_log_Sunday    = PNSD_log[,indices_by_day$Sunday,]

dim(PNSD_log_Monday)[2]    # 409
dim(PNSD_log_Tuesday)[2]   # 409
dim(PNSD_log_Wednesday)[2] # 409
dim(PNSD_log_Thursday)[2]  # 408
dim(PNSD_log_Friday)[2]    # 408
dim(PNSD_log_Saturday)[2]  # 409
dim(PNSD_log_Sunday)[2]    # 409

# because of different days, we consider T = 408 for consistency (Thursday - Wednesday(last day))

PNSD_log_Mon_408 = PNSD_log_Monday[,1:408,]
PNSD_log_Tue_408 = PNSD_log_Tuesday[,1:408,]
PNSD_log_Wed_408 = PNSD_log_Wednesday[,1:408,]
PNSD_log_Thu_408 = PNSD_log_Thursday[,1:408,]
PNSD_log_Fri_408 = PNSD_log_Friday[,1:408,]
PNSD_log_Sat_408 = PNSD_log_Saturday[,1:408,]
PNSD_log_Sun_408 = PNSD_log_Sunday[,1:408,]

# split training (75%) and testing (25%) sets

PNSD_Mon_test = 10^(PNSD_log_Mon_408[,(1 + (408 * 3/4)):408,]) - 1
PNSD_Tue_test = 10^(PNSD_log_Tue_408[,(1 + (408 * 3/4)):408,]) - 1
PNSD_Wed_test = 10^(PNSD_log_Wed_408[,(1 + (408 * 3/4)):408,]) - 1
PNSD_Thu_test = 10^(PNSD_log_Thu_408[,(1 + (408 * 3/4)):408,]) - 1
PNSD_Fri_test = 10^(PNSD_log_Fri_408[,(1 + (408 * 3/4)):408,]) - 1
PNSD_Sat_test = 10^(PNSD_log_Sat_408[,(1 + (408 * 3/4)):408,]) - 1
PNSD_Sun_test = 10^(PNSD_log_Sun_408[,(1 + (408 * 3/4)):408,]) - 1

