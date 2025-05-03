###########
# rankings
###########

n_methods = 7
MAPE_dynamic_hour_array = array(NA, dim = c(7, 22, n_methods), 
                dimnames = list(1:7, 3:24, c("TS\nMLFTS", "BM\nMLFTS", "TS\nFactor + MLFTS", "BM\nFactor + MLFTS", "OLS", "Ridge", "PLS")))
MAPE_dynamic_hour_array[,,1] = TS_BM_period_days_mean_hour
MAPE_dynamic_hour_array[,,2] = BM_days_mean_hour
MAPE_dynamic_hour_array[,,3] = factor_TS_BM_period_days_mean_hour
MAPE_dynamic_hour_array[,,4] = factor_BM_days_mean_hour
MAPE_dynamic_hour_array[,,5] = ols_days_array_mean_hour
MAPE_dynamic_hour_array[,,6] = ridge_days_array_mean_hour
MAPE_dynamic_hour_array[,,7] = PLS_days_array_mean_hour

MAPE_dynamic_hour_ranking_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(MAPE_dynamic_hour_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(MAPE_dynamic_hour_ranking_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(MAPE_dynamic_hour_array[ih,,],1, which.min))
    if(length(temp_rank) == n_methods)
    {
        MAPE_dynamic_hour_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, n_methods)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        MAPE_dynamic_hour_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

###################
# create a heatmap
###################

MAPE_dynamic_hour_ranking_mat = data.frame(Horizon = 1:7, MAPE_dynamic_hour_ranking_mat)

MAPE_dynamic = MAPE_dynamic_hour_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

ggsave("Fig_2a.png")
ggplot(MAPE_dynamic, aes(Model, ordered(Horizon, levels = 7:1))) +
  geom_tile(aes(fill = count)) +
  geom_text(aes(label = count)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none") + 
  scale_fill_gradient2(high="darkgray",mid="lightgray",low="white", midpoint = ceiling(dim(MAPE_dynamic_hour_array)[2]/n_methods)) + 
  ylab("Updating period") +
  scale_x_discrete(labels = c("TS\nMLFTS", "BM\nMLFTS", "TS\nFactor\nMLFTS", "BM\nFactor\nMLFTS", "OLS", "Ridge", "PLS")) + 
  scale_y_discrete(labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")) + 
  theme_bw()
dev.off()

