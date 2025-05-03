################
### alpha = 0.05
################

## CPD

# rankings

n_methods = 6
CPD_0.05_dynamic_array = array(NA, dim = c(7, 22, n_methods), 
                          dimnames = list(1:7, 3:24, c("MLFTS", "factor + MLFTS", "Ridge", "MLFTS", "factor + MLFTS", "Ridge")))
CPD_0.05_dynamic_array[,,1] = PNSD_days_sd_test_TS_CPD_0.05
CPD_0.05_dynamic_array[,,2] = factor_PNSD_days_sd_test_TS_0.05_CPD
CPD_0.05_dynamic_array[,,3] = PNSD_days_sd_ridge_PI_test_array_CPD_0.05

CPD_0.05_dynamic_array[,,4] = PNSD_days_TS_conformal_dynamic_0.05_CPD
CPD_0.05_dynamic_array[,,5] = factor_PNSD_days_TS_conformal_dynamic_0.05_CPD
CPD_0.05_dynamic_array[,,6] = PNSD_days_ridge_conformal_array_CPD_0.05

CPD_0.05_dynamic_ranking_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(CPD_0.05_dynamic_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(CPD_0.05_dynamic_ranking_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(CPD_0.05_dynamic_array[ih,,],1, which.min))
    if(length(temp_rank) == n_methods)
    {
        CPD_0.05_dynamic_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, n_methods)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        CPD_0.05_dynamic_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

# create a heatmap

CPD_0.05_dynamic_ranking_mat = data.frame(Horizon = 1:7, CPD_0.05_dynamic_ranking_mat)

CPD_0.05_dynamic = CPD_0.05_dynamic_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

ggsave("Fig_6c.png")
ggplot(CPD_0.05_dynamic, aes(Model, ordered(Horizon, levels = 7:1))) +
  geom_tile(aes(fill = count)) +
  geom_text(aes(label = count)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none") + 
  scale_fill_gradient2(high="darkgray",mid="lightgray",low="white", midpoint = 5) + 
  ylab("Updating period") +
  scale_x_discrete(labels = c("MLFTS\nsd", "factor + MLFTS\nsd", "Ridge\nsd", "MLFTS\nconformal", "factor + MLFTS\nconformal", "Ridge\nconformal")) + 
  scale_y_discrete(labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")) + 
  theme_bw()
dev.off()

## interval score

n_methods = 6
score_0.05_dynamic_array = array(NA, dim = c(7, 22, n_methods), 
                               dimnames = list(1:7, 3:24, c("MLFTS", "factor + MLFTS", "Ridge", "MLFTS", "factor + MLFTS", "Ridge")))
score_0.05_dynamic_array[,,1] = PNSD_days_sd_test_TS_score_0.05
score_0.05_dynamic_array[,,2] = factor_PNSD_days_sd_test_TS_0.05_score
score_0.05_dynamic_array[,,3] = PNSD_days_sd_ridge_PI_test_array_score_0.05

score_0.05_dynamic_array[,,4] = PNSD_days_TS_conformal_dynamic_0.05_score
score_0.05_dynamic_array[,,5] = factor_PNSD_days_TS_conformal_dynamic_0.05_score
score_0.05_dynamic_array[,,6] = PNSD_days_ridge_conformal_array_score_0.05


score_0.05_dynamic_ranking_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(score_0.05_dynamic_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(score_0.05_dynamic_ranking_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(score_0.05_dynamic_array[ih,,],1, which.min))
    if(length(temp_rank) == n_methods)
    {
        score_0.05_dynamic_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, n_methods)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        score_0.05_dynamic_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

# create a heatmap

score_0.05_dynamic_ranking_mat = data.frame(Horizon = 1:7, score_0.05_dynamic_ranking_mat)

score_0.05_dynamic = score_0.05_dynamic_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

ggsave("Fig_6d.png")
ggplot(score_0.05_dynamic, aes(Model, ordered(Horizon, levels = 7:1))) +
  geom_tile(aes(fill = count)) +
  geom_text(aes(label = count)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none") + 
  scale_fill_gradient2(high="darkgray",mid="lightgray",low="white", midpoint = ceiling(dim(CPD_0.05_dynamic_array)[2]/n_methods)) + 
  ylab("Updating period") +
  scale_x_discrete(labels = c("MLFTS\nsd", "factor + MLFTS\nsd", "Ridge\nsd", "MLFTS\nconformal", "factor + MLFTS\nconformal", "Ridge\nconformal")) + 
  scale_y_discrete(labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")) + 
  theme_bw()
dev.off()

