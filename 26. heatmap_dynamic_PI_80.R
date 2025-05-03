################
### alpha = 0.2
################

## CPD

# rankings

n_methods = 6
CPD_dynamic_array = array(NA, dim = c(7, 22, n_methods), 
    dimnames = list(1:7, 3:24, c("MLFTS", "factor + MLFTS", "Ridge", "MLFTS", "factor + MLFTS", "Ridge")))
CPD_dynamic_array[,,1] = PNSD_days_sd_test_TS_CPD
CPD_dynamic_array[,,2] = factor_PNSD_days_sd_test_TS_CPD
CPD_dynamic_array[,,3] = PNSD_days_sd_ridge_PI_test_array_CPD

CPD_dynamic_array[,,4] = PNSD_days_TS_conformal_dynamic_CPD
CPD_dynamic_array[,,5] = factor_PNSD_days_TS_conformal_dynamic_CPD
CPD_dynamic_array[,,6] = PNSD_days_ridge_conformal_array_CPD


CPD_dynamic_ranking_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(CPD_dynamic_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(CPD_dynamic_ranking_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(CPD_dynamic_array[ih,,],1, which.min))
    if(length(temp_rank) == n_methods)
    {
        CPD_dynamic_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, n_methods)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        CPD_dynamic_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

# create a heatmap

CPD_dynamic_ranking_mat = data.frame(Horizon = 1:7, CPD_dynamic_ranking_mat)

CPD_dynamic = CPD_dynamic_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

ggsave("Fig_6a.png")
ggplot(CPD_dynamic, aes(Model, ordered(Horizon, levels = 7:1))) +
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

# rankings

n_methods = 6
score_dynamic_array = array(NA, dim = c(7, 22, n_methods), 
                          dimnames = list(1:7, 3:24, c("MLFTS", "factor + MLFTS", "Ridge", "MLFTS", "factor + MLFTS", "Ridge")))
score_dynamic_array[,,1] = PNSD_days_sd_test_TS_score
score_dynamic_array[,,2] = factor_PNSD_days_sd_test_TS_score
score_dynamic_array[,,3] = PNSD_days_sd_ridge_PI_test_array_score

score_dynamic_array[,,4] = PNSD_days_TS_conformal_dynamic_score
score_dynamic_array[,,5] = factor_PNSD_days_TS_conformal_dynamic_score
score_dynamic_array[,,6] = PNSD_days_ridge_conformal_array_score


score_dynamic_ranking_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(score_dynamic_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(score_dynamic_ranking_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(score_dynamic_array[ih,,],1, which.min))
    if(length(temp_rank) == n_methods)
    {
        score_dynamic_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, n_methods)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        score_dynamic_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

# create a heatmap

score_dynamic_ranking_mat = data.frame(Horizon = 1:7, score_dynamic_ranking_mat)

score_dynamic = score_dynamic_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

ggsave("Fig_6b.png")
ggplot(score_dynamic, aes(Model, ordered(Horizon, levels = 7:1))) +
  geom_tile(aes(fill = count)) +
  geom_text(aes(label = count)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none") + 
  scale_fill_gradient2(high="darkgray",mid="lightgray",low="white", midpoint = ceiling(dim(CPD_dynamic_array)[2]/n_methods)) + 
  ylab("Updating period") +
  scale_x_discrete(labels = c("MLFTS\nsd", "factor + MLFTS\nsd", "Ridge\nsd", "MLFTS\nconformal", "factor + MLFTS\nconformal", "Ridge\nconformal")) + 
  scale_y_discrete(labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")) + 
  theme_bw()
dev.off()

