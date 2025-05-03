###########
# rankings
###########

n_methods = 4
score_95_TS_array = array(NA, dim = c(7, 51, n_methods), 
                        dimnames = list(1:7, x_grid, c("MLFTS\nsd", "MLFTS\nconformal", "Factor+MLFTS\nsd", "Factor+MLFTS\nconformal")))
score_95_TS_array[,,1] = t(PNSD_score_sd_test_0.95)
score_95_TS_array[,,2] = t(PNSD_score_conformal_test_0.95)
score_95_TS_array[,,3] = t(factor_PNSD_score_sd_test_0.95)
score_95_TS_array[,,4] = t(factor_PNSD_score_conformal_test_0.95)


score_95_ranking_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(score_95_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(score_95_ranking_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(score_95_TS_array[ih,,],1, which.min))
    if(length(temp_rank) == 4)
    {
        score_95_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, 4)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        score_95_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

###################
# create a heatmap
###################

score_95_ranking_mat = data.frame(Horizon = 1:7, score_95_ranking_mat)

score_95_days = score_95_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

ggsave("Fig_5d.png")
ggplot(score_95_days, aes(Model, ordered(Horizon, levels = 7:1))) +
  geom_tile(aes(fill = count)) +
  geom_text(aes(label = count)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none") + 
  scale_fill_gradient2(high="darkgray",mid="lightgray",low="white", 
                       na.value="yellow", midpoint = ceiling(dim(score_95_TS_array)[2]/n_methods)) + 
  ylab("Day of the week") +
  scale_x_discrete(labels = c("sd\nMLFTS", "conformal\nMLFTS", "sd\nFactor+MLFTS", "conformal\nFactor+MLFTS")) + 
  scale_y_discrete(labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")) + 
  theme_bw()
dev.off()

