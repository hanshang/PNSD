###########
# rankings
###########

n_methods = 2
MAPE_TS_size_array = array(NA, dim = c(7, 51, n_methods), 
                      dimnames = list(1:7, 1:51, c("MLFTS", "Factor+MLFTS")))
MAPE_TS_size_array[,,1] = mape_h1_days_TS_size_mean
MAPE_TS_size_array[,,2] = mape_h1_days_factor_size_mean

MAPE_ranking_size_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(MAPE_ranking_size_mat) = paste("M", 1:n_methods, sep = "")
rownames(MAPE_ranking_size_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(MAPE_TS_size_array[ih,,],1, which.min))
    if(length(temp_rank) == n_methods)
    {
      MAPE_ranking_size_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, n_methods)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        MAPE_ranking_size_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

###################
# create a heatmap
###################

MAPE_ranking_size_mat = data.frame(Horizon = 1:7, MAPE_ranking_size_mat)

MAPE_days_size = MAPE_ranking_size_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

# save figure

ggsave("Fig_1b.png")
ggplot(MAPE_days_size, aes(Model, ordered(Horizon, levels = 7:1))) +
  geom_tile(aes(fill = count)) +
  geom_text(aes(label = count)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none") + 
  scale_fill_gradient2(high="darkgray",mid="lightgray",low="white", 
                       na.value="yellow", midpoint = 12) + 
  ylab("Day of the week") +
  scale_x_discrete(labels = c("MLFTS", "Factor+MLFTS")) + 
  scale_y_discrete(labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")) + 
  theme_bw()
dev.off()

