###########
# rankings
###########

n_methods = 4
CPD_95_TS_array = array(NA, dim = c(7, 51, n_methods), 
                     dimnames = list(1:7, x_grid, c("MLFTS\nsd", "MLFTS\nconformal", "Factor+MLFTS\nsd", "Factor+MLFTS\nconformal")))
CPD_95_TS_array[,,1] = t(PNSD_CPD_sd_test_0.95)
CPD_95_TS_array[,,2] = t(PNSD_CPD_conformal_test_0.95)
CPD_95_TS_array[,,3] = t(factor_PNSD_CPD_sd_test_0.95)
CPD_95_TS_array[,,4] = t(factor_PNSD_CPD_conformal_test_0.95)

CPD_95_ranking_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(CPD_95_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(CPD_95_ranking_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(CPD_95_TS_array[ih,,],1, which.min))
    if(length(temp_rank) == 4)
    {
        CPD_95_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, 4)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        CPD_95_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

###################
# create a heatmap
###################

CPD_95_ranking_mat = data.frame(Horizon = 1:7, CPD_95_ranking_mat)

CPD_95_days = CPD_95_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

ggsave("Fig_5b.png")
ggplot(CPD_95_days, aes(Model, ordered(Horizon, levels = 7:1))) +
  geom_tile(aes(fill = count)) +
  geom_text(aes(label = count)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none") + 
  scale_fill_gradient2(high="darkgray",mid="lightgray",low="white", 
                       na.value="yellow", midpoint = ceiling(dim(CPD_95_TS_array)[2]/n_methods)) + 
  ylab("Day of the week") +
  scale_x_discrete(labels = c("sd\nMLFTS", "conformal\nMLFTS", "sd\nFactor+MLFTS", "conformal\nFactor+MLFTS")) + 
  scale_y_discrete(labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")) + 
  ggtitle("Nominal coverage probability of 95%") + 
  theme_bw()
dev.off()

