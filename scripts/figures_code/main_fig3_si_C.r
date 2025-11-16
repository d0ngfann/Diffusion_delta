# 메인 Figure 3 및 SI 섹션 C 그림 생성 (네트워크 유형별 확산 차이 KS 검정)
###############################
# Figure 3
#Figure S3 from SI Section C - KS Test
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(cowplot)
library(stats)

## FIGURE SPECIFIC FUNCTIONS
get_sum_ks <- function(file_list){
  df_long <- data.frame()
  for(file in file_list){
    print(file)
    file2 <- gsub("0.0perc", "1.0perc", file)
    df_temp <- rbind(read_csv(file), read_csv(file2))%>%
      filter(seed_strat == "seed_strat_one_nbrs")%>%
      select(-starts_with(c("spread_time", "time_to" , "n_", "full_")))%>%
      group_by(k, perc_rewire, sig_thresh_sd, 
               time_of_influence, p1, p2, thrshld, seed_strat, G_name, randomization_type)  %>% 
      summarise(value = list(max_spread_norm), mean = mean(max_spread_norm))%>%
      pivot_wider(names_from = "perc_rewire", values_from = c("value", "mean"), names_prefix = "ms") %>% 
      mutate(ks_test = map2(value_ms0, value_ms1, stats::ks.test)) %>% 
      mutate(ks_stat = map_dbl(ks_test, "statistic"),
             ks_pval = map_dbl(ks_test, "p.value"),
             mdif = mean_ms0-mean_ms1,
             mavg = (mean_ms0 + mean_ms1)/2,) %>%
      select(-starts_with("value"), -ks_test)
    
    df_long <- rbind(df_long, df_temp)}
  return(df_long)
}

get_ks_bar_plot<- function(df, group, threshold, title,subtitle, x_axis){
  group_var <- rlang::sym(group)
  
  df_sum_plot_ks <- df %>%
    mutate(mdif = mean_ms0 - mean_ms1,
           sig = ifelse(ks_pval < 0.05, 1, 0),
           mavg =  (mean_ms0 +mean_ms1)/2,
           l_win = ifelse(mdif > 0 & sig == 1, 1,0),
           r_win = ifelse(mdif < 0 & sig == 1, 1,0),
           equal_g = ifelse(l_win == 0 & r_win == 0, 1,0),
           n_eq_over = ifelse(equal_g ==1 & mean_ms0 >= threshold & mean_ms1 >= threshold, 1,0),
           n_eq_under = ifelse(equal_g ==1 & n_eq_over ==0, 1,0) ) %>%
    group_by(!!group_var)%>% 
    summarise(n_l = sum(l_win), n_r = sum(r_win), n_eq_over = sum(n_eq_over), n_eq_under= sum(n_eq_under))%>%
    pivot_longer(starts_with("n_"), names_to = "type", values_to = "count") %>% 
    group_by(!!group_var)%>% 
    mutate(total = sum(count),
           prop = count/total,
           prop_label = ifelse(type == "n_l", as.character(round(prop, digits = 2)), NA),
           type = case_when( type == "n_l" ~ "Greater Spread on Clustered Network",
                             type == "n_r"~ "Greater Spread on Random Network",
                             type == "n_eq_over" ~ "Full Equal Spread (Random Network Faster)",
                             TRUE  ~ "Minimal Equal Spread"),
           type = factor(type, ordered = TRUE, 
                         levels = c("Full Equal Spread (Random Network Faster)", "Minimal Equal Spread", 
                                    "Greater Spread on Clustered Network", 
                                    "Greater Spread on Random Network")))
  
  plot<- ggplot(df_sum_plot_ks, aes(x = !!group_var, y =prop, fill = type , group = type))+geom_col()+
    scale_fill_manual(values = c("gray80","gray30","dodgerblue", "goldenrod1"))+theme_bw(base_size = 18)+
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.00))+
    geom_label(aes(label = prop_label),fill = "dodgerblue",color = "white", position = position_stack(vjust = 0.5))+
    labs(title = title , 
         subtitle = subtitle,
         y = expression(paste( "Proportion of ", p[1] <= p[2], " space")), x = x_axis, fill = "")
         
  return(plot)
}

get_reg_bar_plot<- function(df, group, margin, threshold,title, subtitle, x_axis){
  group_var <- rlang::sym(group)
  neg_margin = -1*margin
  
  df_sum_plot_ks <- df %>%
    mutate(mdif = mean_ms0 - mean_ms1,
           #sig = ifelse(ks_pval < 0.05, 1, 0),
           mavg =  (mean_ms0 +mean_ms1)/2,
           l_win = ifelse(mdif > margin, 1,0),
           r_win = ifelse(mdif < neg_margin, 1,0),
           equal_g = ifelse(l_win == 0 & r_win == 0, 1,0),
           n_eq_over = ifelse(equal_g ==1 & mean_ms0 >= threshold & mean_ms1 >= threshold, 1,0),
           n_eq_under = ifelse(equal_g ==1 & n_eq_over ==0, 1,0)) %>%
    group_by(!!group_var)%>% 
    summarise(n_l = sum(l_win), n_r = sum(r_win), n_eq_over = sum(n_eq_over), n_eq_under= sum(n_eq_under))%>%
    pivot_longer(starts_with("n_"), names_to = "type", values_to = "count") %>% 
    group_by(!!group_var)%>% 
    mutate(total = sum(count),
           prop = count/total,
           prop_label = ifelse(type == "n_l", as.character(round(prop, digits = 2)), NA),
           type = case_when( type == "n_l" ~ "Greater Spread on Clustered Network",
                             type == "n_r"~ "Greater Spread on Random Network",
                             type == "n_eq_over" ~ "Full Equal Spread (Random Network Faster)",
                             TRUE  ~ "Minimal Equal Spread"),
           type = factor(type, ordered = TRUE, 
                         levels = c("Full Equal Spread (Random Network Faster)", "Minimal Equal Spread", 
                                    "Greater Spread on Clustered Network", 
                                    "Greater Spread on Random Network")))
  
  plot<- ggplot(df_sum_plot_ks, aes(x = !!group_var, y =prop, fill = type , group = type))+geom_col()+
    scale_fill_manual(values = c("gray80","gray30","dodgerblue", "goldenrod1"))+theme_bw(base_size = 18)+
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.00))+
    geom_label(aes(label = prop_label),fill = "dodgerblue",color = "white", position = position_stack(vjust = 0.5),
               size = 3)+
    labs(subtitle = paste(title, "\n",subtitle, sep = ""),
         y = expression(paste( "Proportion of ", p[1] <= p[2], " space")), x = x_axis, fill = "")
  
  return(plot)
}


## CREATE AND SAVE AGGREGATE FILES
dir = "output_data/cf"

# By Degree, k
file_list_k <-  list.files(dir, pattern = "WS_.+_1T.+2td_.+0.0perc_0.0sd_ES", full.names = TRUE)

df_sum_k <- get_sum_ks(file_list_k)

write.csv(df_sum_k,"output_data/cf_sum_k_ks_df.csv", row.names = FALSE )

## By Social Reinforcement Threshold, i
file_list_i <- list.files(dir, pattern = "WS_20k_.+_1T_.+0.0perc_0.0sd_ES", full.names = TRUE)
df_sum_i <- get_sum_ks(file_list_i)
write.csv(df_sum_i,"output_data/cf_sum_i_ks_df.csv", row.names = FALSE )

## By Time of Influence, T
file_list_T<- c() 
for(filename in list.files(dir, pattern = "WS_8k_.+_2td_.+0.0perc_0.0sd_ES", full.names = TRUE)){
  if(grepl("_0T_", filename) == FALSE){
    file_list_T <- c(file_list_T, filename)
  }
}

df_sum_T <- get_sum_ks(file_list_T)
write.csv(df_sum_T,"output_data/cf_sum_T_ks_df.csv", row.names = FALSE )

### MAKE FIGURES

df_sum_k <- read.csv("output_data/cf_sum_k_ks_df.csv")

df_sum_i<- read.csv("output_data/cf_sum_i_ks_df.csv")

df_sum_T <- read.csv("output_data/cf_sum_T_ks_df.csv")

#### Figure 3
##### Figure 3A
k_reg_plot_0.5 <- get_reg_bar_plot(df_sum_k, "k", 0.05, 0.6,"By Degree", "Example parameters i = 2, T = 1",
                                   "Degree (k)")+
  scale_x_continuous(breaks = c(4,6,8,10,12,14,16,18,20), expand = c(0.01, 0.01)) #+coord_fixed(9/0.8)
  
##### Figure 3B
T_reg_plot_0.5 <- get_reg_bar_plot(df_sum_T, "time_of_influence", 0.05, 0.6,
                                   "By Time of Influence","Example parameters k = 8, i = 2", 
                                   "Time steps influential (T)")+
  scale_x_continuous(breaks = 1:10, expand = c(0.01, 0.01)) #+ coord_fixed(5/0.8)


##### Figure 3C
i_reg_plot_0.5 <- get_reg_bar_plot(df_sum_i, "thrshld", 0.05, 0.6,
                                   "By Social Reinforcement Threshold", "Example parameters k = 20, T = 1", 
                                   "Social Reinforcement Threshold (i)")+
  scale_x_continuous(breaks = 2:10, expand = c(0.01, 0.01)) #+coord_fixed(4.5/0.8)


# Combine in wide format
# size = 3
legend <- get_legend(T_reg_plot_0.5 +guides(fill=guide_legend(ncol=2, override.aes = list(size = 0.5))))
plots<- plot_grid(k_reg_plot_0.5+theme_bw(base_size = 16) + theme(legend.position="none"), 
          T_reg_plot_0.5+theme_bw(base_size = 16) + theme(legend.position="none"),
          i_reg_plot_0.5+theme_bw(base_size = 16) + theme(legend.position="none"), 
        labels = c("A.", "B.", "C."), label_size = 24,
          ncol = 3)
plot_grid(plots, legend, ncol = 1, rel_heights = c(1, 0.2))
ggsave("output_plots/main_fig3_wide.pdf", width = 14, height = 5, units = "in")

#### Figure S3 - KS test bar plot with 60 percent_threshold
thresh = 0.6
k_ks_plot_0.5 <- get_ks_bar_plot(df_sum_k, "k", thresh,"By Degree", "Example parameters i = 2, T = 1",
                                   "Degree (k)")+
  scale_x_continuous(breaks = c(4,6,8,10,12,14,16,18,20))+
  coord_fixed(9/0.8)

i_ks_plot_0.5 <- get_ks_bar_plot(df_sum_i, "thrshld", thresh,
                                 "By Social Reinforcement Threshold", "Example parameters k = 20, T = 1", 
                                 "Social Reinforcement Threshold (i)")+
  scale_x_continuous(breaks = 2:10)+
  coord_fixed(4.5/0.8)

T_ks_plot_0.5 <- get_ks_bar_plot(df_sum_T, "time_of_influence", thresh,
                                   "By Time of Influence","Example parameters k = 8, i = 2", 
                                 "Time steps influential (T)")+
  scale_x_continuous(breaks = 1:10)+
  coord_fixed(5/0.8)

legend <- get_legend(k_ks_plot_0.5+theme_bw(base_size = 20))
top_row <- plot_grid(k_ks_plot_0.5 + theme(legend.position="none"), 
                     T_ks_plot_0.5 + theme(legend.position="none"), 
                     labels = c("A.", "B."), label_size = 25)
bottom_row <- plot_grid(i_ks_plot_0.5 + theme(legend.position="none"), legend, labels = c("C.", ""),
                        label_size = 25)

plot_grid(
  top_row,
  bottom_row,
  ncol =1,
  align = 'v'
) +theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave("output_plots/si_ks_bar_plot_60.pdf",
       width = 13, height = 10.5, units = "in")