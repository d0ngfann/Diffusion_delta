# SI Figure 12-15 생성 (SI 섹션 L - Centola 2010 실험의 컴퓨터 시뮬레이션 복제)
###############################
# Figure S12 to S15 in SI Section L - Replication of Centola 2010
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(ggplot2)
library(stringr)
library(stats)
library(cowplot)

### FUNCTIONS SPECIFIC TO CENTOLA 2010 REPLICATION
get_expanded_ts <- function(df){
  names_vec <- c()
  for(i in 0:max(df$time_to_spread)){
    names_vec <- c(names_vec, paste0("t_", i))
  }
  
  df_time_series_expand <- df %>%
    mutate(full_timeseries = gsub("\\[|\\]" ,"", full_timeseries)) %>%
    separate_wider_delim(full_timeseries, 
                         delim = ", ",
                         names = names_vec,
                         too_few = "align_start")
  return(df_time_series_expand)
}

get_pairs_df <- function(df_ts, p_1, p_2, margin){
  df_temp <- df_ts %>% filter(p1 == p_1) %>%
    filter(p2 == p_2)%>% 
    pivot_longer(starts_with("t_"), names_to = "t", values_to = "spread")%>%
    group_by(perc_rewire, seed) %>%
    mutate(t = as.numeric(gsub("t_", "", t)),
           spread = as.numeric(spread)/128,
           spread = ifelse(is.na(spread) == TRUE, max_spread_norm, spread))%>%
    select(seed, perc_rewire, spread, t)%>%
    pivot_wider(names_from = perc_rewire, values_from = spread, names_prefix = "perc_")
  
  all_pairs_big <- data.frame()
  for(time in 1:25){
    df_temp2 <- df_temp %>% filter(t == time) 
    vec0 <-df_temp2 %>% pull(perc_0)
    vec1 <-df_temp2 %>% pull(perc_1)
    all_pairs<- expand.grid(perc_0 = vec0, 
                            perc_1 = vec1)%>%
      mutate(dif = perc_0 - perc_1, 
             mean1 = mean(dif),
             median = median(dif),
             lci = quantile(dif, 0.025),
             uci = quantile(dif, 0.975),
             #mean2 = mean(vec0 ,na.rm = T) - mean(vec1, na.rm = T),
             t = time, 
             #means_dif = mean1 - mean2,
             win0= ifelse(dif> margin, 1,0),
             win1= ifelse(dif< margin, 1,0),
             percent_win0 =sum(win0)/10000,
             percent_win1 = sum(win1)/10000)
    all_pairs_big <-rbind(all_pairs_big, all_pairs)
  }
  return(all_pairs_big)
}

### CREATE AND SAVE AGGREGATE FILES
df_big <- get_long_file(dir = "output_data/centola_2010/", 
                          match_string = "cf_sim_.+_[0.0|0.1]sd")
df_big_sum <- df_big %>%
    filter(p2 < 1.01)%>%
    group_by(k, perc_rewire, sig_thresh_sd, n_nodes,
             time_of_influence, p1, p2, thrshld, seed_strat, G_name) %>%
    summarise(maxspread_m = mean(max_spread_norm, na.rm = T),
              maxspread_lci = quantile(max_spread_norm, 0.025),
              maxspread_uci = quantile(max_spread_norm, 0.975),
              t_total = mean(time_to_spread, na.rm = T),
              t_total_lci = quantile(time_to_spread, 0.025),
              t_total_uci = quantile(time_to_spread, 0.975),
              t_60 = mean(time_to_60_spread, na.rm = T),
              t_60_lci = quantile(time_to_60_spread, 0.025, na.rm = T),
              t_60_uci = quantile(time_to_60_spread, 0.975, na.rm = T),
              t_75 = mean(time_to_75_spread, na.rm = T),
              t_75_lci = quantile(time_to_75_spread, 0.025, na.rm = T),
              t_75_uci = quantile(time_to_75_spread, 0.975, na.rm = T),
              t_90 = mean(time_to_90_spread, na.rm = T),
              t_90_lci = quantile(time_to_90_spread, 0.025, na.rm = T),
              t_90_uci = quantile(time_to_90_spread, 0.975, na.rm = T),
              sum_over_90 = sum(max_spread_norm>0.9),
              sum_over_75 = sum(max_spread_norm>0.75),
              sum_over_60 = sum(max_spread_norm>0.6))

write.csv(df_big_sum, "output_data/centola_2010_cf_sum.csv")
write.csv(df_big %>% filter(p1 %in% seq(0, 1, 0.1)) %>% filter(p2 %in% seq(0, 1, 0.1)),
          "output_data/centola_2010_df_ts_sum.csv")

### FIGURES 
#### Figure S12

df_big_sum = read.csv("output_data/centola_2010_cf_sum.csv")

dif_plot_df <- get_dif_plot_df(df_big_sum %>% 
                             select(-starts_with(c("t_60", "t_75", "t_90", "sum_over")), -ends_with(c("lci", "uci"))))%>%
  mutate(k_name = paste0("k = ", k),
         T_name = paste0("T = ", time_of_influence))

strip_margin = 2
strip_size = 18
plot1 <- get_dif_plot_spread(dif_plot_df %>% filter(sig_thresh_sd == 0))+
  labs(subtitle = expression(paste("Homogeneous Individuals (", sigma == 0, ", i = 2)")),
       fill = "")+ 
  facet_grid(vars(T_name), vars(k_name))+
  theme(plot.subtitle = element_text(size=14))+
scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))

plot2 <-  get_dif_plot_spread(dif_plot_df %>% filter(sig_thresh_sd == 0.1))+
  labs(subtitle = expression(paste("Heterogeneous Individuals (", sigma == 0.1, ", i = 2)")),
    fill = "Difference in Spread\n(Random - Clustered)")+ 
  facet_grid(vars(T_name), vars(k_name)) + 
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.direction="horizontal",legend.key.width = unit(1, "cm"),legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))

legend <- get_legend(plot2)
plot_grid(
  plot1 + theme(legend.position="none"),
  plot2+ theme(legend.position="none"),
  legend,
  align = 'vh',
  hjust = -1,
  nrow = 1, axis = "b",
  rel_widths = c(1,1,.5),
  labels = c("A.", "B.", ""), label_size = 25
)

ggsave("output_plots/si_c2010_main_compare.pdf", 
       width = 13, height = 6.5, units = "in")


#### Figure S13 - Full time series for sd = 0
full_ts_df <- read_csv("output_data/centola_2010_df_ts_sum.csv") %>%
  filter(time_of_influence == 1) %>%
  filter(k ==6)%>%
  filter(sig_thresh_sd %in% c(0, 0.1))
  
df_ts_plot_big <-  get_expanded_ts(full_ts_df) %>% 
 #filter(p1 %in% c(0, 0.2, 0.4, 0.6, 0.8, 1))%>% 
  #filter(p2 %in% c(0, 0.2, 0.4, 0.6, 0.8, 1))%>% 
  pivot_longer(starts_with("t_"), names_to = "t", values_to = "spread")%>%
  mutate(t = as.numeric(gsub("t_", "", t)),
         spread = as.numeric(spread)/128) %>% filter(is.na(spread) == FALSE)

strip_size = 16
strip_margin = 2
ggplot(df_ts_plot_big %>% #filter(t<26)%>%
         filter(sig_thresh_sd == 0)%>%
         mutate(perc_rewire = ifelse(perc_rewire == 0,
                                     "Clustered", "Random"),
                p1 = paste0("p1 = ",p1),
                p2 = paste0("p2 = ",p2)), 
       aes(x = t, y = spread, group = interaction(seed, perc_rewire), color = factor(perc_rewire)))+
  geom_line(alpha = 0.3, linewidth = 0.5)+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  facet_grid(rows = vars(p2), cols = vars(p1),
             as.table = FALSE)+
  scale_color_manual(values = c("dodgerblue", "goldenrod1"))+
  labs(color = "Network Type", y = "Fraction Adopted", x = "Time steps",
       title = expression(paste("Homogeneous Individuals, ", sigma == 0, ", k = 6, n = 128, i = 2, T = 1")))+
         
         #"Homogenous Individuals, k = 6, n = 128, i = 2, time of influence = 1, het = 0")+
  theme_bw(base_size = 16)+
  theme(strip.text.x = element_text(size = strip_size, margin = margin( b = strip_margin, t = strip_margin) ),
        strip.text.y = element_text(size = strip_size, margin = margin( l = strip_margin, r = strip_margin) ),
        legend.position = "bottom", legend.key.width = unit(1, "cm"))+ 
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))

ggsave("output_plots/si_c2010_ts_0sd.pdf", 
       width = 18, height = 13, units = "in")

#### Figure S14 - Full time series for sd = 0.1
ggplot(df_ts_plot_big %>% #filter(t<26)%>%
         filter(sig_thresh_sd == 0.1)%>%
         mutate(perc_rewire = ifelse(perc_rewire == 0,
                                     "Clustered", "Random"),
                p1 = paste0("p1 = ",p1),
                p2 = paste0("p2 = ",p2)), 
       aes(x = t, y = spread, group = interaction(seed, perc_rewire), color = factor(perc_rewire)))+
  geom_line(alpha = 0.3, linewidth = 0.5)+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  facet_grid(rows = vars(p2), cols = vars(p1),
             as.table = FALSE)+
  scale_color_manual(values = c("dodgerblue", "goldenrod1"))+
  labs(color = "Network Type", y = "Fraction Adopted", x = "Time steps",
       title = expression(paste("Heterogeneous Individuals, ", sigma == 0.1, ", k = 6, n = 128, i = 2, T = 1")))+
  theme_bw(base_size = 16)+
  theme(strip.text.x = element_text(size = strip_size, margin = margin( b = strip_margin, t = strip_margin) ),
        strip.text.y = element_text(size = strip_size, margin = margin( l = strip_margin, r = strip_margin) ),
        legend.position = "bottom", legend.key.width = unit(1, "cm"))+ 
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))

ggsave("output_plots/si_c2010_ts_0.1sd.pdf", 
       width = 18, height = 13, units = "in")

#### Figure S15 - Same network same seed test
dir = "output_data/centola_2010/"
df_rand <- read_csv(paste0(dir, "cf_sameR_6k_1T_2td_0.0sd0.2p1.csv"))%>%
  filter(row_number() <1001) %>%
  mutate(G_name = factor(paste("Trial",G_name+1),
                         ordered = TRUE, levels = paste("Trial", as.character(1:10))),
         group = "Same Network, Different Seed")

df_seed <- read_csv(paste0(dir, "cf_sameR_sameS_6k_1T_2td_0.0sd0.2p1.csv"))%>%
  mutate(G_name = factor(paste("Trial",G_name+1),
                         ordered = TRUE, levels = paste("Trial", as.character(1:10))),
         group = "Same Network, Same Seed")

df_difnet_difseed<- full_ts_df%>%
  select(- "...1")%>%
  filter(p1 == 0.2)%>%
  filter(p2 == 0.6)%>%
  filter(sig_thresh_sd == 0) %>%
  filter(perc_rewire == 1)%>%
  mutate(group = "Different Network, Different Seed")

df_ts_plot <- get_expanded_ts(rbind( df_rand, df_seed, df_difnet_difseed)) %>%
  pivot_longer(starts_with("t_"), names_to = "t", values_to = "spread")%>%
  mutate(t = as.numeric(gsub("t_", "", t)),
         spread = as.numeric(spread)/128)%>% filter(is.na(spread) == FALSE)

strip_size = 16
strip_margin = 2
plot1 <- ggplot(df_ts_plot %>% filter(group == "Different Network, Different Seed"), 
       aes(x = t, y = spread, group = seed))+
  geom_line(alpha = 0.3, linewidth = 0.5, color = "orange")+
  labs( y = "Fraction Adopted", x = "Time steps",
       subtitle = "Different Network, Different Seed",
       #subtitle = "p1 = 0.2, p2 = 0.6, k = 6, n = 128, i = 2, T = 1, het = 0"
       )+
  theme_bw(base_size = 16) +
  scale_x_continuous(limits = c(0,25))+
  coord_fixed(25)

plot2 <- ggplot(df_ts_plot %>%filter( group == "Same Network, Different Seed"), 
         aes(x = t, y = spread, group = seed))+
    geom_line(alpha = 0.3, linewidth = 0.5, color = "orange")+
    facet_wrap(~G_name, ncol = 2)+
    labs( y = "Fraction Adopted", x = "Time steps",
         subtitle = "Same Network, Different Seed",
         #subtitle = "p1 = 0.2, p2 = 0.6, k = 6, n = 128, i = 2, T = 1, het = 0"
         )+
    theme_bw(base_size = 16) + scale_x_continuous(limits = c(0,25))+
  coord_fixed(25)+
  theme(strip.text.x = element_text(size = strip_size, margin = margin( b = strip_margin, t = strip_margin) ),
        strip.text.y = element_text(size = strip_size, margin = margin( l = strip_margin, r = strip_margin) ))


plot3 <-  ggplot(df_ts_plot %>%filter( group == "Same Network, Same Seed"), 
       aes(x = t, y = spread, group = seed))+
  geom_line(alpha = 0.3, linewidth = 0.5, color = "orange")+
  facet_wrap(~G_name, ncol = 2)+
  labs(y = "Fraction Adopted", x = "Time steps",
       subtitle = "Same Network, Same Seed",
       #subtitle = "p1 = 0.2, p2 = 0.6, k = 6, n = 128, i = 2, T = 1, het = 0"
       )+
  theme_bw(base_size = 16) + scale_x_continuous(limits = c(0,25))+
  coord_fixed(25) +
  theme(strip.text.x = element_text(size = strip_size, margin = margin( b = strip_margin, t = strip_margin) ),
        strip.text.y = element_text(size = strip_size, margin = margin( l = strip_margin, r = strip_margin) ))

plot_grid(plot1, plot2, plot3, align = 'h',
          #hjust = -1,
          nrow = 1, axis = "b",
          rel_widths = c(1,1,1),
          labels = c("A.", "B.", "C."), label_size = 25)

ggsave("output_plots/si_c2010_rand_seed.pdf", 
       width = 13, height = 12, units = "in")