# SI Figure 1, 2 생성 (SI 섹션 B - 클러스터 네트워크에서의 분석적 경계 조건)
###############################
#Figure S1 and S2 from SI Section B - Analytical Boundary on Clustered Networks
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)

### CREATE AGGREGATE FILES
files1 = list.files(dir, pattern = "WS_+[_4|_8]k_.+_1T_.+_0.0sd_ES", full.names = TRUE)
files1 = list.files(dir, pattern = "WS_12k_.+_1T_.+_0.0sd_ES", full.names = TRUE)

file_list_main <-  c(files1, files2)

save_file <- "output_data/sum_main_df_long.csv"
i=1
for(file in file_list_main){
  df_temp <- read_csv(file) %>% 
    filter(p2 < 1.01)%>% 
    filter(seed_strat == "seed_strat_one_nbrs")%>%
    group_by(k, perc_rewire, sig_thresh_sd, n_nodes,
               time_of_influence, p1, p2, thrshld, seed_strat, G_name)%>%
    summarise(maxspread_m = mean(max_spread_norm),
              i_spread_m = mean(spread_time_1-spread_time_0))%>%
    mutate(over_thresh = ifelse(i_spread_m>=thrshld,  1,0))%>%
    ungroup()
    if(file.exists(save_file) == TRUE){
      existing_file <- read_csv(save_file)
      df_temp <- rbind(existing_file, df_temp)}
    write.csv(df_temp, save_file, row.names = FALSE) 
    print(i)
    i = i+1
}

### MAKE FIGURES
df_long_main<- read_csv("output_data/sum_main_df_long.csv")

df_sum_main <- df_long_main %>%
  pivot_wider(names_from = perc, 
              values_from = c("maxspread_m","i_spread_m", "over_thresh"))

#### Figure S1A
plot<-ggplot(data =df_sum_main%>%filter(i_k %in% c(0.25,0.5) ) %>%filter(`T`== 1),
       aes(x = p1, y = p2)) +
  geom_tile(aes(fill = factor(over_thresh)))+
  scale_fill_manual( values = c("gray95", "lightskyblue"))+
  geom_line(data = lbound_th  %>%filter(i_k %in% c(0.25,0.5)) %>% filter(t == 1),
aes(x = p1, y= p2), color = "blue" , size = 1)+
  scale_y_continuous(limits = c(-0.01,1.01))+
  labs(title = "Initial spread on clustered networks",
       fill = "Spread")+ facet_grid(vars(i_k), vars(k))+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  theme_bw(base_size=18) +coord_fixed()
ggsave("output_plots/si_cluster_net1.pdf",plot, width = 14, height = 7, units = "in")

#### Figure S1B
plot<-ggplot(data =df_sum_main%>%filter(i_k %in% c(0.25,0.5) ) %>%filter(`T`== 1),
             aes(x = p1, y = p2)) +
  geom_tile(aes(fill = maxspread_m_0))+
  scale_fill_gradient( low = "gray95", high = "dodgerblue")+
  geom_line(data = big_df2 %>%filter(i_k %in% c(0.25,0.5) ),aes(x = p1vec, y= p2vec), color = "blue" , size = 1)+
  scale_y_continuous(limits = c(-0.01,1.01))+
  labs(title = "Final spread on clustered networks",
       fill = "Spread")+ facet_grid(vars(i_k), vars(k))+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  theme_bw(base_size=18) +coord_fixed()
ggsave("output_plots/si_cluster_net2.pdf",plot, width = 14, height = 7, units = "in")

#### Figure S2
plot <- ggplot(df_long_main %>% filter(k == 8)%>%filter(`T`== 1), aes(x = i_spread_m, y= maxspread_m, color = factor(perc)))+
  scale_color_manual(values = c("dodgerblue", "goldenrod1"))+
  geom_point(alpha = 0.3)+
  scale_x_continuous(breaks = c(0,5,10,15,20,25))+
  geom_vline(aes(xintercept = thrshld), size = 0.75, linetype = "dashed")+
  labs(x = "Number of adopers in first time step", y = "Final Proportion of Spread", color = "Network Type" )+
  facet_wrap(~thrshld)+theme_bw(base_size=18)+coord_fixed(ratio = 25)

ggsave("output_plots/si_i_plot.pdf",plot, width = 18, height = 8, units = "in")