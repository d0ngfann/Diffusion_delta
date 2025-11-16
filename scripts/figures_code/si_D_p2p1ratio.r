# SI Figure 4 생성 (SI 섹션 D - 확산을 위한 최소 p2/p1 비율 분석)
###############################
# Figure S4 in SI Section D - Minimum p2/p1 ratio
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(ggplot2)
library(arrow)
library(cowplot)

### GET AGGREGATE FILES
dir = "output_data/cf/"
# k
df_k_big <- get_long_file(dir, "WS_.+_1T.+2td_.+_0.0sd_ES")
final_k <- df_k_big %>% filter(seed_strat == "seed_strat_one_nbrs")

write_parquet(final_k,"output_data/df_k_long_final.parquet")

# i
df_i_big <- get_long_file(dir, "WS_20k_.+_1T_.+_0.0sd_ES")
final_i <- df_i_big %>% filter(seed_strat == "seed_strat_one_nbrs")

write_parquet(final_i,"output_data/df_i_long_final.parquet")

#T
df_T_big <- get_long_file(dir, "WS_8k_.+_2td_.+_0.0sd_ES")
final_T <- df_T_big %>% filter(seed_strat == "seed_strat_one_nbrs") %>% filter(time_of_influence != 0)

write_parquet(final_T,"output_data/df_T_long_final.parquet")

### MAKE FIGURES
# load in data from bootstrap
sum_k <- read_csv("output_data/FigA4-BootstrappedCI-Degree.csv")
sum_T <- read_csv("output_data/FigA4-BootstrappedCI-Time.csv")
sum_i <- read_csv("output_data/FigA4-BootstrappedCI-Thrshld.csv")

#### Figure S4A
plot_k <- ggplot(sum_k %>% mutate(Margin = factor(Margin)), 
       aes(x = Degree, y = P2P1, color = Margin, group = Margin, ymin = P2P1 - ci, ymax = P2P1+ci, fill = Margin))+
  geom_ribbon(alpha = 0.3, color = NA)+
  #geom_point(size = 1)+
  geom_line()+
  labs(subtitle = "By Degree\nExample parameters i = 2, T = 1", 
       x = "Degree (k)", y = expression(paste( "Minimum ratio of ", p[2]," to ", p[1])))+
  theme_bw(base_size = 16)+
  scale_y_continuous(limits = c(1,18), breaks = 1:18)+
  scale_x_continuous(limits = c(4,20), breaks = seq(4,20,2)) +
  scale_color_manual(values = c(  "maroon1", "gray40","firebrick1","green3" ))+
  scale_fill_manual(values = c(  "maroon1", "gray40","firebrick1","green3" ))

#### Figure S4B
plot_T <- ggplot(sum_T %>% mutate(Margin = factor(Margin)), 
       aes(x = Time, y = P2P1, color = Margin, group = Margin, ymin = P2P1 - ci, ymax = P2P1+ci, fill = Margin))+
  geom_ribbon(alpha = 0.3, color = NA)+
  #geom_pointrange(size = 0.1)+
  geom_line()+
  labs(subtitle = "By Time of Influence\nExample parameters k = 8, i = 2", 
       x = "Time steps influential (T)", y = expression(paste( "Minimum ratio of ", p[2]," to ", p[1])))+
  theme_bw(base_size = 16)+
  scale_y_continuous(limits = c(1,18), breaks = 1:18)+
  scale_x_continuous(limits = c(1,10), breaks = seq(1,10,1)) +
  scale_color_manual(values = c(  "maroon1", "gray40","firebrick1","green3" ))+
  scale_fill_manual(values = c(  "maroon1", "gray40","firebrick1","green3" ))

#### Figure S4C
  plot_i <- ggplot(sum_i %>% mutate(Margin = factor(Margin)), 
      aes(x = Thrshld, y = P2P1, color = Margin, group = Margin, ymin = P2P1 - ci, ymax = P2P1 +ci, fill = Margin))+
  geom_ribbon(alpha = 0.3, color = NA)+
  geom_line()+
  labs(subtitle = "By Social Reinforcement Threshold\nExample parameters k = 20, T = 1", 
       x = "Social Reinforcement Threshold (i)", y = expression(paste( "Minimum ratio of ", p[2]," to ", p[1])))+
  theme_bw(base_size = 16)+
  scale_y_continuous(limits = c(1,18), breaks = 1:18)+
  scale_x_continuous(limits = c(2,10), breaks = seq(2,10,1)) +
  scale_color_manual(values = c(  "maroon1", "gray40","firebrick1","green3" ))+
  scale_fill_manual(values = c(  "maroon1", "gray40","firebrick1","green3" ))
  
grobs <- ggplotGrob(plot_i  + theme(legend.position = "bottom"))$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


plots <- plot_grid(plot_k + theme(legend.position = "none",
                         panel.grid.minor = element_blank()), 
          plot_T+ theme(legend.position = "none",
                        panel.grid.minor = element_blank()), 
          plot_i+ theme(legend.position = "none",
                        panel.grid.minor = element_blank()), ncol = 3, 
          labels = c("A.", "B.", "C."), label_size = 25)
plot_grid(plots, legend, ncol = 1, rel_heights = c(1, 0.2))
ggsave("/output_plots/p2_p1_ratio.pdf",
       width = 13, height = 7.5, units = "in")
