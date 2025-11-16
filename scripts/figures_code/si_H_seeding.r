# SI Figure 8 생성 (SI 섹션 H - 초기 시딩 전략 비교 분석)
###############################
# Figure S8 from SI Section H - Seeding Strategies
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(ggplot2)
library(cowplot)

### MAKE FIGURES
#### seeding is the same as 1d vs. 2d summary file
df_big_sum <- read_csv( "output_data/cf_1D_2D_sum.csv")

df_seed_compare <- df_big_sum %>%
  filter(!(perc_rewire == 1 & seed_strat == "seed_strat_adj"))%>%
  filter(!(G_name == "MR" & seed_strat == "seed_strat_adj"))%>%
  filter(!(G_name == "MR" & perc_rewire == 1))%>%
  mutate(G_name = ifelse(perc_rewire ==1, "ES", G_name))%>%
 pivot_wider(names_from= "seed_strat", 
                            values_from = c("maxspread_m", "t_60", "t_75", "t_90", "t_total")) %>%
  mutate(ss_nbrs_rand =  maxspread_m_seed_strat_one_nbrs - maxspread_m_seed_strat_random ,
         ss_rand_adj = maxspread_m_seed_strat_random - maxspread_m_seed_strat_adj,
         ss_nbrs_adj = maxspread_m_seed_strat_one_nbrs - maxspread_m_seed_strat_adj)%>%
  mutate(G_name = factor(case_when(G_name == "ES" ~ "Rewired Ring (Random)",
                            G_name == "WS" ~ "Ring (1D Clustered)",
                            G_name == "MR" ~ "Moore (2D Clustered)"),
                         ordered = TRUE, 
                         levels = c("Rewired Ring (Random)", 
                                    "Ring (1D Clustered)", 
                                    "Moore (2D Clustered)")),
         thrshld_name = paste0("i = ", thrshld))

#### Figure S8A
plot1 <- ggplot(data =df_seed_compare,
       aes(x = p1, y = p2)) +
  geom_tile(aes(fill = ss_nbrs_rand))+
  scale_fill_gradient2(limits = c(-1, 1),high = "green3", mid = "gray95", low = "firebrick1")+
  scale_y_continuous(limits = c(-0.01,1.01))+
  labs( title = "“Neighbor” vs. Random Seeding",
       subtitle = "Example Parameters k = 8, T = 1",
       fill = "Difference in Spread by Seeding Strategy\n(Neighbor - Random)")+ 
  facet_grid(rows = vars(thrshld_name), cols = vars(G_name))+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  theme_bw(base_size=18) +coord_fixed()+ 
  guides(fill = guide_colorbar(title.position = "top", ))+
  theme(legend.position= "bottom", legend.key.width = unit(1, "cm"), legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))

#### Figure S8B
plot2<- ggplot(data =df_seed_compare %>%
         filter(G_name== "Ring (1D Clustered)"),
       aes(x = p1, y = p2)) +
  geom_tile(aes(fill = ss_nbrs_adj))+
  scale_fill_gradient2(limits = c(-0.25,0.25),high = "green3", mid = "gray95", low = "firebrick1", 
                       breaks = c(-0.2, 0, 0.2))+
  scale_y_continuous(limits = c(-0.01,1.01))+
  labs( title = "“Neighbor” vs. Adjacent Seeding",
       subtitle = "Example Parameters k = 8, T = 1",
       fill = "Difference in Spread by Seeding Strategy\n(Neighbor - Adjacent)")+ 
  facet_grid(rows = vars(thrshld_name), cols = vars(G_name))+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  theme_bw(base_size=18) +coord_fixed()+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position= "bottom", legend.key.width = unit(1, "cm"), legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))

plot_grid(plot1, plot2, rel_widths = c(1, 0.6),ncol = 2, align = 'h', axis = "t", 
          labels = c("A.", "B."), label_size = 25)
ggsave("output_plots/si_seeding.pdf", 
       width = 13, height = 8, units = "in")