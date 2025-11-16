# SI Figure 7 생성 (SI 섹션 G - 1D vs 2D 클러스터 네트워크 비교)
###############################
# Figure S7 in SI Section G - Comparing 2D vs. 1D Clustered Networks
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(ggplot2)
library(cowplot)

### SAVE TO AGGREGATE FILE
save_long_file(dir = "output_data/cf/",
                   match_string = "_8k_.+_1T.+[2|4]td_.+_0.0sd_ES",
               summary_file =  "output_data/cf_1D_2D_sum.csv")

### READ AGGREGATE FILE
df_big_sum <- read_csv( "output_data/cf_1D_2D_sum.csv")

df_dim_compare <- df_big_sum %>% filter(seed_strat == "seed_strat_one_nbrs")%>%
  mutate(G_name = factor(case_when(G_name == "WS" ~ "Ring (1D)",
                                   G_name == "MR" ~ "Moore (2D)"),
                         ordered = TRUE, 
                         levels = c("Ring (1D)", 
                                    "Moore (2D)")),
         thrshld_name = paste0("i = ", thrshld))

### MAKE FIGURES
basesize = 16

#### Figure S7A - Spread on clustered networks
plot_c <- get_rand_clust_plot(df_dim_compare %>%
                      filter(perc_rewire== 0),
                    color = "dodgerblue")+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  labs(fill = "Proportion of Adopters")+ facet_grid(vars(thrshld_name), vars(G_name))+
  theme_bw(base_size = basesize)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position = "bottom",legend.key.width = unit(1, "cm"),legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))

#### Figure S7B - Spread on random rewired networks
plot_r <- get_rand_clust_plot(df_dim_compare %>%
                      filter(perc_rewire== 1),
                    color = "goldenrod1")+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  labs(fill = "Proportion of Adopters")+ facet_grid(vars(thrshld_name), vars(G_name))+
  theme_bw(base_size = basesize)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position = "bottom",legend.key.width = unit(1, "cm"),legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))

#### Figure S7C - Difference in spread in random vs. clustered networks
dif_plot_df <- get_dif_plot_df(df_dim_compare%>%
                                 select(-starts_with(c("t_60", "t_75", "t_90")), -ends_with(c("uci", "lci"))))

plot_dif <- get_dif_plot_spread(dif_plot_df)+ facet_grid(vars(thrshld_name), vars(G_name))+
  theme_bw(base_size = basesize)+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+labs(fill = "Difference in Spread\n(Random - Clustered)")+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position = "bottom",legend.key.width = unit(1, "cm"),legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))


#### Figure S7D - Speed on clustered networks
plot_t1 <- ggplot(data =df_dim_compare %>%
         filter(perc_rewire== 0),
       aes(x = p1, y = p2)) +
  geom_tile(aes(fill = t_total))+
  scale_fill_gradient( high = "dodgerblue", low = "gray95")+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  labs(fill = "Simulation Run Time")+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  theme_bw(base_size=basesize) +coord_fixed()+ facet_grid(vars(thrshld_name), vars(G_name))+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position = "bottom",legend.key.width = unit(1, "cm"),legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))

#### Figure S7D - Difference in speed
dif_plot_df_t <- df_dim_compare %>% 
  select(-ends_with(c("lci", "uci")))%>%
pivot_wider(names_from= perc_rewire, 
            values_from = c("maxspread_m","t_total", "t_90", "t_60", "t_75"))%>%
  mutate(mdif = maxspread_m_0 -maxspread_m_1
         t60dif = t_60_0 - t_60_1,
         t90dif = t_90_0 - t_90_1) 

plot_t2 <- ggplot(data =dif_plot_df_t,
       aes(x = p1, y = p2)) +
  geom_tile(aes(fill = t60dif))+
  scale_fill_gradient2(limits = c(-705, 705), high = "goldenrod1", mid = "gray95", low = "dodgerblue",
                       breaks = c(-700, -350, 0, 350, 700))+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  geom_hline(data= dif_plot_df_t %>% filter(thrshld == 2), aes(yintercept=1), colour="green3") + 
  geom_hline(data= dif_plot_df_t %>% filter(thrshld == 2), aes(yintercept=0.6), colour="red") + 
  geom_hline(data= dif_plot_df_t %>% filter(thrshld == 4), aes(yintercept=1), colour="purple1") + 
  geom_hline(data= dif_plot_df_t %>% filter(thrshld == 4), aes(yintercept=0.6), colour="cyan2") + 
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  labs(fill = "Difference in time to 60% spread\n(Random - Clustered)")+
  theme_bw(base_size=basesize) +coord_fixed()+ facet_grid(vars(thrshld_name), vars(G_name))+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position = "bottom",legend.key.width = unit(1, "cm"),legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))

#### Figure S7E - Difference in speed at specific points
df_plot_ts <- df_big_sum %>% 
  filter(p2 %in% c( 0.6, 1))%>%
  mutate(G_name_perc = case_when(
    G_name == "WS" & perc_rewire == 0 ~ "Ring (1D Lattice)",
    G_name == "WS" & perc_rewire == 1 ~ "Rewired Ring (Random)",
    G_name == "MR" & perc_rewire == 0 ~ "Moore (2D Lattice)",
    G_name == "MR" & perc_rewire == 1 ~ "Rewired Moore (Random)"),
    G_name_perc = factor(G_name_perc, ordered = TRUE,
                         levels = c("Ring (1D Lattice)", "Moore (2D Lattice)",
                                    "Rewired Ring (Random)", "Rewired Moore (Random)")),
    thrshld = paste0("i = ", thrshld),
    p2  =  paste0("p2 = ", p2))%>% 
  filter(seed_strat == "seed_strat_one_nbrs")%>%
  filter(!(G_name_perc == "Ring (1D Lattice)" & p1 < 0.1 ))
                                                                 
plot_t3 <- ggplot(df_plot_ts, 
       aes(x = p1, y= t_60, color = G_name_perc, group = G_name_perc,
           ymin = t_60_lci, ymax = t_60_uci, fill = G_name_perc))+geom_line()+ 
  geom_ribbon(alpha = 0.3, color = NA)+
  scale_color_manual(values = c("dodgerblue", "blue", "goldenrod1", "sienna1"))+
  scale_fill_manual(values = c("dodgerblue", "blue", "goldenrod1", "sienna1"))+
  scale_x_continuous( breaks = c(0, 0.5, 1))+
  facet_grid(cols =vars(p2), rows = vars(thrshld))+
  theme_bw(base_size = basesize)+
  labs( y = "Time steps to 60% spread",color = "Network Type", fill = "Network Type")+
  xlab(bquote(p[1]))+
  coord_fixed(1/450)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            data = df_plot_ts %>%filter(p2 == "p2 = 0.6")%>% filter(thrshld == "i = 2") %>% select(p2, thrshld)%>%distinct(), 
            colour = "red", size = 2, alpha = 1, fill = NA, inherit.aes = FALSE) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            data = df_plot_ts %>%filter(p2 == "p2 = 1")%>% filter(thrshld == "i = 2") %>% select(p2, thrshld)%>%distinct(), 
            colour = "green3", size = 2, alpha = 1, fill = NA, inherit.aes = FALSE) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            data = df_plot_ts %>%filter(p2 == "p2 = 0.6")%>% filter(thrshld == "i = 4") %>% select(p2, thrshld)%>%distinct(), 
            colour = "purple1", size = 2, alpha = 1, fill = NA, inherit.aes = FALSE) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            data = df_plot_ts %>%filter(p2 == "p2 = 1")%>% filter(thrshld == "i = 4") %>% select(p2, thrshld)%>%distinct(), 
            colour = "cyan2", size = 2, alpha = 1, fill = NA, inherit.aes = FALSE) +
  theme(legend.position = "bottom", legend.text = element_text(size = 10),legend.title=element_text(size=12),
        panel.border = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.key.width = unit(1, "cm"), 
        legend.direction='vertical') 

#### Combine plots and save
spread_plot <- plot_grid(plot_c, plot_r, plot_dif, ncol = 3,
                         labels = c("A.", "B.", "C."), label_size = 25)
time_plot <- plot_grid(plot_t1, plot_t2, plot_t3, ncol = 3, align = "h", axis = "b",
                       labels = c("D.", "E.", "F."), label_size = 25)

plot_grid(spread_plot, time_plot, ncol = 1)+
  theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave("output_plots/si_dim_compare_final.pdf", 
        width = 13, height = 13.5, units = "in")

