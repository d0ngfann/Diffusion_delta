# SI Figure 10, 11 생성 (SI 섹션 J - Rewiring 효과, 섹션 K - Centola & Macy 2007 복제)
###############################
# Figure S10 from SI Section J - Rewiring
# Figure S11 from SI Section K - Centola and Macy Replication
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(ggplot2)
library(scales)
library(cowplot)

### FIGURE SPECIFIC FUNCTIONS
get_mc_sigmoid <- function(x , k , m){
  denom = 1+(exp((k -x)*m))
  p = 1/denom
  return(p)
}

### CREATE AND SAVE AGGREGATE FILES

dir = "output_data/rewiring/"
files <- list.files(dir)
summary_file =  "output_data/rewiring_df_sum.csv"
save_long_file(dir,files, summary_file)

files_mr = list.files(dir, pattern = "MR_.+_2td_.+_0.0sd_ES")
save_long_file(dir, files_mr, summary_file = "output_data/cf_sum_MR_2td_ES.csv")

### MAKE REWIRING FIGURES (Figure S10 from SI Section J)
df_sum <- read_csv("output_data/rewiring_df_sum.csv")%>%
  filter(p2 != 0.75)%>%
  mutate(G_name = factor(ifelse(G_name == "WS", "Ring (1D)", "Moore (2D)"),
                         ordered =TRUE, levels = c("Ring (1D)", "Moore (2D)")),
         p2 = paste("p2 =", p2)) 


#### Figure S10A
plot_t1_s <- ggplot(df_sum%>% filter(time_of_influence == 1) %>% filter(p2>0),
            aes(x = perc_rewire, y = maxspread_m, color = p1, group = p1))+
  geom_line()+
  scale_color_stepsn(limits = c(-0.1, 0.2),colors = c("green3","olivedrab2", "gold", "maroon1", "purple"),
                     values = rescale(x = c(-0.1, 0.001, 0.01, 0.1, 0.2), from = c(-0.1, 0.2)),
                     breaks = c(0, 0.001, 0.01, 0.1, 0.2))+
  facet_grid(rows = vars(p2), cols = vars(G_name), as.table = FALSE) +
  labs(subtitle = "T = 1", x = "Percent edges rewired", y = "Proportion Adopted")+
  scale_y_continuous(limits = c(0,1),breaks = c(0, 0.5, 1))+
  scale_x_continuous(labels = scales::label_number(drop0trailing=TRUE))+
  theme_bw(base_size = 16)

#### Figure S10B
plot_t0_s <- ggplot(df_sum%>% filter(time_of_influence == 0) %>% filter(p2>0),
                    aes(x = perc_rewire, y = maxspread_m, color = p1, group = p1))+
  geom_line()+
  scale_color_stepsn(limits = c(-0.1, 0.2),colors = c("green3","olivedrab2", "gold", "maroon1", "purple"),
                     values = rescale(x = c(-0.1, 0.001, 0.01, 0.1, 0.2), from = c(-0.1, 0.2)),
                     breaks = c(0, 0.001, 0.01, 0.1, 0.2))+
  facet_grid(rows = vars(p2), cols = vars(G_name), as.table = FALSE) +
  labs(subtitle = "Adopters influential for entire simulation", x = "Percent edges rewired", y = "Proportion Adopted")+
  scale_y_continuous(limits = c(0,1),breaks = c(0, 0.5, 1))+
  scale_x_continuous(labels = scales::label_number(drop0trailing=TRUE))+
  theme_bw(base_size = 16)

#### Figure S10C
plot_t1_t <-ggplot(df_sum%>% filter(time_of_influence == 1), 
                   aes(x = perc_rewire, y = t_90, color = p1, group = p1))+
  geom_line()+
  scale_x_log10(labels = scales::label_number(drop0trailing=TRUE))+
  scale_color_stepsn(limits = c(-0.1, 0.2),colors = c("green3","olivedrab2", "gold", "maroon1", "purple"),
                     values = rescale(x = c(-0.1, 0.001, 0.01, 0.1, 0.2), from = c(-0.1, 0.2)),
                     breaks = c(0, 0.001, 0.01, 0.1, 0.2))+
  facet_grid(rows = vars(p2), cols = vars(G_name), as.table = FALSE) +
  labs(subtitle = "T = 1", x = "Percent edges rewired", y = "Time steps 90% to spread")+
  theme_bw(base_size = 16)
  
#### Figure S10D
plot_t0_t <- ggplot(df_sum%>% filter(time_of_influence == 0), 
                    aes(x = perc_rewire, y = t_90, color = p1, group = p1))+
  geom_line()+
  scale_x_log10(labels = scales::label_number(drop0trailing=TRUE))+
  scale_color_stepsn(limits = c(-0.1, 0.2),colors = c("green3","olivedrab2", "gold", "maroon1", "purple"),
                        values = rescale(x = c(-0.1, 0.001, 0.01, 0.1, 0.2), from = c(-0.1, 0.2)),
                     breaks = c(0, 0.001, 0.01, 0.1, 0.2))+
  facet_grid(rows = vars(p2), cols = vars(G_name), as.table = FALSE) +
  labs(subtitle = "Adopters influential for entire simulation", x = "Percent edges rewired", y = "Time steps to 90% spread")+
  theme_bw(base_size = 16)

#### Figure S10E
df_map_raw <- rbind(read_csv("output_data/cf_sum_ESrand.csv") %>%
                      select(-i_k),
                    read_csv("output_data/cf_sum_MR_2td_ES.csv"))%>%
  filter(k ==8) %>% filter(thrshld == 2)%>%
  select(-ends_with(c("uci", "lci", "_60", "_75", "_90", "filename"))) %>% 
  filter(seed_strat == "seed_strat_one_nbrs")%>%
  filter(sig_thresh_sd == 0)%>%
  filter(time_of_influence != 5) %>% select(-n_nodes)

df_map <- rbind(df_map_raw, 
                df_map_raw %>% filter(G_name == "WS")%>% 
                  filter(perc_rewire == 1)%>% filter(time_of_influence == 0) %>% mutate(G_name = "MR") )%>%
  mutate(G_name = factor(ifelse(G_name == "WS", "Ring (1D)", "Moore (2D)"),
                         ordered =TRUE, levels = c("Ring (1D)", "Moore (2D)")))

df_points <- expand.grid(p1 = c(0, 0.001, 0.01, 0.1, 0.2), p2 = c(0.25, 0.5, 1))%>%
  filter(p1 <= p2)

dif_df0 <- get_dif_plot_df(df_map)
plot_map1 <- get_dif_plot_spread(dif_df0 %>% filter(time_of_influence == 1))+
  geom_point(data = df_points, aes(x = p1, y = p2, color = p1))+
  scale_color_stepsn(limits = c(-0.1, 0.2),colors = c("green3","olivedrab2", "gold", "maroon1", "purple"),
                     values = rescale(x = c(-0.1, 0.001, 0.01, 0.1, 0.2), from = c(-0.1, 0.2)),
                     breaks = c(0, 0.001, 0.01, 0.1, 0.2), labels = scales::label_number(drop0trailing=TRUE)) + 
  scale_fill_gradient2( limits = c(-1,1),high = "goldenrod1", low ="dodgerblue", mid = "gray95", guide = "none")+
  labs(title = "T = 1", color =  expression(paste( "Base rate adoption, ", p[1] ))  )+
  facet_grid(~ G_name)+theme_bw(base_size = 16)+
  guides(color = guide_colorsteps(title.position = "top"))+
  theme(legend.position= "bottom", legend.key.width = unit(1.5, "cm"), legend.title.align=0.5,
        legend.title=element_text(size=12), legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

#### Figure S10F
plot_map2 <- get_dif_plot_spread(dif_df0 %>% filter(time_of_influence == 0))+
  geom_point(data = df_points, aes(x = p1, y = p2, color = p1))+
  scale_color_stepsn(limits = c(-0.1, 0.2),colors = c("green3","olivedrab2", "gold", "maroon1", "purple"),
                     values = rescale(x = c(-0.1, 0.001, 0.01, 0.1, 0.2), from = c(-0.1, 0.2)),
                     breaks = c(0, 0.001, 0.01, 0.1, 0.2), labels = scales::label_number(drop0trailing=TRUE),
                     guide = "none") + 
  labs(subtitle = "Adopters influential for entire simulation",
       fill = "Difference in Spread (Random - Clustered)")+
  facet_grid(~ G_name)+theme_bw(base_size = 16)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position= "bottom", legend.key.width = unit(1.5, "cm"), legend.title.align=0.5,
        legend.title=element_text(size=12), legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

#### Combine and save figure
top <- plot_grid(plot_t1_s + theme(legend.position = "none"),
                 plot_t0_s  + theme(legend.position = "none"), 
                 plot_t1_t  + theme(legend.position = "none"),
                 plot_t0_t  + theme(legend.position = "none"), 
                 labels = c("A.", "B.", "C.", "D."),
                 label_size = 25,
                   ncol = 2)

bottom <- plot_grid(plot_map1, plot_map2, labels = c("E.", "F."),
                    label_size = 25)

plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 0.5))

ggsave("output_plots/si_rewiring.pdf", width = 13, height = 16.75, units = "in")

### MAKE FIGURES FOR CENTOLA & MACY REPLICATION (Figure S11 from SI Section K) 

df_mc_raw<- get_long_file(dir = "output_data/macy_centola_07/",
                   match_string = NULL)
df_sum_mc <- df_mc_raw %>%
  group_by(k, perc_rewire, sig_thresh_sd, n_nodes, m,
           time_of_influence, p1, p2, thrshld, seed_strat, G_name) %>%
  summarise(maxspread_m = mean(max_spread_norm, na.rm = T),
            maxspread_lci = quantile(max_spread_norm, 0.025),
            maxspread_uci = quantile(max_spread_norm, 0.975),
            t_90 = mean(time_to_90_spread, na.rm = T),
            t_90_lci = quantile(time_to_90_spread, 0.025, na.rm = T),
            t_90_uci = quantile(time_to_90_spread, 0.975, na.rm = T),
            t_total = mean(time_to_spread, na.rm = T),
            t_total_lci = quantile(time_to_spread, 0.025),
            t_total_uci = quantile(time_to_spread, 0.975))%>%
  mutate(G_name = factor(ifelse(G_name == "WS", "Ring (1D)", "Moore (2D)"),
                         ordered =TRUE, levels = c("Ring (1D)", "Moore (2D)"))) 

#### Figure S11A
df_sigmoid <- expand.grid(x = seq(0,5, 0.1),
                          k = 2,
                          m = 1:10)%>%
  mutate(p = get_mc_sigmoid(x, k, m),
         red = ifelse(m %in% c(2, 4, 10), "use", "not use"),
         text = ifelse(x == 3.5 & m %in% c(1,2,10), paste("m = ", m, sep = ""), NA))

plot_sigmoid <- ggplot(df_sigmoid, aes(x = x, y = p, group = m, color = m, label = text))+
  scale_color_gradient2(low= "maroon1", high = "green3", mid = "goldenrod2", midpoint = 5.5,
                        breaks = 1:10)+
  geom_line()+theme_bw(base_size = 18)+
  coord_fixed(5)+
  guides(color = guide_colorbar(title.position = "top"))+
  labs(x = "Num. exposures to different neighbors", y = "Probability of adoption", 
       color = "Slope parameter, m")+
  theme(legend.direction= "horizontal", legend.key.width = unit(1, "cm"), legend.title.align=0.5,
        legend.title=element_text(size=12))

#### Figure S11B
plot_mc_10 <-ggplot(df_sum_mc %>% filter(time_of_influence == 0) %>% filter(m< 11), 
                 aes(x = perc_rewire, y = t_90, group = m, color = m))+
  scale_color_gradient2(low= "maroon1", high = "green3", mid = "goldenrod2", midpoint = 5.5)+
  scale_x_log10(limits = c(0.001, 6), labels = scales::label_number(drop0trailing=TRUE))+
  geom_line()+
  geom_text(aes(color = m, label = paste("m =", m)), nudge_x = 0.5, nudge_y = -0.5,
            data = df_sum_mc %>% filter(time_of_influence == 0) %>% filter(m > 7)%>% filter(perc_rewire == 1))+
  facet_wrap(~G_name, scales = "free_y")+
  theme_bw(base_size = 18)+
  labs( x = "Percent edges rewired", y = "Time steps 90% to spread")+
  theme(legend.position = "none")

#### Figure S11C
plot_mc_5 <-ggplot(df_sum_mc %>% filter(time_of_influence == 0) %>% filter(m< 6), 
                 aes(x = perc_rewire, y = t_90, group = m, color = m))+
  scale_color_gradient2(limits = c(1,10),low= "maroon1", high = "green3", mid = "goldenrod2", midpoint = 5.5)+
  geom_text(aes(color = m, label = paste("m =", m)), nudge_x = 0.5, nudge_y = -0.5,
            data = df_sum_mc %>% filter(time_of_influence == 0) %>% filter(m %in% c(3,4,5))%>% filter(perc_rewire == 1))+
  scale_x_log10(limits = c(0.001, 6), labels = scales::label_number(drop0trailing=TRUE))+
  geom_line()+
  facet_wrap(~G_name, scales = "free_y")+
  theme_bw(base_size = 18)+
  labs( x = "Percent edges rewired", y = "Time steps 90% to spread")+
  theme(legend.position = "none")

#### Figure S11D
plot_mc_compare <- ggplot()+
  scale_color_gradient2(limits = c(1,10),low= "maroon1", high = "green3", mid = "goldenrod2", midpoint = 5.5)+
  scale_x_log10(limits = c(0.001, 6), labels = scales::label_number(drop0trailing=TRUE))+
  #scale_y_log10()+
  geom_line(data = df_sum_mc %>% filter(time_of_influence == 0) %>% filter(m %in% c(2,4)), 
            aes(x = perc_rewire, y = t_90, group = m, color = m))+
  geom_text(aes(x = perc_rewire, y = t_90, color = m, label = paste("m =", m)), nudge_x = 0.5, nudge_y = -0.5,
            data = df_sum_mc %>% filter(time_of_influence == 0) %>% filter(m %in% c(2,4)) %>% filter(perc_rewire == 1))+
  geom_line(data = df_sum%>% filter(time_of_influence == 0) %>% filter(p1 %in% c(0.01, 0.1)) 
            %>% filter(p2 == "p2 = 0.5"),
            aes(x = perc_rewire, y = t_90, group = p1), color = "black", linetype = "dashed")+
  facet_wrap(~G_name, scales = "free_y")+
  theme_bw(base_size = 18)+
  labs( x = "Percent edges rewired", y = "Time steps 90% to spread")+
  theme(legend.position = "none")

plot_grid(plot_sigmoid, plot_mc_10,plot_mc_5, plot_mc_compare,
          labels = c("A.", "B.", "C.", "D."), label_size = 25, vjust = 1)
ggsave("output_plots/si_mc07.pdf", width = 13, height = 9, units = "in")
