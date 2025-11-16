# 메인 Figure 2 및 SI 섹션 E, F 그림 생성 (확산 속도 및 영향력 시간 분석)
###############################
# Figure 2
#Figure S5 from SI Section E - Speed of Spread
#Figure S6 from SI Section F - Time of Influence
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)

### CREATE AND SAVE AGGREGATE FILES
dir = "output_data/cf/"
files1 = list.files(dir, pattern = "WS_+[_4|_8]k_.+_[0|1|5]T_.+_0.0sd_ES")
files2 = list.files(dir, pattern = "WS_12k_.+_[0|1|5]T_.+_0.0sd_ES")

final_files <- c(files1, files2)
save_long_file(dir, final_files, 
               summary_file = "output_data/cf_sum_ESrand.csv")

### MAKE FIGURES
#### Figure 2A

df_sum1 <- read_csv("output_data/cf_sum_ESrand.csv")%>% 
  distinct()%>%
  filter(grepl("_ES", filename) == TRUE) %>%
  filter(seed_strat == "seed_strat_one_nbrs") %>%
  filter(sig_thresh_sd == 0) %>%
  filter(G_name == "WS") %>%
  filter(k %in% c(4,8,12))%>%
  filter(time_of_influence %in% c(0,5, 1))%>%
  select(-ends_with(c("_90", "_75", "_90_lci", "_75_lci", "_75_uci", "_90_uci")), - filename)
  
df_sum <- df_sum1%>%
  pivot_wider(names_from = perc_rewire, 
              values_from = c("maxspread_m","maxspread_lci","maxspread_uci", "t_total","t_total_lci",
                             "t_total_uci", "t_60", "t_60_lci", "t_60_uci", "sum_over_60"))%>%
  mutate(mdif =maxspread_m_1-maxspread_m_0,
         tdif = t_total_0 -  t_total_1,
         t60dif = t_60_0 - t_60_1,
         i_k = thrshld/k) %>% 
  filter(i_k %in% c(0.25,0.5) ) %>%
  mutate(i_k = paste0("i/k = ",thrshld/k), 
         k = factor(paste0("k = ", k), ordered = TRUE,
                    levels = c("k = 4", "k = 8", "k = 12")))

bounds_df <- get_bounds_df(c(4,8,12), c(2,3,4,6), t_val = c(1,5))%>%
  #filter(i_k %in% c(0.25,0.5)) %>%
  mutate(i_k = paste0("i/k = ",i/k), 
         k = factor(paste0("k = ", k), ordered = TRUE,
                    levels = c("k = 4", "k = 8", "k = 12")))

bounds_df_theory <- bounds_df %>%
  filter(t == 1) %>% 
  filter(k %in% c("k = 4", "k = 8"))%>%
  filter(i %in% c(2,3,4)) %>%
  mutate( i = paste0("i = ", i)) 

p1_val <- 1- ((1-(1/(4-1)))^(1/1))
bounds_df_shaded <-bounds_df %>%
  filter(i == 2)%>%
  filter(k == "k = 4") %>%
  filter(t ==1) %>%
  filter(p1 <= p2)

bounds_df_blue<- bounds_df_shaded %>%
  filter(p1 < p1_val)%>%
  select(p1, p2) %>%
tibble::add_row(p1 = c(p1_val, 0), p2 = c(1, 1))

bounds_df_orange<- bounds_df_shaded %>% 
  filter(p1 > p1_val)%>%
  select(p1, p2) %>%
  tibble::add_row(p1 = c(p1_val), p2 = c(p1_val))

bounds_df_dark_gray<- bounds_df_shaded %>% 
  filter(p1 < p1_val)%>%
  select(p1, p2) %>%
  tibble::add_row(p1 = c( p1_val,0,0), p2 = c(p1_val, 0, 1))

bounds_df_light_gray<- bounds_df_shaded %>% 
  filter(p1 > p1_val)%>%
  select(p1, p2) %>%
  tibble::add_row(p1 = c(1,p1_val), p2 = c(1, 1))

plot1 <- ggplot()+
  #geom_tile( data = bounds_df_shaded, aes(x = p1, y = p2, fill = region_type), alpha = 0.5)+
  #scale_fill_manual(values = c("dodgerblue", "goldenrod1", "gray30","gray80", "cornsilk"))+
  geom_polygon(data = bounds_df_blue,
               aes(x = p1, y = p2),
               fill = "dodgerblue",
               alpha = 0.5)+
  geom_polygon(data = bounds_df_dark_gray,
               aes(x = p1, y = p2),
               fill = "gray30",
               alpha = 0.5)+
  geom_polygon(data = bounds_df_light_gray,
               aes(x = p1, y = p2),
               fill = "gray80",
               alpha = 0.5)+
  geom_polygon(data = bounds_df_orange,
               aes(x = p1, y = p2),
               fill = "goldenrod1",
               alpha = 0.5)+
  geom_polygon(data = data.frame(x = c(0, 1, 1), y = c(0, 1, 0)),
               aes(x = x, y = y), fill = 'cornsilk', alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")+
  geom_line(data = bounds_df_theory %>% filter(k == "k = 4"), 
            aes(x = p1, y= p2), color = "blue", linewidth = 1)+
  geom_vline(data = bounds_df_theory %>% filter(k == "k = 4"),
             aes(xintercept =p1_star ), color = "orange", linewidth = 1)+
  scale_x_continuous(limits = c(-0.01,1.01))+
  scale_y_continuous(limits = c(-0.01,1.01))+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  labs(subtitle = "Analytical Regions of Spread\nExample Parameters k = 4, i = 2, T = 1")+
  theme_bw(base_size=16) +coord_fixed()+
  theme(legend.position = "none")

#### Figure 2B
plot2 <- ggplot()+
  geom_line(data = bounds_df_theory, 
            aes(x = p1, y= p2, color = i, group = i), linewidth = 1)+
  geom_vline(data = bounds_df_theory,
             aes(xintercept =p1_star ), color = "orange", linewidth = 1)+
  scale_x_continuous(limits = c(-0.01,1.01))+
  scale_y_continuous(limits = c(-0.01,1.01))+
  scale_color_manual(values = c("blue", "dodgerblue", "lightskyblue1"))+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  theme_bw(base_size=16) +coord_fixed()+
  facet_wrap(~k)+
  labs(color = "", subtitle = "Analytical Regions by k and i, T = 1")+
  theme(legend.position = "bottom",legend.direction = "vertical",
        legend.text=element_text(size=16),
        strip.text.x = element_text(size =16, margin = margin(0.1,0,0.1,0, "cm")))

#### Figure 2C and Figure S6A
plot_t1_s <- get_dif_plot_spread(df_sum %>%
                      filter(time_of_influence== 1))+
  geom_line(data = bounds_df  %>% filter(i_k %in% c("i/k = 0.25", "i/k = 0.5")) %>%filter(t == 1),
            aes(x = p1, y= p2), color = "blue", linewidth = 1 )+
  geom_vline(data = bounds_df %>% filter(i_k %in% c("i/k = 0.25", "i/k = 0.5")) %>%filter(t == 1),
             aes(xintercept =p1_star ), color = "orange", linewidth = 1)+
  geom_polygon(data = data.frame(x = c(0,0, 1, 1), y = c(0, 1,1, 0), k = rep("k = 4",4), i_k = rep("i/k = 0.25",4) ) %>%
                 mutate(k = factor(k, ordered = TRUE,
                                   levels = c("k = 4", "k = 8", "k = 12"))),
               aes(x = x, y = y), fill = 'white')+
  labs(subtitle = "Empirical Regions of Spread, T = 1",
       fill = "Difference in Spread\n(Random - Clustered)")+ facet_grid(vars(i_k), vars(k))+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme_bw(base_size = 16)+
  theme(legend.position= "bottom", legend.key.width = unit(1, "cm"), 
        #plot.title = element_text(size=16), 
        #plot.subtitle = element_text(size=14),
        legend.title=element_text(size=16),
        legend.title.align=0.5,
        strip.text.x = element_text(size =16, margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(size =16, margin = margin( 0,0.1,0,0.1,"cm")))

#### Combine and save Figure 2
toprow <- plot_grid(plot1, 
                    plot2, 
                    ncol = 1,
                    labels = c("A.", "B."), label_size = 24, align = "v", axis = "l", rel_heights = c(1, 0.86))

bottom_row <- plot_grid(plot_t1_s,
                        labels = c( "C."), label_size = 24)

plot_grid(toprow, bottom_row, ncol = 2, align = "h", axis = "t", rel_widths = c(1, 1.5))
ggsave("output_plots/main_fig2.pdf", width = 14, height = 11, units = "in")


#### Figure S6C
plot_t5_s <- get_dif_plot_spread(df_sum %>%
                                   filter(time_of_influence== 5))+
  geom_line(data = bounds_df%>% filter(i_k %in% c("i/k = 0.25", "i/k = 0.5"))  %>% filter(t == 5),
            aes(x = p1, y= p2), color = "blue", linewidth = 1 )+
  geom_vline(data = bounds_df%>% filter(i_k %in% c("i/k = 0.25", "i/k = 0.5"))  %>% filter(t == 5),
             aes(xintercept =p1_star ), color = "orange", linewidth = 1)+
  labs(title = "",
       fill = "")+ facet_grid(vars(i_k), vars(k))+
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm"))

#### Figure S6E
plot_t0_s <-get_dif_plot_spread(df_sum  %>%
                                  filter(time_of_influence== 0))+
  labs(title = "",
       fill = "Difference in Spread\n(Random - Clustered)")+ facet_grid(vars(i_k), vars(k))+
  guides(fill = guide_colorbar(title.position = "top", ))+
  theme(legend.position= "bottom", legend.key.width = unit(1, "cm"), legend.title=element_text(size=14),
        legend.title.align=0.5)

grobs <- ggplotGrob(plot_t0_s)$grobs
legend_s <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

#### Figure S6B
df_plot_hist <- df_sum1 %>%
  mutate(i_k = thrshld/k )%>%
  filter(i_k %in% c(0.25,0.5)) %>%
  mutate(i_k = paste0("i/k = ",i_k), 
         k = factor(paste0("k = ", k), ordered = TRUE,
                    levels = c("k = 4", "k = 8", "k = 12")))

plot_t1_h <-ggplot(data =df_plot_hist %>%
         filter(time_of_influence== 1),
       aes(x = maxspread_m)) +
  geom_histogram(aes(y=after_stat(count)/sum(..count..)),position = "identity", fill = "gray75", color = "black")+ 
  facet_grid(vars(i_k), vars(k))+theme_bw(base_size =18)+
  labs(y = expression(paste( "Proportion of ", p[1] <= p[2], " space")), x = "Proportion of Adopted")+
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.1, 0.2))+
  coord_fixed(5)

#### Figure S6D
plot_t5_h <- ggplot(data =df_plot_hist %>%
                    filter(time_of_influence== 5),
                  aes(x = maxspread_m)) +
  geom_histogram(aes(y= after_stat(count)/sum(..count..)),position = "identity", fill = "gray75", color = "black")+ 
  facet_grid(vars(i_k), vars(k))+theme_bw(base_size =18)+
  labs(y = expression(paste( "Proportion of ", p[1] <= p[2], " space")), x = "Proportion of Adopted")+
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.1, 0.2))+
  coord_fixed(5)

#### Figure S6F
plot_t0_h <- ggplot(data =df_plot_hist %>%
                      filter(time_of_influence== 0),
                    aes(x = maxspread_m)) +
  geom_histogram(aes(y=after_stat(count)/sum(..count..)),position = "identity", fill = "gray75", color = "black")+ 
  facet_grid(vars(i_k), vars(k))+theme_bw(base_size =18)+
  labs(y = expression(paste( "Proportion of ", p[1] <= p[2], " space")), x = "Proportion of Adopted")+
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.1, 0.2))+
  coord_fixed(5)

#### Combine and save Figure S6
plot_grid(plot_t1_s+
            scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
            scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
            theme(legend.position = "none")+labs(title = "", subtitle = "T = 1"), 
          
          plot_t1_h+
            scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+labs(subtitle = "T = 1"), 
          
          plot_t5_s +scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
            scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
            theme(legend.position = "none")+labs(subtitle = "T = 5"), 
          
          plot_t5_h +
            scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+labs(subtitle = "T = 5"), 
          
          plot_t0_s +scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
            scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
            theme(legend.position = "none")+labs(subtitle = "Adopters influential for entire simulation"), 
          
          plot_t0_h +
            scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))
          +labs(subtitle = "Adopters influential for entire simulation"),
          
          legend_s, 
          rel_heights = c(1,1,1,0.4),
          labels = c("A.", "B.", "C.", "D.", "E.", "F.", "", ""),
          label_size = 25, 
          ncol = 2,  align = "hv", axis = "t")+
  theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave("output_plots/si_time_inf.pdf", width = 13, height = 16, units = "in")


df_sum_t <- df_sum %>% filter(sum_over_60_1 == 100) %>% filter(sum_over_60_0 == 100)
max_tdif = max(df_sum_t$t60dif, na.rm = TRUE)
min_tdif = min(df_sum_t$t60dif, na.rm = TRUE)
basesize = 16

plot_t1_t <- ggplot(data =df_sum %>%
                      filter(time_of_influence== 1),
                    aes(x = p1, y = p2)) +
  geom_tile(aes(fill = t60dif))+
  scale_fill_gradient2(high = "goldenrod1", mid = "gray95", low = "dodgerblue")+
  #scale_fill_gradient2(limits = c(-705, 705), high = "goldenrod1", mid = "gray95", low = "dodgerblue",
  #breaks = c(-700, -350, 0, 350, 700))+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  labs(subtitle = "T = 1", fill = "Difference in time to 60% spread (Random - Clustered)")+
  theme_bw(base_size=basesize) +coord_fixed()+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+ 
  facet_grid(vars(i_k), vars(k))+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position= "bottom", legend.key.width = unit(1, "cm"), legend.title=element_text(size=12),
        legend.title.align=0.5, legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

plot_t5_t <- ggplot(data =df_sum %>%
         filter(time_of_influence== 5),
       aes(x = p1, y = p2)) +
  geom_tile(aes(fill = t60dif))+
  scale_fill_gradient2(high = "goldenrod1", mid = "gray95", low = "dodgerblue")+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  labs(subtitle = "T = 5", fill = "Difference in time to 60% spread (Random - Clustered)")+
  theme_bw(base_size=basesize) +coord_fixed()+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+ 
  facet_grid(vars(i_k), vars(k))+
  guides(fill = guide_colorbar(title.position = "top" ))+
  theme(legend.position= "bottom", legend.key.width = unit(1, "cm"), legend.title=element_text(size=12),
        legend.title.align=0.5, legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

plot_t0_t <-ggplot(data =df_sum %>%
         filter(time_of_influence== 0),
       aes(x = p1, y = p2)) +
  geom_tile(aes(fill = t60dif))+
  scale_fill_gradient2( high = "goldenrod1", mid = "gray95", low = "dodgerblue")+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  labs(subtitle = "Adopters influential for entire simulation",fill = "Difference in time to 60% spread (Random - Clustered)")+
  theme_bw(base_size=basesize) +coord_fixed()+
  scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+
  scale_x_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1))+ 
  facet_grid(vars(i_k), vars(k))+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.position= "bottom", legend.key.width = unit(1, "cm"), legend.title=element_text(size=12),
        legend.title.align=0.5, legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

df_sum_t <- df_sum1  %>% 
  filter(p2 ==0.6)%>% 
  filter(sum_over_60 > 10)%>% 
  mutate(perc_rewire = ifelse(perc_rewire == 0, "Clustered Network", "Random Network"),
         i_k = paste0("i/k = ",thrshld/k),
         k = factor(paste0("k = ", k), ordered = TRUE,
                    levels = c("k = 4", "k = 8", "k = 12")))%>%
  filter(i_k %in% c("i/k = 0.25", "i/k = 0.5"))

max_t <- max(df_sum_t$t_60_uci)
min_t <- min(df_sum_t$t_60_lci)
                                                                          
plot_t1_tl <- ggplot(df_sum_t  %>%  filter(time_of_influence ==1),
       aes(x = p1, y= t_60, color = perc_rewire, group = perc_rewire,
           ymin = t_60_lci, ymax = t_60_uci, fill =perc_rewire))+geom_line()+
  geom_ribbon(alpha = 0.3, color = NA)+ 
  scale_color_manual(values = c("dodgerblue",  "sienna1"))+
  scale_fill_manual(values = c("dodgerblue",  "goldenrod1"))+
  scale_y_continuous(limits = c(min_t, max_t))+
  facet_grid(cols =vars(k), rows = vars(i_k))+
  xlab(bquote(p[1]))+
  labs(subtitle = "T = 1", color = "", fill = "", y = "Time to 60% Spread")+
  coord_fixed(0.6/max_t)+
  theme_bw(base_size = basesize)

plot_t5_tl <-ggplot(df_sum_t  %>%  filter(time_of_influence ==5),
       aes(x = p1, y= t_60, color = perc_rewire, group = perc_rewire,
           ymin = t_60_lci, ymax = t_60_uci, fill =perc_rewire))+geom_line()+
  geom_ribbon(alpha = 0.3, color = NA)+ 
  scale_color_manual(values = c("dodgerblue",  "sienna1"))+
  scale_fill_manual(values = c("dodgerblue",  "goldenrod1"))+
  scale_y_continuous(limits = c(min_t, max_t))+
  facet_grid(cols =vars(k), rows = vars(i_k))+
  xlab(bquote(p[1]))+
  labs(subtitle = "T = 5", color = "", fill = "", y = "Time to 60% Spread")+
  coord_fixed(0.6/max_t)+
  theme_bw(base_size = basesize)

plot_t0_tl <-ggplot(df_sum_t  %>%  filter(time_of_influence ==0),
       aes(x = p1, y= t_60, color = perc_rewire, group = perc_rewire,
           ymin = t_60_lci, ymax = t_60_uci, fill =perc_rewire))+geom_line()+
  geom_ribbon(alpha = 0.3, color = NA)+ 
  scale_color_manual(values = c("dodgerblue",  "sienna1"))+
  scale_fill_manual(values = c("dodgerblue",  "goldenrod1"))+
  scale_y_continuous(limits = c(min_t, max_t))+
  facet_grid(cols =vars(k), rows = vars(i_k))+
  xlab(bquote(p[1]))+
  labs(subtitle = "Adopters influential for entire simulation", color = "", fill = "", y = "Time to 60% Spread")+
  coord_fixed(0.6/max_t)+
  theme_bw(base_size = basesize)

row1 <- plot_grid(plot_t1_t, 
          
          plot_t1_tl+ theme(legend.position = "bottom"), align = "hv", axis = "b", rel_widths = c(1, 1),
          labels = c("A.", "B."), label_size = 25)

row2 <- plot_grid(plot_t5_t, 
                  
                  plot_t5_tl+ theme(legend.position = "bottom"), align = "hv", axis = "b",
                  labels = c("C.", "D."), label_size = 25)

row3 <- plot_grid(plot_t0_t, 
                  
                  plot_t0_tl+ theme(legend.position = "bottom"), align = "hv", axis = "b",
                  labels = c("E.", "F."), label_size = 25)

plot_grid(row1, row2, row3, ncol = 1)
ggsave("output_plots/si_speed.pdf", width = 13, height = 16, units = "in")