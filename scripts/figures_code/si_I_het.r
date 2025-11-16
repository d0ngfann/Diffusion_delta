# SI Figure 9 생성 (SI 섹션 I - 이질적 개인 채택 성향 효과 분석)
###############################
# Figure S9 from SI Section I - Heterogeneous Adoption
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)
library(ggplot2)

## CREATE AND SAVE AGGREGATE FILES
df<- get_long_file(dir = "output_data/cf/",
                   match_string = "_8k_.+1T_2td_.+_0.0sd_ES")

df2<- get_long_file(dir = "output_data/cf/",
                   match_string = "_8k_.+1T_2td_.+_0.1sd_ES")

het_sum <- rbind(df, df2) %>%
  filter(k == 8)%>%
  group_by(k, perc_rewire, sig_thresh_sd, n_nodes, 
           time_of_influence, p1, p2, thrshld, seed_strat, G_name) %>%
  summarise(maxspread_m = mean(max_spread_norm),
            t_total = mean(time_to_spread, na.rm = T)) %>% 
filter(seed_strat == "seed_strat_one_nbrs")%>%
  filter(!(G_name == "MR" & perc_rewire == 1))

write.csv(het_sum, "output_data/cf_het_sum_0.0_0.1sd.csv", row.names = FALSE)

## MAKE FIGURES
het_sum <- read.csv("output_data/cf_het_sum_0.0_0.1sd.csv")

df_plot <- het_sum %>%
  mutate(G_name = ifelse(perc_rewire ==1, "ES", G_name))%>%
  pivot_wider(names_from= sig_thresh_sd, 
              values_from = c("maxspread_m", "t_total"))%>%
  mutate(mdif =maxspread_m_0 - maxspread_m_0.1,
         tdif = t_total_0 - t_total_0.1,
         G_name = factor(case_when(G_name == "ES" ~ "Rewired Ring (Random)",
                                   G_name == "WS" ~ "Ring (1D Clustered)",
                                   G_name == "MR" ~ "Moore (2D Clustered)"),
                         ordered = TRUE, 
                         levels = c("Rewired Ring (Random)", 
                                    "Ring (1D Clustered)", 
                                    "Moore (2D Clustered)")))

ggplot(data =df_plot,
       aes(x = p1, y = p2)) +
  geom_tile(aes(fill = mdif))+
  scale_fill_gradient2( limits = c(-0.27, 0.27) ,high = "green3", mid = "gray95", low = "firebrick1",
                        breaks = c(-0.2, 0, 0.2))+
  scale_y_continuous(limits = c(-0.01,1.01))+
       subtitle = "Example parameters k = 8, i = 2, T = 1",
       fill = "Difference in spread\n(Homogenous \u2212 Heterogeneous Adoption)")+ 
  facet_wrap(~G_name)+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  theme_bw(base_size=18) +coord_fixed()+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(legend.direction= "horizontal", legend.key.width = unit(1, "cm"), legend.title=element_text(size=12),
        legend.title.align=0.5, plot.title = element_text(size=16), plot.subtitle = element_text(size=14))

ggsave("output_plots/si_het.pdf", 
       width = 12, height = 6, units = "in")
