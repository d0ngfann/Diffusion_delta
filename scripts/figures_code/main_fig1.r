# 메인 Figure 1 생성 스크립트 (전파 모델 유형 및 이론적 채택 곡선 시각화)
###############################
# Figure 1
###############################

library(tidyverse)
library(ggplot2)
library(cowplot)
setwd("/home/wan.a/complex_contagion_repo/")
source("scripts/analysis_plotting/helper_functions.r")

### Figure 1B
bounds_df <- get_bounds_df(c(4,8,12), c(2,3,4,6), t_val = c(1,5))%>%
  mutate(i_k = paste0("i/k = ",i/k), 
         k = factor(paste0("k = ", k), ordered = TRUE,
                    levels = c("k = 4", "k = 8", "k = 12")))

bounds_df_theory <- bounds_df %>%
  filter(t == 1) %>% 
  filter(k %in% c("k = 4", "k = 8"))%>%
  filter(i %in% c(2,3,4)) %>%
  mutate( i = paste0("i = ", i)) 


### Theoretical adoption curves plots
alpha = 0.666666667
p1 = 0.25
p2= 0.75
outlist <- list(c(0,0))
for(nc in 1:6){
  prod <- 1
  for(k in 1:nc){
    pk <- 1-((1-p1)*((1-alpha)^(k-1)))
    prod <- prod*(1-pk)
  }
  outlist[[length(outlist)+1]] <- c(k, 1-prod)
}
p_cont_cum <- as.data.frame(do.call(rbind, outlist))
colnames(p_cont_cum) <- c("k", "p_cont_cum")

outlist <- list(c(0,0))
pk_list <- c(0.25, 0.75, 0.75, 0.75, 0.75, 0.75)
for(nc in 1:6){
  prod <- 1
  for(k in 1:nc){
    pk <- pk_list[k]
    prod <- prod*(1-pk)
  }
  outlist[[length(outlist)+1]] <- c(k, 1-prod)
}
p_stoch_cum <- as.data.frame(do.call(rbind, outlist))
colnames(p_stoch_cum) <- c("k", "p_stoch_cum")


adoptdf <- expand.grid( p1 = 0.2, p2 = 0.8, k = seq(0,5, 0.01)) %>%
  mutate(p_det = ifelse(k>=2,1,0),
         p_stoch = case_when(
           k<1 ~0,
           k>=2 ~ p2,
           TRUE ~ p1),
         p_simple = ifelse(k>=1, 1, 0),
         p_sstoch = ifelse(k>=1, p1, 0),
         # cumulative prob
         p_det_cum = ifelse(k>=2,1,0),
         p_simple_cum = ifelse(k>=1,1,0),
         p_sstoch_cum = 1-((1-p1)^(k)))

adoptdf<- left_join(adoptdf, p_stoch_cum)%>% #drop_na()%>%
  pivot_longer(starts_with("p_"), names_to = "type", values_to = "pk")%>%
  mutate(p_type = ordered(ifelse(grepl("_cum", type),"Cumulative\nProb. F(c)", "Per-exposure\nProb. p(c)"),
                          levels = c( "Per-exposure\nProb. p(c)", "Cumulative\nProb. F(c)")),
         official_name = case_when(
           grepl("p_det", type)==TRUE ~ "Deterministic Complex Contagion",
           grepl("p_cont", type)==TRUE ~ "Stochastic Continuous Complex Contagion",
           grepl("p_simple", type) == TRUE ~ "Deterministic Simple Contagion",
           grepl("p_sstoch", type) == TRUE ~ "Stochastic Simple Contagion",
           grepl("p_stoch", type) == TRUE ~ "Stochastic Discrete Complex Contagion"
         ),
         det_stoch = ifelse(grepl("Deterministic", official_name)==TRUE, 
                            "Deterministic\nAdoption", "Probabilistic\nAdoption"),
         sim_com = ifelse(grepl("Simple", official_name)==TRUE, 
                          "Simple", "Complex")) %>% drop_na()


plot1 <- ggplot(adoptdf %>% filter(!(type %in% c("p_cont", "p_cont_cum"))), aes(x = k, y = pk, 
                                                                       color = official_name, linetype = sim_com))+
  geom_line(linewidth = 1)+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6))+
  scale_color_manual(values = c( "firebrick1", "green3", "gray40", "maroon1"))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  labs(x = "Num. different influential neighbors\nindividual is exposed to, c", y= "",
       color = "Contagion Type", linetype ="")+
  coord_fixed(ratio = 4)+
  theme_bw(base_size = 16)+
  facet_grid(cols = vars(det_stoch),
             rows = vars(p_type),  switch = "y")+
  theme(legend.position="bottom")+ 
  theme(legend.key.width = unit(4, "line"), strip.background.y = element_blank(),
        strip.placement.y = "outside")+ guides(color="none")+
  theme(panel.grid = element_blank(),
        strip.text.x = element_text(size =16, margin = margin(0.1,0,0.1,0, "cm")),
         strip.text.y = element_text(size =16, margin = margin( 0,0.1,0,0.1,"cm")))

### Figure 1C
plot2 <- ggplot()+
  geom_polygon(data = data.frame(x = c(0, 0, 1), y = c(0, 1, 1)),
               aes(x = x, y = y, alpha = 0.5), fill = 'gray40')+
  geom_polygon(data = data.frame(x = c(0, 1, 1), y = c(0, 1, 0)),
               aes(x = x, y = y, alpha = 0.5), fill = 'cornsilk')+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linewidth = 1, linetype = "dashed", color = "maroon1")+
  geom_segment(aes(x = 0.01, y = 0.01, xend = 0.01, yend = 1), linewidth = 1, linetype = "dashed", color = "cyan")+
  geom_segment(aes(x = 0, y = 0.99, xend = 1, yend = 0.99), linewidth = 1, linetype = "dashed", color = "gold1")+
  geom_point(aes(x = 0.99,y = 0.99), color = "green3" , size = 3, shape = 15)+
  geom_point(aes(x = 0.01,y = 0.99), color = "red" , size = 3,  shape = 15)+
  scale_x_continuous(limits = c(-0.01,1.01))+
  scale_y_continuous(limits = c(-0.01,1.01))+
  xlab(bquote(p[1]))+
  ylab(bquote(p[2]))+
  theme_bw(base_size=16) +coord_fixed()+
  theme(legend.position = "none",
        legend.justification = "left", legend.direction = "vertical")

plot_grid(plot1,plot2, align = "h", axis = "b", rel_widths = c(1, 1))
ggsave("output_plots/main_fig1.pdf", width = 14, height = 5, units = "in")


