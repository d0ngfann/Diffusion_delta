# 시뮬레이션 데이터 분석 및 플로팅을 위한 헬퍼 함수 모음 (데이터 집계, 히트맵, 이론적 임계값 계산)
library(tidyverse)
library(ggplot2)

get_simulation_conf_file <- function(
    G_name_val, 
    k_val,
    b_val,
    i_val,
    p1_val, 
    perc_val,
    seed_strat_val , #c("random", "adj", "one_nbrs")
    sig_thresh_sd_val,
    rand_val){
  
  constant_string = "python simulation.py --trials '100' --G_name '"
  
  df <- expand.grid( G_name = G_name_val, 
                     k = k_val,
                     b = b_val,
                     i = i_val,
                     p1 = p1_val, 
                     perc = perc_val,
                     seed_strat = seed_strat_val,
                     sig_thresh_sd= sig_thresh_sd_val,
                     rand = rand_val)%>%
    filter(i<=k/2)%>%
    mutate(n = k*250,
           command = paste0(constant_string,
                            G_name,
                            "' --n '", as.character(n),
                            "' --k '", as.character(k),
                            "' --b '", as.character(b),
                            "' --thrshld '", as.character(i),
                            "' --p1 '", as.character(p1),
                            "' --perc '", as.character(perc),
                            "' --seed_strat '", seed_strat,
                            "' --sig_thresh_sd '", as.character(sig_thresh_sd),
                            "' --rand '", rand,
                            "'"))
  return(df)
}

get_long_file <- function(dir, match_string){

  files <- list.files(dir, pattern = match_string)
  
  df_big <- data.frame()
  for(f in files){
    df <- read_csv(paste0(dir,f)) %>%
      filter(p2 < 1.01)
    df_big <- rbind(df_big, df)
  }
  
  return(df_big)
}

save_long_file <- function(dir,files, summary_file){
  
  #files <- list.files(dir, pattern = match_string)
  for(f in files){
    df <- read_csv(paste0(dir,f)) %>%
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
                sum_over_60 = sum(max_spread_norm>0.6))%>%
      mutate(filename = f)
    
    if(file.exists(summary_file)==TRUE){
      existing_cf <- read_csv(summary_file)
      df <- rbind(existing_cf, df)
      write.csv(df, summary_file,row.names=FALSE)}
    else{ write.csv(df, summary_file,row.names=FALSE)}
    
  }
}

check_full_run <- function(files, summary_file){
for(f in files){
  df <- read_csv(paste0(dir,f)) %>%
    filter(p2 < 1.01)%>%
    filter(seed_strat == "seed_strat_one_nbrs")
  
  df_clean <- df %>%
    select(p1, G_name, k, n_nodes, time_of_influence, thrshld, perc_rewire, sig_thresh_sd, randomization_type) %>%
    distinct() %>%
    mutate(filename = f,nrow = nrow(df))
  
  
  if(file.exists(summary_file)==TRUE){
    existing_cf <- read_csv(summary_file)
    df_clean <- rbind(existing_cf, df_clean)
    write.csv(df_clean, summary_file,row.names=FALSE)}
  else{ write.csv(df_clean, summary_file,row.names=FALSE)}
  
}
}
  
get_dif_plot_df <- function(df_sum){
  dif_df <- df_sum %>% 
    pivot_wider(names_from= perc_rewire, 
                values_from = c("maxspread_m","t_total"))%>%
    mutate(mdif = maxspread_m_1 -maxspread_m_0, # puts random network on the right
           tdif = t_total_0 - t_total_1)
  return(dif_df)
}

get_dif_plot_spread <- function(dif_df){
  plot <- ggplot(data =dif_df,
         aes(x = p1, y = p2)) +
    geom_tile(aes(fill = mdif))+
    scale_fill_gradient2( limits = c(-1,1),high = "goldenrod1", low ="dodgerblue", mid = "gray95")+
    scale_y_continuous(limits = c(-0.01,1.01))+
    xlab(bquote(p[1]))+
    ylab(bquote(p[2]))+
    theme_bw(base_size=18) +coord_fixed()
  
  return(plot)
}

get_dif_plot_time <- function(dif_df, lim, colors){
  plot <- ggplot(data =dif_df,
                 aes(x = p1, y = p2)) +
    geom_tile(aes(fill = tdif))+
    scale_fill_gradient2( limits = lim,high = colors[1], low =colors[2], mid = "gray95")+
    scale_y_continuous(limits = c(-0.01,1.01))+
    xlab(bquote(p[1]))+
    ylab(bquote(p[2]))+
    theme_bw(base_size=18) +coord_fixed()
  
  return(plot)
}



get_rand_clust_plot <- function(df, color){
  plot <- ggplot(data =df,
                 aes(x = p1, y = p2)) +
    geom_tile(aes(fill = maxspread_m))+
    scale_fill_gradient( limits = c(0,1),high = color, low = "gray95")+
    #geom_vline(data = big_df2 %>%filter(i_k %in% c(0.25,0.5) ),aes(xintercept =p1_2,), color ="red", linetype = "dashed")+
    scale_y_continuous(limits = c(-0.01,1.01))+
    xlab(bquote(p[1]))+
    ylab(bquote(p[2]))+
    theme_bw(base_size=18) +coord_fixed()
  return(plot)
}

get_p2_T<-function(k,i, p1, t){
  p1_rev <- 1-p1
  
  sum_vec <- c()
  for(j in 1:(i-1)){
    sum_vec <- c(sum_vec, 1-(p1_rev^(t*j) ))
  }
  
  ft <- 2*(sum(sum_vec))
  top <- ft-(2*i)+k+2-(i)
  bottom <- (k-(2*i)+2)*(p1_rev^(i-1))
  
  p2<- 1-((top/bottom)^(1/((t*i)-i+1)))
  return(p2)
}

get_p2_Tv <- Vectorize(get_p2_T)

get_bounds_df <- function(k_val, i_val, t_val){
  lbound_th <- expand.grid(p1 = seq(0,1, 0.001),
                           k = k_val,
                           i = i_val, t = t_val)%>%
    filter(i<= k/2)%>%
    mutate(p2 = get_p2_Tv(k,i,p1,t),
           i_k_group = paste(i,k, sep ="_"),
           i_k = i/k,
           p1_star = 1- ((1-(1/(k-1)))^(1/t)))
  
  return(lbound_th)
}
