# 시뮬레이션 작업 배열 실행을 위한 배치 설정 파일(.conf) 생성 스크립트
###############################
# Script for generating configuration files to run simulations in job array
###############################

setwd("/home/wan.a/complex_contagion_repo/")

source("scripts/analysis_plotting/helper_functions.r")
library(tidyverse)

## Main Fig 2 + SI Section B, E, F
# k = 8, b =1 already covered in seeding
df <-get_simulation_conf_file(
    G_name_val = "WS", 
    k_val = c(4,8,12),
    b_val = c(1, 5, 0),
    i_val = c(2,3,4,6),
    p1_val = seq(0, 1, 0.02), 
    perc_val = c(0,1),
    seed_strat_val = "one_nbrs",
    sig_thresh_sd_val= 0,
    rand_val = "ES") %>% 
    mutate(i_k = i/k)%>%
    filter(i_k %in% c(0.25, 0.5))
write.table(df %>% filter(perc == 0) %>%select(command), 
            "scripts/batch_main_perc0.conf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(df %>% filter(perc == 1) %>%select(command), 
            "scripts/batch_main_perc1.conf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

## Main Fig 3A + SI Section C & D 
# By Degree (k), i = 2, T = 1
# k = 4,8 already covered in main fig 2
df <-get_simulation_conf_file(
    G_name_val = "WS", 
    k_val = c(6,10, 12, 14,16,18, 20),
    b_val = 1,
    i_val = 2,
    p1_val = seq(0, 1, 0.02), 
    perc_val = c(0,1),
    seed_strat_val = "one_nbrs",
    sig_thresh_sd_val= 0,
    rand_val = "ES")
write.table(df %>%select(command), 
            "scripts/batch_k.conf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

## Main Fig 3B + SI Section C & D
# By time of influence (T; b_val in this code), k = 8, i = 2
# b = 1,5 already covered in main fig 2
df <- get_simulation_conf_file(
    G_name_val = "WS", 
    k_val = 8,
    b_val = c(2,3,4,6,7,8,9,10),
    i_val = 2,
    p1_val = seq(0, 1, 0.02), 
    perc_val = c(0,1),
    seed_strat_val = "one_nbrs",
    sig_thresh_sd_val= 0,
    rand_val = "ES")
write.table(df  %>%select(command), 
            "scripts/batch_T.conf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

## Main Fig 3B + SI Section C & D
# By social reinforcement threshold (i), k = 20, T = 1
# k=20, i= 2 already done in Main Fig 3A
df <- get_simulation_conf_file(
    G_name_val = "WS", 
    k_val = 20,
    b_val = 1,
    i_val = 3:10,
    p1_val = seq(0, 1, 0.02), 
    perc_val =c(0, 1),
    seed_strat_val = "one_nbrs",
    sig_thresh_sd_val= 0,
    rand_val = "ES")
write.table(df  %>%select(command), 
            "scripts/batch_i.conf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
            
## SI Section G (Seeding) and H (2D Lattice)
# "one_nbrs" seeding on WS ring lattices already done in main fig 2
# "adj" seeding only needed for WS ring lattice
# perc == 0 only for Moore
df <- get_simulation_conf_file(
    G_name_val = c("WS", "MR"), 
    k_val = 8,
    b_val = 1,
    i_val = c(2,4),
    p1_val = seq(0, 1, 0.02), 
    perc_val = c(0,1),
    seed_strat_val = c("random", "adj", "one_nbrs"),
    sig_thresh_sd_val= 0,
    rand_val = "ES")%>%
  filter(!(G_name == "WS" & seed_strat == "one_nbrs"))%>%
  filter(!(G_name == "MR" & perc == 1 ))%>%
  filter(!(G_name == "MR" & seed_strat == "adj"))

write.table(df  %>%select(command), 
              "scripts/batch_seeding.conf", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
    
## SI Section I (Heterogeneous Adoption)
#sig_thresh_sd = 0 already done in main fig 2
df <- get_simulation_conf_file(
  G_name_val = c("WS", "MR"), 
  k_val = 8,
  b_val = 1,
  i_val = 2,
  p1_val = seq(0, 1, 0.02), 
  perc_val = c(0, 1),
  seed_strat_val = "one_nbrs",
  sig_thresh_sd_val=  0.1,
  rand_val = "ES")%>%
  filter(!(perc == 1 & G_name == "MR"))
write.table(df  %>%select(command), 
            "scripts/batch_het.conf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
  
## SI Section J (Rewiring)
constant_string = "python si_rewiring.py --trials '100' --G_name '"
  
df <- expand.grid( G_name = c("WS", "MR"), 
                     k = 8,
                     b = c(0,1),
                     thrshld = 2,
                     p1 = c(0, 0.001, 0.01, 0.1, 0.2),
                     p2 = c(0.25, 0.5, 0.75, 1),
                     sig_thresh_sd= 0)%>%
    mutate(n = k*250,
           command = paste0(constant_string,
                            G_name,
                            "' --n '", as.character(n),
                            "' --k '", as.character(k),
                            "' --b '", as.character(b),
                            "' --thrshld '", as.character(thrshld),
                            "' --p1 '", as.character(p1),
                            "' --p2 '", as.character(p2),
                            "' --sig_thresh_sd '", as.character(sig_thresh_sd),
                            "'")) %>%
    filter(p1<= p2)
  write.table(df  %>%select(command), 
              "scripts/batch_rewiring.conf", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)


## SI Section K (Centoal and Macy Replication)
constant_string = "python si_macy_centola_rep.py --trials '100' --G_name '"
  
df <- expand.grid( G_name = c("WS", "MR"), 
                     k = 8,
                     b = 0,
                     perc = c(seq(0, 0.01, 0.001), seq(0.02, 0.5, 0.01),seq(0.52,1,0.02)),
                     m = 1:10,
                     seed_strat ="one_nbrs",
                     sig_thresh_sd= 0,
                     rand = "ES")%>%
    mutate(n = k*250,
           command = paste0(constant_string,
                            G_name,
                            "' --n '", as.character(n),
                            "' --k '", as.character(k),
                            "' --b '", as.character(b),
                            "' --perc '", as.character(perc),
                            "' --m '", as.character(m),
                            "' --seed_strat '", seed_strat,
                            "' --sig_thresh_sd '", as.character(sig_thresh_sd),
                            "' --rand '", rand,
                            "'"))
write.table(df  %>%select(command), 
            "scripts/batch_mc07.conf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##  SI Section L (Centola 2010 Replication)
constant_string = "python si_centola_2010_rep.py --r_start '0' --trials '100' --k '"

df <- expand.grid( p1 = seq(0, 1, 0.02), i = 2, b = c(1,2), k = c(6,8), sig_thresh_sd = c(0, 0.1))%>%
  mutate(
    command = paste0(constant_string,
                     as.character(k),
                     "' --b '", as.character(b),
                     "' --p1 '", as.character(p1),
                     "' --sig_thresh_sd '", as.character(sig_thresh_sd),
                     "'"))
write.table(df %>%select(command), 
            "scripts/batch_centola_2010.conf", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
