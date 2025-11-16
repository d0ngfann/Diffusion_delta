# SI Figure 4C를 위한 부트스트랩 신뢰구간 계산 (사회적 강화 임계값(i)별 분석)
###############################
## Bootstrap CI for Figure S4C in SI Section D
###############################

setwd("/home/wan.a/complex_contagion_repo/")

# Load libraries
source("scripts/analysis_plotting/si_D_summarySE.R")
require(ggplot2)
require(arrow)
require(data.table)
require(foreach)
require(doParallel)
detectCores()
registerDoParallel( detectCores() ) 


nSamples <- 10
nReps <- 1000

dat <- read_parquet("output_data/df_i_long_final.parquet")
dat <- dat[dat$seed_strat == "seed_strat_one_nbrs",]
dim(dat)
table(dat$perc_rewire)  

grid <- expand.grid( Thrshld=unique(dat$thrshld), 
                     P1=unique(dat$p1),
                     P2=unique(dat$p2) )
grid <- grid[grid$P1 <= grid$P2, ]
dim(grid)

boostrapOne <- function(thrshld, p1, p2) {
  # print(paste0("k=", degree, ", p1=", p1, ", p2=", p2))
  # Pick the right values
  # pick 10 max_spread_norm from perc_rewire == 0, and 10 from maxspread_norm from perc_rewire ==1
  rnd <- dat$max_spread_norm[dat$perc_rewire==0 & dat$thrshld==thrshld & dat$p1==p1 & dat$p2==p2]
  cls <- dat$max_spread_norm[dat$perc_rewire==1 & dat$thrshld==thrshld & dat$p1==p1 & dat$p2==p2]
  out <- replicate(nReps, {
    rrnd <- sample( rnd, nSamples)
    rcls <- sample( cls, nSamples)
    # Take the average across the 10 sample by network type. And then the difference in the sample average ("mdif")
    mean(rrnd) - mean(rcls)
  })
  data.frame( Thrshld=thrshld, P1=p1, P2=p2, Seed=1:nReps, Mdif=out)
}
boostrapOne(thrshld=4, p1=0.1, p2=0.1)
simOut <- foreach( row=1:nrow(grid) ) %dopar% { 
  boostrapOne( thrshld=grid$Thrshld[row], p1=grid$P1[row], p2=grid$P2[row] )
}
simOut <- rbindlist(simOut)
dim(simOut)

# Calculate p2/p1 
simOut$P2P1 <- simOut$P2 / simOut$P1

# Then pick 10 random rows (or a different number) where perc_rewire = 0 (clustered net) and 10 random rows where perc_rewire = 1 (random net).
# Filter to where mdif>margin, where mdif is in . return the minimum p2/p1 where mdif > margin.
# different cutoffs: c(0.001, 0.01, 0.05, 0.1)
agg <- rbind( cbind( Margin=0.1,   aggregate( P2P1 ~ Thrshld + Seed, simOut[simOut$Mdif > 0.1,  ], min) ),
              cbind( Margin=0.05,  aggregate( P2P1 ~ Thrshld + Seed, simOut[simOut$Mdif > 0.05, ], min) ),
              cbind( Margin=0.01,  aggregate( P2P1 ~ Thrshld + Seed, simOut[simOut$Mdif > 0.01, ], min) ),
              cbind( Margin=0.001, aggregate( P2P1 ~ Thrshld + Seed, simOut[simOut$Mdif > 0.001,], min) ) )

sumSE <- summarySE(agg, measurevar="P2P1", groupvars = c("Margin", "Thrshld"))
sumSE
write.csv(sumSE, "output_data/FigA4-BootstrappedCI-Thrshld.csv", row.names=FALSE)

sumSE$Margin <- factor(sumSE$Margin)

ggplot(sumSE, aes(x=Thrshld, y=P2P1, fill=Margin)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  geom_point(aes(color=Margin)) +
  geom_line(aes(color=Margin))

