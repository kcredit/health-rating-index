### ----------------------------
### Environmental Health Rating Index Calculation
### Author: [Redacted]
### Description: Estimating the population-level health burden from environmental noise and air pollution exposure at the small area level in Dublin, Ireland.
### Required datasets include pre-calculated accessibility to health 'benefits' (primary healthcare and open space).
### ----------------------------

# Clear global environment
rm(list = ls())

# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(sf)
library(spdep)
library(ggplot2)
library(magrittr)
library(randomForest)
library(ranger)
library(mgcv)
library(xgboost)
library(e1071)
library(hydroGOF)
library(ALEPlot)
library(conflicted)

# Set preference for commonly masked functions
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")
conflicts_prefer(pdp::partial)
conflicts_prefer(vip::vi)

# Set working directory to project root (adjust to suit your machine or set relatively)
# setwd("~/your/project/path")

# ----------------------------
# Noise Pollution Health Impact Calculation
# ----------------------------

# noise <- read.csv("data/noise_contours.csv")

noise_d <- noise %>%
  spread(DBs, area) %>%
  select(SA_PUB2022, SHAPE_Area, T1_1AGETT, `55`, `60`, `65`, `70`, `75`)

# Coalesce values grouped by spatial unit
noise_w <- noise_d %>%
  group_by(SA_PUB2022) %>% 
  summarise(across(everything(), ~ reduce(., coalesce)), .groups = 'drop')

# Compute proportion exposed to noise bands
noise_w$D55_P <- noise_w$`55`/noise_w$SHAPE_Area
noise_w$D60_P <- noise_w$`60`/noise_w$SHAPE_Area
noise_w$D65_P <- noise_w$`65`/noise_w$SHAPE_Area
noise_w$D70_P <- noise_w$`70`/noise_w$SHAPE_Area
noise_w$D75_P <- noise_w$`75`/noise_w$SHAPE_Area

noise_w$D55_P[is.na(noise_w$D55_P)] <- 0
noise_w$D60_P[is.na(noise_w$D60_P)] <- 0
noise_w$D65_P[is.na(noise_w$D65_P)] <- 0
noise_w$D70_P[is.na(noise_w$D70_P)] <- 0
noise_w$D75_P[is.na(noise_w$D75_P)] <- 0

# Define odds ratios
noise_w$OR55 <- ifelse(noise_w$D55_P==0,1,1.000195)
noise_w$OR60 <- ifelse(noise_w$D60_P==0,1,1.01296)
noise_w$OR65 <- ifelse(noise_w$D65_P==0,1,1.061315)
noise_w$OR70 <- ifelse(noise_w$D70_P==0,1,1.15078)
noise_w$OR75 <- ifelse(noise_w$D75_P==0,1,1.286875)

# Overall regional mortality rate (I_k)
Ik <- 0.000654

# Calculate estimated risked of mortality (ER) for each band
noise_w$ER55 <- (noise_w$OR55-1)*noise_w$T1_1AGETT*Ik*noise_w$D55_P
noise_w$ER60 <- (noise_w$OR60-1)*noise_w$T1_1AGETT*Ik*noise_w$D60_P
noise_w$ER65 <- (noise_w$OR65-1)*noise_w$T1_1AGETT*Ik*noise_w$D65_P
noise_w$ER70 <- (noise_w$OR70-1)*noise_w$T1_1AGETT*Ik*noise_w$D70_P
noise_w$ER75 <- (noise_w$OR65-1)*noise_w$T1_1AGETT*Ik*noise_w$D75_P

# Total ER from noise exposure per small area
noise_w$ERTOT <- noise_w$ER55+noise_w$ER60+noise_w$ER65+noise_w$ER70+noise_w$ER75

# Overall total ER from noise exposure
sum(noise_w$ERTOT)

# ----------------------------
# Air Pollution Health Impact Calculation
# ----------------------------

# air <- read.csv("data/SA2022_Dublin_data.csv")

air <- merge(air,noise_w,on="SA_PUB2022",all.x=TRUE)

# Calculate population over 30 years
air$POP30P <- air$T1_1AGE3_1 + air$T1_1AGE3_2 + air$T1_1AGE4_1 + air$T1_1AGE4_2 + air$T1_1AGE5_1 + air$T1_1AGE5_2 + air$T1_1AGE6_1 + air$T1_1AGE6_2 + air$T1_1AGE7_1 + air$T1_1AGE7_2 + air$T1_1AGE8_1 + air$T1_1AGEG_1

# Relative risks and population attributable fractions (PAF)
air$RRPM25 <- ifelse(
  !is.na(air$PM25_ugm3_) & air$PM25_ugm3_ >= 5,
  (((air$PM25_ugm3_ - 5) / 10) * 0.08) + 1,
  1
)

air$RRNO2 <- ifelse(
  !is.na(air$NO2_ugm3_m) & air$NO2_ugm3_m >= 10,
  (((air$NO2_ugm3_m - 10) / 10) * 0.02) + 1,
  1
)

air$RRO3 <- ifelse(
  !is.na(air$O3_ugm3_me) & air$O3_ugm3_me >= 71.05,
  (((air$O3_ugm3_me - 71.05) / 10) * 0.0043) + 1,
  1
)

air$PAFPM25 <- (air$RRPM25-1)/air$RRPM25
air$PAFNO2 <- (air$RRNO2-1)/air$RRNO2
air$PAFO3 <- (air$RRO3-1)/air$RRO3

# Crude death rates
CDRAll17 <- 0.005917
CDR3017 <- 0.009650213

# Estimate attributable deaths by pollutant for each small area
air$ERPM25 <- air$PAFPM25*air$POP30P*CDR3017
air$ERNO2 <- air$PAFNO2*air$POP30P*CDR3017
air$ERO3 <- air$PAFO3*air$T1_1AGETT*CDRAll17
air$ERAQ <- air$ERPM25 + air$ERNO2 + air$ERO3

# Overall total premature deaths from pollution exposure
sum(air$ERAQ)

# Clean missing values
air <- air %>%
  mutate(across(c(ERAQ, ERTOT, GP_AI, grn_ndx_me, GS_LogitT_), ~ replace_na(., 0)))

# Standardized death rate per 100k
air$ERAQPC <- (air$ERAQ / air$T1_1AGETT) * 100000
air$ERTOTPC <- (air$ERTOT / air$T1_1AGETT) * 100000

# Z-scores for indicators
air$ERAQPCZ <- scale(air$ERAQPC)
air$ERTOTPCZ <- scale(air$ERTOTPC)
air$greenRZ <- scale(air$grn_ndx_me)
air$greensZ <- scale(air$GS_LogitT_)
air$GPZ <- scale(air$GP_AI)

# Health Risk Index (HRI): HRI_1 = 2/3 risk weighted, HRI_2 = evenly-weighted
air$HRI_1 <- (1/6)*air$GPZ + (1/12)*air$greenRZ + (1/12)*air$greensZ - (1/3)*air$ERTOTPCZ - (1/3)*air$ERAQPCZ
air$HRI_2 <- (1/4)*air$GPZ + (1/8)*air$greenRZ + (1/8)*air$greensZ - (1/4)*air$ERTOTPCZ - (1/4)*air$ERAQPCZ

# Additional variables for modeling
air$Bach_p     <- with(air, (T10_4_OD_2 + T10_4_HD_2 + T10_4_PDT + T10_4_DT) / T10_4_TT)
air$NoAuto_p   <- with(air, (T11_1_FT + T11_1_BIT + T11_1_BUT + T11_1_TDLT) / T11_1_TT)
air$BVBHth_p   <- with(air, (T12_3_BT + T12_3_VBT) / T12_3_TT)
air$log_POP    <- log(air$T1_1AGETT)
air$log_dist   <- log(air$HubDist)
air$Area       <- air$SHAPE_Area
air$POP        <- air$T1_1AGETT

air <- air %>%
  mutate(across(c(log_dist, Bach_p), ~ replace_na(., 0)))

# Final selection
air <- air %>%
  select(SA_PUB2022, RoutingKey, grn_ndx_me, greenRZ, greensZ, GPZ,
         starts_with("D"), starts_with("OR"), starts_with("ER"), starts_with("PAF"),
         ERAQ, ERAQPC, ERTOTPC, ERAQPCZ, ERTOTPCZ,
         POP30P, POP, Bach_p, NoAuto_p, BVBHth_p, log_POP, Area, log_dist,
         HRI_1, HRI_2)

# ----------------------------
# Spatial Lags and Modeling
# ----------------------------

# SA <- st_read("data/SA2022_Dublin_AllData3.shp")

coords <- st_centroid(st_geometry(SA), of_largest_polygon = TRUE)
sar.nb20 <- knearneigh(coords, k = 20)
sar.nb20 <- knn2nb(sar.nb20)
sar.wt <- nb2listw(sar.nb20, style = "W")

SA$w_HRI_1 <- lag.listw(sar.wt, SA$HRI_1)
SA$w_HRI_2 <- lag.listw(sar.wt, SA$HRI_2)
SA$w_Bach_p <- lag.listw(sar.wt, SA$Bach_p)
SA$w_NoAuto_p <- lag.listw(sar.wt, SA$NoAuto_p)
SA$w_BVBHth_p <- lag.listw(sar.wt, SA$BVBHth_p)
SA$POPD <- SA$POP / SA$Area
SA$w_POPD <- lag.listw(sar.wt, SA$POPD)

# Convert to data.frame for modeling
SA <- as.data.frame(SA)

# Create dummy variables for Dublin postal districts
SA$D01 <- ifelse(SA$RoutingKey=="DUBLIN 1",1,0)
SA$D02 <- ifelse(SA$RoutingKey=="DUBLIN 2",1,0)
SA$D03 <- ifelse(SA$RoutingKey=="DUBLIN 3",1,0)
SA$D04 <- ifelse(SA$RoutingKey=="DUBLIN 4",1,0)
SA$D05 <- ifelse(SA$RoutingKey=="DUBLIN 5",1,0)
SA$D06 <- ifelse(SA$RoutingKey=="DUBLIN 6",1,0)
SA$D6W <- ifelse(SA$RoutingKey=="DUBLIN 6W",1,0)
SA$D07 <- ifelse(SA$RoutingKey=="DUBLIN 7",1,0)
SA$D08 <- ifelse(SA$RoutingKey=="DUBLIN 8",1,0)
SA$D09 <- ifelse(SA$RoutingKey=="DUBLIN 9",1,0)
SA$D10 <- ifelse(SA$RoutingKey=="DUBLIN 10",1,0)
SA$D11 <- ifelse(SA$RoutingKey=="DUBLIN 11",1,0)
SA$D12 <- ifelse(SA$RoutingKey=="DUBLIN 12",1,0)
SA$D13 <- ifelse(SA$RoutingKey=="DUBLIN 13",1,0)
SA$D17 <- ifelse(SA$RoutingKey=="DUBLIN 17",1,0)
SA$D20 <- ifelse(SA$RoutingKey=="DUBLIN 20",1,0)

# Split into train/test
set.seed(1111)
train.index <- sample(seq_len(nrow(SA)), size = floor(0.7 * nrow(SA)))
train.SA <- SA[train.index, ]
test.SA <- SA[-train.index, ]

# Ready data for final modeling with Random Forest - to see results with XGBoost or for dependent variables HRI_1 or BVBHth_p, code will need to be adjusted
Xw <- train.SA %>%
  select(In22_ED,NoAuto_p,BVBHth_p,log_dist,POPD,SOUTH,w_HRI_2,w_NoAuto_p,w_BVBHth_p,w_POPD,
         D01,D02,D03,D04,D05,D06,D6W,D07,D08,D09,D10,D11,D12,D13,D17,D20)

my_rfw <- ranger(HRI_2 ~ In22_ED+NoAuto_p+BVBHth_p+log_dist+POPD+SOUTH+w_HRI_2+w_NoAuto_p+w_BVBHth_p+w_POPD+
                   D01+D02+D03+D04+D05+D06+D6W+D07+D08+D09+D10+D11+D12+D13+D17+D20,
                 data = train.SA, seed=1111,
                 importance = "permutation")

# Plot permutation importances
impw <- as.data.frame(rbind(importance(my_rfw))) %>%
  pivot_longer(
    cols = "In22_ED":"D20", 
    names_to = "variable",
    values_to = "importance"
  )

hrim <- ggplot(impw, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Mean Decrease in Accuracy")+
  xlab("")+
  ggtitle("HRI: Permutation Importance")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue") 

plot(hrim)

# Produce model RMSE
yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata)$predictions)

yhat(my_rfw, 
     test.SA
)

RMSE_rfw=rmse(yhat(my_rfw, test.SA),test.SA$HRI_2)
RMSE_rfw/(max(test.SA$HRI_2)-min(test.SA$HRI_2))

# Accumulated Local Effects (ALE) plots
pdjw1 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=1, K=20, NA.plot = TRUE)
pdjw2 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=2, K=20, NA.plot = TRUE)
pdjw3 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=3, K=20, NA.plot = TRUE)
pdjw4 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=4, K=20, NA.plot = TRUE)
pdjw5 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=5, K=20, NA.plot = TRUE)
pdjw6 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=6, K=20, NA.plot = TRUE)

pdjw7 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=7, K=20, NA.plot = TRUE)
pdjw8 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=8, K=20, NA.plot = TRUE)
pdjw9 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=9, K=20, NA.plot = TRUE)
pdjw10 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=10, K=20, NA.plot = TRUE)
pdjw11 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=11, K=20, NA.plot = TRUE)

pdjw12 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=12, K=20, NA.plot = TRUE)
pdjw13 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=13, K=20, NA.plot = TRUE)
pdjw14 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=14, K=20, NA.plot = TRUE)
pdjw15 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=15, K=20, NA.plot = TRUE)
pdjw16 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=16, K=20, NA.plot = TRUE)
pdjw17 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=17, K=20, NA.plot = TRUE)
pdjw18 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=18, K=20, NA.plot = TRUE)
pdjw19 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=19, K=20, NA.plot = TRUE)
pdjw20 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=20, K=20, NA.plot = TRUE)
pdjw21 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=21, K=20, NA.plot = TRUE)
pdjw22 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=22, K=20, NA.plot = TRUE)
pdjw23 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=23, K=20, NA.plot = TRUE)
pdjw24 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=24, K=20, NA.plot = TRUE)
pdjw25 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=25, K=20, NA.plot = TRUE)
pdjw26 = ALEPlot(Xw, my_rfw, pred.fun=yhat, J=26, K=20, NA.plot = TRUE)

# Produce labelled ALE plots used in manuscript
plot(pdjw8$x.values, pdjw8$f.values, type="l", cex.lab = 1.5, cex.axis = 1.5,  xlab="HRI of neighbours", ylab="ALE of y")
plot(pdjw4$x.values, pdjw4$f.values, type="l", cex.lab = 1.5, cex.axis = 1.5, xlab="Distance to nearest primary or secondary road", ylab="ALE of y")
plot(pdjw1$x.values, pdjw1$f.values, type="l", cex.lab = 1.5, cex.axis = 1.5,  xlab="Pobal Deprivation Index", ylab="ALE of y")
plot(pdjw5$x.values, pdjw5$f.values, type="l", cex.lab = 1.5, cex.axis = 1.5,  xlab="Population density", ylab="ALE of y")
plot(pdjw2$x.values, pdjw2$f.values, type="l", cex.lab = 1.5, cex.axis = 1.5,  xlab="% Non-auto commuting", ylab="ALE of y")
plot(pdjw6$x.values, pdjw6$f.values, type="l", cex.lab = 1.5, cex.axis = 1.5,  xlab="South of the Liffey", ylab="ALE of y")



