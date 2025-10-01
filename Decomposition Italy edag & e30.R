# ---------------------------------------------------------------------------- #
# Title:   Example of decomposition stepwise-replacement considering population composition,
#          to life expectancy  and lifespan variation changes
# Data:    Turin Longitudinal Study
# Author:  J.Daniel Zazueta-Borboa
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
#     0. Working directory, package, objects and life table function
# ---------------------------------------------------------------------------- #

# To clear everything in R, before start the analysis
rm(list = ls())

# Working directory from safe environment
setwd("/Users/DanielZ/Desktop/Third paper/R Code 1x1 - MPIDR")

# Packages
library(tidyverse) 
library(data.table) 
library(foreign)
library(magrittr)
library(patchwork)
library(hrbrthemes)
library(here)
library(ggdark)
library(segmented)
library(lmtest)
library(scales)
library(ggpubr)
library(DemoDecomp)

# Important labels by country and SES group
Sex <- c("Males", "Females")
SES <- c("Low","Middle","High", "Overall")
Census_Italy <- c(1971,1976,1981,1986,
                  1991,1996,2001,2006,
                  2011,2016)
Color_SES <- c("#D9D9D9", "#737373", "#252525")
Color_SES_Total <- c("#D9D9D9", "#737373", "#252525", "Black")

SES_Edu <-c("1"="Low",
            "2"="Medium",
            "3"="HIgh")

Labels.age  <- c('30-34','35-39','40-44','45-49',
                 '50-54','55-59','60-64','65-69',
                 "70-74","75-79","80-84","85-89",
                 "90-94","95-99")


# ---------------
#   Life table by single age groups
# ---------------

SLT_1 <- function (nmx, age, sex){
  n <- 1
  x <- age
  #n <- c(diff(x), NA)
  nax <- rep(0, length(x))
  m <- length(x)
  for(i in 1:length(nax)){
    nax[i] <- 0.5*n
    nax[m] <- 0 # Last nax with 0
  }# formulta to compute last nax 
  
  nqx <- (n * nmx)/(1 + (n - nax) * nmx)
  nqx <- c(nqx[-(length(nqx))], 1)
  for (i in 1:length(nqx)) {
    if (nqx[i] > 1) 
      nqx[i] <- 1
  }
  npx <- 1 - nqx
  l0 = 100000 #standard
  lx <- round(cumprod(c(l0, npx)))
  ndx <- -diff(lx)
  lxpn <- lx[-1]
  nLx <- n * lxpn + ndx * nax
  nLx[m] <- ndx[m]/nmx[m] # We compute the last nLx
  
  if (nmx[m]>=0) {
    nLx[m] <- 0
  } else {
    nLx[m] <- ndx[m]/nmx[m]
  }
  Tx <- rev(cumsum(rev(nLx)))
  lx <- lx[1:length(age)]
  ex <- Tx/lx
  e30 <- ex[1]
  # Adjusting the ax value fo the last age groups
  nax[m] <- ex[m] # we replace the nax in the last age group with the ex
  
  v        <- ((nax*c(ex[-1L],0) + (1-nax)*ex))
  v[length(ex)] <- ex[length(ex)]
  v <- (ndx*v)/l0
  e.dagger <- rev(cumsum(rev(v)))
  e.dagger[1]
  
  # SD
  vx <- sum(ndx*(((age-30)+nax)-e30)^2)/l0
  sd <- sqrt(vx)
  sd <- sd
  
  
  # Adjusting the ax value for the last age groups
  nax[m] <- ex[m] # we replace the nax in the last age group with the ex
  
  lt <- data.frame(age = age, 
                   n,
                   nmx = round(nmx, 4),
                   nqx = round(nqx, 4),
                   nax = round(nax, 4),
                   npx = round(npx,4),
                   lx, ndx,
                   lxpn = round(lxpn,4),
                   Lx = round(nLx,3),
                   Tx = round(Tx,3),
                   ex = round(ex,3),
                   e.dagger = round(e.dagger,3),
                   sd = round(sd,3))
  
  return(lt)
}



LTweighted <- function(mxwx =  mxwx){
  
  lengthvec <- length(mxwx)
  mx <- mxwx[1:(lengthvec / 2)]
  wx <- mxwx[(lengthvec / 2 + 1):lengthvec]
  
  
  Mx <- matrix(mx, nrow = 66) # number of age groups
  wx <- matrix(wx, nrow = 66) # number of age groups
  
  
  weighted_mx <- rowSums(Mx*wx)
  
  LT <- SLT_1(nmx = weighted_mx, age = c(seq(30, 95, 1)))
  e0 <- LT$ex[1]
  
  
  return(e0)
}



LTweighted_edag <- function(mxwx =  mxwx){
  
  lengthvec <- length(mxwx)
  mx <- mxwx[1:(lengthvec / 2)]
  wx <- mxwx[(lengthvec / 2 + 1):lengthvec]
  
  
  Mx <- matrix(mx, nrow = 66) # number of age groups
  wx <- matrix(wx, nrow = 66) # number of age groups
  
  
  weighted_mx <- rowSums(Mx*wx)
  
  LT <- SLT_1(nmx = weighted_mx, age = c(seq(30, 95, 1)))
  edag <- LT$e.dagger[1]
  
  
  return(edag)
}

# ---------------------------------------------------------------------------- #
#     1. Read and prepare data
# ---------------------------------------------------------------------------- #


Italy_Data_example <- read.csv("Data/Data_example.csv", header = T)

# Compute the national weighted mx
Data_for_LT <- Italy_Data_example %>% 
  mutate(Mx = (mx_1*Share_1) + (mx_2*Share_2) + (mx_3*Share_3)) %>% 
  data.table()

# ---------------------------------------------------------------------------- #
#     2. Compute Life tables 
# ---------------------------------------------------------------------------- #

# Life table National weighted
Italy_Lifetables_Total_weighted <- Data_for_LT[,cbind(SLT_1(nmx = Mx,
                                                                 sex = sex,
                                                                 age = age)), by = list(year, sex)]

Italy_e30_edagger_Total_weighted <- Italy_Lifetables_Total_weighted  %>% 
  filter(age==30) %>% 
  dplyr::select(year, sex, ex, e.dagger)


# ---------------------------------------------------------------------------- #
#  3. Decomposition analysis 
# ---------------------------------------------------------------------------- #

# 1. We need to create one single dataset where mx goes first and then wx (share)
#.   for that reason we aggregate the data

# -----------------------------------------
#  3.0 Prepare data for decomposition
# -----------------------------------------

Italy_mx <- Data_for_LT %>% 
  dplyr::select(-c(Mx, Share_1, Share_2, Share_3)) %>% 
  pivot_longer(!c(year, sex, age),
               names_to = "SES",
               values_to = "mxwx") %>% 
  mutate(Edu=case_when(SES=="mx_1" ~ 1,
                       SES=="mx_2" ~ 2,
                       SES=="mx_3" ~ 3)) %>% 
  dplyr::select(-SES) %>% 
  arrange(sex, year, Edu,age) # This decomposition is sensitive to the order, for that reason we sort the data by
# sex, period, education and age

# ---------------------------------------
#   Share
# ---------------------------------------

Italy_share <- Data_for_LT %>% 
  dplyr::select(-c(Mx, mx_1, mx_2, mx_3)) %>% 
  pivot_longer(!c(year, sex, age),
               names_to = "SES",
               values_to = "mxwx") %>% 
  mutate(Edu=case_when(SES=="Share_1" ~ 1,
                       SES=="Share_2" ~ 2,
                       SES=="Share_3" ~ 3)) %>% 
  dplyr::select(-SES) %>% 
  arrange(sex, year, Edu, age) %>% 
  arrange(sex, year, Edu,age) # This decomposition is sensitive to the order, for that reason we sort the data by
# sex, period, education and age


# ------
# We combine both datasets
# -----

Data_for_decomposition <- rbind(Italy_mx, Italy_share)

# ---------------------------------------
#   Females
# ---------------------------------------

# 1975
Italy_F_mxwx_1975 <- Data_for_decomposition %>% 
  filter(sex==2 & year==1975) %>% 
  dplyr::select(mxwx) %>% 
  as.matrix()


# 2015
Italy_F_mxwx_2015 <- Data_for_decomposition %>% 
  filter(sex==2 & year==2015) %>% 
  dplyr::select(mxwx) %>% 
  as.matrix()


# ----------------------------------------------------------------------------- #
#  Females - Life expectancy decomposition  1975-2015
# ----------------------------------------------------------------------------- #


Ita_F_75_15 <- stepwise_replacement(func = LTweighted,
                                    pars1 = Italy_F_mxwx_1975[,1],
                                    pars2 = Italy_F_mxwx_2015[,1])

sum(Ita_F_75_15)

Decomp_Italy_F_75_15 <- matrix(Ita_F_75_15, nrow = 66)
p_effects <- rowSums(Decomp_Italy_F_75_15[,4:6])
tot_effects <- rowSums(Decomp_Italy_F_75_15[,1:6])
Decomp_Italy_F_75_15 <- cbind(Decomp_Italy_F_75_15, p_effects, tot_effects)
Decomp_Italy_F_75_15 <- cbind(seq(30,95,1), Decomp_Italy_F_75_15)
colnames(Decomp_Italy_F_75_15) <- c("Age",
                                    "M1", "M2","M3", 
                                    "P1", "P2","P3", 
                                    "Ptotal", "TotEff")
Decomp_Italy_F_75_15 <- as.data.frame(Decomp_Italy_F_75_15)
Decomp_Italy_F_75_15$Sex <- "Females"
Decomp_Italy_F_75_15$Country <- "Italy"
Decomp_Italy_F_75_15$Period <- "1975-2015"


# ----------------------------------------------------------------------------- #
#  Females - Lifespan inequality decomposition 1975-1985
# ----------------------------------------------------------------------------- #

# ---------------
#   1975-2015
# ---------------

Ita_F_75_15 <- stepwise_replacement(func = LTweighted_edag,
                                    pars1 = Italy_F_mxwx_1975[,1],
                                    pars2 = Italy_F_mxwx_2015[,1])

sum(Ita_F_75_15)

Decomp_Italy_F_75_15 <- matrix(Ita_F_75_15, nrow = 66)
p_effects <- rowSums(Decomp_Italy_F_75_15[,4:6])
tot_effects <- rowSums(Decomp_Italy_F_75_15[,1:6])
Decomp_Italy_F_75_15 <- cbind(Decomp_Italy_F_75_15, p_effects, tot_effects)
Decomp_Italy_F_75_15 <- cbind(seq(30,95,1), Decomp_Italy_F_75_15)
colnames(Decomp_Italy_F_75_15) <- c("Age",
                                    "M1", "M2","M3", 
                                    "P1", "P2","P3", 
                                    "Ptotal", "TotEff")
Decomp_Italy_F_75_15 <- as.data.frame(Decomp_Italy_F_75_15)
Decomp_Italy_F_75_15$Sex <- "Females"
Decomp_Italy_F_75_15$Country <- "Italy"
Decomp_Italy_F_75_15$Period <- "1975-2015"


