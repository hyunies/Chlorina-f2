# HPLC Analysis ####

## Clear R brain
rm(list= ls())

## Import libraries 
library(dplyr)
library(ggplot2)
library(car)
library(lsmeans)

## set working directory
setwd("~/School Notes/R studio/HPLC_2022.csv") 

## Import HPLC Data
hplc <- read_csv("School Notes/R studio/module 1/HPLC_2022.csv")

#####################/
#####################/
## Xanthophyll 1 ####
#####################/
#####################/

## Plot 
ggplot(hplc, aes(x = treatment, y = x1, fill = genotype)) + 
  geom_boxplot() + 
  theme_bw() + 
  ylab(expression(paste("Xanthophyll 1 (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment")

## Need to transform x1
hplc <- mutate(hplc, log.x1 = log(1 + x1))

## ANOVA to compare treatments and genotypes
lm.x1 <- lm(log.x1 ~ genotype*treatment, data = hplc)

## ANOVA table
summary(aov(lm.x1))

#                       Df Sum Sq Mean Sq F value   Pr(>F)    
#   genotype             1  17.99  17.994  22.897 4.19e-06 ***
#   treatment            2   4.82   2.410   3.067  0.04963 *  
#   genotype:treatment   2  11.02   5.508   7.008  0.00125 ** 
#   Residuals          144 113.17   0.786                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Run lsmeans to compare treatments within each genotype
lsmeans(lm.x1, pairwise ~ treatment | genotype, adjust = "tukey") 

# $lsmeans
# genotype = f2:
# treatment lsmean    SE  df lower.CL upper.CL
# CNT        1.142 0.177 144    0.792     1.49
# PI         1.482 0.177 144    1.132     1.83
# REC        0.871 0.177 144    0.520     1.22
# 
# genotype = WT:
# treatment lsmean    SE  df lower.CL upper.CL
# CNT        2.364 0.177 144    2.014     2.71
# PI         1.430 0.177 144    1.080     1.78
# REC        1.779 0.177 144    1.429     2.13
# 
# Confidence level used: 0.95 
# 
# $contrasts
# genotype = f2:
# contrast  estimate    SE  df t.ratio p.value
# CNT - PI    -0.340 0.251 144  -1.355  0.3676
# CNT - REC    0.272 0.251 144   1.083  0.5260
# PI - REC     0.611 0.251 144   2.438  0.0421
# 
# genotype = WT:
# contrast  estimate    SE  df t.ratio p.value
# CNT - PI     0.934 0.251 144   3.726  0.0008
# CNT - REC    0.585 0.251 144   2.335  0.0543
# PI - REC    -0.349 0.251 144  -1.392  0.3478
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 


## Run lsmeans to compare genotypes within each treatment
lsmeans(lm.x1, pairwise ~ genotype | treatment, adjust = "tukey") 

# $lsmeans
# treatment = CNT:
#   genotype lsmean    SE  df lower.CL upper.CL
# f2        1.142 0.177 144    0.792     1.49
# WT        2.364 0.177 144    2.014     2.71
# 
# treatment = PI:
#   genotype lsmean    SE  df lower.CL upper.CL
# f2        1.482 0.177 144    1.132     1.83
# WT        1.430 0.177 144    1.080     1.78
# 
# treatment = REC:
#   genotype lsmean    SE  df lower.CL upper.CL
# f2        0.871 0.177 144    0.520     1.22
# WT        1.779 0.177 144    1.429     2.13
# 
# Confidence level used: 0.95 
# 
# $contrasts
# treatment = CNT:
#   contrast estimate    SE  df t.ratio p.value
# f2 - WT    -1.222 0.251 144  -4.873  <.0001
# 
# treatment = PI:
#   contrast estimate    SE  df t.ratio p.value
# f2 - WT     0.052 0.251 144   0.208  0.8359
# 
# treatment = REC:
#   contrast estimate    SE  df t.ratio p.value
# f2 - WT    -0.908 0.251 144  -3.622  0.0004


## Make an interaction plot

ggplot(hplc, aes(x = genotype, y = log.x1, color = treatment, group = treatment)) + 
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun = mean, geom = "line") + 
  geom_smooth(formula = y ~ x, method = "lm") + theme_classic() + ggtitle("Interaction plot:NOT FOR REPORTING")

## Make plot with indicators of significance

hplc.x1 <- hplc %>%
  group_by(genotype, treatment) %>%
  summarise(x1 = 35)

## Use * to indicate significant differences between genotypes within treatments
## Use lowercase letters to indicate significant differences between treatments in the f2 genotype
## Use UPPERCASE LETTERS to indicate significant differences between treatments in the WT genotype

## CNT  PI  REC
## *        *
## ab   a   b 
## A    B   AB

hplc.x1$sig.label <- c("* ab A", "a B", "* b AB", "", "", "")


ggplot(hplc, aes(x = treatment, y = x1, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Xanthophyll 1 (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment") + 
  geom_text(data = hplc.x1, aes(x = treatment, y = x1, label = sig.label) )

hplc_summary <- hplc %>%
  group_by(genotype, treatment) %>%
  summarise(n_genotype = length(genotype), n_treatment = length(treatment))


#####################/
#####################/
## Xanthophyll 2 ####
#####################/
#####################/

## Plot 
ggplot(hplc, aes(x = treatment, y = x2, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Xanthophyll 2 (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment")

## Log transform x2
hplc <- mutate(hplc, log.x2 = log(1 + x2))

lm.x2 <- lm(log.x2 ~ genotype*treatment, data = hplc)

summary(aov(lm.x2))

#                      Df Sum Sq Mean Sq F value  Pr(>F)    
# genotype             1 118.96  118.96 209.296 < 2e-16 ***
# treatment            2   2.11    1.05   1.855 0.16016    
# genotype:treatment   2   7.60    3.80   6.685 0.00167 ** 
#   Residuals          144  81.85    0.57                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Run lsmeans to compare treatments within each genotype
lsmeans(lm.x2, pairwise ~ treatment | genotype, adjust = "tukey")

# $lsmeans
# genotype = f2:
#   treatment lsmean    SE  df lower.CL upper.CL
# CNT        0.810 0.151 144    0.512    1.108
# PI         0.709 0.151 144    0.411    1.007
# REC        0.467 0.151 144    0.169    0.766
# 
# genotype = WT:
#   treatment lsmean    SE  df lower.CL upper.CL
# CNT        1.978 0.151 144    1.680    2.276
# PI         2.649 0.151 144    2.351    2.947
# REC        2.703 0.151 144    2.405    3.001
# 
# Confidence level used: 0.95 
# 
# $contrasts
# genotype = f2:
#   contrast  estimate    SE  df t.ratio p.value
# CNT - PI     0.101 0.213 144   0.474  0.8834
# CNT - REC    0.343 0.213 144   1.607  0.2459
# PI - REC     0.242 0.213 144   1.133  0.4956
# 
# genotype = WT:
#   contrast  estimate    SE  df t.ratio p.value
# CNT - PI    -0.671 0.213 144  -3.147  0.0057
# CNT - REC   -0.725 0.213 144  -3.400  0.0025
# PI - REC    -0.054 0.213 144  -0.253  0.9653
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

## Run lsmeans to compare genotype within treatments
lsmeans(lm.x2, pairwise ~ genotype | treatment, adjust = "tukey")

# $lsmeans
# treatment = CNT:
#   genotype lsmean    SE  df lower.CL upper.CL
# f2        0.810 0.151 144    0.512    1.108
# WT        1.978 0.151 144    1.680    2.276
# 
# treatment = PI:
#   genotype lsmean    SE  df lower.CL upper.CL
# f2        0.709 0.151 144    0.411    1.007
# WT        2.649 0.151 144    2.351    2.947
# 
# treatment = REC:
#   genotype lsmean    SE  df lower.CL upper.CL
# f2        0.467 0.151 144    0.169    0.766
# WT        2.703 0.151 144    2.405    3.001
# 
# Confidence level used: 0.95 
# 
# $contrasts
# treatment = CNT:
#   contrast estimate    SE  df t.ratio p.value
# f2 - WT     -1.17 0.213 144  -5.476  <.0001
# 
# treatment = PI:
#   contrast estimate    SE  df t.ratio p.value
# f2 - WT     -1.94 0.213 144  -9.098  <.0001
# 
# treatment = REC:
#   contrast estimate    SE  df t.ratio p.value
# f2 - WT     -2.24 0.213 144 -10.484  <.0001


## Make plot with indicators of significance

hplc.x2 <- hplc %>%
  group_by(genotype, treatment) %>%
  summarise(x2 = 50)

## Use * to indicate significant differences between genotypes within treatments
## Use lowercase letters to indicate significant differences between treatments in the f2 genotype
## Use UPPERCASE LETTERS to indicate significant differences between treatments in the WT genotype

## CNT  PI  REC
## *    *   *
## State in figure caption no significant differences for f2
## A    B   B

hplc.x2$sig.label <- c("* A", "* B", "* B", "", "", "")


ggplot(hplc, aes(x = treatment, y = x2, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Xanthophyll 2 (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment") + 
  geom_text(data = hplc.x2, aes(x = treatment, y = x2, label = sig.label) )



#####################/
#####################/
## Xanthophyll 3 ####
#####################/
#####################/

## Plot 
ggplot(hplc, aes(x = treatment, y = x3, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Xanthophyll 3 (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment")

## ANOVA
lm.x3 <- lm(x3 ~ genotype*treatment, data = hplc)

summary(aov(lm.x3))

#                      Df Sum Sq Mean Sq F value  Pr(>F)    
# genotype             1   1036    1036   2.987  0.0861 .  
# treatment            2  11380    5690  16.409 3.8e-07 ***
# genotype:treatment   2   2287    1143   3.297  0.0398 *  
# Residuals          144  49934     347                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


## Run lsmeans to compare treatments within each genotype
lsmeans(lm.x3, pairwise ~ treatment | genotype, adjust = "tukey")

# $lsmeans
# genotype = f2:
#   treatment lsmean   SE  df lower.CL upper.CL
# CNT         28.9 3.72 144    21.51     36.2
# PI          13.4 3.72 144     6.01     20.7
# REC         23.1 3.72 144    15.70     30.4
# 
# genotype = WT:
#   treatment lsmean   SE  df lower.CL upper.CL
# CNT         31.3 3.72 144    23.91     38.6
# PI          10.8 3.72 144     3.45     18.2
# REC         39.0 3.72 144    31.63     46.3
# 
# Confidence level used: 0.95 
# 
# $contrasts
# genotype = f2:
#   contrast  estimate   SE  df t.ratio p.value
# CNT - PI     15.51 5.27 144   2.944  0.0105
# CNT - REC     5.81 5.27 144   1.103  0.5136
# PI - REC     -9.70 5.27 144  -1.841  0.1600
# 
# genotype = WT:
#   contrast  estimate   SE  df t.ratio p.value
# CNT - PI     20.46 5.27 144   3.884  0.0005
# CNT - REC    -7.71 5.27 144  -1.465  0.3110
# PI - REC    -28.17 5.27 144  -5.349  <.0001
# 
# P value adjustment: tukey method for comparing a family of 3 estimates

## Run lsmeans to compare genotype within each treatment
lsmeans(lm.x3, pairwise ~ genotype | treatment, adjust = "tukey")

# $lsmeans
# treatment = CNT:
#   genotype lsmean   SE  df lower.CL upper.CL
# f2         28.9 3.72 144    21.51     36.2
# WT         31.3 3.72 144    23.91     38.6
# 
# treatment = PI:
#   genotype lsmean   SE  df lower.CL upper.CL
# f2         13.4 3.72 144     6.01     20.7
# WT         10.8 3.72 144     3.45     18.2
# 
# treatment = REC:
#   genotype lsmean   SE  df lower.CL upper.CL
# f2         23.1 3.72 144    15.70     30.4
# WT         39.0 3.72 144    31.63     46.3
# 
# Confidence level used: 0.95 
# 
# $contrasts
# treatment = CNT:
#   contrast estimate   SE  df t.ratio p.value
# f2 - WT     -2.40 5.27 144  -0.455  0.6496
# 
# treatment = PI:
#   contrast estimate   SE  df t.ratio p.value
# f2 - WT      2.55 5.27 144   0.485  0.6286
# 
# treatment = REC:
#   contrast estimate   SE  df t.ratio p.value
# f2 - WT    -15.92 5.27 144  -3.023  0.0030


## Make plot with indicators of significance

hplc.x3 <- hplc %>%
  group_by(genotype, treatment) %>%
  summarise(x3 = 120)

## Use * to indicate significant differences between genotypes within treatments
## Use lowercase letters to indicate significant differences between treatments in the f2 genotype
## Use UPPERCASE LETTERS to indicate significant differences between treatments in the WT genotype

## CNT  PI  REC
##          *
## a    b   ab
## A    B   A

hplc.x3$sig.label <- c("a A", "b B", "* ab A", "", "", "")


ggplot(hplc, aes(x = treatment, y = x3, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Xanthophyll 3 (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment") + 
  geom_text(data = hplc.x3, aes(x = treatment, y = x3, label = sig.label) )




#####################/
#####################/
## Lutein ###########
#####################/
#####################/

## Plot 
ggplot(hplc, aes(x = treatment, y = lutein, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Lutein (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment")

## ANOVA
lm.lutein <- lm(lutein ~ genotype*treatment, data = hplc)

summary(aov(lm.lutein))

#                     Df Sum Sq Mean Sq F value   Pr(>F)    
# genotype             1  57618   57618  52.830 2.13e-11 ***
# treatment            2  12459    6229   5.712   0.0041 ** 
# genotype:treatment   2   6429    3214   2.947   0.0557 .  
# Residuals          144 157052    1091                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Interaction between genotype:treatment is not significant. Remove and rerun ANOVA

lm.lutein2 <- lm(lutein ~ genotype + treatment, data = hplc)

summary(aov(lm.lutein2))

#                 Df Sum Sq Mean Sq F value  Pr(>F)    
#   genotype      1  57618   57618  51.457 3.4e-11 ***
#   treatment     2  12459    6229   5.563 0.00469 ** 
#   Residuals   146 163480    1120                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Run Tukey test to find differences
TukeyHSD(aov(lm.lutein2))

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = lm.lutein2)
# 
# $genotype
# diff      lwr      upr p adj
# WT-f2 39.19797 28.39847 49.99747     0
# 
# $treatment
#            diff        lwr       upr     p adj
# PI-CNT  22.171553   6.324731 38.018375 0.0033192
# REC-CNT 13.339174  -2.507648 29.185996 0.1174955
# REC-PI  -8.832379 -24.679201  7.014443 0.3865107

## State difference in genotype in text (body and in figure caption)
## Use uppercase lettering to indicate significant differences among treatments OVERALL - be clear in figure caption that this is the case (different from when we have an interaction).




## Make plot with indicators of significance
hplc.lutein <- hplc %>%
  group_by(genotype, treatment) %>%
  summarise(lutein = 225)

## Use uppercase to indicate overall differences among treatments

## CNT  PI  REC
## A    B   AB

hplc.lutein$sig.label <- c("A", "B", "AB", "", "", "")


ggplot(hplc, aes(x = treatment, y = lutein, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Lutein (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment") + 
  geom_text(data = hplc.lutein, aes(x = treatment, y = lutein, label = sig.label) )



#####################/
#####################/
## Chlorophyll b ####
#####################/
#####################/

## Plot 
ggplot(hplc, aes(x = treatment, y = chl.b, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Chlorophyll b (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment")

hplc <- mutate(hplc, log.chl.b = log(1 + chl.b))

## ANOVA
lm.chl.b <- lm(log.chl.b ~ genotype*treatment, data = hplc)

summary(aov(lm.chl.b))

#                     Df Sum Sq Mean Sq  F value Pr(>F)    
# genotype             1  917.2   917.2 5604.205 <2e-16 ***
# treatment            2    0.3     0.2    0.921  0.401    
# genotype:treatment   2    0.3     0.1    0.880  0.417    
# Residuals          144   23.6     0.2                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Remove non-sig interaction

lm.chl.b2 <- lm(log.chl.b ~ genotype + treatment, data = hplc)

summary(aov(lm.chl.b2))

#              Df Sum Sq Mean Sq  F value Pr(>F)    
# genotype      1  917.2   917.2 5613.462 <2e-16 ***
# treatment     2    0.3     0.2    0.922    0.4    
# Residuals   146   23.9     0.2                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#####################/
#####################/
## Chlorophyll a ####
#####################/
#####################/

## Plot 
ggplot(hplc, aes(x = treatment, y = chl.a, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste("Chlorophyll a (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment")

## ANOVA
lm.chl.a <- lm(chl.a ~ genotype*treatment, data = hplc)

summary(aov(lm.chl.a))

#                      Df  Sum Sq Mean Sq F value   Pr(>F)    
# genotype             1  837835  837835  21.113 9.38e-06 ***
# treatment            2   26187   13094   0.330   0.7195    
# genotype:treatment   2  241318  120659   3.041   0.0509 .  
# Residuals          144 5714363   39683                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Remove non-sig interaction
lm.chl.a2 <- lm(chl.a ~ genotype + treatment, data = hplc)

summary(aov(lm.chl.a2))

#               Df  Sum Sq Mean Sq F value   Pr(>F)    
# genotype      1  837835  837835  20.539 1.21e-05 ***
# treatment     2   26187   13094   0.321    0.726    
# Residuals   146 5955681   40792                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#####################/
#####################/
## B carotene    ####
#####################/
#####################/

## Plot 
ggplot(hplc, aes(x = treatment, y = b.carotene, fill = genotype)) + 
  geom_boxplot() + 
  theme_classic() + 
  ylab(expression(paste(Beta, "Carotene (", mu, "g/g)", sep = ""))) + 
  xlab("Treatment")

## ANOVA
lm.b.carotene <- lm(b.carotene ~ genotype*treatment, data = hplc)

summary(aov(lm.b.carotene))

#                     Df Sum Sq Mean Sq F value  Pr(>F)   
# genotype             1   3419    3419   8.618 0.00388 **
# treatment            2    608     304   0.766 0.46681   
# genotype:treatment   2   2258    1129   2.846 0.06137 . 
# Residuals          144  57123     397                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Remove non-significant interaction
lm.b.carotene2 <- lm(b.carotene ~ genotype + treatment, data = hplc)

summary(aov(lm.b.carotene2))

#               Df Sum Sq Mean Sq F value  Pr(>F)   
# genotype      1   3419    3419   8.405 0.00432 **
# treatment     2    608     304   0.747 0.47559   
# Residuals   146  59381     407                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
