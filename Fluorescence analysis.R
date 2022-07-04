library(ggplot2)
library(dplyr)
library(car)
library(ggpattern)

rm(list = ls())

# Set working directory 
setwd("~/School Notes/R studio/fluorescence_2022.csv")

# Import data
fluor <- read.csv("~/School Notes/R studio/module 1/fluorescence_2022.csv")

# Examine an initial plot of fv.fm with treatment and genotype
ggplot(data = fluor, aes(x = treatment, y = fv.fm, fill = genotype)) + 
  geom_boxplot()

## make the anova object to test for differences in fv.fm among treatments and genotypes
fv.fm.lm <- lm(fv.fm ~ treatment*genotype, data = fluor)

## These do the same - nesting & piping
summary(aov(fv.fm.lm))

fv.fm.lm %>% aov() %>% summary()

#                      Df Sum Sq Mean Sq F value   Pr(>F)
# treatment            2 0.1032 0.05162  12.813 7.57e-06 ***
# genotype             1 0.0041 0.00411   1.021    0.314
# treatment:genotype   2 0.0083 0.00413   1.024    0.362
# Residuals          144 0.5801 0.00403
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Remove non-sig interaction
fv.fm.lm2 <- lm(fv.fm ~ treatment + genotype, data = fluor) %>% 
  aov() 

summary(fv.fm.lm2)

#                Df Sum Sq Mean Sq F value  Pr(>F)    
#   treatment     2 0.1032 0.05162   12.81 7.5e-06 ***
#   genotype      1 0.0041 0.00411    1.02   0.314    
#   Residuals   146 0.5884 0.00403                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

TukeyHSD(fv.fm.lm2)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = .)
# 
# $treatment
#               diff          lwr           upr     p adj
# PI-CNT  -0.06423692 -0.094300757 -0.0341730843 0.0000037
# REC-CNT -0.03058961 -0.060653444 -0.0005257709 0.0451306
# REC-PI   0.03364731  0.003583477  0.0637111498 0.0241170
# 
# $genotype
#             diff        lwr       upr     p adj
# WT-f2 -0.0104722 -0.03096049 0.0100161 0.3140846


##########################################/
##########################################/
## qp across treatment and genotype ####
##########################################/
##########################################/

fluor2 <- filter(fluor, qp <= 1)

## Examine an initial plot of qp with treatment and genotype
ggplot(data = fluor2, aes(x = treatment, y = qp, color = genotype)) + 
  geom_boxplot() +
  theme_bw () +
  scale_color_manual(values = c("black", "black")) + 
  labs(x = "Light Treatment",
       y = expression(paste("qP")))

## make the anova object to test for differences in qp among treatments and genotypes
qp.lm <- lm(qp ~ treatment*genotype, data = fluor2)

summary(aov(qp.lm))

#                      Df  Sum Sq  Mean Sq F value  Pr(>F)   
# treatment            2 0.00191 0.000957   2.992 0.05351 . 
# genotype             1 0.00350 0.003498  10.938 0.00121 **
# treatment:genotype   2 0.00060 0.000298   0.932 0.39631   
# Residuals          135 0.04317 0.000320                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Remove non-significant interaction

qp.lm2 <- lm(qp ~ treatment + genotype, data = fluor2)

summary(aov(qp.lm2))

#               Df  Sum Sq  Mean Sq F value Pr(>F)   
# treatment     2 0.00191 0.000957   2.995 0.0533 . 
# genotype      1 0.00350 0.003498  10.949 0.0012 **
#   Residuals   137 0.04377 0.000319                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1






##########################################/
##########################################/
## npq across treatment and genotype ####
##########################################/
##########################################/

## Examine an initial plot of npq with treatment and genotype
ggplot(data = fluor, aes(x = treatment, y = npq, color = genotype)) + 
  geom_boxplot() +
  theme_bw () +
  labs(x = "Light Treatment",
       y = expression(paste("Fv/Fm")))


## make the anova object to test for differences in npq among treatments and genotypes
npq.lm <- lm(npq ~ treatment*genotype, data = fluor)

summary(aov(npq.lm))

#                         Df Sum Sq Mean Sq F value   Pr(>F)
# treatment            2 0.1949 0.09745   8.920 0.000223 ***
# genotype             1 0.1605 0.16047  14.690 0.000189 ***
# treatment:genotype   2 0.0152 0.00760   0.696 0.500360
# Residuals          144 1.5731 0.01092
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Remove non-significant interaction

npq.lm2 <- lm(npq ~ treatment + genotype, data = fluor) %>% aov()

summary(npq.lm2)

#               Df Sum Sq Mean Sq F value   Pr(>F)    
# treatment     2 0.1949 0.09745   8.958 0.000214 ***
# genotype      1 0.1605 0.16047  14.752 0.000182 ***
# Residuals   146 1.5883 0.01088                     
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Run TukeyHSD test to find differences among treatments

TukeyHSD(npq.lm2)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = .)
# 
# $treatment
#             diff        lwr         upr     p adj
# PI-CNT  -0.08213186 -0.13152551 -0.03273821 0.0003723
# REC-CNT -0.06912826 -0.11852191 -0.01973461 0.0033085
# REC-PI   0.01300360 -0.03639005  0.06239725 0.8075477
# 
# $genotype
#           diff       lwr         upr     p adj
# WT-f2 -0.06541658 -0.099078 -0.03175515 0.0001823






