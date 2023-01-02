% ************************************************************************
% The code was developed at the Neurosceince Lab, Constructor University (Jacobs University)

% Contributors to MADE pipeline:
% Amir Jahanian Najafabadi (a.jahaniannajafabadi@jacobs-university.de)

% Please cite the following references for in any manuscripts produced utilizing MADE code:

% EEGLAB: A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
% analysis of single-trial EEG dynamics. Journal of Neuroscience Methods, 134, 9?21.

% This code is released under the GNU General Public License version 3.

% ************************************************************************

PALMS R Code for Statistical Analayses and Plotting  

#install.packages("readxl")
library("readxl")

my_data <- read_excel("directory.fileformat")

head(my_data)

#install.packages("tidyr")
library(tidyr)

#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2) 

#install.packages("sjPlot")
library(sjPlot) 

#install.packages("sjlabelled")
library(sjlabelled) 

#install.packages("sjmisc")
library(sjmisc) 

#install.packages("jtools")
library(jtools) 

#install.packages("ggstance")
library(ggstance) 

#install.packages("huxtable")
library(huxtable) 

#install.packages("gamlj")
#library(gamlj)

#install.packages("jmv")
library(jmv) 

#install.packages("jmvcore")
library(jmvcore)

#install.packages("car")
library(car)

#install.packages("report")
library(report)

#install.packages("ggcorrplot")
library(ggcorrplot)

#install.packages("corrplot")
library(corrplot)


setwd("~/R/PALMS")
#load libraries (make sure, packages are installed)
library(tidyverse)
## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
## v ggplot2 3.3.5     v purrr   0.3.4
## v tibble  3.1.6     v dplyr   1.0.7
## v tidyr   1.1.4     v stringr 1.4.0
## v readr   2.0.0     v forcats 0.5.1
## Warning: Paket 'tibble' wurde unter R Version 4.1.2 erstellt
## Warning: Paket 'tidyr' wurde unter R Version 4.1.2 erstellt
## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
library(ggplot2) 
library(sjPlot) 
## Learn more about sjPlot with 'browseVignettes("sjPlot")'.
library(sjlabelled) 
## 
## Attache Paket: 'sjlabelled'
## Das folgende Objekt ist maskiert 'package:forcats':
## 
##     as_factor
## Das folgende Objekt ist maskiert 'package:dplyr':
## 
##     as_label
## Das folgende Objekt ist maskiert 'package:ggplot2':
## 
##     as_label
library(sjmisc) 
## 
## Attache Paket: 'sjmisc'
## Das folgende Objekt ist maskiert 'package:purrr':
## 
##     is_empty
## Das folgende Objekt ist maskiert 'package:tidyr':
## 
##     replace_na
## Das folgende Objekt ist maskiert 'package:tibble':
## 
##     add_case
library(jtools) 
## 
## Attache Paket: 'jtools'
## Die folgenden Objekte sind maskiert von 'package:sjmisc':
## 
##     %nin%, center
library(ggstance) 
## 
## Attache Paket: 'ggstance'
## Die folgenden Objekte sind maskiert von 'package:ggplot2':
## 
##     geom_errorbarh, GeomErrorbarh
library(huxtable) 
## 
## Attache Paket: 'huxtable'
## Die folgenden Objekte sind maskiert von 'package:sjmisc':
## 
##     add_columns, add_rows
## Das folgende Objekt ist maskiert 'package:sjlabelled':
## 
##     set_label
## Das folgende Objekt ist maskiert 'package:sjPlot':
## 
##     font_size
## Das folgende Objekt ist maskiert 'package:dplyr':
## 
##     add_rownames
## Das folgende Objekt ist maskiert 'package:ggplot2':
## 
##     theme_grey
library(gamlj)
library(jmv)
library(jmvcore)
## Warning: Paket 'jmvcore' wurde unter R Version 4.1.2 erstellt
## 
## Attache Paket: 'jmvcore'
## Das folgende Objekt ist maskiert 'package:dplyr':
## 
##     select
## Die folgenden Objekte sind maskiert von 'package:base':
## 
##     endsWith, format, startsWith
library(car)
## Lade nötiges Paket: carData
## 
## Attache Paket: 'car'
## Das folgende Objekt ist maskiert 'package:dplyr':
## 
##     recode
## Das folgende Objekt ist maskiert 'package:purrr':
## 
##     some
library(report)
## Warning: Paket 'report' wurde unter R Version 4.1.2 erstellt
library(ggcorrplot)
library(corrplot)
## corrplot 0.90 loaded
library(readxl) 
#define contrasts
options(contrasts=c("contr.sum", "contr.poly"))
#loading data: 
dataset <- read_excel("Directory")
datasetV <- read_excel("Directory.xlsx")
datasetW <- read_excel("Directory.xlsx")

# factorize categorical variables
dataset$Orientation <- factor(dataset$Orientation)
dataset$Session <- factor(dataset$Session)
dataset$Lable <- factor(dataset$Lable, levels = c("Green","Blue","Yellow","Red"), labels=c("20.2","30.0","40.0","57.6"))
dataset$Feedback <- factor(dataset$Feedback)
datasetV$Session <- factor(datasetV$Session)
datasetV$Feedback <- factor(datasetV$Feedback)
datasetW$Session <- factor(datasetV$Session)
datasetW$Feedback <- factor(datasetV$Feedback)

#set theme
set_theme(
  geom.outline.color = "black", 
  geom.outline.size = 1, 
  geom.label.size = 2,
  geom.label.color = "black",
  title.color = "red", 
  title.size = 1.5, 
  axis.angle.x = 0, 
  axis.textcolor = "black", 
  base = theme_bw()
)

### JAMOVI
jmv::descriptives(
  formula = ErrorBase_mm ~ Orientation:Distance_mm,
  data = dataset,
  box = TRUE,
  dot = TRUE,
  n = FALSE,
  missing = FALSE,
  mean = FALSE,
  median = FALSE,
  sd = FALSE,
  min = FALSE,
  max = FALSE)

##
###R - lm
asym<-lm(ErrorBase_mm ~ Orientation*Distance_mm,data = dataset)
summary(report(asym))

tab_model(asym, p.style="stars", auto.label = FALSE, show.std = TRUE, show.stat = TRUE)
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.



datasetMean <- read_excel("Data/_2-3_TDJ_02May2022residuals_wide.xlsx", sheet = "Tabelle1")
ggplot(data=datasetMean, aes(x=Res_Error_H,y=Res_Error_V)) +
  geom_point(alpha=0.5) +
  labs(x= "Mean Residual Error Mediolateral", y="Mean Residual Error Proximodistal") +
  geom_abline(intercept = 0, slope = 1)


#
jmv::ttestPS(
  data = datasetMean,
  pairs = list(
    list(
      i1="Res_Error_H",
      i2="Res_Error_V")),
  wilcoxon = TRUE,
  norm = TRUE,
  meanDiff = TRUE,
  effectSize = TRUE,
  desc = TRUE)


#
jmv::ttestOneS(
  data = datasetW,
  vars = vars(RES_Error_mm_V),
  wilcoxon = TRUE,
  norm = TRUE,
  desc = TRUE)

###R - lm
TDJ<-lm(RES_Error_mm ~ Orientation*Session*Feedback,data = dataset)
summary(report(TDJ))


tab_model(TDJ, p.style="stars", auto.label = FALSE, show.std = TRUE, show.stat = TRUE)
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.


## One case with an outstanding residual error of -61.6 was excluded from further analysis.
datasetV2<-datasetV[datasetV$GLM_RES_Error_mm > (-60.),]
# Descriptives
BO_BA_desc <- descriptives(
  formula = BO + BOr + BA + BAr ~ Session:Feedback,
  data = datasetV,
  bar = FALSE)
BO_BA_desc


# bar plot with ggplot2
# convert BO/BA to long
data_long<-gather(datasetV,BOA, measure, BO:BAr, factor_key = TRUE)

# plot it
ggplot(data=data_long) +
  geom_boxplot(data=data_long, mapping=aes(x=Feedback, y=measure, fill=Session))+
  labs(y="Ratings", x="Feedback order")+
  geom_point(position=position_dodge(width=.75),aes(x=Feedback, y=measure, group=Session)) +
  facet_wrap(~BOA, nrow=2)



# MANOVA
R_manova<-manova(cbind(BO, BOr, BA, BAr) ~Session*Feedback, data=datasetV)
summary(R_manova, test="Pillai")


##
mBO <- lm(BO ~ Session*Feedback, data=datasetV)
mBOr <- lm(BOr ~ Session*Feedback, data=datasetV)
mBA <- lm(BA ~ Session*Feedback, data=datasetV)
mBAr <- lm(BAr ~ Session*Feedback, data=datasetV)
#
summary(report(mBO))

tab_model(mBO, p.style="stars", auto.label = FALSE, show.std = TRUE, show.stat = TRUE)
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.


summary(report(mBOr))

tab_model(mBOr, p.style="stars", auto.label = FALSE, show.std = TRUE, show.stat = TRUE)
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.


tab_model(mBA, p.style="stars", auto.label = FALSE, show.std = TRUE, show.stat = TRUE)
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.


summary(report(mBAr))


tab_model(mBAr, p.style="stars", auto.label = FALSE, show.std = TRUE, show.stat = TRUE)
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.
## Argument 'df_method' is deprecated. Please use 'ci_method' instead.

jmv::corrMatrix(
    data = datasetV2,
    vars = vars(BO, BOr, BA, BAr, GLM_RES_Error_mm),
    spearman = TRUE,
    flag = TRUE)


data_long<-data_long[data_long$GLM_RES_Error_mm > (-60.),]
ggplot(aes(x=GLM_RES_Error_mm,y=measure), data=data_long) +
  geom_point(alpha=0.5, group=data_long$BOA) +
  labs(x= "Residual Estimation Error", y="BO rating")+
  geom_smooth(method=lm)+
  facet_wrap(~data_long$BOA, nrow=2)
## `geom_smooth()` using formula 'y ~ x'

ggplot(data=data_long) +
  geom_point(data=data_long, mapping=aes(x=RES_Error_mm, y=measure, color=Feedback, shape=Session))+
  facet_wrap(~BOA, nrow=2)+
  geom_smooth(data=data_long, mapping=aes(x=RES_Error_mm, y=measure), method='lm')+
  xlim(-40,40)



