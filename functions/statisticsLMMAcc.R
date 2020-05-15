### statistical analysis of accuracy
#     Copyright (C) 2019  Leonardo Jost
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


library(lme4)
library(optimx)
source("functions/helpers.R")

#load dataset
datasetForLMM=datasetAnalysis
#scaling
datasetForLMM$deg=datasetForLMM$deg/100
datasetForLMM$degY=datasetForLMM$degY/100
datasetForLMM$degZ=datasetForLMM$degZ/100
datasetForLMM$endTime=datasetForLMM$endTime/30 #30 minutes (time is already in minutes)
#prepare dataset
dataset.noOutlier=datasetForLMM[which(!datasetForLMM$outlier),]
dataset.acc=dataset.noOutlier
dataset.rt=dataset.noOutlier[which(dataset.noOutlier$typeOutlier=="hit"),]
#normalizing time is necessary to analyze main effects of partial interaction (block*group) when higher-
#order interactions are present (time*block*group). Main effects are calculated for value 0
#0 of time: difference between blocks
#degree is not centered, as 0 is a meaningful intercept and is necessary to distinguish between axes

##accuracy
#base model
aBase=glmer((type=="hit")~fixedEffects+(randomSlopes|ID)+(randomSlopes|modelNumber),family=binomial(),data=dataset.acc.axis,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
aBase.summary=modelSummary(aBase)
save(aBase,aBase.summary,file="AccModelaBase.RData")

#split interactions


#nonsignificant effects

