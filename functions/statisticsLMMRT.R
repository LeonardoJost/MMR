### statistical analysis of reaction time
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

dir.create("statmodels")

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

##reaction time
#base model
mBase=lmer(reactionTime~degY*endTime*block*group+degZ*endTime+degZ*block+degY*correct_response+degY*endTime+degZ*correct_response+Gender*block+Experience*block+(deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase.summary=modelSummary(mBase)
save(mBase,mBase.summary,file="statmodels/RTmBase.RData")
##partial effects of interactions
#degY*endTime*block*group
mdegYTimeBlock=lmer(reactionTime~degY*endTime*block*group-degY:endTime:block+
                           degZ*endTime+degZ*block+degY*correct_response+degY*endTime+degZ*correct_response+Gender*block+Experience*block+(deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mdegYTimeBlock.summary=modelSummary(mdegYTimeBlock)
save(mdegYTimeBlock,mdegYTimeBlock.summary,file="statmodels/RTmdegYTimeBlock.RData")
#degY*endTime*group
mdegYTimeGroup=lmer(reactionTime~degY*endTime*block*group-degY:endTime:group+
                           degZ*endTime+degZ*block+degY*correct_response+degY*endTime+degZ*correct_response+Gender*block+Experience*block+(deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mdegYTimeGroup.summary=modelSummary(mdegYTimeGroup)
save(mdegYTimeGroup,mdegYTimeGroup.summary,file="statmodels/RTmdegYTimeGroup.RData")
#degY*block*group
mdegYBlockGroup=lmer(reactionTime~degY*endTime*block*group-degY:block:group+
                           degZ*endTime+degZ*block+degY*correct_response+degY*endTime+degZ*correct_response+Gender*block+Experience*block+(deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mdegYBlockGroup.summary=modelSummary(mdegYBlockGroup)
save(mdegYBlockGroup,mdegYBlockGroup.summary,file="statmodels/RTmdegYBlockGroup.RData")
#endTime*block*group
mTimeBlockGroup=lmer(reactionTime~degY*endTime*block*group-endTime:block:group+
                           degZ*endTime+degZ*block+degY*correct_response+degY*endTime+degZ*correct_response+Gender*block+Experience*block+(deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeBlockGroup.summary=modelSummary(mTimeBlockGroup)
save(mTimeBlockGroup,mTimeBlockGroup.summary,file="statmodels/RTmTimeBlockGroup.RData")
#nonsignificant effects






