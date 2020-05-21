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

##reaction time
#base model
mBase=lmer(reactionTime~degY*endTime*block*group+
           degZ*block+
           deg*correct_response+deg*endTime+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase.summary=modelSummary(mBase)
save(mBase,mBase.summary,file="statmodels/RTmBase.RData")
##partial effects of interactions
#degY:endTime:block (same model due to splitting of factor parameters in numerical interaction)
#degY*endTime*group (same model due to splitting of factor parameters in numerical interaction)
# -> split interaction
mDegYTimeBlockGroup=lmer(reactionTime~endTime*block*group+degY*block*group+degY*endTime*group+degY*endTime*block+
             degZ*block+
             deg*correct_response+deg*endTime+
             Gender*block+Experience*block+
             (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mDegYTimeBlockGroup.summary=modelSummary(mDegYTimeBlockGroup)
save(mDegYTimeBlockGroup,mDegYTimeBlockGroup.summary,file="statmodels/RTmDegYTimeBlockGroup.RData")
#degY*block*group
mdegYBlockGroup=lmer(reactionTime~degY*endTime*block*group-degY:block:group+
                       degZ*block+
                       deg*correct_response+deg*endTime+
                       Gender*block+Experience*block+
                       (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mdegYBlockGroup,mBase)
#endTime*block*group
mTimeBlockGroup=lmer(reactionTime~degY*endTime*block*group-endTime:block:group+
                       degZ*block+
                       deg*correct_response+deg*endTime+
                       Gender*block+Experience*block+
                       (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mTimeBlockGroup,mBase)
#block*group
mBlockGroup=lmer(reactionTime~degY*endTime*block*group-block:group+
                       degZ*block+
                       deg*correct_response+deg*endTime+
                       Gender*block+Experience*block+
                       (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBlockGroup,mBase)
#nonsignificant effects






