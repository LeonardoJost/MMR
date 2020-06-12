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
#comparing if models were previously trained or not
mTrainedModels=lmer(reactionTime~degY*endTime*block*group+
                      degZ*block+
                      deg*correct_response+deg*endTime+
                      Gender*block+Experience*block+
                      trainedModel*block*group+
                      (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBase,mTrainedModels)
summary(mTrainedModels)
#comparing training performance
#firstDeviationTime for each participant
mfirstDeviationTime=lmer(reactionTime~degY*endTime*block*group+
                      degZ*block+
                      deg*correct_response+deg*endTime+
                      Gender*block+Experience*block+
                      firstDeviationTimeAvgByID*block*group+
                      (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBase,mfirstDeviationTime)
summary(mfirstDeviationTime)
