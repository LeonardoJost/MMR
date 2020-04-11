### analysis incorporating time as fixed effect
#     Copyright (C) 2020  Leonardo Jost
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
dataset.rt=datasetAnalysis

#model generation

#preliminary
mBaseTime=lmer(reactionTime~deg*endTime*block*group+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTime.summary=modelSummary(mBaseTime,0)
#stepwise remove nonsignificant effects
mTime2=lmer(reactionTime~endTime*block*group+deg*block*group+deg*endTime*group+deg*endTime*block+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime2.summary=modelSummary(mTime2,0)
#split deg*block*group
mTime3=lmer(reactionTime~endTime*block*group+deg*endTime*group+deg*endTime*block+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime3.summary=modelSummary(mTime3,0)
#split endTime*block*group
mTime4=lmer(reactionTime~block*group+deg*endTime*group+deg*endTime*block+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime4.summary=modelSummary(mTime4,0)
#split deg*block*time
mTime5=lmer(reactionTime~block*group+deg*endTime*group+deg*block+endTime*block+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime5.summary=modelSummary(mTime5,0)
#split deg*block
mTime6=lmer(reactionTime~block*group+deg*endTime*group+endTime*block+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime6.summary=modelSummary(mTime6,0)
#split deg*endTime*group
mTime7=lmer(reactionTime~block*group+deg*endTime+deg*group+endTime*group+endTime*block+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime7.summary=modelSummary(mTime7,0)
#split block*group
mTime8=lmer(reactionTime~deg*endTime+deg*group+endTime*group+endTime*block+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime8.summary=modelSummary(mTime8,0)
#split deg*group
mTime9=lmer(reactionTime~deg*endTime+endTime*group+endTime*block+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime9.summary=modelSummary(mTime9,0)
#all effects significant

