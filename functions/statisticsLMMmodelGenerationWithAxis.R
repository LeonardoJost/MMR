### Model generation for statistical analysis excluding 0° cases
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
datasetForLMM0=datasetAnalysis[which(datasetAnalysis$deg!=0),]
#scaling
datasetForLMM0$deg=datasetForLMM0$deg/100
datasetForLMM0$endTime=datasetForLMM0$endTime/30 #30 minutes (time is already in minutes)
#prepare dataset
dataset.noOutlier0=datasetForLMM0[which(!datasetForLMM0$outlier),]
dataset.noOutlier0$deg=dataset.noOutlier0$deg-mean(dataset.noOutlier0$deg) #center degree
dataset.acc0=dataset.noOutlier0
dataset.rt0=dataset.noOutlier0[which(dataset.noOutlier0$typeOutlier=="hit"),]
dataset.rt0$deg=dataset.rt0$deg-mean(dataset.rt0$deg) #center degree
#normalizing time and centering degree are necessary to analyze main effects of partial interaction (block*group) when higher-
#order interactions are present (deg*block*group+time*block*group). Main effects are calculated for value 0
#0 of degree: average effect due to centering (this is "standard" main effect of removing higher order interaction)
#0 of time: difference between blocks
dataset.posttest0=dataset.rt0[which(dataset.rt0$block=="postTest"),]

#model generation
#remove random correlation
m1=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+deg*correct_response*axis+
          (deg+endTime+block+correct_response+axis+
             deg*endTime+deg*block+deg*correct_response+
             endTime*block+endTime*correct_response+
             block*correct_response+
             deg*axis+correct_response*axis||ID)+
          (deg+endTime+block+correct_response+Gender+Experience+group+axis||modelNumber),
        data=dataset.rt0,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m1),corr=FALSE)
summary(rePCA(m1))
VarCorr(m1)
#remove all effects with correlation <-.95,>.95 or NaN
m2=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+deg*correct_response*axis+
          (deg+endTime+block+correct_response+
             deg*endTime+deg*block+
             endTime*block+endTime*correct_response+
             block*correct_response||ID)+
          (deg+endTime+group||modelNumber),
        data=dataset.rt0,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m2),corr=FALSE)
summary(rePCA(m2))
VarCorr(m2)
#remove all effects with correlation <-.95,>.95 or NaN
m3=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+deg*correct_response*axis+
          (deg+endTime+block+correct_response+
             deg*endTime+deg*block+
             endTime*block||ID)+
          (deg||modelNumber),
        data=dataset.rt0,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m3),corr=FALSE)
summary(rePCA(m3))
VarCorr(m3)
#remove all effects with correlation <-.95,>.95 or NaN
m4=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+deg*correct_response*axis+
          (deg+endTime+block+
             deg*endTime+deg*block+
             endTime*block||ID)+
          (deg||modelNumber),
        data=dataset.rt0,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m4),corr=FALSE)
summary(rePCA(m4))
VarCorr(m4)
#readd random correaltion
m5=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+deg*correct_response*axis+
          (deg+endTime+block+
             deg*endTime+deg*block+
             endTime*block|ID)+
          (deg|modelNumber),
        data=dataset.rt0,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m5),corr=FALSE)
summary(rePCA(m5))
VarCorr(m5)
#model shows pca value 0, id shows very low values <.001
#reduce random slopes by lowest sd in m4
m6=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+deg*correct_response*axis+
          (deg+endTime+block+
             deg*endTime+
             endTime*block|ID)+
          (1|modelNumber),
        data=dataset.rt0,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m6),corr=FALSE)
summary(rePCA(m6))
VarCorr(m6)
#stepwise remove nonsignificant effects
m6.summary=modelSummary(m6,0)
#split deg*side*axis
m7=lmer(reactionTime~deg*endTime*block*group+
          deg*endTime*correct_response+
          Gender*block+Experience*block+
          correct_response*axis+deg*axis+
          (deg+endTime+block+deg*endTime+endTime*block|ID)+(1|modelNumber),data=dataset.rt0,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m7.summary=modelSummary(m7,0)
#all significant
m15.summary=modelSummary(m15)
plot(m15)


