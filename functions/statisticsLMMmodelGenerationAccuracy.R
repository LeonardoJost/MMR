### Model generation for statistical analysis including 0° cases
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
a0=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+block+correct_response+
              deg*endTime+deg*block+deg*correct_response+
              endTime*block+endTime*correct_response+
              block*correct_response|ID)+
           (deg+endTime+block+correct_response+Gender+Experience+group|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(a0),corr=FALSE)
summary(rePCA(a0))
#remove correlation
a1=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+block+correct_response+
              deg*endTime+deg*block+deg*correct_response+
              endTime*block+endTime*correct_response+
              block*correct_response||ID)+
           (deg+endTime+block+correct_response+Gender+Experience+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a1)
anova(a1, a0)
summary(rePCA(a1))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a2=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+
              deg*endTime+
              block*correct_response||ID)+
           (deg+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a2)
summary(rePCA(a2))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a3=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+
              deg*endTime+
              block+correct_response||ID)+
           (deg||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a3)
summary(rePCA(a3))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a4=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+
              deg*endTime+correct_response||ID)+
           (deg||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a4)
summary(rePCA(a4))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a5=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+
              deg*endTime||ID)+
           (deg||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a5)
summary(rePCA(a5))
#add random slope correlation again
a6=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+
              deg*endTime|ID)+
           (deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a6)
summary(rePCA(a6))
#is ok
##remove nonsignificant fixed effects
a6.summary=modelSummary(a6,0)
#split block*experience
a7=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience+
           (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a7.summary=modelSummary(a7,0)
#remove experience
a8=glmer((type=="hit")~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+
           (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a8.summary=modelSummary(a8,0)
#split deg*time*block*group
a9=glmer((type=="hit")~endTime*block*group+deg*block*group+deg*endTime*group+deg*endTime*block+
           deg*endTime*correct_response+Gender*block+
           (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a9.summary=modelSummary(a9,0)
#split deg*time*block
a10=glmer((type=="hit")~endTime*block*group+deg*block*group+deg*endTime*group+
            deg*endTime*correct_response+Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a10.summary=modelSummary(a10,0)
#split deg*group*block
a11=glmer((type=="hit")~endTime*block*group+deg*block+deg*endTime*group+
            deg*endTime*correct_response+Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a11.summary=modelSummary(a11,0)
#split deg*time*correct_response
a12=glmer((type=="hit")~endTime*block*group+deg*block+deg*endTime*group+
            endTime*correct_response+deg*correct_response+deg*endTime+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a12.summary=modelSummary(a12,0)
#split deg*correct_response
a13=glmer((type=="hit")~endTime*block*group+deg*block+deg*endTime*group+
            endTime*correct_response+deg*endTime+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a13.summary=modelSummary(a13,0)
#split time*correct_response
a14=glmer((type=="hit")~endTime*block*group+deg*block+deg*endTime*group+
            correct_response+deg*endTime+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a14.summary=modelSummary(a14,0)
#remove correct_response
a15=glmer((type=="hit")~endTime*block*group+deg*block+deg*endTime*group+
            deg*endTime+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a15.summary=modelSummary(a15,0)
#split time*block*group
a16=glmer((type=="hit")~block*group+endTime*block+
            deg*block+deg*endTime*group+
            deg*endTime+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a16.summary=modelSummary(a16,0)
#split group*time*deg
a17=glmer((type=="hit")~block*group+endTime*block+
            deg*block+
            endTime*group+deg*group+
            deg*endTime+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a17.summary=modelSummary(a17,0)
#split group*deg
a18=glmer((type=="hit")~block*group+endTime*block+
            deg*block+
            endTime*group+
            deg*endTime+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a18.summary=modelSummary(a18,0)
#split deg*time
a19=glmer((type=="hit")~block*group+endTime*block+
            deg*block+
            endTime*group+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a19.summary=modelSummary(a19,0)
#split group*time
a20=glmer((type=="hit")~block*group+endTime*block+
            deg*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a20.summary=modelSummary(a20,0)
#split group*block
a21=glmer((type=="hit")~group+endTime*block+
            deg*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a21.summary=modelSummary(a21,0)
#remove group
a22=glmer((type=="hit")~endTime*block+
            deg*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a22.summary=modelSummary(a22,0)
#split block*gender
a23=glmer((type=="hit")~endTime*block+
            deg*block+
            Gender+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a23.summary=modelSummary(a23,0)
#remove gender
a24=glmer((type=="hit")~endTime*block+
            deg*block+
            (deg+endTime+deg*endTime|ID)+(deg|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a24.summary=modelSummary(a24,0)
#all values significant
#visual inspection of normality
plot(a24)