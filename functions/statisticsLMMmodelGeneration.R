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
datasetForLMM$endTime=datasetForLMM$endTime/30 #30 minutes (time is already in minutes)
#prepare dataset
dataset.noOutlier=datasetForLMM[which(!datasetForLMM$outlier),]
dataset.noOutlier$deg=dataset.noOutlier$deg-mean(dataset.noOutlier$deg) #center degree
dataset.acc=dataset.noOutlier
dataset.rt=dataset.noOutlier[which(dataset.noOutlier$typeOutlier=="hit"),]
dataset.rt$deg=dataset.rt$deg-mean(dataset.rt$deg) #center degree
#normalizing time and centering degree are necessary to analyze main effects of partial interaction (block*group) when higher-
#order interactions are present (deg*block*group+time*block*group). Main effects are calculated for value 0
#0 of degree: average effect due to centering (this is "standard" main effect of removing higher order interaction)
#0 of time: difference between blocks
dataset.posttest=dataset.rt[which(dataset.rt$block=="postTest"),]

#model generation
#remove random correlation
m1=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block+correct_response+
             deg*endTime+deg*block+deg*correct_response+
             endTime*block+endTime*correct_response+
             block*correct_response||ID)+
          (deg+endTime+block+correct_response+Gender+Experience+group||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m1),corr=FALSE)
summary(rePCA(m1))
VarCorr(m1)
#remove all effects with correlation <-.95,>.95 or NaN
m2=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block+
             deg*endTime+deg*block+
             endTime*block+
             block*correct_response||ID)+
          (deg||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m2),corr=FALSE)
summary(rePCA(m2))
VarCorr(m2)
#remove all effects with correlation <-.95,>.95 or NaN
m3=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block+
             deg*endTime+
             endTime*block||ID)+
          (deg||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m3),corr=FALSE)
summary(rePCA(m3))
VarCorr(m3)
anova(m2,m3)
#remove all effects with correlation <-.95,>.95 or NaN
m4=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block+
             deg*endTime||ID)+
          (deg||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m4),corr=FALSE)
summary(rePCA(m4))
VarCorr(m4)
#readd random correaltion
m5=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block+
             deg*endTime|ID)+
          (deg|modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m5),corr=FALSE)
summary(rePCA(m5))
VarCorr(m5)
#deg|modelNumber shows corr 1
m6=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block+
             deg*endTime|ID)+
          (1|modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m6),corr=FALSE)
summary(rePCA(m6))
VarCorr(m6)
anova(m1,m2,m3,m4,m5,m6)

m6a=lmer(reactionTime~deg*endTime*block*group+deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+block+
              deg*endTime|ID)+
           (0+deg|modelNumber),
         data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m6a),corr=FALSE)
summary(rePCA(m6a))
VarCorr(m6a)
anova(m5,m6,m6a)
#m6 is better than m6a. m5 is not used although it is best because pca shows component variance 0

#stepwise remove nonsignificant effects
m6.summary=modelSummary(m6,0)
#split deg*time*side
m7=lmer(reactionTime~deg*endTime*block*group+Gender*block+Experience*block+
          endTime*correct_response+deg*correct_response+
          (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m7.summary=modelSummary(m7,0)
#split time*side
m8=lmer(reactionTime~deg*endTime*block*group+Gender*block+Experience*block+
          deg*correct_response+
          (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m8.summary=modelSummary(m8,0)
#split deg*endTime*block*group
m9=lmer(reactionTime~endTime*block*group+deg*block*group+deg*endTime*group+deg*endTime*block+
          Gender*block+Experience*block+deg*correct_response+
          (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m9.summary=modelSummary(m9,0)
#split deg*block*group
m10=lmer(reactionTime~endTime*block*group+deg*endTime*group+deg*endTime*block+
           Gender*block+Experience*block+deg*correct_response+
           (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m10.summary=modelSummary(m10,0)
#split deg*block*time
m11=lmer(reactionTime~endTime*block*group+deg*endTime*group+deg*block+
           Gender*block+Experience*block+deg*correct_response+
           (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m11.summary=modelSummary(m11,0)
#split group*block*time
m12=lmer(reactionTime~endTime*block+deg*endTime*group+deg*block+
           Gender*block+Experience*block+deg*correct_response+
           (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m12.summary=modelSummary(m12,0)
#split group*block*time
m13=lmer(reactionTime~endTime*block+deg*endTime*group+
           Gender*block+Experience*block+deg*correct_response+
           (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m13.summary=modelSummary(m13,0)
#split group*deg*time
m14=lmer(reactionTime~endTime*block+endTime*group+deg*endTime+deg*group+
           Gender*block+Experience*block+deg*correct_response+
           (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m14.summary=modelSummary(m14,0)
#split group*deg
m15=lmer(reactionTime~endTime*block+endTime*group+deg*endTime+
           Gender*block+Experience*block+deg*correct_response+
           (deg+endTime+block+deg*endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m15.summary=modelSummary(m15,0)
#all significant
m15.summary=modelSummary(m15)
plot(m15)



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