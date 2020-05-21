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

#model generation
#remove random correlation
m1=lmer(reactionTime~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
          ((degY+degZ)+endTime+block+correct_response+
             (degY+degZ)*endTime+(degY+degZ)*block+(degY+degZ)*correct_response+
             endTime*block+endTime*correct_response+
             block*correct_response||ID)+
          ((degY+degZ)+endTime+block+correct_response+Gender+Experience+group||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m1),corr=FALSE)
summary(rePCA(m1))
VarCorr(m1)
#remove all effects with correlation <-.95,>.95 or NaN
m2=lmer(reactionTime~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
          ((degY+degZ)+endTime+block+
             (degY+degZ)*endTime+(degY+degZ)*block+
             endTime*block+
             block*correct_response||ID)+
          ((degY+degZ)+endTime+group||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m2),corr=FALSE)
summary(rePCA(m2))
VarCorr(m2)
#remove all effects with correlation <-.95,>.95 or NaN or sd<avgsd/100
m3=lmer(reactionTime~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
          ((degY+degZ)+endTime+block+
             (degY+degZ)*endTime+
             endTime*block||ID)+
          ((degY+degZ)+group||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m3),corr=FALSE)
summary(rePCA(m3))
VarCorr(m3)
#remove all effects with correlation <-.95,>.95 or NaN or sd<avgsd/100
m4=lmer(reactionTime~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
          ((degY+degZ)+endTime+block+
             endTime*block||ID)+
          ((degY+degZ)||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m4),corr=FALSE)
summary(rePCA(m4))
VarCorr(m4)
#test removal of degree separation
m5=lmer(reactionTime~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block+
             endTime*block||ID)+
          (deg||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m5),corr=FALSE)
summary(rePCA(m5))
VarCorr(m5)
anova(m5,m4)
#m5 is the better model, but id shows very low variance component -> remove interaction
m6=lmer(reactionTime~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block||ID)+
          (deg||modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m6),corr=FALSE)
summary(rePCA(m6))
VarCorr(m6)
anova(m1,m2,m3,m4,m5,m6)
#readd random correlation
m7=lmer(reactionTime~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block|ID)+
          (deg|modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m7),corr=FALSE)
summary(rePCA(m7))
VarCorr(m7)
#0 component |model -> remove deg
m8=lmer(reactionTime~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
          (deg+endTime+block|ID)+
          (1|modelNumber),
        data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(m8),corr=FALSE)
summary(rePCA(m8))
VarCorr(m8)
#for information, if everything looks ok
anova(m1,m2,m3,m4,m5,m6,m7,m8)

#stepwise remove nonsignificant effects
m8.summary=modelSummary(m8,0)
m8a=lmer(reactionTime~deg*endTime*block*group+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+
           (1|modelNumber),
         data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m8,m8a)
m8b=lmer(reactionTime~(degY+degZ)*endTime*block*group+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+
           (1|modelNumber),
         data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m8,m8b)
#use m8b
m8b.summary=modelSummary(m8b,0)
#split degZ:endTime:block:group
m9=lmer(reactionTime~degY*endTime*block*group+
          degZ*block*group+degZ*endTime*group+degZ*endTime*block+
          deg*endTime*correct_response+
          Gender*block+Experience*block+
          (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m9.summary=modelSummary(m9,0)
#merge degZ and degY (merging only in the triple interactions is redundant, as degZ+degY=deg and they are included in the four-way interaction)
m9a=lmer(reactionTime~deg*endTime*block*group+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m9,m9a)
#split degZ:block:group
m10=lmer(reactionTime~degY*endTime*block*group+
           degZ*endTime*group+degZ*endTime*block+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m10a=lmer(reactionTime~deg*endTime*block*group+
            deg*endTime*correct_response+
            Gender*block+Experience*block+
            (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m10,m10a)
m10.summary=modelSummary(m10,0)
#split endTime:block:degZ 
m11=lmer(reactionTime~degY*endTime*block*group+
           degZ*endTime*group+degZ*block+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m11a=lmer(reactionTime~deg*endTime*block*group+
            deg*endTime*correct_response+
            Gender*block+Experience*block+
            (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m11,m11a)
m11.summary=modelSummary(m11,0)
#split endTime:degZ:group (keep degZ*endTime for now, in case deg*time*correct is split)
m12=lmer(reactionTime~degY*endTime*block*group+
           #degZ*endTime+
           degZ*group+
           degZ*block+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m12a=lmer(reactionTime~deg*endTime*block*group+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m12,m12a)
m12.summary=modelSummary(m12,0)
#split degZ:group
m13=lmer(reactionTime~degY*endTime*block*group+
           #degZ*endTime+
           degZ*block+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m13a=lmer(reactionTime~deg*endTime*block*group+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m13,m13a)
m13.summary=modelSummary(m13,0)
#split endTime:deg:correct_response
m14=lmer(reactionTime~degY*endTime*block*group+
           #degZ*endTime+
           degZ*block+
           endTime*correct_response+deg*correct_response+deg*endTime+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m14a=lmer(reactionTime~deg*endTime*block*group+
           endTime*correct_response+deg*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m14,m14a)
m14.summary=modelSummary(m14,0)
#split endTime:correct_response
m15=lmer(reactionTime~degY*endTime*block*group+
           #degZ*endTime+
           degZ*block+
           deg*correct_response+deg*endTime+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m15a=lmer(reactionTime~deg*endTime*block*group+
           deg*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m15,m15a)
m15b=lmer(reactionTime~degY*endTime*block*group+
           degZ*endTime+
           degZ*block+
           deg*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m15,m15b)
m15.summary=modelSummary(m15,0)
#all significant (using degZ*endTime or deg*endTime does not change the model only coefficients of degY*endTime)
plot(m15)
