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

##accuracy
a0=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           ((degY+degZ)+endTime+block+correct_response+
              (degY+degZ)*endTime+(degY+degZ)*block+(degY+degZ)*correct_response+
              endTime*block+endTime*correct_response+
              block*correct_response||ID)+
           ((degY+degZ)+endTime+block+correct_response+Gender+Experience+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(a0),corr=FALSE)
summary(rePCA(a0))
#remove correlation
a1=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           ((degY+degZ)+endTime+block+correct_response+
              (degY+degZ)*endTime+(degY+degZ)*block+(degY+degZ)*correct_response+
              endTime*block+endTime*correct_response+
              block*correct_response||ID)+
           ((degY+degZ)+endTime+block+correct_response+Gender+Experience+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a1)
anova(a1, a0)
summary(rePCA(a1))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a2=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           ((degY+degZ)+endTime+block+correct_response+
              (degY+degZ)*endTime+
              endTime*block+endTime*correct_response||ID)+
           ((degY+degZ)+endTime+Experience+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a2)
summary(rePCA(a2))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a3=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           ((degY+degZ)+endTime+block+correct_response+
              (degY+degZ)*endTime||ID)+
           ((degY+degZ)+endTime+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a3)
anova(a3, a2, a1)
summary(rePCA(a3))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a4=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           (degY+degZ+endTime+
              deg:endTime||ID)+
           (degY+degZ+endTime||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a4)
anova(a4,a3, a2, a1)
summary(rePCA(a4))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a5=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+
              deg*endTime||ID)+
           (deg+endTime||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a5)
anova(a5,a4,a3, a2, a1)
summary(rePCA(a5))
#readd random correlation
a6=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+
              deg*endTime|ID)+
           (deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a6)
anova(a6,a5,a4,a3, a2, a1)
summary(rePCA(a6))
##remove nonsignificant fixed effects
a6.summary=modelSummary(a6,0)
#split degZ*endTime*block*group
a7=glmer((type=="hit")~degY*endTime*block*group+
           degZ*block*group+degZ*endTime*group+degZ*endTime*block+
           (degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a7.summary=modelSummary(a7,0)
a7a=glmer((type=="hit")~degY*endTime*block*group+
           degZ*block*group+degZ*endTime*group+degZ*endTime*block+
           deg*endTime*correct_response+Gender*block+Experience*block+
           (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(a7,a7a)
#split block*group*degZ
a8=glmer((type=="hit")~degY*endTime*block*group+
           degZ*endTime*group+degZ*endTime*block+
           (degY+degZ)*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a8.summary=modelSummary(a8,0)
a8a=glmer((type=="hit")~degY*endTime*block*group+
           degZ*endTime*group+degZ*endTime*block+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(a8,a8a)
#split block*endTime*degZ
a9=glmer((type=="hit")~degY*endTime*block*group+
           degZ*endTime*group+degZ*block+
           (degY+degZ)*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a9.summary=modelSummary(a9,0)
a9a=glmer((type=="hit")~degY*endTime*block*group+
           degZ*endTime*group+degZ*block+
           deg*endTime*correct_response+
           Gender*block+Experience*block+
           (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(a9,a9a)
#split degY*endTime*correct_response
a10=glmer((type=="hit")~degY*endTime*block*group+
            degZ*endTime*group+degZ*block+
            degY*correct_response+
            degZ*endTime*correct_response+
            Gender*block+Experience*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a10.summary=modelSummary(a10,0)
#split Experience*block
a11=glmer((type=="hit")~degY*endTime*block*group+
            degZ*endTime*group+degZ*block+
            degY*correct_response+
            degZ*endTime*correct_response+
            Gender*block+Experience+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a11.summary=modelSummary(a11,0)
#remove Experience
a12=glmer((type=="hit")~degY*endTime*block*group+
            degZ*endTime*group+degZ*block+
            degY*correct_response+
            degZ*endTime*correct_response+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a12.summary=modelSummary(a12,0)
#split degY*correct_response
a13=glmer((type=="hit")~degY*endTime*block*group+
            degZ*endTime*group+degZ*block+
            degZ*endTime*correct_response+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a13.summary=modelSummary(a13,0)
#split degZ*endTime*correct_response
a14=glmer((type=="hit")~degY*endTime*block*group+
            degZ*endTime*group+degZ*block+
            degZ*correct_response+endTime*correct_response+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a14.summary=modelSummary(a14,0)
#split degZ*endTime*group
a15=glmer((type=="hit")~degY*endTime*block*group+
            degZ*group+degZ*endTime+
            degZ*block+
            degZ*correct_response+endTime*correct_response+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a15.summary=modelSummary(a15,0)
#split degZ*endTime
a16=glmer((type=="hit")~degY*endTime*block*group+
            degZ*group+
            degZ*block+
            degZ*correct_response+endTime*correct_response+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a16.summary=modelSummary(a16,0)
#split degZ*group
a17=glmer((type=="hit")~degY*endTime*block*group+
            degZ*block+
            degZ*correct_response+endTime*correct_response+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a17.summary=modelSummary(a17,0)
#split endTime*correct_response
a18=glmer((type=="hit")~degY*endTime*block*group+
            degZ*block+
            degZ*correct_response+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a18.summary=modelSummary(a18,0)
#split degZ*correct_response
a19=glmer((type=="hit")~degY*endTime*block*group+
            degZ*block+
            correct_response+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a19.summary=modelSummary(a19,0)
#remove correct_response
a20=glmer((type=="hit")~degY*endTime*block*group+
            degZ*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a20.summary=modelSummary(a20,0)
#split degY*endTime*block*group
a21=glmer((type=="hit")~endTime*block*group+degY*block*group+degY*endTime*group+degY*endTime*block+
            degZ*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a21.summary=modelSummary(a21,0)
#split degY*block*group
a22=glmer((type=="hit")~endTime*block*group+degY*endTime*group+degY*endTime*block+
            degZ*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a22.summary=modelSummary(a22,0)
#split degY*block*endTime
a23=glmer((type=="hit")~endTime*block*group+degY*endTime*group+degY*block+
            degZ*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a23.summary=modelSummary(a23,0)
#summarize degY*block and degZ*block
a23a=glmer((type=="hit")~endTime*block*group+degY*endTime*group+
            deg*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(a23,a23a)
#split degY*block
a24=glmer((type=="hit")~endTime*block*group+degY*endTime*group+
            degZ*block+
            Gender*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
a24.summary=modelSummary(a24,0)




#all values significant
#visual inspection of normality
plot(a24)