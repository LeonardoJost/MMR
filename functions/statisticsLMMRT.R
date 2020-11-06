### statistical analysis of reaction time
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
#reduce to get degY*endTime
mDegYTime1=lmer(reactionTime~endTime*block*group+degY*endTime+
                 degZ*block+
                 deg*correct_response+deg*endTime+
                 Gender*block+Experience*block+
                 (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mDegYTime1.summary=modelSummary(mDegYTime1,0)
mDegYTime2=lmer(reactionTime~block*group+endTime*group+endTime*block+
                  degY*endTime+
                  degZ*block+
                  deg*correct_response+deg*endTime+
                  Gender*block+Experience*block+
                  (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mDegYTime2.summary=modelSummary(mDegYTime2,0)
mDegYTime3=lmer(reactionTime~endTime*group+endTime*block+
                  degY*endTime+
                  degZ*block+
                  deg*correct_response+deg*endTime+
                  Gender*block+Experience*block+
                  (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mDegYTime3.summary=modelSummary(mDegYTime3)
#degree*time*side
mDegTimeSide=lmer(reactionTime~degY*endTime*block*group+
                    degZ*block+
                    deg*endTime+deg*correct_response+deg:endTime:correct_response+
                    Gender*block+Experience*block+
                    (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBase,mDegTimeSide)
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
anova(mBase,mBlockGroup)
summary(mBlockGroup)
#main effect of block
mBlock=lmer(reactionTime~degY*endTime*block*group+
              degZ*block+
              deg*correct_response+deg*endTime+
              Gender*block+Experience*block-block+
              (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBlock,mBase)
#difference between degZ*block and degY*block
mDegBlock=lmer(reactionTime~degY*endTime*block*group-
                 degY:block+
                 degZ*block+
                 deg*correct_response+deg*endTime+
                 Gender*block+Experience*block+
                 (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBase,mDegBlock)
save(mdegYBlockGroup,mTimeBlockGroup,mBlockGroup,mBlock,mDegBlock,file="statmodels/RTmBlockEffects.RData")
#main effects of gender and experience
#gender
mGender=lmer(reactionTime~degY*endTime*block*group+
              degZ*block+
              deg*correct_response+deg*endTime+
              Gender+Experience*block+
              (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mGenderNull=update(mGender, formula = ~ . -Gender)
anova(mGender,mGenderNull)
mGender.summary=modelSummary(mGender)
#experience
mExperience=lmer(reactionTime~degY*endTime*block*group+
               degZ*block+
               deg*correct_response+deg*endTime+
               Gender*block+Experience+
               (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mExperienceNull=update(mExperience, formula = ~ . -Experience)
anova(mExperience,mExperienceNull)
#average effect size of block (remove main effects of gender and experience due to nonsignificance)
mBlockE=lmer(reactionTime~degY*endTime*block*group+
               degZ*block+
               deg*correct_response+deg*endTime+
              (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
summary(mBlockE)
mBlockE.summary=modelSummary(mBlockE)



#split pre- and posttest
mPostTest1=lmer(reactionTime~degY*endTime*group+
                  deg*correct_response+deg*endTime+
                  Gender+Experience+
                  (deg+endTime|ID)+(1|modelNumber),
                data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPostTest1.summary=modelSummary(mPostTest1,0)
#remove gender
mPostTest2=lmer(reactionTime~degY*endTime*group+
                  deg*correct_response+deg*endTime+
                  Experience+
                  (deg+endTime|ID)+(1|modelNumber),
                data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPostTest2.summary=modelSummary(mPostTest2,0)
#split degY*endTime*group
mPostTest3=lmer(reactionTime~endTime*group+degY*group+degY*endTime+
                  deg*correct_response+deg*endTime+
                  Experience+
                  (deg+endTime|ID)+(1|modelNumber),
                data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPostTest3.summary=modelSummary(mPostTest3,0)
#split degY*endTime
mPostTest4=lmer(reactionTime~endTime*group+degY*group+
                  deg*correct_response+deg*endTime+
                  Experience+
                  (deg+endTime|ID)+(1|modelNumber),
                data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPostTest4.summary=modelSummary(mPostTest4,0)
#remove Experience
mPostTest=lmer(reactionTime~endTime*group+degY*group+
                  deg*correct_response+deg*endTime+
                  (deg+endTime|ID)+(1|modelNumber),
                data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPostTest.summary=modelSummary(mPostTest)
#all significant
plot(mPostTest)
save(mPostTest,mPostTest.summary,file="statmodels/RTmPostTest.RData")
#gender effects
mPostTestGender=lmer(reactionTime~endTime*group+degY*group+
                       deg*correct_response+deg*endTime+Gender+
                       (deg+endTime|ID)+(1|modelNumber),
                     data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mPostTest,mPostTestGender)
mPostTestGender.summary=modelSummary(mPostTestGender)
#experience effects
mPostTestExperience=lmer(reactionTime~endTime*group+degY*group+
                       deg*correct_response+deg*endTime+Experience+
                       (deg+endTime|ID)+(1|modelNumber),
                     data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mPostTest,mPostTestExperience)
mPostTestExperience.summary=modelSummary(mPostTestExperience)
save(mPostTestGender.summary,mPostTestExperience.summary,file="statmodels/RTmPostTestGenderExperience.RData")

#pretest
mPreTest1=lmer(reactionTime~degY*endTime*group+
                 deg*correct_response+deg*endTime+
                 Gender+Experience+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPreTest1.summary=modelSummary(mPreTest1,0)
#split deg*endTime
mPreTest2=lmer(reactionTime~degY*endTime*group+
                 deg*correct_response+
                 Gender+Experience+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPreTest2.summary=modelSummary(mPreTest2,0)
#remove gender
mPreTest3=lmer(reactionTime~degY*endTime*group+
                 deg*correct_response+
                 Experience+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPreTest3.summary=modelSummary(mPreTest3,0)
#split dgY*endTime*group
mPreTest4=lmer(reactionTime~endTime*group+degY*group+degY*endTime+
                 deg*correct_response+
                 Experience+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPreTest4.summary=modelSummary(mPreTest4,0)
#split degY*group
mPreTest5=lmer(reactionTime~endTime*group+degY*endTime+
                 deg*correct_response+
                 Experience+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPreTest5.summary=modelSummary(mPreTest5,0)
#split endTime*group
mPreTest6=lmer(reactionTime~group+degY*endTime+
                 deg*correct_response+
                 Experience+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPreTest6.summary=modelSummary(mPreTest6,0)
#remove group
mPreTest7=lmer(reactionTime~degY*endTime+
                 deg*correct_response+
                 Experience+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPreTest7.summary=modelSummary(mPreTest7,0)
#split degY*endTime
mPreTest=lmer(reactionTime~degY+endTime+
                 deg*correct_response+
                 Experience+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mPreTest.summary=modelSummary(mPreTest)
#all significant
plot(mPreTest)
save(mPreTest,mPreTest.summary,file="statmodels/RTmPreTest.RData")
#gender effects
mPreTestGender=lmer(reactionTime~degY+endTime+
                deg*correct_response+
                Experience+Gender+
                (deg+endTime|ID)+(1|modelNumber),
              data=dataset.rt.preTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mPreTest,mPreTestGender)
mPreTestGender.summary=modelSummary(mPreTestGender)
save(mPreTestGender.summary,file="statmodels/RTmPreTestGender.RData")





