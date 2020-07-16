### exploratory statistical analysis of reaction time
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

###base models

mBase=lmer(reactionTime~degY*endTime*block*group+
             degZ*block+
             deg*correct_response+deg*endTime+
             Gender*block+Experience*block+
             (deg+endTime+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#posttest base model
mPostTest=lmer(reactionTime~endTime*group+degY*group+
                 deg*correct_response+deg*endTime+
                 (deg+endTime|ID)+(1|modelNumber),
               data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))

###comparing if models were previously trained or not in posttest performance

mTrainedModels=lmer(reactionTime~endTime*group+degY*group+
                      deg*correct_response+deg*endTime+
                      trainedModel*group+
                      (deg+endTime|ID)+(1|modelNumber),
                    data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainedModels2=lmer(reactionTime~endTime*group+degY*group+
                       deg*correct_response+deg*endTime+
                       trainedModel+
                       (deg+endTime|ID)+(1|modelNumber),
                     data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#combining models prePost and traiPost as trained models
mTrainedModels3=lmer(reactionTime~endTime*group+degY*group+
                       deg*correct_response+deg*endTime+
                       trainedModelPost*group+
                       (deg+endTime|ID)+(1|modelNumber),
                     data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainedModels4=lmer(reactionTime~endTime*group+degY*group+
                       deg*correct_response+deg*endTime+
                       trainedModelPost+
                       (deg+endTime|ID)+(1|modelNumber),
                     data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mPostTest,mTrainedModels,mTrainedModels2,mTrainedModels3,mTrainedModels4)
mTrainedModels4.summary=modelSummary(mTrainedModels4)
save(mTrainedModels4,mTrainedModels4.summary,file="statmodels/RTmTrainedModels.RData")

###comparing training performance

#all possible training effects
#careful: training parameters are related
#-> leave out number of training trials, possible comparison at end
mTrainingEffects=lmer(reactionTime~endTime*group+degY*group+
                        deg*correct_response+deg*endTime+
                        firstDeviationTimeAvgByID*group+
                        rotationSpeedAbsAvgByID*group+
                        comparisonTimeAvgByID*group+
                        shortDirectionPropByID*group+
                        numberOfSwitchesByID*group+
                        #numberOfTrainingTrialsByID*group+
                        numberOfPretestTrialsByID*group+
                        (deg+endTime|ID)+(1|modelNumber),
                      data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects.summary=modelSummary(mTrainingEffects,0)
#split firstDeviationTimeAvgByID*group
mTrainingEffects2=lmer(reactionTime~endTime*group+degY*group+
                         deg*correct_response+deg*endTime+
                         firstDeviationTimeAvgByID+
                         rotationSpeedAbsAvgByID*group+
                         comparisonTimeAvgByID*group+
                         shortDirectionPropByID*group+
                         numberOfSwitchesByID*group+
                         #numberOfTrainingTrialsByID*group+
                         numberOfPretestTrialsByID*group+
                         (deg+endTime|ID)+(1|modelNumber),
                       data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects2.summary=modelSummary(mTrainingEffects2,0)
#split rotationSpeedAbsAvgByID*group
mTrainingEffects3=lmer(reactionTime~endTime*group+degY*group+
                         deg*correct_response+deg*endTime+
                         firstDeviationTimeAvgByID+
                         rotationSpeedAbsAvgByID+
                         comparisonTimeAvgByID*group+
                         shortDirectionPropByID*group+
                         numberOfSwitchesByID*group+
                         #numberOfTrainingTrialsByID*group+
                         numberOfPretestTrialsByID*group+
                         (deg+endTime|ID)+(1|modelNumber),
                       data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects3.summary=modelSummary(mTrainingEffects3,0)
#remove rotationSpeedAbsAvgByID
mTrainingEffects4=lmer(reactionTime~endTime*group+degY*group+
                         deg*correct_response+deg*endTime+
                         firstDeviationTimeAvgByID+
                         comparisonTimeAvgByID*group+
                         shortDirectionPropByID*group+
                         numberOfSwitchesByID*group+
                         #numberOfTrainingTrialsByID*group+
                         numberOfPretestTrialsByID*group+
                         (deg+endTime|ID)+(1|modelNumber),
                       data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects4.summary=modelSummary(mTrainingEffects4,0)
#split numberOfSwitchesByID*group
mTrainingEffects5=lmer(reactionTime~endTime*group+degY*group+
                         deg*correct_response+deg*endTime+
                         firstDeviationTimeAvgByID+
                         comparisonTimeAvgByID*group+
                         shortDirectionPropByID*group+
                         numberOfSwitchesByID+
                         #numberOfTrainingTrialsByID*group+
                         numberOfPretestTrialsByID*group+
                         (deg+endTime|ID)+(1|modelNumber),
                       data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects5.summary=modelSummary(mTrainingEffects5,0)
#remove numberOfSwitchesByID
mTrainingEffects6=lmer(reactionTime~endTime*group+degY*group+
                         deg*correct_response+deg*endTime+
                         firstDeviationTimeAvgByID+
                         comparisonTimeAvgByID*group+
                         shortDirectionPropByID*group+
                         #numberOfTrainingTrialsByID*group+
                         numberOfPretestTrialsByID*group+
                         (deg+endTime|ID)+(1|modelNumber),
                       data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects6.summary=modelSummary(mTrainingEffects6,0)
#remove firstDeviationTimeAvgByID
mTrainingEffects7=lmer(reactionTime~endTime*group+degY*group+
                         deg*correct_response+deg*endTime+
                         comparisonTimeAvgByID*group+
                         shortDirectionPropByID*group+
                         #numberOfTrainingTrialsByID*group+
                         numberOfPretestTrialsByID*group+
                         (deg+endTime|ID)+(1|modelNumber),
                       data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects7.summary=modelSummary(mTrainingEffects7,0)
#split numberOfPretestTrialsByID*group
mTrainingEffects8=lmer(reactionTime~endTime*group+degY*group+
                         deg*correct_response+deg*endTime+
                         comparisonTimeAvgByID*group+
                         shortDirectionPropByID*group+
                         #numberOfTrainingTrialsByID*group+
                         numberOfPretestTrialsByID+
                         (deg+endTime|ID)+(1|modelNumber),
                       data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects8.summary=modelSummary(mTrainingEffects8,0)
#split comparisonTimeAvgByID*group+
mTrainingEffects9=lmer(reactionTime~endTime*group+degY*group+
                         deg*correct_response+deg*endTime+
                         comparisonTimeAvgByID+
                         shortDirectionPropByID*group+
                         #numberOfTrainingTrialsByID*group+
                         numberOfPretestTrialsByID+
                         (deg+endTime|ID)+(1|modelNumber),
                       data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects9.summary=modelSummary(mTrainingEffects9,0)
#all significant
#compare with numberOfTrainingTrialsByID
mTrainingEffects10=lmer(reactionTime~endTime*group+degY*group+
                          deg*correct_response+deg*endTime+
                          comparisonTimeAvgByID+
                          shortDirectionPropByID*group+
                          numberOfTrainingTrialsByID*group+
                          numberOfPretestTrialsByID+
                          (deg+endTime|ID)+(1|modelNumber),
                        data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects10.summary=modelSummary(mTrainingEffects10,0)
#remove comparisonTimeAvgByID
mTrainingEffects11=lmer(reactionTime~endTime*group+degY*group+
                          deg*correct_response+deg*endTime+
                          shortDirectionPropByID*group+
                          numberOfTrainingTrialsByID*group+
                          numberOfPretestTrialsByID+
                          (deg+endTime|ID)+(1|modelNumber),
                        data=dataset.rt.postTest,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTrainingEffects11.summary=modelSummary(mTrainingEffects11)
#allsignificant
plot(mTrainingEffects11)
save(mTrainingEffects11,mTrainingEffects11.summary,file="statmodels/RTmTrainingEffects.RData")
