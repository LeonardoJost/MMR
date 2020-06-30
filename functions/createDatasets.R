### create and plot dataset according to time analysis
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

### functions
source("functions/helpers.R")
source("functions/generateGraphsAndTables.R", encoding="utf-8")

#load full dataset
myData=read.csv(file="dataset\\dataset.csv",sep=";")
#normalize time to minutes
myData$endTime=myData$endTime/60000 
#get pre and posttest data and separate training
myData$block=toChar(myData$block)
myDataTraining=myData[which(myData$block=="training"),]
myData=myData[which(myData$block!="training"),]
#inspection
myData[which(myData$endTime>10.1),]
#subtract 10 minutes from preTest times to set transition between tests to 0
myData[which(myData$block=="preTest"),"endTime"]=myData[which(myData$block=="preTest"),"endTime"]-10
#separate rotation angle by axis
myData$degY=myData$deg*(myData$axis=="y")
myData$degZ=myData$deg*(myData$axis=="z")
#0,1-> pretest and training
#2,3-> pretest, posttest
#4,5-> only pretest
#6,7-> training, posttest
#8,9-> only training
#10,11-> only posttest
myData$trainedModel=ifelse(myData$modelIndex<2,"preTrai",
                    ifelse(myData$modelIndex<4,"prePost",
                    ifelse(myData$modelIndex<6,"pre",
                    ifelse(myData$modelIndex<8,"traiPost",
                    ifelse(myData$modelIndex<10,"trai","post")))))
myData$trainedModelPost=ifelse(myData$trainedModel=="post",FALSE,TRUE)

#dataset for analysis
datasetAnalysis=myData
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
#center degree
dataset.rt$deg=dataset.rt$deg-mean(dataset.rt$deg) 
dataset.rt$degY=dataset.rt$degY-mean(dataset.rt$degY) 
dataset.rt$degZ=dataset.rt$degZ-mean(dataset.rt$degZ) 
dataset.acc$deg=dataset.acc$deg-mean(dataset.acc$deg) 
dataset.acc$degY=dataset.acc$degY-mean(dataset.acc$degY) 
dataset.acc$degZ=dataset.acc$degZ-mean(dataset.acc$degZ) 
#normalizing time is necessary to analyze main effects of partial interaction (block*group) when higher-
#order interactions are present (time*block*group). Main effects are calculated for value 0
#0 of time: difference between blocks
#degree is centered to the mean for every axis (correct main effects for average but not meaningful intercept)
#split pre- and posttest data
dataset.rt.preTest=dataset.rt[which(dataset.rt$block=="preTest"),]
dataset.rt.postTest=dataset.rt[which(dataset.rt$block=="postTest"),]

#plot block*group interaction over time for pre/posttest
myData$cond=paste(myData$group,myData$block,sep="*")
myData$condForLineTypes=myData$group
generateTableAndGraphsForCondition(myData,"datasetTime",TRUE,TRUE,"Group",TRUE)
#plot axis differences for pre/posttest
myData$cond=paste(myData$axis,myData$block,sep="*")
myData$condForLineTypes=myData$axis
generateTableAndGraphsForCondition(myData,"AxisBlock",TRUE,TRUE,"Axis",TRUE)

#plot groups over time for training
myDataTraining$cond=myDataTraining$group
generateTableAndGraphsForCondition(myDataTraining,"TrainingGroup",TRUE,TRUE,"Group")
#plot firstDeviationTime instead of reaction Time
myDataTraining$reactionTime2=myDataTraining$reactionTime
myDataTraining$reactionTime=myDataTraining$firstDeviationTime
generateTableAndGraphsForCondition(myDataTraining,"TrainingGroupFirstDeviationTime",TRUE,TRUE,"Group")
#plot firstAllowedAnswerTime instead of reaction Time
myDataTraining$reactionTime=myDataTraining$firstAllowedAnswerTime
generateTableAndGraphsForCondition(myDataTraining,"TrainingGroupFirstAllowedAnswerTime",TRUE,TRUE,"Group")
#plot time from first allowed answer until actual reaction
myDataTraining$reactionTime=myDataTraining$reactionTime2-myDataTraining$firstAllowedAnswerTime
generateTableAndGraphsForCondition(myDataTraining,"TrainingGroupFirstAllowedAnswerTimeDifference",TRUE,TRUE,"Group")
#plot rotationSpeed instead of reaction Time
myDataTraining$reactionTime=abs(myDataTraining$rotationSpeed)
generateTableAndGraphsForCondition(myDataTraining,"TrainingGroupRotationSpeed",TRUE,TRUE,"Group")
#plot numberOfSwitches instead of reaction Time
myDataTraining$reactionTime=myDataTraining$numberOfSwitches
generateTableAndGraphsForCondition(myDataTraining,"TrainingGroupNumberOfSwitches",TRUE,TRUE,"Group")

#plot side differences
myData$cond=myData$correct_response
generateTableAndGraphsForCondition(myData,"side")

#plot axis differences

#plot four-way-interaction
