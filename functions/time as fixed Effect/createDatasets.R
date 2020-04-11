### plot dataset according to time analysis
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
getReactionTimeDataset=function(myData){
  datasetForLMM=myData
  #scaling
  datasetForLMM$deg=datasetForLMM$deg/100
  datasetForLMM$endTime=datasetForLMM$endTime/30 #30 minutes (time is already in minutes)
  #prepare dataset
  dataset.noOutlier=datasetForLMM[which(!datasetForLMM$outlier),]
  dataset.rt=dataset.noOutlier[which(dataset.noOutlier$typeOutlier=="hit"),]
  dataset.rt$deg=dataset.rt$deg-mean(dataset.rt$deg) #center degree
  #normalizing time and centering degree are necessary to analyze main effects of partial interaction (block*group) when higher-
  #order interactions are present (deg*block*group+time*block*group). Main effects are calculated for value 0
  #0 of degree: average effect due to centering (this is "standard" main effect of removing higher order interaction)
  #0 of time: difference between blocks
  return(dataset.rt)
}


### main script
#load full dataset
myData=read.csv(file="output\\dataset.csv",sep=";")
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

#dataset for analysis
datasetAnalysis=getReactionTimeDataset(myData)

#plot block*group interaction over time for pre/posttest
myData$cond=paste(myData$group,myData$block,sep="*")
myData$condForLineTypes=myData$group
generateTableAndGraphsForCondition(myData,"datasetTime",FALSE,TRUE,"Group",TRUE)

#plot groups over time for training
myDataTraining$cond=myDataTraining$group
generateTableAndGraphsForCondition(myDataTraining,"TrainingGroup",TRUE,TRUE,"Block")
