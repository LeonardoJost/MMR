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
datasetAnalysis=myData

#plot block*group interaction over time for pre/posttest
myData$cond=paste(myData$group,myData$block,sep="*")
myData$condForLineTypes=myData$group
generateTableAndGraphsForCondition(myData,"datasetTime",TRUE,TRUE,"Group",TRUE)

#plot groups over time for training
myDataTraining$cond=myDataTraining$group
generateTableAndGraphsForCondition(myDataTraining,"TrainingGroup",TRUE,TRUE,"Block")
