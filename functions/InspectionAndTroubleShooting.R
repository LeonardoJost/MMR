### Script for inspection of dataset
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

#number of each trainingType
trainingTypes=unique(datasetAnalysis[,c("ID","trainingType","group","Gender")])
nrow(trainingTypes)
n_occur=data.frame(table(paste(trainingTypes$trainingType,paste(trainingTypes$group,trainingTypes$Gender))))
print(n_occur)

#number of outliers, incorrect answers, unrotated stimuli for each part
dataset=read.csv(file="dataset\\dataset.csv",sep=";")
#overall outliers
n_occur=data.frame(table(dataset$outlier))
print(n_occur)
#overall incorrect answers and no rotation
datasetNO=dataset[which(!dataset$outlier),]
n_occur=data.frame(table(datasetNO$type))
print(n_occur)
n_occur=data.frame(table(datasetNO$deg))
print(n_occur)
#values for each block
for(thisblock in levels(as.factor(dataset$block))){
  datasetPerBlock=dataset[which(dataset$block==thisblock),]
  print(thisblock)
  print(nrow(datasetPerBlock))
  #outliers
  n_occur=data.frame(table(datasetPerBlock$outlier))
  print(n_occur)
  #incorrect answers and no rotation
  datasetPerBlock=datasetPerBlock[which(!datasetPerBlock$outlier),]
  n_occur=data.frame(table(datasetPerBlock$type))
  print(n_occur)
  n_occur=data.frame(table(datasetPerBlock$deg))
  print(n_occur)
}

#average number of training trials per group
trainingTrialsByID=unique(datasetAnalysis[,c("ID","numberOfTrainingTrialsByID","group")])
for(group in levels(as.factor(trainingTrialsByID$group))) {
  print(group)
  print(meanMode(trainingTrialsByID$numberOfTrainingTrialsByID[which(trainingTrialsByID$group==group)],TRUE,5))
}
  