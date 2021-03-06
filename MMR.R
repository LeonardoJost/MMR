### main program
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

source("functions/helpers.R")
source("functions/readData.R", encoding="utf-8")
source("functions/generateGraphsAndTables.R", encoding="utf-8")
source("functions/calculateTrainingData.R", encoding="utf-8")

##create output directories, if they don't exist (outputs warnings otherwise)
dir.create("figs")
dir.create("output")
dir.create("figs/MR")
dir.create("figs/MR/allData")
dir.create("figs/MR/meanData")
dir.create("figs/MR/Timed/")
dir.create("figs/MR/accData/")

##options, parameters
options(digits=6)
#set data folder
folder="data\\Logfiles\\"
verbose=3 #detail of output
experimentalSoftware="OpenSesame" #"OpenSesame", no support for presentation atm
questionaireOutFile="output\\questionaire" #.csv added at end, leave empty if no output desired
handednessGraphFile="figs\\HandednessMW.png" #leave empty if no output desired
outlierFactor=3 #factor of sd to define outliers in MR
block=c("preTest","training","postTest")#name of intersting block of data
questionaireDataCols=c("ID","Gender","Experience") #which questionaire columns shall be kept for statistical analysis

##read and write data
#read data
questionaireData=getQuestionaireData(experimentalSoftware,verbose,folder)
MRData=getMRData(experimentalSoftware,verbose,folder,block)
MRDataTraining=getMRData(experimentalSoftware,verbose,folder,"trainingShowStimulus")
#modify data 
questionaireData=modifyQuestionaireData(experimentalSoftware,questionaireData)
MRData=modifyMRData(experimentalSoftware,verbose,MRData,outlierFactor)
#calculate means from questionaire (and save to csv)
calculateMeansQuestionaire(verbose,questionaireData,questionaireOutFile,handednessGraphFile)
#remove not analyzed questionaire data to protect participant identity
questionaireData=subset(questionaireData,select=questionaireDataCols)
#summarize training data (too much data otherwise)
MRDataTrainingSummary=summarizeTrainingData(MRDataTraining,verbose)
#unify data
dataset=merge(MRData,questionaireData,by="ID")
dataset=merge(dataset,MRDataTrainingSummary,by=c("ID","startTime"),all.x=TRUE)
#save average training data to each ID
dataset=summarizeTrainingDataByID(dataset)
#anonymise IDs to protect participant identity
dataset$ID=as.factor(dataset$ID)
levels(dataset$ID)=paste("id",sample.int(length(levels(dataset$ID))),sep="")

#set groupnames by trainingtype
dataset$group=ifelse(dataset$trainingType==1,"wheel",ifelse(dataset$trainingType==2,"buttons","visual"))

#save full dataset to csv
write.table(dataset,file="output\\dataset.csv",sep=";", col.names=NA)

##create datasets for analysis and plot reaction time and accuracy by interesting conditions
source("functions/createDatasets.R", encoding="utf-8")
