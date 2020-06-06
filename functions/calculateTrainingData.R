### Summarize training data
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

#summarize trainingData
#dataset: dataset containing data while stimuli are shown
#verbose: detail of output
summarizeTrainingData=function(dataset,verbose){
  #rename variables
  dataset$deg=toNumeric(dataset$angle)
  dataset$reactionTime=dataset$response_time
  dataset$angle=NULL
  dataset$response_time=NULL
  #get times
  #name end for each stimulus
  dataset$endTime=dataset$duration+dataset$reactionTime #use time_Stimulus instead of duration to account for framerate of monitor?
  #name startTime for each stimulus
  dataset$startTime=dataset$endTime-dataset$reactionTime
  #create dataset for saving/returning data
  returnset=data.frame(matrix(nrow=length(unique(paste(dataset$ID,dataset$startTime))), ncol=6), stringsAsFactors = FALSE)
  names(returnset)=c("ID","startTime","firstDeviationTime","rotationSpeed","numberOfSwitches","shortDirection")
  if (verbose>1) {
    print(paste("Summarizing training data of", nrow(returnset), "stimuli..."))
  }
  #loop over stimuli
  i=1
  for(thisID in unique(dataset$ID)) {
    dataId=dataset[which(dataset$ID==thisID),]
    if (verbose>2) {
      print(paste("Summarizing training data of participant ID", thisID, "with", length(unique(dataId$startTime)), "stimuli"))
    }
    for(thisStartTime in unique(dataId$startTime)) {
      #get data to this stimulus
      dataIdTime=dataId[which(dataId$startTime==thisStartTime),]
      dataIdTime=dataIdTime[order(dataIdTime$reactionTime),]
      #get first entry -> starting degrees
      startTime=min(dataIdTime$reactionTime)
      startDegrees=dataIdTime$deg[which(dataIdTime$reactionTime==startTime)]
      #get last time where degree==startDegree (start of first deviation -1)
      firstDeviationTime=min(dataIdTime$reactionTime[which(abs(dataIdTime$deg-startDegrees)>1)-1])-startTime
      #get finish time -> angle small enough that answer is allowed (<10 or >350) or jump over 360
      firstAllowedAnswerTime=min(dataIdTime$reactionTime[
        min(which(dataIdTime$deg<10),which(dataIdTime$deg>350),which(abs(diff(dataIdTime$deg))>300)+1)])-startTime
      firstAllowedAnswerDegrees=dataIdTime$deg[which(dataIdTime$reactionTime==firstAllowedAnswerTime+startTime)]
      #calculate the rotation in the "correct" direction and rotation in the "wrong" direction
      #get data for rotation
      dataIdTimeRotation=dataIdTime[which(dataIdTime$reactionTime>=firstDeviationTime+startTime & dataIdTime$reactionTime<=firstAllowedAnswerTime+startTime),]
      #get stepwise rotations
      rotationSteps=diff(dataIdTimeRotation$deg)
      #in case of rotation acroos 0/360°, correct by 360 for analysis (this can only happen in the last step)
      if(rotationSteps[length(rotationSteps)]>300){
        rotationSteps[length(rotationSteps)]=-1
        firstAllowedAnswerDegrees=firstAllowedAnswerDegrees-360
      }
      if(rotationSteps[length(rotationSteps)]<(-300)){
        rotationSteps[length(rotationSteps)]=1
        firstAllowedAnswerDegrees=firstAllowedAnswerDegrees+360
      }
      #select only actual rotations
      rotationSteps=rotationSteps[which(rotationSteps!=0)]
      #get number of switches of rotation direction
      numberOfSwitches=sum(diff(sign(rotationSteps))!=0)
      #rotation speed from first deviation to finish (in °/s)
      rotationSpeed=1000*(firstAllowedAnswerDegrees-startDegrees)/(firstAllowedAnswerTime-firstDeviationTime)
      #was the rotation in the short direction? (180 undefined -> considered later)
      shortDirection=FALSE
      if(sign(startDegrees-180)==sign(rotationSpeed))
        shortDirection=TRUE
      returnset[i,]=list(thisID,thisStartTime,firstDeviationTime,rotationSpeed,numberOfSwitches,shortDirection)
      if (verbose>3) {
        print(paste("Stimulus", i, "of", nrow(returnset), "finished."))
      }
      i=i+1
    }
  }
  if (verbose>1) {
    print("Summarizing training data finished.")
  }
  return(returnset)
}

#summarizes training data by ID
#dataset: dataset containing all summarized training stimuli for one id
summarizeDataByID=function(dataset){
  dataset$firstDeviationTimeAvgByID=mean(dataset$firstDeviationTime)
  dataset$rotationSpeedAvgByID=mean(dataset$rotationSpeed)
  dataset$numberOfSwitchesByID=mean(dataset$numberOfSwitches)
  return(dataset)
}

#dataset=MRDataTraining
#dataIdTime=dataset[which(dataset$ID==1581076965.8099999 & dataset$startTime==873317),]
#print(paste(startTime,startDegrees,firstDeviationTime,firstAllowedAnswerTime,firstAllowedAnswerDegrees,rotationSpeed,numberOfSwitches,sep=", "))
