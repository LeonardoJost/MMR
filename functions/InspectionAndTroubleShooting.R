
trainingTypes=unique(MRData[,c("ID","trainingType")])
nrow(trainingTypes)
n_occur=data.frame(table(trainingTypes$trainingType))
