#number of each trainingType
trainingTypes=unique(dataset[,c("ID","trainingType","group")])
nrow(trainingTypes)
n_occur=data.frame(table(paste(trainingTypes$trainingType,trainingTypes$group)))
print(n_occur)
