### statistical analysis of accuracy
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

##accuracy
#base model
aBase=glmer((type=="hit")~endTime*block+
            degY*endTime*group+
            degZ*block+
            (deg+endTime+deg*endTime|ID)+(deg+endTime|modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
aBase.summary=modelSummary(aBase)
save(aBase,aBase.summary,file="statmodels/AccModelaBase.RData")

#split interactions


#nonsignificant effects

