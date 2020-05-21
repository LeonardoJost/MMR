### Model generation for statistical analysis including 0° cases
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

library(lme4)
library(optimx)
source("functions/helpers.R")

##accuracy
a0=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           ((degY+degZ)+endTime+block+correct_response+
              (degY+degZ)*endTime+(degY+degZ)*block+(degY+degZ)*correct_response+
              endTime*block+endTime*correct_response+
              block*correct_response||ID)+
           ((degY+degZ)+endTime+block+correct_response+Gender+Experience+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
print(summary(a0),corr=FALSE)
summary(rePCA(a0))
#remove correlation
a1=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           ((degY+degZ)+endTime+block+correct_response+
              (degY+degZ)*endTime+(degY+degZ)*block+(degY+degZ)*correct_response+
              endTime*block+endTime*correct_response+
              block*correct_response||ID)+
           ((degY+degZ)+endTime+block+correct_response+Gender+Experience+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a1)
anova(a1, a0)
summary(rePCA(a1))
#remove parameter with correlation 1 and sd 0 or very close or NaN
a2=glmer((type=="hit")~(degY+degZ)*endTime*block*group+(degY+degZ)*endTime*correct_response+Gender*block+Experience*block+
           ((degY+degZ)+endTime+block+correct_response+
              (degY+degZ)*endTime+
              endTime*block+endTime*correct_response||ID)+
           ((degY+degZ)+endTime+Experience+group||modelNumber),family=binomial(),data=dataset.acc,control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
VarCorr(a2)
anova(a2, a1)
summary(rePCA(a2))
#remove parameter with correlation 1 and sd 0 or very close or NaN


##remove nonsignificant fixed effects

#split

#all values significant
#visual inspection of normality
plot(a24)