### Approximation of Bayes factors according to Wagenmakers (2007)
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

#compute BIC difference from LRT
#chiSq=2*(logLik(m1)-logLik(m2))
#BIC(m1)=-2*logLik(m1)+log(n)*df(m1)
#delBIC=BIC(m1)-BIC(m2)
#      =-2*(logLik(m1)-logLik(m2))+log(n)*(df(m1)-df(m2))
#      =-chiSq(m1,m2)+log(n)*delDf
computeDelBICfromLRT=function(chiSq,n,delDf){
  return(-chiSq+log(n)*delDf)
}
#returns bayes factor computed from the difference in BIC
computeBFfromDelBIC=function(delBIC){
  return(exp(delBIC/2))
}
#returns bayes factor in favor of the less complex model (null hypothesis)
computeBFfromLRT=function(chiSq,n,delDf){
  return(computeBFfromDelBIC(computeDelBICfromLRT(chiSq,n,delDf)))
}
#returns bayes factor in favor of the less complex model (null hypothesis)
#m1 should be the alternative, m2 should be the null hypothesis
computeBFfromTwoModels=function(m1,m2,n,delDf){
  return(computeBFfromLRT(-2*(logLik(m2)-logLik(m1)),n,delDf))
}
