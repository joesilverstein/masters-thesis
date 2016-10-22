library("dlm")
library("tseries")
library("stats")
library("forecast")
library("sandwich")
library("lmtest")
library("MTS") # multivariate time series analysis (cite Ruey S. Tsay)
library("splines")
library("boot")
library("ridge")
library("MASS")
library("glmnet")
library("MESS")
library("parcor")

source('C:/Users/Joe/Google Drive/Masters Thesis/adalassomae.R')
source('C:/Users/Joe/Google Drive/Masters Thesis/mylarsmae.R')

setwd("C:/Users/Joe/Google Drive/Masters Thesis")
rm(list=ls()) # clears all environment variables

hpi1q08 = read.csv('historical_data/1q08hpiregpo.txt', sep="\t")
hpi2q08 = read.csv('historical_data/2q08hpiregpo.txt', sep="\t")
hpi3q08 = read.csv('historical_data/3q08hpiregpo.txt', sep="\t")
hpi4q08 = read.csv('historical_data/4q08hpiregpo.txt', sep="\t")
hpi1q09 = read.csv('historical_data/1q09hpiregpo.txt', sep="\t")
hpi2q09 = read.csv('historical_data/2q09hpiregpo.txt', sep="\t")
hpi3q09 = read.csv('historical_data/3q09hpiregpo.txt', sep="\t")
hpi4q09 = read.csv('historical_data/4q09hpiregpo.txt', sep="\t") 
hpi1q10 = read.csv('historical_data/1q10hpiregpo.txt', sep="\t")
hpi2q10 = read.csv('historical_data/2q10hpiregpo.txt', sep="\t")
hpi3q10 = read.csv('historical_data/3q10hpiregpo.txt', sep="\t")
hpi4q10 = read.csv('historical_data/4q10hpiregpo.txt', sep="\t")
hpi1q11 = read.csv('historical_data/1q11hpiregpo.txt', sep="\t")
hpi2q11 = read.csv('historical_data/2q11hpiregpo.txt', sep="\t")
hpi3q11 = read.csv('historical_data/3q11hpiregpo.txt', sep="\t")
hpi4q11 = read.csv('historical_data/4q11hpiregpo.txt', sep="\t")
hpi1q12 = read.csv('historical_data/1q12hpiregpo.txt', sep="\t")
hpi2q12 = read.csv('historical_data/2q12hpiregpo.txt', sep="\t")
hpi3q12 = read.csv('historical_data/3q12hpiregpo.txt', sep="\t")
hpi4q12 = read.csv('historical_data/4q12hpiregpo.txt', sep="\t")
hpi1q13 = read.csv('historical_data/1q13hpiregpo.txt', sep="\t")
hpi2q13 = read.csv('historical_data/2q13hpiregpo.txt', sep="\t")
hpi3q13 = read.csv('historical_data/3q13hpiregpo.txt', sep="\t")
hpi4q13 = read.csv('historical_data/4q13hpiregpo.txt', sep="\t")
hpi1q14 = read.csv('historical_data/1q14hpiregpo.txt', sep="\t")
hpi2q14 = read.csv('historical_data/2q14hpiregpo.txt', sep="\t")
hpi3q14 = read.csv('historical_data/3q14hpiregpo.txt', sep="\t")
hpi4q14 = read.csv('historical_data/4q14hpiregpo.txt', sep="\t")
hpi1q15 = read.csv('historical_data/1q15hpiregpo.txt', sep="\t")
hpi2q15 = read.csv('historical_data/2q15hpiregpo.txt', sep="\t")

# States
hpiS4q09 = read.csv('historical_data/4q09hpistspo.txt', sep="\t") 
hpiS1q10 = read.csv('historical_data/1q10hpistspo.txt', sep="\t")
hpiS2q10 = read.csv('historical_data/2q10hpistspo.txt', sep="\t")
hpiS3q10 = read.csv('historical_data/3q10hpistspo.txt', sep="\t")
hpiS4q10 = read.csv('historical_data/4q10hpistspo.txt', sep="\t")
hpiS1q11 = read.csv('historical_data/1q11hpistspo.txt', sep="\t")
hpiS2q11 = read.csv('historical_data/2q11hpistspo.txt', sep="\t")
hpiS3q11 = read.csv('historical_data/3q11hpistspo.txt', sep="\t")
hpiS4q11 = read.csv('historical_data/4q11hpistspo.txt', sep="\t")
hpiS1q12 = read.csv('historical_data/1q12hpistspo.txt', sep="\t")
hpiS2q12 = read.csv('historical_data/2q12hpistspo.txt', sep="\t")
hpiS3q12 = read.csv('historical_data/3q12hpistspo.txt', sep="\t")
hpiS4q12 = read.csv('historical_data/4q12hpistspo.txt', sep="\t")
hpiS1q13 = read.csv('historical_data/1q13hpistspo.txt', sep="\t")

names(hpi1q08)[5] = "hpi1q08"
names(hpi2q08)[5] = "hpi2q08"
names(hpi3q08)[5] = "hpi3q08"
names(hpi4q08)[5] = "hpi4q08"
names(hpi1q09)[5] = "hpi1q09"
names(hpi2q09)[5] = "hpi2q09"
names(hpi3q09)[5] = "hpi3q09"
names(hpi4q09)[5] = "hpi4q09"
names(hpi1q10)[5] = "hpi1q10"
names(hpi2q10)[5] = "hpi2q10"
names(hpi3q10)[5] = "hpi3q10"
names(hpi4q10)[5] = "hpi4q10"
names(hpi1q11)[5] = "hpi1q11"
names(hpi2q11)[5] = "hpi2q11"
names(hpi3q11)[5] = "hpi3q11"
names(hpi4q11)[5] = "hpi4q11"
names(hpi1q12)[5] = "hpi1q12"
names(hpi2q12)[5] = "hpi2q12"
names(hpi3q12)[5] = "hpi3q12"
names(hpi4q12)[5] = "hpi4q12"
names(hpi1q13)[5] = "hpi1q13"
names(hpi2q13)[5] = "hpi2q13"
names(hpi3q13)[5] = "hpi3q13"
names(hpi4q13)[5] = "hpi4q13"
names(hpi1q14)[5] = "hpi1q14"
names(hpi2q14)[5] = "hpi2q14"
names(hpi3q14)[5] = "hpi3q14"
names(hpi4q14)[5] = "hpi4q14"
names(hpi1q15)[5] = "hpi1q15"
names(hpi2q15)[5] = "hpi2q15"

names(hpiS4q09)[5] = "hpi4q09"
names(hpiS1q10)[5] = "hpi1q10"
names(hpiS2q10)[5] = "hpi2q10"
names(hpiS3q10)[5] = "hpi3q10"
names(hpiS4q10)[5] = "hpi4q10"
names(hpiS1q11)[5] = "hpi1q11"
names(hpiS2q11)[5] = "hpi2q11"
names(hpiS3q11)[5] = "hpi3q11"
names(hpiS4q11)[5] = "hpi4q11"
names(hpiS1q12)[5] = "hpi1q12"
names(hpiS2q12)[5] = "hpi2q12"
names(hpiS3q12)[5] = "hpi3q12"
names(hpiS4q12)[5] = "hpi4q12"
names(hpiS1q13)[5] = "hpi1q13"

hpi1q08$index_po_nsa = NULL
hpi2q08$index_po_nsa = NULL
hpi3q08$index_po_nsa = NULL
hpi4q08$index_po_nsa = NULL
hpi1q09$index_po_nsa = NULL
hpi2q09$index_po_not_seasonally_adjusted = NULL
hpi3q09$index_po_not_seasonally_adjusted = NULL
hpi4q09$index_po_not_seasonally_adjusted = NULL
hpi1q10$index_po_not_seasonally_adjusted = NULL
hpi2q10$index_po_not_seasonally_adjusted = NULL
hpi3q10$index_po_not_seasonally_adjusted = NULL
hpi4q10$index_po_not_seasonally_adjusted = NULL
hpi1q11$index_po_not_seasonally_adjusted = NULL
hpi2q11$index_po_not_seasonally_adjusted = NULL
hpi3q11$index_po_not_seasonally_adjusted = NULL
hpi4q11$index_po_not_seasonally_adjusted = NULL
hpi1q12$index_po_not_seasonally_adjusted = NULL
hpi2q12$index_po_not_seasonally_adjusted = NULL
hpi3q12$index_po_not_seasonally_adjusted = NULL
hpi4q12$index_po_not_seasonally_adjusted = NULL
hpi1q13$index_po_not_seasonally_adjusted = NULL
hpi2q13$index_po_not_seasonally_adjusted = NULL
hpi3q13$index_po_not_seasonally_adjusted = NULL
hpi4q13$index_po_not_seasonally_adjusted = NULL
hpi1q14$index_po_not_seasonally_adjusted = NULL
hpi2q14$index_po_not_seasonally_adjusted = NULL
hpi3q14$index_po_not_seasonally_adjusted = NULL
hpi4q14$index_po_not_seasonally_adjusted = NULL
hpi1q15$index_po_not_seasonally_adjusted = NULL
hpi2q15$index_po_not_seasonally_adjusted = NULL

hpiS4q09$index_nsa = NULL
hpiS1q10$index_nsa = NULL
hpiS2q10$index_nsa = NULL
hpiS3q10$index_nsa = NULL
hpiS4q10$index_nsa = NULL
hpiS1q11$index_nsa = NULL
hpiS2q11$index_nsa = NULL
hpiS3q11$index_nsa = NULL
hpiS4q11$index_nsa = NULL
hpiS1q12$index_nsa = NULL
hpiS2q12$index_nsa = NULL
hpiS3q12$index_nsa = NULL
hpiS4q12$index_nsa = NULL
hpiS1q13$index_nsa = NULL

hpiList = list(hpi1q08, hpi2q08, hpi3q08, hpi4q08,
               hpi1q09, hpi2q09, hpi3q09, hpi4q09,
               hpi1q10, hpi2q10, hpi3q10, hpi4q10,
               hpi1q11, hpi2q11, hpi3q11, hpi4q11,
               hpi1q12, hpi2q12, hpi3q12, hpi4q12,
               hpi1q13, hpi2q13, hpi3q13, hpi4q13,
               hpi1q14, hpi2q14, hpi3q14, hpi4q14,
               hpi1q15, hpi2q15)

hpiSList = list(hpiS4q09,
                hpiS1q10, hpiS2q10, hpiS3q10, hpiS4q10,
                hpiS1q11, hpiS2q11, hpiS3q11, hpiS4q11,
                hpiS1q12, hpiS2q12, hpiS3q12, hpiS4q12,
                hpiS1q13)

hpi = Reduce(function(x, y) merge(x, y, by=c("division","year","qtr"), all.y=TRUE), hpiList)

hpiENC = subset(hpi,division=="DV_ENC")
rownames(hpiENC) = NULL

hpiESC = subset(hpi,division=="DV_ESC")
rownames(hpiESC) = NULL

hpiMA = subset(hpi,division=="DV_MA")
rownames(hpiMA) = NULL

hpiMT = subset(hpi,division=="DV_MT")
rownames(hpiMT) = NULL

hpiNE = subset(hpi,division=="DV_NE")
rownames(hpiNE) = NULL

hpiPAC = subset(hpi,division=="DV_PAC")
rownames(hpiPAC) = NULL

hpiSA = subset(hpi,division=="DV_SA")
rownames(hpiSA) = NULL

hpiWNC = subset(hpi,division=="DV_WNC")
rownames(hpiWNC) = NULL

hpiWSC = subset(hpi,division=="DV_WSC")
rownames(hpiWSC) = NULL

hpiUSA = subset(hpi, division=="USA")
rownames(hpiUSA) = NULL

# # Don't run the following merge for now because it is too slow
# hpiS = Reduce(function(x, y) merge(x, y, by=c("state","yr","qtr"), all.y=TRUE), hpiSList)
# hpiS = subset(hpiS, division=="USA")
# rownames(hpiS) = NULL

hpiMergedDivList = list(hpiENC, hpiESC, hpiMA, hpiMT, hpiNE, hpiPAC, hpiSA, hpiWNC, hpiWSC, hpiUSA)

for (i in 1:length(hpiMergedDivList)) {
  tmp = data.frame(hpiMergedDivList[i])
  
  # Create variables for next first revisions and first revisions one, two, three, and four periods back
  tmp[,"rev0"] = NA
  tmp[,"rev1"] = NA
  tmp[,"rev2"] = NA
  tmp[,"rev3"] = NA
  tmp[,"rev4"] = NA
  x = 8
  for (j in 73:98) {
    tmp[j-1,34] = tmp[j-1,x] / tmp[j-1,x-1] - 1
    tmp[j-1,35] = tmp[j-2,x-1] / tmp[j-2,x-2] - 1
    tmp[j-1,36] = tmp[j-3,x-2] / tmp[j-3,x-3] - 1
    tmp[j-1,37] = tmp[j-4,x-3] / tmp[j-4,x-4] - 1
    tmp[j-1,38] = tmp[j-5,x-4] / tmp[j-5,x-5] - 1
    x = x + 1
  }
  
  # Create variables for first, second, third, and fourth lagged house price appreciation:
  tmp[,"app1"] = NA
  tmp[,"app2"] = NA
  tmp[,"app3"] = NA
  tmp[,"app4"] = NA
  tmp[,"app5"] = NA
  tmp[,"app6"] = NA
  tmp[,"app7"] = NA
  tmp[,"app8"] = NA
  x = 7
  for (j in 72:97) {
    tmp[j,39] = tmp[j,x] / tmp[j-1,x] - 1
    tmp[j,40] = tmp[j-1,x] / tmp[j-2,x] - 1
    tmp[j,41] = tmp[j-2,x] / tmp[j-3,x] - 1
    tmp[j,42] = tmp[j-3,x] / tmp[j-4,x] - 1
    tmp[j,43] = tmp[j-4,x] / tmp[j-5,x] - 1
    tmp[j,44] = tmp[j-5,x] / tmp[j-6,x] - 1
    tmp[j,45] = tmp[j-6,x] / tmp[j-7,x] - 1
    tmp[j,46] = tmp[j-7,x] / tmp[j-8,x] - 1
    x = x + 1
  }
  
  hpiMergedDivList[i] = list(tmp)
}

tmp1 = data.frame(hpiMergedDivList[1])
tmp2 = data.frame(hpiMergedDivList[2])
cor(na.omit(tmp1$app1), na.omit(tmp2$app1))
cor(na.omit(tmp1$app2), na.omit(tmp2$app2))
cor(na.omit(tmp1$app3), na.omit(tmp2$app3))
cor(na.omit(tmp1$app4), na.omit(tmp2$app4))
cor(na.omit(tmp1$rev1), na.omit(tmp2$rev1))

# Want to see if the state revisions are correlated, to see if we can get more predictability by running a VAR
# A: They are correlated, so should run a VAR

### Decide on the appropriate model ###

# Test to make sure the series of revisions is stationary
# Start with first division in the dataset
div1 = data.frame(hpiMergedDivList[1])
rev1 = na.omit(div1$rev1)
adf.test(rev1)
# since p > 0.05, I accept the null hypothesis of non-stationarity.
# Need to take first difference:
diffRev1 = diff(rev1)
adf.test(diffRev1)
# Still non-stationary. Need to take 2nd difference.
diff2Rev1 = diff(diffRev1)
adf.test(na.omit(diff2Rev1))

app1 = na.omit(div1$app1)

# Start by looking at the autocorrelation and partial autocorrelation functions
acf(rev1)$acf # Very little autocorrelation. This leads me to believe revisions should not be modeled as an autoregression.
pacf(rev1)
cor(rev1,app1) # Somewhat correlated. I think there is enough correlation to use it for prediction.
acf(app1)$acf # very little autocorrelation - no discernable pattern
pacf(app1)

# Formally test for autocorrelation in rev1 and app1 using the Ljung-Box test:
Box.test(rev1, type="Ljung-Box") # Since p > 0.05, I accept the null hypothesis of independence
Box.test(app1, type="Ljung-Box") # Since p > 0.05, I accept the null hypothesis of independence

# Therefore, I conclude that both time series are in fact a series of independent observations.

# Since the series rev1 of first revisions exhibits little autocorrelation, should not run an autoregression.
# However, since rev1 is correlated with app1, we can just regress rev1 on app1. 
# Since app1 exhibits no autocorrelation and we want to avoid overfitting, just use the most recent lag of app1 in predicting rev1

### Estimate the Model ###

# Since app1 is highly correlated across divisions, we should run a vector regression.

# first run OLS on just first division's data to just take an initial glance at the results:

reg = lm(rev1 ~ app1 + app2 + app3 + app4, data=div1)
summary(reg)
newdata = data.frame(app1=0.01, app2=0.005, app3=0.013, app4=0.02)
predict(reg, newdata, interval="confidence")
predict(reg, newdata, interval="prediction")

# FIRST, TEST FOR THE DEGREE OF HETEROSKEDASTICITY JUST TO GET AN IDEA.
# THEN, SHOULD MODIFY THE CODE TO USE HETEROSKEDASTICITY-CONSISTENT STANDARD ERRORS!!!

# Now run a multivariate linear regression using data for all the divisions to hopefully get a more accurate prediction:

# Construct the vectors to be used in the regression:

div1 = data.frame(hpiMergedDivList[1])
rev0div1 = na.omit(div1$rev0)
rev1div1 = na.omit(div1$rev1)
rev2div1 = na.omit(div1$rev2)
rev3div1 = na.omit(div1$rev3)
rev4div1 = na.omit(div1$rev4)
app1div1 = na.omit(div1$app1)
app2div1 = na.omit(div1$app2)
app3div1 = na.omit(div1$app3)
app4div1 = na.omit(div1$app4)
app5div1 = na.omit(div1$app5)
app6div1 = na.omit(div1$app6)
app7div1 = na.omit(div1$app7)
app8div1 = na.omit(div1$app8)

div2 = data.frame(hpiMergedDivList[2])
rev0div2 = na.omit(div2$rev0)
rev1div2 = na.omit(div2$rev1)
rev2div2 = na.omit(div2$rev2)
rev3div2 = na.omit(div2$rev3)
rev4div2 = na.omit(div2$rev4)
app1div2 = na.omit(div2$app1)
app2div2 = na.omit(div2$app2)
app3div2 = na.omit(div2$app3)
app4div2 = na.omit(div2$app4)
app5div2 = na.omit(div2$app5)
app6div2 = na.omit(div2$app6)
app7div2 = na.omit(div2$app7)
app8div2 = na.omit(div2$app8)

div3 = data.frame(hpiMergedDivList[3])
rev0div3 = na.omit(div3$rev0)
rev1div3 = na.omit(div3$rev1)
rev2div3 = na.omit(div3$rev2)
rev3div3 = na.omit(div3$rev3)
rev4div3 = na.omit(div3$rev4)
app1div3 = na.omit(div3$app1)
app2div3 = na.omit(div3$app2)
app3div3 = na.omit(div3$app3)
app4div3 = na.omit(div3$app4)
app5div3 = na.omit(div3$app5)
app6div3 = na.omit(div3$app6)
app7div3 = na.omit(div3$app7)
app8div3 = na.omit(div3$app8)

div4 = data.frame(hpiMergedDivList[4])
rev0div4 = na.omit(div4$rev0)
rev1div4 = na.omit(div4$rev1)
rev2div4 = na.omit(div4$rev2)
rev3div4 = na.omit(div4$rev3)
rev4div4 = na.omit(div4$rev4)
app1div4 = na.omit(div4$app1)
app2div4 = na.omit(div4$app2)
app3div4 = na.omit(div4$app3)
app4div4 = na.omit(div4$app4)
app5div4 = na.omit(div4$app5)
app6div4 = na.omit(div4$app6)
app7div4 = na.omit(div4$app7)
app8div4 = na.omit(div4$app8)

div5 = data.frame(hpiMergedDivList[5])
rev0div5 = na.omit(div5$rev0)
rev1div5 = na.omit(div5$rev1)
rev2div5 = na.omit(div5$rev2)
rev3div5 = na.omit(div5$rev3)
rev4div5 = na.omit(div5$rev4)
app1div5 = na.omit(div5$app1)
app2div5 = na.omit(div5$app2)
app3div5 = na.omit(div5$app3)
app4div5 = na.omit(div5$app4)
app5div5 = na.omit(div5$app5)
app6div5 = na.omit(div5$app6)
app7div5 = na.omit(div5$app7)
app8div5 = na.omit(div5$app8)

div6 = data.frame(hpiMergedDivList[6])
rev0div6 = na.omit(div6$rev0)
rev1div6 = na.omit(div6$rev1)
rev2div6 = na.omit(div6$rev2)
rev3div6 = na.omit(div6$rev3)
rev4div6 = na.omit(div6$rev4)
app1div6 = na.omit(div6$app1)
app2div6 = na.omit(div6$app2)
app3div6 = na.omit(div6$app3)
app4div6 = na.omit(div6$app4)
app5div6 = na.omit(div6$app5)
app6div6 = na.omit(div6$app6)
app7div6 = na.omit(div6$app7)
app8div6 = na.omit(div6$app8)

div7 = data.frame(hpiMergedDivList[7])
rev0div7 = na.omit(div7$rev0)
rev1div7 = na.omit(div7$rev1)
rev2div7 = na.omit(div7$rev2)
rev3div7 = na.omit(div7$rev3)
rev4div7 = na.omit(div7$rev4)
app1div7 = na.omit(div7$app1)
app2div7 = na.omit(div7$app2)
app3div7 = na.omit(div7$app3)
app4div7 = na.omit(div7$app4)
app5div7 = na.omit(div7$app5)
app6div7 = na.omit(div7$app6)
app7div7 = na.omit(div7$app7)
app8div7 = na.omit(div7$app8)

div8 = data.frame(hpiMergedDivList[8])
rev0div8 = na.omit(div8$rev0)
rev1div8 = na.omit(div8$rev1)
rev2div8 = na.omit(div8$rev2)
rev3div8 = na.omit(div8$rev3)
rev4div8 = na.omit(div8$rev4)
app1div8 = na.omit(div8$app1)
app2div8 = na.omit(div8$app2)
app3div8 = na.omit(div8$app3)
app4div8 = na.omit(div8$app4)
app5div8 = na.omit(div8$app5)
app6div8 = na.omit(div8$app6)
app7div8 = na.omit(div8$app7)
app8div8 = na.omit(div8$app8)

div9 = data.frame(hpiMergedDivList[9])
rev0div9 = na.omit(div9$rev0)
rev1div9 = na.omit(div9$rev1)
rev2div9 = na.omit(div9$rev2)
rev3div9 = na.omit(div9$rev3)
rev4div9 = na.omit(div9$rev4)
app1div9 = na.omit(div9$app1)
app2div9 = na.omit(div9$app2)
app3div9 = na.omit(div9$app3)
app4div9 = na.omit(div9$app4)
app5div9 = na.omit(div9$app5)
app6div9 = na.omit(div9$app6)
app7div9 = na.omit(div9$app7)
app8div9 = na.omit(div9$app8)

cor(rev1div1, rev2div1) # significant correlation, so lagged first revisions may be worth including
cor(rev1div1, rev3div1)
cor(rev1div1, rev4div1)

DIV = data.frame(rev0div1, rev0div2, rev0div3, rev0div4, rev0div5, rev0div6, rev0div7, rev0div8, rev0div9,
                 rev1div1, rev1div2, rev1div3, rev1div4, rev1div5, rev1div6, rev1div7, rev1div8, rev1div9,
                 rev2div1, rev2div2, rev2div3, rev2div4, rev2div5, rev2div6, rev2div7, rev2div8, rev2div9,
                 rev3div1, rev3div2, rev3div3, rev3div4, rev3div5, rev3div6, rev3div7, rev3div8, rev3div9,
                 rev4div1, rev4div2, rev4div3, rev4div4, rev4div5, rev4div6, rev4div7, rev4div8, rev4div9,
                 app1div1, app1div2, app1div3, app1div4, app1div5, app1div6, app1div7, app1div8, app1div9,
                 app2div1, app2div2, app2div3, app2div4, app2div5, app2div6, app2div7, app2div8, app2div9,
                 app3div1, app3div2, app3div3, app3div4, app3div5, app3div6, app3div7, app3div8, app3div9,
                 app4div1, app4div2, app4div3, app4div4, app4div5, app4div6, app4div7, app4div8, app4div9,
                 app5div1, app5div2, app5div3, app5div4, app5div5, app5div6, app5div7, app5div8, app5div9,
                 app6div1, app6div2, app6div3, app6div4, app6div5, app6div6, app6div7, app6div8, app6div9,
                 app7div1, app7div2, app7div3, app7div4, app7div5, app7div6, app7div7, app7div8, app7div9,
                 app8div1, app8div2, app8div3, app8div4, app8div5, app8div6, app8div7, app8div8, app8div9)

### First do cross-validation on the Elul et. al. model to show it doesn't work:

costFun = function(y, yhat) {
  return(abs(y-yhat))
}

reg1 = glm(rev0div1 ~ rev1div1 + rev2div1 + rev3div1 + rev4div1 + app1div1 + app2div1 + app3div1 + app4div1)
cv.glm(data=DIV, glmfit=reg1, K=length(rev0div1), cost=costFun)$delta[1]

reg2 = glm(rev0div2 ~ rev1div2 + rev2div2 + rev3div2 + rev4div2 + app1div2 + app2div2 + app3div2 + app4div2)
cv.glm(data=DIV, glmfit=reg2, K=length(rev0div1), cost=costFun)$delta[1]

reg3 = glm(rev0div3 ~ rev1div3 + rev2div3 + rev3div3 + rev4div3 + app1div3 + app2div3 + app3div3 + app4div3)
cv.glm(data=DIV, glmfit=reg3, K=length(rev0div1), cost=costFun)$delta[1]

reg4 = glm(rev0div4 ~ rev1div4 + rev2div4 + rev3div4 + rev4div4 + app1div4 + app2div4 + app3div4 + app4div4)
cv.glm(data=DIV, glmfit=reg4, K=length(rev0div1), cost=costFun)$delta[1]

reg5 = glm(rev0div5 ~ rev1div5 + rev2div5 + rev3div5 + rev4div5 + app1div5 + app2div5 + app3div5 + app4div5)
cv.glm(data=DIV, glmfit=reg5, K=length(rev0div1), cost=costFun)$delta[1]

reg6 = glm(rev0div6 ~ rev1div6 + rev2div6 + rev3div6 + rev4div6 + app1div6 + app2div6 + app3div6 + app4div6)
cv.glm(data=DIV, glmfit=reg6, K=length(rev0div1), cost=costFun)$delta[1]

reg7 = glm(rev0div7 ~ rev1div7 + rev2div7 + rev3div7 + rev4div7 + app1div7 + app2div7 + app3div7 + app4div7)
cv.glm(data=DIV, glmfit=reg7, K=length(rev0div1), cost=costFun)$delta[1]

reg8 = glm(rev0div8 ~ rev1div8 + rev2div8 + rev3div8 + rev4div8 + app1div8 + app2div8 + app3div8 + app4div8)
cv.glm(data=DIV, glmfit=reg8, K=length(rev0div1), cost=costFun)$delta[1]

reg9 = glm(rev0div9 ~ rev1div9 + rev2div9 + rev3div9 + rev4div9 + app1div9 + app2div9 + app3div9 + app4div9)
cv.glm(data=DIV, glmfit=reg9, K=length(rev0div1), cost=costFun)$delta[1]

###

# Don't run a multivariate regression since the predictions will be the same.
# Instead, regress the revision in each division on the appreciation in all the divisions, to take advantage of the correlation in appreciation accross divisions.
# It makes sense that appreciation is correlated, since it is driven by the same forces across divisions.
# Question: Should we use confidence or prediction intervals?
reg1 = glm(rev1div1 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9
                    + app2div1 + app2div2 + app2div3 + app2div4 + app2div5 + app2div6 + app2div7 + app2div8 + app2div9
                    + app3div1 + app3div2 + app3div3 + app3div4 + app3div5 + app3div6 + app3div7 + app3div8 + app3div9
                    + app4div1 + app4div2 + app4div3 + app4div4 + app4div5 + app4div6 + app4div7 + app4div8 + app4div9,
                    data=DIV)
reg1 = glm(rev1div1 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
newdata = DIV1[5,]
# There might be a bug in the ridge:::predict.ridgeLinear code. A potential fix is below:
# http://stackoverflow.com/questions/31663440/predict-with-linearridge-error-in-as-matrixmm-beta-non-conformable-arg
# This bug seems to be pretty serious. Can instead use lm.ridge in MASS package or use rms package
predict(reg1, newdata, se.fit=TRUE)
# Now the prediction is statistically significant!
# REMEMBER TO GET PREDICT TO PRODUCE ROBUST STANDARD ERRORS, SINCE THIS IS NOT THE DEFAULT!!!
reg2 = glm(rev1div2 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
predict(reg2, newdata, se.fit=TRUE)
reg3 = glm(rev1div3 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
predict(reg3, newdata, se.fit=TRUE)
reg4 = glm(rev1div4 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
predict(reg4, newdata, se.fit=TRUE)
reg5 = glm(rev1div5 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
predict(reg5, newdata, se.fit=TRUE)
reg6 = glm(rev1div6 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
predict(reg6, newdata, se.fit=TRUE)
reg7 = glm(rev1div7 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
predict(reg7, newdata, se.fit=TRUE)
reg8 = glm(rev1div8 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
predict(reg8, newdata, se.fit=TRUE)
reg9 = glm(rev1div9 ~ app1div1 + app1div2 + app1div3 + app1div4 + app1div5 + app1div6 + app1div7 + app1div8 + app1div9, data=DIV)
predict(reg9, newdata, se.fit=TRUE)

# Conclusion: the prediction is usually statistically significant (for most of the divisions)

# Now do exhaustive cross validation.
# cv.glm does "leave-one-out cross validation" (LOOCV), estimating the regression model on n-1 obs and then testing on the remaining ob in all n=22 ways
# http://robjhyndman.com/hyndsight/crossvalidation/
# In the above webpage, Hyndman explains that there is a closed-form solution for the LOOCV statistic for a linear model with independent observations. Since we tested and saw that the obs are independent, we can use this closed form solution with cv.glm. Thus, it doesn't need to be re-estimated K times.
# Note that since the model is estimated on a smaller sample each time, the estimates of prediction error will tend to be biased slightly upward. They are an upper bound.
# Have to run the regression using glm instead of lm, otherwise $delta=NaN
# cv.glm uses the average squared error function as its cost function
# cv.glm: The data is divided randomly into K=22 groups. For each group the generalized linear model is fit to data omitting that group, then the function cost is applied to the observed responses in the group that was omitted from the fit and the prediction made by the fitted models for those observations.
# The cost function is mean-squared-error: cost = function(y, yhat) mean((y - yhat)^2)
# The average mean-squared prediction errors are reported below for each model, along with the average error:

# Should also display average absolute revision size to get an idea of how far-off the predictions tend to be
mse1 = cv.glm(DIV, reg1)$delta[1]
(e1 = sqrt(mse1)) # RMSE of prediction during cross-validation (not so useful)
(absRev1 = mean(abs(rev0div1))) # mean absolute revision size in division 1 (compare with MAE in a table side-by-side later)
mse2 = cv.glm(DIV, reg2)$delta[1]
(e2 = sqrt(mse2))
(absRev2 = mean(abs(rev0div2)))
mse3 = cv.glm(DIV, reg3)$delta[1]
(e3 = sqrt(mse3))
(absRev3 = mean(abs(rev0div3)))
mse4 = cv.glm(DIV, reg4)$delta[1]
(e4 = sqrt(mse4))
(absRev4 = mean(abs(rev0div4)))
mse5 = cv.glm(DIV, reg5)$delta[1]
(e5 = sqrt(mse5))
(absRev5 = mean(abs(rev0div5)))
mse6 = cv.glm(DIV, reg6)$delta[1]
(e6 = sqrt(mse6))
(absRev6 = mean(abs(rev0div6)))
mse7 = cv.glm(DIV, reg7)$delta[1]
(e7 = sqrt(mse7))
(absRev7 = mean(abs(rev0div7)))
mse8 = cv.glm(DIV, reg8)$delta[1]
(e8 = sqrt(mse8))
(absRev8 = mean(abs(rev0div8)))
mse9 = cv.glm(DIV, reg9)$delta[1]
(e9 = sqrt(mse9))
(absRev9 = mean(abs(rev0div9)))
# Note: RMSE is not the same as average absolute error, since the squaring overweights outliers:
# http://people.duke.edu/~rnau/compare.htm

# Note: model no longer identified -- should use PCR or something like it to reduce the dimensionality.
# Description of dimensionality-reduction techniques linked below:
# Diebold's blog (PCR vs. Ridge Regression)
# Diebold says to use ridge: "ridge effectively includes all PC's and shrinks according to sizes of associated eigenvalues, whereas PCR effectively shrinks some PCs completely to zero (those not included) and doesn't shrink others at all (those included). "
# http://fxdiebold.blogspot.com/2015/10/whither-econometric-principal.html
# https://computation.llnl.gov/casc/sapphire/pubs/148494.pdf
# Elements of Statistical Inference (Hastie et. al.): see section on Ridge regression and equation (3.4.1), (3.4.4)
# As discussed in the book, there are 3 main choices for dimensionality reduction: Subset selection using PCA, Ridge regression, or the Lasso.
# In the ridge package, the ridge parameter is selected using the method of Cule et. al. (2012)
# If it is computationally feasible, the ridge parameter can be selected to minimize the LOOCV MSE using the method of Golub et. al. (1979)
# Maybe should instead use the LASSO (cite Tibshirani) because it does both shrinkage and selection, whereas Ridge only does shrinkage
# There is also a Bayesian LASSO (Park and Casella 2008)
# There is also the Elastic Net, which is a compromise between ridge and lasso
# The elastic net can yield more that N nonzero coefficients when p > N, a potential advantage over the lasso.
# If p>>N, the elastic net is better (see Hui and Zou 2005)
# There is also a Bayesian elastic net and an R package for it (2010 paper)
# See R package EBglmnet
# Note: Empirical Bayes is Bayes in which the prior distribution is estimated from the data
# Start with the Lasso, and then move to elastic net and empirical Bayes if time
# http://users.stat.umn.edu/~zouxx019/software.html
# "The Lasso shrinks AND selects: (Diebold Forecasting book)
# Should be using the elastic net for the following reason (from Diebold forecasting book):
# "Unlike Lasso, it moves strongly correlated predictors in or out of the model together, hopefully producing improving prediction accuracy relative to Lasso."
# There is also the "adaptive Lasso" and the "adaptive elastic net", which weight the term penaliztions differently and have the "oracle property"
# Def. of Oracle Property: same as consistency. See page 112 of Diebold's Forecasting book
# Conclusion: Should use the adaptive elastic net (Diebold Forecasting p. 120)

# reg1 = linearRidge(rev1div1 ~ ., data=DIV1)
# reg1 = lm.ridge(rev1div1 ~ ., data=DIV1) # lambda=0 by default, making it the same as OLS

# From Chernozukhov JEP paper (2014):
# LASSO-type estimators have also been shown to have appealing properties
# under plausible assumptions that allow for approximation errors, heteroskedasticity,
# clustering and fixed effects, and non-normality

### Now I do the same using the HPI data ###

x = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,1])

#note alpha=1 for lasso only and can blend with ridge penalty down to alpha=0 ridge only
mod1 = glmnet(x,y,alpha=1,family='gaussian') # for continuous non-censored response variable, use family='gaussian'

#plot variable coefficients vs. shrinkage parameter lambda.
plot(mod1,xvar="lambda")
grid()

#coefficents can be extracted from the glmmod. Here shown with 3 variables selected.
coef(mod1)[,10]

# Lastly, cross validation can also be used to select lambda.
cv.mod1 = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV)) # should extend the lambda vector to larger values of lambda to ensure there is a local min
plot(cv.mod1)
(best_lambda = cv.mod1$lambda.min) # MSE is only slightly lower than before (when using fewer variables)
# want to extend the lambda vector to see whether it is really a local min
logLambda = log(cv.mod1$lambda)
diff = logLambda[2] - logLambda[1]
tmp = seq(from = logLambda[1] - diff*50, to = logLambda[1] - diff, by = diff)
logLambda2 = c(tmp, logLambda) # the larger values of lambda don't seem to be getting used correctly (the CV plot is flat on the RHS)
cv.mod1 = cv.glmnet(x, y, alpha=1, lambda=exp(logLambda2), nfolds=nrow(DIV))
plot(cv.mod1)
(minRMSE = sqrt(min(cv.mod1$cvm))) # CV RMSE

# Now do CV using mean absolute error (MAE), which is more easily comparable to revision size
weights <- adaptive.weights(x, y)
cv.mod1_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor=weights$weights)
plot(cv.mod1_mae)
(minMAE1 = min(cv.mod1_mae$cvm)) # This is much lower, implying that there are a lot of outliers that it does not predict well that are getting overweighted when using RMSE
(absRev1 = mean(abs(rev1div1))) # For comparison

# Including all of the lagged first house price revisions lowers the number of obs down from 22 to 19, so it actually makes things worse. Really need to get the missing quarter's observations from Will.
# NOTE THAT WE CAN USE ALMOST ALL THE DATA AND DON'T HAVE TO WORRY ABOUT THE MISSING QUARTER IF WE DON'T INCLUDE THE LAGGED FIRST REVISIONS AS PREDICTORS
# MAYBE WANT TO JUST USE 1 LAGGED REVISION AND 1 LAGGED APPRECIATION TO KEEP THE MODEL PARSIMONIOUS

# My prior beliefs are that all of the regressiors should have some predictive value, so none of the coefficients should be shrunk completely to 0. However, this would not result in a sparse model.
coef(mod1)[,45] # why are all the coefficients 0 except for app1div9, even when running ridge?
# Also, why are most of the lasso coefficients the same as the ridge coefficients?

# Also, see "MODEL SELECTION WHEN THE NUMBER OF VARIABLES EXCEEDS THE NUMBER OF OBSERVATIONS" (Stodden 2006)

# should loop over various values of alpha and compare the MSE's to determine the optimal elastic net

# NOW DO THE SAME FOR DIVISIONS 2 THROUGH 9:

y = as.matrix(DIV[,2])
weights <- adaptive.weights(x, y)
cv.mod2_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor = weights$weights)
plot(cv.mod2_mae)
(minMAE2 = min(cv.mod2_mae$cvm)) # This is much lower, implying that there are a lot of outliers that it does not predict well that are getting overweighted when using RMSE
(absRev2 = mean(abs(rev1div2)))

y = as.matrix(DIV[,3])
weights <- adaptive.weights(x, y)
cv.mod3_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor = weights$weights)
plot(cv.mod3_mae)
(minMAE3 = min(cv.mod3_mae$cvm)) # This is much lower, implying that there are a lot of outliers that it does not predict well that are getting overweighted when using RMSE
(absRev3 = mean(abs(rev1div3)))

y = as.matrix(DIV[,4])
weights <- adaptive.weights(x, y)
cv.mod4_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor = weights$weights)
plot(cv.mod4_mae)
(minMAE4 = min(cv.mod4_mae$cvm)) # This is much lower, implying that there are a lot of outliers that it does not predict well that are getting overweighted when using RMSE
(absRev4 = mean(abs(rev1div4)))

y = as.matrix(DIV[,5])
weights <- adaptive.weights(x, y)
cv.mod5_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor = weights$weights)
plot(cv.mod5_mae)
(minMAE5 = min(cv.mod5_mae$cvm)) # This is much lower, implying that there are a lot of outliers that it does not predict well that are getting overweighted when using RMSE
(absRev5 = mean(abs(rev1div5)))

y = as.matrix(DIV[,6])
weights <- adaptive.weights(x, y)
cv.mod6_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor = weights$weights)
plot(cv.mod6_mae)
(minMAE6 = min(cv.mod6_mae$cvm)) # This is much lower, implying that there are a lot of outliers that it does not predict well that are getting overweighted when using RMSE
(absRev6 = mean(abs(rev1div6)))

y = as.matrix(DIV[,7])
weights <- adaptive.weights(x, y)
cv.mod7_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor = weights$weights)
plot(cv.mod7_mae)
(minMAE7 = min(cv.mod7_mae$cvm)) # This is much lower, implying that there are a lot of outliers that it does not predict well that are getting overweighted when using RMSE
(absRev7 = mean(abs(rev1div7)))

y = as.matrix(DIV[,8])
weights <- adaptive.weights(x, y)
cv.mod8_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor = weights$weights)
plot(cv.mod8_mae)
(minMAE8 = min(cv.mod8_mae$cvm)) # This is much lower, implying that there are a lot of outliers that it does not predict well that are getting overweighted when using RMSE
(absRev8 = mean(abs(rev1div8)))

y = as.matrix(DIV[,9])
weights <- adaptive.weights(x, y)
cv.mod9_mae = cv.glmnet(x, y, alpha=1, nfolds=nrow(DIV), type.measure = "mae", penalty.factor = weights$weights)
plot(cv.mod9_mae)
(minMAE9 = min(cv.mod9_mae$cvm)) 
(absRev9 = mean(abs(rev1div9)))

# Mention in paper that the Bayesian elastic net does not offer much improvement in prediction, and that cross-validation is not well grounded from a Bayesian perspective.

# First run adalasso to make sure it works, then run adalassomae.R (modified version using MAE cost function)
X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,1])
lassoFit1 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit1$cv.adalasso
(absRev1 = mean(abs(rev0div1))) # For comparison

X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,2])
lassoFit2 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit2$cv.adalasso
(absRev2 = mean(abs(rev0div2))) # For comparison

X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,3])
lassoFit3 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit3$cv.adalasso
(absRev3 = mean(abs(rev0div3))) # For comparison

X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,4])
lassoFit4 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit4$cv.adalasso
(absRev4 = mean(abs(rev0div4))) # For comparison

X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,5])
lassoFit5 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit5$cv.adalasso
(absRev5 = mean(abs(rev0div5))) # For comparison

X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,6])
lassoFit6 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit6$cv.adalasso
(absRev6 = mean(abs(rev0div6))) # For comparison

X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,7])
lassoFit7 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit7$cv.adalasso
(absRev7 = mean(abs(rev0div7))) # For comparison

X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,8])
lassoFit8 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit8$cv.adalasso
(absRev8 = mean(abs(rev0div8))) # For comparison

X = as.matrix(DIV[,10:ncol(DIV)])
y = as.matrix(DIV[,9])
lassoFit9 = adalassomae(X, y, k=nrow(DIV)) 
lassoFit9$cv.adalasso
(absRev9 = mean(abs(rev0div9))) # For comparison


# Plot the most recent release of HPI from 4q09 through 2q15 to see whether we are analyzing both upward and downward trends in the market (which means we don't have to worry as much about overfitting):


# Now adjust the first revision of the most recent (2q15) release based on the prediction, and plot the revised series on top of the original series to see the difference:


# Now do this with the states. We would have even more information so supposedly the prediction will be even more significant.

# Now do it Bayesian (first look at papers to see whether Bayesian prediction performs better in practice):

# bind the vectors together
REV1 = cbind(rev1div1, rev1div2, rev1div3, rev1div4, rev1div5, rev1div6, rev1div7, rev1div8, rev1div9)
APP1 = cbind(app1div1, app1div2, app1div3, app1div4, app1div5, app1div6, app1div7, app1div8, app1div9)
DATA = cbind(REV1, APP1)

# Run the multivariate regression:
REG = Mlm(REV1, APP1, output=FALSE) # mvreg REV1 APP1
summary(REG)

# Do prediction (will need to write the predict function since it isn't in the MTS package)
coef = REG$beta
se = REG$se.beta
newdata = as.matrix(APP1[4, , drop=FALSE]) # took random row of the old data (vector of first appreciations) to use for prediction

# Each row is a division. 1st column is the prediction, 2nd and 3rd columns are the lower and upper bounds of the 95% confidence interval
pred = matrix(, nrow=9, ncol=3) 
# coef[2,1]*newdata[,1] + coef[3,1]*newdata[,2] + ... + coef[10,1]*newdata[,9]
for (i in 1:9) {
  tmp = coef[2:10,i]
  pred[i,1] = coef[1,i] + tmp %*% newdata[,1:9]
  
  # Need to do use a different formula to calculate the "prediction interval" (assuming that the regression errors are normally distributed)
  # There is the question of whether to use "prediction" or "confidence" intervals. 
  # http://robjhyndman.com/hyndsight/intervals/
  # http://stats.stackexchange.com/questions/9131/obtaining-a-formula-for-prediction-limits-in-a-linear-model
  tmp = coef[2:10,i] - 1.96 * se[2:10,i]
  pred[i,2] = (coef[1,i] - 1.96 * se[1,i]) + tmp %*% newdata[,1:9]
  tmp = coef[2:10,i] + 1.96 * se[2:10,i]
  pred[i,3] = (coef[1,i] + 1.96 * se[1,i]) + tmp %*% newdata[,1:9]
}
  


# Now run a Bayesian multivariate linear regression, since this approach is theoretically better:



### Now, focus on predicting the first revision using lagged revisions with AR(1) state-space model. ###

# Now it is stationary. We need this to guarantee that the Kalman Filter converges to the Matrix-Ricatti equation. 
# See Quant Econ page about Kalman filter for more details.

# Start with the following parsimonious model (to avoid overfitting):
# rev1_t = c + app1_{t} + app1_{t-1} + eps_t
# (In this case, it is just OLS, so don't need to use Kalman filter to get optimal forecast.)

# Now do this for every state to make sure the above parsimonious model makes sense:



# Also, use BIC to determine the optimal number of AR lags. Take average across states.
auto.arima(app1,ic="bic")

# Need to also make sure the series of HPI appreciation is stationary? 

# Check whether the state series are correlated. If so, then run a VAR to do prediction.


### Time-Domain MLE ###

# Model:
# Regular Representation:
# 
# State Space Representation:
# y_t = theta_t
# theta_t = mu*(1-rho) + rho*theta_{t-1} + eps_t, eps_t ~ N(0,sigma^2)

# Legend:
# x[1] = mu
# x[2] = rho
# x[3] = sigma
# Note that the model is time-varying
buildMod = function(x) {
  print(x[2])
  # need to ensure that W and C0 do not become negative
  mod = dlm(FF = as.matrix(c(1,app1), nrow=1, ncol=2), V=as.matrix(0), GG=as.matrix(x[2]), W=as.matrix(x[3]^2), m0=as.matrix(x[1]), C0=as.matrix(abs(x[3]^2/(1-x[2]^2))))
  return(mod)
}

# regression to determine initial parameter guesses:
regPrior = ar(diff2Rev1, order.max=1, method="ols")
print(regPrior)

# perform the maximum likelihood estimation
# L-BFGS-B is gradient-based
fitML = dlmMLE(diff2Rev1, 
                parm=c(regPrior$x.mean, 0.5, sqrt(regPrior$var.pred[1])), 
                build=buildMod, hessian=TRUE, method="L-BFGS-B", 
                control = list(maxit = 500));

modML = buildMod(fitML$par)

drop(GG(modML))
drop(W(modML))

# optim minimizes the negative log-likelihood, so it returns the negative Hessian.
# The covariance matrix is the inverse of the population expectation of the 
# negative Hessian (i.e., the inverse of the information matrix).
sqrt(diag(solve(fitML$hessian)))

# Need to run Kalman smoother first before forecasting?

forecast = dlmForecast(modML, nAhead=4)
predictions = forecast$f # forecasts
predVar = forecast$Q # variances of forecasts
# The above forecasts are not statistically significant. Need to make the model richer to fix this. Start by including all 4 AR lags.

# Linear regression with time-varying coefficients + time-varying constant 
# http://lalas.github.io/quantitativeThoughts/r/2014/09/01/dlmTutorial.html
# Will combine the two models later.


orig = diffinv(diff2Rev1, xi=rev1[1:2], differences = 2) # use this to recover the original series
# Note that when we integrate the series of forecasted revisions, we do not know the initial two integrated revision forecasts, so we have to just assume they are 0

# Remember to do cross-validation and diagnostic testing (Hamilton Ch. 5 and 14)
# http://www.rinfinance.com/agenda/2012/workshop/Zivot+Yollin.pdf

# Next, Focus on predicting the first revision using lagged house price appreciation.
# Just use OLS for this?


