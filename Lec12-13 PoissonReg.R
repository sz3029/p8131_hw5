### Wave Damage Data (Page 205 of GLM)
#The data concern a type of damage caused by waves to the forward section of certain cargo-carrying vessels. For the purpose of setting standards for hull construction, one needs to know the risk of damage associated with the three classifying factors shown below:
#
#           Ship Type: A-E
#Year of Construction: 1960-64, 65-69, 70-74, 75-79
# Period of operation: 1960-74, 75-79
#
#Two other variables are Aggregated months of service and Number of damage accidents.
#
#The question of interest is how the number of damage accidents depends on the other four variables.

## formulate the data into a data.frame
month=c(127,63, 1095,1095,1512, 3353, 0,2244, 44882, 17176, 28609, 20370,  7064, 13099,0, 7117,1179,552,781,676,783,1948,0,274,251,105,
        288,192,349,1208, 0,2051,45,0,789,437,1157, 2161, 0,542)
damage=c(0 , 0  ,3  ,4  ,6, 18,  0, 11, 39, 29, 58, 53, 12, 44,  0, 18,  1,  1,  0,  1,  6,  2,  0,  1,  0,  0,  0,  0,
         2, 11,  0,  4,  0,  0,  7,  7,  5, 12,  0,  1)
ship <- rep(LETTERS[1:5], rep(8,5)) # 5 groups
year <- rep(rep(c("60-64", "65-69", "70-74","75-79"), rep(2,4)), 5) # 4 groups
period <- rep(c("60-74", "75-79"), 20) # 2 groups
wave <- data.frame(ship, year, period, month, damage)

#   ship  year period month damage
#1     A 60-64  60-74   127      0
#2     A 60-64  75-79    63      0
#3     A 65-69  60-74  1095      3
#4     A 65-69  75-79  1095      4
#5     A 70-74  60-74  1512      6
#...

# fit Poisson log linear model
wave.glm1 <- glm(damage~ship+year+period+offset(log(month)), family=poisson, data=wave, subset=month>0)
summary(wave.glm1) # n=34, p=5-1 + 4-1 + 2-1 + 1=9

### check interactions 
wave.glm2 <- glm(damage~ship*year+period+offset(log(month)), family=poisson(link=log), data=wave, subset=month>0)
summary(wave.glm2) # how to calc df?  


## deviance analysis (ignoring the over dispersion)
test.stat=wave.glm1$deviance-wave.glm2$deviance
df=25-13
pval=1-pchisq(test.stat,df=df) # chisq test
pval # rej, go with the bigger model






##############################################################################
### estimate the dispersion parameter (from the additive model)
# the traditional way of calc constant dispersion parameter
res.p1=residuals(wave.glm1,type='pearson',data=wave, subset=month>0)  # exactly the same as pearson residual for wave.glm3
G1=sum(res.p1^2) # calc dispersion param based on full model
pval=1-pchisq(G1,df=25) # lack of fit
phi=G1/(34-9)
phi # 1.69
wave.glm1$deviance/wave.glm1$df.residual # 1.55

summary(wave.glm1,dispersion=phi)
# an equivalent way of estimating dispersion parameter and fit over-dispersed poisson regression
wave.glm3 <- glm(damage~ship+year+period+offset(log(month)), family=quasi(link=log,variance=mu), data=wave, subset=month>0)
summary(wave.glm3) # check ?quasi 


# test over-dispersion (half normal plot)
plot(qnorm((34+1:34+0.5)/(2*34+1.125)),sort(abs(res.p1)),xlab='Expected Half-Normal Order Stats',ylab='Ordered Abs Pearson Residuals')
abline(a=0,b=1)
abline(a=0,b=sqrt(phi),lty=2)  # controversial?



# deviance analysis
test.stat=wave.glm1$deviance-wave.glm2$deviance # deviance (from original model fitting)
df=25-13
res.p=residuals(wave.glm2,type='pearson')  
res.p 
G=sum(res.p^2) # calc dispersion param based on larger model
phi=G/13
F.stat=test.stat/(df*phi)
pval=1-pf(F.stat,df,13)
pval # not rej, go with the smaller model


######################################
library(MASS)
newwave=wave[wave$month>0,]
# negative binomial regression
wave.nb=glm.nb(damage~ship+year+period+offset(log(month)),data=newwave)
summary(wave.nb) # dispersion param theta will be estimated simultaneously 
# the result (huge theta value) indicates no over-dispersion 



##################################################
# Analyze fish data for Zero-inflated Poisson regression
# install.packages("pscl")

zinb <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv")
zinb=zinb[,c(3,4,5,8)]
head(zinb)

library(pscl)
m1 <- zeroinfl(count ~ child + camper | persons, data = zinb) # child and camper for poisson, persons for binary (event is true zero, aka not fishing)
summary(m1)
# interpretation:
# more people in the group, more likely to fish
# more children, fewer fish; those who camp tend to catch more fish

pr <- predict(m1,type="zero")  # pi(do not fish)
mu <- predict(m1,type="count") # given fish, how many caught
zip <- pr + (1-pr)*exp(-mu) # total prob of catch 0 fish 
cbind(zip,zinb$count)
