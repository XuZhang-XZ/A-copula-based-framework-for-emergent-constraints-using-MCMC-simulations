
empirical_bivariate = function(x, y) {
  probss = 0
  for (i in 1:length(x)) {
    probss[i] = length(which(x <= x[i] & y <= y[i]))
  }
  probss =  probss / (length(x) + 0)
  return(probss)
  
}


# x <- c(1,3,2,2,8,2,1,3,1,1,3,3,1,1,2,1,2,1,1,3,4,1,1,3,1,1,2,1,3,7,1,4,6,1,2,1,1,3,1,2,2,3,4,1,1,1,1,2,2,12,1,1,2,1,1,1,3,4)
# y <- c(1.42,5.15,2.52,2.29,12.36,2.82,1.49,3.53,1.17,1.03,4.03,5.26,1.65,1.41,3.75,1.09,3.44,1.36,1.19,4.76,5.58,1.23,2.29,7.71,1.12,1.26,2.78,1.13,3.87,15.43,1.19,4.95,7.69,1.17,3.27,1.44,1.05,3.94,1.58,2.29,2.73,3.75,6.80,1.16,1.01,1.00,1.02,2.32,2.86,22.90,1.42,1.10,2.78,1.23,1.61,1.33,3.53,10.44)
# 
# library(mltools)
# library(data.table)
# 
# # set data in a data.table
# dt <- data.table(x = x, y = y)
# 
# empirical_cdf(dt, ubounds = data.table(x = 3, y = 5))
# mean(x <= 3 & y <= 5) # same result

# qqnorm(c(55,57,51,96,88))
# 
# y <- rt(200, df = 5)
# qqnorm(y); qqline(y, col = 2)
# qqplot(y, rt(300, df = 5))
# 
# ?rt
# 
# 
# set.seed(1)
# x <- rnorm(10)
# y <- rnorm(500, mean = 2.5, sd = .95)
# ex <- TRUE
# ### p = 0.112
# (pval <- ks.test(x, y, exact = ex)$p.value)
# ## 88.8% confidence band with bisecting line
# ## touching the lower bound
# qqplot(x, y, pch = 19, conf.level = 1 - pval, 
#        conf.args = list(exact = ex, col = "lightgrey"))
# abline(a = 0, b = 1)
# 
