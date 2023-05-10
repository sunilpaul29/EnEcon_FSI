# Load the packaes required for replicating the figures
library("tidyverse")
library("psych")
library("Hmisc")
library("dplyr")

# load the data
FSIdata <- read.csv("FSI_data.csv", header=TRUE)
FSI_data1 <- as.matrix(FSIdata[,c(2,3,4,5,6)])
# function to calculate population stddev
sd.p=function(x){sd(x)*sqrt((length(x)-1)/length(x))}

# formula to calculate pearson's kurtosis
kurtosis_Pearson <- function(x) {
  n <- length(x)
  m <- mean(x)
  v <- sd.p(x)
 mean(((x - m)/v)^4)
}

####  Table 1  ######
# following code replicates Table 1 in the paper. The significance level is marked with in paper usisnf *,**,*** based on the p values
p_values <- rcorr(as.matrix(FSI_data1))
print("Table 1: Cross correlations between financial market sub-indices")
p_values

#### Table 2 #####
#codes to replicate Table 2
df_summary <- psych::describe(FSIdata)
df_summary <- as.data.frame(df_summary)
df_summary <- df_summary %>% select(median, max, min, sd, skew, kurtosis)
Table2 <- t(df_summary[c(-1,-7,-8),])
Table2[6,] <- apply(FSI_data1, 2, kurtosis_Pearson)
Table2 <- round(Table2,3)
rownames(Table2)<- c("Median", "Maximum","Minimum", "Std.Dev.","Skewness", "Kurtosis")
print("Table 2: Descriptive statistics")

Table2
