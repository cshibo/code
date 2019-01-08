# CS544 Final Project
# Shibo Cao

# dataset: Ramen rating


### Part 1: Data processing
# import statistical packages
library(dplyr)
library(sampling)
library(stringr)

# set working directory
setwd("~/Desktop/544/CS544Final_Cao")

# read file
ramen_data <- read.csv("ramen-ratings.csv")

# eliminate all data with no information for Stars and Style
ramen_data <- ramen_data[ramen_data$Stars != "Unrated",]
ramen_data <- ramen_data[ramen_data$Style != "",]

# set all Star rating with two decimal places
ramen_data$Stars <- as.character(ramen_data$Stars)
for (i in 1: nrow(ramen_data)) {
  if (nchar(ramen_data$Stars[i]) == 1){
    ramen_data$Stars[i] <- paste(ramen_data$Stars[i], ".", "00", sep = "")
  }
  if (nchar(ramen_data$Stars[i]) == 3){
    ramen_data$Stars[i] <- paste(ramen_data$Star[i], "0", sep = "")
  }
}
ramen_data$Stars <- as.numeric(ramen_data$Stars)
# In "Conutry" column, convert all "United States" to "USA"
row_USA_num <- which(ramen_data$Country == "United States")
row_USA <- ramen_data[row_USA_num,]
row_USA$Country <- "USA"
ramen_data[row_USA_num,] <- row_USA

### Part 2A: Categorical data analysis: origin country of product
# frequency table:
ramen_freq <- table(ramen_data$Country)

# sort out the table based on frequency with descending order
ramen_freq_ordered <- ramen_data %>% 
  count(Country) %>%
  arrange(desc(n))
ramen_freq_ordered <- as.data.frame(ramen_freq_ordered)

# Combine those categories with a frequency less than 100 into "Others"
# locate the index and calculate the sum
less_than_hundred_index <- which(ramen_freq_ordered$n <= 100)
total_freq_others <- sum(ramen_freq_ordered$n[less_than_hundred_index])
ramen_freq_ordered <- subset(ramen_freq_ordered, n > 100)
Others_row <- data.frame(Country = "Others", n = total_freq_others)
ramen_freq_ordered <- rbind(ramen_freq_ordered, Others_row)

# bar plot:
barplot(as.vector(ramen_freq_ordered$n), horiz = TRUE, xlim = c(0,400), col = rainbow(nrow(ramen_freq_ordered)), xlab = "Frequency", names.arg = ramen_freq_ordered$Country, las = 1, cex.names=0.6)



### Part 2B: Numerical data analysis: Stars rating
# mean, median and mode of stars rating across all data
mean_rating <- mean(ramen_data$Stars, trim = 0.1)
median_rating <- median(ramen_data$Stars)
mode_rating <- which(table(ramen_data$Stars) == max(table(ramen_data$Stars)))

# distribution of data
range_rating <- range(ramen_data$Stars)
variance_rating <- var(ramen_data$Stars)
stdev_rating <- sd(ramen_data$Stars)
f <- fivenum(ramen_data$Stars)
Quantile_range <- quantile(ramen_data$Stars, c(0, 0.25, 0.5, 0.75, 1))
Interquar_range <- IQR(ramen_data$Stars)

# histogram and boxplot of Stars rating
attach(ramen_data)
hist_stars <- hist(Stars, col = hcl(0), xlim = c(0,5), ylim = c(0, 800))
boxplot(Stars, horizontal = TRUE, xaxt = "n")
axis(side = 1, at= fivenum(Stars), labels = TRUE)

# outliers
outlier_range <- c(f[2] - 1.5*(f[4] - f[2]))
outlier_num <- Stars[Stars <= outlier_range]



### Part 3: Bivariate analysis: Country and Style
country_style_table <- rbind(table(Country, Style))

# distribution of percentage of each style of ramen different countries produce
country_style_prop_table <- prop.table(country_style_table,1)
barplot(t(country_style_table), horiz = TRUE, col = rainbow(ncol(country_style_table)), xlab = "Frequency", xlim = c(0, 450), ylim = c(0,38), las = 1, cex.names=0.5, legend.text=TRUE)



### Part 4: Distribution of numerical data
# assume the Star rating follows normal distribution
# PDF
x <- seq(0,5, by = 0.1)
pdf <- dnorm(x, mean = mean_rating, sd = stdev_rating)
plot(x, pdf, type = "l", col = "red", xlim = c(0,5), ylim = c(0,0.5), main = "Star Rating PDF" )

# CDF
x <- seq(0,5, by = 0.1)
cdf <- pnorm(x, mean = mean_rating, sd = stdev_rating)
plot(x, cdf, type = "l", col = "red", xlim = c(0,5), ylim = c(0,1), main = "Star Rating CDF" )
abline(h = 0)



### Part 5: Central Limit theorem
# if sample size = 50, sampling 1000 times
samples <- 1000
sample.size <- 50

xbar <- numeric(samples)
for (i in 1: samples){
  xbar[i] <- mean(rnorm(sample.size, mean = mean_rating, sd = stdev_rating))
}

hist(xbar, prob = TRUE, xlim = c(3.25,4.25), breaks = 10, ylim = c(0,6), xlab = "Star rating", main = "Sample size = 50 / 1000 times")

# if sample size = 100, sampling 1000 times
samples <- 1000
sample.size <- 100

xbar <- numeric(samples)
for (i in 1: samples){
  xbar[i] <- mean(rnorm(sample.size, mean = mean_rating, sd = stdev_rating))
}

hist(xbar, prob = TRUE, xlim = c(3.25,4.25), breaks = 10, ylim = c(0,6), xlab = "Star rating", main = "Sample size = 100 / 1000 times")

# if sample size = 200, sampling 1000 times
samples <- 1000
sample.size <- 200

xbar <- numeric(samples)
for (i in 1: samples){
  xbar[i] <- mean(rnorm(sample.size, mean = mean_rating, sd = stdev_rating))
}

hist(xbar, prob = TRUE, xlim = c(3.25,4.25), breaks = 10, ylim = c(0,6), xlab = "Star rating", main = "Sample size = 200 / 1000 times")



### Part 6: Sampling method
# Simple random sampling (100 sample without replacement)
s <- srswor(100, nrow(ramen_data))
sample_1 <- ramen_data[s != 0,]

# Systematic sampling (100)
N <- nrow(ramen_data)
n <- 100
k <- ceiling(N / n)
r <- sample(k,1)
s <- seq(r, by = k, length = n)
sample_2 <- ramen_data[s,]

# Unequal stratified sampling based on Style
# order data
index <- order(ramen_data$Style)
style_ordered_data <- ramen_data[index,]
# Stratified sample
freq_style <- table(style_ordered_data$Style)
str_size <- round(100 * freq_style / sum(freq_style))
# Based on the str_size, we do not draw any sample from some of the categories
# therefore we drop them by create a new subset
subset <- ramen_data[!ramen_data$Style %in% c("Bar", "Box", "Can"), ]
subset <- subset[order(subset$Style),]
freq_style_subset <- table(subset$Style)
freq_style_subset <- freq_style_subset[freq_style_subset != 0]
str_size_subset <- str_size[str_size != 0]
sample_3 <- getdata(style_ordered_data,strata(subset, stratanames = c("Style"), size = str_size_subset, method = "srswor", description = TRUE))

# Cluster sampling
# Choose 3 out of 6 styles
cl <- cluster(ramen_data, c("Style"), size = 3, method = "srswor")
first_cluster <- getdata(ramen_data, cl)
# 2nd step sampling: simple random sampling of 100 samples
s <- srswor(100, nrow(first_cluster))
sample_4 <- first_cluster[s != 0,]

# Compare mean & sd of each sample
mean(sample_1$Stars, na.rm = TRUE)
sd(sample_1$Stars, na.rm = TRUE)
mean(sample_2$Stars, na.rm = TRUE)
sd(sample_2$Stars, na.rm = TRUE)
mean(sample_3$Stars, na.rm = TRUE)
sd(sample_3$Stars, na.rm = TRUE)
mean(sample_4$Stars, na.rm = TRUE)
sd(sample_4$Stars, na.rm = TRUE)

### Part 7: Implementing features not mentioned in the specification
# create a subset with data of top_ten rating 2012-2016
all_top_ten <- ramen_data[(ramen_data$Top.Ten != "" & ramen_data$Top.Ten != "\n"),]
all_top_ten <- all_top_ten[order(all_top_ten$Top.Ten),]
# add year to the subset
yr <- rep(c(2012,2013,2014,2015,2016), c(9,7,8,7,6))
all_top_ten <- cbind(all_top_ten, yr)
# make the subset as a tibble
all_top_ten <- as_tibble(all_top_ten)

# 2012-2016 Top 10 ramen country of origin
yr_country_count <- all_top_ten %>% group_by(yr,Country) %>% summarise(count = n())
sum_2012 <- filter(yr_country_count, yr == 2012)
arrange(sum_2012, desc(count))
sum_2013 <- filter(yr_country_count, yr == 2013)
arrange(sum_2013, desc(count))
sum_2014 <- filter(yr_country_count, yr == 2014)
arrange(sum_2014, desc(count))
sum_2015 <- filter(yr_country_count, yr == 2015)
arrange(sum_2015, desc(count))
sum_2016 <- filter(yr_country_count, yr == 2016)
arrange(sum_2016, desc(count))

# all country of origin
Country_count <- all_top_ten %>% group_by(Country) %>% summarise(count = n())
arrange(Country_count, desc(count))


