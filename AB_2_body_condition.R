## Body condition analysis for WCS territory size study
## Date: June 2022
## Author: Ruth Simberloff (rsimberloff@gmail.com)


# Method
# Scaled Mass Index of body condition, following Peig & Green 2009
# This index is based on residuals from an OLS regression of mass (y) and wing length (x)

# Scaled Mass Index = (body mass measurement)(mean population length / length measurement)^(bSMA)
# bSMA = (slope from OLS regression / Pearson's correlation coefficient r)


### Set up
library(tidyverse)


### Import data
sparrows <- read_csv("condition_data.csv")
glimpse(sparrows)

### visualize how mass scales with wing length

ggplot(sparrows, aes(x = wing, y = weight)) +
  geom_point()


### Calculate bSMA
###### OLS regression
OLSregression <- lm(weight ~ wing, data = sparrows)
OLSslope <- OLSregression$coefficients[2]

###### Pearson's correlation coefficient
cor.coef <- cor(sparrows$weight, sparrows$wing)

###### bSMA calculation
bSMA <- OLSslope/cor.coef

### Calculate L0 (mean population length)
L0 <- mean(sparrows$wing)


### Scaled Mass Index

###### add a "condition" column
sparrows <- sparrows %>%
  mutate(condition = (weight)*(L0/wing)^bSMA)

glimpse(sparrows)


###### visualize
ggplot(sparrows, aes(x = wing, y = weight, color = condition, shape = habitat)) +
  geom_point()
