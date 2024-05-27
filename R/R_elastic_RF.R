# Install and load necessary packages
if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret")
}
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}

if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}
library(caret)
library(glmnet)
library(dplyr)

library(tidyr)

###################

library(randomForest)
library(datasets)
library(caret)

file_path <- "C:/Users/46705/Documents/SpiderSilk/code/random_forest/test_in_r.csv"

# Read the CSV file into a data frame
data <- read.csv(file_path)


# Split the data into predictors (X) and target variable (y)
X <- subset(data, select = -tensile_strength)
y <- data$tensile_strength

set.seed(123)  # Set seed for reproducibility
# Generate random indices for sampling
sample_indices <- sample(c(TRUE, FALSE), nrow(data), replace = TRUE, prob = c(0.8, 0.2))

# Subset the dataframe into training and testing datasets
train <- data[sample_indices, ]
test <- data[!sample_indices, ]



print(head(data))

ts.rf <- randomForest(tensile_strength ~ ., data = train, mtry = 3, 
                         importance = TRUE, na.action = na.omit)

print(ts.rf)

plot(ts.rf)

# Predict the tensile_strength values for the test dataset
predicted_values <- predict(ts.rf, newdata = test)

# Extract the actual tensile_strength values from the test dataset
actual_values <- test$tensile_strength

# Calculate R^2 value
r_squared <- cor(predicted_values, actual_values)^2

# Plot predicted vs. actual values
plot(actual_values, predicted_values, xlab = "Actual Tensile Strength", ylab = "Predicted Tensile Strength", main = "Predicted vs. Actual Tensile Strength")

# Add a diagonal line for reference
abline(0, 1, col = "red")

# Print R^2 value
cat("R-squared value:", r_squared, "\n")

###########################






# Specify the local path to your CSV file
file_path <- "C:/Users/46705/Documents/SpiderSilk/generalized/this.csv"

# Read the CSV file into a data frame
data <- read.csv(file_path)


# Split the data into predictors (X) and target variable (y)
X <- subset(data, select = -tensile_strength)
y <- data$tensile_strength

set.seed(42)
cv_5 = trainControl(method = "cv", number = 5)

hit_elnet = train(
  y ~ ., data = X,
  method = "glmnet",
  trControl = cv_5
)
