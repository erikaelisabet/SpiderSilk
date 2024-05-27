library(caret)
library(glmnet)

file_path <- "C:/Users/46705/Documents/SpiderSilk/data/r_analysis/this.csv"
file_path <- "C:/Users/46705/Documents/SpiderSilk/data/post_filtering/filtered_full_ara_log.csv"

# Load the CSV file
data <- read.csv(file_path)
data

# Remove first 3 rows
data_ <- data[-c(1:3), ]
# Create a vector of indices for the columns to keep
selected_columns <- c(8:(ncol(data_) - 29), which(names(data_) == "tensile_strength"))
# Subset the dataframe with selected columns
df_selected <- data_[, selected_columns, drop = FALSE]
df_selected[is.na(df_selected)] <- 0



X <- df_selected[, !names(df_selected) %in% c("tensile_strength", "tensile_strength.1")]
y <- df_selected[["tensile_strength"]]



fit <- glmnet(X, y, alpha = 0.5)

plot(fit)

print(fit)

coef(fit, s = 0.1)




# Set seed for reproducibility
set.seed(150)

# Pre-define the folds
foldid <- createFolds(y, k = 5)

# Create a list to store cross-validation results for different alpha values
cv_results <- list()

# Define alpha values to test
alpha_values <- c(0.1, 0.25, 0.5, 0.75, 0.9)

# Iterate over alpha values
for (alpha in alpha_values) {
  # Perform cross-validation with cv.glmnet using pre-defined folds and current alpha
  cvfit <- cv.glmnet(x = X, y = y, alpha = alpha, foldid = unlist(foldid))
  
  # Store the cross-validation results
  cv_results[[as.character(alpha)]] <- cvfit
}

# Print or access cross-validation results for each alpha value
for (alpha in alpha_values) {
  cat("Alpha:", alpha, "\n")
  print(cv_results[[as.character(alpha)]])
}



#####################################
# Create an empty plot with the desired y-axis limits
# Create an empty plot with the desired y-axis limits
plot(NULL, xlim = range(unlist(lapply(cv_results, function(cvfit) cvfit$lambda))), 
     ylim = c(y_min, y_max), xlab = "Lambda", ylab = "MSE", 
     main = "Cross-Validation Results", type = "n")

# Set up colors for lines
line_colors <- rainbow(length(alpha_values))

# Print or access cross-validation results for each alpha value
for (i in seq_along(alpha_values)) {
  alpha <- alpha_values[i]
  cvfit <- cv_results[[as.character(alpha)]]
  cvm_values <- cvfit$cvm
  lambda_values <- cvfit$lambda
  
  # Plot cross-validation results
  lines(lambda_values, cvm_values, col = line_colors[i], 
        main = paste("Alpha =", alpha), type = "l")
  
  # Get the minimum value of the y-axis
  min_y <- min(cvm_values)
  
  # Add a point indicating the minimum point
  points(lambda_values[which.min(cvm_values)], min_y, col = "black", pch = 19)
  
  # Add text indicating the minimum value of the y-axis
  text(lambda_values[which.min(cvm_values)], min_y, paste("", round(min_y, digits = 3)), col = "black", pos = 3, cex = 0.8)
}

# Add a legend
legend("topright", legend = alpha_values, col = line_colors, lty = 1, title = "Alpha")

# Get the index of the alpha value with the lowest minimum error
min_alpha_index <- which.min(sapply(cv_results, function(cvfit) min(cvfit$cvm)))
min_alpha <- alpha_values[min_alpha_index]
min_cvfit <- cv_results[[as.character(min_alpha)]]

# Get the lambda value with the lowest MSE
min_lambda <- min_cvfit$lambda[which.min(min_cvfit$cvm)]
min_mse <- min(min_cvfit$cvm)

# Add vertical line from the lowest MSE to the x-axis
abline(v = min_lambda, col = "red", lty = 2)

# Add horizontal line from the lowest lambda to the y-axis
abline(h = min_mse, col = "red", lty = 2)

# Add point indicating the minimum MSE
points(min_lambda, min_mse, col = "red", pch = 19)

# Add text indicating the minimum MSE and corresponding lambda
text(min_lambda, min_mse, paste("MSE =", round(min_mse, digits = 3), "\nLambda =", round(min_lambda, digits = 3)), col = "red", pos = 1, cex = 0.8)


# Add a legend
legend("topright", legend = alpha_values, col = line_colors, lty = 1, title = "Alpha")

# Get the index of the alpha value with the lowest minimum error
min_alpha_index <- which.min(sapply(cv_results, function(cvfit) min(cvfit$cvm)))
min_alpha <- alpha_values[min_alpha_index]
min_cvfit <- cv_results[[as.character(min_alpha)]]

# Get the lambda value with the lowest MSE
min_lambda <- min_cvfit$lambda[which.min(min_cvfit$cvm)]
min_mse <- min(min_cvfit$cvm)

# Add vertical line from the lowest MSE to the x-axis
abline(v = min_lambda, col = "red", lty = 2)

# Add horizontal line from the lowest lambda to the y-axis
abline(h = min_mse, col = "red", lty = 2)

# Add point indicating the minimum MSE
points(min_lambda, min_mse, col = "red", pch = 19)

# Add text indicating the minimum MSE and corresponding lambda
text(min_lambda, min_mse, paste("MSE =", round(min_mse, digits = 3), "\nLambda =", round(min_lambda, digits = 3)), col = "red", pos = 1, cex = 0.8)
##########################

plot(cv_results[[as.character(0.5)]], xlab = "Lambda")

###################################################


#######################################################
cvfit <- cv.glmnet(X, y, alpha = 0.5)
plot(cvfit)
cvfit$lambda.min
coefficients <- coef(cvfit, s = "lambda.min")

plot(cvfit)


# Select non-zero coefficients
non_zero_indices <- which(coefficients != 0)
non_zero_coefficients <- coefficients[non_zero_indices]

# Print non-zero coefficients along with their corresponding feature names
for (i in 1:length(non_zero_indices)) {
  feature_index <- non_zero_indices[i]
  feature_name <- feature_names[feature_index]
  coefficient <- non_zero_coefficients[i]
  cat("Feature:", feature_name, "Coefficient:", coefficient, "\n")
} 

# Initialize an empty dataframe to store non-zero coefficients and feature names
non_zero_df <- data.frame(
  Feature = character(),
  Coefficient = numeric(),
  stringsAsFactors = FALSE
)

# Iterate over non-zero indices and populate the dataframe
for (i in 1:length(non_zero_indices)) {
  feature_index <- non_zero_indices[i]
  feature_name <- feature_names[feature_index]
  coefficient <- non_zero_coefficients[i]
  
  # Append to the dataframe
  non_zero_df <- rbind(non_zero_df, data.frame(Feature = feature_name, Coefficient = coefficient))
}

# Save the dataframe as a CSV file
write.csv(non_zero_df, file = "C:/Users/46705/Documents/SpiderSilk/data/r_analysis/EN_non_zero_features.csv", row.names = FALSE)


######################################################3


set.seed(42)
cv_5 = trainControl(method = "cv", number = 5)


hit_elnet <- train(
  tensile_strength ~ ., data = df_selected,  # Corrected formula
  method = "glmnet",
  trControl = cv_5, 
  tuneLength = 15
  
)

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

get_best_result(hit_elnet)

fit <- cv.glmnet(as.matrix(X), y, alpha = 0.5)

# Print the cross-validated mean squared error for each lambda
print(fit$cvm)

# Choose the lambda that minimizes mean squared error
best_lambda <- fit$lambda.min

# Fit the final elastic net model using the selected lambda
final_model <- glmnet(as.matrix(X), y, alpha = 0.5, lambda = best_lambda)

# Print the coefficients of the final model
print(final_model$beta)
# Get the coefficients of the final model
coefficients <- final_model$beta

# Get the names of the predictor variables (features)
feature_names <- colnames(X)

# Select non-zero coefficients
non_zero_indices <- which(coefficients != 0)
non_zero_coefficients <- coefficients[non_zero_indices]

# Print non-zero coefficients along with their corresponding feature names
for (i in 1:length(non_zero_indices)) {
  feature_index <- non_zero_indices[i]
  feature_name <- feature_names[feature_index]
  coefficient <- non_zero_coefficients[i]
  cat("Feature:", feature_name, "Coefficient:", coefficient, "\n")
}

