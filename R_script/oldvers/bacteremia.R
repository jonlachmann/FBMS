#devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)

library(Metrics)# For C-index
library(FBMS)
set.seed(123)

# Function to calculate C-index manually for binary classification
cindex_manual <- function(predictions, labels) {
  # Ensure predictions and labels are numeric
  n <- length(labels)
  
  # Initialize counters for concordant, discordant, and ties
  concordant <- 0
  discordant <- 0
  ties <- 0
  
  # Loop through all possible pairs
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (labels[i] != labels[j]) {  # Only consider pairs with different labels
        if (predictions[i] == predictions[j]) {
          ties <- ties + 1
        } else if ((predictions[i] > predictions[j] && labels[i] > labels[j]) ||
                   (predictions[i] < predictions[j] && labels[i] < labels[j])) {
          concordant <- concordant + 1
        } else {
          discordant <- discordant + 1
        }
      }
    }
  }
  
  # Calculate the C-index
  total_pairs <- concordant + discordant + ties
  c_index <- (concordant + 0.5 * ties) / total_pairs
  return(c_index)
}


df = read.csv("https://zenodo.org/records/7554815/files/Bacteremia_public_S2.csv", 
              header = TRUE, sep = ",", dec = ".")

df = df[,!(names(df) %in% c("MONOR", "LYMR", "NEUR", "EOSR", "BASOR", "WBC", "MCV", "HCT"))]
df$BloodCulture = ifelse(df$BloodCulture == "yes", 1, 0)

#df = na.omit(df)
trid = sample.int(dim(df)[1],round(dim(df)[1]*2/3))

df.train = (df[trid,]) 
df.test = (df[-trid,]) 


# Number of bootstrap iterations
n_bootstrap <- 100

# Store results for each bootstrap iteration
accuracy_oob <- numeric(n_bootstrap)
accuracy_boot <- numeric(n_bootstrap)
cindex_oob <- numeric(n_bootstrap)
cindex_boot <- numeric(n_bootstrap)

# Full model performance
result.nonlinear = fbms(formula = BloodCulture ~ 1 + ., family = "binomial", data = df.train,beta_prior = list(type = "Jeffreys-BIC"), impute = T,  method = "gmjmcmc.parallel",P = 10,cores = 6, runs = 6, transforms = c("sigmoid","sin","cos","exp_dbl"))
summary(result.nonlinear)

# Full model performance
result = fbms(formula = BloodCulture ~ 1 + ., family = "binomial", data = df.train,beta_prior = list(type = "Jeffreys-BIC"), impute = T, method = "mjmcmc")
summary(result)
preds = predict(result.nonlinear, df.test[,-45])
prob_full <- sigmoid(preds$aggr$mean)

# Accuracy on full model
accuracy_full <- mean((prob_full > 0.5) == df.test$BloodCulture)

# AUC on full model
auc_full <- auc(df.test$BloodCulture, prob_full)

# C-index on full model
cindex_full <- cindex_manual(prob_full, df.test$BloodCulture)

# Theoretical performance (no effect) for AUC and C-index
p_M0_auc <- 0.5  # For AUC
p_M0_cindex <- 0.5  # For C-index
p_M0_accuracy <- 0.5  # For accuracy (random guessing)

# Bootstrap procedure
for (i in 1:n_bootstrap) {
  
  # Bootstrap resample from training data
  boot_indices <- sample(1:nrow(df.train), replace = TRUE)
  df_boot <- df.train[boot_indices, ]
  
  # Fit model on bootstrap sample
  result_boot <- fbms(formula = BloodCulture ~ 1 + ., family = "binomial", data = df_boot, impute = TRUE, method = "mjmcmc")
  
  # Predictions on bootstrap sample (in-sample performance)
  preds_boot <- predict(result_boot, df_boot[,-45])
  prob_boot <- sigmoid(preds_boot$mean)
  
  # Predictions on the original data (out-of-bag performance)
  preds_oob <- predict(result_boot, df.test[,-45])
  prob_oob <- sigmoid(preds_oob$mean)
  
  # Accuracy
  accuracy_boot[i] <- mean((prob_boot > 0.5) == df_boot$BloodCulture)
  accuracy_oob[i] <- mean((prob_oob > 0.5) == df.test$BloodCulture)
  
  # AUC
  auc_boot[i] <- auc(df_boot$BloodCulture, prob_boot)
  auc_oob[i] <- auc(df.test$BloodCulture, prob_oob)
  
  # C-index
  cindex_boot[i] <- cindex_manual(prob_boot, df_boot$BloodCulture)
  cindex_oob[i] <- cindex_manual(prob_oob, df.test$BloodCulture)
}

# Calculate overfitting rate R for accuracy
R_accuracy <- mean(accuracy_boot - accuracy_oob) / (p_M0_accuracy - accuracy_full)

# Calculate overfitting rate R for AUC
R_auc <- mean(auc_boot - auc_oob) / (p_M0_auc - auc_full)

# Calculate overfitting rate R for C-index
R_cindex <- mean(cindex_boot - cindex_oob) / (p_M0_cindex - cindex_full)

# .632+ weight for accuracy
w_accuracy <- 0.632 / (1 - 0.368 * R_accuracy)

# .632+ weight for AUC
w_auc <- 0.632 / (1 - 0.368 * R_auc)

# .632+ weight for C-index
w_cindex <- 0.632 / (1 - 0.368 * R_cindex)
# .632+ estimate
cindex_632plus <- (1 - w_cindex) * cindex_full + w_cindex * mean(cindex_boot)
accuracy_632plus <- (1 - w_accuracy) * accuracy_full + w_accuracy * mean(accuracy_boot)
auc_632plus <- (1 - w_auc) * auc_full + w_auc * mean(auc_boot)

# Confidence intervals (95% CI using bootstrap percentiles)
cindex_ci <- quantile(cindex_boot, probs = c(0.025, 0.975))
accuracy_ci <- quantile(accuracy_boot, probs = c(0.025, 0.975))
auc_ci <- quantile(auc_boot, probs = c(0.025, 0.975))

# Print results
cat("Full model C-index:", cindex_full, "\n")
cat(".632+ C-index estimate:", cindex_632plus, "\n")
cat("C-index 95% CI from bootstrap:", cindex_ci, "\n\n")

cat("Full model Accuracy:", accuracy_full, "\n")
cat(".632+ Accuracy estimate:", accuracy_632plus, "\n")
cat("Accuracy 95% CI from bootstrap:", accuracy_ci, "\n\n")

cat("Full model AUC:", auc_full, "\n")
cat(".632+ AUC estimate:", auc_632plus, "\n")
cat("AUC 95% CI from bootstrap:", auc_ci, "\n")