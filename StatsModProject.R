#Loading the required libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(survival)
library(survminer)
library(mice)
library(gridExtra)
#Loading the Dataset
data <- read.csv("/Users/manasarthak/Downloads/cirrhosis.csv")
head(data)
# Convert days to years for age
data$Age <- data$Age / 365.25

# Calculate the percentage of missing data for each column
missing_data_summary<-data %>%
  summarise(across(everything(), ~sum(is.na(.))/n()*100)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "MissingPercentage")%>%
  mutate(Variable = factor(Variable, levels = Variable[order(MissingPercentage)]))

ggplot(missing_data_summary, aes(x = Variable, y = MissingPercentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Percentage of Missing Data by Variable", y = "Percentage Missing", x = "Variables")


md.pattern(data)
#beautification
install.packages("VIM")
library(VIM)
# Aggregated missing data plot
aggr(data, col=c('navyblue', 'yellow'), numbers=TRUE, sortVars=TRUE,
     labels=names(data), cex.axis=.7,
     gap=3, ylab=c("Missing data","Pattern"))


# Convert 'Drug' to factor with appropriate levels
data$Drug <- factor(data$Drug, levels = c("D-penicillamine", "Placebo"), labels = c("D-penicillamine", "Placebo"))
data$Status <- factor(data$Status, levels = c("C", "D","CL"), labels = c("C", "D","CL"))
data$Sex <- factor(data$Sex, levels = c("F", "M"), labels = c("F", "M"))

# List of true binary variables
binary_vars <- c("Edema","Ascites", "Hepatomegaly", "Spiders")

# Convert binary variables to factors with levels "N" and "Y"
data[binary_vars] <- lapply(data[binary_vars], function(x) {
  factor(x, levels = c("N", "Y"), labels = c("No", "Yes"))
})

# Check the structure to confirm changes
str(data[c("Drug", binary_vars)])

# Continuous Variables
continuous_vars <- c("Bilirubin", "Cholesterol", "Albumin", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin")
df_select <- data[,continuous_vars]
library(naniar)
mcar_test(data=df_select)

continuous_plots <- lapply(continuous_vars, function(var) {
  ggplot(data, aes_string(x = var)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", alpha = 0.7) +
    geom_density(color = "red") +
    labs(title = paste("Distribution of", var), x = var, y = "Density") +
    theme_minimal()
})

# Categorical Variables
categorical_vars <- c("Status", "Drug", "Sex", "Ascites", "Hepatomegaly", "Spiders", "Edema", "Stage")

# Create plots for each categorical variable
categorical_plots <- lapply(categorical_vars, function(var) {
  ggplot(data, aes_string(x = var)) +
    geom_bar(fill = "steelblue") +
    labs(title = paste("Distribution of", var), x = var, y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
})
# Combine all plots
all_plots <- c(continuous_plots, categorical_plots)

# Use grid.arrange to display all plots together
do.call(grid.arrange, c(all_plots, ncol = 3))

# Initialize the MICE model for PMM
init_data <- mice(data, method = 'pmm', m = 5, seed = 123, printFlag = FALSE)
summary(init_data)

#visualize the "similarity" between 
densityplot(init_data)


#For surival curves we will use one of the imputed datasets
#One of the 5 imputed datasets
ex_data<-complete(init_data,4)

# Convert Status to a survival object, assuming 'D' stands for the event (death)
surv_obj <- Surv(ex_data$N_Days, ex_data$Status == "D")

# Fit Kaplan-Meier survival curves
km_fit <- survfit(surv_obj ~ Sex, data = ex_data)

# Plot the Kaplan-Meier survival curves
g <- ggsurvplot(
  km_fit, 
  data = data,
  pval = TRUE,                 # Show p-value of the log-rank test
  conf.int = TRUE,             # Show confidence intervals
  palette = "Dark2",           # Color palette
  xlab = "Days",               # Label for the x-axis
  ylab = "Survival probability", # Label for the y-axis
  title = "Kaplan-Meier Survival Curve by Sex" # Title of the plot
)

print(g)

# Fit Kaplan-Meier survival curve startified by Drug(primary variable for analysis)
km_fit <- survfit(surv_obj ~ Drug, data = ex_data)

# Plot the Kaplan-Meier survival curves
g <- ggsurvplot(
  km_fit, 
  data = ex_data,
  pval = TRUE,                 # Show p-value of the log-rank test
  conf.int = TRUE,             # Show confidence intervals
  palette = "Dark2",           # Color palette
  xlab = "Days",               # Label for the x-axis
  ylab = "Survival probability", # Label for the y-axis
  title = "Kaplan-Meier Survival Curve by Drug" # Title of the plot
)

print(g)

# Fit Kaplan-Meier survival curve startified by Ascites
km_fit <- survfit(surv_obj ~ Ascites, data = ex_data)

# Plot the Kaplan-Meier survival curves
g <- ggsurvplot(
  km_fit, 
  data = ex_data,
  pval = TRUE,                 # Show p-value of the log-rank test
  conf.int = TRUE,             # Show confidence intervals
  palette = "Dark2",           # Color palette
  xlab = "Days",               # Label for the x-axis
  ylab = "Survival probability", # Label for the y-axis
  title = "Kaplan-Meier Survival Curve by Ascites" # Title of the plot
)

print(g)

# Fit Kaplan-Meier survival curve startified by Hepatomegaly
km_fit <- survfit(surv_obj ~ Hepatomegaly, data = ex_data)

# Plot the Kaplan-Meier survival curves
g <- ggsurvplot(
  km_fit, 
  data = ex_data,
  pval = TRUE,                 # Show p-value of the log-rank test
  conf.int = TRUE,             # Show confidence intervals
  palette = "Dark2",           # Color palette
  xlab = "Days",               # Label for the x-axis
  ylab = "Survival probability", # Label for the y-axis
  title = "Kaplan-Meier Survival Curve by Hepatomegaly" # Title of the plot
)

print(g)

# Fit Kaplan-Meier survival curve startified by Spider
km_fit <- survfit(surv_obj ~ Spiders, data = ex_data)

# Plot the Kaplan-Meier survival curves
g <- ggsurvplot(
  km_fit, 
  data = ex_data,
  pval = TRUE,                 # Show p-value of the log-rank test
  conf.int = TRUE,             # Show confidence intervals
  palette = "Dark2",           # Color palette
  xlab = "Days",               # Label for the x-axis
  ylab = "Survival probability", # Label for the y-axis
  title = "Kaplan-Meier Survival Curve by Spiders" # Title of the plot
)

print(g)

# Fit Kaplan-Meier survival curve startified by Spider
km_fit <- survfit(surv_obj ~ Edema, data = ex_data)

# Plot the Kaplan-Meier survival curves
g <- ggsurvplot(
  km_fit, 
  data = ex_data,
  pval = TRUE,                 # Show p-value of the log-rank test
  conf.int = TRUE,             # Show confidence intervals
  palette = "Dark2",           # Color palette
  xlab = "Days",               # Label for the x-axis
  ylab = "Survival probability", # Label for the y-axis
  title = "Kaplan-Meier Survival Curve by Edema" # Title of the plot
)

print(g)

#Cox Regression

# Fit a Cox proportional hazards model to each imputed dataset
cox_models <- with(data = init_data, exp = {
  # Ensure the event is correctly coded within the function scope
  event <- ifelse(Status == "D", 1, 0)
  surv_obj <- Surv(time = N_Days, event = event)
  coxph(surv_obj ~ Drug+Ascites+Hepatomegaly+Spiders+Age+Sex+Bilirubin+Albumin +Cholesterol+ Copper+ Alk_Phos+ SGOT+ Tryglicerides+Platelets+Prothrombin)
})
# Pool the results of fitting the model to each imputed dataset
pooled_results <- pool(cox_models)

# Print the summary of the pooled results
df<-as_tibble(summary(pooled_results))
df$significance <- ifelse(df$p.value < 0.001, "***",
                                  ifelse(df$p.value < 0.01, "**",
                                         ifelse(df$p.value < 0.05, "*", "")))
df


# Assuming 'ex_data' is your dataset with imputed values
# Ensure the Status column is correctly encoded for the event (1 for event occurred, 0 for censored)
ex_data$event <- ifelse(ex_data$Status == "D", 1, 0)

# Create the survival object
surv_obj <- Surv(time = ex_data$N_Days, event = ex_data$event)

# Fit the Cox proportional hazards model
cox_model <- coxph(surv_obj ~ Drug+Ascites+Hepatomegaly+Spiders+Age+Sex+Bilirubin+Albumin +Cholesterol+ Copper+ Alk_Phos+ SGOT+ Tryglicerides+Platelets+Prothrombin, data = ex_data)

# Summary of the Cox model
summary(cox_model)

library(corrplot)
continuous_vars <- ex_data[, c("Age", "Bilirubin", "Albumin", "Prothrombin", "Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets")]

# Compute correlation matrix
cor_matrix <- cor(continuous_vars, use = "complete.obs")  # Handling missing values

# Plot the correlation matrix
corrplot(cor_matrix, method = "circle", type = "upper", 
         order = "hclust", 
         tl.col = "black", tl.srt = 30,  # Text label color and rotation
         title = "Correlation Matrix of Continuous Predictors")