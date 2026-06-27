library(dplyr)
library(tidyr)
library(glmnet)



set.seed(1)

seroprev_est_df <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/mixture_model_df.csv")

long_luminex <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/lavstsen/long_luminex.csv")

# LASSO: use antigen-specific log_mfi to predict symp_prop_12_24
# Assumes df has at least: id (or id_cat), antigen, log_mfi, symp_prop_12_24
set.seed(1)

# 1) Keep needed cols and make one row per person (wide: one column per antigen)
dat_wide <- df %>%
  filter( timepoint=="52 weeks")%>%
  transmute(
    id = if ("id" %in% names(df)) id else id_cat,
    antigen,
    log_mfi,
    y = symp_prop_12_24
  ) %>%
  filter(!is.na(id), !is.na(antigen),) %>%
  pivot_wider(names_from = antigen, values_from = log_mfi)

# 2) Build model matrix X and outcome y
y <- dat_wide$y
X_df <- dat_wide %>% select(-id, -y)

# Drop rows with missing outcome
keep <- !is.na(y)
y <- y[keep]
X_df <- X_df[keep, , drop = FALSE]

# glmnet can't take NA in X; impute missing log_mfi as 0 (common with sparse panels)
# Alternatives: column means, KNN, etc.
X_df <- X_df %>% mutate(across(everything(), ~ replace_na(.x, 0)))

# Convert to numeric matrix for glmnet
X <- as.matrix(X_df)

# 3) Fit LASSO with cross-validation (Gaussian since y is continuous proportion)
cvfit <- cv.glmnet(
  x = X, y = y,
  alpha = 0.5,              # lasso
  family = "gaussian",
  standardize = TRUE,
  nfolds = 10
)

# 4) Extract selected antigens (non-zero coefficients)
coef_min <- coef(cvfit, s = "lambda.min")
coef_1se <- coef(cvfit, s = "lambda.1se")

selected_min <- rownames(coef_min)[as.vector(coef_min != 0)]
selected_1se <- rownames(coef_1se)[as.vector(coef_1se != 0)]

selected_min <- setdiff(selected_min, "(Intercept)")
selected_1se <- setdiff(selected_1se, "(Intercept)")

cat("Selected antigens @ lambda.min:\n")
print(selected_min)


# 5) Save a tidy coefficient table (lambda.min)
coef_tbl <- data.frame(
  antigen = rownames(coef_min),
  beta = as.numeric(coef_min)
) %>%
  filter(antigen != "(Intercept)", beta != 0) %>%
  arrange(desc(abs(beta)))

print(coef_tbl)

# Optional: quick plot of CV curve
plot(cvfit)

# Optional: simple in-sample predictions at lambda.1se
pred <- as.numeric(predict(cvfit, newx = X, s = "lambda.1se"))
# e.g., correlation between predicted and observed
cat("\nCorrelation(y, pred) @ lambda.1se:", cor(y, pred), "\n")