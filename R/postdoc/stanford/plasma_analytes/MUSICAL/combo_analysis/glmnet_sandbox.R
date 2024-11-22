library(glmnet)

targetNames <- unique(big_df2$targetName)

nulisa_for_ridge <- big_df2 %>%
  filter(timepoint=="baseline", class %in% c("A", "S"))%>%
  ungroup()%>%
  pivot_wider(names_from = targetName, values_from = concentration, id_cols = c(id, class))%>%
  mutate(class=if_else(class=="S", 1, 0))

write.csv(nulisa_for_ridge, "~/Downloads/nulisa_for_ridge.csv")
# Getting the independent variable
x_var <- as.matrix(nulisa_for_ridge[,c( targetNames)])
# Getting the dependent variable
y_var <- as.matrix(nulisa_for_ridge[,c("class")])
# Setting the range of lambda values
lambda_seq <- 10^seq(4, -2, by = -.1)
# Using glmnet function to build the ridge regression in r
fit <- glmnet(x_var, y_var, alpha = 0, lambda  = lambda_seq)
# Checking the model
summary(fit)

# Using cross validation glmnet
lambda_list <- list(vector(length=100))
for(i in 1:100){
  ridge_cv <- cv.glmnet(x_var, y_var, alpha = 0, lambda = lambda_seq)
  lambda_list[[i]]=ridge_cv$lambda.min
}
#0.1 most common value
table(unlist(lambda_list))
# Best lambda value


MSEs <- NULL
for (i in 1:100){
  cv <- cv.glmnet(x = x_var, y = y_var, alpha=1,)  
  MSEs <- cbind(MSEs, cv$cvm)
}
rownames(MSEs) <- cv$lambda
lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))


best_ridge <- glmnet(x_var, y_var, alpha = 1, lambda = 0.1061993, family = "binomial")
hist(coef(best_ridge)[1:250])

names(subset(coef(best_ridge)[,1], coef(best_ridge)[,1]!=0))
