# UNIVERSITY CARLOS III, MADRID
# PREDICTIVE MODELIN 18/19 - SI 
# Mar�a del Roc�o Roldán de la Rosa.

## Predictive models about the chemical influence on the liver ##
# Necessary libraries
library(psych)
library(GGally)
library(ggplot2)
library(gridExtra)
library(nortest)
library(lmtest)
library(corrplot)
library(car)
library(factoextra)
library(pls)
library(leaps)
library(glmnet)
library(ISLR)

########################################################
# DATA
dataliver <- read.csv("dataliver.csv", header = TRUE, sep = ";",dec = "," )
View(dataliver)
attach(dataliver)
dataliver
########################################################
# Preliminary analysis
round(cor(x = dataliver, method = "pearson"), 3)
multi.hist(x = dataliver, dcol = c("blue", "red"), dlty = c("dotted", "solid"),
           main = "")

#library(GGally)
ggpairs(dataliver, lower = list(continuous = "smooth"),
        diag = list(continuous = "bar"), axisLabels = "none")
########################################################
# SECTION: 3 STATISTICAL ANALYSIS 
## 3.1 MULTIPLE LINEAR MODEL 

# Response variable: Direct_Bilirubin
# Predictors: Age, Gender, Total_Bilirubin, Alkaline_Phosphotase, 
#Alamine_Aminotransferase, Aspartate_Aminotransferase, 
#Total_Protiens + Albumin, Albumin_and_Globulin_Ratio, Dataset


# 3.1.1. Generate the model (Mixed:All components)
model1 <- lm(Direct_Bilirubin ~ Age+Gender+Total_Bilirubin
             +Alkaline_Phosphotase+Alamine_Aminotransferase
             +Aspartate_Aminotransferase+Total_Protiens+ Albumin
             +Albumin_and_Globulin_Ratio + Dataset, data = dataliver )
summary(model1)
#Residual standard error: 1.333 on 568 degrees of freedom
#(4 observations deleted due to missingness)
#Multiple R-squared:   0.78,	Adjusted R-squared:  0.7762 
#F-statistic: 201.4 on 10 and 568 DF,  p-value: < 2.2e-16

### Selection of the best predictors 
step(object = model1, direction = "both", trace = 1)

########################################################

### 3.1.3 The best model (With the non-significant components)
model2 <- (lm(formula = Direct_Bilirubin ~ Total_Bilirubin + Albumin
              + Total_Protiens + Albumin_and_Globulin_Ratio +
                Alkaline_Phosphotase + Dataset, data = dataliver))
summary(model2)
#Residual standard error: 1336000 on 572 degrees of freedom
#(4 observations deleted due to missingness)
#Multiple R-squared:  0.7773,	Adjusted R-squared: 0.775
#F-statistic: 332.8 on 6 and 572 DF,  p-value: < 2.2e-16

#IC 
# It is important to obtain the CI for each of the predictors for 
#each of the partial regression coefficients. Being, indicative of 
#those who contribute most to the model.
confint(lm(formula = Direct_Bilirubin ~ Total_Bilirubin + Albumin
           + Total_Protiens + Albumin_and_Globulin_Ratio +
             Alkaline_Phosphotase + Dataset, data = dataliver))

########################################################
### 3.1.2.	Model diagnosis
# Validation of conditions for the model
#library(ggplot2)
#library(gridExtra)

### Linearity
par(mfrow=c(1,1))
plot(model2, 3, col = "black", lty = "dotted")

### Normalitity of residuals
# Normality Q-Q 
qqnorm(model2$residuals, col = "black", lty = "dotted")
qqline(model2$residuals, col = "red", lty = "dotted")

shapiro.test(model2$residuals)
#Shapiro-Wilk normality test
#data:  model2$residuals
#W = 0.34063, p-value < 2.2e-16
#It is significant, therefore it complies with the model.

# Lilliefors test - an adaptation of Kolmogorov-Smirnov test
#install.packages("nortest")
#library(nortest)
lillie.test(model2$residuals)
#Lilliefors (Kolmogorov-Smirnov) normality test
#data:  model2$residuals
#D = 0.28592, p-value < 2.2e-16

### Homoscedasticity (Constant variability of waste)
#library(lmtest)
bptest(model2)
#studentized Breusch-Pagan test
#data:  model2
#BP = 180.35, df = 6, p-value < 2.2e-16


### Independence of the errors
plot(model2$residuals, type = "l", col="darkslategray4")
abline(h=0, col= "Red")

### Multicollinearity 
#library(corrplot)
corrplot(cor(dplyr::select(dataliver, Total_Bilirubin, Albumin, 
                           Total_Protiens, Albumin_and_Globulin_Ratio,
                          Alkaline_Phosphotase, Dataset)),
         method = "number", tl.col = "black")

modMultiCo <- lm(Direct_Bilirubin ~ Total_Bilirubin + Albumin
                 + Total_Protiens + Albumin_and_Globulin_Ratio +
                   Alkaline_Phosphotase + Dataset, data = dataliver)
summary(modMultiCo)
#Residual standard error: 1336000 on 572 degrees of freedom
#(4 observations deleted due to missingness)
#Multiple R-squared:  0.7773,	Adjusted R-squared:  0.775 
#F-statistic: 332.8 on 6 and 572 DF,  p-value: < 2.2e-16
# It is significant, therefore it complies with the model.


# Autocorrelación 
#library(car)
dwt(model2, alternative = "two.sided")
#lag Autocorrelation D-W Statistic p-value
#1       0.06717512      1.865502   0.054
#Alternative hypothesis: rho != 0
# The p-value is not significant in this case, but not by a high
#value so it can be said that it can be said that there are indications 
#of almost colinearity.

#Extra

# VIF____ quantifies the intensity of multicollinearity in a normal
#least squares regression analysis. In this case, you get:
#library(car)
vif(model2)
#Total_Bilirubin  1.215992
#Albumin 9.682426 
#Total_Protiens  5.304848
#Albumin_and_Globulin_Ratio 3.656397
#Alkaline_Phosphotase 1.105232 
#Dataset 1.093449

#### OUTLIERS____
plot(model2, 5, col = "black")
########################################################

## **Other techniques** ##
# Influence of the predictors of the model in id.
summary(influence.measures(model2))

# Visualicemos la influencia
influencePlot(model2)
# Hay exitencia de valores que son muy influyentes
#StudRes        Hat       CookD
#167 -51.4580352 0.267900949 24.59786115
#247   6.2405582 0.003932087  0.02059627
#532   4.5847472 0.042871600  0.12995469
#549  -0.7844221 0.277444997  0.03377534
#576   0.8180896 0.318860656  0.044783698

#It shows that there are several influential observations
# (position 167 and 247) that exceed the limits of concern
# for the values of Leverages or Cook Distance.


########################################################
#3.3.	 DIMENSION REDUCTION TECHNIQUES 
## 3.3.1.  PCA & PCR
dataRo <- subset(dataliver, select = -c(Dataset))

dataliver <- read.csv("dataliver.csv", header = TRUE, sep = ";",dec = "," )
Sample.scaled <- data.frame(apply(dataRo,2,scale))
Sample.scaled.2 <- data.frame(t(na.omit(t(Sample.scaled))))
pca.Sample.2 <- prcomp(Sample.scaled.2, retx=TRUE)

apply(X = dataliver, MARGIN = 2, FUN = mean)
apply(X = dataliver, MARGIN = 2, FUN = var)

## use cor = TRUE to standarize variables
pca <- prcomp(~ ., data=dataRo, na.action=na.omit, scale=TRUE)
summary(pca)

pca$center
pca$scale
pca$rotation
head(pca$x)

### Dimension of the main components
dim(pca$x)

# Composition of each component
#library(corrplot)
corrplot::corrplot(cor(dataRo,pca$scores), addCoef.col = "darkslategray3")
plot(pca, type = "l")

#install.packages("factoextra")
#library("factoextra")
# Plot of the percentage explained variance by each PC  
biplot(x = pca, scale = 0, cex = 0.6, col = c("blue4", "brown3"))

pca$rotation <- -pca$rotation
pca$x        <- -pca$x
biplot(x = pca, scale = 0, cex = 0.6, col = c("blue4", "brown3"))


########################################################
#3.4.	BEST MODEL SELECTION
#3.4.1.	ANOVA
# An ANOVA can be used to perform a hypothesis test and
#get a p-value that evaluates the null hypothesis that both models
# they fit the data just as well.
model3 <- update(model2, .~. -Dataset -Alkaline_Phosphotase)
summary(model3)

#Call:
#  lm(formula = Direct_Bilirubin ~ Total_Bilirubin + Albumin + Total_Protiens + 
#       Albumin_and_Globulin_Ratio, data = dataliver)
#Residual standard error: 1346000 on 574 degrees of freedom
#(4 observations deleted due to missingness)
#Multiple R-squared:  0.7732,	Adjusted R-squared:  0.7716 
#F-statistic: 489.2 on 4 and 574 DF,  p-value: < 2.2e-16

anova(model2, model3)
#Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
#1    572 1.0211e+15                                  
#2    574 1.0400e+15 -2 -1.8914e+13 5.2978 0.005251 **

# How you can influence between Direct_Bilirubin and Albumin_and_Globulin_Ratio
modelo_lineal <- lm(formula = Direct_Bilirubin ~ Albumin_and_Globulin_Ratio, data = dataliver)
summary(modelo_lineal)

#Residual standard error: 2762000 on 577 degrees of freedom
#Multiple R-squared:  0.04005,	Adjusted R-squared:  0.03839 
#F-statistic: 24.07 on 1 and 577 DF,  p-value: 1.208e-06

plot(x = Albumin_and_Globulin_Ratio, y = Direct_Bilirubin, main = "Influencia en la enfermedad", pch = 20,
     col = "grey")
abline(modelo_lineal, lwd = 3, col = "red")

########################################################
## 3.2. SHRINKAGE METHODS 
coef(lm(Direct_Bilirubin ~ Total_Bilirubin + Albumin + Total_Protiens + 
          Albumin_and_Globulin_Ratio, data = dataliver))

#According to the model, for each unit that increases Direct_Bilirubin of suffering
#the disease, varying the other predictors, Albumin_and_Globulin_Ratio, Total_Bilirubin
#and Total_Protiens increases on average 1032132.3, 379443.1 and 505984.8

## 3.2.1. RIDGE REGRESSION 
x <- model.matrix(Direct_Bilirubin~., data = dataliver)[, -1]
head(x)
y <- dataliver$Direct_Bilirubin

#library(glmnet)
# Para obtener un ajuste mediante ridge regression se indica argumento alpha=0.
modelos_ridge <- glmnet(x = x, y = y, alpha = 0)

set.seed(1)
cv_error_ridge <- cv.glmnet(x = x, y = y, alpha = 0, nfolds = 10,
                            type.measure = "mse")
plot(cv_error_ridge)

# Lambda value with which the minimum test-error is achieved
cv_error_ridge$lambda.min
#[1] 246084.7

# Optimal lambda value: higher lambda value with which the test-error is not
# moves more than 1 sd away from the minimum possible error-test.
cv_error_ridge$lambda.1se
#[1] 4401592

# The value of the coefficients for the optimum lambda value is displayed
modelo_final_ridge <- glmnet(x = x, y = y, alpha = 0, lambda = 4401592)
coef(modelo_final_ridge)

########################################################
# 3.2.2.LASSO
#library(glmnet)
modelos_lasso <- glmnet(x = x, y = y, alpha = 1)
plot(modelos_lasso, xvar = "lambda", label = TRUE)

set.seed(1)
cv_error_lasso <- cv.glmnet(x = x, y = y, alpha = 1, nfolds = 10)
plot(cv_error_lasso)

cv_error_lasso$lambda.min
#[1] 219068.3
cv_error_lasso$lambda.1se
#[1] 1408187

# The model with all observations is readjusted using the value of
# optimum lambda
modelo_final_lasso <- glmnet(x = x, y = y, alpha = 1, lambda = cv_error_lasso$lambda.1se)
coef(modelo_final_lasso)

par(mfrow = c(1,2))
plot(cv_error_ridge,ylab = "Mean Square Error ridge regression" )
abline(h = 120000)
plot(cv_error_lasso,ylab = "Mean Square Error lasso")
abline(h = 120000)
#par(mfrow = c(1,1))

####################################################
#3.3.	 Dimension reduction techniques  +++
## 3.3.2. Principal Components Regression PCR
library(pls)
set.seed(2)
# Important to standardize the variables indicating it with the argument scale
# Indicating validation = CV, 10-fold-cross-validation is used to
# Identify the optimal number of components.
modelo_pcr <- pcr(Direct_Bilirubin~., data = dataliver, scale = TRUE, validation = "CV")
summary(modelo_pcr)

#Express the cumulative explained variance and with 7 you get
#more than 90%.
validationplot(modelo_pcr, val.type = "RMSEP")

## PRESS (Other)
# PRESS es el Predicted Sum of Squares
plot(as.numeric(sqrt(modelo_pcr$validation$PRESS)), type = "b", pch = 19,
     ylab = expression(sqrt("PRESS")))
axis(side = 1, at = 1:19)

#To know the number of components with which the error is minimized
which.min(x = modelo_pcr$validation$PRESS)
#[1] 8
# It is observed that the smallest error occurs with 6 components.

###################################################################3
# Best models with other techniques
#library(leaps)
mejores_modelos <- regsubsets(Direct_Bilirubin~., data = dataliver, nvmax = 19)
# The nvmax argument determines the maximum size of the models to inspect.
# If you want to perform best subset selection evaluating all possible
# models, nvmax has to be equal to the number of available variables
summary(mejores_modelos)
names(summary(mejores_modelos))
summary(mejores_modelos)$adjr2
summary(mejores_modelos)$bic
summary(mejores_modelos)$rss
summary(mejores_modelos)$rsq
summary(mejores_modelos)$which
summary(mejores_modelos)$cp
#it is identified which model has the maximum value of R adjusted
which.max(summary(mejores_modelos)$adjr2)#[1] 7
which.max(summary(mejores_modelos)$bic)#[1] 10
which.max(summary(mejores_modelos)$rss)#[1] 1
which.max(summary(mejores_modelos)$rsq)#[1] 10
which.max(summary(mejores_modelos)$which)#[1] 1
which.max(summary(mejores_modelos)$cp)#[1] 1

#The model that occupies position 7 is the one that has the highest R2adjusted
#achieves. Therefore, the best model is the one that contains 7 predictors.

coef(object = mejores_modelos, id = 7)

#(Intercept)                  Total_Bilirubin 
#-974615.2165                 368125.4581 
#Alkaline_Phosphotase         Alamine_Aminotransferase 
#512.3619                     755.3909 
#Total_Protiens               Albumin 
#511529.6084                  -978478.4404 
#Albumin_and_Globulin_Ratio   Dataset 
#1127335.0419                 -213672.3463 

summary(mejores_modelos)$adjr2[7]
#[1] 0.7768209 #better model
summary(mejores_modelos)$adjr2[6]
#[1] 0.7761427
summary(mejores_modelos)$adjr2[5]
#[1] 0.7744075

######################################################
# 3.4.2 AIC & BIC
# BICs
BIC(model1) #18039.16
BIC(model2) #18020.86 #smaller = better 

#AICs
AIC(model1) #17986.83 
AIC(model2) #17985.97 #smaller = better 

#Check the summaries
summary(model1) #R=0.78   R2= 0.7762 #better
summary(model2) #R=0.7773 R2= 0.775 

# Full model = model2 = called DATA
Data <- subset(dataliver, select = -c(Gender, Age ))
mod <- lm(Direct_Bilirubin ~ ., data = Data)
mod

######################################################
# Comparison between methods
# 4. CONCLUSIONS
library(ISLR)
dataliver <- na.omit(dataliver)
dim(dataliver)
#[1] 579  11

set.seed(1)
indices_entrenamiento <- sample(x = 1:nrow(dataliver),
                                size = round(nrow(dataliver) * (2/3)))
# 2/3 of the observations
indices_test <- (1:nrow(dataliver))[-indices_entrenamiento]
dataliver_1 <- dataliver[indices_entrenamiento,]
dataliver_2 <- dataliver[indices_test,]

### Regression by least squares
modelo_OLS   <- lm(formula = Direct_Bilirubin ~ ., data = dataliver_1)
test_MSE_OLS <- mean((predict(modelo_OLS, dataliver_2) - dataliver_2$Direct_Bilirubin)^2)
test_MSE_OLS
#[1] 1.715697e+12

### Ridge regression
# The glmnet () function requires passing the predictors as the matrix and the variable
# dependent as a vector.
x_dataliver_1 <- model.matrix(Direct_Bilirubin~., data = dataliver_1)[, -1]
y_dataliver_1 <- dataliver_1$Direct_Bilirubin

x_dataliver_2 <- model.matrix(Direct_Bilirubin~., data = dataliver_2)[, -1]
y_dataliver_2 <- dataliver_2$Direct_Bilirubin


set.seed(1)
# The best lambda value is identified by k-cross-validation
# ridge regression
cv_error_ridge <- cv.glmnet(x = x_dataliver_1, y = y_dataliver_1, alpha = 0,
                            nfolds = 10,
                            type.measure = "mse")
# To obtain a setting using * ridge regression * argument is indicated
# alpha = 0
modelo_ridge <- glmnet(x = x_dataliver_1, y = y_dataliver_1, alpha = 0,
                       lambda = cv_error_ridge$lambda.1se)
# Predictions are stored in a separate variable so as not to concatenate
# so much code
predicciones <- predict(object = modelo_ridge, newx = x_dataliver_2,
                        s = cv_error_ridge$lambda.1se, exact = TRUE)
test_MSE_ridge <- mean((predicciones - dataliver_2$Direct_Bilirubin)^2)
test_MSE_ridge
#[1] 6.62969e+12

###Lasso
set.seed(1)
# The best lambda value is identified by k-cross-validation
# lasso regresion
cv_error_lasso <- cv.glmnet(x = x_dataliver_1, y = y_dataliver_1, alpha = 1,
                            nfolds = 10,
                            type.measure = "mse")
# To obtain an adjustment using * lasso regression * argument is indicated
# alpha = 1
modelo_lasso <- glmnet(x = x_dataliver_1, y = y_dataliver_1, alpha = 1,
                       lambda = cv_error_ridge$lambda.1se)
# Predictions are stored in a separate variable so as not to concatenate
# so much code
predicciones <- predict(object = modelo_ridge, newx = x_dataliver_2,
                        s = cv_error_ridge$lambda.1se, exact = TRUE)

test_MSE_lasso <- mean((predicciones - dataliver_2$Direct_Bilirubin)^2)
test_MSE_lasso
#[1] 6.62969e+12

### PCR
#library(pls)
set.seed(233)
# Important to standardize the variables indicating it with the argument
# scale = TRUE Indicating validation = CV, 10-fold-cross-validation is used
# to identify the optimal number of components.
modelo_pcr <- pcr(Direct_Bilirubin~., data = dataliver_1, scale = TRUE, validation = "CV")
validationplot(modelo_pcr, val.type = "RMSEP")

# See in more detail from component 1
# PRESS is the Predicted Sum of Squares
plot(as.numeric(sqrt(modelo_pcr$validation$PRESS)), type = "b", pch = 19,
     ylab = expression(sqrt("PRESS")))
axis(side = 1, at = 1:19)

#To know the number of components with which the error is minimized
which.min(x = modelo_pcr$validation$PRESS)
#[1] 7

predicciones <- predict(object = modelo_pcr, newdata = dataliver_2, ncomp = 6)
test_MSE_PCR <- mean((predicciones - dataliver_2$Direct_Bilirubin)^2)
test_MSE_PCR
#[1] 3.143411e+12

# Even though the minimum PRESS is reached
#with 7 components, using only 7 you get a very valuable
#considerable.

# Conclusion
metodo <- c("OLS", "ridge regresion", "lasso", "PCR")
test_MSE <- c(test_MSE_OLS, test_MSE_ridge, test_MSE_lasso,
              test_MSE_PCR)
resultados <- data.frame(metodo, test_MSE)
resultados
#             metodo            test_MSE
#1           OLS                1.715697e+12
#2           ridge regresion    6.629690e+12
#3           lasso              6.629690e+12
#4           PCR                3.143411e+12

ggplot(data = resultados, aes(x = reorder(metodo, test_MSE),
                              y = sqrt(test_MSE))) +
  geom_bar(stat = "identity") +
  labs(x = "método de regresión", y = expression(sqrt("tes-MSR"))) +
  theme_bw()
