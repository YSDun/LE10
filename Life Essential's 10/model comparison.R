#install.packages("BiocManager")
#BiocManager::install("survcomp")
#install.packages("compareC")
library(survcomp)
library(compareC)
library(survival)
library(rms) 
library(pec)

#LE8
model_LE8 <- coxph(Surv(time,status)~age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score, data=data_cox, x=TRUE)
Cindex_LE8 <- summary(model_LE8)$concordance
Cindex_LE8
#LE8+alcohol consumption score
model_LE9 <- coxph(Surv(time,status)~age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+alcohol_score, data=data_cox, x=TRUE)
Cindex_LE9 <- summary(model_LE9)$concordance
Cindex_LE9
#LE8+mental health score
model_LC9 <- coxph(Surv(time,status)~age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+mental_health_score, data=data_cox, x=TRUE)
Cindex_LC9 <- summary(model_LC9)$concordance
Cindex_LC9
#LE10
model_LE10 <- coxph(Surv(time,status)~age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+alcohol_score+mental_health_score, data=data_cox, x=TRUE)
Cindex_LE10 <- summary(model_LE10)$concordance
Cindex_LE10

#####################################################################
# Calculate C-indexes with survcomp package
#####################################################################
Cindex_LE8 <- concordance.index(x=predict(model_LE8),
                                surv.time=data_cox$time,
                                surv.event=data_cox$status,
                                method="noether")
Cindex_LE8

#LE8+alcohol consumption score
Cindex_LE9 <- concordance.index(x=predict(model_LE9),
                                surv.time=data_cox$time,
                                surv.event=data_cox$status,
                                method="noether")
Cindex_LE9

#LE8+alcohol consumption score
Cindex_LC9 <- concordance.index(x=predict(model_LC9),
                                surv.time=data_cox$time,
                                surv.event=data_cox$status,
                                method="noether")
Cindex_LC9

#LE8+mental health score
Cindex_LE10 <- concordance.index(x=predict(model_LE10),
                                 surv.time=data_cox$time,
                                 surv.event=data_cox$status,
                                 method="noether")
Cindex_LE10

cindex.comp1 <- cindex.comp(Cindex_LC9, Cindex_LE8)
cindex.comp2 <- cindex.comp(Cindex_LE9, Cindex_LE8)
cindex.comp3 <- cindex.comp(Cindex_LE10, Cindex_LE8)

p.value_decimal1 <- format(cindex.comp1$p.value, scientific = FALSE)
p.value_decimal2 <- format(cindex.comp2$p.value, scientific = FALSE)
p.value_decimal3 <- format(cindex.comp3$p.value, scientific = FALSE)

#####################################################################
# Calculate time-dependent AUC
#####################################################################
pred_LE8 <- predict(model_LE8, type = "lp")
pred_LE9 <- predict(model_LE9, type = "lp")
pred_LC9 <- predict(model_LC9, type = "lp")
pred_LE10 <- predict(model_LE10, type = "lp")

time_points <- seq(0, 15.6, by = 0.1)
time_roc_LE8 <- timeROC(T = data_cox$time, delta = data_cox$status, marker = pred_LE8, cause = 1, times = time_points)
time_roc_LC9 <- timeROC(T = data_cox$time, delta = data_cox$status, marker = pred_LC9, cause = 1, times = time_points)
time_roc_LE9 <- timeROC(T = data_cox$time, delta = data_cox$status, marker = pred_LE9, cause = 1, times = time_points)
time_roc_LE10 <- timeROC(T = data_cox$time, delta = data_cox$status, marker = pred_LE10, cause = 1, times = time_points)


# Fit the logistic model for KDM-BA acceleration
model1 <- glm(kdm_AA_binary ~ age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score, data = kdm_UK, family = binomial, x = TRUE,y = TRUE)
model2 <- glm(kdm_AA_binary ~ age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+mental_health_score, data = kdm_UK, family = binomial, x = TRUE,y = TRUE)
model3 <- glm(kdm_AA_binary ~ age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+alcohol_score, data = kdm_UK, family = binomial, x = TRUE,y = TRUE)
model4 <- glm(kdm_AA_binary ~ age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+alcohol_score+mental_health_score, data = kdm_UK, family = binomial, x = TRUE,y = TRUE)

# Calculate prediction values
kdm_UK$predvalue1 <- predict(model1)
kdm_UK$predvalue2 <- predict(model2)
kdm_UK$predvalue3 <- predict(model3)
kdm_UK$predvalue4 <- predict(model4)

# Fit ROC
modelROC1 <- roc(kdm_UK$kdm_AA_binary, kdm_UK$predvalue1)
auc1 <- auc(modelROC1)
ci1 <- ci(auc1)
modelROC2 <- roc(kdm_UK$kdm_AA_binary, kdm_UK$predvalue2)
auc2 <- auc(modelROC2)
ci2 <- ci(auc2)
modelROC3 <- roc(kdm_UK$kdm_AA_binary, kdm_UK$predvalue3)
auc3 <- auc(modelROC3)
ci3 <- ci(auc3)
modelROC4 <- roc(kdm_UK$kdm_AA_binary, kdm_UK$predvalue4)
auc4 <- auc(modelROC4)
ci4 <- ci(auc4)

# Calculate the difference between two models using the roc.test function in the pROC package
roc_test_result1 <- roc.test(modelROC1, modelROC2)
roc_test_result2 <- roc.test(modelROC1, modelROC3)
roc_test_result3 <- roc.test(modelROC1, modelROC4)


# Fit the logistic model for PhenoAge acceleration
model1 <- glm(pheno_AA_binary ~ age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score, data = pheno_UK, family = binomial, x = TRUE,y = TRUE)
model2 <- glm(pheno_AA_binary ~ age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+mental_health_score, data = pheno_UK, family = binomial, x = TRUE,y = TRUE)
model3 <- glm(pheno_AA_binary ~ age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+alcohol_score, data = pheno_UK, family = binomial, x = TRUE,y = TRUE)
model4 <- glm(pheno_AA_binary ~ age+sex+TDI+ethnicity+education+income+family_history+physical_activity_score+smoke_score+HDI_score+sleep_score+BMI_score+BL_score+BG_score+BP_score+alcohol_score+mental_health_score, data = pheno_UK, family = binomial, x = TRUE,y = TRUE)

# Calculate prediction values
pheno_UK$predvalue1 <- predict(model1)
pheno_UK$predvalue2 <- predict(model2)
pheno_UK$predvalue3 <- predict(model3)
pheno_UK$predvalue4 <- predict(model4)

# Fit ROC
modelROC1 <- roc(pheno_UK$pheno_AA_binary, pheno_UK$predvalue1)
auc1 <- auc(modelROC1)
ci1 <- ci(auc1)
modelROC2 <- roc(pheno_UK$pheno_AA_binary, pheno_UK$predvalue2)
auc2 <- auc(modelROC2)
ci2 <- ci(auc2)
modelROC3 <- roc(pheno_UK$pheno_AA_binary, pheno_UK$predvalue3)
auc3 <- auc(modelROC3)
ci3 <- ci(auc3)
modelROC4 <- roc(pheno_UK$pheno_AA_binary, pheno_UK$predvalue4)
auc4 <- auc(modelROC4)
ci4 <- ci(auc4)

# Calculate the difference between two models using the roc.test function in the pROC package
roc_test_result1 <- roc.test(modelROC1, modelROC2)
roc_test_result2 <- roc.test(modelROC1, modelROC3)
roc_test_result3 <- roc.test(modelROC1, modelROC4)
