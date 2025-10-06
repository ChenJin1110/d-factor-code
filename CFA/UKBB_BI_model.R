library(lavaan)

UKBBdata<- read.csv("ukbb.csv")

Bi.model <- ' 
M =~
+Conduct
+Depression
+Anxiety
+Phobia
+Pervasive_developmental_disorders
+OCD
+Schizophrenia
+Bipolar
+Psychosis
+Mania
+Adjustment_disorders

P =~ 
+Stroke
+Asthma
+Dermatitis
+Arthritis
+Dorsalgia
+Ulcer
+Colitis
+IBS
+Kidney_ureter_stones
+Hearing
+Visual
+Obesity
+Diabetes1
+Diabetes2
+High_blood_pressure1
+High_blood_pressure2
+Hernia
+Bronchitis
+Movement

D =~
+Conduct
+Depression
+Anxiety
+Phobia
+Pervasive_developmental_disorders
+OCD
+Schizophrenia
+Bipolar
+Psychosis
+Mania
+Adjustment_disorders
+Stroke
+Asthma
+Dermatitis
+Arthritis
+Dorsalgia
+Ulcer
+Colitis
+IBS
+Kidney_ureter_stones
+Hearing
+Visual
+Obesity
+Diabetes1
+Diabetes2
+High_blood_pressure1
+High_blood_pressure2
+Hernia
+Bronchitis
+Movement

# restrictions
M ~~ 0 * D
P ~~ 0 * D
M ~~ 0 * P
'

Uni.model <- ' 
D =~
+Conduct
+Depression
+Anxiety
+Phobia
+Pervasive_developmental_disorders
+OCD
+Schizophrenia
+Bipolar
+Psychosis
+Mania
+Adjustment_disorders
+Stroke
+Asthma
+Dermatitis
+Arthritis
+Dorsalgia
+Ulcer
+Colitis
+IBS
+Kidney_ureter_stones
+Hearing
+Visual
+Obesity
+Diabetes1
+Diabetes2
+High_blood_pressure1
+High_blood_pressure2
+Hernia
+Bronchitis
+Movement
'
Cor.model <- ' 
M =~
+Conduct
+Depression
+Anxiety
+Phobia
+Pervasive_developmental_disorders
+OCD
+Schizophrenia
+Bipolar
+Psychosis
+Mania
+Adjustment_disorders

P =~ 
+Stroke
+Asthma
+Dermatitis
+Arthritis
+Dorsalgia
+Ulcer
+Colitis
+IBS
+Kidney_ureter_stones
+Hearing
+Visual
+Obesity
+Diabetes1
+Diabetes2
+High_blood_pressure1
+High_blood_pressure2
+Hernia
+Bronchitis
+Movement

# restrictions
M ~~ P
'

fit <- cfa(Cor.model, data = UKBBdata, estimator='WLSMV', ordered = TRUE)
fit_uni <- cfa(Uni.model, data = UKBBdata, estimator='WLSMV', ordered = TRUE)
bifactorIndices(fit)

summary(fit,fit.measures=TRUE, rsquare = TRUE,standardized = TRUE)

standardizedSolution(fit)

"------------------------------------------------------------------"
library(psych)

efa_result <- fa(UKBBdata, nfactors = 2, rotate = "varimax", fm = "minres")

print(efa_result)

loadings <- efa_result$loadings
print(loadings)

fa.diagram(efa_result)

"------------------------------------------------------------------"
omega_results <- omega(UKBBdata, nfactors = 3, sl = FALSE, fm = "regression")
print(omega_results)

"------------------------------------------------------------------"
library(semTools)
reliability_results <- reliability(fit)
print(reliability_results)

"------------------------------------------------------------------"
loadings <- inspect(fit, what = "std")$lambda
ave <- colMeans(loadings^2)
print(ave)

"------------------------------------------------------------------"
loadings <- inspect(fit, what = "std")$lambda
theta <- diag(inspect(fit, what = "theta"))
H <- colSums(loadings)^2 / (colSums(loadings)^2 + colSums(theta))
print(H)

"------------------------------------------------------------------"
CR <- colSums(loadings)^2 / (colSums(loadings)^2 + colSums(1 - loadings^2))
print(CR)

"------------------------------------------------------------------"
unstd_loadings <- parameterEstimates(fit, standardized = FALSE)
unstd_loadings <- unstd_loadings[unstd_loadings$op == "=~", ]
print(unstd_loadings)
write.csv(unstd_loadings, file = "unstandardized_loadings.csv", row.names = FALSE)

std_loadings <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
std_loadings <- std_loadings[std_loadings$op == "=~", ]
print(std_loadings)
write.csv(std_loadings, file = "standardized_loadings_with_CI.csv", row.names = FALSE)