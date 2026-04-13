# In this file we compare the survival and phenotypic differences of the clusters


library(survival)
library(survminer)


# 1.- Survival object -----------------------------------------------------

# 1.1 Create the survival object

s_obj <- Surv(time = metadata_her2_class$os_months, event = metadata_her2_class$os_status)

metadata_her2_class$surv_obj <- s_obj


# 2.- Survival based on original cluster ----------------------------------

# 2.1 Fit the curves based on the initial clusters

fit_surv <- survfit(s_obj ~ Mol_Status, data = metadata_her2_class)

# 2.2 Cox survival on initial clusters

fit_cox <- coxph(s_obj ~ Mol_Status, data = metadata_her2_class)

summary(fit_cox)


# 2.3 Kaplan-Meier Curve

ggsurvplot(fit_surv, 
           data = metadata_her2_class,
           pval = TRUE,             # Does the difference matter statistically?
           risk.table = TRUE,       # Show how many patients are left over time
           conf.int = FALSE,        # Keeps the plot clean for 4 groups
           palette = c("red", "darkblue", "green", "lightblue"),
           title = "Survival: True Zero vs. High-RNA Zero vs. Low")



# 3.- Phenotype data based on initial clusters ---------------------------

# 3.1 Create a table of clusters with its count of PAM50 sub types

subtype_check <- table(metadata_her2_class$Mol_Status, metadata_her2_class$pam50)

# 3.1.2 Create a table of clusters with its count of PAM50 sub types and the ER status

subtype_check <- table(metadata_her2_class$Mol_Status, metadata_her2_class$pam50, metadata_her2_class$er_pred_sgc)


# 3.2 Convert to proportions and modify to make it easier to read

# 3.2.1 Independently of ER status

prop.table(subtype_check, 1) %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  mutate(across(-Var1, ~ . * 100)) %>% 
  mutate(Luminal = LumA + LumB) # Sum proportions of luminal subtype

# 3.2.2 Dependance on ER status

metadata_her2_class %>% 
  dplyr::count(er_pred_sgc, 
               Mol_Status, 
               pam50) %>%
  group_by(er_pred_sgc, Mol_Status) %>% # Group so as to anbalyze phenotype based on this 2 factors
  mutate(prop = 100 * n / sum(n)) %>% 
  select(-n) %>% 
  pivot_wider(names_from = pam50, values_from = prop, values_fill = 0) %>% 
  mutate(Luminal = LumA + LumB)


###############################################################################
###############################################################################
#---------------------------//Cluster corrected//------------------------------
###############################################################################
###############################################################################


# 4.- Survival analysis based on corrected cluster ------------------------


# 4.1 Fit the curves based on the corrected clusters

fit_surv <- survfit(s_obj ~ Extended_Status, data = metadata_her2_class)

# 4.2 Cox survival on initial clusters

fit_cox <- coxph(s_obj ~ Extended_Status, data = metadata_her2_class)

summary(fit_cox)

# 4.3 Kaplan-Meier Curves

ggsurvplot(fit_surv, 
           data = metadata_her2_class,
           pval = TRUE,             
           risk.table = TRUE,       
           conf.int = FALSE,        
           palette = c("red", "darkblue", "green", "lightblue"),
           title = "Survival: True Zero vs. High-RNA Zero vs. Low")



# 5.- Phenotype data based on initial clusters ---------------------------

# 5.1 Create a table of clusters with its count of PAM50 sub types

subtype_check <- table(metadata_her2_class$Extended_Status, metadata_her2_class$pam50) 

# 5.2 Convert to proportions and modify to make it easier to read

# 5.2.1 Independently of ER status

prop.table(subtype_check, 1) %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  mutate(across(-Var1, ~ . * 100)) %>% 
  mutate(Luminal = LumA + LumB)


# 5.2.3 Dependance on ER status

metadata_her2_class %>% 
  dplyr::count(er_pred_sgc, 
               Extended_Status, 
               pam50) %>%
  group_by(er_pred_sgc, Extended_Status) %>% 
  mutate(prop = 100 * n / sum(n)) %>% 
  select(-n) %>% 
  pivot_wider(names_from = pam50, values_from = prop, values_fill = 0) %>% 
  mutate(Luminal = LumA + LumB)



##############################################################################
#--------------------Survival dependency on ER status ------------------------
##############################################################################

# 6.1 Create objects with only one ER status 

# 6.1.1 Only ER+

metadata_her2_class_er_pos <- 
  metadata_her2_class %>% 
  filter(er_pred_sgc == 1)

# 6.1.2 Only ER-

metadata_her2_class_er_neg <- 
  metadata_her2_class %>% 
  filter(er_pred_sgc == 0)

# Object to use for cox

cox_er_dep <- metadata_her2_class_er_neg

# 6.2 Survival object

s_obj <- Surv(time = cox_er_dep$os_months, event = cox_er_dep$os_status)

# 6.3 Fit the curves based on corrected clusters and in only one ER status

fit_surv <- survfit(s_obj ~ Extended_Status, data = cox_er_dep)

fit_cox <- coxph(s_obj ~ Extended_Status, data = cox_er_dep)

fit_cox_sum <- summary(fit_cox)


# 5. Kaplan-Meier Curves

ggsurvplot(fit_surv, 
           data = cox_er_dep,
           pval = TRUE,             # Does the difference matter statistically?
           risk.table = TRUE,       # Show how many patients are left over time
           conf.int = FALSE,        # Keeps the plot clean for 4 groups
           palette = c("red", "darkblue", "green", "lightblue"),
           title = "Survival: True Zero vs. High-RNA Zero vs. Low")

# Results 

for (i in 1:3) {
  
  hr <- fit_cox_sum$coefficients[i,2]
  pval <- fit_cox_sum$coefficients[i,5]
  low_ci <- fit_cox_sum$conf.int[i, 3]
  high_ci <- fit_cox_sum$conf.int[i, 4]
  obj <- fit_cox$xlevels$Extended_Status[i + 1]
  
  cat((paste0("HER2-", obj, " had a HR of ", round(hr, 2), " (CI 95%: ", round(low_ci, 2), " - ", round(high_ci, 2), ", pval = ", round(pval, 5), ")")), ", ")
  
}


