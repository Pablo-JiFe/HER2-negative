# In this file we run the Gaussian mixture model and analyze the clusters


library(mclust)
library(fpc)
library(cluster)


# 1.- GMM -----------------------------------------------------------------

# 1.1 Extract ERBB2 values

erbb2_vals <- counts_her2_low["ERBB2", ] %>% 
  as.numeric()

# 1.1.2 Assign names to the counts  

names(erbb2_vals) <- colnames(counts_her2_low)


# 1.2 Gaussian Mixture Model

set.seed(123)
fit <- Mclust(erbb2_vals, G = 3)

fit$BIC # Measure the difference between clusters so as to decide the best
plot(fit, what = "BIC")

# 1.3 Check for separation with the silhouette score

sil <- silhouette(fit$classification, dist(erbb2_vals))
mean(sil[, 3])   # average silhouette score

# 1.4 Intersection and stability

clusterboot(erbb2_vals, clustermethod = kmeansCBI, k = fit$G)



# 2.- Obvserve the clusters -----------------------------------------------


# 2.1 Assign labels to metadata corresponding to thee clusters identified in the GSM

# 2.1.1 Make sure that the order of patients in metadata and the fit values are thse same so as to assign the clusters correctly

metadata_her2_class <- metadata_her2_low[match(rownames(fit$data), metadata_her2_low$title), ]

# 2.1.2 Assing the labels to its respective objects

metadata_her2_class <- metadata_her2_class %>%
  mutate(Mol_Status = case_when(
    fit$classification == 1 ~ "Mol_zero",
    fit$classification == 2 ~ "Mol_ultralow",
    fit$classification == 3 ~ "Mol_low"
  ))



# 2.2 Visualize the split

ggplot(metadata_her2_class, aes(x = erbb2_vals, fill = Mol_Status)) +
  geom_density(alpha = 0.5) +
  labs(title = "GGM clustering based on ERBB2 expression",
       x = "ERBB2 Expression (Log2 FPKM)", y = "Density") +
  theme(text = element_text(size = 30))

# 2.3 Create a plot to see if ER status explains the three peaks

ggplot(metadata_her2_class, aes(x = erbb2_vals, fill = er_status)) + 
  geom_density(alpha = 0.5) +
  facet_wrap(~Mol_Status) +
  labs(title = "Is the Low/Zero split just Estrogen Status?") +
  theme(text = element_text(size = 30))


# 2.4 Create a plot to see if ER status explains the two peaks

ggplot(metadata_her2_class, aes(x = erbb2_vals, fill = Mol_Status)) + 
  geom_density(alpha = 0.5) +
  facet_wrap(~er_status) +
  labs(title = "HER2 split dependence on clinical ER status") +
  theme(text = element_text(size = 30))

# 2.4.2 Create a plot to see if ER status explains the two peaks, similñar but with the SCAN-B classification

ggplot(metadata_her2_class, aes(x = erbb2_vals, fill = Mol_Status)) + 
  geom_density(alpha = 0.5) +
  facet_wrap(~er_pred_sgc) +
  labs(title = "HER2 split dependence on single gene classifier prediciton of ER status") +
  theme(text = element_text(size = 30))


# 3.- Local minimum -------------------------------------------------------

# 3.1 Get the expression values for only the Mol_zero group sinceits the one wiuth the bimodal structure

zero_vals <- erbb2_vals[fit$classification == 1]

# 3.2 Compute the kernel density estimate

dens <- density(zero_vals)

# 3.3 Find the local minimum between the two peaks

# 3.3.1 We want to avoid the vminimum being selected in the extremes of the curve since we want the division that splits it into 2

valley_region <- which(dens$x > 3.0 & dens$x < 7.5)

# 3.3.2 We look for the point where density is lowest 

min_idx <- valley_region[which.min(dens$y[valley_region])]

# 3.3.3 Obtain the value

valley_x <- dens$x[min_idx]

# 3.4 Metadata with the HER2-zero division

metadata_her2_class <- metadata_her2_class %>%
  mutate(Extended_Status = case_when(
    Mol_Status == "Mol_zero" & erbb2_vals > valley_x ~ "Zero_High",
    Mol_Status == "Mol_zero" ~ "True_Zero",
    Mol_Status == "Mol_ultralow" ~ "Ultralow",
    Mol_Status == "Mol_low" ~ "Low"
  ))

# 3.4.2 Check how many 'Mismatched' patients you found
table(metadata_her2_class$Extended_Status)


#> //Dictionary//#######################################################
#> 
#> metadata_her2_class <- Here metadata has now columns with the clusters
#> 
