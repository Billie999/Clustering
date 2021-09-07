# Assignment 2
# title: Clustering gene-expression data
# name: Biljana Simonovikj

library(foreign)
library(dplyr)
library(stats)
library(graphics)
library(factoextra)
library(dbscan)
library(reshape)
library(reshape2)
library(dendextend)
library(tidyverse)


# 1. Read the dataset directly from the ARFF file into a data frame.
golub <- read.arff("/Users/Biljana/Data Mining/Ass 2/golub-1999-v1_database.arff") # read data from Weka Attribute-Relation File Format (ARFF) files into a data frame (golub)
dim(golub) # retrieve dimensions of the data frame
sum(is.na(golub)) # overall number of missing values in the data frame
levels(golub$Classe) # factor levels of the categorical variable Classe

# 2. Set aside the rightmost column (containing the class labels) from the data, storing it separately from the remaining data frame (with the 1868 predictors).
golub_classe <- golub %>% select(-c(1:1868)) # allocate the rightmost variable "Classe" from the data frame by dropping a sequence of variables (1:1868)
golub_predictors <- golub %>% select(c(1:1868)) # store the remaining sequence of predictors variables in a data frame (golub_predictors)

# 3 Use the 72 X 1868 data frame to compute a matrix containing all the pairwise Euclidean distances between observations, that is, a 72 x 72 matrix with distances
# between tissue samples according to their 1868 expression levels. This matrix must be of type dist, which can be achieved either by using the function dist() from
# the base R package stats or by coercion using the function as.dist().
golub_matrix <- dist(golub_predictors, method = "euclidean", diag = FALSE, upper = FALSE, p = 2) # distance matrix computation of class "dist" object by using the
# specified distance measure (euclidean) that computes the dissimilarity distances between the rows of a data matrix
head(as.matrix(golub_matrix)) # conversion to conventional distance matrix

# 4. Use the distance matrix as input to call the Single-Linkage clustering algorithm available from the base R package stats and plot the resulting dendrogram.
# Do not use any class labels to perform this step.
golub_sl <- hclust(golub_matrix, method = "single") # hierarchical cluster analysis with the Single-Linkage clustering algorithm

plot(golub_sl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "red3", labels = FALSE,
     main = "Cluster Dendrogram with Single Linkage Method") # plot the dendogram with Single-Linkage method (see Fig 1)

# 5. # Use the distance matrix as input to call the Complete-Linkage clustering algorithm available from the base R package stats and plot the resulting dendrogram.
# Do not use any class labels to perform this step.
golub_cl <- hclust(golub_matrix, method = "complete") # hierarchical cluster analysis with the Complete-Linkage clustering algorithm

plot(golub_cl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "turquoise3", labels = FALSE,
     main = "Cluster Dendrogram with Complete Linkage Method") # plot the dendrogram with Complete-Linkage method (see Fig 2)

# 6 Use the distance matrix as input to call the Average-Linkage clustering algorithm available from the base R package stats and plot the resulting dendrogram.
# Do not use any class labels to perform this step.
golub_al <- hclust(golub_matrix, method = "average") # hierarchical cluster analysis with the Average-Linkage clustering algorithm

plot(golub_al, xlab = "", sub = "", cex = 0.6, hang = -1, col = "slateblue2", labels = FALSE,
     main = "Cluster Dendrogram with Average Linkage Method") # plot the dendrogram with Average-Linkage method (see Fig 3)

# 7. Use the distance matrix as input to call Ward’s clustering algorithm available from the base R package stats and plot the resulting dendrogram.
# Do not use any class labels to perform this step.
golub_wl <- hclust(golub_matrix, method = "ward.D2") # hierarchical cluster analysis with Ward’s clustering algorithm

plot(golub_wl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "maroon3", labels = FALSE,
     main = "Cluster Dendrogram with Ward 2 Linkage Method") # plot the dendrogram with Ward's 2 Linkage method (see Fig 4)

# 8. Compare the dendrograms plotted in Items 4 to 7. Visually, the dendrograms suggest that some clustering algorithm(s) generate more clear clusters
# than the others.In your opinion, which algorithm(s) may we be referring to and why? In particular, in which aspects do the results produced by this/these
# algorithm(s) look more clear? Perform Item 9 below only for this/those algorithm(s).
# Plot dendrograms with Single, Complete, Average and Ward's 2 Linkage method together:
opar <- par(mfrow = c(2, 2))
plot(golub_sl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "red3", labels = FALSE,
     main = "Cluster Dendrogram with Single Linkage Method")
plot(golub_cl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "turquoise3", labels = FALSE,
     main = "Cluster Dendrogram with Complete Linkage Method")
plot(golub_al, xlab = "", sub = "", cex = 0.6, hang = -1, col = "slateblue2", labels = FALSE,
     main = "Cluster Dendrogram with Average Linkage Method")
plot(golub_wl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "maroon3", labels = FALSE,
     main = "Cluster Dendrogram with Ward 2 Linkage Method")
par(opar)

# 9. Redraw the dendrogram(s) for the selected algorithm(s) in Item 8, now using the class labels that you stored separately in Item 2 to label the observations
# (as disposed along the horizontal axis of the dendrogram). Do some prominent clusters in the dendrogram(s) correspond approximately to the classes
# (that is, the two subtypes of leukemia)?
Classe <- golub_classe$Classe # assigning variable as a vector

# Plot cluster dendrogram with Complete-Linkage method along with class labels:
dendrogram_cl <- as.dendrogram(golub_cl)
par(mar = c(12,4,1,1))
dendrogram_cl %>%
  set("labels_col", value = c("blue", "maroon3"), k = 2) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color", value = c("blue", "maroon3"), k = 2) %>%
  place_labels(paste(golub_classe$Classe, sep = "_")) %>%
  plot( main = "Cluster Dendrogram with Complete Linkage Method")
dendrogram_cl_rect <- rect.dendrogram(dendrogram_cl, k = 2, lty = 5, lwd = 0, x = 1, col = rgb(0.1, 0.2, 0.4, 0.1))

# Plot cluster dendrograms with Ward's 2 Linkage method along with class labels:
dendrogram_wl <- as.dendrogram(golub_wl)
par(mar = c(12,4,1,1))
dendrogram_wl %>%
  set("labels_col", value = c("blue", "maroon3"), k = 2) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color", value = c("blue", "maroon3"), k = 2) %>%
  place_labels(paste(golub_classe$Classe, sep = "_")) %>%
  plot( main = "Cluster Dendrogram with Ward 2 Linkage Method")
dendrogram_wl_rect <- rect.dendrogram(dendrogram_wl, k = 2, lty = 5, lwd = 0, x = 1, col = rgb(0.1, 0.2, 0.4, 0.1))

# 10. Repeat the analysis, now using normalised data. The 1868 predictors have not been normalised before computing the distance matrix in Item 3.
# Normalisation is a non-trivial aspect in unsupervised clustering, as there is no ground truth to assess whether or not it improves performance.
# On the one hand, it may, prevent variables with wider value ranges to dominate distance computations, but on the other hand it, may distort clusters
# by removing natural differences in variance that help characterise them as clusters. Normalisation thus becomes an aspect of Exploratory Data Analysis when it
# comes to clustering: the analyst will usually generate and try to interpret results both with normalised and non-normalised versions of the data. The type of
# normalisation depends on the application in hand. Here, we are computing Euclidean distance between rows of the dataset, so the type of normalisation that applies
# is typically the so-called  z-score normalisation of columns, where each column is rescaled to have zero mean and standard deviation of 1. In this item, you are
# first asked to normalise the data this way before computing the distance matrix in Item 3.Then, repeat Items 4to 9. Does normalisation improve or worsen the results in this dataset?

# Step 3 - Normalize the data and compute the distance matrix with Euclidean method of measure the distances between the rows of a data matrix:
golub_scaled <- scale(golub_predictors) # standardize the variables to have zero mean and standard deviation 1

round(sd(golub_scaled[,45], 0)) # check the function on random selected variables
round(mean(golub_scaled[,1345], 0)) # check the function on random selected variables

golub_matrix_scaled <- dist(golub_scaled, method = "euclidean", diag = FALSE, upper = FALSE, p = 2) # compute distance dissimilarity matrix with Euclidean distance
# measure that computes the distances between the rows of the normalized matrix as an object of class "dist"

# Step 4 to 8 - Hierarchical clustering with Single, Complete, Average and  Ward's 2 Linkage method on normalized variables:
golub_sl_scaled <- hclust(golub_matrix_scaled, method = "single")
golub_cl_scaled <- hclust(golub_matrix_scaled, method = "complete")
golub_al_scaled <- hclust(golub_matrix_scaled, method = "average")
golub_wl_scaled <- hclust(golub_matrix_scaled, method = "ward.D2")

# Plot dendograms with Single, Complete, Average and Ward's 2 Linkage methods:
plot(golub_sl_scaled, xlab = "", sub = "", cex = 0.6, hang = -1, col = "red3", labels = FALSE,
     main = "Single Linkage Method on Normalized Data")
plot(golub_cl_scaled, xlab = "", sub = "", cex = 0.6, hang = -1, col = "turquoise3", labels = FALSE,
     main = "Complete Linkage Method on Normalized Data")
plot(golub_al_scaled, xlab = "", sub = "", cex = 0.6, hang = -1, col = "slateblue2", labels = FALSE,
     main = "Average Linkage Method on Normalized Data")
plot(golub_wl_scaled, xlab = "", sub = "", cex = 0.6, hang = -1, col = "maroon3", labels = FALSE,
     main = "Ward 2 Linkage Method on Normalized Data")

# Step 9 - Plot cluster dendrograms with Complete and Ward's 2 Linkage algorithms applied with class labels that generate the most clear clusters:
# Dendrogram with Complete Linkage algorithm on normalized dataset:
dendrogram_cl_scaled_golub <- as.dendrogram(golub_cl_scaled)
par(mar = c(12,4,1,1))
dendrogram_cl_scaled_golub %>%
  set("labels_col", value = c("blue", "maroon3"), k = 2) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color", value = c("blue", "maroon3"), k = 2) %>%
  place_labels(paste(golub_classe$Classe, sep = "_")) %>%
  plot( main = "Complete Linkage Method on Normalized Data")
dendrogram_cl_rect_scaled <- rect.dendrogram(dendrogram_cl_scaled_golub, k = 2, lty = 5, lwd = 0, x = 1, col = rgb(0.1, 0.2, 0.4, 0.1))

# Dendrograms with Ward 2 Linkage algorithm on normalized dataset:
dendrogram_wl_scaled_golub <- as.dendrogram(golub_wl_scaled)
par(mar = c(12,4,1,1))
dendrogram_wl_scaled_golub %>%
  set("labels_col", value = c("blue", "maroon3"), k = 2) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color", value = c("blue", "maroon3"), k = 2) %>%
  place_labels(paste(golub_classe$Classe, sep = "_")) %>%
  plot( main = "Ward's 2 Linkage Method on Normalized Data")
dendrogram_wl_rect_scaled <- rect.dendrogram(dendrogram_wl_scaled_golub, k = 2, lty = 5, lwd = 0, x = 1, col = rgb(0.1, 0.2, 0.4, 0.1))

# Activity 2: Clustering genes (Part A)
# 11. Read the dataset directly from the ARFF file into a data frame.
yeast <- read.arff("/Users/Biljana/Data Mining/Ass 2/yeast.arff") # read data from Weka Attribute-Relation File Format (ARFF) files into a dataframe
dim(yeast) # retrieve dimensions of the data frame

# 12. Set aside the rightmost column (containing the class labels) from the data, storing it separately from the remaining data frame (with the 20 predictors).
yeast_classe <- yeast %>% select(-c(1:20)) # allocate the variable "Classe" in a yeast_classe data frame by dropping a sequence of variables (1:20)
yeast_predictors <- yeast %>% select(c(1:20)) # storing the remaining sequence of predictors variables (1:20) in data frame
summary(yeast_predictors) # summarize the data set
# 13. Use the 205 X 20 data frame to compute a matrix containing all the pairwise Pearson-based dissimilarities between observations, that is, a 205 X 205
# matrix with dissimilarities between genes according to their 20 expression measurements. Important Note: It is well-known that co-regulated genes are better
# characterised by similar trends in their gene expression profiles, rather than similar expression levels in terms of their absolute values. In other words,
# the similarity between genes in terms of their expression profiles for different measurements is better captured by a correlation measure, such as Pearson correlation
# (James, Witten, Hastie, & Tibshirani, 2013), which is the most widely adopted similarity measure for practical applications of gene clustering. For this reason, in this
# activity we will use Pearson correlation instead of Euclidean distance. However, recall that Pearson is a similarity measure that ranges from -1(lowest similarity) to +1
# (highest similarity). After computing the 205 X 205 Pearson similarity matrix, you have to convert it to a dissimilarity matrix whose values range from 0 (lowest dissimilarity)
# to +1(highest dissimilarity). Once you have this Pearson-based dissimilarity matrix, you can coerce it into type dist as required by the hierarchical clustering methods in the
# base R package stats.
yeast_predictors_cm <- cor(t(yeast_predictors), method = "pearson") # compute Pearson's correlation matrix (similarity matrix)
yeast_predictors_dd <- as.dist((1 - yeast_predictors_cm)/2) # compute Pearson's based dissimilarity distance matrix as object type "dist"
head(as.matrix(yeast_predictors_dd)) # displays first 6 rows and 6 columns of the matrix

# 14. Repeat the clustering analysis in Items 4 to 9 of Activity 1, now using the dissimilarity matrix for the YeastGalactose data computed in Item 13
# (and, when applicable, the class labels that you stored separately in Item 12 to label observations as disposed along the horizontal axis of the relevant dendrograms).
yeast_sl <- hclust(yeast_predictors_dd, method = "single") # hierachical cluster analysis with Single-Linkage method
plot(yeast_sl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "red2", labels = FALSE,
     main = "Single Linkage Method with Correlation Based Distance") # plot the dendogram

yeast_al <- hclust(yeast_predictors_dd, method = "average") # hierachical cluster analysis with Average-Linkage method
plot(yeast_al, xlab = "", sub = "", cex = 0.6, hang = -1, col = "slateblue2", labels = FALSE,
     main = "Average Linkage Method with Correlation Based Distance") # plot the dendogram

yeast_cl <- hclust(yeast_predictors_dd, method = "complete") # hierachical cluster analysis with Complete-Linkage method
plot(yeast_cl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "turquoise3", labels = FALSE,
     main = "Complete Linkage Method with Correlation Based Distance") # plot the dendogram

yeast_wl <- hclust(yeast_predictors_dd, method = "ward.D2") # hierachical cluster analysis with Ward's 2 Linkage method
plot(yeast_wl, xlab = "", sub = "", cex = 0.6, hang = -1, col = "maroon3", labels = FALSE,
     main = "Ward 2 Linkage Method with Correlation Based Distance") # plot the dendogram

# Plot cluster dendograms with Complete, Average and Ward 2 Linkage Method and class labels that generate clear clusters:
Classe <- yeast_classe$Classe # assigning avariable Classe as vector

# Plot cluster dendrogram with Complete-Linkage algorithm and class labels:
dendrogram_cl_yeast <- as.dendrogram(yeast_cl) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color", value = c("black", "blue", "green", "violet"), k = 4)
  colours_to_use_yeast_cl <- as.numeric(yeast_classe$Classe)
  colours_to_use_yeast_cl <- colours_to_use_yeast_cl[order.dendrogram(dendrogram_cl_yeast)]
  labels_colors(dendrogram_cl_yeast) <- colours_to_use_yeast_cl
  dend_list_yeast_cl <- as.character(yeast_classe$Classe)
  labels(dendrogram_cl_yeast) <- dend_list_yeast_cl[order.dendrogram(dendrogram_cl_yeast)]
plot(dendrogram_cl_yeast, main = "Complete Linkage Method with Correlation Based Distance", ylab = "Height")
dendrogram_cl_rect_yeast <- rect.dendrogram(dendrogram_cl_yeast, k = 4, lty = 5, lwd = 0, x = 1, col = rgb(0.1, 0.2, 0.4, 0.1))
legend("topright",
       legend = c("Claster1","Cluster2","Cluster3","Claster4"),
       col = c("black", "red", "green", "blue"),
       title = "Cluster Labels",
       pch = c(20,20), bty = "n", pt.cex = 1.5, cex = 0.8,
       text.col = c("black"), horiz = F, inset = c(0,0.1))

# Plot cluster dendrogram with Ward's 2 Linkage algorithm and class labels:
dendrogram_wl_yeast <- as.dendrogram(yeast_wl) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color", value = c("green", "black", "blue", "red"), k = 4)
colours_to_use_yeast_wl <- as.numeric(yeast_classe$Classe)
colours_to_use_yeast_wl <- colours_to_use_yeast_wl[order.dendrogram(dendrogram_wl_yeast)]
labels_colors(dendrogram_wl_yeast) <- colours_to_use_yeast_wl
dend_list_yeast_wl <- as.character(yeast_classe$Classe)
labels(dendrogram_wl_yeast) <- dend_list_yeast_wl[order.dendrogram(dendrogram_wl_yeast)]
plot(dendrogram_wl_yeast, main = "Ward 2 Linkage Method with Correlation Based Distance", ylab = "Height")
dendrogram_wl_rect_yeast <- rect.dendrogram(dendrogram_wl_yeast, k = 4, lty = 5, lwd = 0, x = 1, col = rgb(0.1, 0.2, 0.4, 0.1))
legend("topright",
       legend = c("Claster1","Cluster2","Cluster3","Claster4"),
       col = c("black", "red", "green", "blue"),
       title = "Cluster Labels",
       pch = c(20,20), bty = "n", pt.cex = 1.5, cex = 0.8,
       text.col = c("black"), horiz = F, inset = c(0,0.1))

# 15. Rescale the data frame in a row-wise fashion so that each rescaled row has magnitude 1.You can achieve this by dividing each element of a row by the magnitude of
#the row.

# Apply min-max normalization to all variables:
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

yeast_predictors_norm = t(apply(yeast_predictors, 1, normalize)) # apply the function to the data frame

# Apply the Euclidean distance normalization so that each row is a unit vector with magnitude one:
magnitude_unit = function(x){
  x/sqrt(sum(x^2))
}

yeast_predictors_scaled = t(apply(yeast_predictors_norm, 1, magnitude_unit)) # apply the function to the data frame

sqrt(sum(yeast_predictors_scaled[100,]^2)) == 1 # check the function with random chosen rows
sqrt(sum(yeast_predictors_scaled[180,]^2)) == 1 # check the function with random chosen rows
sqrt(sum(yeast_predictors_scaled[90,]^2)) == 1 # check the function with random chosen rows

# Another way of checking to see if all rows add up to 1:
result <- sqrt(rowSums(yeast_predictors_scaled^2))
if (mean(result) == 1) {
  print("All rows have a magnitude of 1.")
} else {
  print("Error! Rows are not standardised.")
}

# 16. Run HDBSCAN* (with Euclidean distance) on the rescaled version of the data frame obtained in Item 15. You can (optionally) try different values for the parameter
# MinPts, but MinPts = 5 is required. Plot the resulting HDBSCAN* dendrograms with and without the class labels along the horizontal axis, just like in Items 4–9
# (Activity 1) and Item 14 (Activity 2).
yeast_predictors_hdbs <- hdbscan(yeast_predictors_scaled, minPts = 5) # using single parameter minPts = 5, HDBSCAN finds 4 clusters and
# 24 noise points
yeast_predictors_hdbs # print results of HDBSCAN algorithm

# Condense the complicated cluster hierachy into the simplified cluster tree which shows cluster-wide changes over an infinite number of eps thresholds:
plot(yeast_predictors_hdbs, gradient = c("purple", "blue", "green", "yellow"), show_flat = T) # plot the simplified cluster_tree
color_1to8 <- function(x) ifelse(x == 0,1,((x - 1) %% 7) + 2) # function to apply colors to clusters
plot(yeast_predictors_scaled, pch = 19, col = color_1to8(yeast_predictors_hdbs$cluster), main = "HDBSCAN* Data Set (minpts = 5): 4 clusters", xlab = "x", ylab = "y")
legend("topright",
       legend = c("Outlier","1","2","3","4"),
       col = c("black", "red", "green", "blue", "cyan"),
       pch = c(20,20), pt.cex = 1.5, cex = 0.8,
       text.col = c("black"), horiz = F, inset = c(0.05,0.1)) # plot the extracted partition with clusters in colours and noise in black

# Presenting the labeled 24 outliers(black points) as data points and colored clusters:
outliers <- order(yeast_predictors_hdbs$outlier_scores, decreasing = T)[1:24]
colors <- mapply(function(col, i) adjustcolor(col, alpha.f = yeast_predictors_hdbs$outlier_scores[i]),
                 palette()[yeast_predictors_hdbs$cluster + 1], seq_along(yeast_predictors_hdbs$cluster))
plot(yeast_predictors_scaled, col = colors, pch = 19, main = "HDBSCAN* Data Set (minpts = 5): 4 Clusters with Labelled Outliers", xlab = "x", ylab = "y")
text(yeast_predictors_scaled[outliers, ], labels = outliers, pos = 3)
legend("topright",
       legend = c("Outlier","1","2","3","4"),
       col = c("black", "red", "green", "darkviolet", "cyan"),
       pch = c(20,20), pt.cex = 1.5, cex = 0.8,
       text.col = c("black"), horiz = F, inset = c(0.05,0.1))

# Plot  HDBSCAN* dendrograms with and without the class labels
plot(yeast_predictors_hdbs$hc, main = "HDBSCAN* Hierarchy", xlab = "", sub = "", hang = -1, col = "darkorchid3", labels = FALSE, cex = 0.6) # dendogram without class labels

dend_hdbs <- as.dendrogram(yeast_predictors_hdbs$hc)
colours_to_use_hdbs <- as.numeric(yeast_classe$Classe)
colours_to_use_hdbs <- colours_to_use_hdbs[order.dendrogram(dend_hdbs)]
labels_colors(dend_hdbs) <- colours_to_use_hdbs
dend_list_hdbs <- as.character(yeast_classe$Classe)
labels(dend_hdbs) <- dend_list_hdbs[order.dendrogram(dend_hdbs)]
plot(dend_hdbs, main = "HDBSCAN* Hierarchical Clustering", ylab = "Height")
dendrogram_hdbs_rect_yeast_scaled <- rect.dendrogram(dend_hdbs, k = 4, lty = 5, lwd = 0,
                                                   x = 1, col = rgb(0.1, 0.2, 0.4, 0.1))
legend("topright",
       legend = c("Claster1","Cluster2","Cluster3","Claster4"),
       col = c("black", "red", "green", "blue"),
       title = "Cluster Labels",
       pch = c(20,20), bty = "n", pt.cex = 1.5, cex = 0.8,
       text.col = c("black"), horiz = F, inset = c(0,0.1)) # dendogram with class labels

# From Part 2 activity, I performed out of interest, dendrogram with Single Linkage method on normalized dataset with z-score transformation before Pearson dissimilarity matrix computation:
yeast_predictors_scaled <- scale(yeast_predictors) # normalize the data set to have mean 0 and standard deviation 1
yeast_predictors_cm <- cor(t(yeast_predictors_scaled), method = "pearson") # compute Pearson's correlation matrix (similarity matrix)
yeast_predictors_dd <- as.dist((1 - yeast_predictors_cm)/2) # compute Pearson's based dissimilarity distance matrix as object type "dist"
yeast_sl_norm <- hclust(yeast_predictors_dd, method = "single") # hierachical cluster analysis with Single-Linkage method
plot(yeast_sl_norm, xlab = "", sub = "", cex = 0.6, hang = -1, col = "red2", labels = FALSE,
     main = "Single Linkage Method with Correlation Based Distance on Normalized Dataset") # plot the dendrogram

# 17. Plot a contingency table. By setting MinPts = 5, the automatic cluster extraction method provided by HDBSCAN* extracts four clusters from the resulting hierarchy.
# Plot a contingency table of these clusters (labelled ‘0’, ‘1’, ‘2’, ‘3’ and ‘4’, where ‘0’ means objects left unclustered as noise/outliers) against the
# ground truth class labels that you stored separately in Item 12 (a factor with levels ‘cluster1’, ‘cluster2’,‘cluster3’, ‘cluster4’).
yeast_table <- table(yeast_predictors_hdbs$cluster, labels = yeast_classe$Classe) # plot the contingency table
yeast_table

#19.Plot the genes grouped by their class labels (that is, functional categories ‘cluster1’, ‘cluster2’, ‘cluster3’and ‘cluster4’), in such a way that all the genes
# belonging to the same class are plotted in a separate sub-figure (four sub-figures in total, each one in a different colour). Plot each gene as a time-series with
# 20 data points (where each point is connected by lines to its adjacent points in the series).
Classe <- yeast_classe$Classe

yeast_plot <- cbind(yeast_predictors, cluster = yeast_classe$Classe) # bind the Classe variable to the original yeast data frame
yeast_plot_final <- cbind(yeast_plot, row_number = seq(1, nrow(yeast_plot))) # bind row_number column to the data frame
yeast_plot_melt <- melt(yeast_plot_final, id.vars = c("row_number", "cluster" )) # convert the variables into a moletn data frame with 4 columns identified with two id.variables
colnames(yeast_plot_melt) <- c("Row_number", "Cluster", "Parameter", "Expression") # define the column names of the 4 variables

# Assign the variables as vectors:
Cluster <- yeast_plot_melt$Cluster
Parameter <- yeast_plot_melt$Parameter
Expression <- yeast_plot_melt$Expression
Row_number <- yeast_plot_melt$Row_number

# Create the plot (ggplot2) of gene expressions levels of subset of 205 selected genes of S. cerevisiae from 20 different expreimental conditions as time_series: time series of expression of 4 clusters of genes
ggplot(yeast_plot_melt, aes(x = Parameter, y = Expression, colour = factor(Cluster))) +
  geom_line(aes(group = Row_number)) + facet_grid(Cluster~., space = "free") + ggtitle("Genes Expression Patterns Grouped in Four Functional Categories")

# 20. Plot a figure analogous to the one in Item 19, but now with genes grouped in separate sub-figures according to their cluster as assigned by HDBSCAN* (‘1’, ‘2’, ‘3’ and ‘4’),
# rather than by class labels. Do not plot genes that were left unclustered as noise by HDBSCAN* (labelled ‘0’). Use the best class-tocluster association, as in your answer to Item 18,
# in order to assign each sub-figure of a cluster the same colour used in the sub-figure of the corresponding class in Item 19. For instance, supposing that the best association of class
# ‘clusterX’ in the ground truth is with HDBSCAN* cluster ‘Y’, according to the contingency table in Item 18, then if the genes belonging to class ‘clusterX’ have been plotted in red in
# Item 19, then the genes belonging to HDBSCAN* cluster ‘Y’ should also be plotted in red.

yeast_plot <- cbind(yeast_predictors, cluster = yeast_predictors_hdbs$cluster) # adding in the clusters from HDBSCAN* solution
yeast_plot_final <- cbind(yeast_plot, row_number = seq(1, nrow(yeast_plot))) # bind row_number column to the data frame
yeast_plot_melt <- melt(yeast_plot_final, id.vars = c("row_number", "cluster" )) # convert the variables into a moletn data frame with 4 columns identified with two id.variables
colnames(yeast_plot_melt) <- c("Row_number", "HDBS_Cluster", "Parameter", "Expression") # define the column names of the 4 variables
HDBS_Cluster <- yeast_plot_melt$HDBS_Cluster#  assign the variable as a vectors
sub_yeast_plot <- yeast_plot_melt %>% filter(HDBS_Cluster > 0 ) %>% droplevels() # remove the cluster"0"
sub_yeast_plot$HDBS_Cluster = factor(sub_yeast_plot$HDBS_Cluster, ordered = F,
                                levels = c(4,3,1,2)) # re-ordering the factor levels of HDBSCAN* solution according to the
# ground truth class labels
# Assign the variables as vectors:
Row_number <- sub_yeast_plot$Row_number
HDBS_Cluster <- sub_yeast_plot$HDBS_Cluster
Parameter <- sub_yeast_plot$Parameter
Expression <- sub_yeast_plot$Expression

# Create the plot (ggplot2) of gene expressions levels of subset of 205 selected genes of S. cerevisiae from 20 different expreimental conditions as time_series according to HDBSCAN*
# cluster labels ordered according to the ground truth class labels:
ggplot(sub_yeast_plot, aes(x = Parameter, y = Expression, colour = factor(HDBS_Cluster))) +
  geom_line(aes(group = Row_number)) + facet_grid(HDBS_Cluster~.,space = "free") + ggtitle("Genes Expression Patterns in YeastGalatcose Dataset by HDBSCAN* Algorithm",
          subtitle = "HDBSCAN* generated clusters arranged according to the order of ground truth class labels\n where 4 = cluster1, 3 = cluster2, 1 = cluster3, 2 = cluster4")

