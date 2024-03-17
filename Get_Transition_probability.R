#filter the grids by number of cells
library(ks)
kde(freq_grid$Freq)
median(freq_grid$Freq)
dim(freq_grid[freq_grid$Freq > 10, ])
filtered_freq_grid <- freq_grid[freq_grid$Freq > 10, ]
filtered_grid_cell <- grid_cell[grid_cell$gridGroup %in% filtered_freq_grid$Var1, ]

grid_expmat <- list()
for (z in 1:length(filtered_freq_grid$Var1)){
  print(z)
  temp_grid <- filtered_grid_cell[filtered_grid_cell$gridGroup == as.character(filtered_freq_grid$Var1[z]), ]
  subset_genecellmat <- t_genecell_mat[rownames(t_genecell_mat) %in% temp_grid$cell_type,]
  grid_expmat[[z]] <- subset_genecellmat
}

#find average expression per grid
avg_expmat_pergrid <- data.frame(matrix(nrow = 313, ncol = length(grid_expmat)))
for (i in 1:length(grid_expmat)){
  avg_expmat_pergrid[,i] <- unname(unlist(colSums(grid_expmat[[i]])))
}
colnames(avg_expmat_pergrid) <- filtered_freq_grid$Var1
rownames(avg_expmat_pergrid) <- rownames(genecell_mat)

#implementing KNN and calinski harabasz score
pca <- prcomp(t(avg_expmat_pergrid), scale=F, center=T)
pca_importance <- summary(pca)
plot(pca_importance$importance[2,])
coord_pca <- data.frame(pca$x)
coord_pca_pcs <- dplyr::select(coord_pca, c(paste0("PC", as.character(seq(1, 10)))))
dist_coord_pca_pcs <- as.matrix(dist(coord_pca_pcs))

#calculate the distance between grids in KNN space and divided it by actual distance (average centroids from cells in grids and apply dist function)

coords_arranged <- coords[match(grid_cell$cell_type,
                                coords$cell_id), ]
all(coords_arranged$cell_id == grid_cell$cell_type)
coords_arranged$grid <- grid_cell$gridGroup
coords_arranged <- coords_arranged[coords_arranged$grid %in% colnames(avg_expmat_pergrid), ]

mean_x <- aggregate(coords_arranged$x_centroid, list(coords_arranged$grid), FUN=mean)
mean_y <- aggregate(coords_arranged$y_centroid, list(coords_arranged$grid), FUN=mean)[,2]
mean_total <- data.frame("mean_x" = mean_x, "mean_y" = mean_y)
rownames(mean_total) <- mean_total$mean_x.Group.1
mean_total <- dplyr::select(mean_total, -"mean_x.Group.1")
dist_total <- as.matrix(dist(mean_total))
dim(dist_total)

#calculate 
dist_coord_pca_pcs <- replace(dist_coord_pca_pcs, dist_coord_pca_pcs == 0 , 1)
reci_dist_coord_pca_pcs <- 1/dist_coord_pca_pcs
diag(reci_dist_coord_pca_pcs) <- 0
dist_coord_pca_pcs <- dist_coord_pca_pcs[rownames(dist_coord_pca_pcs) %in% names(emission_probabilities),
                                         colnames(dist_coord_pca_pcs) %in% names(emission_probabilities)]
dist_total <-  replace(dist_total, dist_total == 0 , 1)
reci_dist_total <- 1/dist_total
diag(reci_dist_total) <- 0
reci_dist_total <- reci_dist_total[rownames(reci_dist_total) %in% names(emission_probabilities),
                                   colnames(reci_dist_total) %in% names(emission_probabilities)]
#addition of individual similarity indices
scaled <- scale(reci_dist_coord_pca_pcs) + scale(reci_dist_total)
minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}


#kmeans (calculate centroid distances within kmeans, 
BiocManager::install("fpc")
library(fpc)
library(ggfortify)

normalizeddata 

BiocManager::install("clValid")
library(clValid)
scaled #addition of individual similarity indices (calculate centroid distances before min max normalization)
highest_ratio <- c() #highest ratio desired
k_cluster_assignments <- list()
centroid_distances <- list()
for (i in 2: 10){
  print(i)
  km.res <- kmeans(scaled, i, nstart = 25)
  cluster_assignments <- km.res$cluster
  k_cluster_assignments[[i]] <- cluster_assignments
  temp_mat <- as.matrix(dist(km.res$centers))
  centroid_distances[[i]] <- temp_mat
  highest_ratio <- c(highest_ratio, km.res$betweenss/km.res$totss)
}

#getting cluster assignments and identifying ideal centroid distances

df <- as.data.frame(k_cluster_assignments[[10]]) #preliminary index calculation based on wss
min_max_scaled <- as.data.frame(apply(centroid_distances[[10]], MARGIN = 2, minMax))
min_max_scaled <- 1 - min_max_scaled
write.csv(min_max_scaled, "1732023_preliminary_transition_probabilities.csv")
write.csv(df, "1732023_cluster_assignments_transition_probabilities.csv")