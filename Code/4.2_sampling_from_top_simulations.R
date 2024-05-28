##################################################################
##                                                              ##
##            SAMPLING FROM TOP SIMULATIONS - ABC               ##
##                                                              ##
##                                                              ##
## Last update : 28/06/2023                                     ##
## Oscar Cortes Azuero                                          ##
##################################################################

# Load required packages
library(sf)
library(doParallel)
library(foreach)

sampling_circle <- function(i, radius_sample_circle, width_polygon, circles_){
  # Creates a circle of radius radius_sample_circle as a proxy of the area of circulation of
  # Pteropus bats in a country, and counts how many circles from the circles_ object are 
  # intersected by this sampled circle. The circles_ object represents the genetic clusters
  # in a specific area
  polygon =
    list(
      matrix(c(0, 0, width_polygon-radius_sample_circle, 0, width_polygon-radius_sample_circle, width_polygon-radius_sample_circle, 0, width_polygon-radius_sample_circle, 0, 0), ncol=2, byrow=T)
    )
  polygon = sf::st_polygon(polygon)
  center <- st_sample(polygon, size=1)
  sample_circle <- st_buffer(center, radius_sample_circle)
  inters <- st_intersects(sample_circle, circles_)
  return(length(inters[[1]]))
}

rep_sampling_circle <- function(radius_sample_circle, n, width_polygon, circles_){
  # Repeats sampling_circle over n iterations
  return(sapply(1:n, sampling_circle, radius_sample_circle, width_polygon, circles_))
}

n_circs_inters <- function(x, areas=areas_to_sample){
  # Within a simulated square, simulate n_clust genetic clusters of Nipah with a standard deviation
  # of stdev, and for the areas of countries with Pteropus bat circulation, bootstrap the number
  # of clusters that could circulate in each specific area over 1000 iterations.
  width <- 11000
  radii_sample_circles <- sqrt(areas/pi)
  
  n_clust <- x[1]
  stdev <- x[2]
  
  polygon =
    list(
      matrix(
        c(0, 0, width, 0, width, width, 0, width, 0, 0),
        ncol=2, byrow=T
      )
    ) 
  polygon = sf::st_polygon(polygon)
  
  centers <- st_sample(polygon, size=n_clust)
  centers <- st_sf(centers)
  centers$cluster <- c(1:n_clust)
  centers$cluster <- as.factor(centers$cluster)
  
  circles <- st_buffer(centers, dist=2*stdev)
  
  n_intersections_df <- as.data.frame(sapply(radii_sample_circles, rep_sampling_circle, 1000, width, circles))
  colnames(n_intersections_df) <- paste0("a", areas)
  
  return(n_intersections_df)
}

top_sim_params <- as.matrix(read.csv("top_2pc_sim_params.csv")[,1:2])
coverage <- read.csv("../Data/tree_pteropus_mean_coverage__thresholds_10_30_50_all_countries.csv")
areas_to_sample <- coverage$extent

registerDoParallel(52)
results <- foreach (i=1:(nrow(top_sim_params)), .combine=rbind) %dopar% {n_circs_inters(top_sim_params[i,])}
write.csv(results, "sampling_results_all_countries.csv", row.names=F)
stopImplicitCluster()




