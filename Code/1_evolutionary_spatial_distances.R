##################################################################
##                                                              ##
##                 NIPAH VIRUS GENETIC DIVERSITY                ##
##                    SPATIAL SPREAD ANALYSIS                   ##
##                                                              ##
## Last update : 16/03/2021                                     ##
## Oscar Cortes Azuero                                          ##
##                                                              ##  
##################################################################

## Load necessary packages ######################################
library(ape)
library(geodist)
library(ggplot2)

## Load tree and metadata ######################################
tree <- read.nexus("../Data/NiV_172_final.mcc.tre")
metadata <- read.csv("../Data/seq_metadata.csv")

## Figure 2: SPATIAL DISTANCE IN FUNCTION OF EVOLUTIONARY DISTANCE ###########

# Extracting patristic distance
pat_dist <- cophenetic.phylo(tree)
cases <- colnames(pat_dist)
case_id <- stringr::str_split_i(tree$tip.label, "_", 1)
colnames(pat_dist) <- rownames(pat_dist) <- case_id

duplicates <- metadata$accession[!is.na(metadata$equal_to)] #work directly with sequences with loc data
repeated_loc <- metadata$equal_to[!is.na(metadata$equal_to)]

# Replicating data for duplicates + keeping the duplicates' ids
pat_dist_2 <- matrix(data=NA, nrow=length(case_id) + length(duplicates), 
                     ncol=length(case_id) + length(duplicates))
pat_dist_2[1:length(case_id), 1:length(case_id)] <- pat_dist
for (i in 1:length(repeated_loc)){
  pat_dist_2[1:length(case_id),length(case_id)+i] <- pat_dist[,which(case_id == repeated_loc[i])]
}
for (i in 1:length(repeated_loc)){
  pat_dist_2[length(case_id)+i,] <- pat_dist_2[which(case_id == repeated_loc[i]),]
}
colnames(pat_dist_2) <- rownames(pat_dist_2) <- c(case_id, duplicates)
diag(pat_dist_2) <- NA

# Keep 1 seq per roost per year or per outbreak
metadata <- metadata[metadata$unique_sampling_event =="x",] 

# Generate distance matrix using most precise location data available 
loc_data <- data.frame(x=metadata$E, y=metadata$N)
spat_dist <- geodist(loc_data, measure="geodesic")
colnames(spat_dist) <- rownames(spat_dist) <- metadata$accession
diag(spat_dist) <- NA

# Subsetting sequences to include only sequences for which we have both loc AND evol metadata
seq_tre_loc <- intersect(colnames(spat_dist), colnames(pat_dist_2))
spat_dist_subset <- spat_dist[seq_tre_loc, seq_tre_loc]
diag(spat_dist_subset) <- NA
pat_dist_subset <- pat_dist_2[seq_tre_loc, seq_tre_loc]
diag(pat_dist_subset) <- NA

# Function to compute spatial distance in function of evolutionary distance
window_dist <- function(max_window, pat_dist=pat_dist_subset, spat_dist=spat_dist_subset){
  return(mean(spat_dist[which(pat_dist <= max_window, arr.ind = T)]))
}

# Wrapper function to handle subsampling dataset in bootstrapping
bootstrap_dist <- function(seq_list, max_win, pat_dist=pat_dist_subset, spat_dist=spat_dist_subset){
  pat_dist_2 = pat_dist[seq_list, seq_list]
  spat_dist_2 = spat_dist[seq_list, seq_list]
  return(sapply(max_win, FUN=window_dist, pat_dist=pat_dist_2, spat_dist=spat_dist_2))
}

# Bootstrapping
set.seed(199)
seq_lists <- list()
for (i in 1:1000){
  seq_lists[[i]] <- sample(seq_tre_loc, size=length(seq_tre_loc), replace=T)
}

max_windows <- seq(0,120, by=10)

results <- do.call(rbind, lapply(seq_lists, FUN=bootstrap_dist, max_win=max_windows))

res_plot_df <- data.frame(t(apply(results, 2, quantile, c(0.5, 0.025, 0.975), na.rm=T)))
colnames(res_plot_df) <- c("med", "lwr", "upr")
res_plot_df$window <- max_windows

ggplot(res_plot_df, aes(x=window, y=med/1000)) +
  geom_ribbon(aes(ymin=lwr/1000, ymax=upr/1000), color=NA, fill="#5BD1D7", alpha=0.65) +
  geom_line(linewidth=0.75) +  theme_bw(base_size = 14) +
  labs(x='Evolutionary distance inferior to (years)', y='Mean pairwise spatial distance (km)') +
  scale_y_continuous(breaks=c(0, 250, 500, 750, 1000, 1250)) +
  theme(axis.text=element_text(colour="black"), axis.ticks = element_line(colour="black"))
