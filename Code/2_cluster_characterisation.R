##################################################################
##                                                              ##
##       REGRESSIONS FOR CORRELATION OF GENETIC CLUSTERS        ##
##                                                              ##
## Last update : 08/02/2022                                     ##
## Oscar Cortes Azuero                                          ##
##                                                              ##  
##################################################################

library(ape)
library(dplyr)
library(geodist)
library(ggplot2)


metadata <- read.csv("../Data/seq_metadata.csv")
metadata <- metadata %>%
  filter(!is.na(N), !(short_host %in% c("Human", "Pig", "?Human?", "Dog")))

tree <- ape::read.nexus("../Data/NiV_172_final.mcc.tre")
rownames(metadata) <- metadata$seq_name

# Extract distance, and genetic cluster comparison data for sequences #################
prepare_df <- function(metadata, tree){
  # Belong to the same genetic cluster
  clust_pair <- outer(metadata$cluster, metadata$cluster, "==")
  colnames(clust_pair) <- rownames(clust_pair) <- metadata$accession
  
  # Pairwise spatial distance
  loc_data <- data.frame(x=metadata$E, y=metadata$N)
  spat_dist <- geodist(loc_data, measure="geodesic")
  colnames(spat_dist) <- rownames(spat_dist) <- metadata$accession
  
  # Bat host species
  bat_pair <- outer(metadata$short_host, metadata$short_host, "==")
  colnames(bat_pair) <- rownames(bat_pair) <- metadata$accession
  
  # Sampling time difference
  time_diff <- abs(outer(metadata$dec_date, metadata$dec_date, "-"))
  colnames(time_diff) <- rownames(time_diff) <- metadata$accession
  
  # Extracting patristic distance
  pat_dist <- cophenetic.phylo(tree)
  cases <- colnames(pat_dist)
  case_id <- (do.call(rbind, strsplit(cases, split="_"))[,1])
  colnames(pat_dist) <- rownames(pat_dist) <- case_id
  
  duplicates <- metadata$accession[!is.na(metadata$equal_to)] 
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
  pat_dist_2 <- pat_dist_2[metadata$accession, metadata$accession]
  
  # Names
  pairs <- expand.grid(metadata$accession, metadata$accession)
  pair_names <- matrix(data=paste(pairs[,1], pairs[,2],sep="_"), 
                       nrow=nrow(clust_pair), ncol=ncol(clust_pair))
  colnames(pair_names) <- rownames(pair_names) <- metadata$accession
  
  # Select permutations
  distances <- spat_dist[lower.tri(spat_dist, diag=FALSE)] #data read column-wise
  clusters <- clust_pair[lower.tri(clust_pair, diag = FALSE)]
  hosts <- bat_pair[lower.tri(bat_pair, diag = FALSE)]
  pairs <- pair_names[lower.tri(clust_pair, diag = FALSE)]
  sample_time_diff <- time_diff[lower.tri(time_diff, diag = FALSE)]
  pat_distance <- pat_dist_2[lower.tri(pat_dist_2, diag = FALSE)]
  
  # Prepare output
  seq_data <- data.frame(pair=pairs, dist=distances, clust=clusters, host=hosts, 
                         pat_dist=pat_distance, time_diff=sample_time_diff)
  
  return(seq_data)
}

seq_data <- prepare_df(metadata, tree)
seq_data$same_loc <- seq_data$dist == 0
seq_data$diff_loc <- 1 - seq_data$same_loc

## FIGURE 3C: SAME CLUSTER ACCORDRING TO SAME LOCATION AND DISTANCE ###################
set.seed(199)
clust_model_loc_dist <- glm(clust ~ same_loc + diff_loc:dist, data=seq_data, family = binomial) 

sim_seq_data <- data.frame(dist=seq(0, 2e6, length.out=10000), host=FALSE, same_loc=FALSE)
sim_seq_data$same_loc[1] <- T
sim_seq_data$diff_loc <- 1 - sim_seq_data$same_loc

ilink <- family(clust_model_loc_dist)$linkinv
sim_seq_data$med <- ilink(predict(clust_model_loc_dist, newdata=sim_seq_data))

# Visualise model estimates
coeffs_conf_int <- as.data.frame(t(confint(clust_model_loc_dist, method="boot")))
eta <- coeffs_conf_int$`(Intercept)`[1] + coeffs_conf_int$same_locTRUE[1]*sim_seq_data$same_loc + coeffs_conf_int$`diff_loc:dist`[1]*sim_seq_data$diff_loc*sim_seq_data$dist
sim_seq_data$lwr <- ilink(eta)
eta <- coeffs_conf_int$`(Intercept)`[2] + coeffs_conf_int$same_locTRUE[2]*sim_seq_data$same_loc + coeffs_conf_int$`diff_loc:dist`[2]*sim_seq_data$diff_loc*sim_seq_data$dist
sim_seq_data$lwr <- ilink(eta)

# Visualise data
loc_data <- data.frame(x=metadata$E, y=metadata$N)
spat_dist <- geodist(loc_data, measure="geodesic")
colnames(spat_dist) <- rownames(spat_dist) <- metadata$accession
diag(spat_dist) <- NA

cluster_pair <- outer(metadata$cluster, metadata$cluster, "==")
colnames(cluster_pair) <- rownames(cluster_pair) <- metadata$accession
diag(cluster_pair) <- NA

bootstrap_dist <- function(seq_list, max_win, spat_dist_bs=spat_dist, cluster_pa_bs=cluster_pair){
  cluster_pa_2 = cluster_pa_bs[seq_list, seq_list]
  spat_dist_2 = spat_dist_bs[seq_list, seq_list]
  return(sapply(max_win, FUN=window_dist_cl, spat_dist_=spat_dist_2, cluster_pa_=cluster_pa_2))
}

window_dist_cl <- function(window, spat_dist_=spat_dist, cluster_pa_=cluster_pair){ #, bat_pa_=bat_pair){
  return(mean(cluster_pa_[which((((spat_dist_ >= window[1]) & (spat_dist_ <= window[2]))), arr.ind = T)], na.rm=TRUE))
}

set.seed(199)
seq_lists_2 <- list()
for (i in 1:1000){
  seq_lists_2[[i]] <- sample(metadata$accession, size=length(metadata$accession), replace=T)
}

Pmax <- seq(100, 2250, 250)*1e3
Pmin <- Pmax - 1000*1e3 
Pmin[Pmin <=0] <- 50e3
Pmid<-(Pmin+Pmax)/2
win_lab <- paste(round(Pmin/1000), round(Pmax/1000), sep="-")
windows <- list()
for (i in 1:length(Pmin)){
  windows[[i]] <- c(Pmin[i], Pmax[i])
}

results <- do.call(rbind, lapply(seq_lists_2, FUN=bootstrap_dist, max_win=windows))

res_plot_1 <- data.frame(window=c(1:length(windows)), inf=NA, med=NA, max=NA)
for (i in 1:length(windows)){
  a <- quantile(results[,i], probs = c(0.025, 0.975), na.rm = TRUE)
  res_plot_1[i, 2] <- a[[1]]
  res_plot_1[i, 3] <- median(results[,i], na.rm=T)
  res_plot_1[i, 4] <- a[[2]]
}
res_plot_1 <- cbind(res_plot_1, lab=win_lab, x=Pmid)

result <- do.call(rbind, lapply(seq_lists_2, FUN=bootstrap_dist, max_win=list(c(0, 0))))
roost_result <- data.frame(x=0, min=NA, med=NA, max=NA)
a <- quantile(result, probs = c(0.025, 0.975), na.rm = TRUE)
roost_result[1, 2] <- a[[1]]
roost_result[1, 3] <- median(result, na.rm=T)
roost_result[1, 4] <- a[[2]]

## Prob same cluster based on same_loc and distance
p <- ggplot(sim_seq_data, aes(x=dist/1000)) +
  geom_ribbon(data = sim_seq_data[sim_seq_data$dist > 0,], aes(ymin=lwr, ymax=upr), alpha=1, fill="#80cdc1") +
  geom_line(data = sim_seq_data[sim_seq_data$dist > 0,], aes(x=dist/1000, y=med), 
            size=0.5, color="#018571") +
  geom_point(inherit.aes = F, data = sim_seq_data[sim_seq_data$same_loc == T,], 
             aes(x=dist/1000, y=med), color="#80cdc1", size=4, shape=15) + 
  geom_errorbar(inherit.aes = F, data = sim_seq_data[(sim_seq_data$same_loc == T),], 
                aes(x=dist/1000, ymin=lwr, ymax=upr), color="#80cdc1",width=0, size=1.5) +
  geom_point(data=res_plot_1, aes(x=x/1e3, y=med), size=2, color="#018571") +
  geom_errorbar(data=res_plot_1, aes(x=x/1e3, ymin=inf, ymax=max), width=0, size=0.75, color="#018571") +
  geom_point(data=roost_result, aes(x=x, y=med), color="#018571", size=2) +
  geom_errorbar(data=roost_result, aes(x=x, ymin=min, ymax=max), color="#018571", width=0, size=0.75 ) +
  theme_classic(base_size = 12) + labs(x="Distance (km)", y="Probabilty same cluster") +
  scale_x_continuous(breaks = seq(0, 2000, 500), labels = c("Same loc.", "500", "1000", "1500", "2000")) +
  theme(legend.position = c(0.8, 0.9))
p

## FIGURE 3D: SAME CLUSTER ACCORDING TO HOST SPECIES ####################
clust_model_host <- glm(clust ~ host, data=seq_data, family = binomial)

sim_seq_data <- data.frame(host=c(TRUE, FALSE))
ilink <- family(clust_model_host)$linkinv
sim_seq_data$med <- ilink(predict(clust_model_host, newdata=sim_seq_data))

coeffs_conf_int <- as.data.frame(t(confint(clust_model_host, method="boot")))
eta <- coeffs_conf_int$`(Intercept)`[2] + coeffs_conf_int$hostTRUE[2]*sim_seq_data$host
sim_seq_data$upr <- ilink(eta)
eta <- coeffs_conf_int$`(Intercept)`[1] + coeffs_conf_int$hostTRUE[1]*sim_seq_data$host
sim_seq_data$lwr <- ilink(eta)

## Computing odds ratio for cluster~host
odds <- function(p){
  return(p/(1-p))
}

odds_model <- apply(sim_seq_data[,c(2:4)], MARGIN=2, FUN=odds)
odds_ratio <- odds_model[1,]/odds_model[2,]

p <- ggplot(as.data.frame(t(odds_ratio)), aes(x=1, y=med)) +
  geom_point(size=2, colour="#201887") +
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0, size=0.75, colour="#201887") +
  scale_x_continuous(breaks=c(1)) +
  theme_classic(base_size = 12) + 
  labs(x="Same host species", y="OR")
p

## FIGURE S5: SAME HOST ACCORDING TO DISTANCE ####
host_model_dist <- glm(host ~ dist, data=seq_data, family = binomial)

sim_seq_data <- data.frame(dist=seq(0, 2e6, length.out=10000))
ilink <- family(host_model_dist)$linkinv
sim_seq_data$med <- ilink(predict(host_model_dist, newdata=sim_seq_data))

coeffs_conf_int <- as.data.frame(t(confint(host_model_dist, method="boot")))
eta <- coeffs_conf_int$`(Intercept)`[1] +  coeffs_conf_int$dist[1]*sim_seq_data$dist
sim_seq_data$lwr <- ilink(eta)
eta <- coeffs_conf_int$`(Intercept)`[2] +  coeffs_conf_int$dist[2]*sim_seq_data$dist
sim_seq_data$upr <- ilink(eta)

bootstrap_dist <- function(seq_list, max_win, spat_dist_bs=spat_dist, bat_pa_bs=bat_pair){
  bat_pa_2 = bat_pa_bs[seq_list, seq_list]
  spat_dist_2 = spat_dist_bs[seq_list, seq_list]
  return(sapply(max_win, FUN=window_dist_cl, spat_dist_=spat_dist_2, bat_pa_=bat_pa_2))
}

window_dist_cl <- function(window, spat_dist_=spat_dist, bat_pa_=bat_pair){
  return(mean(bat_pa_[which((((spat_dist_ >= window[1]) & (spat_dist_ <= window[2]))), arr.ind = T)], na.rm=TRUE))
}

windows[[1]][1] <- 0
results <- do.call(rbind, lapply(seq_lists_2, FUN=bootstrap_dist, max_win=windows))

results_plot <- data.frame(window=c(1:length(windows)), lwr=NA, med=NA, upr=NA)
for (i in 1:length(windows)){
  a <- quantile(results[,i], probs = c(0.025, 0.975), na.rm = TRUE)
  results_plot[i, 2] <- a[[1]]
  results_plot[i, 3] <- median(results[,i], na.rm=T)
  results_plot[i, 4] <- a[[2]]
}

p <- ggplot(sim_seq_data, aes(x=dist/1000)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3, fill="#a6611a") +
  geom_line(aes(y=med), size=1, color="#a6611a") +
  geom_point(data=results_plot, aes(x=Pmid/1e3, y=med), size=3, color="#a6611a") +
  geom_linerange(data=results_plot, aes(x=Pmid/1e3, ymin=lwr, ymax=upr), linewidth=0.8, color="#a6611a") +
  theme_bw(base_size = 14) + labs(x="Distance (km)", y="Probabilty same host species") +
  theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour="black"))
p
