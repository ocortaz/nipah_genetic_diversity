# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 08:15:25 2023

@author: Oscar
"""

import numpy as np
from scipy.stats import truncnorm
from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
import os
import time

def mooreNeighbourhood(coords, a=0, b=7):
    """
    Takes in the column and the row of a cell, returns the Moore Neighbourhood of the cell
        
    Args:
    -----------
    coords: list
        List of length 2, containing the column and the row of a cell within the gridspace, respectively
    a: int
        Value of the smallest possible coordinate for each cell
    b: int
        Value of the largest possible coordinate for each cell

    Returns:
    --------
    Numpy array containing the numbers of all cells in the Moore Neighbourhood of the input cell, including the cell itself.
    A cell in row m and column n will have the number n*8 + m. Column and row numbers start at 0.
    """
    [x, y] = coords
    neighbs_x = np.unique([max(a, x-1), x, min((x+1), b)]).astype(int)
    neighbs_y = np.unique([max(a, y-1), y, min((y+1), b)]).astype(int)
    return(np.array([(b+1)*w+z for w in neighbs_x for z in neighbs_y]))

def searchNNClusts(sel_samples):
    """
    For each point within a given cell of the gridspace, this function finds all points within a radius of 510km of said point, 
    attributing them to their corresponding spatial window, and comparing their genetic clusters.
    Args:
    -----
    sel_samples: tuple
        Tuple of length 2. The first element is a subset of the simulated samples: each sample is an array of length 4 (exact coordinates, genetic cluster, cell number).
        The second element is the cell of interest. The subset samples in the first element are all points in the cell or in its Moore Neighbourhood.

    Returns:
    --------
    Numby array of shape (6,2). Each row corresponds to a spatial window (0-10km, 10-110km, ..., 410-510km). For each row, the first element corresponds to the number
    of pairs of points within that spatial window that belong to the same cluster, and the second element corresponds to the total number of pairs of points within
    that window. 
    """
    samples_neighbours, cell = sel_samples[0], sel_samples[1]
    samples_cell = samples_neighbours[np.where(samples_neighbours[:,3]==cell)[0],:2] # Get samples within the cell of interest only
    knn_obj = NearestNeighbors(radius=510, algorithm='kd_tree').fit(samples_neighbours[:,:2])
    results = np.array(knn_obj.radius_neighbors(samples_cell, sort_results=True), dtype=object) # Compare points from cell to all points in neighbourhood
    clust_pair_tmp = np.concatenate([samples_neighbours[:,2][results[1][i][1:]] == samples_neighbours[:,2][results[1][i][0]] for i in range(len(results[1]))]) # Compare clusters of all points
    dist_pair_tmp = np.concatenate([np.delete(results[0][i], 0) for i in range(len(results[0]))]) # For each point, delete comparisons to itself
    dist_pair_tmp = (dist_pair_tmp - 10)//100
    return(np.array([[np.nansum(clust_pair_tmp, where=(dist_pair_tmp == win)),np.count_nonzero((dist_pair_tmp==win))] for win in range(-1,5)]))

def simClustSample(seed, min_n_clust, max_n_clust, min_r, max_r, n_obs, n_cpus):
    """
    Runs a full simulation for ABC framework.
    Args:
    -----
        seed: int
            Seed
        min_n_clust: int
            Minimum number of clusters that can be simulated.
        max_n_clust: int
            Maximum number of clusters that can be simulated.
        min_r: int
            Minimum possible value for the standard deviation of the truncated normal distribution that will determine the distance between simulated points 
            and the "center" of their genetic cluster.
        max_r: int
            Maximum possible value for the standard deviation of the truncated normal distribution that will determine the distance between simulated points 
            and the "center" of their genetic cluster.
        n_obs: int
            Number of points that will be simulated for each genetic cluster.
        n_cpus: int
            Number of CPUs that will be used to parallelize the computation of summary statistics for the simulation.

    Returns:
    --------
        Numpy array containing:
            For each spatial window, the proportion of pairs of simulated points within that window that belong to the same cluster
            The number of clusters that were simulated
            The standard deviation used to sample the distance between simulated points and the center of their genetic cluster
            The number of clusters for which at least one simulated point was contained within the sampling area.
    """
    rng = np.random.default_rng(seed)
    [n_clust] = rng.integers(min_n_clust, max_n_clust, 1)
    [r] = rng.uniform(min_r, max_r, 1)
    buffer = 3e3 # Padding surrounding the area for which summary statistics will be computed (to avoid edge effects)
    grid_width = 625 # Cell width and height for the grid that will be used to distribute the computation of summary statistics 
    sampling_width = 5e3 # Total width of the area for which summary statistics will be computed
    area_width = 2*buffer + sampling_width

    # Simulating centers and points for each genetic cluster
    centers = rng.uniform(0, area_width, (n_clust,2))
    samples = np.repeat(centers, n_obs, axis=0)
    radii = truncnorm.rvs(-2, 2, loc=0, scale=r, size = n_clust*n_obs, random_state=seed)
    angles = rng.uniform(0, 2*np.pi, n_clust*n_obs)
    samples[:,0] = samples[:,0] + radii*np.cos(angles)
    samples[:,1] = samples[:,1] + radii*np.sin(angles)
    samples = np.c_[samples, np.repeat(range(1,n_clust+1), n_obs)] # Genetic cluster of each simulated point
    inds = np.where(((samples[:,0] >= buffer) & (samples[:,0] <= (buffer+sampling_width))) & ((samples[:,1] >= buffer) & (samples[:,1] <= (buffer+sampling_width))))
    
    if len(inds[0]) == 0:
        return(np.repeat(5e3, 7))
    
    samples = samples[inds[0], :] # Keep only samples within the sampling area
    samples = np.c_[samples, 8*((samples[:,1] - buffer)//grid_width) + (samples[:,0] - buffer)//grid_width]
    observed_coords = np.unique(samples[:,3])
    
    # Compute results for each cell in the sampling area
    pool = mp.Pool(n_cpus)
    res = pool.map_async(searchNNClusts, [(samples[np.isin(samples[:,3], mooreNeighbourhood(divmod(cell, 8))).nonzero()[0]], cell) for cell in observed_coords])
    res.wait()
    results = sum(res.get())
    pool.close()
    
    # Aggregate results and compute final summary statistics
    tau_obs = []
    for i in range(len(results)):
        if results[i,1] == 0:
            tau_obs.append(5e3)
        else:
            tau_obs.append(results[i,0]/results[i,1])

    return(np.concatenate([tau_obs, [n_clust, r, len(np.unique(samples[:,2]))]]))
