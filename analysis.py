import numpy as np
import matplotlib.pyplot as plt
import random
import sys

import skbio
from skbio.stats.distance import mantel
from igraph import *

############################### I/O functions for Ising model predictions ###############################

def select_few_OTUs():
	data, cluster_names = read_data('Data/uclust.otus.97.mtx', normalize = True, threshold = 10000, binary_threshold = True)
	with open('Data/bacteria_data_pared.csv', 'w') as f:
		for i in xrange(data.shape[0]):
			f.write(','.join([str(j) for j in data[i]]) + '\n')

def plot_pred_vs_act_probs():
	pred_probs = [float(val) for val in open('./sample_probabilities.csv', 'r').readlines()]
	actual_data = open('Data/bacteria_data_pared.csv', 'r').readlines()
	bit_seqs = {} # str sequence : (count, first index seen)
	for i, line in enumerate(actual_data):
		curr_val = ''.join([str(bit.split('.')[0]) for bit in line.split(',')])
		if curr_val not in bit_seqs:
			bit_seqs[curr_val] = [0.0, i]
		bit_seqs[curr_val][0] += (1.0 / len(actual_data))
	print bit_seqs

	scatter_data = np.zeros((len(bit_seqs),2))
	for i, bits in enumerate(bit_seqs):
		scatter_data[i, 0] = bit_seqs[bits][0]
		scatter_data[i, 1] = pred_probs[bit_seqs[bits][1]]
	plt.scatter(scatter_data[:,0], scatter_data[:,1])
	plt.title('Predicted versus Actual Colony Patterns')
	plt.xlabel('Actual')
	plt.ylabel('Predicted')
	plt.show()

################################ Mantel Test (and related analysis) ###################

def mantel_test(m1, m2):
	return mantel(m1, m2)

def plot_corrs_vs_dists(corrs, dists):
	''' Plots correlation coefficients versus OTU distances '''
	# Can't plot all points, there are too many.  Take a random sample
	random_indices = random.sample(range(corrs.shape[0] * corrs.shape[1]), 10000)
	corrs = np.ndarray.flatten(corrs)[random_indices]
	dists = np.ndarray.flatten(dists)[random_indices]
	plt.scatter(dists, corrs)
	plt.title('Correlation vs OTU Distance')
	plt.xlabel('Distance')
	plt.ylabel('OTU Distance')
	plt.show()

################################ Correlation Analysis Graphs ###########################

def generate_graph(strong_corrs, cluster_names, filename):
	''' Displays correlation graph using igraph library '''
	edges = []
	edge_weights = []
	for i in xrange(strong_corrs.shape[0]):
		for j in xrange(i + 1, strong_corrs.shape[0]):
			if strong_corrs[i,j] > 0:
				edges.append((i,j))
				edge_weights.append(strong_corrs[i,j])
	g = Graph(edges)

	# Name vertices with OTU cluster names
	g.vs['cluster'] = cluster_names
	g.es['weight'] = edge_weights

	layout = g.layout('kk')
	visual_style = {}
	visual_style["vertex_size"] = 10
	visual_style["vertex_label"] = g.vs['cluster']
	visual_style["edge_width"] = [1 + int(2 * ew) for ew in g.es['weight']]
	visual_style["layout"] = layout
	visual_style["bbox"] = (3000, 3000)
	visual_style["margin"] = 20

	plot(g, filename, **visual_style)

def hist_corrs(corrs):
	# Get only half the matrix, not including the center
	flattened_corrs = []
	for i in xrange(corrs.shape[0]):
		for j in xrange(i + 1, corrs.shape[0]):
			flattened_corrs.append(corrs[i,j])
	plt.hist(flattened_corrs, 100)
	plt.title('Histogram of Pairwise OTU Correlations')
	plt.xlabel("Pearson's Correlation Coefficient")
	plt.ylabel('Frequency')
	plt.show()

def hist_OTU_counts(data, write_to_file = False):
	# Get sums of OTUs (over all 111 samples)
	sums = np.log10(np.sum(data, axis = 0))

	if write_to_file:
		with open('OTU_summary.txt', 'w') as f:
			for i in xrange(len(sums)):
				f.write(str(sums[i]) + '\n')

	plt.hist(sums, 150)
	plt.title('Histogram of OTU Counts')
	plt.xlabel('log(OTU count) before normalization')
	plt.ylabel('Frequency')
	plt.show()


##################################### Data Fetching Functions ##################################

def get_correlation_matrix(data, save_data = False):
	''' Returns the correlation coefficient matrix between OTUs.
	The data is assumed to be organized as columns representing OTUs and rows
	representing participant samples. '''
	corrs = np.corrcoef(data.transpose())
	if save_data:
		np.savetxt('OTU_correlations.txt', corrs, fmt='%1.2f', delimiter = ',')
	return corrs

def get_strong_corrs(corrs, threshold = 0.3):
	''' Thresholds correlation matrix for strong correlations. '''
	strong_corrs = np.array(corrs)
	strong_corrs[np.where(np.abs(strong_corrs) < threshold)] = 0.0
	return strong_corrs

def check_phylo_class(dist_matrix):
	raw_data = open('Data/uclust.otus.97.tax', 'r').readlines()
	clusters = {}
	tax_names = {}
	tax_count = 0
	i = 0
	level = 5
	for line in raw_data[2:]:
		taxonomy = line.split()[1]
		taxonomy = taxonomy.lower()
		if taxonomy != 'Unknown' and taxonomy.count(';') > level - 1:
			taxonomy = taxonomy.split(';')[level] # Get class
			if taxonomy not in clusters:
				clusters[taxonomy] = []
				tax_names[tax_count] = taxonomy
				tax_count += 1
			clusters[taxonomy].append(i)
		i += 1

	pruned_clusters = {}
	pruned_tax_names = {}
	pruned_count = 0
	for name in clusters:
		if len(clusters[name]) > 30:
			print len(clusters[name])
			pruned_clusters[name] = clusters[name]
			pruned_tax_names[pruned_count] = name
			pruned_count += 1

	# Sample and print indices within each cluster and between clusters
	num_intra_samples = 1000
	num_inter_samples = 1000
	intra_cluster_dists = np.zeros((num_intra_samples, 1))
	inter_cluster_dists = np.zeros((num_inter_samples, 1))
	for i in xrange(num_intra_samples):
		# Choose cluster, than choose two random OTUs in that cluster
		cluster_name = pruned_tax_names[random.randint(0, len(pruned_clusters) - 1)]
		rn1 = random.randint(0, len(pruned_clusters[cluster_name]) - 1)
		rn2 = random.randint(0, len(pruned_clusters[cluster_name]) - 1)
		while rn2 == rn1:
			rn2 = random.randint(0, len(pruned_clusters[cluster_name]) - 1)
		cluster_i1 = pruned_clusters[cluster_name][rn1]
		cluster_i2 = pruned_clusters[cluster_name][rn2]
		intra_cluster_dists[i] = dist_matrix[cluster_i1, cluster_i2]
		
	for i in xrange(num_inter_samples):
		# Choose two random clusters, then choose one random OTU within each cluster
		cluster_name1 = tax_names[random.randint(0, len(clusters) - 1)]
		cluster_name2 = tax_names[random.randint(0, len(clusters) - 1)]
		while cluster_name2 == cluster_name1:
			cluster_name2 = tax_names[random.randint(0, len(clusters) - 1)]
		rn1 = random.randint(0, len(clusters[cluster_name1]) - 1)
		rn2 = random.randint(0, len(clusters[cluster_name2]) - 1)
		cluster_i1 = clusters[cluster_name1][rn1]
		cluster_i2 = clusters[cluster_name2][rn2]
		inter_cluster_dists[i] = dist_matrix[cluster_i1, cluster_i2]

	fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
	axes[0].hist(intra_cluster_dists, bins=[i * 0.05 for i in xrange(8)], normed = True)
	axes[1].hist(inter_cluster_dists, bins=[i * 0.05 for i in xrange(8)], normed = True)
	plt.show()
	print np.mean(intra_cluster_dists), np.std(intra_cluster_dists), np.mean(inter_cluster_dists), np.std(inter_cluster_dists)


def check_phylo_phylum(dist_matrix):
	raw_data = open('Data/uclust.otus.97.tax', 'r').readlines()
	clusters = {}
	tax_names = {}
	tax_count = 0
	i = 0
	only_one_bac = set(['Crenarchaeota', 'OP3', 'Gemmatimonadetes', 'TM6', 'Chlamydiae', 'Mitochondria'])
	for line in raw_data[2:]:
		taxonomy = line.split()[1]
		if taxonomy != 'Unknown':
			print taxonomy, only_one_bac
			if ';' in taxonomy:
				taxonomy = taxonomy.split(';')[1]
			if taxonomy not in clusters and taxonomy not in only_one_bac:
				clusters[taxonomy] = []
				tax_names[tax_count] = taxonomy
				tax_count += 1
			if taxonomy not in only_one_bac:
				clusters[taxonomy].append(i)
		i += 1

	# Sample and print indices within each cluster and between clusters
	num_intra_samples = 1000
	num_inter_samples = 4000
	intra_cluster_dists = np.zeros((num_intra_samples, 1))
	inter_cluster_dists = np.zeros((num_inter_samples, 1))
	for i in xrange(num_intra_samples):
		# Choose cluster, than choose two random OTUs in that cluster
		cluster_name = tax_names[random.randint(0, len(clusters) - 1)]
		rn1 = random.randint(0, len(clusters[cluster_name]) - 1)
		rn2 = random.randint(0, len(clusters[cluster_name]) - 1)
		while rn2 == rn1:
			rn2 = random.randint(0, len(clusters[cluster_name]) - 1)
		cluster_i1 = clusters[cluster_name][rn1]
		cluster_i2 = clusters[cluster_name][rn2]
		intra_cluster_dists[i] = dist_matrix[cluster_i1, cluster_i2]
		
	for i in xrange(num_inter_samples):
		# Choose two random clusters, then choose one random OTU within each cluster
		cluster_name1 = tax_names[random.randint(0, len(clusters) - 1)]
		cluster_name2 = tax_names[random.randint(0, len(clusters) - 1)]
		while cluster_name2 == cluster_name1:
			cluster_name2 = tax_names[random.randint(0, len(clusters) - 1)]
		rn1 = random.randint(0, len(clusters[cluster_name1]) - 1)
		rn2 = random.randint(0, len(clusters[cluster_name2]) - 1)
		cluster_i1 = clusters[cluster_name1][rn1]
		cluster_i2 = clusters[cluster_name2][rn2]
		print cluster_i1, cluster_i2, dist_matrix[cluster_i1, cluster_i2]
		inter_cluster_dists[i] = dist_matrix[cluster_i1, cluster_i2]

	print np.mean(intra_cluster_dists), np.std(intra_cluster_dists), np.mean(inter_cluster_dists), np.std(inter_cluster_dists)
	
def get_OTU_dists(cluster_names):
	''' Returns a 2D matrix of OTU distances in the same order as the clusters in the correlation matrix '''
	dists = np.zeros((len(cluster_names),len(cluster_names)))
	raw_data = open('Data/uclust.seeds.97.dist', 'r').readlines()
	lt = index_lookup_table(cluster_names)
	cluster_names = set(cluster_names)
	for line in raw_data:
		curr_data = line.split()
		cluster1 = int(curr_data[0].split(';')[0][7:])
		cluster2 = int(curr_data[1].split(';')[0][7:])
		i1 = None
		i2 = None
		if cluster1 in cluster_names and cluster2 in cluster_names:
			i1 = lt[cluster1]
			i2 = lt[cluster2]
			dist = float(curr_data[2])
			dists[i1, i2] = dist
			dists[i2, i1] = dist
	return dists

def index_lookup_table(cluster_names):
	''' Used to look up the index for each cluster in the order it appears in the correlation matrix ''' 
	lookup_table = {}
	for i in xrange(len(cluster_names)):
		lookup_table[cluster_names[i]] = i
	return lookup_table

def read_data(filename, normalize = True, threshold = 0, binary_threshold = True, binary_threshold_value = 0.1):
	''' Reads the normalized data into a numpy array, thresholding as desired.

	    threshold = the raw threshold.  If the OTU value for a particular sample is less than this threshold,
		we threshold to 0.
        binary threshold = a threshold that represents presence or lack of presence of an OTU in a particular sample.
		This threshold is based on the normalized values, not the raw data values  '''
	data = []
	lines = open(filename, 'r').readlines()

	cluster_names = []
	for l in lines[1:]:
		string_data = l.split()
		curr_data = [float(val) for val in string_data[1:]]
		if sum(curr_data) >= threshold:
			cluster_names.append(int(string_data[0][7:]))
			data.append(curr_data)

	# Transpose and then normalize
	# Rows now represent each sample, with columns representing OTUs
	data = np.array(data).transpose()

	if normalize:
		sums = np.sum(data, axis = 1)
		for row_i in xrange(data.shape[0]):
			data[row_i] /= sums[row_i]

	if binary_threshold:
		for i in xrange(data.shape[0]):
			for j in xrange(data.shape[1]):
				if data[i][j] >= binary_threshold_value:
					data[i][j] = 1.0
				else:
					data[i][j] = 0.0

		# Remove columns with few '1's after binary thresholding
		nonzero_cols = np.where(np.sum(data, axis = 0) > 15)[0]
		for i in xrange(data.shape[1] - 1, -1, -1):
			if i not in nonzero_cols:
				del cluster_names[-1 * i]
		data = data[:, nonzero_cols]
	return data, cluster_names		

if __name__ == "__main__":
	job = None
	if len(sys.argv) > 1:
		job = int(sys.argv[1])

	# Check the OTU distances
	if job == 1:
		data, cluster_names = read_data('Data/uclust.otus.97.mtx', normalize = False, threshold = -1, binary_threshold = False)
		dist_matrix = get_OTU_dists(cluster_names)
		check_phylo_class(dist_matrix)

	# Write limited set of clusters for Ising model
	if job == 2:
		select_few_OTUs()

	# Show graph for Ising model results
	if job == 3:
		# Reads from file, so make sure files are correct
		plot_pred_vs_act_probs()

	# Mantel test
	if job == 4:
		thresh = 1000
		data, cluster_names = read_data('Data/uclust.otus.97.mtx', normalize = True, threshold = thresh, binary_threshold = False)
		dist_matrix = get_OTU_dists(cluster_names)
		corrs = get_correlation_matrix(data)
		print "Mantel test: ", mantel_test(1-corrs, dist_matrix)
		plot_corrs_vs_dists(corrs, dist_matrix)

	# Draw correlation graph
	if job == 5:
		data, cluster_names = read_data('Data/uclust.otus.97.mtx', normalize = False, threshold = -1, binary_threshold = False)
		corrs = get_correlation_matrix(data, save_data = False)
		strong_corrs = get_strong_corrs(corrs, threshold = 0.6)
		generate_graph(strong_corrs, cluster_names, 'temp.pdf')
