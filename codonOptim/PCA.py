
from sklearn.decomposition import PCA
from sklearn import mixture
import matplotlib.pyplot as plt
import util
import numpy as np, pandas as pd


"""Code for performing PCA analysis"""


class PrincipalComponentAnalysis:

	def __init__(self, data, n_components=5, prior_weight = 1.0, labels=None):
		"""Initialise and perform the PCA using n_components
				labels is a dictionary mapping a string label to a list of data indexes to which they refer"""

		self._prior_weight = prior_weight
		self._labels = labels
		
		#normalise the data
		self._prior = data.sum(0)
		self._mean = None
		self._data = self._normalise(data.copy())

		#perform the PCA
		self._pca = PCA(n_components=n_components)
		self._pca.fit(self._data)

		#get PCA scores
		self._scores = self._get_scores(self._data)



	@classmethod
	def from_GenomeStats(cls, gs, n_components=5, prior_weight=1.0):
		return cls(gs.fo(), 
							 n_components = n_components, 
							 prior_weight = prior_weight, 
							 labels = {gs.name(): gs.fo().index,},)

	@classmethod
	def from_GMM(cls, 
							 data, 
							 K=3, 
							 PCA_components=3, 
							 prior_weight=1., 
							 rename_classes=True, 
							 covariance_type='diag'):
		"""Perform PCA and classify resulting scores using GMM/EM"""

		#initiate the object and perform PCA as normal
		ret = cls(data, PCA_components, prior_weight)

		#perform EM
		gmm = mixture.GMM(n_components=K, covariance_type=covariance_type)
		gmm.fit(ret._scores)

		#keep hold of the GMM object
		ret._gmm = gmm

		#predict classes
		classes = gmm.predict(ret._scores)

		#generate labels
		labels = {}
		for i in range(K):
			labels[i] = []
		for i,cls in enumerate(classes):
			labels[cls].append(data.index[i])
		for i in range(K):
			labels["Class {}".format(i+1)] = labels.pop(i)

		#set the labels
		ret._labels = labels

		if rename_classes:
			ax = ret.plot()
			ax.figure.show()
			for i in range(K):
				n = input("Enter new name for class {}: ".format(i+1))
				ret._labels[n] = ret._labels.pop("Class {}".format(i+1))

		return ret

	def add_sequence(self, label, seq):
		if type(seq) is str:
			bias = util.get_bias(seq)
		else:
			bias = util.get_bias(str(seq.seq))
		original = self._data.loc[112]
		
		#make sure the ordering is the same as the original
		bias = bias.reindex(self._data.columns)

		#normalise
		bias = self._normalise(bias)

		#get score
		self._scores.loc[label] = np.dot(self._pca.components_, bias)

		#set label
		self._labels[label] = [label,]
		

	def labels(self):
		return self._labels

	def _get_scores(self, data):

		scores = pd.DataFrame(np.zeros((len(data), 
																	 self._pca.n_components)),
													index=data.index)

		for i,r in data.iterrows():
			scores.loc[i,:] = np.dot(self._pca.components_, r)

		return scores


	def _normalise(self, data):
		
		#convert to dataframe
		itype = type(data)
		if itype is pd.Series:
			data = data.to_frame().transpose()


		#weight the prior according to prior_weight
		prior = self._prior.copy()
		for aa, codon_list in util.codon_table.items():
			prior[codon_list] = self._prior_weight*self._prior[codon_list] / float(prior[codon_list].sum())

		#normalise by codon and codon box number
		for aa, codon_list in util.codon_table.items():
			if aa == '*':
				continue
			n = data[codon_list].add(prior[codon_list], axis=1) 
			m = n.sum(axis=1) / float(len(codon_list))

			data[codon_list] = n.div(m, axis=0)

		if self._mean is None:
			self._mean = data.mean(0)

		#if I was given a series, return a series
		if itype is pd.Series:
			return (data - self._mean).iloc[0]
		return data - self._mean

	def indexes_from_circle(self, x, y, radius):
		"""return indexes of all scores that are within the circle"""
		distances = np.sqrt(np.square(self._scores[0] - x) + 
												np.square(self._scores[1] - y))
		return distances[distances <= radius].index

	def label_from_circle(self, name, x, y, radius):
		self._labels[name] = self.indexes_from_circle(x,y,radius)
		return self._labels[name]

	def plot(self, x=0, y=1, ax=None, order=None, colors=None):
		
		if not ax:
			ax = plt.figure().gca()

		if not order:
			order = list(self._labels.keys())

		cmap = plt.get_cmap()
		if not colors:
			colors = [cmap(i/float(len(order))) for i in range(len(order))]

		handles = ([ax.scatter(self._scores.loc[self._labels[k], x], 
													 self._scores.loc[self._labels[k], y], 
													 color=colors[i%len(colors)])
								for i,k in enumerate(order)])

		# Shrink current axis by 20%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

		# Put a legend to the right of the current axis
		ax.legend(handles, 
							order, 
							loc='center left', 
							bbox_to_anchor=(1, 0.5))

		ax.set_xlabel('$Z_{{{}}}$'.format(x+1))
		ax.set_ylabel('$Z_{{{}}}$'.format(y+1))

		return ax
