
from sklearn.decomposition import PCA
from sklearn import mixture
import matplotlib.pyplot as plt


"""Code for performing PCA analysis"""


class PrincipalComponentAnalysis:

	def __init__(self, data, n_components=5, prior_weight = 1.0, labels=None):
		"""Initialise and perform the PCA using n_components
				labels is a dictionary mapping a string label to a list of data indexes to which they refer"""

		self._prior_weight = prior_weight
		self._labels = labels
		
		#normalise the data
		self._prior = data.sum(0)
		self._data = self._normalise(data)

		#perform the PCA
		pca = PCA(n_components=n_components)
		pca.fit(self._data)

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
			labels[cls].append(ret._scores.index[i])
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

	def _get_scores(self, data):

		scores = pd.DataFrame(np.zeros((len(data), 
																	 self._pca.n_components)),
													index=data.index)

		for i,r in data.iterrows():
			scores.loc[i,:] = np.dot(pca.components_, r)

		return scores


	def _normalise(self, data):
		#weight the prior according to prior_weight
		prior = self._prior.copy()
		for aa, codon_list in util.codon_table.iteritems():
			prior[codon_list] = self._prior_weight*self._prior[codon_list] / float(prior[codon_list].sum())

		for aa, codon_list in util.codon_table.iteritems():
			if aa == '*':
				continue
			n = data[codon_list].add(prior[codon_list], axis=1) 
			m = n.sum(axis=1) / float(len(codon_list))
			data[codon_list] = n.div(m, axis=0)

		mean = data.mean(0)

		return data - mean

	def plot(self, x=0, y=1, ax=None):

		cmap = plt.get_cmap()
		colors = [cmap(i/float(len(self._labels))) for i in range(len(self._labels))]
		
		if not ax:
			ax = plt.figure().gca()

		handles = ([ax.scatter(scores.loc[n, x], 
													 scores.loc[n, y], 
													 color=colors[i%len(colors)])
								for i,n in enumerate(self._labels.values())])

		# Shrink current axis by 20%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

		# Put a legend to the right of the current axis
		ax.legend(handles, self._labels.keys(), 
							loc='center left', 
							bbox_to_anchor=(1, 0.5))

		ax.set_xlabel('$Z_{{{}}}$'.format(x+1))
		ax.set_ylabel('$Z_{{{}}}$'.format(y+1))

		return ax
