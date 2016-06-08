import numpy as np, pandas as pd, itertools, random
import Bio.SeqIO as SeqIO
import matplotlib.pyplot as plt
import os.path
from sklearn.decomposition import PCA
import util

random.seed()



def _extract(genome, f):
	"""a zillion times faster than biopython"""
	seq = str(genome.seq[int(f.location.start):int(f.location.end)]).upper()
	if f.location.strand < 0:
		complement = {'A':'T','T':'A','C':'G','G':'C'}
		seq = ''.join(list(reversed([complement[b] for b in seq])))
	return seq

def list_codons():
	return [''.join(x) for x in itertools.product('ATGC', repeat=3)]

def list_non_stop_codons():
	return [c for c in list_codons() if c not in util.codon_table['*']]

def get_bias(seq):
	d = {}
	for cdn in list_codons():
		d[cdn] = 0
	seq = str(seq).upper()
	for cdn in (seq[i:i+3] for i in range(0, len(seq), 3)):
		if len(cdn) == 3:
			d[cdn] = d[cdn] + 1

	return pd.Series(d)

def add_second_order(so, seq):
	seq = str(seq).upper()
	for cdnA,cdnB in ((seq[i-3:i], seq[i:i+3]) for i in range(3, len(seq), 3)):
		if len(cdnB) == 3:
			so.at[cdnA, cdnB] = so.at[cdnA, cdnB] + 1

def normalise(bias):
		out = pd.Series(np.zeros(len(bias)), index=bias.index)
		for aa, cdn_list in util.codon_table.iteritems():
			out[cdn_list] = bias[cdn_list] / float(bias[cdn_list].sum())
		return out

def so_normalise(so):
	out = pd.DataFrame(np.zeros((64,64)),
										 index=so.index,
										 columns=so.columns)
	for i in out.index:
		out.loc[i,:] = normalise(so.loc[i,:])

	return out

class GenomeStats:
	"""Represent statistical information about a genome or collection of genes"""
	def __init__(self, _name, _data, _second_order, _scores):
		"""Calculate codon bias in CDS annotations.
		sr: genome seqrecord"""

		self._name = _name
		self._data = _data
		self._second_order = _second_order
		self._scores = _scores
		self._bias = self._data.sum(0)
		self._normed = normalise(self._bias)
		self._so_normed = so_normalise(self._second_order)
		
	@classmethod
	def from_seqrecord(cls, sr, featuretype='CDS', name=None):

		if not name:
			name = sr.name

		CDS = [f for f in sr.features if f.type == featuretype]

		print(CDS[3589])

		_data = pd.DataFrame(np.zeros((len(CDS), 64), dtype=int), 
															columns = list_codons())

		_second_order = pd.DataFrame(np.zeros((64,64), dtype=int),
																			index = list_codons(),
																			columns = list_codons())

		_scores = pd.DataFrame(np.zeros((len(CDS), 2)), 
															columns = ['first', 'second',])
		
		_seqs = [_extract(sr, cds) for cds in CDS]
		for i,seq in enumerate(_seqs):
			_data.loc[i,:] = get_bias(seq)
			add_second_order(_second_order, seq)

		ret = cls(name, _data, _second_order, _scores)
		#calculate scores
		for i,seq in enumerate(_seqs):
			_scores.at[i,'first'] = ret.score(seq)
			_scores.at[i,'second'] = ret.so_score(seq)

		ret._scores = _scores
		return ret

	def save_stats(self, folder, name=None):
		if not name:
			name = self._name

		self._data.to_csv(os.path.join(folder, name + ".1.csv"))
		self._second_order.to_csv(os.path.join(folder, name + ".2.csv"))
		self._scores.to_csv(os.path.join(folder, name + ".scores.csv"))

	
	@classmethod
	def from_file(cls, folder, name = None):
		if not name:
			folder, name = os.path.split(folder)

		_data   = pd.read_csv(os.path.join(folder, name + ".1.csv"), index_col=0)
		_so     = pd.read_csv(os.path.join(folder, name + ".2.csv"), index_col=0)
		_scores = pd.read_csv(os.path.join(folder, name + ".scores.csv"), index_col=0)

		return cls(name, _data, _so, _scores)


	def raw_bias(self):
		return self._bias

	def norm_bias(self):
		return self._normed

	def so(self):
		return self._second_order

	def so_normed(self):
		return self._so_normed

	def name(self):
		return self._name

	def emit(self, AAseq, cutoff=0.):
		"""Emit a codon to code for the amino acid 'aa'. Ignore codons with 
		frequency below cutoff (e.g. 0.1 = 10%)
		"""
		values = {}
		for aa in util.AA:
			cdn_list = [c for c in util.codon_table[aa] if self._normed[c] > cutoff]
			if not cdn_list:
				raise ValueError("No possible codons for aa=\'{}\'. Try a lower cutoff".format(aa))

			#emit a codon with probability equal to the remains
			values[aa] = self._bias[cdn_list] / float(self._bias[cdn_list].sum())

		oseq = []
		for aa in AAseq:
			r = random.random()
			for cdn,s in values[aa].cumsum().iterkv():
				if r < s:
					oseq.append(cdn)
					break
		return ''.join(oseq)

	def generate_codons(self, AAseq, bias=None, cutoff=0.):
		"""Generate unordered codons to use for each amino acid present in AAseq 
		such that the codon usage is as close as possible to bias."""
		if not bias:
			bias = self._normed
		else:
			bias = normalise(bias)

		out = {}
		for aa in util.AA:
			#list all codons which could be used for this aa
			cdn_list = [c for c in util.codon_table[aa] if bias[c] > cutoff]
			#how many codons do we need for this aa?
			count = len([1 for aas in AAseq if aas == aa])
			#what number of each codon should we have?
			counts = (bias[cdn_list] / np.sum(bias[cdn_list]))*count
			#sort by smallest residual
			counts = pd.DataFrame({'c':counts, 'r':np.abs(counts-np.around(counts))}).sort('r')['c']
			#assign integers
			overflow = 0.
			icounts = pd.Series(np.zeros(len(counts), dtype=int), index=counts.index)
			for i in range(len(counts)):
				icounts[i] = int(np.round(counts[i]+overflow))
				overflow = counts[i] - icounts[i]
			#list of codons
			out[aa] = []
			for cdn,count in icounts.iterkv():
				out[aa] = out[aa] + [cdn,]*count
			#shuffle the list (in some schemes, the codons are taken in list order
			#when the genome lacks information)
			random.shuffle(out[aa])

		return out

	def score(self, seq, bias=None):
		if bias is None:
			bias = self._normed

		r = 0
		seq = str(seq).upper()
		for cdn in [seq[i:i+3] for i in range(0,len(seq),3)]:
			if len(cdn) == 3:
				r = r + np.log(bias[cdn])

		return r / (len(seq)/3)

	def so_score(self, seq):
		r = 0
		seq = str(seq).upper()
		for cdnA,cdnB in [(seq[i-3:i],seq[i:i+3]) for i in range(3,len(seq),3)]:
			if len(cdnB) == 3:
				r = r + np.log(self._so_normed.at[cdnA,cdnB])

		return r / (len(seq)/3)

	def plot_score(self, names=[], seqs=[], colors=None, order=1):
		if order == 1:
			scores = self._scores['first']
		elif order == 2:
			scores = self._scores['second']
		else:
			raise ValueError("order must be 1 or 2, not {}".format(order))

		ax = plt.figure().gca()

		if colors is None:
			cmap = plt.get_cmap()
			colors = [cmap(i/float(len(names))) for i in range(len(names))]

		#histogram of the scores
		ax.hist(scores, bins=20, label = 'background', color='gray')

		for name,seq,color in zip(names,seqs,colors):
			if order == 1:
				score = self.score(seq)
			elif order == 2:
				score = self.so_score(seq)
			ax.axvline(score, color=color, label=name)

		ax.legend(loc=0)

		if order == 1:
			ax.set_xlabel("First Order Average log-probability")
		elif order == 2:
			ax.set_xlabel("Second Order Average log-probability")

		return ax

	def get_pca_scores(self, prior_weight, components):
		data = self._data.copy()
		pca = do_pca(data, components=components, prior_weight=prior_weight)
		return get_pca_scores(pca, data)

	def plot_pca(self, names=[], sequences=[], colors=None, x=0, y=1, prior_weight=0.5, ax=None):
		if len(names) != len(sequences):
			raise ValueError("Should have the same number of names and sequences")
		data = self._data.copy()
		for n,s in zip(names, sequences):
			data.loc[n,:] = get_bias(s.seq)

		if not colors:
			cmap = plt.get_cmap()
			colors = [cmap(i/float(len(names))) for i in range(len(names))]
		
		pca = do_pca(data, components=max(x,y)+1, prior_weight=prior_weight)
		scores = get_pca_scores(pca, data)

		if not ax:
			ax = plt.figure().gca()
		handles = ([ax.scatter(scores[x], scores[y], color='grey'),] +
						   [ax.scatter(scores.loc[n, x], 
													 scores.loc[n, y], 
													 color=colors[i%len(colors)])
								for i,n in enumerate(names)])

		# Shrink current axis by 20%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

		# Put a legend to the right of the current axis
		ax.legend(handles, [self._name,] + names, 
							loc='center left', 
							bbox_to_anchor=(1, 0.5))

		ax.set_xlabel('$Z_{{{}}}$'.format(x+1))
		ax.set_ylabel('$Z_{{{}}}$'.format(y+1))

		return ax



	def __str__(self):
		s = ["{}: {:,} CDSs ({:,} codons)".format(self._name, 
																							len(self._data),
																							np.sum(self._bias)),
				 "fields: [triplet] ([amino-acid]) [normalised frequency] ([count])",]
		cols = int(np.ceil(np.log10(np.max(self._bias))))
		cols = cols + cols/3
		fmt = ("{} ({}) {:2.2f} ({:"+str(cols)+",d})")

		for a in ['T', 'C', 'A', 'G',]:
			for c in ['T', 'C', 'A', 'G',]:
				line = []
				for b in ['T', 'C', 'A', 'G',]:
					cdn = a+b+c
					aa = util.inv_codon_table[cdn]
					line.append(fmt.format(cdn, aa, self._normed[cdn], self._bias[cdn]))

				s.append('  '.join(line)) 
			s.append('')
		return '\n'.join(s[:-1])

def do_pca(data, components=5, prior_weight=20):

	data = pca_normalise(data, prior_weight)
	pca = PCA(n_components=components)
	pca.fit(data)

	return pca

def get_pca_scores(pca, data):

	scores = pd.DataFrame(np.zeros((len(data), pca.n_components)),
												index=data.index)

	for i,r in data.iterrows():
		scores.loc[i,:] = np.dot(pca.components_, r)

	return scores

def pca_normalise(data, prior_weight=20.):

	prior = data.sum(0)
	for aa, codon_list in util.codon_table.iteritems():
		prior[codon_list] = prior_weight*prior[codon_list] / float(prior[codon_list].sum())

	for aa, codon_list in util.codon_table.iteritems():
		if aa == '*':
			continue
		n = data[codon_list].add(prior[codon_list], axis=1) 
		m = n.sum(axis=1) / float(len(codon_list))
		data[codon_list] = n.div(m, axis=0)

	mean = data.mean(0)

	return data - mean
