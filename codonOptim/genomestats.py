import numpy as np, pandas as pd, itertools, random
import Bio.SeqIO as SeqIO
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

random.seed()

AA = ['F', 'L', 'S', 'Y', '*', 'C', 'W', 'P', 'H', 'Q', 'R', 'I', 'M', 'T', 'N', 'K', 'V', 'A', 'D', 'E', 'G']
codon_table = {
	'A': ['GCT', 'GCC', 'GCA', 'GCG'], 
	'C': ['TGT', 'TGC'], 
	'E': ['GAA', 'GAG'], 
	'D': ['GAT', 'GAC'], 
	'G': ['GGT', 'GGC', 'GGA', 'GGG'], 
	'F': ['TTT', 'TTC'], 
	'I': ['ATT', 'ATC', 'ATA'], 
	'H': ['CAT', 'CAC'], 
	'K': ['AAA', 'AAG'], 
	'*': ['TAA', 'TAG', 'TGA'], 
	'M': ['ATG'], 
	'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 
	'N': ['AAT', 'AAC'], 
	'Q': ['CAA', 'CAG'], 
	'P': ['CCT', 'CCC', 'CCA', 'CCG'], 
	'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 
	'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 
	'T': ['ACT', 'ACC', 'ACA', 'ACG'], 
	'W': ['TGG'], 
	'V': ['GTT', 'GTC', 'GTA', 'GTG'], 
	'Y': ['TAT', 'TAC']
}
inv_codon_table = {}
for aa,cdn_list in codon_table.iteritems():
	for cdn in cdn_list:
		inv_codon_table[cdn] = aa


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
	return [c for c in list_codons() if c not in codon_table['*']]

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
		for aa, cdn_list in codon_table.iteritems():
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
	def __init__(self, sr):
		"""Calculate codon bias in CDS annotations.
		sr: genome seqrecord"""

		CDS = [f for f in sr.features if f.type == 'CDS']

		self._data = pd.DataFrame(np.zeros((len(CDS), 64), dtype=int), 
															columns = list_codons())

		self._second_order = pd.DataFrame(np.zeros((64,64), dtype=int),
																			index = list_codons(),
																			columns = list_codons())
		
		self._seqs = [_extract(sr, cds) for cds in CDS]
		for i,seq in enumerate(self._seqs):
			self._data.loc[i,:] = get_bias(seq)
			add_second_order(self._second_order, seq)

		self._name = sr.name

		self._bias = self._data.sum(0)

		self._normed = normalise(self._bias)
		self._so_normed = so_normalise(self._second_order)

	def raw_bias(self):
		return self._bias

	def norm_bias(self):
		return self._normed

	def so(self):
		return self._second_order

	def so_normed(self):
		return self._so_normed

	def emit(self, AAseq, cutoff=0.):
		"""Emit a codon to code for the amino acid 'aa'. Ignore codons with 
		frequency below cutoff (e.g. 0.1 = 10%)
		"""
		values = {}
		for aa in AA:
			cdn_list = [c for c in codon_table[aa] if self._normed[c] > cutoff]
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
		for aa in AA:
			#list all codons which could be used for this aa
			cdn_list = [c for c in codon_table[aa] if bias[c] > cutoff]
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

	def plot_score(self, names, seqs, colors=None, order=1):
		if order not in [1,2]:
			raise ValueError("order must be 1 or 2, not {}".format(order))

		if order == 1:
			scores = [self.score(iseq) for iseq in self._seqs]
		elif order == 2:
			scores = [self.so_score(iseq) for iseq in self._seqs]

		ax = plt.figure().gca()

		if colors is None:
			cmap = plt.get_cmap()
			colors = [cmap(i/float(len(names))) for i in range(len(names))]

		ax.hist(scores, bins=20)
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
		

	def plot_pca(self, names=[], sequences=[], colors=None, x=0, y=1, prior_weight=0.5, ax=None):
		if len(names) != len(sequences):
			raise ValueError("Should have the same number of names and sequences")
		data = self._data.copy()
		for n,s in zip(names, sequences):
			data.loc[n,:] = get_bias(s.seq)

		if not colors:
			cmap = plt.get_cmap()
			colors = [cmap(i/float(len(names))) for i in range(len(names))]
		
		scores = do_pca(data, x, y, prior_weight)

		if not ax:
			ax = plt.figure().gca()
		handles = ([ax.scatter(scores['x'], scores['y'], color='grey'),] +
						   [ax.scatter(scores.loc[n,'x'], 
													 scores.loc[n, 'y'], 
													 color=colors[i%len(colors)])
								for i,n in enumerate(names)])

		# Shrink current axis by 20%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

		# Put a legend to the right of the current axis
		ax.legend(handles, ['genome',] + names, 
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
					aa = inv_codon_table[cdn]
					line.append(fmt.format(cdn, aa, self._normed[cdn], self._bias[cdn]))

				s.append('  '.join(line)) 
			s.append('')
		return '\n'.join(s[:-1])

def do_pca(data, x=0, y=1, prior_weight=20):

	data = pca_normalise(data, prior_weight)
	pca = PCA(n_components=max(x,y)+1)
	pca.fit(data)

	scores = pd.DataFrame(np.zeros((len(data), 2)),
												index=data.index,
												columns=['x','y'])

	for i,r in data.iterrows():
		scores.loc[i,'x'] = np.dot(r, pca.components_[x])
		scores.loc[i,'y'] = np.dot(r, pca.components_[y])

	return scores

def pca_normalise(data, prior_weight=20.):
	out = pd.DataFrame(np.zeros((len(data), 61)),
										 columns=list_non_stop_codons())
	prior = data.sum(0)
	for aa, codon_list in codon_table.iteritems():
		prior[codon_list] = prior[codon_list] / float(prior[codon_list].sum())

	for i,r in data.iterrows():
		for aa, codon_list in codon_table.iteritems():
			if aa == '*':
				continue
			total = (prior_weight*prior[codon_list].sum()+r[codon_list].sum()) / float(len(codon_list))
			out.loc[i,codon_list] = (prior_weight*prior[codon_list] + r[codon_list]) / total

	mean = out.mean(0)

	return out - mean
