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

class Bias:
	"""Represent the codon bias in a genome or collection of genes"""
	def __init__(self, sr):
		"""Calculate codon bias.
		sr: genome seqrecord"""

		self._from_genome(sr)

	def raw_bias(self):
		return self._bias

	def norm_bias(self):
		return self._normed

	def emit(self, aa, cutoff=0.):
		"""Emit a codon to code for the amino acid 'aa'. Ignore codons with 
		frequency below cutoff (e.g. 0.1 = 10%)
		"""
		cdn_list = [c for c in codon_table[aa] if self._normed[c] > cutoff]
		if not cdn_list:
			raise ValueError("No possible codons for aa=\'{}\'. Try a lower cutoff".format(aa))

		#emit a codon with probability equal to the remains
		values = self._bias[cdn_list] / float(self._bias[cdn_list].sum())
		r = random.random()
		for cdn,s in values.cumsum().iterkv():
			if r < s:
				return cdn

	def plot_pca(self, names=[], sequences=[], colors=None, x=0, y=1):
		if len(names) != len(sequences):
			raise ValueError("Should have the same number of names and sequences")
		data = self._data.copy()
		for n,s in zip(names, sequences):
			data.loc[n,:] = get_bias(s.seq)

		if not colors:
			cmap = plt.get_cmap()
			colors = [cmap(i/float(len(names))) for i in range(len(names))]
		
		scores = do_pca(data, x, y, 0.5)

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


	def _from_genome(self, sr):
		CDS = [f for f in sr.features if f.type == 'CDS']

		self._data = pd.DataFrame(np.zeros((len(CDS), 64), dtype=int), 
															columns = list_codons())

		for i,cds in enumerate(CDS):
			seq = _extract(sr, cds)
			self._data.loc[i,:] = get_bias(seq)

		self._name = sr.name

		self._bias = self._data.sum(0)
		self._normed = self._get_normed()

	def _get_normed(self):
		d = {}
		for aa, cdn_list in codon_table.iteritems():
			total = float(self._bias[cdn_list].sum())
			for cdn in cdn_list:
				d[cdn] = self._bias[cdn] / total
		return pd.Series(d)

	def __str__(self):
		s = ["{}: {:,} CDSs ({:,} codons)".format(self._name, 
																							len(self._data),
																							np.sum(self._bias)),
				 "fields: [triplet] [amino-acid] [normalised frequency] ([count])",]
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
