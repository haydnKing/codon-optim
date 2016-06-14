import numpy as np, pandas as pd, itertools, random
import Bio.SeqIO as SeqIO
import matplotlib.pyplot as plt
import os.path
import util

random.seed()


_complement = {'A':'T','T':'A','C':'G','G':'C'}

def _extract(genome, f):
	"""a zillion times faster than biopython"""
	if f.location_operator:
		if f.location_operator != 'join':
			print(("Unsupported location_operator \'{}\', "+
						 "ignoring feature \'{}\'").format(f.location_operator, 
																							 f.qualifiers['gene']))
		locs = f.location.parts
	else:
		locs = [f.location,]
	seq = ''
	for l in locs:
		s = str(genome.seq[int(l.start):int(l.end)]).upper()
		if l.strand < 0:
			s = ''.join(list(reversed([_complement[b] for b in s])))
		seq = seq + s


	if len(seq) == 0:
		print("len(seq) = 0")
		print("locs = {}".format(locs))
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
		for aa, cdn_list in util.codon_table.items():
			out[cdn_list] = bias[cdn_list] / float(bias[cdn_list].sum())
		return out

def so_normalise(so):
	out = pd.DataFrame(np.zeros((64,64)),
										 index=so.index,
										 columns=so.columns)
	for i in out.index:
		out.loc[i,:] = normalise(so.loc[i,:])

	return out

def score(bias, seq):
	r = 0
	seq = str(seq).upper()
	for cdn in [seq[i:i+3] for i in range(0,len(seq),3)]:
		if len(cdn) == 3:
			r = r + np.log(bias[cdn])

	return r / (len(seq)/3)

def so_score(so_weight, seq):
	r = 0
	seq = str(seq).upper()
	for cdnA,cdnB in [(seq[i-3:i],seq[i:i+3]) for i in range(3,len(seq),3)]:
		if len(cdnB) == 3:
			r = r + np.log(so_weight.at[cdnA,cdnB])

	return r / (len(seq)/3)

class GenomeStats:
	"""Represent statistical information about a genome or collection of genes"""
	def __init__(self, _name, _data, _second_order, _scores, _ndata=None, _nso=None):
		"""Calculate codon bias in CDS annotations.
		sr: genome seqrecord"""

		self._name = _name
		self._data = _data
		self._second_order = _second_order
		self._scores = _scores
		self._bias = self._data.sum(0)
		if _ndata is None:
			self._normed = normalise(self._bias)
		else:
			self._normed = _ndata
		if _nso is None:
			self._so_normed = so_normalise(self._second_order)
		else:
			self._so_normed = _nso
		
	@classmethod
	def from_seqrecord(cls, sr, featuretype='CDS', name=None):

		if not name:
			name = sr.name

		CDS = [f for f in sr.features if f.type == featuretype]

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

		#calculate scores
		_nd = normalise(_data.sum(0))
		_nso= so_normalise(_second_order)
		for i,seq in enumerate(_seqs):
			_scores.at[i,'first'] = score(_nd, seq)
			_scores.at[i,'second'] = so_score(_nso, seq)

		return cls(name, _data, _second_order, _scores, _nd, _nso)

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

	def fo(self):
		return self._data

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

	def score(self, seq):
		return score(self._normed, seq)

	def so_score(self, seq):
		return so_score(self._so_normed, seq)

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

	def __str__(self):
		s = ["{}: {:,} CDSs ({:,} codons)".format(self._name, 
																							len(self._data),
																							np.sum(self._bias)),
				 "fields: [triplet] ([amino-acid]) [normalised frequency] ([count])",]
		cols = int(np.ceil(np.log10(np.max(self._bias))))
		#extra for commas
		cols = cols + int(cols/3)
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

