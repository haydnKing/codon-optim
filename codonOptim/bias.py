import numpy as np, pandas as pd, itertools, random
import Bio.SeqIO as SeqIO

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

class Bias:
	"""Represent the codon bias in a genome or collection of genes"""
	def __init__(self, sr):
		"""Calculate codon bias.
		sr: genome seqrecord"""

		self._stat_string = ''
		self._bias = self._from_genome(sr)
		self._normed = self._normalise()

	def raw_bias(self):
		return self._bias

	def norm_bias(self):
		return self._normed

	def stats(self):
		return self._stat_string

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

	def _from_genome(self, sr):
		CDS = [f for f in sr.features if f.type == 'CDS']

		d = {}
		for cdn in self._list_codons():
			d[cdn] = 0

		for cds in CDS:
			seq = _extract(sr, cds)
			for cdn in (seq[i:i+3] for i in range(0, len(seq), 3)):
				if len(cdn) == 3:
					d[cdn] = d[cdn] + 1

		self._stat_string = ("Codon bias table build from {} CDSs, {} codons. " +
												 "Average CDS length {:.1f} codons").format(
														len(CDS),
														sum(d.values()),
														sum(d.values())/float(len(CDS)))

		return pd.Series(d)

	def _list_codons(self):
		return [''.join(x) for x in itertools.product('ATGC', repeat=3)]

	def _normalise(self):
		d = {}
		for aa, cdn_list in codon_table.iteritems():
			total = float(self._bias[cdn_list].sum())
			for cdn in cdn_list:
				d[cdn] = self._bias[cdn] / total
		return pd.Series(d)


	def __str__(self):
		s = ""
		for aa in AA:
			s += "{}:\n".format(aa)
			for cdn in codon_table[aa]:
				s += "\t\'{}\': {} ({:.1f}%)\n".format(cdn, self._bias[cdn], self._normed[cdn]*100.)
		return s


