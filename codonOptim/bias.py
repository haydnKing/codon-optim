import numpy as np, pandas as pd
import Bio.SeqIO as SeqIO

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


def _extract(genome, feature):
	"""a zillion times faster than biopython"""
	complement = {'A':'T','T':'A','C':'G','G':'C'}
	seq = str(genome.Seq[int(gene.location.start):int(gene.location.end)]).upper()
	if gene.location.strand < 0:
		seq = ''.join(list(reversed([complement[b] for b in seq])))
	return seq

class Bias:
	"""Represent the codon bias in a genome or collection of genes"""
	def __init__(self, sr):
		"""Calculate codon bias.
		sr: genome seqrecord"""

		self._bias = self._from_genome(sr)
		self._normed = self._normalise()

	def raw_bias(self):
		return self._bias

	def norm_bias(self):
		return self._normed

	def _from_genome(self, sr):
		genes = [f for f in sr.features if 
				(f.type == 'CDS' and f.qualifiers.has_key('gene'))]
		
		d = {}
		for cdn in self.list_codons():
			d[cdn] = 0

		for gene in genes:
			seq = _extract(sr, gene)
			for cdn in (seq[i:i+3] for i in range(0, len(seq), 3)):
				d[cdn] = d[cdn] + 1

		return pd.Series(d)

	def _list_codons(self):
		return [''.join(x) for x in itertools.product('ATGC', repeat=3)]

	def _normalised(self):
		d = {}
		for aa, cdn_list in codon_table.iteritems():
			total = float(self._bias[codon_list].sum())
			for cdn in cdn_list:
				d[cdn] = self._bias[cdn] / total
		return pd.Series(d)
