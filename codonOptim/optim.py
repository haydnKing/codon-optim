"""Library for optimisations"""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.Alphabet as Alphabet
import numpy as np, pandas as pd
import util, PCA

def _verify(in_record, out_str, desc_tag=". codon optimised"):
	in_aa = util.translate(in_record.seq)
	out_aa= util.translate(out_str)

	if out_aa != in_aa:
		raise ValueError("Optimisation failed: translations don't match")

	return SeqRecord(Seq(out_str, Alphabet.generic_dna),
									id = in_record.id,
									name=in_record.name,
									description=in_record.description + desc_tag)

def exact(gs, sr, rare_codon_cutoff=0.):
	AA = util.translate(str(sr.seq))
	codons = _generate_codons(AA, gs.norm_bias(), cutoff=rare_codon_cutoff)

	oseq = []
	for aa in AA:
		cdn = codons[aa].pop(np.random.randint(0, len(codons[aa])))
		oseq.append(cdn)
	oseq = ''.join(oseq)

	return ('exact', _verify(sr, ret))

def simple(gs, sr, rare_codon_cutoff=0.):

	AA = util.translate(str(sr.seq))
	oseq = gs.emit(AA, rare_codon_cutoff)

	return ('simple', _verify(sr, oseq))

def second(gs, sr, rare_codon_cutoff=0., mode='rand'):

	if mode not in ['rand','maximum','minimum','irand']:
		raise ValueError("Unknown mode \'{}\'".format(mode))
	
	AAseq = util.translate(str(sr.seq))
	codons = _generate_codons(AAseq, gs.norm_bias(), cutoff=rare_codon_cutoff)
	oseq = _second(gs.so(), AAseq, codons, mode)
	return ('second.{}'.format(mode), _verify(sr, oseq))


def _second(so, AAseq, codons, mode='rand'):

	oseq = []
	for aa in AAseq:
		if not oseq:
			oseq.append(codons[aa].pop(np.random.randint(0, len(codons[aa]))))
		else:
			p = so.loc[oseq[-1], codons[aa]]
			#if there are no instances, make the distribution uniform
			if p.sum() == 0:
				p[:] = np.ones(len(p))
			p = p / float(p.sum())
			if mode in ['rand', 'irand']:
				r = np.random.random()
				if mode == 'irand' and len(p) > 1:
					p = 1-p
					p = p / p.sum()
				for i,s in enumerate(p.cumsum()):
					if r < s:
						oseq.append(codons[aa].pop(i))
						break
			elif mode == 'maximum':
				p.index = range(len(p))
				oseq.append(codons[aa].pop(np.argmax(p)))
			elif mode == 'minimum':
				p.index = range(len(p))
				oseq.append(codons[aa].pop(np.argmin(p)))

			
	return ''.join(oseq)


def second_PCA(gs, 
							 sr, 
							 rare_codon_cutoff=0., 
							 GMM_components=3, 
							 prior_weight=1., 
							 PCA_components=3,
							 mode='rand'):

	ret = []
	#Perform PCA and GMM/EM clustering
	pca = PCA.PrincipalComponentAnalysis.from_GMM(gs.fo(),
																								K=GMM_components,
																								PCA_components=PCA_components,
																								prior_weight=prior_weight)

	#Amino acid sequence
	AAseq = util.translate(str(sr.seq))

	#for each cluster
	for name, indexes in pca.labels().items():

		#calculate first order bias
		data = gs.fo().loc[indexes]

		print(util.fo_to_string(name, data))
		bias = util.normalise(data.sum(0))

		#generate codon lists given fo bias	
		codons = _generate_codons(AAseq, bias, cutoff=rare_codon_cutoff)

		#order codons according to whole genome so preference
		oseq = _second(gs.so(), AAseq, codons, mode)
		
		seq = _verify(sr, oseq)

		ret.append(('s.class_{}'.format(name),seq))
	
	return ret


def _generate_codons(AAseq, bias, cutoff=0.):
	"""Generate unordered codons to use for each amino acid present in AAseq 
	such that the codon usage is as close as possible to bias."""
	bias = util.normalise(bias)

	out = {}
	for aa in util.AA:
		#list all codons which could be used for this aa
		cdn_list = [c for c in util.codon_table[aa] if bias[c] > cutoff]
		#how many codons do we need for this aa?
		count = len([1 for aas in AAseq if aas == aa])
		#what number of each codon should we have?
		counts = (bias[cdn_list] / np.sum(bias[cdn_list]))*count
		#sort by smallest residual
		counts = pd.DataFrame({'c':counts, 
													 'r':np.abs(counts-np.around(counts))
													 }).sort_values(by='r')['c']
		#assign integers
		overflow = 0.
		icounts = pd.Series(np.zeros(len(counts), dtype=int), index=counts.index)
		for i in range(len(counts)):
			icounts[i] = int(np.round(counts[i]+overflow))
			overflow = counts[i] - icounts[i]
		#list of codons
		out[aa] = []
		for cdn,count in icounts.iteritems():
			out[aa] = out[aa] + [cdn,]*count
		#shuffle the list (in some schemes, the codons are taken in list order
		#when the genome lacks information)
		np.random.shuffle(out[aa])

	return out
