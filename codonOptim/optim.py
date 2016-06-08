"""Library for optimisations"""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.Alphabet as Alphabet
import util
import numpy as np

def _translate(seq):
	return [util.inv_codon_table[str(seq[i:i+3]).upper()] for i in range(0, len(seq), 3)]

def _verify(in_record, out_str, desc_tag=". codon optimised"):
	in_aa = _translate(in_record.seq)
	out_aa= _translate(out_str)

	if out_aa != in_aa:
		raise ValueError("Optimisation failed: translations don't match")

	return SeqRecord(Seq(out_str, Alphabet.generic_dna),
									id = in_record.id,
									name=in_record.name,
									description=in_record.description + desc_tag)

def exact(gs, sr, rare_codon_cutoff=0.):
	AA = _translate(str(sr.seq))
	codons = gs.generate_codons(AA, cutoff=rare_codon_cutoff)

	oseq = []
	for aa in AA:
		cdn = codons[aa].pop(np.random.randint(0, len(codons[aa])))
		oseq.append(cdn)
	oseq = ''.join(oseq)

	return _verify(sr, ret)

def simple(gs, sr, rare_codon_cutoff=0.):

	AA = _translate(str(sr.seq))
	oseq = gs.emit(AA, rare_codon_cutoff)

	return _verify(sr, oseq)

def second(gs, sr, rare_codon_cutoff=0., mode='rand'):

	if mode not in ['rand','maximum','minimum','irand']:
		raise ValueError("Unknown mode \'{}\'".format(mode))
	
	AA = _translate(str(sr.seq))
	codons = gs.generate_codons(AA, cutoff=rare_codon_cutoff)

	oseq = []
	for aa in AA:
		if not oseq:
			oseq.append(codons[aa].pop(np.random.randint(0, len(codons[aa]))))
		else:
			p = gs.so().loc[oseq[-1], codons[aa]]
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

			
	oseq = ''.join(oseq)

	return _verify(sr, oseq)
