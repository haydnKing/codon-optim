import Bio.SeqIO as SeqIO
import pandas as pd, numpy as np
import itertools

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
for aa,cdn_list in codon_table.items():
	for cdn in cdn_list:
		inv_codon_table[cdn] = aa

def load_sequence(filename):

	try:
		return SeqIO.read(filename, "fasta")
	except ValueError:
		pass

	try:
		return SeqIO.read(filename, "genbank")
	except ValueError:
		pass
	
	raise ValueError("Could not load file \'{}\'.".format(filename))

def check_folder_exists(name):
	if not os.path.isdir(name):
		resp = ''
		while resp.upper() not in ['Y','N',]:
			resp = raw_input(("Output directory \'{}\' doesn't exist. "+
												"Create it? [Y/N]: ").format(args.output_folder))
		if resp.upper() == 'Y':
			os.mkdir(args.output_folder)
		else:
			return False
	return True

def normalise(codon_counts):
		out = pd.Series(np.zeros(len(codon_counts)), index=codon_counts.index)
		for aa, cdn_list in codon_table.items():
			out[cdn_list] = codon_counts[cdn_list] / float(codon_counts[cdn_list].sum())
		return out

def so_normalise(so):
	out = pd.DataFrame(np.zeros((64,64)),
										 index=so.index,
										 columns=so.columns)
	for i in out.index:
		out.loc[i,:] = normalise(so.loc[i,:])

	return out

def translate(seq):
	return [inv_codon_table[str(seq[i:i+3]).upper()] for i in range(0, len(seq), 3)]

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

def fo_to_string(name, data):
	totals = data.sum(0)
	normed = normalise(totals)
	s = ["{}: {:,} CDSs ({:,} codons)".format(name, 
																						len(data),
																						np.sum(totals)),
			 "fields: [triplet] ([amino-acid]) [normalised frequency] ([count])",]
	cols = int(np.ceil(np.log10(np.max(totals))))
	#extra for commas
	cols = cols + int(cols/3)
	fmt = ("{} ({}) {:2.2f} ({:"+str(cols)+",d})")

	for a in ['T', 'C', 'A', 'G',]:
		for c in ['T', 'C', 'A', 'G',]:
			line = []
			for b in ['T', 'C', 'A', 'G',]:
				cdn = a+b+c
				aa = inv_codon_table[cdn]
				line.append(fmt.format(cdn, aa, normed[cdn], totals[cdn]))

			s.append('  '.join(line)) 
		s.append('')
	return '\n'.join(s[:-1])
